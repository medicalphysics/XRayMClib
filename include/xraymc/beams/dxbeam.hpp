/*This file is part of XRayMClib.

XRayMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

XRayMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with XRayMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <mutex>
#include <span>

namespace xraymc {

/**
 * @brief A single exposure of a diagnostic X-ray (DX) beam.
 *
 * Returned by `DXBeam::exposure()`. Emits photons from a fixed source position with
 * directions sampled uniformly over a symmetric rectangular collimation field using
 * `SphereSamplingRectangularField`. The local-frame direction is transformed to world
 * space via the supplied direction cosine basis. Energies are drawn from an owned copy
 * of a `SpecterDistribution`.
 *
 * Unlike `CBCTBeamExposure`, the spectrum is owned by value (copied at construction)
 * rather than referenced by pointer, so `DXBeamExposure` is self-contained.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class DXBeamExposure {
public:
    /**
     * @brief Constructs a DX exposure with full geometry and spectrum.
     *
     * Derives the central beam direction as `dircosines[0] × dircosines[1]` and builds a
     * `SphereSamplingRectangularField` for the supplied symmetric collimation half-angles.
     *
     * @param pos                   Source position [cm].
     * @param dircosines            Two orthonormal vectors `{cos_x, cos_y}` spanning the plane
     *                              perpendicular to the beam axis.
     * @param N                     Number of photon histories.
     * @param weight                Statistical weight of each sampled photon.
     * @param collimationHalfAngles `{half_x, half_y}` symmetric collimation half-angles [rad].
     * @param specter               Energy spectrum sampler (copied by value).
     */
    DXBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double> specter)
        : m_directionSampler(collimationHalfAngles)
        , m_pos(pos)
        , m_dirCosines(dircosines)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_NParticles(N)
        , m_weight(weight)
        , m_specter(specter)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const { return m_NParticles; }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * Draws a direction from `SphereSamplingRectangularField` in the local frame,
     * transforms it to world space via `particleDirection()`, then draws an energy
     * from the spectrum. Returns a photon with the configured weight.
     *
     * @param state  Per-thread PRNG state.
     * @return `Particle` or `ParticleTrack` depending on `ENABLETRACKING`.
     */
    auto sampleParticle(RandomState& state) const noexcept
    {
        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specter.sampleValue(state),
                .weight = m_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specter.sampleValue(state),
                .weight = m_weight
            };
            return p;
        }
    }

protected:
    /**
     * @brief Samples a world-space photon direction using the current direction cosine basis.
     *
     * Draws a local-frame direction from `SphereSamplingRectangularField` and transforms
     * it to world space via `changeBasis(cos_x, cos_y, dir, local_dir)`.
     *
     * @param state  Per-thread PRNG state.
     * @return World-space unit direction vector.
     */
    std::array<double, 3> particleDirection(RandomState& state) const
    {
        const auto dir = m_directionSampler(state);
        return vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dir);
    }

private:
    SphereSamplingRectangularField m_directionSampler;                                        ///< Samples directions within the symmetric rectangular collimation field.
    std::array<double, 3> m_pos = { 0, 0, 0 };                                               ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 };                                               ///< Central beam direction (derived from direction cosines).
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };   ///< Orthonormal basis perpendicular to m_dir.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };                                ///< Symmetric collimation half-angles {half_x, half_y} [rad].
    std::uint64_t m_NParticles = 100;                                                          ///< Number of photon histories.
    double m_weight = 1;                                                                       ///< Statistical weight of each photon.
    SpecterDistribution<double> m_specter;                                                     ///< Owned copy of the energy spectrum sampler.
};

/**
 * @brief A diagnostic X-ray (DX / radiography) beam, satisfying `BeamType`.
 *
 * Models a single fixed-geometry X-ray tube: photons are emitted from a fixed
 * source position with directions sampled over a symmetric rectangular collimation field
 * `{half_x, half_y}`. The beam orientation is specified by two direction cosines; the
 * central beam axis is their cross product. Photon energies are drawn from a
 * `SpecterDistribution` built from the internal `Tube` model.
 *
 * All exposures are identical (the index @p i passed to `exposure()` is ignored).
 *
 * `calibrationFactor()` converts accumulated energy scores to dose by dividing a
 * user-supplied measured DAP value by the expected air KERMA from the simulated fluence,
 * using the NIST dry-air mass energy-transfer coefficient.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class DXBeam {
public:
    /**
     * @brief Constructs a DX beam at the given position and orientation.
     *
     * @param pos                 Source position [cm]. Default: origin.
     * @param dircosines          Two orthonormal vectors `{cos_x, cos_y}`; normalised internally.
     *                            Default: `{{1,0,0}, {0,1,0}}` (beam along +z axis).
     * @param filtrationMaterials Map of `{atomic number Z → thickness [mm]}` added as filtration
     *                            to the internal `Tube`. Default: empty.
     * @param updateSpecter       If true (default), rebuild the spectrum cache after applying
     *                            filtration. Pass false to defer the rebuild.
     */
    DXBeam(
        const std::array<double, 3>& pos = { 0, 0, 0 },
        const std::array<std::array<double, 3>, 2>& dircosines = { { { 1, 0, 0 }, { 0, 1, 0 } } },
        const std::map<std::uint8_t, double>& filtrationMaterials = { }, bool updateSpecter = true)
        : m_pos(pos)
    {
        setDirectionCosines(dircosines);

        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);

        if (updateSpecter)
            tubeChanged();
    }

    /// @brief Returns the number of exposures.
    std::uint64_t numberOfExposures() const { return m_Nexposures; }

    /**
     * @brief Sets the number of exposures (clamped to ≥ 1).
     * @param n  Desired exposure count.
     */
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    /// @brief Returns the total number of photon histories across all exposures.
    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }

    /// @brief Returns the number of photon histories per exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    /**
     * @brief Sets the number of photon histories per exposure.
     * @param n  Histories per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /// @brief Returns the source position [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the source position [cm].
     * @param pos  3-D source position in world space [cm].
     */
    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }

    /// @brief Returns the two normalised direction cosines `{cos_x, cos_y}`.
    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

    /**
     * @brief Sets the beam orientation from a pair of direction cosine vectors; both are normalised.
     * @param dir  `{cos_x, cos_y}` — two orthogonal vectors spanning the beam-perpendicular plane.
     */
    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    /**
     * @brief Sets the beam orientation from two explicit direction-cosine vectors; both are normalised.
     * @param xdir  First in-plane direction cosine vector.
     * @param ydir  Second in-plane direction cosine vector, perpendicular to @p xdir.
     */
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dirCosines[0] = xdir;
        m_dirCosines[1] = ydir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    /// @brief Returns the measured dose-area product (DAP) used for dose calibration.
    double DAPvalue() const { return m_measuredDAP; }

    /**
     * @brief Sets the measured DAP value used for dose calibration; the absolute value is used.
     * @param dap  Measured dose-area product [Gy·cm²].
     */
    void setDAPvalue(double dap) { m_measuredDAP = std::abs(dap); }

    /// @brief Returns a const reference to the internal X-ray tube model.
    const Tube& tube() const { return m_tube; }

    /**
     * @brief Replaces the internal tube model; always rebuilds the spectrum cache.
     * @param tube  New `Tube` instance (copied).
     */
    void setTube(const Tube& tube)
    {
        m_tube = tube;
        tubeChanged();
    }

    /**
     * @brief Sets the tube voltage [kV].
     * @param voltage       Tube voltage [kV].
     * @param updateSpecter If true (default), rebuild the spectrum cache immediately.
     */
    void setTubeVoltage(double voltage, bool updateSpecter = true)
    {
        m_tube.setVoltage(voltage);
        if (updateSpecter)
            tubeChanged();
    }

    /**
     * @brief Sets the anode angle [rad].
     * @param ang           Anode angle [rad].
     * @param updateSpecter If true (default), rebuild the spectrum cache immediately.
     */
    void setTubeAnodeAngle(double ang, bool updateSpecter = true)
    {
        m_tube.setAnodeAngle(ang);
        if (updateSpecter)
            tubeChanged();
    }

    /**
     * @brief Sets the anode angle [deg].
     * @param ang           Anode angle [deg].
     * @param updateSpecter If true (default), rebuild the spectrum cache immediately.
     */
    void setTubeAnodeAngleDeg(double ang, bool updateSpecter = true)
    {
        m_tube.setAnodeAngleDeg(ang);
        if (updateSpecter)
            tubeChanged();
    }

    /**
     * @brief Adds a filtration material to the tube.
     * @param Z             Atomic number of the filter material.
     * @param mm            Filter thickness [mm].
     * @param updateSpecter If true (default), rebuild the spectrum cache immediately.
     */
    void addTubeFiltrationMaterial(std::size_t Z, double mm, bool updateSpecter = true)
    {
        m_tube.addFiltrationMaterial(Z, mm);
        if (updateSpecter)
            tubeChanged();
    }

    /**
     * @brief Returns the filtration thickness [mm] of element @p Z.
     * @param Z  Atomic number of the filter material.
     * @return Filter thickness [mm], or 0 if not present.
     */
    double tubeFiltration(std::size_t Z) const
    {
        return m_tube.filtration(Z);
    }

    /**
     * @brief Removes all filtration materials from the tube.
     * @param updateSpecter If true (default), rebuild the spectrum cache immediately.
     */
    void clearTubeFiltrationMaterials(bool updateSpecter = true)
    {
        m_tube.clearFiltrationMaterials();
        if (updateSpecter)
            tubeChanged();
    }

    /**
     * @brief Sets the energy resolution of the tube spectrum.
     * @param energyResolution  Energy bin width [keV].
     * @param updateSpecter     If true (default), rebuild the spectrum cache immediately.
     */
    void setTubeEnergyResolution(double energyResolution, bool updateSpecter = true)
    {
        m_tube.setEnergyResolution(energyResolution);
        if (updateSpecter)
            tubeChanged();
    }

    /// @brief Returns the aluminium half-value layer of the tube spectrum [mm Al].
    double tubeAlHalfValueLayer()
    {
        return m_tube.mmAlHalfValueLayer();
    }

    /// @brief Returns the mean photon energy of the tube spectrum [keV].
    double tubeMeanSpecterEnergy()
    {
        return m_tube.meanSpecterEnergy();
    }

    /// @brief Returns the symmetric collimation half-angles `{half_x, half_y}` [rad].
    const std::array<double, 2>& collimationHalfAngles() const
    {
        return m_collimationHalfAngles;
    }

    /**
     * @brief Sets the symmetric collimation half-angles from an array; each clamped to [0, π].
     * @param angles  `{half_x, half_y}` [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 2>& angles)
    {
        setCollimationHalfAngles(angles[0], angles[1]);
    }

    /**
     * @brief Sets the symmetric collimation half-angles from two scalars; each clamped to [0, π].
     * @param X  Half-angle in the x-direction [rad].
     * @param Y  Half-angle in the y-direction [rad].
     */
    void setCollimationHalfAngles(double X, double Y)
    {
        m_collimationHalfAngles[0] = std::clamp(X, 0.0, PI_VAL());
        m_collimationHalfAngles[1] = std::clamp(Y, 0.0, PI_VAL());
    }

    /// @brief Returns the symmetric collimation half-angles `{half_x, half_y}` [deg].
    std::array<double, 2> collimationHalfAnglesDeg() const
    {
        auto d = m_collimationHalfAngles;
        d[0] *= RAD_TO_DEG();
        d[1] *= RAD_TO_DEG();
        return d;
    }

    /**
     * @brief Sets the symmetric collimation half-angles from an array [deg].
     * @param angles  `{half_x, half_y}` [deg].
     */
    void setCollimationHalfAnglesDeg(const std::array<double, 2>& angles)
    {
        setCollimationHalfAnglesDeg(angles[0], angles[1]);
    }

    /**
     * @brief Sets the symmetric collimation half-angles from two scalars [deg].
     * @param X  Half-angle in the x-direction [deg].
     * @param Y  Half-angle in the y-direction [deg].
     */
    void setCollimationHalfAnglesDeg(double X, double Y)
    {
        setCollimationHalfAngles(DEG_TO_RAD() * X, DEG_TO_RAD() * Y);
    }

    /**
     * @brief Sets the collimation half-angles from a physical beam size at a known SDD.
     *
     * Converts beam dimensions at the detector plane to half-angles using:
     *   half_angle = atan(beamSize / (2 × sourceDetectorDistance))
     *
     * Has no effect if @p sourceDetectorDistance ≤ 0.
     *
     * @param beamSizeX              Full beam width at the detector [cm].
     * @param beamSizeY              Full beam height at the detector [cm].
     * @param sourceDetectorDistance Source-to-detector distance [cm].
     */
    void setBeamSize(double beamSizeX, double beamSizeY, double sourceDetectorDistance)
    {
        if (sourceDetectorDistance > 0) {
            setCollimationHalfAngles(
                m_collimationHalfAngles[0] = std::atan(std::abs(beamSizeX) / (2 * sourceDetectorDistance)),
                m_collimationHalfAngles[1] = std::atan(std::abs(beamSizeY) / (2 * sourceDetectorDistance)));
        }
    }

    /**
     * @brief Returns a `DXBeamExposure` for the given index.
     *
     * All exposures are identical — the index @p i is accepted for `BeamType`
     * interface compatibility but is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A fully configured `DXBeamExposure` with an owned copy of the spectrum.
     */
    DXBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {

        DXBeamExposure<ENABLETRACKING> exp(m_pos, m_dirCosines, m_particlesPerExposure, m_weight, m_collimationHalfAngles, m_specter);
        return exp;
    }

    /**
     * @brief Returns the dose calibration factor for this beam.
     *
     * Computes the ratio of the measured DAP to the expected air KERMA from the
     * simulated fluence:
     *
     *   factor = measuredDAP / (numberOfParticles × Σ w_i · E_i · μ_tr/ρ(E_i))
     *
     * where the sum is over tube spectrum bins and μ_tr/ρ is the NIST dry-air
     * mass energy-transfer coefficient. Returns 0 if the NIST air material cannot
     * be loaded.
     *
     * @param progress  Unused; accepted for `BeamType` interface compatibility.
     * @return Calibration factor to multiply accumulated energy scores by.
     */
    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, true);
        auto air_cand = Material<5>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();

        const auto kerma_per_history = std::transform_reduce(std::execution::par_unseq, energies.cbegin(), energies.cend(), weights.cbegin(), 0.0, std::plus<>(), [&](const auto e, const auto w) -> double {
            const auto uen = air.massEnergyTransferAttenuation(e);
            return w * e * uen;
        });

        const auto kerma_total = kerma_per_history * numberOfParticles();
        return m_measuredDAP / kerma_total;
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMDXBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMDXBeam";
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether @p data begins with the expected magic identifier.
     * @param data  Byte span to inspect; must be at least 32 bytes.
     * @return True if the first 32 bytes match `magicID()`, false otherwise.
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the beam configuration to a byte buffer.
     *
     * Writes position, both direction cosine vectors, collimation half-angles, exposure
     * count, particles per exposure, weight, measured DAP, and the serialized `Tube`
     * using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_dirCosines[0], buffer);
        Serializer::serialize(m_dirCosines[1], buffer);
        Serializer::serialize(m_collimationHalfAngles, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_measuredDAP, buffer);

        Serializer::serializeItem(m_tube, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a `DXBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`, deserializes the embedded `Tube`
     * sub-block, and calls `tubeChanged()` to rebuild the spectrum cache. Returns
     * `std::nullopt` if the `Tube` sub-block cannot be deserialized.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<DXBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<DXBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> pos;
        buffer = Serializer::deserialize(pos, buffer);

        std::array<std::array<double, 3>, 2> cosines;
        buffer = Serializer::deserialize(cosines[0], buffer);
        buffer = Serializer::deserialize(cosines[1], buffer);

        DXBeam<ENABLETRACKING> item(pos, cosines, { }, false);
        buffer = Serializer::deserialize(item.m_collimationHalfAngles, buffer);
        buffer = Serializer::deserialize(item.m_Nexposures, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);
        buffer = Serializer::deserialize(item.m_weight, buffer);
        buffer = Serializer::deserialize(item.m_measuredDAP, buffer);

        auto name = Serializer::getNameIDTemplate();
        auto tube_buffer = Serializer::getEmptyBuffer();
        buffer = Serializer::deserializeItem(name, tube_buffer, buffer);

        auto tube_opt = Tube::deserialize(tube_buffer);
        if (!tube_opt)
            return std::nullopt;

        item.m_tube = tube_opt.value();
        item.tubeChanged();
        return std::make_optional(item);
    }

protected:
    /**
     * @brief Rebuilds the `SpecterDistribution` cache from the current `Tube` state.
     *
     * Called automatically after any tube parameter change. Queries the tube for its
     * energy bins and unnormalised weights, then constructs a new `SpecterDistribution`
     * for use in `exposure()`.
     */
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };                                               ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };   ///< Normalised direction cosines {cos_x, cos_y} perpendicular to the beam axis.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };                                ///< Symmetric collimation half-angles {half_x, half_y} [rad].
    std::uint64_t m_Nexposures = 100;                                                          ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100;                                                ///< Photon histories per exposure.
    double m_weight = 1;                                                                       ///< Base photon weight.
    double m_measuredDAP = 1;                                                                  ///< Measured dose-area product for calibration [Gy·cm²].
    Tube m_tube;                                                                               ///< X-ray tube model (owned copy).
    SpecterDistribution<double> m_specter;                                                     ///< Spectrum cache rebuilt on tube changes.
};
}
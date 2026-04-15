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
#include <span>

namespace xraymc {

/**
 * @brief A single gantry exposure used internally by `CBCTBeam`.
 *
 * Represents the X-ray tube at one gantry angle of a CBCT arc scan. The source position
 * is SDD/2 from the isocentre along the beam direction derived from the supplied direction
 * cosines. Photon energies are drawn from a `SpecterDistribution` and directions are
 * sampled from a symmetric `SphereSamplingRectangularField` over the collimated field.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class CBCTBeamExposure {
public:
    /**
     * @brief Constructs a single CBCT exposure.
     *
     * Derives the central beam direction as `dircosines[0] × dircosines[1]` and builds a
     * `SphereSamplingRectangularField` for the supplied symmetric collimation half-angles.
     *
     * @param pos                  Source position in world space [cm].
     * @param dircosines           Two orthonormal vectors `{cos_x, cos_y}` spanning the plane
     *                             perpendicular to the beam axis.
     * @param N                    Number of photon histories for this exposure.
     * @param weight               Statistical weight of each sampled photon.
     * @param collimationHalfAngles `{half_x, half_y}` symmetric collimation half-angles [rad].
     * @param specter              Non-owning pointer to the energy spectrum sampler.
     */
    CBCTBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double>* specter)
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

    /// @brief Default constructor is deleted — an exposure must always have geometry and a spectrum.
    CBCTBeamExposure() = delete;

    /// @brief Returns the source position for this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the two direction cosines `{cos_x, cos_y}` for this exposure.
    const std::array<std::array<double, 3>, 2>& directionCosines() const { return m_dirCosines; }

    /// @brief Returns the symmetric collimation half-angles `{half_x, half_y}` [rad].
    const std::array<double, 2>& collimationHalfAngles() const
    {
        return m_collimationHalfAngles;
    }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

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
                .energy = m_specter->sampleValue(state),
                .weight = m_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specter->sampleValue(state),
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
    const SpecterDistribution<double>* m_specter = nullptr;                                    ///< Non-owning pointer to the energy spectrum.
};

/**
 * @brief A cone-beam CT (CBCT) arc beam, satisfying `BeamType`.
 *
 * Simulates an X-ray tube rotating around a patient isocentre over a configurable angular
 * arc. Exposures are evenly spaced by `m_angleStep` between `m_angleStart` and `m_angleStop`.
 * At each gantry angle the source is placed at SDD/2 from the isocentre; the in-plane normal
 * and the rotation axis form the direction cosine frame. Photon energies are drawn from a
 * `SpecterDistribution` built from the internal `Tube` model; directions are sampled over a
 * symmetric rectangular collimation field `{half_x, half_y}`.
 *
 * `calibrationFactor()` converts accumulated energy scores to dose by dividing a
 * user-supplied measured DAP value by the expected air KERMA from the simulated fluence,
 * using the dry-air mass energy-transfer coefficient.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CBCTBeam {
public:
    /**
     * @brief Constructs a CBCT beam centred at @p isocenter rotating around @p rotationAxis.
     *
     * @param isocenter           Isocentre position [cm]. Default: origin.
     * @param rotationAxis        Gantry rotation axis; normalised internally. Default: +z axis.
     * @param filtrationMaterials Map of `{atomic number Z → thickness [mm]}` added as filtration
     *                            to the internal `Tube`. Default: empty (no extra filtration).
     */
    CBCTBeam(
        const std::array<double, 3>& isocenter = { 0, 0, 0 },
        const std::array<double, 3>& rotationAxis = { 0, 0, 1 },
        const std::map<std::size_t, double>& filtrationMaterials = { })
        : m_isocenter(isocenter)
    {
        setRotationAxis(rotationAxis);
        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }

    /**
     * @brief Returns the number of gantry exposures in the arc.
     *
     * Computed as max(1, |stopAngle − startAngle| / stepAngle).
     */
    std::uint64_t numberOfExposures() const
    {
        const auto dAngle = std::abs(m_angleStop - m_angleStart);
        const auto steps = std::max(1.0, dAngle / m_angleStep);
        return static_cast<std::uint64_t>(steps);
    }

    /// @brief Returns the total number of photon histories across all exposures.
    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }

    /// @brief Returns the number of photon histories per exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    /**
     * @brief Sets the number of photon histories per exposure (clamped to ≥ 1).
     * @param n  Histories per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = std::max(n, std::uint64_t { 1 }); }

    /// @brief Returns the isocentre position [cm].
    const std::array<double, 3>& isocenter() const { return m_isocenter; }

    /**
     * @brief Sets the isocentre position [cm].
     * @param pos  3-D isocentre position in world space [cm].
     */
    void setIsocenter(const std::array<double, 3>& pos) { m_isocenter = pos; }

    /// @brief Returns the normalised gantry rotation axis.
    const std::array<double, 3>& rotationAxis() const
    {
        return m_direction;
    }

    /**
     * @brief Sets the gantry rotation axis; normalised internally.
     * @param dir  Rotation axis vector (need not be a unit vector).
     */
    void setRotationAxis(const std::array<double, 3>& dir)
    {
        m_direction = vectormath::normalized(dir);
    }

    /// @brief Returns the arc start angle [rad].
    double startAngle() const { return m_angleStart; }
    /// @brief Sets the arc start angle [rad].
    void setStartAngle(double angle) { m_angleStart = angle; }
    /// @brief Returns the arc start angle [deg].
    double startAngleDeg() const { return m_angleStart * RAD_TO_DEG(); }
    /// @brief Sets the arc start angle [deg].
    void setStartAngleDeg(double angle) { m_angleStart = angle * DEG_TO_RAD(); }

    /// @brief Returns the angular step between consecutive exposures [rad].
    double stepAngle() const { return m_angleStep; }

    /**
     * @brief Sets the angular step between consecutive exposures [rad].
     *
     * Clamped to a minimum of 0.1° to avoid degenerate arcs.
     *
     * @param angle  Desired step size [rad]; the absolute value is used.
     */
    void setStepAngle(double angle)
    {
        m_angleStep = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }
    /// @brief Returns the angular step between consecutive exposures [deg].
    double stepAngleDeg() const { return m_angleStep * RAD_TO_DEG(); }
    /// @brief Sets the angular step between consecutive exposures [deg].
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    /// @brief Returns the arc stop angle [rad].
    double stopAngle() const { return m_angleStop; }
    /// @brief Sets the arc stop angle [rad].
    void setStopAngle(double angle) { m_angleStop = angle; }
    /// @brief Returns the arc stop angle [deg].
    double stopAngleDeg() const { return m_angleStop * RAD_TO_DEG(); }
    /// @brief Sets the arc stop angle [deg].
    void setStopAngleDeg(double angle) { m_angleStop = angle * DEG_TO_RAD(); }

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
     * @brief Replaces the internal tube model; rebuilds the spectrum cache.
     * @param tube  New `Tube` instance (moved in).
     */
    void setTube(const Tube&& tube)
    {
        m_tube = tube;
        tubeChanged();
    }

    /**
     * @brief Sets the tube voltage [kV]; rebuilds the spectrum cache.
     * @param voltage  Tube voltage [kV].
     */
    void setTubeVoltage(double voltage)
    {
        m_tube.setVoltage(voltage);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle [rad]; rebuilds the spectrum cache.
     * @param ang  Anode angle [rad].
     */
    void setTubeAnodeAngle(double ang)
    {
        m_tube.setAnodeAngle(ang);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle [deg]; rebuilds the spectrum cache.
     * @param ang  Anode angle [deg].
     */
    void setTubeAnodeAngleDeg(double ang)
    {
        m_tube.setAnodeAngleDeg(ang);
        tubeChanged();
    }

    /**
     * @brief Adds a filtration material to the tube; rebuilds the spectrum cache.
     * @param Z   Atomic number of the filter material.
     * @param mm  Filter thickness [mm].
     */
    void addTubeFiltrationMaterial(std::size_t Z, double mm)
    {
        m_tube.addFiltrationMaterial(Z, mm);
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

    /// @brief Removes all filtration materials from the tube; rebuilds the spectrum cache.
    void clearTubeFiltrationMaterials()
    {
        m_tube.clearFiltrationMaterials();
        tubeChanged();
    }

    /**
     * @brief Sets the energy resolution of the tube spectrum; rebuilds the spectrum cache.
     * @param energyResolution  Energy bin width [keV].
     */
    void setTubeEnergyResolution(double energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
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
     * @brief Sets the symmetric collimation half-angles from an array; each clamped to [0, π/2].
     * @param angles  `{half_x, half_y}` [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 2>& angles)
    {
        setCollimationHalfAngles(angles[0], angles[1]);
    }

    /**
     * @brief Sets the symmetric collimation half-angles from two scalars; each clamped to [0, π/2].
     * @param X  Half-angle in the x-direction [rad].
     * @param Y  Half-angle in the y-direction [rad].
     */
    void setCollimationHalfAngles(double X, double Y)
    {
        m_collimationHalfAngles[0] = std::clamp(X, 0.0, PI_VAL() * 0.5);
        m_collimationHalfAngles[1] = std::clamp(Y, 0.0, PI_VAL() * 0.5);
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

    /// @brief Returns the source-to-detector distance [cm].
    double sourceDetectorDistance() const
    {
        return m_SDD;
    }

    /**
     * @brief Sets the source-to-detector distance [cm]; clamped to ≥ 1 cm.
     * @param SDD_cm  Source-to-detector distance [cm].
     */
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    /**
     * @brief Returns the `CBCTBeamExposure` for gantry position @p i.
     *
     * Computes the gantry angle as `startAngle + i × stepAngle` (or
     * `startAngle − i × stepAngle` when startAngle > stopAngle). Derives the in-plane
     * normal by rotating a seed vector around the rotation axis by that angle, then
     * places the source at isocentre − SDD/2 × (normal × rotationAxis).
     *
     * @param i  Zero-based exposure index; should be in [0, numberOfExposures()).
     * @return A fully configured `CBCTBeamExposure` for the requested gantry position.
     */
    CBCTBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        auto angle = i * m_angleStep;
        if (m_angleStart > m_angleStop)
            angle = -angle;

        // finding normal vector to direction
        const auto normal_ind = vectormath::argmin3(m_direction);
        std::array<double, 3> normal = { 0, 0, 0 };
        normal[normal_ind] = 1;
        normal = vectormath::normalized(vectormath::cross(normal, m_direction));
        normal = vectormath::rotate(normal, m_direction, m_angleStart + angle);

        const auto beamdir = vectormath::cross(normal, m_direction);

        const std::array<std::array<double, 3>, 2> cosines = { normal, m_direction };

        auto pos = vectormath::add(m_isocenter, vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis

        CBCTBeamExposure<ENABLETRACKING> exp(pos, cosines, m_particlesPerExposure, m_weight, m_collimationHalfAngles, &m_specter);
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
     * where the sum is over the tube spectrum bins and μ_tr/ρ is the NIST dry-air
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
     * @return Fixed-length tag "BEAMCBCTBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMCBCTBeam";
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
     * Writes isocentre, rotation axis, collimation half-angles, start/step/stop angles,
     * SDD, particles per exposure, weight, measured DAP, and the serialized `Tube` using
     * the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        Serializer::serialize(m_isocenter, buffer);
        Serializer::serialize(m_direction, buffer);
        Serializer::serialize(m_collimationHalfAngles, buffer);
        Serializer::serialize(m_angleStart, buffer);
        Serializer::serialize(m_angleStep, buffer);
        Serializer::serialize(m_angleStop, buffer);
        Serializer::serialize(m_SDD, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_measuredDAP, buffer);
        Serializer::serializeItem(m_tube, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a `CBCTBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`, deserializes the embedded `Tube`, and
     * calls `tubeChanged()` to rebuild the spectrum cache. Returns `std::nullopt` if
     * the `Tube` sub-block cannot be deserialized.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CBCTBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<CBCTBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {

        std::array<double, 3> pos, dir;
        buffer = Serializer::deserialize(pos, buffer);
        buffer = Serializer::deserialize(dir, buffer);

        CBCTBeam<ENABLETRACKING> item(pos, dir);
        buffer = Serializer::deserialize(item.m_collimationHalfAngles, buffer);
        buffer = Serializer::deserialize(item.m_angleStart, buffer);
        buffer = Serializer::deserialize(item.m_angleStep, buffer);
        buffer = Serializer::deserialize(item.m_angleStop, buffer);
        buffer = Serializer::deserialize(item.m_SDD, buffer);
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
     * Called automatically after any tube parameter change (voltage, filtration, etc.).
     * Queries the tube for its energy bins and unnormalised weights, then constructs a
     * new `SpecterDistribution` for use in `exposure()`.
     */
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_isocenter = { 0, 0, 0 };         ///< Isocentre position [cm].
    std::array<double, 3> m_direction = { 0, 0, 1 };          ///< Normalised gantry rotation axis.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 }; ///< Symmetric collimation half-angles {half_x, half_y} [rad].
    double m_angleStart = 0;                                    ///< Arc start angle [rad].
    double m_angleStop = PI_VAL();                             ///< Arc stop angle [rad]. Default: π (half rotation).
    double m_angleStep = PI_VAL() / 180;                      ///< Angular step between exposures [rad]. Default: 1°.
    double m_SDD = 50;                                         ///< Source-to-detector distance [cm].
    std::uint64_t m_particlesPerExposure = 100000;             ///< Photon histories per gantry position.
    double m_weight = 1;                                       ///< Base photon weight.
    double m_measuredDAP = 1;                                  ///< Measured dose-area product for calibration [Gy·cm²].
    Tube m_tube;                                               ///< X-ray tube model (owned copy).
    SpecterDistribution<double> m_specter;                     ///< Spectrum cache rebuilt on tube changes.
};
}
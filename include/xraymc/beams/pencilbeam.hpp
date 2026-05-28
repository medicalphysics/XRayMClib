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

Copyright 2022 Erlend Andersen
*/

#pragma once

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
 * @brief A single pencil-beam exposure: a monoenergetic, mono-directional particle source.
 *
 * Returned by `PencilBeam::exposure()` and consumed by the transport worker loop.
 * All particles from one exposure share the same position, direction, energy, and weight.
 * When `ENABLETRACKING` is true, `sampleParticle` returns a `ParticleTrack` with the
 * initial position already registered; otherwise it returns a plain `Particle`.
 *
 * @tparam ENABLETRACKING  If true, returns `ParticleTrack` from `sampleParticle()`
 *                         so that transport history can be recorded. Default: false.
 */
template <bool ENABLETRACKING = false>
class PencilBeamExposure {
public:
    /**
     * @brief Constructs an exposure with the given geometry and particle parameters.
     *
     * @param pos     Source position [cm].
     * @param dir     Beam direction unit vector (should be normalised by the caller).
     * @param energy  Photon energy [keV].
     * @param weight  Statistical weight of each sampled particle.
     * @param N       Number of particles to emit in this exposure.
     */
    PencilBeamExposure(const std::array<double, 3>& pos, const std::array<double, 3>& dir, double energy, double weight, std::uint64_t N)
        : m_energy(energy)
        , m_weight(weight)
        , m_pos(pos)
        , m_dir(dir)
        , m_NParticles(N)
    {
    }

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the number of particles to be emitted in this exposure.
    std::uint64_t numberOfParticles() const { return m_NParticles; }

    /**
     * @brief Samples (constructs) a particle for this exposure.
     *
     * All particles are identical — position, direction, energy, and weight are fixed.
     * When `ENABLETRACKING` is true, returns a `ParticleTrack` with the start position
     * already registered in the circular history buffer; otherwise returns a `Particle`.
     * The `state` parameter is accepted for interface compatibility but is not used.
     *
     * @param state  Per-thread PRNG state (unused by this exposure type).
     * @return `Particle` or `ParticleTrack` depending on `ENABLETRACKING`.
     */
    auto sampleParticle(RandomState& state) const noexcept
    {
        if constexpr (ENABLETRACKING) {
            ParticleTrack p = { .pos = m_pos,
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            p.registerPosition();
            return p;
        } else {
            Particle p = { .pos = m_pos,
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            return p;
        }
    }

private:
    double m_energy = 60; ///< Photon energy [keV].
    double m_weight = 1; ///< Statistical weight of each particle.
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Beam direction unit vector.
    std::uint64_t m_NParticles = 100; ///< Number of particles in this exposure.
};

/**
 * @brief A monoenergetic, mono-directional pencil beam satisfying `BeamType`.
 *
 * Emits `numberOfExposures() × numberOfParticlesPerExposure()` photons all from the
 * same position, in the same direction, at a single energy. Useful for beam-quality
 * studies, point-kernel calculations, and unit tests.
 *
 * The beam calibration is expressed as an air-KERMA value (`setAirKerma()`). The
 * `calibrationFactor()` converts accumulated energy scores to dose by dividing the
 * prescribed air KERMA by the expected KERMA from the simulated fluence, using the
 * NIST dry-air mass energy-transfer coefficient at the beam energy.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, `exposure().sampleParticle()` returns `ParticleTrack`
 *                         instead of `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class PencilBeam {
public:
    /**
     * @brief Constructs a pencil beam with the given geometry and energy.
     *
     * @param pos     Source position [cm]. Default: origin.
     * @param dir     Beam direction; normalised internally. Default: +z.
     * @param energy  Photon energy [keV]. Default: 60 keV.
     */
    PencilBeam(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 }, double energy = 60)
        : m_energy(energy)
        , m_pos(pos)
        , m_dir(dir)
    {
        vectormath::normalize(m_dir);
    }

    /**
     * @brief Sets the photon energy [keV].
     * @param energy  Photon energy; the absolute value is used.
     */
    void setEnergy(double energy)
    {
        m_energy = std::abs(energy);
    }

    /// @brief Returns the photon energy [keV].
    double energy() const { return m_energy; }

    /// @brief Returns the total number of exposures.
    std::uint64_t numberOfExposures() const { return m_Nexposures; }

    /**
     * @brief Sets the number of exposures (clamped to ≥ 1).
     * @param n  Desired exposure count.
     */
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    /**
     * @brief Sets the number of particles emitted per exposure.
     * @param n  Particles per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /// @brief Returns the number of particles emitted per exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    /// @brief Returns the total number of particles across all exposures.
    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }

    /**
     * @brief Sets the source position [cm].
     * @param pos  3-D source position in world space [cm].
     */
    void setPosition(const std::array<double, 3>& pos)
    {
        m_pos = pos;
    }

    /// @brief Returns the source position [cm].
    const std::array<double, 3>& position() const
    {
        return m_pos;
    }

    /**
     * @brief Sets the beam direction; normalised internally.
     * @param dir  Beam direction vector (need not be a unit vector).
     */
    void setDirection(const std::array<double, 3>& dir)
    {
        m_dir = dir;
        vectormath::normalize(m_dir);
    }

    /// @brief Returns the normalised beam direction unit vector.
    const std::array<double, 3>& direction() const
    {
        return m_dir;
    }

    /**
     * @brief Returns two orthogonal unit vectors spanning the plane perpendicular to the beam.
     *
     * Constructs an arbitrary but consistent pair of direction cosines by finding the
     * component of `m_dir` with the smallest magnitude, placing a unit vector there, and
     * computing two successive cross products.
     *
     * @return Array `{cos_x, cos_y}` where both vectors are perpendicular to `direction()`.
     */
    std::array<std::array<double, 3>, 2> directionCosines() const
    {
        std::array<double, 3> cand = { 0, 0, 0 };
        const auto minIdx = vectormath::argmin3(m_dir);
        cand[minIdx] = 1;

        std::array<std::array<double, 3>, 2> cos;
        cos[0] = vectormath::cross(m_dir, cand);
        vectormath::normalize(cos[0]);
        cos[1] = vectormath::cross(m_dir, cos[0]);
        return cos;
    }

    /**
     * @brief Sets the beam direction from a pair of direction cosines.
     *
     * The beam direction is computed as the cross product `dir[0] × dir[1]` and then
     * normalised. `dir[0]` and `dir[1]` should be orthogonal unit vectors spanning the
     * plane perpendicular to the desired beam axis.
     *
     * @param dir  `{cos_x, cos_y}` — two orthogonal vectors in the beam-perpendicular plane.
     */
    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dir = vectormath::cross(dir[0], dir[1]);
        vectormath::normalize(m_dir);
    }

    /**
     * @brief Sets the beam direction from two explicit direction-cosine vectors.
     *
     * Equivalent to `setDirectionCosines({xdir, ydir})`.
     *
     * @param xdir  First in-plane direction cosine vector.
     * @param ydir  Second in-plane direction cosine vector, perpendicular to @p xdir.
     */
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dir = vectormath::cross(xdir, ydir);
        vectormath::normalize(m_dir);
    }

    /**
     * @brief Returns the collimation angles for this beam.
     *
     * A pencil beam has no angular divergence, so all four angles are zero.
     *
     * @return `{0, 0, 0, 0}` — x_min, y_min, x_max, y_max collimation angles [rad].
     */
    const std::array<double, 4> collimationAngles() const
    {
        return std::array { 0.0, 0.0, 0.0, 0.0 };
    }

    /**
     * @brief Sets the statistical weight of each sampled particle.
     * @param weight  Particle weight. Default: 1.
     */
    void setParticleWeight(double weight = 1)
    {
        m_weight = weight;
    }

    /**
     * @brief Returns a `PencilBeamExposure` for the given index.
     *
     * All exposures are identical (same position, direction, energy, weight, and particle
     * count), so the index @p i is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A `PencilBeamExposure` configured with the current beam parameters.
     */
    PencilBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        PencilBeamExposure<ENABLETRACKING> exp(m_pos, m_dir, m_energy, m_weight, m_particlesPerExposure);
        return exp;
    }

    /**
     * @brief Sets the prescribed air-KERMA used for dose calibration [arbitrary units].
     *
     * The absolute value of @p k is stored. Used by `calibrationFactor()` to scale
     * dose scores to the desired output unit.
     *
     * @param k  Air KERMA value; the absolute value is used.
     */
    void setAirKerma(double k)
    {
        m_airKerma = std::abs(k);
    }

    /// @brief Returns the prescribed air-KERMA value used for dose calibration.
    double airKerma() const { return m_airKerma; }

    /**
     * @brief Returns the dose calibration factor for this beam.
     *
     * Computes the ratio of the prescribed air KERMA to the expected air KERMA from
     * the simulated fluence:
     *
     *   factor = airKerma / (numberOfParticles · μ_tr/ρ(E) · E)
     *
     * where μ_tr/ρ is the NIST mass energy-transfer coefficient of dry air at the
     * beam energy. Returns 0 if the NIST air material cannot be loaded.
     *
     * @param progress  Unused; accepted for `BeamType` interface compatibility.
     * @return Calibration factor to multiply accumulated energy scores by.
     */
    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        auto air_cand = Material<5>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();
        const double kerma = m_energy * numberOfParticles() * air.massEnergyTransferAttenuation(m_energy); // keV/g
        return m_airKerma / kerma; //(mGy/kg)(g/keV)
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMPencilBeam" padded with spaces, used by the
     *         serializer to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMPencilBeam";
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
     * Writes energy, weight, air KERMA, position, direction, exposure count, and
     * particles-per-exposure using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_energy, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_airKerma, buffer);
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_dir, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a `PencilBeam` from a serialized byte buffer.
     *
     * Reads the seven fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<PencilBeam>` containing the restored beam.
     */
    static std::optional<PencilBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        PencilBeam<ENABLETRACKING> item;
        buffer = Serializer::deserialize(item.m_energy, buffer);
        buffer = Serializer::deserialize(item.m_weight, buffer);
        buffer = Serializer::deserialize(item.m_airKerma, buffer);
        buffer = Serializer::deserialize(item.m_pos, buffer);
        buffer = Serializer::deserialize(item.m_dir, buffer);
        buffer = Serializer::deserialize(item.m_Nexposures, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        return std::make_optional(item);
    }

private:
    double m_energy = 60; ///< Photon energy [keV].
    double m_weight = 1; ///< Statistical weight of each particle.
    double m_airKerma = 1; ///< Prescribed air KERMA for dose calibration.
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Normalised beam direction unit vector.
    std::uint64_t m_Nexposures = 100; ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100; ///< Particles emitted per exposure.
};
}
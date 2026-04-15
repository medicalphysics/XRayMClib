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

#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <span>

namespace xraymc {

/**
 * @brief A single exposure of a monoenergetic isotropic point source with rectangular collimation.
 *
 * Returned by `IsotropicMonoEnergyBeam::exposure()`. Emits photons from a fixed position
 * with directions sampled uniformly over an asymmetric rectangular solid angle defined by
 * four half-angles `{minX, minY, maxX, maxY}`. The local-frame direction is transformed to
 * world space using the orthonormal basis `{m_dirCosines[0], m_dirCosines[1], m_dir}` derived
 * from the caller-supplied direction cosines. All photons have a fixed energy and weight 1.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicMonoEnergyBeamExposure {
public:
    /**
     * @brief Constructs an exposure at the given position with the given orientation and energy.
     *
     * Derives the central beam direction as the cross product of the two supplied direction
     * cosines: `m_dir = dircosines[0] × dircosines[1]`.
     *
     * @param pos         Source position [cm].
     * @param dircosines  Two orthonormal vectors `{cos_x, cos_y}` spanning the plane
     *                    perpendicular to the central beam direction.
     * @param energy      Fixed photon energy [keV].
     * @param N           Number of photon histories.
     */
    IsotropicMonoEnergyBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, double energy, std::uint64_t N)
        : m_energy(energy)
        , m_pos(pos)
        , m_dirCosines(dircosines)
        , m_NParticles(N)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the asymmetric rectangular collimation half-angles and rebuilds the direction sampler.
     * @param angles  `{minX, minY, maxX, maxY}` half-angles [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 4>& angles)
    {
        m_collimationHalfAngles = angles;
        m_directionSampler.setData(m_collimationHalfAngles);
    }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const { return m_NParticles; }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * Draws a direction from `SphereSamplingRectangularField` in the local frame and
     * transforms it to world space via `particleDirection()`. Assigns the fixed energy
     * `m_energy` and weight 1.
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
                .energy = m_energy,
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_energy,
                .weight = 1
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
    SphereSamplingRectangularField m_directionSampler;                                        ///< Samples directions within the rectangular collimation field.
    double m_energy = 60;                                                                      ///< Fixed photon energy [keV].
    std::array<double, 3> m_pos = { 0, 0, 0 };                                               ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 };                                               ///< Central beam direction (derived from direction cosines).
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };   ///< Orthonormal basis perpendicular to m_dir.
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 };                          ///< Rectangular collimation half-angles {minX, minY, maxX, maxY} [rad].
    std::uint64_t m_NParticles = 100;                                                          ///< Number of photon histories.
};

/**
 * @brief A monoenergetic isotropic point source with rectangular collimation, satisfying `BeamType`.
 *
 * Emits photons from a fixed position with directions sampled uniformly over an asymmetric
 * rectangular solid angle. The beam orientation is specified via two direction cosines
 * `{cos_x, cos_y}` (normalised internally); the central beam axis is their cross product.
 * Energy is set directly and clamped to [MIN_ENERGY(), MAX_ENERGY()]. All photons have
 * weight 1; dose calibration is left to the caller (`calibrationFactor()` returns 1).
 * All exposures are identical.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicMonoEnergyBeam {
public:
    /**
     * @brief Constructs a monoenergetic isotropic beam at the given position and orientation.
     *
     * @param pos         Source position [cm]. Default: origin.
     * @param dircosines  Two orthonormal vectors `{cos_x, cos_y}`; normalised internally.
     *                    Default: `{{1,0,0}, {0,1,0}}` (beam along +z axis).
     * @param energy      Photon energy [keV]. Default: 60 keV.
     */
    IsotropicMonoEnergyBeam(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<std::array<double, 3>, 2>& dircosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }, double energy = 60)
        : m_energy(energy)
        , m_pos(pos)
    {
        setDirectionCosines(dircosines);
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

    /**
     * @brief Sets the number of photon histories per exposure.
     * @param n  Histories per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /**
     * @brief Sets the source position [cm].
     * @param pos  3-D source position in world space [cm].
     */
    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }

    /// @brief Returns the source position [cm].
    const std::array<double, 3>& position() const { return m_pos; }

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

    /**
     * @brief Sets the photon energy [keV], clamped to [MIN_ENERGY(), MAX_ENERGY()].
     * @param energy  Desired photon energy [keV].
     */
    void setEnergy(double energy) { m_energy = std::min(std::max(MIN_ENERGY(), energy), MAX_ENERGY()); }

    /// @brief Returns the photon energy [keV].
    double energy() const { return m_energy; }

    /// @brief Returns the asymmetric rectangular collimation half-angles `{minX, minY, maxX, maxY}` [rad].
    const std::array<double, 4>& collimationHalfAngles() const { return m_collimationHalfAngles; }

    /**
     * @brief Sets the rectangular collimation half-angles from an array.
     * @param angles  `{minX, minY, maxX, maxY}` half-angles [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 4>& angles) { m_collimationHalfAngles = angles; }

    /**
     * @brief Sets the rectangular collimation half-angles from four individual values.
     * @param minX  Negative-x half-angle [rad].
     * @param minY  Negative-y half-angle [rad].
     * @param maxX  Positive-x half-angle [rad].
     * @param maxY  Positive-y half-angle [rad].
     */
    void setCollimationHalfAngles(double minX, double minY, double maxX, double maxY)
    {
        m_collimationHalfAngles[0] = minX;
        m_collimationHalfAngles[1] = minY;
        m_collimationHalfAngles[2] = maxX;
        m_collimationHalfAngles[3] = maxY;
    }

    /**
     * @brief Sets symmetric rectangular collimation half-angles from two magnitudes.
     *
     * The absolute values of @p X and @p Y are used; the field is centred on the beam axis:
     * `{−|X|, −|Y|, +|X|, +|Y|}`.
     *
     * @param X  Half-angle magnitude in the x-direction [rad].
     * @param Y  Half-angle magnitude in the y-direction [rad].
     */
    void setCollimationHalfAngles(double X, double Y)
    {
        X = std::abs(X);
        Y = std::abs(Y);
        m_collimationHalfAngles[0] = -X;
        m_collimationHalfAngles[1] = -Y;
        m_collimationHalfAngles[2] = X;
        m_collimationHalfAngles[3] = Y;
    }

    /**
     * @brief Returns an `IsotropicMonoEnergyBeamExposure` for the given index.
     *
     * All exposures are identical — the index @p i is accepted for `BeamType`
     * interface compatibility but is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A fully configured `IsotropicMonoEnergyBeamExposure`.
     */
    IsotropicMonoEnergyBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicMonoEnergyBeamExposure<ENABLETRACKING> exp(m_pos, m_dirCosines, m_energy, m_particlesPerExposure);
        exp.setCollimationHalfAngles(m_collimationHalfAngles);
        return exp;
    }

    /**
     * @brief Returns the dose calibration factor (always 1 for monoenergetic isotropic beams).
     *
     * Dose calibration is left to the caller. The `progress` parameter is accepted
     * for `BeamType` interface compatibility but is not used.
     *
     * @param progress  Unused.
     * @return 1.0.
     */
    double calibrationFactor(TransportProgress* progress = nullptr) const noexcept
    {
        return 1;
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMIsotropicBeamMono" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMIsotropicBeamMono";
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
     * Writes energy, position, both direction cosine vectors, collimation half-angles,
     * exposure count, and particles per exposure using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_energy, buffer);
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_dirCosines[0], buffer);
        Serializer::serialize(m_dirCosines[1], buffer);
        Serializer::serialize(m_collimationHalfAngles, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs an `IsotropicMonoEnergyBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<IsotropicMonoEnergyBeam>` containing the restored beam.
     */
    static std::optional<IsotropicMonoEnergyBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        IsotropicMonoEnergyBeam<ENABLETRACKING> item;
        buffer = Serializer::deserialize(item.m_energy, buffer);
        buffer = Serializer::deserialize(item.m_pos, buffer);
        buffer = Serializer::deserialize(item.m_dirCosines[0], buffer);
        buffer = Serializer::deserialize(item.m_dirCosines[1], buffer);
        buffer = Serializer::deserialize(item.m_collimationHalfAngles, buffer);
        buffer = Serializer::deserialize(item.m_Nexposures, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        return std::make_optional(item);
    }

private:
    double m_energy = 60;                                                                      ///< Fixed photon energy [keV].
    std::array<double, 3> m_pos = { 0, 0, 0 };                                               ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };   ///< Normalised direction cosines {cos_x, cos_y} perpendicular to the beam axis.
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 };                          ///< Rectangular collimation half-angles {minX, minY, maxX, maxY} [rad].
    std::uint64_t m_Nexposures = 100;                                                          ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100;                                                ///< Photon histories per exposure.
};

}
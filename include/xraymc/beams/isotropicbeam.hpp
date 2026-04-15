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
 * @brief A single exposure of an isotropic point source with a rectangular collimated field.
 *
 * Returned by `IsotropicBeam::exposure()`. Samples photon directions uniformly over the
 * solid angle defined by the asymmetric collimation half-angles using
 * `SphereSamplingRectangularField`, transforms them into world space via the stored
 * direction cosines, and draws energies from a `SpecterDistribution`. All photons have
 * weight 1.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicBeamExposure {
public:
    /**
     * @brief Constructs an exposure at the given position with the given beam orientation.
     *
     * Derives the central beam direction as the cross product of the two direction cosines.
     *
     * @param pos         Source position [cm].
     * @param dircosines  `{cos_x, cos_y}` — two orthogonal unit vectors defining the beam
     *                    plane. The central beam direction is `cos_x × cos_y`.
     * @param N           Number of photon histories. Default: 10⁶.
     */
    IsotropicBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_dirCosines(dircosines)
        , m_NParticles(N)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the asymmetric rectangular collimation half-angles [rad] and
     *        rebuilds the direction sampler.
     *
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 4>& angles)
    {
        m_collimationHalfAngles = angles;
        m_directionSampler.setData(m_collimationHalfAngles);
    }

    /**
     * @brief Sets the energy spectrum sampler for this exposure.
     * @param s  `SpecterDistribution` to copy into the exposure.
     */
    void setSpecterDistribution(const SpecterDistribution<double>& s)
    {
        m_specterDist = s;
    }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const { return m_NParticles; }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * Draws a direction from the rectangular-field sampler, transforms it to world
     * space, and draws an energy from the spectrum. Weight is always 1.
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
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            return p;
        }
    }

protected:
    /**
     * @brief Samples a direction in world space within the collimated field.
     *
     * Draws a direction from the rectangular-field sampler (in local beam coordinates)
     * and transforms it to world space using `changeBasis`.
     *
     * @param state  Per-thread PRNG state.
     * @return Unit direction vector in world space.
     */
    std::array<double, 3> particleDirection(RandomState& state) const
    {
        const auto dir = m_directionSampler(state);
        return vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dir);
    }

private:
    SphereSamplingRectangularField m_directionSampler;                              ///< Samples directions within the collimated field.
    std::array<double, 3> m_pos = { 0, 0, 0 };                                     ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 };                                     ///< Central beam direction (cos_x × cos_y).
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction cosines {cos_x, cos_y}.
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 };                ///< Asymmetric collimation bounds {x_min, y_min, x_max, y_max} [rad].
    std::uint64_t m_NParticles = 100;                                               ///< Number of photon histories.
    SpecterDistribution<double> m_specterDist;                                      ///< Energy spectrum sampler.
};

/**
 * @brief An isotropic point source with a rectangular collimated field, satisfying `BeamType`.
 *
 * Emits photons from a fixed point in all directions within an asymmetric rectangular
 * solid angle defined by four collimation half-angles `{x_min, y_min, x_max, y_max}`.
 * Photon energies are drawn from a user-supplied energy spectrum. All photons have
 * weight 1; dose calibration is left to the caller (`calibrationFactor()` returns 1).
 *
 * The beam orientation is set via two orthogonal direction cosines `{cos_x, cos_y}`;
 * the central beam direction is their cross product. All exposures are identical.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicBeam {
public:
    /**
     * @brief Constructs an isotropic beam at the given position and orientation.
     *
     * @param pos        Source position [cm]. Default: origin.
     * @param dircosines `{cos_x, cos_y}` — two orthogonal direction cosines defining the
     *                   beam plane. Normalised internally. Default: x- and y-axes.
     */
    IsotropicBeam(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<std::array<double, 3>, 2>& dircosines = { { { 1, 0, 0 }, { 0, 1, 0 } } })
        : m_pos(pos)
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

    /// @brief Returns the two direction cosines `{cos_x, cos_y}` defining the beam plane.
    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

    /**
     * @brief Sets the beam orientation from a pair of direction cosines.
     *
     * Both vectors are normalised internally. The central beam direction is
     * computed as `cos_x × cos_y`.
     *
     * @param dir  `{cos_x, cos_y}` — two orthogonal direction vectors.
     */
    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    /**
     * @brief Sets the beam orientation from two explicit direction vectors.
     *
     * Both vectors are normalised internally. The central beam direction is
     * computed as `xdir × ydir`.
     *
     * @param xdir  First in-plane direction cosine.
     * @param ydir  Second in-plane direction cosine; should be perpendicular to @p xdir.
     */
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dirCosines[0] = xdir;
        m_dirCosines[1] = ydir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    /**
     * @brief Sets the energy spectrum from a vector of `{energy [keV], weight}` pairs.
     *
     * Constructs a `SpecterDistribution` from the supplied histogram; weights are
     * normalised internally.
     *
     * @param specter  Vector of `{energy [keV], relative weight}` pairs.
     */
    void setEnergySpecter(const std::vector<std::pair<double, double>>& specter)
    {
        m_specter = SpecterDistribution(specter);
    }

    /// @brief Returns the current collimation half-angles `{x_min, y_min, x_max, y_max}` [rad].
    const std::array<double, 4>& collimationHalfAngles() const { return m_collimationHalfAngles; }

    /**
     * @brief Sets the asymmetric collimation half-angles from a four-element array [rad].
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds [rad].
     */
    void setCollimationHalfAngles(const std::array<double, 4>& angles) { m_collimationHalfAngles = angles; }

    /**
     * @brief Sets the asymmetric collimation half-angles from four scalar values [rad].
     * @param minX  Minimum x angular bound [rad].
     * @param minY  Minimum y angular bound [rad].
     * @param maxX  Maximum x angular bound [rad].
     * @param maxY  Maximum y angular bound [rad].
     */
    void setCollimationHalfAngles(double minX, double minY, double maxX, double maxY)
    {
        m_collimationHalfAngles[0] = minX;
        m_collimationHalfAngles[1] = minY;
        m_collimationHalfAngles[2] = maxX;
        m_collimationHalfAngles[3] = maxY;
    }

    /**
     * @brief Returns an `IsotropicBeamExposure` for the given index.
     *
     * All exposures are identical — position, orientation, collimation, and spectrum
     * are fixed. The index @p i is accepted for `BeamType` interface compatibility
     * but is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A fully configured `IsotropicBeamExposure`.
     */
    IsotropicBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicBeamExposure<ENABLETRACKING> exp(m_pos, m_dirCosines, m_particlesPerExposure);
        exp.setCollimationHalfAngles(m_collimationHalfAngles);
        exp.setSpecterDistribution(m_specter);
        return exp;
    }

    /**
     * @brief Returns the dose calibration factor (always 1 for isotropic beams).
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
     * @return Fixed-length tag "BEAMIsotropicBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMIsotropicBeam";
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
     * Writes position, direction cosines, collimation half-angles, exposure count,
     * particles per exposure, and the internal spectrum data using the `Serializer` format.
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
        Serializer::serialize(m_specter.copyInteralData(), buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs an `IsotropicBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`. Returns `std::nullopt` if the
     * spectrum internal data is invalid.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<IsotropicBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<IsotropicBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> pos;
        buffer = Serializer::deserialize(pos, buffer);

        std::array<std::array<double, 3>, 2> cosines;
        buffer = Serializer::deserialize(cosines[0], buffer);
        buffer = Serializer::deserialize(cosines[1], buffer);

        IsotropicBeam<ENABLETRACKING> item(pos, cosines);

        buffer = Serializer::deserialize(item.m_collimationHalfAngles, buffer);
        buffer = Serializer::deserialize(item.m_Nexposures, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        std::vector<double> specter_data;
        buffer = Serializer::deserialize(specter_data, buffer);
        auto specter_opt = item.m_specter.fromInternalData(specter_data);
        if (!specter_opt)
            return std::nullopt;
        else
            item.m_specter = specter_opt.value();

        return std::make_optional(item);
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };                                      ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction cosines {cos_x, cos_y}, normalised.
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 };                 ///< Collimation bounds {x_min, y_min, x_max, y_max} [rad].
    std::uint64_t m_Nexposures = 100;                                                ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100;                                      ///< Photon histories per exposure.
    SpecterDistribution<double> m_specter;                                           ///< Energy spectrum sampler (owned copy).
};
}
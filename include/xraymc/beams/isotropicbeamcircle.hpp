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

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "xraymc/beams/utilities/spheresamplingcircularfield.hpp"
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
 * @brief A single exposure of an isotropic point source collimated to a circular cone.
 *
 * Returned by `IsotropicBeamCircle::exposure()`. Samples photon directions uniformly
 * over the spherical cap defined by the cone half-angle around @p direction using
 * `SphereSamplingCircularField`, transforms them to world space, and draws energies
 * from a `SpecterDistribution`. All photons have weight 1.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicBeamCircleExposure {
public:
    /**
     * @brief Constructs an exposure at the given position pointing in the given direction.
     *
     * Computes an orthonormal basis (two direction cosines) perpendicular to @p direction
     * for use in the world-space transformation inside `sampleParticle`.
     *
     * @param pos        Source position [cm].
     * @param direction  Central beam direction unit vector (normalised by the caller).
     * @param N          Number of photon histories. Default: 10⁶.
     */
    IsotropicBeamCircleExposure(const std::array<double, 3>& pos, const std::array<double, 3>& direction, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_dir(direction)
        , m_NParticles(N)
    {

        m_dirCosines = calculateDirectionCosines(m_dir);
    }

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the cone half-angle and rebuilds the direction sampler.
     * @param angle  Cone half-angle [rad].
     */
    void setCollimationHalfAngle(double angle)
    {
        m_directionSampler.setData(angle);
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
     * Draws a direction from `SphereSamplingCircularField` (in the local +z frame),
     * transforms it to world space via `changeBasis`, then draws an energy from the
     * spectrum. Weight is always 1.
     *
     * @param state  Per-thread PRNG state.
     * @return `Particle` or `ParticleTrack` depending on `ENABLETRACKING`.
     */
    auto sampleParticle(RandomState& state) const noexcept
    {

        const auto dirz = m_directionSampler(state);
        const auto dir = vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dirz);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = dir,
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = dir,
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            return p;
        }
    }

protected:
    /**
     * @brief Computes an orthonormal basis perpendicular to @p dir.
     *
     * Finds the component of @p dir with the smallest magnitude, places a unit vector
     * there, then builds two orthogonal vectors via successive cross products:
     *   cos[0] = normalise(dir × k)
     *   cos[1] = cos[0] × dir
     *
     * @param dir  The central beam direction (assumed normalised).
     * @return `{cos_x, cos_y}` — two unit vectors spanning the plane perpendicular to @p dir.
     */
    static std::array<std::array<double, 3>, 2> calculateDirectionCosines(const std::array<double, 3>& dir)
    {
        std::array<std::array<double, 3>, 2> cos = { { { 0, 0, 0 }, { 0, 0, 0 } } };
        const auto minInd = vectormath::argmin3<std::uint_fast32_t, double>(dir);
        std::array<double, 3> k { 0, 0, 0 };
        k[minInd] = 1;

        const auto vec_xy_raw = vectormath::cross(dir, k);
        cos[0] = vectormath::normalized(vec_xy_raw);
        cos[1] = vectormath::cross(cos[0], dir);
        return cos;
    }

private:
    SphereSamplingCircularField m_directionSampler;                                      ///< Samples directions within the circular cone.
    std::array<double, 3> m_pos = { 0, 0, 0 };                                          ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 };                                          ///< Central beam direction unit vector.
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< Orthonormal basis perpendicular to m_dir.
    std::uint64_t m_NParticles = 100;                                                    ///< Number of photon histories.
    SpecterDistribution<double> m_specterDist;                                           ///< Energy spectrum sampler.
};

/**
 * @brief An isotropic point source collimated to a circular cone, satisfying `BeamType`.
 *
 * Emits photons from a fixed point within the spherical cap defined by a single cone
 * half-angle around a central beam direction. Photon energies are drawn from a
 * user-supplied energy spectrum. All photons have weight 1; dose calibration is left
 * to the caller (`calibrationFactor()` returns 1). All exposures are identical.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicBeamCircle {
public:
    /**
     * @brief Constructs an isotropic circular-cone beam at the given position and direction.
     *
     * @param pos  Source position [cm]. Default: origin.
     * @param dir  Central beam direction; normalised internally. Default: +x axis.
     */
    IsotropicBeamCircle(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 1, 0, 0 })
        : m_pos(pos)
        , m_dir(dir)
    {
        vectormath::normalize(m_dir);
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

    /// @brief Returns the central beam direction unit vector.
    const std::array<double, 3>& direction() const
    {
        return m_dir;
    }

    /**
     * @brief Sets the central beam direction; normalised internally.
     * @param dir  Beam direction vector (need not be a unit vector).
     */
    void setDirection(const std::array<double, 3>& dir)
    {
        m_dir = vectormath::normalized(dir);
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

    /**
     * @brief Sets the cone half-angle for the circular collimation field.
     *
     * Clamped to [0, π]. A value of 0 produces a pencil beam; π covers the full
     * hemisphere.
     *
     * @param ang  Cone half-angle [rad]; the absolute value is clamped to [0, π].
     */
    void setCollimationHalfAngle(double ang)
    {
        m_collimationHalfAngle = std::clamp(std::abs(ang), 0.0, PI_VAL());
    }

    /// @brief Returns the cone half-angle [rad].
    double collimationHalfAngle() const { return m_collimationHalfAngle; }

    /**
     * @brief Returns an `IsotropicBeamCircleExposure` for the given index.
     *
     * All exposures are identical — the index @p i is accepted for `BeamType`
     * interface compatibility but is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A fully configured `IsotropicBeamCircleExposure`.
     */
    IsotropicBeamCircleExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicBeamCircleExposure<ENABLETRACKING> exp(m_pos, m_dir, m_particlesPerExposure);
        exp.setCollimationHalfAngle(m_collimationHalfAngle);
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
     * @return Fixed-length tag "BEAMIsotropicBeamCircle" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMIsotropicBeamCircle";
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
     * Writes position, direction, cone half-angle, exposure count, particles per
     * exposure, and spectrum internal data using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_dir, buffer);
        Serializer::serialize(m_collimationHalfAngle, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serialize(m_specter.copyInteralData(), buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs an `IsotropicBeamCircle` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`. Returns `std::nullopt` if the
     * spectrum internal data is invalid.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<IsotropicBeamCircle>` on success, or `std::nullopt` on failure.
     */
    static std::optional<IsotropicBeamCircle<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> pos;
        buffer = Serializer::deserialize(pos, buffer);

        std::array<double, 3> dir;
        buffer = Serializer::deserialize(dir, buffer);

        IsotropicBeamCircle<ENABLETRACKING> item(pos, dir);
        buffer = Serializer::deserialize(item.m_collimationHalfAngle, buffer);
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
    std::array<double, 3> m_pos = { 0, 0, 0 };  ///< Source position [cm].
    std::array<double, 3> m_dir = { 1, 0, 0 };  ///< Normalised central beam direction.
    double m_collimationHalfAngle = 0;           ///< Cone half-angle [rad], clamped to [0, π].
    std::uint64_t m_Nexposures = 100;            ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100;  ///< Photon histories per exposure.
    SpecterDistribution<double> m_specter;        ///< Energy spectrum sampler (owned copy).
};
}
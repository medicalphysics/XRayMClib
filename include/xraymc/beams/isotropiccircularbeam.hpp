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
 * @brief A single exposure of an isotropic source distributed uniformly on a circular ring.
 *
 * Returned by `IsotropicCircularBeam::exposure()`. Each call to `sampleParticle()` places
 * the photon source at a uniformly random position on a circle of radius @p radius centred
 * at @p pos. The circle lies in the XY plane (z = 0 in world space): a random azimuthal
 * angle φ is drawn, the local X axis is rotated by φ around the Z axis, and the source is
 * placed at `pos + Z × (−radius)` where Z = X × Y (pointing inward toward the centre).
 * Photon directions are sampled from a `SphereSamplingRectangularField` and transformed to
 * world space via `changeBasis`. Energies are drawn from a `SpecterDistribution`. All
 * photons have weight 1.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicCircularBeamExposure {
public:
    /**
     * @brief Constructs an exposure centred at @p pos with the given ring radius.
     *
     * @param pos     Centre of the source ring in world space [cm].
     * @param radius  Radius of the source ring [cm]. Default: 0 (point source at centre).
     * @param N       Number of photon histories. Default: 10⁶.
     */
    IsotropicCircularBeamExposure(const std::array<double, 3>& pos, double radius = 0, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_radius(radius)
        , m_NParticles(N)
    {
    }

    /// @brief Returns the centre position of the source ring [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the asymmetric rectangular collimation half-angles and rebuilds the direction sampler.
     * @param angles  `{minX, minY, maxX, maxY}` half-angles [rad] passed to `SphereSamplingRectangularField`.
     */
    void setCollimationHalfAngles(const std::array<double, 4>& angles)
    {
        m_directionSampler.setData(angles);
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
     * 1. Draws a uniform random azimuthal angle φ ∈ [0, 2π).
     * 2. Rotates the default Y-axis direction `{0,1,0}` by φ around `{0,0,1}` to obtain X.
     * 3. Sets Y = `{0,0,1}` and computes Z = X × Y (the inward radial direction).
     * 4. Places the source at `m_pos + Z × (−m_radius)`, i.e. on the ring at angle φ.
     * 5. Draws a direction from `SphereSamplingRectangularField` in the local frame and
     *    transforms it to world space via `changeBasis(X, Y, Z, dirz)`.
     * 6. Draws an energy from the spectrum. Weight is always 1.
     *
     * @param state  Per-thread PRNG state.
     * @return `Particle` or `ParticleTrack` depending on `ENABLETRACKING`.
     */
    auto sampleParticle(RandomState& state) const noexcept
    {
        const auto angle = state.randomUniform(2 * PI_VAL());
        const auto sang = std::sin(angle);
        const auto cang = std::cos(angle);
        const auto X = vectormath::rotate({ 0, 1, 0 }, { 0, 0, 1 }, sang, cang);
        constexpr std::array<double, 3> Y = { 0, 0, 1 };
        const auto Z = vectormath::cross(X, Y);
        const auto pos = vectormath::add(vectormath::scale(Z, -m_radius), m_pos);

        const auto dirz = m_directionSampler(state);
        const auto dir = vectormath::changeBasis(X, Y, Z, dirz);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = pos,
                .dir = dir,
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = pos,
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
     * @param dir  A direction vector (assumed normalised).
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
    SphereSamplingRectangularField m_directionSampler;        ///< Samples directions within the rectangular collimation field.
    std::array<double, 3> m_pos = { 0, 0, 0 };               ///< Centre of the source ring [cm].
    std::uint64_t m_NParticles = 100;                         ///< Number of photon histories.
    double m_radius = 0;                                       ///< Radius of the source ring [cm].
    SpecterDistribution<double> m_specterDist;                 ///< Energy spectrum sampler.
};

/**
 * @brief An isotropic source distributed on a circular ring, satisfying `BeamType`.
 *
 * Emits photons from positions uniformly distributed on a circle of configurable radius
 * centred at @p pos. Each call to `exposure(i).sampleParticle()` selects a fresh random
 * azimuthal position on the ring, making this beam suitable for simulating annular or
 * rotationally symmetric source geometries (e.g. a CT ring source). Photon directions
 * are sampled from an asymmetric rectangular collimation field defined by four half-angles
 * `{minX, minY, maxX, maxY}`. Energies are drawn from a user-supplied spectrum. All
 * photons have weight 1; dose calibration is left to the caller (`calibrationFactor()`
 * returns 1). All exposures are identical.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class IsotropicCircularBeam {
public:
    /**
     * @brief Constructs an isotropic circular beam centred at @p pos.
     *
     * @param pos  Centre of the source ring in world space [cm]. Default: origin.
     */
    IsotropicCircularBeam(const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_pos(pos)
    {
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
     * @brief Sets the centre position of the source ring [cm].
     * @param pos  3-D centre position in world space [cm].
     */
    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }

    /// @brief Returns the centre position of the source ring [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /**
     * @brief Sets the radius of the source ring; the absolute value is used.
     * @param r  Ring radius [cm].
     */
    void setRadius(double r)
    {
        m_radius = std::abs(r);
    }

    /// @brief Returns the ring radius [cm].
    double radius() const { return m_radius; }

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
     * @brief Returns an `IsotropicCircularBeamExposure` for the given index.
     *
     * All exposures are identical — the index @p i is accepted for `BeamType`
     * interface compatibility but is ignored.
     *
     * @param i  Exposure index (ignored).
     * @return A fully configured `IsotropicCircularBeamExposure`.
     */
    IsotropicCircularBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicCircularBeamExposure<ENABLETRACKING> exp(m_pos, m_radius, m_particlesPerExposure);
        exp.setCollimationHalfAngles(m_collimationHalfAngles);
        exp.setSpecterDistribution(m_specter);
        return exp;
    }

    /**
     * @brief Returns the dose calibration factor (always 1 for isotropic circular beams).
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
     * @return Fixed-length tag "BEAMIsotropicCircularBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMIsotropicCircularBeam";
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
     * Writes position, ring radius, collimation half-angles, exposure count, particles
     * per exposure, and spectrum internal data using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_radius, buffer);
        Serializer::serialize(m_collimationHalfAngles, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serialize(m_specter.copyInteralData(), buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs an `IsotropicCircularBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`. Returns `std::nullopt` if the
     * spectrum internal data is invalid.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<IsotropicCircularBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<IsotropicCircularBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> pos;
        buffer = Serializer::deserialize(pos, buffer);

        IsotropicCircularBeam<ENABLETRACKING> item(pos);
        buffer = Serializer::deserialize(item.m_radius, buffer);
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
    std::array<double, 3> m_pos = { 0, 0, 0 };                   ///< Centre of the source ring [cm].
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 }; ///< Rectangular collimation half-angles {minX, minY, maxX, maxY} [rad].
    double m_radius = 0;                                           ///< Ring radius [cm].
    std::uint64_t m_Nexposures = 100;                             ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100;                   ///< Photon histories per exposure.
    SpecterDistribution<double> m_specter;                         ///< Energy spectrum sampler (owned copy).
};
}
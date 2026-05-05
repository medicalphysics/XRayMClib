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

Copyright 2026 Erlend Andersen
*/

#pragma once

#include "xraymc/interactions.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/particletracker.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

namespace xraymc {

/**
 * @brief A thin, planar personal dosimeter modelling an Hp(10) dose-equivalent detector.
 *
 * The dosimeter is represented as a flat rectangular slab (width × height × depth) with
 * a configurable orientation. When a photon crosses the sensitive plane it is fully
 * absorbed, and its contribution is scored both as air-KERMA and as Hp(10) dose
 * equivalent using tabulated conversion coefficients (Šolc et al., 2025).
 *
 * Optionally (controlled by `USEDIRECTIONALDEPENDENCE`) an empirical angular-response
 * model modelled on the RaySafe i3 personal dosimeter is applied to correct the scored
 * quantity for off-axis incidence.
 *
 * @tparam USEDIRECTIONALDEPENDENCE  When true (default), the scored energy is weighted
 *                                   by the angular-response model. When false, every
 *                                   photon contributes with its full weight regardless
 *                                   of incidence angle.
 */
template <bool USEDIRECTIONALDEPENDENCE = true>
class PersonalDosimeter {
public:
    /**
     * @brief Default constructor. Loads dry-air material data and sets up the
     *        Hp(10) conversion-factor spline and the initial AABB.
     */
    PersonalDosimeter()
        : m_air(Material<12>::byNistName("Air, Dry (near sea level)").value())
    {
        setupHp10ConvertionFactor();
        m_aabb = calculateAABB();
    }

    /// @brief Returns the axis-aligned bounding box as {xmin, ymin, zmin, xmax, ymax, zmax} [cm].
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Translates the dosimeter by @p dist, updating the AABB in place.
     * @param dist  Translation vector [cm].
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    /// @brief Returns the centre of the dosimeter in world coordinates [cm].
    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    /**
     * @brief Sets the centre of the dosimeter and recomputes the AABB.
     * @param c  New centre position in world coordinates [cm].
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        m_aabb = calculateAABB();
    }

    /**
     * @brief Sets the orientation via two in-plane direction cosines and recomputes the AABB.
     *
     * The normal vector (beam-facing direction) is the cross product of the two cosines.
     *
     * @param dir_cosX  Direction cosine along the dosimeter x-axis (starboard→port); will be normalised.
     * @param dir_cosY  Direction cosine along the dosimeter y-axis (bottom→top); will be normalised.
     */
    void setDirectionCosines(const std::array<double, 3>& dir_cosX, const std::array<double, 3>& dir_cosY)
    {
        m_direction_cosines[0] = vectormath::normalized(dir_cosX);
        m_direction_cosines[1] = vectormath::normalized(dir_cosY);
        m_aabb = calculateAABB();
    }

    /// @brief Returns the two normalised in-plane direction cosines {x-axis, y-axis}.
    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_direction_cosines;
    }

    /// @brief Returns the outward normal unit vector (cross product of the two direction cosines).
    std::array<double, 3> normalVector() const
    {
        return vectormath::cross(m_direction_cosines[0], m_direction_cosines[1]);
    }

    /**
     * @brief Rotates the dosimeter about an arbitrary axis and recomputes the AABB.
     * @param angle  Rotation angle [radians].
     * @param axis   Unit vector defining the rotation axis.
     */
    void rotate(double angle, const std::array<double, 3>& axis)
    {
        m_direction_cosines[0] = vectormath::rotate(m_direction_cosines[0], axis, angle);
        m_direction_cosines[1] = vectormath::rotate(m_direction_cosines[1], axis, angle);
        m_aabb = calculateAABB();
    }

    /**
     * @brief Tests a particle ray against the dosimeter geometry.
     *
     * First performs a fast AABB test, then refines to an exact plane intersection.
     * Particles that miss the rectangular sensitive area are rejected.
     *
     * @param p  Particle whose position and direction define the ray.
     * @return   Intersection result; `valid()` is false if the particle misses.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        auto res = basicshape::AABB::intersect(p, m_aabb);
        if (res.valid()) {
            const auto t = intersectPlane(p);
            if (t) {
                res.intersection = t.value();
                res.rayOriginIsInsideItem = false;
            } else {
                res.intersectionValid = false;
            }
        }
        return res;
    }

    /**
     * @brief Visualization variant of the intersection test.
     *
     * Behaves like `intersect()` but additionally stores the surface normal and
     * the current accumulated air-KERMA dose in the result for rendering purposes.
     *
     * @tparam U  Scalar type used for the visualization value.
     * @param p   Particle whose position and direction define the ray.
     * @return    Visualization intersection result with normal and dose value set on hit.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        auto res = basicshape::AABB::intersectVisualization<U>(p, m_aabb);
        if (res.valid()) {
            const auto t = intersectPlane(p);
            if (t) {
                res.intersection = t.value();
                res.rayOriginIsInsideItem = false;
                res.normal = normalVector();
                res.value = m_dose_airKermaScored.dose();
            } else {
                res.intersectionValid = false;
            }
        }
        return res;
    }

    /**
     * @brief Transports a particle through the dosimeter, scoring its contribution and absorbing it.
     *
     * Computes the air-KERMA and Hp(10) contributions from the photon energy. When
     * `USEDIRECTIONALDEPENDENCE` is true the scored values are multiplied by the
     * angular-response weight and a cosine area-correction factor. The particle is fully
     * absorbed (energy set to zero) after scoring.
     *
     * If @p P is `ParticleTrack`, the particle is also registered with the internal tracker.
     *
     * @param p      Particle to transport; energy is set to zero on return.
     * @param state  Random-number generator state (unused here but required by the interface).
     */
    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }

        const auto airKerma = p.energy * m_air.massEnergyTransferAttenuation(p.energy);
        const auto hp10 = airKerma * hp10ConvertionFactor(p.energy);

        if constexpr (USEDIRECTIONALDEPENDENCE) {
            const auto angular_weight = angularResponseWeight(p);

            // If we correct for the detector model we must cancel the area effect
            const auto areaScaling = -1.0 / vectormath::dot(normalVector(), p.dir);

            m_energyScored.scoreEnergy(hp10 * p.weight * angular_weight * areaScaling);
            m_energy_airKermaScored.scoreEnergy(airKerma * p.weight * angular_weight * areaScaling);
        } else {
            m_energyScored.scoreEnergy(hp10 * p.weight);
            m_energy_airKermaScored.scoreEnergy(airKerma * p.weight);
        }

        p.energy = 0; // we kill the particle, it is absorbed by the dosimeter
    }

    /**
     * @brief Returns the Hp(10) conversion factor h_pk(10, E, 0) for a given photon energy.
     *
     * The factor converts air-KERMA to Hp(10) dose equivalent [Sv/Gy]. Data are taken
     * from Šolc et al. (2025) and interpolated with a log-linear cubic spline.
     * Returns 0 for energies outside the tabulated range [7.5, MAX_ENERGY()] keV.
     *
     * @param energy  Photon energy [keV].
     * @return        Hp(10) conversion factor [Sv/Gy], or 0 if out of range.
     */
    double hp10ConvertionFactor(double energy) const
    {
        if (energy < 7.5 || energy > MAX_ENERGY()) {
            return 0;
        } else {
            // Unitless convertionfactor from air KERMA to Hp10 dose
            // Data is log-lin interpolated
            const auto energy_log = std::log(energy);
            return m_hp10factor_interpolator(energy_log);
        }
    }

    /**
     * @brief Returns the energy-score accumulator at the given index.
     *
     * Index 0 returns the Hp(10)-weighted energy scorer; index 1 returns the
     * air-KERMA energy scorer.
     *
     * @param index  Scorer index (0 = Hp(10), 1 = air-KERMA). Defaults to 0.
     * @return       Reference to the requested `EnergyScore`.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_energyScored;

        return m_energy_airKermaScored;
    }

    /// @brief Resets both energy-score accumulators to zero.
    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_energy_airKermaScored.clear();
    }

    /**
     * @brief Accumulates the current energy scores into the dose-score accumulators.
     *
     * Uses an effective area and density of 1.0 (the dosimeter's dose has already been
     * computed in the energy-score step). Optionally applies a beam calibration factor.
     *
     * @param calibration_factor  Multiplicative calibration factor (default 1).
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        constexpr double area = 1.0; // width() * height();
        constexpr double density = 1.0;

        // We have already calculated the dose in energy score.
        m_dose.addScoredEnergy(m_energyScored, area, density, calibration_factor);
        m_dose_airKermaScored.addScoredEnergy(m_energy_airKermaScored, area, density, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator at the given index.
     *
     * Index 0 returns the Hp(10) dose scorer; index 1 returns the air-KERMA dose scorer.
     *
     * @param index  Scorer index (0 = Hp(10), 1 = air-KERMA). Defaults to 0.
     * @return       Reference to the requested `DoseScore`.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_dose;

        return m_dose_airKermaScored;
    }

    /// @brief Resets both dose-score accumulators to zero.
    void clearDoseScored()
    {
        m_dose.clear();
        m_dose_airKermaScored.clear();
    }

    /// @brief Returns the accumulated air-KERMA dose score (convenience alias for `doseScored(1)`).
    const DoseScore& airKermaScored() const
    {
        return doseScored(1);
    }

    /// @brief Returns the accumulated Hp(10) dose score (convenience alias for `doseScored(0)`).
    const DoseScore& hp10DoseScored() const
    {
        return doseScored(0);
    }

    /**
     * @brief Computes the angular-response weight for a given particle direction.
     *
     * Models the directional sensitivity of the RaySafe i3 personal dosimeter using
     * an empirical cosine-power law (exponent 0.2) fitted separately in the x and y
     * in-plane directions. Angles are clamped to ±87° in x and an offset range in y
     * to account for the dosimeter being worn tilted 15° downward.
     *
     * Returns 0 for particles traveling away from the front face (dot product with
     * normal ≥ 0), preventing double-counting from behind.
     *
     * @param p  Particle whose direction is evaluated.
     * @return   Dimensionless angular-response weight in [0, 1].
     */
    double angularResponseWeight(const ParticleType auto& p) const
    {

        const auto normal = normalVector();
        const auto dz = vectormath::dot(p.dir, normal);

        if (dz >= 0) {
            return 0; // No response to particles traveling away from the front face
        }

        const auto dx = vectormath::dot(p.dir, m_direction_cosines[0]);
        const auto dy = vectormath::dot(p.dir, m_direction_cosines[1]);

        constexpr auto max_angle = 87.0 * std::numbers::pi_v<double> / 180.0;
        constexpr auto Y_corr = -15.0 * std::numbers::pi_v<double> / 180.0;

        const auto angX = std::clamp(std::atan2(dx, -dz), -max_angle, max_angle);
        const auto angY = std::clamp(std::atan2(dy, -dz), -max_angle + Y_corr, max_angle + Y_corr);

        // The model is a simple cosine to the power of 0.2
        // determined by cognitive fit (aka eyeballing) to data
        // from the RaySafe i3 personal dosimeter.
        const auto weightX = std::pow(std::cos(angX), 0.2);

        // 15 degrees pointed downwards for y weight at unity
        const auto weightY = std::pow(std::cos(angY - Y_corr), 0.2);

        return weightX * weightY;
    }

    /// @brief Returns the physical height of the dosimeter [cm] (58 mm).
    constexpr static double height()
    {
        return 58.0 / 10.0; // cm
    }

    /// @brief Returns the physical width of the dosimeter [cm] (40 mm).
    constexpr static double width()
    {
        return 40.0 / 10.0; // cm
    }

    /// @brief Returns the effective depth (thickness) of the sensitive layer [cm].
    constexpr static double depth()
    {
        return 1.0E-3;
        // return 17.0 / 10.0; // cm
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "PersonalDosimeter1" padded with spaces, used by the
     *         serializer to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "PersonalDosimeter1" + std::to_string(USEDIRECTIONALDEPENDENCE);
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether a raw data buffer begins with the expected magic identifier.
     * @param data Buffer to inspect; must be at least 32 bytes for a positive result.
     * @return true if the first 32 bytes match magicID().
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the dosimeter state (center, direction cosines, dose, and air kerma scores)
     *        to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_center, buffer);
        Serializer::serialize(m_direction_cosines[0], buffer);
        Serializer::serialize(m_direction_cosines[1], buffer);
        Serializer::serializeDoseScore(m_dose, buffer);
        Serializer::serializeDoseScore(m_dose_airKermaScored, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a scorer from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed scorer (always valid for well-formed input).
     */
    static std::optional<PersonalDosimeter<USEDIRECTIONALDEPENDENCE>> deserialize(std::span<const char> buffer)
    {
        PersonalDosimeter<USEDIRECTIONALDEPENDENCE> item;
        buffer = Serializer::deserialize(item.m_center, buffer);
        buffer = Serializer::deserialize(item.m_direction_cosines[0], buffer);
        buffer = Serializer::deserialize(item.m_direction_cosines[1], buffer);
        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);
        buffer = Serializer::deserializeDoseScore(item.m_dose_airKermaScored, buffer);
        item.m_aabb = item.calculateAABB();
        return std::make_optional(item);
    }

protected:
    /**
     * @brief Initialises the Hp(10) conversion-factor cubic spline interpolator.
     *
     * Loads tabulated h_pk(10, E, 0) data [Sv/Gy] from Šolc et al. (2025) and fits a
     * cubic spline on log-energy vs. conversion-factor coordinates for efficient lookup.
     */
    void setupHp10ConvertionFactor()
    {
        // Data from Šolc, J. et. al, 2025. Re-calculation of air kerma to dose-equivalent conversion coefficients
        // for mono-energetic photons. Journal of Radiological Protection 45
        // https://doi.org/10.1088/1361-6498/adda56
        // h_pk(10, E, 0) [Sv/Gy]
        constexpr std::array<std::pair<double, double>, 71> data = { { { 7, 0.0000009 },
            { 7.5, 0.0000094 },
            { 8, 0.0000737 },
            { 8.5, 0.0003775 },
            { 9, 0.0013479 },
            { 9.5, 0.0036878 },
            { 10, 0.0083370 },
            { 10.5, 0.0162043 },
            { 11, 0.0279772 },
            { 11.5, 0.0442240 },
            { 12, 0.0648614 },
            { 12.5, 0.0895749 },
            { 13, 0.1179816 },
            { 14, 0.1831417 },
            { 15, 0.2547105 },
            { 16, 0.3292878 },
            { 17, 0.4037368 },
            { 18, 0.4747351 },
            { 19, 0.5428724 },
            { 20, 0.6087834 },
            { 22, 0.7305480 },
            { 24, 0.8388087 },
            { 26, 0.9393330 },
            { 28, 1.0340408 },
            { 30, 1.1195263 },
            { 32, 1.2060748 },
            { 34, 1.2854513 },
            { 36, 1.3594318 },
            { 38, 1.4330483 },
            { 40, 1.4997648 },
            { 42, 1.5614685 },
            { 44, 1.6219570 },
            { 46, 1.6764643 },
            { 48, 1.7223829 },
            { 50, 1.7660656 },
            { 52, 1.8008878 },
            { 54, 1.8359704 },
            { 56, 1.8588267 },
            { 58, 1.8823515 },
            { 60, 1.9001095 },
            { 65, 1.9231631 },
            { 70, 1.9260116 },
            { 75, 1.9235539 },
            { 80, 1.9075632 },
            { 85, 1.8869303 },
            { 90, 1.8639250 },
            { 95, 1.8382537 },
            { 100, 1.8140434 },
            { 110, 1.7645023 },
            { 120, 1.7213558 },
            { 130, 1.6821228 },
            { 140, 1.6423432 },
            { 150, 1.6103531 },
            { 160, 1.5834705 },
            { 170, 1.5575980 },
            { 180, 1.5357554 },
            { 190, 1.5140899 },
            { 200, 1.4953936 },
            { 225, 1.4526090 },
            { 240, 1.4349183 },
            { 250, 1.4221813 },
            { 275, 1.3931355 },
            { 300, 1.3696481 },
            { 325, 1.3492310 },
            { 350, 1.3329021 },
            { 375, 1.3168172 },
            { 400, 1.3028021 },
            { 425, 1.2896989 },
            { 450, 1.2783272 },
            { 500, 1.2598832 },
            { 511, 1.2541257 } } };

        // log lin data
        std::vector<std::pair<double, double>> l(data.size());
        std::transform(std::execution::par_unseq, data.cbegin(), data.cend(), l.begin(), [](const auto& p) {
            return std::make_pair(std::log(p.first), p.second);
        });
        m_hp10factor_interpolator = CubicSplineInterpolatorStatic<double, 15>(l);
    }

protected:
    /// @brief Returns {min, max} of two values as a pair.
    static std::pair<double, double> minmax(const double a, const double b)
    {
        return a < b ? std::pair { a, b } : std::pair { b, a };
    }

    /**
     * @brief Tests whether a particle ray intersects the finite sensitive plane.
     *
     * Computes the ray-plane intersection distance @p t, then projects the hit
     * position onto the local x and y axes to verify it falls within the
     * width × height aperture.
     *
     * @param p  Particle defining the ray (position and direction).
     * @return   Distance @p t along the ray to the intersection, or `std::nullopt`
     *           if the ray is parallel to the plane, originates behind it, or
     *           misses the rectangular aperture.
     */
    std::optional<double> intersectPlane(const ParticleType auto& p) const
    {
        const auto n = normalVector();
        const auto dn = vectormath::dot(n, p.dir);
        if (std::abs(dn) < GEOMETRIC_ERROR())
            return std::nullopt;

        auto t = vectormath::dot(n, vectormath::subtract(m_center, p.pos)) / dn;
        if (t < 0.0)
            return std::nullopt;

        auto hitpos = vectormath::add(p.pos, vectormath::scale(p.dir, t));

        auto dist = vectormath::subtract(hitpos, m_center);

        // X
        const auto x = vectormath::dot(dist, m_direction_cosines[0]);
        if (std::abs(x) * 2 > width())
            return std::nullopt;
        // y
        const auto y = vectormath::dot(dist, m_direction_cosines[1]);
        if (std::abs(y) * 2 > height())
            return std::nullopt;
        return t;
    }

    /**
     * @brief Computes the tightest axis-aligned bounding box for the current orientation.
     *
     * Projects the three half-span vectors (x, y, and normal directions scaled by
     * half-width, half-height, and half-depth) onto each world axis and accumulates
     * the min/max extents around the centre.
     *
     * @return AABB as {xmin, ymin, zmin, xmax, ymax, zmax} [cm].
     */
    std::array<double, 6> calculateAABB() const
    {
        const auto normal_vec = normalVector();
        const auto half_span_x = vectormath::scale(m_direction_cosines[0], width() / 2);
        const auto half_span_y = vectormath::scale(m_direction_cosines[1], height() / 2);
        const auto half_span_z = vectormath::scale(normal_vec, depth() / 2);

        const auto span_vectors = std::array { half_span_x, half_span_y, half_span_z };

        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };

        for (const auto& s : span_vectors) {
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                aabb[i] = std::min(aabb[i], mi);
                aabb[i + 3] = std::max(aabb[i + 3], ma);
            }
        }
        return aabb;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 }; ///< Centre of the dosimeter in world coordinates [cm].
    /// In-plane direction cosines {v_x, v_y}; v_x is starboard→port, v_y is bottom→top, normal is back→front.
    std::array<std::array<double, 3>, 2> m_direction_cosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 }; ///< Cached axis-aligned bounding box {xmin,ymin,zmin,xmax,ymax,zmax} [cm].
    Material<12> m_air; ///< Dry-air material used for air-KERMA calculations.
    ParticleTracker m_tracker; ///< Records particle tracks when tracking mode is active.
    EnergyScore m_energyScored; ///< Accumulates Hp(10)-weighted energy per history.
    EnergyScore m_energy_airKermaScored; ///< Accumulates air-KERMA-weighted energy per history.
    DoseScore m_dose; ///< Accumulated Hp(10) dose over all histories.
    DoseScore m_dose_airKermaScored; ///< Accumulated air-KERMA dose over all histories.
    CubicSplineInterpolatorStatic<double, 15> m_hp10factor_interpolator; ///< Log-linear cubic spline for Hp(10) conversion factors.
};
}

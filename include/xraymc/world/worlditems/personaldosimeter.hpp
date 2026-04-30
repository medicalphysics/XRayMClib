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

template <bool USEDIRECTIONALDEPENDENCE = true>
class PersonalDosimeter {
public:
    PersonalDosimeter()
        : m_air(Material<12>::byNistName("Air, Dry (near sea level)").value())
    {
        setupHp10ConvertionFactor();
        m_aabb = calculateAABB();
    }

    /// @brief Returns the axis-aligned bounding box as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        m_aabb = calculateAABB();
    }

    void setDirectionCosines(const std::array<double, 3>& dir_cosX, const std::array<double, 3>& dir_cosY)
    {
        if (vectormath::dot(dir_cosX, dir_cosY) < 2 * GEOMETRIC_ERROR<>()) {
            m_direction_cosines[0] = vectormath::normalized(dir_cosX);
            m_direction_cosines[1] = vectormath::normalized(dir_cosY);
        }
        m_aabb = calculateAABB();
    }

    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_direction_cosines;
    }

    std::array<double, 3> normalVector() const
    {
        return vectormath::cross(m_direction_cosines[0], m_direction_cosines[1]);
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return basicshape::AABB::intersect(p, m_aabb);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        return basicshape::AABB::intersectVisualization<U>(p, m_aabb);
    }

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
            m_energyScored.scoreEnergy(hp10 * p.weight * angular_weight);
            m_energy_airKermaScored.scoreEnergy(airKerma * p.weight * angular_weight);
        } else {
            m_energyScored.scoreEnergy(hp10 * p.weight);
            m_energy_airKermaScored.scoreEnergy(airKerma * p.weight);
        }
        p.energy = 0; // we kill the particle, it is absorbed by the dosimeter
    }

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

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_energyScored;

        return m_energy_airKermaScored;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_energy_airKermaScored.clear();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        constexpr double area = 1.0; // width() * height();
        constexpr double density = 1.0;

        // We have already calculated the dose in energy score.
        m_dose.addScoredEnergy(m_energyScored, area, density, calibration_factor);
        m_dose_airKermaScored.addScoredEnergy(m_energy_airKermaScored, 1.0, 1.0, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_dose;

        return m_dose_airKermaScored;
    }

    void clearDoseScored()
    {
        m_dose.clear();
    }

    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "PersonalDosimeter1";
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

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

    constexpr static double height()
    {
        return 58.0 / 10.0; // cm
    }
    constexpr static double width()
    {
        return 40.0 / 10.0; // cm
    }
    constexpr static double depth()
    {
        return 17.0 / 10.0; // cm
    }

protected:
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
    std::array<double, 3> m_center = { 0, 0, 0 };
    // cosines = {v_x, v_y}, v_x is startboard to port and v_y is bottom to top, normal is back to front
    std::array<std::array<double, 3>, 2> m_direction_cosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 6> m_aabb = { 0, 0, 0 };
    Material<12> m_air;
    ParticleTracker m_tracker;
    EnergyScore m_energyScored;
    EnergyScore m_energy_airKermaScored;
    DoseScore m_dose;
    DoseScore m_dose_airKermaScored;
    CubicSplineInterpolatorStatic<double, 15> m_hp10factor_interpolator;
};
}
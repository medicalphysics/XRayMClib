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

class PersonalDosimeter {
public:
    PersonalDosimeter()
    {
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

        // The model is a simple cosine to the power of 0.2, determined by cognitive fit to data from the RaySafe i3 personal dosimeter.
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
    std::array<std::array<double, 3>, 2> m_direction_cosines = { 1, 0, 0, 0, 1, 0 };
    std::array<double, 6> m_aabb = { 0, 0, 0 };
    ParticleTracker m_tracker;
};
}
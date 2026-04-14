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

#include "xraymc/constants.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"

#include <array>
#include <limits>
#include <optional>

namespace xraymc {
namespace visualizationprops {
    /**
     * @brief Abstract base class for all visualization overlay primitives.
     *
     * Subclasses implement `intersect()` to report the ray parameter and surface
     * normal at the closest intersection point so that VisualizeWorld can composite
     * them on top of the world geometry.
     */
    class VisualizationProp {
    public:
        /**
         * @brief Tests a particle ray against the primitive.
         * @param p Particle whose position and direction define the ray.
         * @return A pair of (t, normal) at the hit point, or std::nullopt on miss.
         *         t is the non-negative ray parameter; normal is the unit outward normal.
         */
        virtual std::optional<std::pair<double, std::array<double, 3>>> intersect(const Particle& p) const noexcept = 0;
    };

    /**
     * @brief A finite cylindrical line prop used to overlay beam paths or particle tracks.
     *
     * Represents a solid cylinder of radius @p radii centered on a ray from @p start
     * along @p dir for at most @p length cm. A negative or zero length is treated as
     * infinite. Intersection is computed analytically as the closest-approach distance
     * between the particle ray and the cylinder axis; the hit is rejected if the
     * closest point falls outside [0, length] along the axis, or if the transverse
     * distance exceeds @p radii.
     */
    class Line : public VisualizationProp {
    public:
        /**
         * @brief Constructs a line prop.
         * @param start  Start position of the cylinder axis in cm.
         * @param dir    Direction of the cylinder axis (normalized internally).
         * @param length Length of the cylinder in cm; ≤ 0 means infinite (default -1).
         * @param radii  Cylinder radius in cm; absolute value is used (default 1).
         */
        Line(const std::array<double, 3> start, const std::array<double, 3> dir, double length = -1, double radii = 1)
            : VisualizationProp()
            , m_pos(start)
            , m_dir(dir)
            , m_radii(std::abs(radii))
        {
            vectormath::normalize(m_dir);
            if (length <= 0)
                m_length = std::numeric_limits<double>::max();
            else
                m_length = length;
        }

        /**
         * @brief Tests a particle ray against the cylinder using closest-approach geometry.
         *
         * Finds the parameter t1 along the particle ray and t2 along the cylinder axis
         * at which the two lines are closest. The hit is rejected when:
         * - the lines are parallel (|dot(ray_dir, axis_dir) − 1| < GEOMETRIC_ERROR),
         * - t1 < 0 (hit is behind the ray origin),
         * - t2 is outside [0, length], or
         * - the closest-approach distance exceeds the cylinder radius.
         * @param p Particle whose position and direction define the ray.
         * @return Pair of (t1, outward_normal) at the hit, or std::nullopt on miss.
         */
        std::optional<std::pair<double, std::array<double, 3>>> intersect(const Particle& p) const noexcept override
        {
            const auto v1v2 = vectormath::dot(p.dir, m_dir);
            if (std::abs(v1v2 - 1) < GEOMETRIC_ERROR()) // parallell lines
                return std::nullopt;

            const auto s = vectormath::subtract(m_pos, p.pos);
            const auto v1s = vectormath::dot(p.dir, s);
            const auto v2s = vectormath::dot(m_dir, s);

            const auto t1 = (v1s - v1v2 * v2s) / (1 - v1v2 * v1v2);

            if (t1 < 0.0)
                return std::nullopt;

            const auto t2 = -v2s + t1 * v1v2;

            if ((t2 < 0.0) || (t2 > m_length))
                return std::nullopt;

            const auto c1 = vectormath::add(p.pos, vectormath::scale(p.dir, t1));
            const auto c2 = vectormath::add(m_pos, vectormath::scale(m_dir, t2));
            const auto d = vectormath::subtract(c1, c2);

            const auto ll = vectormath::length_sqr(d);
            if (ll > m_radii * m_radii)
                return std::nullopt;

            return std::make_optional(std::make_pair(t1, vectormath::normalized(d)));
        }

        /**
         * @brief Returns the axis-aligned bounding box of the cylinder axis segment.
         *
         * The AABB is expanded by GEOMETRIC_ERROR() on each side to ensure the cylinder
         * body is fully contained. For infinite-length lines the extent is unbounded.
         * @return {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
         */
        std::array<double, 6> AABB() const
        {
            const auto end = vectormath::add(m_pos, vectormath::scale(m_dir, m_length));
            std::array<double, 6> aabb;
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(m_pos[i], end[i]) - GEOMETRIC_ERROR();
                aabb[i + 3] = std::max(m_pos[i], end[i]) + GEOMETRIC_ERROR();
            }
            return aabb;
        }

    private:
        std::array<double, 3> m_pos = { 0, 0, 0 };
        std::array<double, 3> m_dir = { 0, 0, 1 };
        double m_radii = 1;
        double m_length = 1;
    };

}

}
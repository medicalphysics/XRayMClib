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

#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"

#include <array>
#include <optional>

namespace xraymc {
namespace basicshape {
    namespace sphere {

        static constexpr bool pointInside(const std::array<double, 3>& pos, const std::array<double, 3>& center, const double radii)
        {
            const std::array<double, 3> dp = { pos[0] - center[0], pos[1] - center[1], pos[2] - center[2] };
            return dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] < radii * radii;
        }

        static constexpr WorldIntersectionResult intersect(const ParticleType auto& p, const std::array<double, 3>& center, const double radii)
        {
            WorldIntersectionResult res;
            // nummeric stable ray sphere intersection
            const auto r2 = radii * radii;
            const auto f = vectormath::subtract(p.pos, center);

            // positive b mean center of sphere is in front of ray
            const auto b = -vectormath::dot(f, p.dir);

            // positive c means ray starts outside of sphere
            const auto c = vectormath::length_sqr(f) - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return res;
            }

            const auto delta1 = vectormath::length_sqr(vectormath::add(f, vectormath::scale(p.dir, b)));

            const auto delta = r2 - delta1;

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return res;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(delta);

            if (c < 0 && b > 0) {
                // inside sphere
                res.rayOriginIsInsideItem = true;
                res.intersection = q;
            } else {
                res.intersection = c / q;
            }
            res.intersectionValid = true;
            return res;
        }

        static constexpr std::optional<std::array<double, 2>> intersectForwardInterval(const ParticleType auto& p, const std::array<double, 3>& center, const double radii)
        {
            // nummeric stable ray sphere intersection
            const auto r2 = radii * radii;
            const auto f = vectormath::subtract(p.pos, center);

            // positive b mean center of sphere is in front of ray
            const auto b = -vectormath::dot(f, p.dir);

            // positive c means ray starts outside of sphere
            const auto c = vectormath::length_sqr(f) - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const auto delta1 = vectormath::length_sqr(vectormath::add(f, vectormath::scale(p.dir, b)));

            const auto delta = r2 - delta1;

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(delta);
            if (c < 0 && b < 0) {
                std::array t = { q, c / q };
                return std::make_optional(t);
            } else {
                std::array t = { c / q, q };
                return std::make_optional(t);
            }
        }

        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p, const std::array<double, 3>& center, const double radii)
        {
            VisualizationIntersectionResult<U> res;
            const auto t_opt = intersectForwardInterval(p, center, radii);
            if (t_opt) {
                const auto& t = t_opt.value();
                res.intersectionValid = true;
                res.rayOriginIsInsideItem = t[0] < 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];

                const auto pos = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                res.normal = vectormath::subtract(center, pos);
                vectormath::normalize(res.normal);
                if (res.rayOriginIsInsideItem) {
                    res.normal = vectormath::scale(res.normal, -1.0);
                }
            }
            return res;
        }
    }
}
}
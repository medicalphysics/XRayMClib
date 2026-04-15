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
    /**
     * @brief Geometric primitives for a sphere.
     *
     * All functions use a numerically stable ray-sphere intersection algorithm
     * based on the `b`/`c`/`q` formulation that avoids catastrophic cancellation.
     * All coordinates and distances are in [cm].
     */
    namespace sphere {

        /**
         * @brief Returns true if @p pos lies strictly inside the sphere.
         *
         * @param pos     3-D point [cm].
         * @param center  Sphere centre [cm].
         * @param radii   Sphere radius [cm].
         * @return True if the squared distance from @p pos to @p center is less than radii².
         */
        static constexpr bool pointInside(const std::array<double, 3>& pos, const std::array<double, 3>& center, const double radii)
        {
            const std::array<double, 3> dp = { pos[0] - center[0], pos[1] - center[1], pos[2] - center[2] };
            return dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] < radii * radii;
        }

        /**
         * @brief Computes the forward ray–sphere intersection for particle transport.
         *
         * Uses a numerically stable quadratic formulation. Returns the nearest forward
         * hit distance; if the ray origin is inside the sphere the exit distance is
         * returned instead.
         *
         * Returns an invalid `WorldIntersectionResult` (`.intersectionValid = false`) if:
         * - The ray origin is outside and the sphere centre is behind the ray.
         * - The discriminant is negative (ray misses).
         *
         * @param p       Particle satisfying `ParticleType` (`.pos`, `.dir`).
         * @param center  Sphere centre [cm].
         * @param radii   Sphere radius [cm].
         * @return `WorldIntersectionResult` with `intersection` = nearest forward hit [cm].
         */
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

        /**
         * @brief Computes both ray–sphere intersection t-values as an interval [t_enter, t_exit].
         *
         * Uses the same numerically stable algorithm as `intersect` but returns both
         * roots so the caller can determine entry and exit distances. t_enter may be
         * negative if the ray origin is inside the sphere.
         *
         * Returns `nullopt` under the same conditions as `intersect` (sphere behind ray,
         * or discriminant < 0).
         *
         * @param p       Particle satisfying `ParticleType`.
         * @param center  Sphere centre [cm].
         * @param radii   Sphere radius [cm].
         * @return `{t_enter, t_exit}` [cm], or `nullopt` if no intersection.
         */
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

        /**
         * @brief Computes a ray–sphere intersection for visualisation, including the surface normal.
         *
         * Delegates to `intersectForwardInterval` for the hit distances, then computes
         * the outward unit normal at the hit point as `normalize(center − hit_pos)`.
         * If the ray origin is inside the sphere the normal is negated (pointing inward,
         * toward the origin).
         *
         * @tparam U      Type tag for the `VisualizationIntersectionResult`.
         * @param p       Particle satisfying `ParticleType`.
         * @param center  Sphere centre [cm].
         * @param radii   Sphere radius [cm].
         * @return `VisualizationIntersectionResult<U>` with intersection distance and outward normal.
         */
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
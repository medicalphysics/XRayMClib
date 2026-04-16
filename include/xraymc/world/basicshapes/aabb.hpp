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

#include <algorithm>
#include <array>
#include <optional>

namespace xraymc {

/**
 * @namespace xraymc::basicshape
 * @brief Geometric primitives.
 *
 * Collection of geometric primitives that are used by items in the World.
 * Implements intersection tests and simple geometric propeties.
 */
namespace basicshape {
    /**
     * @brief Axis-Aligned Bounding Box (AABB) geometry primitives.
     *
     * All functions in this namespace operate on AABBs represented as a 6-element
     * array `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
     */
    namespace AABB {

        /**
         * @brief Returns true if point @p p lies inside or on the boundary of @p aabb.
         *
         * @param p     3-D point [cm].
         * @param aabb  AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @return True if @p p is inside or on the surface of @p aabb.
         */
        static constexpr bool pointInside(const std::array<double, 3>& p, const std::array<double, 6>& aabb)
        {
            return aabb[0] <= p[0] && p[0] <= aabb[3] && aabb[1] <= p[1] && p[1] <= aabb[4] && aabb[2] <= p[2] && p[2] <= aabb[5];
        }

        /**
         * @brief Returns true if two AABBs overlap (share any volume, including touching faces).
         *
         * @param a  First AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @param b  Second AABB.
         * @return True if @p a and @p b overlap or touch.
         */
        static constexpr bool overlap(const std::array<double, 6>& a, const std::array<double, 6>& b)
        {
            const bool x = a[0] <= b[3] && a[3] >= b[0];
            const bool y = a[1] <= b[4] && a[4] >= b[1];
            const bool z = a[2] <= b[5] && a[5] >= b[2];
            return x && y && z;
        }

        /**
         * @brief Returns true if two AABBs collide (overlap in all three dimensions, strictly).
         *
         * Equivalent to `overlap` but uses strict inequality so that touching faces do
         * not count as a collision.
         *
         * @param a  First AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @param b  Second AABB.
         * @return True if @p a and @p b overlap (no gap in any dimension).
         */
        static constexpr bool collide(const std::array<double, 6>& a, const std::array<double, 6>& b)
        {
            // Check if the AABBs are overlapping in all three dimensions
            return !(a[0] > b[3] || a[3] < b[0] || a[1] > b[4] || a[4] < b[1] || a[2] > b[5] || a[5] < b[2]);
        }

        /**
         * @brief Computes the slab-method ray–AABB intersection interval [t_enter, t_exit].
         *
         * Uses the Amy Williams slab method. For each axis where the direction component
         * is non-negligible, the entry and exit t-values are computed; the interval is
         * narrowed to their intersection across all three axes.
         *
         * When `FORWARD = true` (default), t_enter is clamped to ≥ 0 so that only
         * forward intersections are returned. When `FORWARD = false`, both negative
         * (behind the ray origin) and positive t-values are returned — used for
         * visualisation rays that may originate inside the box.
         *
         * @tparam FORWARD  If true, clamp t_enter to 0 (forward-only). Default: true.
         * @param p     Particle (requires `.pos` and `.dir` members satisfying `ParticleType`).
         * @param aabb  AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @return `{t_enter, t_exit}` [cm] if the ray intersects the AABB; `nullopt` otherwise.
         */
        template <bool FORWARD = true>
        std::optional<std::array<double, 2>> intersectForwardInterval(const ParticleType auto& p, const std::array<double, 6>& aabb)
        {
            const std::array<bool, 3> valid = {
                std::abs(p.dir[0]) > std::numeric_limits<double>::epsilon(),
                std::abs(p.dir[1]) > std::numeric_limits<double>::epsilon(),
                std::abs(p.dir[2]) > std::numeric_limits<double>::epsilon(),
            };

            const std::array<double, 3> pdir_inv = {
                valid[0] ? 1.0 / p.dir[0] : 0.0,
                valid[1] ? 1.0 / p.dir[1] : 0.0,
                valid[2] ? 1.0 / p.dir[2] : 0.0
            };

            std::array<double, 2> tm = { std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max() };

            for (int i = 0; i < 3; ++i) {
                if (valid[i]) {
                    const auto t1 = (aabb[i] - p.pos[i]) * pdir_inv[i];
                    const auto t2 = (aabb[i + 3] - p.pos[i]) * pdir_inv[i];

                    tm[0] = std::max(tm[0], std::min(t1, t2));
                    tm[1] = std::min(tm[1], std::max(t1, t2));
                }
            }
            if constexpr (FORWARD)
                tm[0] = std::max(tm[0], 0.0);
            return tm[1] > tm[0] ? std::make_optional(tm) : std::nullopt;
        }

        /**
         * @brief Computes the forward ray–AABB intersection for particle transport.
         *
         * Returns a `WorldIntersectionResult` with:
         * - `intersectionValid = true` if the ray hits the AABB.
         * - `rayOriginIsInsideItem = true` if the particle origin is inside the AABB
         *   (t_enter ≤ 0); in this case `intersection` is the exit distance t_exit.
         * - `intersection` is t_enter if the origin is outside, t_exit if inside [cm].
         *
         * @param p     Particle satisfying `ParticleType`.
         * @param aabb  AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @return `WorldIntersectionResult` describing the intersection.
         */
        static constexpr WorldIntersectionResult intersect(const ParticleType auto& p, const std::array<double, 6>& aabb)
        {
            WorldIntersectionResult res;
            if (const auto t_cand = intersectForwardInterval<true>(p, aabb); t_cand) {
                const auto& t = *t_cand;
                res.rayOriginIsInsideItem = t[0] <= 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
                res.intersectionValid = true;
            }
            return res;
        }

        /**
         * @brief Computes a ray–AABB intersection for visualisation, including the surface normal.
         *
         * Uses `intersectForwardInterval<false>` so both entry and exit are considered
         * regardless of ray origin position. After finding the hit distance, determines
         * the outward face normal by finding which of the six AABB planes is closest to
         * the hit point: the normal component for that face is set to −1 (min face) or
         * +1 (max face), all others remain 0.
         *
         * @tparam U    Type tag for the `VisualizationIntersectionResult` (e.g. material ID).
         * @param p     Particle satisfying `ParticleType`.
         * @param aabb  AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         * @return `VisualizationIntersectionResult<U>` with intersection distance and outward normal.
         */
        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p, const std::array<double, 6>& aabb)
        {
            VisualizationIntersectionResult<U> res;
            const auto t_cand = intersectForwardInterval<false>(p, aabb);
            if (t_cand) {
                const auto& t = *t_cand;
                res.rayOriginIsInsideItem = t[0] <= 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
                res.intersectionValid = true;
                std::array<double, 6> hit_merit;
                for (int i = 0; i < 3; ++i) {
                    const auto hit = p.pos[i] + res.intersection * p.dir[i];
                    hit_merit[i] = std::abs(aabb[i] - hit);
                    hit_merit[i + 3] = std::abs(aabb[i + 3] - hit);
                }
                const auto pos = std::distance(hit_merit.cbegin(), std::min_element(hit_merit.cbegin(), hit_merit.cend()));
                if (pos < 3)
                    res.normal[pos] = -1;
                else
                    res.normal[pos - 3] = 1;
            }
            return res;
        }
    }
}
}
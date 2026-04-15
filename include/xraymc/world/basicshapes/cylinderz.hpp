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
     * @brief Geometric primitives for a z-axis-aligned right-circular cylinder.
     *
     * A specialised, faster alternative to the general `cylinder` namespace for
     * cylinders whose axis is exactly the world z-axis. Parameters are passed
     * directly (centre, radius, half-height) rather than through a struct.
     * All coordinates are in [cm].
     */
    namespace cylinderZ {

        /**
         * @brief Returns true if @p pos lies strictly inside the z-aligned cylinder.
         *
         * The test uses strict inequalities on the z-axis ends and a strict radial
         * inequality, so points exactly on the surface return false.
         *
         * @param pos         3-D point [cm].
         * @param center      Cylinder centre [cm].
         * @param radii       Cylinder radius [cm].
         * @param half_height Half-length along the z-axis [cm].
         * @return True if @p pos is strictly inside the cylinder.
         */
        static constexpr bool pointInside(const std::array<double, 3>& pos, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            const std::array<double, 2> dp = { pos[0] - center[0], pos[1] - center[1] };
            return center[0] - radii <= pos[0] && pos[0] <= center[0] + radii && center[1] - radii <= pos[1] && pos[1] <= center[1] + radii
                && (center[2] - half_height < pos[2]) && (pos[2] < center[2] + half_height) && ((dp[0] * dp[0] + dp[1] * dp[1]) < radii * radii);
        }

        /**
         * @brief Returns the nearest forward ray–cylinder-wall intersection distance [cm].
         *
         * Uses a numerically stable 2-D ray-circle intersection in the xy-plane
         * (ignoring z). Returns the smaller positive t, or the single positive t if
         * the ray origin is inside the circle. Does **not** test the end-caps; the
         * caller must clamp the result to the cylinder's z-extent.
         *
         * Returns `nullopt` if:
         * - The ray direction has no xy-component (parallel to z).
         * - The ray origin is outside and the circle is fully behind the ray.
         * - The ray misses the circle entirely (discriminant < 0).
         *
         * @param p       Particle (only `.pos` and `.dir` are used).
         * @param center  Cylinder centre [cm].
         * @param radii   Cylinder radius [cm].
         * @return Forward intersection distance [cm], or `nullopt` if no hit.
         */
        static constexpr std::optional<double> intersectCylinderWall(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);

            if (c < 0 && b > 0) {
                // inside sphere
                return std::make_optional(q / a);
            } else {
                return std::make_optional(c / q);
            }
        }

        /**
         * @brief Returns the full [t_enter, t_exit] ray–cylinder-wall interval [cm].
         *
         * Uses the same numerically stable 2-D ray-circle formulation as
         * `intersectCylinderWall` but returns both intersection t-values so that the
         * caller can intersect with the z-slab for the end-caps. Both values may be
         * negative (hits behind the ray origin).
         *
         * Returns `nullopt` under the same conditions as `intersectCylinderWall`.
         *
         * @param p       Particle (only `.pos` and `.dir` are used).
         * @param center  Cylinder centre [cm].
         * @param radii   Cylinder radius [cm].
         * @return `{t_enter, t_exit}` [cm], or `nullopt` if no hit.
         */
        static constexpr std::optional<std::array<double, 2>> intersectCylinderWallInterval(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);

            if (c < 0 && b > 0) {
                // inside sphere
                std::array t { c / q, q / a };
                return std::make_optional(t);
            } else {
                std::array t { c / q, q / a };
                return std::make_optional(t);
            }
        }

        /**
         * @brief Returns the ray–disc intersection distance without any forward-direction guard [cm].
         *
         * Tests whether the ray hits a z-plane disc of given @p radii centred at @p center,
         * regardless of whether the hit is in front of or behind the ray origin and
         * regardless of which side the ray approaches from. Used internally when both
         * t-values of a disc are needed for interval computation.
         *
         * Returns `nullopt` if the ray direction has no z-component or the hit point
         * lies outside the disc radius.
         *
         * @param p       Particle (`.pos`, `.dir`).
         * @param center  Disc centre [cm]; only `center[2]` is used for the plane position.
         * @param radii   Disc radius [cm].
         * @return Intersection distance t [cm], or `nullopt` if no hit.
         */
        static std::optional<double> intersectCylinderDiscIntervalZ(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto tz = (center[2] - p.pos[2]) / p.dir[2];
            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz <= radii * radii)
                return std::make_optional(tz);
            return std::nullopt;
        }

        /**
         * @brief Returns the forward ray–disc intersection distance for a z-plane disc [cm].
         *
         * Like `intersectCylinderDiscIntervalZ` but additionally rejects hits where the
         * ray is moving away from the disc plane (i.e. the ray origin is already past the
         * disc in the direction of travel), ensuring only geometrically meaningful forward
         * hits are returned.
         *
         * Returns `nullopt` if:
         * - The ray has no z-component.
         * - The ray is moving away from the disc plane.
         * - The hit point lies outside @p radii.
         *
         * @param p       Particle (`.pos`, `.dir`).
         * @param center  Disc centre [cm]; only `center[2]` sets the plane position.
         * @param radii   Disc radius [cm].
         * @return Forward intersection distance t [cm], or `nullopt` if no hit.
         */
        static std::optional<double> intersectCylinderDiscZ(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<double>::epsilon())
                return std::nullopt;

            if (p.pos[2] > center[2] && p.dir[2] >= 0)
                return std::nullopt;
            if (p.pos[2] < center[2] && p.dir[2] <= 0)
                return std::nullopt;

            const auto tz = (center[2] - p.pos[2]) / p.dir[2];
            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz <= radii * radii)
                return std::make_optional(tz);
            return std::nullopt;
        }

        /**
         * @brief Computes the ray–cylinder intersection for particle transport.
         *
         * Tests the lateral wall (via `intersectCylinderWall`), the bottom cap, and the
         * top cap (via `intersectCylinderDiscZ`) and returns the nearest valid hit.
         * The lateral-wall result is clamped conservatively to the z-extent to avoid
         * spurious misses at the plane–cylinder seam.
         *
         * Returns a `WorldIntersectionResult` with:
         * - `intersectionValid = true` if any of the three tests produced a hit.
         * - `intersection` = the nearest hit distance [cm].
         * - `rayOriginIsInsideItem = true` if the particle starts inside the cylinder.
         *
         * @param p           Particle satisfying `ParticleType`.
         * @param center      Cylinder centre [cm].
         * @param radii       Cylinder radius [cm].
         * @param half_height Half-length along z [cm].
         * @return `WorldIntersectionResult` describing the intersection.
         */
        static WorldIntersectionResult intersect(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            if (t_cand) {
                // we need to be conservative else we will miss intersections on plane cylinder intersection
                const auto tz = std::nextafter(p.pos[2] + p.dir[2] * *t_cand, 0.0);
                if (!(center[2] - half_height <= tz && tz <= center[2] + half_height))
                    t_cand.reset();
            }

            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdisc1_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            centerDisc[2] = center[2] + half_height;
            auto tdisc2_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            const auto t_cand_min = std::min({ t_cand, tdisc1_cand, tdisc2_cand }, [](const auto& lh, const auto& rh) -> bool { return lh.value_or(std::numeric_limits<double>::max()) < rh.value_or(std::numeric_limits<double>::max()); });
            WorldIntersectionResult res;
            if (t_cand_min) {
                res.intersection = *t_cand_min;
                res.intersectionValid = true;
                res.rayOriginIsInsideItem = pointInside(p.pos, center, radii, half_height);
            }

            return res;
        }

        /**
         * @brief Computes the forward ray–cylinder intersection interval [t_enter, t_exit] [cm].
         *
         * Intersects the lateral-wall interval (`intersectCylinderWallInterval`) with the
         * z-slab defined by the two end-caps, then clamps t_enter to ≥ 0. Handles the
         * degenerate case where the ray is exactly aligned with the z-axis by testing
         * radial containment in the xy-plane directly.
         *
         * Returns `nullopt` if the ray misses the cylinder or the entire intersection
         * interval is behind the ray origin.
         *
         * @param p           Particle (`.pos`, `.dir`).
         * @param center      Cylinder centre [cm].
         * @param radii       Cylinder radius [cm].
         * @param half_height Half-length along z [cm].
         * @return `{t_enter, t_exit}` with t_enter ≥ 0 if hit; `nullopt` otherwise.
         */
        static constexpr std::optional<std::array<double, 2>> intersectForwardInterval(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWallInterval(p, center, radii);
            if (!t_cand) {
                if (std::abs(p.dir[2]) == 1) {
                    const auto dx = p.pos[0] - center[0];
                    const auto dy = p.pos[1] - center[1];
                    if (dx * dx + dy * dy < radii * radii) {
                        const auto t1 = std::max(center[2] + half_height - p.pos[2], 0.0);
                        const auto t2 = std::max(center[2] - half_height - p.pos[2], 0.0);
                        std::array<double, 2> tz;
                        if (t1 < t2) {
                            tz = { t1, t2 };
                        } else {
                            tz = { t2, t1 };
                        }
                        if (tz[0] < tz[1]) {
                            return std::make_optional(tz);
                        }
                    }
                }
                return std::nullopt;
            }

            const auto tz1 = (center[2] + half_height - p.pos[2]) / p.dir[2];
            const auto tz2 = (center[2] - half_height - p.pos[2]) / p.dir[2];

            const auto tz_min = tz1 < tz2 ? tz1 : tz2;
            const auto tz_max = tz1 < tz2 ? tz2 : tz1;
            auto& t = t_cand.value();

            t[0] = std::max(t[0], tz_min);
            t[1] = std::min(t[1], tz_max);
            if (t[0] < 0)
                t[0] = 0;
            if (t[0] < t[1])
                return t_cand;
            return std::nullopt;
        }

        /**
         * @brief Computes a ray–cylinder intersection for visualisation, including the surface normal.
         *
         * Tests the lateral wall and both end-caps; picks the nearest hit and sets the
         * outward surface normal:
         * - **Lateral wall**: radial unit vector in the xy-plane pointing away from the axis.
         * - **Top cap** (z = center[2] + half_height): normal = (0, 0, −1).
         * - **Bottom cap** (z = center[2] − half_height): normal = (0, 0, +1).
         *
         * If the ray origin is inside the cylinder the normal is negated (inward-facing).
         *
         * @tparam U      Type tag for the `VisualizationIntersectionResult`.
         * @param p           Particle (`.pos`, `.dir`).
         * @param center      Cylinder centre [cm].
         * @param radii       Cylinder radius [cm].
         * @param half_height Half-length along z [cm].
         * @return `VisualizationIntersectionResult<U>` with intersection distance and outward normal.
         */
        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            if (t_cand) {
                // we need to be conservative else we will miss intersections on plane cylinder intersection
                const auto tz = std::nextafter(p.pos[2] + p.dir[2] * *t_cand, 0.0);
                if (!(center[2] - half_height <= tz && tz <= center[2] + half_height))
                    t_cand.reset();
            }

            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdisc1_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            centerDisc[2] = center[2] + half_height;
            auto tdisc2_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            constexpr auto m = std::numeric_limits<double>::max();

            VisualizationIntersectionResult<U> res;

            if (t_cand.value_or(m) < tdisc1_cand.value_or(m)) {
                if (t_cand.value_or(m) < tdisc2_cand.value_or(m)) {
                    if (t_cand) {
                        res.intersection = *t_cand;
                        res.intersectionValid = true;
                        const auto x = center[0] - (p.pos[0] + res.intersection * p.dir[0]);
                        const auto y = center[1] - (p.pos[1] + res.intersection * p.dir[1]);
                        const auto ll = 1 / std::sqrt(x * x + y * y);
                        res.normal[0] = x * ll;
                        res.normal[1] = y * ll;
                    }
                } else {
                    if (tdisc2_cand) {
                        res.intersection = *tdisc2_cand;
                        res.intersectionValid = true;
                        res.normal[2] = -1;
                    }
                }
            } else {
                if (tdisc2_cand.value_or(m) < tdisc1_cand.value_or(m)) {
                    if (tdisc2_cand) {
                        res.intersection = *tdisc2_cand;
                        res.intersectionValid = true;
                        res.normal[2] = -1;
                    }
                } else {
                    if (tdisc1_cand) {
                        res.intersection = *tdisc1_cand;
                        res.intersectionValid = true;
                        res.normal[2] = 1;
                    }
                }
            }

            if (res.intersectionValid) {
                res.rayOriginIsInsideItem = pointInside(p.pos, center, radii, half_height);
            }
            if (res.rayOriginIsInsideItem) {
                res.normal = vectormath::scale(res.normal, -1.0);
            }
            return res;
        }
    }
}
}
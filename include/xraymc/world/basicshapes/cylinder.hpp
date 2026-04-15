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
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"

#include <array>
#include <optional>

namespace xraymc {
namespace basicshape {
    /**
     * @brief Geometric primitives for a finite right-circular cylinder.
     *
     * The cylinder is parameterised by a centre point, a unit axis direction,
     * a radius, and a half-height. All coordinates are in [cm].
     */
    namespace cylinder {

        /**
         * @brief Geometric description of a finite right-circular cylinder.
         *
         * The cylinder extends from `center − direction × half_height` to
         * `center + direction × half_height` along its axis, with the given radius.
         * The axis direction is always normalised at construction time.
         */
        struct Cylinder {
            std::array<double, 3> center = { 0, 0, 0 };    ///< Centre point of the cylinder [cm].
            std::array<double, 3> direction = { 0, 0, 1 };  ///< Unit vector along the cylinder axis.
            double radius = 1;                               ///< Cylinder radius [cm].
            double half_height = 1;                          ///< Half-length along the axis [cm].

            /// @brief Default-constructs a unit cylinder centred at the origin, aligned with the z-axis.
            Cylinder() = default;

            /**
             * @brief Constructs a cylinder with explicit parameters.
             *
             * @p direction_arr is normalised internally.
             *
             * @param center_arr        Centre of the cylinder [cm].
             * @param direction_arr     Axis direction (normalised internally).
             * @param radii             Cylinder radius [cm].
             * @param half_height_wall  Half-length along the axis [cm].
             */
            Cylinder(const std::array<double, 3>& center_arr, const std::array<double, 3>& direction_arr, double radii, double half_height_wall)
                : center(center_arr)
                , radius(radii)
                , half_height(half_height_wall)
            {
                direction = vectormath::normalized(direction_arr);
            }

            /// @brief Equality comparison (all fields).
            bool operator==(const Cylinder&) const = default;

            /**
             * @brief Returns the volume of the cylinder [cm³].
             * @return π × radius² × 2 × half_height [cm³].
             */
            double volume() const
            {
                return std::numbers::pi_v<double> * radius * radius * half_height * 2;
            }
        };

        /**
         * @brief Computes the axis-aligned bounding box of a cylinder.
         *
         * The AABB is computed from the two end-disc extents and the disc radii
         * projected onto each world axis.
         *
         * @param cyl  The cylinder to bound.
         * @return AABB as `{x_min, y_min, z_min, x_max, y_max, z_max}` [cm].
         */
        static std::array<double, 6> cylinderAABB(const Cylinder& cyl)
        {
            // calculating disc extents
            const std::array<double, 3> e = {
                cyl.radius * std::sqrt(1 - cyl.direction[0] * cyl.direction[0]),
                cyl.radius * std::sqrt(1 - cyl.direction[1] * cyl.direction[1]),
                cyl.radius * std::sqrt(1 - cyl.direction[2] * cyl.direction[2])
            };
            const auto l = vectormath::scale(cyl.direction, cyl.half_height);
            const auto p0 = vectormath::subtract(cyl.center, l);
            const auto p1 = vectormath::add(cyl.center, l);
            std::array<double, 6> aabb = {
                std::min(p0[0], p1[0]) - e[0],
                std::min(p0[1], p1[1]) - e[1],
                std::min(p0[2], p1[2]) - e[2],
                std::max(p0[0], p1[0]) + e[0],
                std::max(p0[1], p1[1]) + e[1],
                std::max(p0[2], p1[2]) + e[2]
            };
            return aabb;
        }

        /**
         * @brief Returns true if @p pos lies inside or on the surface of @p cylinder.
         *
         * Tests radial containment within the infinite cylinder first, then checks
         * that the point lies between the two end-caps.
         *
         * @param pos       3-D point [cm].
         * @param cylinder  Cylinder to test against.
         * @return True if @p pos is inside or on the surface of @p cylinder.
         */
        static constexpr bool pointInside(const std::array<double, 3>& pos, const Cylinder& cylinder)
        {
            const auto d = vectormath::cross(cylinder.direction, vectormath::subtract(pos, cylinder.center));
            // test if point inside infinite cylinder
            if (vectormath::length_sqr(d) <= cylinder.radius * cylinder.radius) {
                // test for side of two end planes
                const auto e = vectormath::scale(cylinder.direction, cylinder.half_height);
                const auto p0 = vectormath::subtract(cylinder.center, e);

                if (vectormath::dot(vectormath::subtract(pos, p0), cylinder.direction) >= 0) {
                    // (pos - p)*normal >= 0
                    const auto p1 = vectormath::add(cylinder.center, e);
                    if (vectormath::dot(vectormath::subtract(pos, p1), cylinder.direction) <= 0) {
                        return true;
                    }
                }
            }
            return false;
        }

        /*static std::array<std::optional<double>, 2> intersectDisc(const ParticleType auto& p, const Cylinder& cylinder)
        {
            std::array<std::optional<double>, 2> res;

            const auto denom = vectormath::dot(p.dir, cylinder.direction);
            const auto r2 = cylinder.radius * cylinder.radius;

            const auto c1 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, cylinder.half_height));
            const auto dist1 = vectormath::subtract(c1, p.pos);
            const auto t1 = vectormath::dot(dist1, cylinder.direction) / denom;
            const auto hit_dist1 = vectormath::length_sqr(vectormath::subtract(vectormath::add(p.pos, vectormath::scale(p.dir, t1)), c1));
            res[0] = hit_dist1 <= r2 ? std::make_optional(t1) : std::nullopt;

            const auto c2 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, -cylinder.half_height));
            const auto dist2 = vectormath::subtract(c2, p.pos);
            const auto t2 = vectormath::dot(dist2, cylinder.direction) / denom;
            const auto hit_dist2 = vectormath::length_sqr(vectormath::subtract(vectormath::add(p.pos, vectormath::scale(p.dir, t2)), c2));
            res[1] = hit_dist2 <= r2 ? std::make_optional(t2) : std::nullopt;

            return res;
        }

        static std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            std::array<double, 2> t = { 0, 0 };

            const auto na = vectormath::cross(p.dir, cylinder.direction);

            if (vectormath::length_sqr(na) < GEOMETRIC_ERROR<>()) {
                // only test caps
                const auto res = intersectDisc(p, cylinder);
                if (res[0] && res[1]) {
                    t[0] = std::min(*res[0], *res[1]);
                    t[1] = std::max(*res[0], *res[1]);
                    return t;
                }
                return std::nullopt;
            }

            const auto b = vectormath::subtract(cylinder.center, p.pos);
            const auto ba = vectormath::cross(b, cylinder.direction);
            const auto naba = vectormath::dot(na, ba);
            const auto nana = vectormath::dot(na, na);
            const auto bna = vectormath::dot(b, na);

            const auto r2 = cylinder.radius * cylinder.radius;
            const auto rot_term = nana * r2 - bna * bna;
            if (rot_term <= 0.0)
                return std::nullopt;

            const auto term = std::sqrt(rot_term);
            const auto nana_inv = 1 / nana;

            // cylinder intersects
            const auto ct1 = (naba + term) * nana_inv;
            const auto ct2 = (naba - term) * nana_inv;
            t[0] = ct1 < ct2 ? ct1 : ct2;
            t[1] = ct1 > ct2 ? ct1 : ct2;

            const auto y1 = vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[0]), b));
            const auto y2 = vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[1]), b));
            if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height && -cylinder.half_height <= y2 && y2 <= cylinder.half_height) {
                // we hit only cylinder return
                return t;
            } else if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height) {
                const auto res = intersectDisc(p, cylinder);
                t[1] = res[0] ? *res[0] : *res[1];
                return t;
            } else if (-cylinder.half_height <= y2 && y2 <= cylinder.half_height) {
                const auto res = intersectDisc(p, cylinder);
                t[0] = res[0] ? *res[0] : *res[1];
                return t;
            } else {
                const auto res = intersectDisc(p, cylinder);
                if (res[0] && res[1]) {
                    t[0] = std::min(*res[0], *res[1]);
                    t[1] = std::max(*res[0], *res[1]);
                    return t;
                }
            }
            return std::nullopt;
        }*/

        /**
         * @brief Computes the full (forward and backward) ray–cylinder intersection interval.
         *
         * Handles two cases:
         * - **Parallel ray** (`|na|² < ε`): the ray is parallel to the cylinder axis; only
         *   end-cap intersections are tested.
         * - **General ray**: the lateral surface intersection is computed from the
         *   cross-product formulation, then each intersection t-value is replaced by
         *   the corresponding cap t-value if the lateral hit is outside the end-caps.
         *
         * Returns `nullopt` if the ray misses the cylinder entirely.
         *
         * @param p        Particle satisfying `ParticleType` (`.pos`, `.dir`).
         * @param cylinder Cylinder to intersect.
         * @return `{t_enter, t_exit}` [cm] if the ray hits; `nullopt` otherwise.
         *         Both t-values may be negative (ray hits are behind the origin).
         */
        static std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto b = vectormath::subtract(cylinder.center, p.pos);
            const auto na = vectormath::cross(p.dir, cylinder.direction);
            const auto nana = vectormath::length_sqr(na);
            const auto r2 = cylinder.radius * cylinder.radius;
            if (nana < GEOMETRIC_ERROR<double>()) {
                std::array<double, 2> t;
                const auto nad = vectormath::dot(cylinder.direction, p.dir);
                const auto c0 = vectormath::add(b, vectormath::scale(cylinder.direction, -cylinder.half_height));
                t[0] = vectormath::dot(cylinder.direction, c0) / nad;
                const auto ndc1 = vectormath::subtract(vectormath::scale(p.dir, t[0]), c0);
                if (vectormath::dot(ndc1, ndc1) < r2) {
                    const auto c1 = vectormath::add(b, vectormath::scale(cylinder.direction, cylinder.half_height));
                    t[1] = vectormath::dot(cylinder.direction, c1) / nad;
                    const auto ndc2 = vectormath::subtract(vectormath::scale(p.dir, t[1]), c1);
                    if (vectormath::dot(ndc2, ndc2) < r2) {
                        if (nad < 0.0)
                            std::swap(t[0], t[1]);
                        return t;
                    }
                }
                return std::nullopt;

            } else {
                const auto bna = vectormath::dot(b, na);
                const auto broot = nana * r2 - bna * bna;
                if (broot <= 0.0) {
                    // no intersection at all
                    return std::nullopt;
                }
                const auto naba = vectormath::dot(na, vectormath::cross(b, cylinder.direction));
                const auto broot_r = std::sqrt(broot);
                std::array<double, 2> t = {
                    (naba - broot_r) / nana,
                    (naba + broot_r) / nana
                };

                const std::array<double, 2> y = {
                    vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[0]), b)),
                    vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[1]), b))
                };

                const auto nad = vectormath::dot(p.dir, cylinder.direction);

                if (-cylinder.half_height <= y[0] && y[0] <= cylinder.half_height && -cylinder.half_height <= y[1] && y[1] <= cylinder.half_height) {
                    // we are inside cylinder
                    return t;
                } else if ((y[0] <= -cylinder.half_height && y[1] >= cylinder.half_height) || (y[1] <= -cylinder.half_height && y[0] >= cylinder.half_height)) {
                    // we are outside but hits both caps
                    if (std::abs(nad) > GEOMETRIC_ERROR<double>()) {
                        const auto nd_inv = 1.0 / nad;
                        const auto c0 = vectormath::add(b, vectormath::scale(cylinder.direction, -cylinder.half_height));
                        const auto c1 = vectormath::add(b, vectormath::scale(cylinder.direction, cylinder.half_height));
                        if (nad >= 0) {
                            t[0] = vectormath::dot(cylinder.direction, c0) * nd_inv;
                            t[1] = vectormath::dot(cylinder.direction, c1) * nd_inv;
                        } else {
                            t[1] = vectormath::dot(cylinder.direction, c0) * nd_inv;
                            t[0] = vectormath::dot(cylinder.direction, c1) * nd_inv;
                        }
                        return t;
                    }
                } else if (-cylinder.half_height <= y[0] && y[0] <= cylinder.half_height) {
                    // cap upper
                    if (std::abs(nad) > GEOMETRIC_ERROR<double>()) {
                        const int side = nad >= 0 ? 1 : -1;
                        const auto c = vectormath::add(b, vectormath::scale(cylinder.direction, side * cylinder.half_height));
                        t[1] = vectormath::dot(cylinder.direction, c) / nad;
                        return t;
                    }
                } else if (-cylinder.half_height <= y[1] && y[1] <= cylinder.half_height) {
                    // cap lower
                    if (std::abs(nad) > GEOMETRIC_ERROR<double>()) {
                        const int side = nad >= 0 ? -1 : 1;
                        const auto c = vectormath::add(b, vectormath::scale(cylinder.direction, side * cylinder.half_height));
                        t[0] = vectormath::dot(cylinder.direction, c) / nad;
                        return t;
                    }
                }
                return std::nullopt;
            }
        }

        /**
         * @brief Computes the forward-only ray–cylinder intersection interval.
         *
         * Delegates to `intersectInterval` then clamps t_enter to ≥ 0 and discards
         * the result entirely if t_exit ≤ 0 (cylinder is fully behind the ray origin).
         *
         * @param p        Particle satisfying `ParticleType`.
         * @param cylinder Cylinder to intersect.
         * @return `{t_enter, t_exit}` with t_enter ≥ 0 if the ray hits; `nullopt` otherwise.
         */
        static std::optional<std::array<double, 2>> intersectForwardInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            auto t_cand = intersectInterval(p, cylinder);
            if (t_cand) {
                auto& v = t_cand.value();
                if (v[1] <= 0)
                    return std::nullopt;
                v[0] = std::max(v[0], 0.0);
            }
            return t_cand;
        }

        /**
         * @brief Computes the ray–cylinder intersection for particle transport.
         *
         * Returns a `WorldIntersectionResult` with:
         * - `intersectionValid = true` if the ray hits the cylinder and t_exit > 0.
         * - `rayOriginIsInsideItem = true` if the particle origin is inside the cylinder
         *   (t_enter ≤ 0); in that case `intersection` = t_exit (the exit distance).
         * - `intersection` = t_enter if the origin is outside [cm].
         *
         * @param p        Particle satisfying `ParticleType`.
         * @param cylinder Cylinder to intersect.
         * @return `WorldIntersectionResult` describing the intersection.
         */
        static WorldIntersectionResult intersect(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto t_cand = intersectInterval(p, cylinder);
            WorldIntersectionResult res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > 0) {
                    res.rayOriginIsInsideItem = v[0] <= 0;
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;
                }
            }
            return res;
        }

        /**
         * @brief Returns the outward surface normal at a point on the cylinder surface.
         *
         * Determines whether @p pos is on a flat end-cap or on the lateral surface by
         * comparing the radial distance from the axis to `radius × (1 − ε)`:
         * - **End-cap**: returns `±direction` depending on which cap was hit.
         * - **Lateral surface**: returns the normalised radial vector from the axis to @p pos.
         *
         * @param pos      Point on the cylinder surface [cm].
         * @param cylinder The cylinder.
         * @return Outward unit normal vector at @p pos.
         */
        static std::array<double, 3> normalOnPoint(const std::array<double, 3>& pos, const Cylinder& cylinder)
        {
            // finding normal
            const auto pa = vectormath::subtract(pos, cylinder.center);
            const auto ns = vectormath::dot(pa, cylinder.direction);
            const auto d = vectormath::subtract(pa, vectormath::scale(cylinder.direction, ns));
            const auto r = cylinder.radius * (1 - GEOMETRIC_ERROR());
            std::array<double, 3> normal;
            if (vectormath::length_sqr(d) < r * r) {
                // we hit plane
                if (vectormath::dot(pa, cylinder.direction) < 0) {
                    normal = vectormath::scale(cylinder.direction, -1.0);
                } else {
                    normal = cylinder.direction;
                }
            } else {
                normal = vectormath::normalized(d);
            }
            return normal;
        }

        /**
         * @brief Computes a ray–cylinder intersection for visualisation, including the surface normal.
         *
         * Uses `intersectInterval` (no forward clamping) and returns both the hit distance
         * and the outward surface normal at the hit point. If the ray origin is inside the
         * cylinder the normal is negated (pointing inward, toward the origin).
         *
         * @tparam U    Type tag for the `VisualizationIntersectionResult`.
         * @param p        Particle satisfying `ParticleType`.
         * @param cylinder Cylinder to intersect.
         * @return `VisualizationIntersectionResult<U>` with intersection distance and outward normal.
         */
        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto t_cand = intersectInterval(p, cylinder);
            VisualizationIntersectionResult<U> res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > 0) {
                    res.rayOriginIsInsideItem = v[0] <= 0;
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;

                    // finding normal
                    const auto p0 = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                    res.normal = normalOnPoint(p0, cylinder);
                    if (res.rayOriginIsInsideItem) {
                        res.normal = vectormath::scale(res.normal, -1.0);
                    }
                }
            }
            return res;
        }
    }
}
}
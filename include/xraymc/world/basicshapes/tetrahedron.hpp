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
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"

#include <array>
#include <numeric>
#include <optional>

namespace xraymc {
namespace basicshape {
    namespace tetrahedron {

        /**
         * @brief Computes the normal vector of the triangle defined by three points.
         * @tparam NORMALIZE If true (default), returns a unit normal; otherwise returns the unnormalized cross product.
         * @param p0 First vertex of the triangle.
         * @param p1 Second vertex of the triangle.
         * @param p2 Third vertex of the triangle.
         * @return Normal vector of the triangle, optionally normalized.
         */
        template <bool NORMALIZE = true>
        static std::array<double, 3> normalVector(const std::array<double, 3>& p0, const std::array<double, 3>& p1, const std::array<double, 3>& p2)
        {
            const auto s1 = vectormath::subtract(p1, p0);
            const auto s2 = vectormath::subtract(p2, p0);
            if constexpr (NORMALIZE)
                return vectormath::normalized(vectormath::cross(s2, s1));
            else
                return vectormath::cross(s2, s1);
        }

        /**
         * @brief Computes the axis-aligned bounding box (AABB) of a tetrahedron.
         * @param v0 First vertex.
         * @param v1 Second vertex.
         * @param v2 Third vertex.
         * @param v3 Fourth vertex.
         * @return Array of six values [xmin, ymin, zmin, xmax, ymax, zmax].
         */
        static constexpr std::array<double, 6> AABB(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3)
        {
            std::array<double, 6> aabb;
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min({ v0[i], v1[i], v2[i], v3[i] });
                aabb[i + 3] = std::max({ v0[i], v1[i], v2[i], v3[i] });
            }
            return aabb;
        }

        /**
         * @brief Computes the axis-aligned bounding box (AABB) of a triangle.
         * @param v0 First vertex.
         * @param v1 Second vertex.
         * @param v2 Third vertex.
         * @return Array of six values [xmin, ymin, zmin, xmax, ymax, zmax].
         */
        static constexpr std::array<double, 6> triangleAABB(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2)
        {
            std::array<double, 6> aabb;
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min({ v0[i], v1[i], v2[i] });
                aabb[i + 3] = std::max({ v0[i], v1[i], v2[i] });
            }
            return aabb;
        }

        /**
         * @brief Tests whether a point lies inside a tetrahedron.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @param point The point to test.
         * @return True if the point is inside or on the boundary of the tetrahedron.
         */
        static bool pointInside(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const std::array<double, 3>& point)
        {
            // F3 (V0V1V2), F2 (V1V0V3), F1 (V2V3V0), F0 (V3V2V1)
            const std::array<std::array<double, 3>, 4> normals = {
                normalVector<false>(v0, v1, v2),
                normalVector<false>(v1, v0, v3),
                normalVector<false>(v2, v3, v0),
                normalVector<false>(v3, v2, v1)
            };

            bool inside = true;
            inside = inside && vectormath::dot(vectormath::subtract(v0, point), normals[0]) >= 0;
            inside = inside && vectormath::dot(vectormath::subtract(v1, point), normals[1]) >= 0;
            inside = inside && vectormath::dot(vectormath::subtract(v2, point), normals[2]) >= 0;
            inside = inside && vectormath::dot(vectormath::subtract(v3, point), normals[3]) >= 0;
            return inside;
        }

        /**
         * @brief Tests whether a triangle overlaps an axis-aligned bounding box.
         *
         * Uses the separating axis theorem approach from Schwarz & Seidel (SIGA 2010).
         * @param v0 First vertex of the triangle.
         * @param v1 Second vertex of the triangle.
         * @param v2 Third vertex of the triangle.
         * @param aabb The AABB as [xmin, ymin, zmin, xmax, ymax, zmax].
         * @return True if the triangle overlaps the AABB.
         */
        static bool triangleAABBoverlap(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 6>& aabb)
        {
            // test for a triangle overlapping a AABB
            // https://michael-schwarz.com/research/publ/files/vox-siga10.pdf
            const auto normal = normalVector(v0, v1, v2);

            const auto& [p, pd] = vectormath::splice(aabb);
            const auto d_aabb = vectormath::subtract(pd, p);

            const std::array c = {
                normal[0] > 0 ? d_aabb[0] : 0.0,
                normal[1] > 0 ? d_aabb[1] : 0.0,
                normal[2] > 0 ? d_aabb[2] : 0.0,
            };

            const auto d1 = vectormath::dot(normal, vectormath::subtract(c, v0));
            const auto d2 = vectormath::dot(normal, vectormath::subtract(vectormath::subtract(d_aabb, c), v0));

            const auto np = vectormath::dot(normal, p);
            const bool overlap = (np + d1) * (np + d2) <= 0;
            if (overlap) {
                const std::array<std::array<double, 3>, 3> edges = {
                    vectormath::subtract(v1, v0),
                    vectormath::subtract(v2, v1),
                    vectormath::subtract(v0, v2)
                };
                const std::array<std::array<double, 3>, 3> vertices = {
                    v0,
                    v1,
                    v2
                };

                // XY
                bool overlap_xy = true;
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto& e = edges[i];
                    const auto& v = vertices[i];
                    std::array n_e = { -e[1], e[0] };
                    if (normal[2] < 0) {
                        n_e[0] *= -1;
                        n_e[1] *= -1;
                    }

                    const auto d_e = -(n_e[0] * v[0] + n_e[1] * v[1]) + std::max(0.0, d_aabb[0] * n_e[0]) + std::max(0.0, d_aabb[1] * n_e[1]);
                    const auto res_xy = (n_e[0] * p[0] + n_e[1] * p[1] + d_e);
                    overlap_xy = overlap_xy && res_xy >= 0;
                }

                // XZ
                bool overlap_xz = true;
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto& e = edges[i];
                    const auto& v = vertices[i];
                    std::array n_e = { -e[2], e[0] };
                    if (normal[1] < 0) {
                        n_e[0] *= -1;
                        n_e[1] *= -1;
                    }

                    const auto d_e = -(n_e[0] * v[0] + n_e[1] * v[2]) + std::max(0.0, d_aabb[0] * n_e[0]) + std::max(0.0, d_aabb[2] * n_e[1]);
                    const auto res_xz = (n_e[0] * p[0] + n_e[1] * p[2] + d_e);
                    overlap_xz = overlap_xz && res_xz >= 0;
                }
                // YZ
                bool overlap_yz = true;
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto& e = edges[i];
                    const auto& v = vertices[i];
                    std::array n_e = { -e[2], e[1] };
                    if (normal[1] < 0) {
                        n_e[0] *= -1;
                        n_e[1] *= -1;
                    }

                    const auto d_e = -(n_e[0] * v[1] + n_e[1] * v[2]) + std::max(0.0, d_aabb[1] * n_e[0]) + std::max(0.0, d_aabb[2] * n_e[1]);
                    const auto res_yz = (n_e[0] * p[1] + n_e[1] * p[2] + d_e);
                    overlap_yz = overlap_yz && res_yz >= 0;
                }
                return overlap_xy || overlap_xz || overlap_yz;
            }
            return overlap;
        }

        /**
         * @brief Tests whether a tetrahedron overlaps an axis-aligned bounding box.
         *
         * First performs a coarse AABB-AABB overlap test, then checks each of the four
         * triangular faces against the AABB for a precise result.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @param aabb The AABB as [xmin, ymin, zmin, xmax, ymax, zmax].
         * @return True if any face of the tetrahedron overlaps the AABB.
         */
        static constexpr bool insideAABB(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const std::array<double, 6>& aabb)
        {

            const std::array<double, 6> tetAABB = {
                std::min({ v0[0], v1[0], v2[0], v3[0] }),
                std::min({ v0[1], v1[1], v2[1], v3[1] }),
                std::min({ v0[2], v1[2], v2[2], v3[2] }),
                std::max({ v0[0], v1[0], v2[0], v3[0] }),
                std::max({ v0[1], v1[1], v2[1], v3[1] }),
                std::max({ v0[2], v1[2], v2[2], v3[2] })
            };
            const bool overlap = AABB::overlap(aabb, tetAABB);

            if (overlap) {
                // test if the 4 triangles intersect the AABB
                const bool tri_overlap1 = triangleAABBoverlap(v0, v1, v2, aabb);
                const bool tri_overlap2 = triangleAABBoverlap(v1, v0, v3, aabb);
                const bool tri_overlap3 = triangleAABBoverlap(v2, v3, v0, aabb);
                const bool tri_overlap4 = triangleAABBoverlap(v3, v2, v1, aabb);
                const bool tri_overlap = tri_overlap1 || tri_overlap2 || tri_overlap3 || tri_overlap4;
                return tri_overlap;

                /* Is this needed, will we miss tetrahedrons covering a whole grid element?
                if (!tri_overlap) { // if no intersection test for a point in aabb is inside tetrahedron
                    const std::array point = { aabb[0], aabb[1], aabb[2] };
                    return pointInside(v0, v1, v2, v3, point);
                } else
                    return true;
                */
            }
            return overlap;
        }

        /**
         * @brief Computes the volume of a tetrahedron.
         * @param v0 First vertex.
         * @param v1 Second vertex.
         * @param v2 Third vertex.
         * @param v3 Fourth vertex.
         * @return Volume of the tetrahedron (always non-negative).
         */
        static double volume(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3)
        {
            const auto a = vectormath::subtract(v1, v0);
            const auto b = vectormath::subtract(v2, v0);
            const auto c = vectormath::subtract(v3, v0);
            static constexpr double scale = 1.0 / 6.0;
            return scale * std::abs(vectormath::tripleProduct(a, b, c));
        }

        /**
         * @brief Returns the outward face normal of the tetrahedron face closest to a given point.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @param p The query point.
         * @return Unit normal of the face nearest to @p p.
         */
        static std::array<double, 3> closestNormalToPoint(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const std::array<double, 3>& p)
        {
            const std::array<std::array<double, 3>, 4> normals = {
                normalVector<true>(v0, v1, v2),
                normalVector<true>(v1, v0, v3),
                normalVector<true>(v2, v3, v0),
                normalVector<true>(v3, v2, v1)
            };

            const std::array<double, 4> dist = {
                std::abs(vectormath::dot(normals[0], vectormath::subtract(p, v0))),
                std::abs(vectormath::dot(normals[1], vectormath::subtract(p, v1))),
                std::abs(vectormath::dot(normals[2], vectormath::subtract(p, v2))),
                std::abs(vectormath::dot(normals[3], vectormath::subtract(p, v3))),
            };
            auto ind = std::min_element(dist.cbegin(), dist.cend());
            return normals[std::distance(dist.cbegin(), ind)];
        }

        /**
         * @brief Verifies that all four face normals of a tetrahedron point outward.
         *
         * For each face, checks that the opposite vertex lies on the negative side of
         * the face's normal, which is the expected convention for outward-facing normals.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @return True if all face normals point outward.
         */
        static bool tetrahedronNormalsOutsideTest(
            const std::array<double, 3>& v0,
            const std::array<double, 3>& v1,
            const std::array<double, 3>& v2,
            const std::array<double, 3>& v3)
        {

            const std::array<std::array<double, 3>, 4> normals = {
                normalVector<false>(v0, v1, v2),
                normalVector<false>(v1, v0, v3),
                normalVector<false>(v2, v3, v0),
                normalVector<false>(v3, v2, v1)
            };

            const std::array<double, 4> d = {
                vectormath::dot(normals[0], vectormath::subtract(v3, v0)),
                vectormath::dot(normals[1], vectormath::subtract(v2, v1)),
                vectormath::dot(normals[2], vectormath::subtract(v1, v2)),
                vectormath::dot(normals[3], vectormath::subtract(v0, v3)),
            };

            return d[0] < 0 && d[1] < 0 && d[2] < 0 && d[3] < 0;
        }

        /**
         * @brief Computes the ray-triangle intersection distance using the Möller–Trumbore algorithm.
         * @param v0 First vertex of the triangle.
         * @param v1 Second vertex of the triangle.
         * @param v2 Third vertex of the triangle.
         * @param p The particle (ray), providing position and direction.
         * @return The parametric distance @c t along the ray to the intersection, or @c std::nullopt if there is no intersection.
         */
        static std::optional<double> intersectTriangle(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const ParticleType auto& p)
        {
            const auto E1 = vectormath::subtract(v1, v0);
            const auto TT = vectormath::subtract(p.pos, v0);
            const auto Q = vectormath::cross(TT, E1);

            const auto E2 = vectormath::subtract(v2, v0);
            const auto P = vectormath::cross(p.dir, E2);

            const auto det = vectormath::dot(P, E1);
            if (std::abs(det) < GEOMETRIC_ERROR())
                return std::nullopt;

            const auto det_inv = 1.0 / det;

            const auto v = vectormath::dot(Q, p.dir) * det_inv;
            const auto u = vectormath::dot(P, TT) * det_inv;

            return v >= 0 && u >= 0 && (u + v) <= 1 && u <= 1 ? std::make_optional(vectormath::dot(Q, E2) * det_inv)
                                                              : std::nullopt;

            /*
            const auto edge1 = vectormath::subtract(v1, v0);
            const auto edge2 = vectormath::subtract(v2, v0);

            const auto h = vectormath::cross(p.dir, edge2);
            const auto a = vectormath::dot(edge1, h);

            // If a is near zero, ray is parallel to triangle
            if (std::abs(a) < GEOMETRIC_ERROR())
                return std::nullopt;

            const auto f = 1.0 / a;
            const auto s = vectormath::subtract(p.pos, v0);
            const auto u = f * vectormath::dot(s, h);

            // Check if intersection is outside triangle
            if (u < 0.0 || u > 1.0)
                return std::nullopt;

            const auto q = vectormath::cross(s, edge1);
            const auto v = f * vectormath::dot(p.dir, q);

            if (v < 0.0 || u + v > 1.0)
                return std::nullopt;

            // Compute t to find intersection point
            double t = f * vectormath::dot(edge2, q);

            return t;
            */
        }

        /**
         * @brief Computes the ray-tetrahedron intersection.
         *
         * Tests the ray against all four triangular faces and returns the closest
         * positive-distance hit, along with whether the ray origin is inside the tetrahedron.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @param particle The particle (ray), providing position and direction.
         * @return A @c WorldIntersectionResult describing the intersection.
         */
        static WorldIntersectionResult intersect(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const ParticleType auto& particle)
        {

            const std::array<std::optional<double>, 4> hits = {
                intersectTriangle(v0, v1, v2, particle),
                intersectTriangle(v1, v0, v3, particle),
                intersectTriangle(v2, v3, v0, particle),
                intersectTriangle(v3, v2, v1, particle)
            };

            auto chit = std::numeric_limits<double>::max();
            for (const auto& h : hits) {
                if (h) {
                    if (*h > 0)
                        chit = std::min(*h, chit);
                }
            }

            WorldIntersectionResult res;
            if (chit < std::numeric_limits<double>::max()) {
                res.rayOriginIsInsideItem = pointInside(v0, v1, v2, v3, particle.pos);
                res.intersection = chit;
                res.intersectionValid = true;
            }
            return res;
        }

        /**
         * @brief Computes the ray-tetrahedron intersection for visualization, including the surface normal at the hit.
         *
         * Like @c intersect, but also computes the outward face normal at the closest intersection
         * point for use in rendering and visualization.
         * @tparam U The scalar type used by @c VisualizationIntersectionResult.
         * @param v0 First vertex of the tetrahedron.
         * @param v1 Second vertex of the tetrahedron.
         * @param v2 Third vertex of the tetrahedron.
         * @param v3 Fourth vertex of the tetrahedron.
         * @param particle The particle (ray), providing position and direction.
         * @return A @c VisualizationIntersectionResult describing the intersection and surface normal.
         */
        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const ParticleType auto& particle)
        {

            const std::array<std::optional<double>, 4> hits = {
                intersectTriangle(v0, v1, v2, particle),
                intersectTriangle(v1, v0, v3, particle),
                intersectTriangle(v2, v3, v0, particle),
                intersectTriangle(v3, v2, v1, particle)
            };

            auto chit = std::numeric_limits<double>::max();
            int normalIdx = 0;
            for (int i = 0; i < 4; ++i) {
                if (auto h = *(hits[i]); hits[i]) {
                    if (h > 0) {
                        if (h < chit) {
                            chit = h;
                            normalIdx = i;
                        }
                    }
                }
            }

            VisualizationIntersectionResult<U> res;
            if (chit < std::numeric_limits<double>::max()) {
                res.rayOriginIsInsideItem = pointInside(v0, v1, v2, v3, particle.pos);
                res.intersection = chit;
                res.intersectionValid = true;
                if (normalIdx == 0)
                    res.normal = normalVector<true>(v0, v1, v2);
                else if (normalIdx == 1)
                    res.normal = normalVector<true>(v1, v0, v3);
                else if (normalIdx == 2)
                    res.normal = normalVector<true>(v2, v3, v0);
                else
                    res.normal = normalVector<true>(v3, v2, v1);
            }
            return res;
        }

    } // namespace tetrahedron
}
}
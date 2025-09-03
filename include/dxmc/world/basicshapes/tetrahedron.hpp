/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <array>
#include <numeric>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace tetrahedron {

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

        static constexpr std::array<double, 6> AABB(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3)
        {
            std::array<double, 6> aabb;
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min({ v0[i], v1[i], v2[i], v3[i] });
                aabb[i + 3] = std::max({ v0[i], v1[i], v2[i], v3[i] });
            }
            return aabb;
        }

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

        static double volume(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3)
        {
            const auto a = vectormath::subtract(v1, v0);
            const auto b = vectormath::subtract(v2, v0);
            const auto c = vectormath::subtract(v3, v0);
            static constexpr double scale = 1.0 / 6.0;
            return scale * std::abs(vectormath::tripleProduct(a, b, c));
        }

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

        static std::optional<double> intersectTriangle(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const ParticleType auto& p)
        {
            const auto E1 = vectormath::subtract(v1, v0);
            const auto TT = vectormath::subtract(p.pos, v0);
            const auto Q = vectormath::cross(TT, E1);

            const auto E2 = vectormath::subtract(v2, v0);
            const auto P = vectormath::cross(p.dir, E2);

            const auto det_inv = 1 / vectormath::dot(P, E1);

            const auto v = vectormath::dot(Q, p.dir) * det_inv;
            const auto u = vectormath::dot(P, TT) * det_inv;

            return v >= 0 && u >= 0 && (u + v) <= 1 ? std::make_optional(vectormath::dot(Q, E2) * det_inv)
                                                    : std::nullopt;
        }

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
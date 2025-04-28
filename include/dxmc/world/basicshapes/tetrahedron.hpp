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

        static constexpr bool pointInside(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3, const std::array<double, 3>& point)
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

        static constexpr double volume(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3)
        {
            const auto a = vectormath::subtract(v1, v0);
            const auto b = vectormath::subtract(v2, v0);
            const auto c = vectormath::subtract(v3, v0);
            static constexpr double scale = 1.0 / 6.0;
            return scale * std::abs(vectormath::tripleProduct(a, b, c));
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
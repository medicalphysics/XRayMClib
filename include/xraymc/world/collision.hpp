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

#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

namespace xraymc {
namespace collision {

    inline bool testCollision(const std::array<double, 6>& AABBa, const std::array<double, 6>& AABBb)
    {
        // AABB layout: {min_x, min_y, min_z, max_x, max_y, max_z}
        // Overlap requires overlap on all three axes.
        for (std::size_t i = 0; i < 3; ++i) {
            if (AABBa[i] > AABBb[i + 3] || AABBb[i] > AABBa[i + 3])
                return false;
        }
        return true;
    }

    namespace detail {

        inline bool coplanarTriTri(const std::array<double, 3>& N, const std::array<std::array<double, 3>, 3>& t1, const std::array<std::array<double, 3>, 3>& t2)
        {
            const std::array<double, 3> A { std::abs(N[0]), std::abs(N[1]), std::abs(N[2]) };
            std::size_t i0 = 0;
            std::size_t i1 = 1;
            if (A[0] > A[1]) {
                if (A[0] > A[2]) {
                    i0 = 1;
                    i1 = 2;
                } else {
                    i0 = 0;
                    i1 = 1;
                }
            } else {
                if (A[2] > A[1]) {
                    i0 = 0;
                    i1 = 1;
                } else {
                    i0 = 0;
                    i1 = 2;
                }
            }

            const auto edgeEdge = [&](const std::array<double, 3>& v0, const std::array<double, 3>& v1,
                                      const std::array<double, 3>& u0, const std::array<double, 3>& u1) {
                const auto ax = v1[i0] - v0[i0];
                const auto ay = v1[i1] - v0[i1];
                const auto bx = u0[i0] - u1[i0];
                const auto by = u0[i1] - u1[i1];
                const auto cx = v0[i0] - u0[i0];
                const auto cy = v0[i1] - u0[i1];
                const auto f = ay * bx - ax * by;
                const auto d = by * cx - bx * cy;
                if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
                    const auto e = ax * cy - ay * cx;
                    if (f > 0)
                        return e >= 0 && e <= f;
                    return e <= 0 && e >= f;
                }
                return false;
            };

            const auto edgeAgainstTri = [&](const std::array<double, 3>& v0, const std::array<double, 3>& v1,
                                            const std::array<std::array<double, 3>, 3>& u) {
                return edgeEdge(v0, v1, u[0], u[1]) || edgeEdge(v0, v1, u[1], u[2]) || edgeEdge(v0, v1, u[2], u[0]);
            };

            const auto pointInTri = [&](const std::array<double, 3>& p,
                                        const std::array<std::array<double, 3>, 3>& u) {
                auto a = u[1][i1] - u[0][i1];
                auto b = -(u[1][i0] - u[0][i0]);
                auto c = -a * u[0][i0] - b * u[0][i1];
                const auto d0 = a * p[i0] + b * p[i1] + c;

                a = u[2][i1] - u[1][i1];
                b = -(u[2][i0] - u[1][i0]);
                c = -a * u[1][i0] - b * u[1][i1];
                const auto d1 = a * p[i0] + b * p[i1] + c;

                a = u[0][i1] - u[2][i1];
                b = -(u[0][i0] - u[2][i0]);
                c = -a * u[2][i0] - b * u[2][i1];
                const auto d2 = a * p[i0] + b * p[i1] + c;

                return d0 * d1 > 0 && d0 * d2 > 0;
            };

            if (edgeAgainstTri(t1[0], t1[1], t2) || edgeAgainstTri(t1[1], t1[2], t2) || edgeAgainstTri(t1[2], t1[0], t2))
                return true;
            return pointInTri(t1[0], t2) || pointInTri(t2[0], t1);
        }

        // Computes the interval [isect0, isect1] carved out of the plane-plane
        // intersection line by a triangle, following Möller's "A Fast Triangle-Triangle
        // Intersection Test". `proj` holds the triangle vertices projected onto the
        // dominant axis of the intersection line, `dist` the signed distances of the
        // vertices to the other triangle's plane. Returns std::nullopt when the two
        // triangles are coplanar (all distances zero), in which case the caller must
        // fall back to the coplanar test.
        inline std::optional<std::array<double, 2>> triangleLineInterval(
            const std::array<double, 3>& proj, const std::array<double, 3>& dist,
            double d0d1, double d0d2)
        {
            std::array<double, 2> isect { };
            const auto isectFn = [&](double p0, double p1, double p2, double a0, double a1, double a2) {
                const auto t0 = a0 / (a0 - a1);
                const auto t1 = a0 / (a0 - a2);
                isect[0] = p0 + (p1 - p0) * t0;
                isect[1] = p0 + (p2 - p0) * t1;
            };

            if (d0d1 > 0.0) // dist[2] is the vertex on its own side of the plane
                isectFn(proj[2], proj[0], proj[1], dist[2], dist[0], dist[1]);
            else if (d0d2 > 0.0)
                isectFn(proj[1], proj[0], proj[2], dist[1], dist[0], dist[2]);
            else if (dist[1] * dist[2] > 0.0 || dist[0] != 0.0)
                isectFn(proj[0], proj[1], proj[2], dist[0], dist[1], dist[2]);
            else if (dist[1] != 0.0)
                isectFn(proj[1], proj[0], proj[2], dist[1], dist[0], dist[2]);
            else if (dist[2] != 0.0)
                isectFn(proj[2], proj[0], proj[1], dist[2], dist[0], dist[1]);
            else
                return std::nullopt; // coplanar

            if (isect[0] > isect[1])
                std::swap(isect[0], isect[1]);
            return isect;
        }

        /// @brief Wraps raw vertex triples (e.g. a tetrahedral-mesh outer contour) as Triangle objects.
        inline std::vector<Triangle> toTriangles(const std::vector<std::array<std::array<double, 3>, 3>>& verts)
        {
            std::vector<Triangle> out;
            out.reserve(verts.size());
            for (const auto& v : verts)
                out.emplace_back(v);
            return out;
        }

        // Akenine-Möller plane/box overlap: does the plane through `vert` with the
        // given `normal` cross the box centered at the origin with half-extents
        // `half`? `vert` is given relative to the box center.
        inline bool planeBoxOverlap(const std::array<double, 3>& normal, const std::array<double, 3>& vert, const std::array<double, 3>& half)
        {
            std::array<double, 3> vmin {};
            std::array<double, 3> vmax {};
            for (std::size_t q = 0; q < 3; ++q) {
                if (normal[q] > 0.0) {
                    vmin[q] = -half[q] - vert[q];
                    vmax[q] = half[q] - vert[q];
                } else {
                    vmin[q] = half[q] - vert[q];
                    vmax[q] = -half[q] - vert[q];
                }
            }
            if (vectormath::dot(normal, vmin) > 0.0)
                return false;
            return vectormath::dot(normal, vmax) >= 0.0;
        }

    } // namespace detail

    inline bool testCollision(const Triangle& a, const Triangle& b)
    {
        const auto& Va = a.vertices();
        const auto& Vb = b.vertices();

        // Plane of triangle a
        const auto Na = vectormath::cross(vectormath::subtract(Va[1], Va[0]), vectormath::subtract(Va[2], Va[0]));
        const auto da = -vectormath::dot(Na, Va[0]);

        // Signed distances of b's vertices to a's plane
        std::array<double, 3> db = {
            vectormath::dot(Na, Vb[0]) + da,
            vectormath::dot(Na, Vb[1]) + da,
            vectormath::dot(Na, Vb[2]) + da
        };
        for (auto& v : db)
            if (std::abs(v) < GEOMETRIC_ERROR<>())
                v = 0.0;
        const auto db0db1 = db[0] * db[1];
        const auto db0db2 = db[0] * db[2];
        if (db0db1 > 0.0 && db0db2 > 0.0)
            return false; // b lies entirely on one side of a's plane

        // Plane of triangle b
        const auto Nb = vectormath::cross(vectormath::subtract(Vb[1], Vb[0]), vectormath::subtract(Vb[2], Vb[0]));
        const auto dbp = -vectormath::dot(Nb, Vb[0]);

        // Signed distances of a's vertices to b's plane
        std::array<double, 3> dva = {
            vectormath::dot(Nb, Va[0]) + dbp,
            vectormath::dot(Nb, Va[1]) + dbp,
            vectormath::dot(Nb, Va[2]) + dbp
        };
        for (auto& v : dva)
            if (std::abs(v) < GEOMETRIC_ERROR<>())
                v = 0.0;
        const auto dva0dva1 = dva[0] * dva[1];
        const auto dva0dva2 = dva[0] * dva[2];
        if (dva0dva1 > 0.0 && dva0dva2 > 0.0)
            return false; // a lies entirely on one side of b's plane

        // Direction of the line where the two planes meet; project onto its
        // dominant axis so the 3-D problem collapses to a 1-D interval overlap.
        const auto D = vectormath::cross(Na, Nb);
        std::size_t index = 0;
        const std::array<double, 3> absD { std::abs(D[0]), std::abs(D[1]), std::abs(D[2]) };
        if (absD[1] > absD[index])
            index = 1;
        if (absD[2] > absD[index])
            index = 2;

        const std::array<double, 3> projA { Va[0][index], Va[1][index], Va[2][index] };
        const std::array<double, 3> projB { Vb[0][index], Vb[1][index], Vb[2][index] };

        const auto isectA = detail::triangleLineInterval(projA, dva, dva0dva1, dva0dva2);
        const auto isectB = detail::triangleLineInterval(projB, db, db0db1, db0db2);

        if (!isectA || !isectB) // coplanar triangles
            return detail::coplanarTriTri(Na, Va, Vb);

        // The triangles intersect iff their intervals on the plane-plane line overlap.
        return !((*isectA)[1] < (*isectB)[0] || (*isectB)[1] < (*isectA)[0]);
    }

    // Triangle vs axis-aligned box (AABB layout {min_x, min_y, min_z, max_x, max_y, max_z}),
    // via Akenine-Möller's "Fast 3D Triangle-Box Overlap Test" (separating-axis theorem).
    inline bool testCollision(const Triangle& a, const std::array<double, 6>& AABB)
    {
        const std::array<double, 3> center {
            (AABB[0] + AABB[3]) * 0.5,
            (AABB[1] + AABB[4]) * 0.5,
            (AABB[2] + AABB[5]) * 0.5
        };
        const std::array<double, 3> half {
            (AABB[3] - AABB[0]) * 0.5,
            (AABB[4] - AABB[1]) * 0.5,
            (AABB[5] - AABB[2]) * 0.5
        };

        // Triangle vertices expressed relative to the box center.
        const auto& V = a.vertices();
        const std::array<std::array<double, 3>, 3> v {
            vectormath::subtract(V[0], center),
            vectormath::subtract(V[1], center),
            vectormath::subtract(V[2], center)
        };

        // Bullet 1: the triangle's own AABB against the box (3 axis tests).
        for (std::size_t i = 0; i < 3; ++i) {
            const auto mn = std::min({ v[0][i], v[1][i], v[2][i] });
            const auto mx = std::max({ v[0][i], v[1][i], v[2][i] });
            if (mn > half[i] || mx < -half[i])
                return false;
        }

        const std::array<std::array<double, 3>, 3> edges {
            vectormath::subtract(v[1], v[0]),
            vectormath::subtract(v[2], v[1]),
            vectormath::subtract(v[0], v[2])
        };

        // Bullet 3: 9 axis tests, axis = boxAxis_k x edge_j.
        const auto axisTest = [&](const std::array<double, 3>& axis) {
            auto mn = vectormath::dot(axis, v[0]);
            auto mx = mn;
            for (std::size_t k = 1; k < 3; ++k) {
                const auto d = vectormath::dot(axis, v[k]);
                mn = std::min(mn, d);
                mx = std::max(mx, d);
            }
            const auto rad = std::abs(axis[0]) * half[0] + std::abs(axis[1]) * half[1] + std::abs(axis[2]) * half[2];
            return !(mn > rad || mx < -rad);
        };

        for (const auto& e : edges) {
            if (!axisTest({ 0.0, -e[2], e[1] })) // boxAxis (1,0,0) x e
                return false;
            if (!axisTest({ e[2], 0.0, -e[0] })) // boxAxis (0,1,0) x e
                return false;
            if (!axisTest({ -e[1], e[0], 0.0 })) // boxAxis (0,0,1) x e
                return false;
        }

        // Bullet 2: the triangle's plane against the box.
        const auto normal = vectormath::cross(edges[0], edges[1]);
        return detail::planeBoxOverlap(normal, v[0], half);
    }

    inline bool testCollision(const std::array<double, 6>& AABB, const Triangle& a)
    {
        return testCollision(a, AABB);
    }

    inline bool testCollision(const std::vector<Triangle>& trianglesA, const std::vector<Triangle>& trianglesB)
    {
        for (const auto& ta : trianglesA) {
            for (const auto& tb : trianglesB) {
                if (testCollision(ta, tb))
                    return true;
            }
        }
        return false;
    }

    template <int NMaterialShellsA, int LOWENERGYCORRECTIONA, int NMaterialShellsB, int LOWENERGYCORRECTIONB>
    inline bool testCollision(const TriangulatedMesh<NMaterialShellsA, LOWENERGYCORRECTIONA>& A, const TriangulatedMesh<NMaterialShellsB, LOWENERGYCORRECTIONB>& B)
    {

        if (!testCollision(A.AABB(), B.AABB()))
            return false;

        const auto& trianglesA = A.triangles();
        const auto& trianglesB = B.triangles();

        return testCollision(trianglesA, trianglesB);
    }

    template <int NMaterialShellsA, int LOWENERGYCORRECTIONA, bool FORCEDINTERACTIONA, int NMaterialShellsB, int LOWENERGYCORRECTIONB, bool FORCEDINTERACTIONB>
    inline bool testCollision(const TetrahedalMesh<NMaterialShellsA, LOWENERGYCORRECTIONA, FORCEDINTERACTIONA>& A, const TetrahedalMesh<NMaterialShellsB, LOWENERGYCORRECTIONB, FORCEDINTERACTIONB>& B)
    {
        if (!testCollision(A.AABB(), B.AABB()))
            return false;

        const auto triA = detail::toTriangles(A.constructOuterContourTriangles());
        const auto triB = detail::toTriangles(B.constructOuterContourTriangles());

        return testCollision(triA, triB);
    }

    template <int NMaterialShellsA, int LOWENERGYCORRECTIONA, bool FORCEDINTERACTIONA, int NMaterialShellsB, int LOWENERGYCORRECTIONB>
    inline bool testCollision(const TetrahedalMesh<NMaterialShellsA, LOWENERGYCORRECTIONA, FORCEDINTERACTIONA>& A, const TriangulatedMesh<NMaterialShellsB, LOWENERGYCORRECTIONB>& B)
    {
        if (!testCollision(A.AABB(), B.AABB()))
            return false;

        const auto triA = detail::toTriangles(A.constructOuterContourTriangles());

        return testCollision(triA, B.triangles());
    }

    template <int NMaterialShellsA, int LOWENERGYCORRECTIONA, bool FORCEDINTERACTIONA, int NMaterialShellsB, int LOWENERGYCORRECTIONB>
    inline bool testCollision(const TriangulatedMesh<NMaterialShellsB, LOWENERGYCORRECTIONB>& A, const TetrahedalMesh<NMaterialShellsA, LOWENERGYCORRECTIONA, FORCEDINTERACTIONA>& B)
    {
        return testCollision(B, A);
    }
}
} // namespace xraymc

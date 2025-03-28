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

Copyright 2025 Erlend Andersen
*/

#pragma once
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"

namespace dxmc {
namespace collision {

    bool collideBasicShape(const dxmc::Triangle& t1, const dxmc::Triangle& t2)
    {
        // Thomas Muller fast triangle intersection test
        const auto& v1 = t1.vertices();
        const auto& v2 = t2.vertices();

        const auto e20 = vectormath::subtract(v2[1], v2[0]);
        const auto e21 = vectormath::subtract(v2[2], v2[0]);
        const auto n2 = vectormath::cross(e20, e21);
        const auto d2 = -vectormath::dot(n2, v2[0]);

        int side = 0;
        for (int i = 0; i < 3; ++i) {
            const auto dv = vectormath::dot(n2, v1[i]) + d2;
            if (dv > 0)
                side++;
            else
                side--;
        }
        if (side == 3 || side == -3)
            return false;

        const auto e10 = vectormath::subtract(v1[1], v1[0]);
        const auto e11 = vectormath::subtract(v1[2], v1[0]);
        const auto n1 = vectormath::cross(e10, e11);
        const auto d1 = -vectormath::dot(n1, v1[0]);

        side = 0;
        for (int i = 0; i < 3; ++i) {
            const auto dv = vectormath::dot(n1, v2[i]) + d1;
            if (dv > GEOMETRIC_ERROR<>())
                side++;
            else
                side--;
        }
        return side != 3 && side != -3;
    }

    bool collideBasicShape(const dxmc::Triangle& tri, const dxmc::Tetrahedron& tet)
    {
        const auto& vertices = tet.vertices();
        const auto t1 = dxmc::Triangle(vertices[0], vertices[1], vertices[2]);
        const auto t2 = dxmc::Triangle(vertices[1], vertices[0], vertices[3]);
        const auto t3 = dxmc::Triangle(vertices[2], vertices[3], vertices[0]);
        const auto t4 = dxmc::Triangle(vertices[3], vertices[2], vertices[1]);

        return collideBasicShape(tri, t1) || collideBasicShape(tri, t2) || collideBasicShape(tri, t3) || collideBasicShape(tri, t4) || tet.pointInside(tri.vertices()[0]);
    }

    bool collideBasicShape(const dxmc::Tetrahedron& tet, const dxmc::Triangle& tri)
    {
        return collideBasicShape(tri, tet);
    }

    bool collideBasicShape(const dxmc::Tetrahedron& tet1, const dxmc::Tetrahedron& tet2)
    {
        const auto& vertices = tet1.vertices();
        const auto t1 = dxmc::Triangle(vertices[0], vertices[1], vertices[2]);
        const auto t2 = dxmc::Triangle(vertices[1], vertices[0], vertices[3]);
        const auto t3 = dxmc::Triangle(vertices[2], vertices[3], vertices[0]);
        const auto t4 = dxmc::Triangle(vertices[3], vertices[2], vertices[1]);
        return collideBasicShape(tet2, t1) || collideBasicShape(tet2, t2) || collideBasicShape(tet2, t3) || collideBasicShape(tet2, t4);
    }

    bool collideBasicShape(const std::vector<dxmc::Triangle>& tri1, const std::vector<dxmc::Triangle>& tri2)
    {
        auto it = std::find_first_of(std::execution::par_unseq, tri1.cbegin(), tri1.cend(), tri2.cbegin(), tri2.cend(), [](const auto& t1, const auto& t2) {
            return collideBasicShape(t1, t2);
        });
        return it != tri1.cend();
    }

    bool collideBasicShape(const std::vector<dxmc::Triangle>& tri1, const std::vector<dxmc::Tetrahedron>& tet2)
    {
        auto it = std::find_first_of(std::execution::par_unseq, tri1.cbegin(), tri1.cend(), tet2.cbegin(), tet2.cend(), [](const auto& t1, const auto& t2) {
            return collideBasicShape(t1, t2);
        });
        return it != tri1.cend();
    }

    bool collideBasicShape(const std::vector<dxmc::Tetrahedron>& tet2, const std::vector<dxmc::Triangle>& tri1)
    {
        return collideBasicShape(tri1, tet2);
    }

    bool collideBasicShape(const std::vector<dxmc::Tetrahedron>& tet1, const std::vector<dxmc::Tetrahedron>& tet2)
    {
        auto it = std::find_first_of(std::execution::par_unseq, tet1.cbegin(), tet1.cend(), tet2.cbegin(), tet2.cend(), [](const auto& t1, const auto& t2) {
            return collideBasicShape(t1, t2);
        });
        return it != tet1.cend();
    }

    template <typename T>
    concept TriangleMeshType = requires(T a) {
        { a.triangles() } -> std::convertible_to<std::vector<dxmc::Triangle>>;
    };

    template <typename T>
    concept TetrahedronMeshType = requires(T a) {
        { a.tetrahedrons() } -> std::convertible_to<std::vector<dxmc::Tetrahedron>>;
    };

    template <typename T1, typename T2>
    bool collide(const T1& m1, const T2& m2)
    {
        if constexpr (TriangleMeshType<T1> && TriangleMeshType<T2>) {
            const auto& tri1 = m1.triangles();
            const auto& tri2 = m2.triangles();
            return collideBasicShape(tri1, tri2);
        } else if constexpr (TriangleMeshType<T1> && TetrahedronMeshType<T2>) {
            const auto& tri1 = m1.triangles();
            const auto& tet2 = m2.tetrahedrons();
            return collideBasicShape(tri1, tet2);
        } else if constexpr (TriangleMeshType<T2> && TetrahedronMeshType<T1>) {
            const auto& tri1 = m2.triangles();
            const auto& tet2 = m1.tetrahedrons();
            return collideBasicShape(tri1, tet2);
        } else if constexpr (TetrahedronMeshType<T1> && TetrahedronMeshType<T2>) {
            const auto& tet1 = m1.tetrahedrons();
            const auto& tet2 = m2.tetrahedrons();
            return collideBasicShape(tet1, tet2);
        }
        return true;
    }

}
}

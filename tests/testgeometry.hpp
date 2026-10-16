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

// Simple test geometries shared by several test executables.

#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangle.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace xraymctest {

// Closed triangulated cube [-scale, scale]^3 with outward facing normals
inline std::vector<xraymc::Triangle> boxTriangles(double scale = 1)
{
    std::vector<std::array<double, 3>> p;
    p.push_back({ 1, 1, 1 }); // 0
    p.push_back({ 1, 1, -1 }); // 1
    p.push_back({ 1, -1, 1 }); // 2
    p.push_back({ -1, 1, 1 }); // 3
    p.push_back({ -1, -1, 1 }); // 4
    p.push_back({ -1, 1, -1 }); // 5
    p.push_back({ 1, -1, -1 }); // 6
    p.push_back({ -1, -1, -1 }); // 7
    for (auto& i : p)
        for (auto& j : i)
            j *= scale;

    std::vector<xraymc::Triangle> t;
    t.push_back({ p[0], p[3], p[4] });
    t.push_back({ p[0], p[4], p[2] });
    t.push_back({ p[6], p[2], p[4] });
    t.push_back({ p[6], p[4], p[7] });
    t.push_back({ p[7], p[4], p[3] });
    t.push_back({ p[7], p[3], p[5] });
    t.push_back({ p[5], p[1], p[6] });
    t.push_back({ p[5], p[6], p[7] });
    t.push_back({ p[1], p[0], p[2] });
    t.push_back({ p[1], p[2], p[6] });
    t.push_back({ p[5], p[3], p[0] });
    t.push_back({ p[5], p[0], p[1] });
    return t;
}

// Closed square pyramid with base [-d, d]^2 in z = 0 and apex at (0, 0, d)
inline std::vector<xraymc::Triangle> pyramidTriangles(double d = 30)
{
    std::vector<std::array<double, 3>> p;
    p.push_back({ 1, 1, 0 }); // 0
    p.push_back({ 1, -1, 0 }); // 1
    p.push_back({ -1, -1, 0 }); // 2
    p.push_back({ -1, 1, 0 }); // 3
    p.push_back({ 0, 0, 1 }); // 4
    for (auto& i : p)
        for (auto& j : i)
            j *= d;

    std::vector<xraymc::Triangle> t;
    t.push_back({ p[0], p[1], p[4] });
    t.push_back({ p[1], p[2], p[4] });
    t.push_back({ p[2], p[3], p[4] });
    t.push_back({ p[3], p[0], p[4] });
    // underside
    t.push_back({ p[0], p[3], p[2] });
    t.push_back({ p[2], p[1], p[0] });
    return t;
}

// Square [-scale, scale]^2 in the z = 0 plane made of two triangles
inline std::vector<xraymc::Triangle> planeTriangles(double scale = 1)
{
    const std::array<double, 3> a = { -scale, -scale, 0 };
    const std::array<double, 3> b = { scale, -scale, 0 };
    const std::array<double, 3> c = { scale, scale, 0 };
    const std::array<double, 3> d = { -scale, scale, 0 };
    return { { a, b, c }, { a, c, d } };
}

// Cube [-1, 1]^3 split into 12 tetrahedra (two per face, meeting at the two
// face centers on z = +/-1), filled with water of density 1.
inline xraymc::TetrahedalMeshData tetCube(double scale = 1)
{
    std::vector<std::array<double, 3>> v(10);
    v[0] = { 1, -1, 1 };
    v[1] = { -1, 1, -1 };
    v[2] = { -1, -1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, 1, 1 };
    v[6] = { 1, 1, 1 };
    v[7] = { 1, 1, -1 };
    v[8] = { 0, 0, 1 };
    v[9] = { 0, 0, -1 };
    for (auto& i : v)
        for (auto& j : i)
            j *= scale;

    std::vector<std::array<std::uint32_t, 4>> t(12);
    t[0] = { 0, 7, 4, 9 };
    t[1] = { 1, 6, 5, 8 };
    t[2] = { 1, 6, 8, 9 };
    t[3] = { 6, 1, 7, 9 };
    t[4] = { 1, 5, 2, 8 };
    t[5] = { 3, 0, 4, 9 };
    t[6] = { 3, 1, 2, 9 };
    t[7] = { 2, 1, 8, 9 };
    t[8] = { 7, 0, 6, 9 };
    t[9] = { 6, 0, 8, 9 };
    t[10] = { 0, 3, 2, 9 };
    t[11] = { 0, 2, 8, 9 };

    xraymc::TetrahedalMeshData data;
    data.nodes = v;
    data.elements = t;
    data.collectionIndices.resize(t.size(), 0);
    data.collectionMaterialComposition.resize(1);
    data.collectionMaterialComposition[0][1] = 0.111894;
    data.collectionMaterialComposition[0][8] = 0.888106;
    data.collectionDensities.resize(1, 1.0);
    data.makeGenericCollectionNames();
    return data;
}

}

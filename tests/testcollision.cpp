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

Copyright 2025 Erlend Andersen
*/

#include "xraymc/world/collisions/collisiontest.hpp"

#include <iostream>

std::vector<xraymc::Tetrahedron> tetCube()
{
    std::vector<std::array<double, 3>> v(16);
    v[0] = { -1, 1, 1 };
    v[1] = { 1, -1, 1 };
    v[2] = { 1, 1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, -1, 1 };
    v[6] = { -1, 1, -1 };
    v[7] = { 1, 1, -1 };

    v[8] = { -1, 1, 3 };
    v[9] = { 1, -1, 3 };
    v[10] = { 1, 1, 3 };
    v[11] = { -1, -1, 1 };
    v[12] = { 1, -1, 1 };
    v[13] = { -1, -1, 3 };
    v[14] = { -1, 1, 1 };
    v[15] = { 1, 1, 1 };

    for (auto& i : v)
        for (auto& n : i)
            n *= 1;

    std::vector<xraymc::Tetrahedron> t(12);

    t[0] = { v[1], v[7], v[0], v[2], 0, 0 }; //*
    t[1] = { v[7], v[3], v[0], v[6], 0, 0 }; //*
    t[2] = { v[1], v[3], v[0], v[4], 0, 0 }; //*
    t[3] = { v[1], v[7], v[4], v[0], 0, 0 }; //*
    t[4] = { v[7], v[3], v[4], v[0], 0, 0 }; //*
    t[5] = { v[1], v[3], v[5], v[0], 0, 0 }; //*

    t[6] = { v[1 + 8], v[7 + 8], v[0 + 8], v[2 + 8], 0, 0 }; //*
    t[7] = { v[7 + 8], v[3 + 8], v[0 + 8], v[6 + 8], 0, 0 }; //*
    t[8] = { v[1 + 8], v[3 + 8], v[0 + 8], v[4 + 8], 0, 0 }; //*
    t[9] = { v[1 + 8], v[7 + 8], v[4 + 8], v[0 + 8], 0, 0 }; //*
    t[10] = { v[7 + 8], v[3 + 8], v[4 + 8], v[0 + 8], 0, 0 }; //*
    t[11] = { v[1 + 8], v[3 + 8], v[5 + 8], v[0 + 8], 0, 0 }; //*

    bool success = true;
    for (auto& tet : t)
        success = success && tet.validVerticeOrientation();

    if (!success)
        t.clear();
    return t;
}
template <std::size_t N = 5, int L = 2, bool Fluence = true>
xraymc::TetrahedalMesh<N, L, Fluence> cubetetrahedron()
{
    auto tets = tetCube();

    std::vector<xraymc::Material<N>> mats;
    mats.push_back(xraymc::Material<N>::byNistName("Water, Liquid").value());
    std::vector<double> dens(tets.size(), 1);

    std::vector<std::string> names(1);
    xraymc::TetrahedalMesh<N, L, Fluence> mesh(tets, dens, mats, names, 8);

    return mesh;
}

bool testTetCube()
{
    auto c1 = cubetetrahedron<5, 1, false>();
    auto c2 = cubetetrahedron<5, 1, false>();
    c1.translate({ 0, 0, 6 });
    auto a1 = c1.AABB();
    auto a2 = c2.AABB();

    auto hit = xraymc::CollisionTests::test(c1, c2);
    return hit;
}

int main()
{
    std::cout << "Testing collision of shapes\n";

    bool success = true;
    success = success && testTetCube();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
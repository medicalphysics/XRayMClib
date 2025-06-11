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

#include "dxmc/world/worlditems/tetrahedalmesh2/tetrahedalmeshdata.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh2/tetrahedalmeshgrid2.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh2/tetrahedalmeshreader2.hpp"

#include <iostream>
#include <string>
#include <vector>

using index_t = std::uint32_t;

dxmc::TetrahedalMeshData tetCube()
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

    std::vector<std::array<index_t, 4>> t(12);

    t[0] = { 1, 7, 0, 2 }; //*
    t[1] = { 7, 3, 0, 6 }; //*
    t[2] = { 1, 3, 0, 4 }; //*
    t[3] = { 1, 7, 4, 0 }; //*
    t[4] = { 7, 3, 4, 0 }; //*
    t[5] = { 1, 3, 5, 0 }; //*

    t[6] = { 1 + 8, 7 + 8, 0 + 8, 2 + 8 }; //*
    t[7] = { 7 + 8, 3 + 8, 0 + 8, 6 + 8 }; //*
    t[8] = { 1 + 8, 3 + 8, 0 + 8, 4 + 8 }; //*
    t[9] = { 1 + 8, 7 + 8, 4 + 8, 0 + 8 }; //*
    t[10] = { 7 + 8, 3 + 8, 4 + 8, 0 + 8 }; //*
    t[11] = { 1 + 8, 3 + 8, 5 + 8, 0 + 8 }; //*

    dxmc::TetrahedalMeshData data;
    data.nodes = v;
    data.elements = t;
    data.collectionIndices.resize(t.size(), 0);
    return data;
}

bool testTetrahedalMeshGrid()
{

    auto data = tetCube();

    // dxmc::TetrahedalMeshGrid2 grid(data);
    return false;
}

int main()
{

    std::string material_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF_media.dat";
    std::string organ_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\icrp145organs.csv";
    std::string node_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.node";
    std::string element_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.ele";

    // dxmc::TetrahedalmeshReader2 testreader(node_file, element_file, material_file, organ_file);

    node_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.node";
    element_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.ele";
    dxmc::TetrahedalmeshReader2 testreader(node_file, element_file);
    dxmc::TetrahedalMeshGrid2 grid(testreader.data(), { 32, 32, 16 });

    auto test = grid.pointInside({ -1, 0, 0 });

    dxmc::Particle p { .pos = { 0, -100, 0 }, .dir = { 0, 1, 0 } };
    auto res = grid.intersect(p);

    bool success = true;

    // auto data = tetCube();
    success = success && testTetrahedalMeshGrid();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
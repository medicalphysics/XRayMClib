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

Copyright 2022 Erlend Andersen
*/

#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh2.hpp"

#include <iostream>
#include <string>
#include <vector>

void test()
{

    dxmc::TetrahedalmeshReader reader;

    std::string node = R"(C:\Users\ander\OneDrive\phantomsMNCP\Pregnant_MRCPs\Pregnant_MRCPs\2. TM-version MRCPs\38wM.node)";
    std::string ele =  R"(C:\Users\ander\OneDrive\phantomsMNCP\Pregnant_MRCPs\Pregnant_MRCPs\2. TM-version MRCPs\38wM.ele)";
    std::string matorg = R"(C:\Users\ander\OneDrive\phantomsMNCP\Pregnant_MRCPs\Pregnant_MRCPs\3. Material Files (for TM-version MRCPs)\38wM.material)";

    reader.readICRPPregnantPhantom(node, ele, matorg);
}

bool testtetmesh2()
{
    using TetMesh = dxmc::TetrahedalMesh2<5, 1, false>;
    using Data = dxmc::TetrahedalMeshData;
    

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

    std::vector<std::array<std::size_t, 4>> t(12);
    t[0] = { 1, 7, 0, 2 }; //*
    t[1] = { 7, 3, 0, 6 }; //*
    t[2] = { 1, 3, 0, 4 }; //*
    t[3] = { 1, 7, 4, 0 }; //*
    t[4] = { 1, 3, 4, 0 }; //*
    t[5] = { 1, 3, 5, 0 }; //*

    t[6] = { 1 + 8, 7 + 8, 0 + 8, 2 + 8 }; //*
    t[7] = { 7 + 8, 3 + 8, 0 + 8, 6 + 8 }; //*
    t[8] = { 1 + 8, 3 + 8, 0 + 8, 4 + 8 }; //*
    t[9] = { 1 + 8, 7 + 8, 4 + 8, 0 + 8 }; //*
    t[10] = { 7 + 8, 3 + 8, 4 + 8, 0 + 8 }; //*
    t[11] = { 1 + 8, 3 + 8, 5 + 8, 0 + 8 }; //*

    Data data;
    data.nodes = v;
    data.elements = t;
    data.collectionIndices.resize(t.size(), 0);
    data.collectionDensities.resize(t.size(), 1);
    data.collectionMaterialComposition.resize(1, { { { 6, 1 } } });
    data.collectionNames.resize(1, "Carbon");
    TetMesh mesh(data);

    
    return false;
}

int main()
{
    testtetmesh2();
    test();
    return EXIT_SUCCESS;
}
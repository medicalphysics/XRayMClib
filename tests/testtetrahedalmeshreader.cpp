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

int main()
{
    test();
    return EXIT_SUCCESS;
}
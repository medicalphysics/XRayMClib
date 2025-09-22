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

Copyright 2022 Erlend Andersen
*/
#include "epicsparser.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

int main()
{

    const bool always_rebuild = XRAYMCLIB_REBUILD_PHYSICSLISTS;

    const std::string eadl(XRAYMCLIB_EADLPATH);
    const std::string epdl(XRAYMCLIB_EPDLPATH);
    const std::string outpath(XRAYMCLIB_PHYSICSLISTSPATH);
    const std::string hartreefock("ComptonProfiles.dat");
    const std::string standardensities("standarddensities.csv");

    auto file_exists = std::filesystem::exists(outpath);
    if (!file_exists || always_rebuild) {
        EPICSparser parser(eadl);
        parser.read(epdl);
        parser.readHartreeFockProfiles(hartreefock);
        parser.readStandardDensities(standardensities);

        auto data = parser.serializeElements();

        std::ofstream of;
        of.open(outpath, std::ios::binary);
        of.write(data.data(), data.size());
        of.close();
        std::cout << "Done generating physicslists.bin\n";
    } else {
        std::cout << "Used previously generated physicslists.bin\n";
    }
    return EXIT_SUCCESS;
}

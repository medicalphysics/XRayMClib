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
#include "epicsparser.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

int main()
{
    const std::string eadl(DXMCLIB_EADLPATH);
    const std::string epdl(DXMCLIB_EPDLPATH);
    const std::string outpath(DXMCLIB_PHYSICSLISTSPATH);
    const std::string hartreefock("hartreeFockProfiles_0.csv");
    const std::string standardensities("standarddensities.csv");
    const std::string comptonprofiles("ComptonProfiles.dat");

    auto file_exists = std::filesystem::exists(outpath);
    if (!file_exists) {
        EPICSparser parser(eadl);
        parser.read(epdl);
        parser.readHartreeFockProfiles(hartreefock);
        parser.readComptonProfiles(comptonprofiles);
        parser.readStandardDensities(standardensities);

        auto data = parser.serializeElements();

        std::ofstream of;
        of.open(outpath, std::ios::binary);
        of.write(data.data(), data.size());
        of.close();
    }

    return EXIT_SUCCESS;
}

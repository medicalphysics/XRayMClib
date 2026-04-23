

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

#include "xraymc/world/worlditems/personaldosimeter.hpp"

#include <iostream>

bool test_angles()
{
    xraymc::PersonalDosimeter dos;
    xraymc::Particle p;
    for (std::size_t i = 1; i < 100; ++i) {

        double ang = -std::numbers::pi_v<double> / 2 + i * std::numbers::pi_v<double> / 100.0;
        p.dir = xraymc::vectormath::rotate({ 0, 0, -1 },  { 1, 0, 0 }, ang);

        auto w = dos.angularResponseWeight(p);
        std::cout << ang << "," << w << std::endl;
    }
    return false;
}

int main()
{
    bool success = test_angles();

    if (success)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}
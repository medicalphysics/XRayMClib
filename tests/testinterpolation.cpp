
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

#include "xraymc/interpolation.hpp"

#include <iostream>
#include <string>
#include <vector>

bool testAkimaSpline()
{

    std::vector<std::pair<double, double>> flat = {
        { 0, 1 },
        { 1, -1 }
    };

    xraymc::AkimaSpline sp;
    sp.setup(flat);

    std::array<double, 4> test = { 0, 0.25, 0.33, 0.5 };
    for (auto a : test) {
        std::cout << a << " " << sp(a) << std::endl;
    }

    // generic filter from a Siemens Definition Flash
    std::vector<std::pair<double, double>> flash = {
        { 0.166511074, 3.53208 },
        { 0.000000000, 13.9167 },
        { 0.041992107, 12.5868 },
        { 0.083836642, 9.41943 },
        { 0.246954945, 1.96665 },
        { 0.324269441, 1.27605 },
        { 0.390607044, 0.947716 }
    };

    std::cout << "Flash filter:" << std::endl;
    sp.setup(flash);
    for (int i = 0; i <= 20; ++i) {
        const auto a = 0.4 * i / 20.0;
        std::cout << a << " " << sp(a) << std::endl;
    }

    xraymc::AkimaSplineStatic<double, 5> sp2;
    sp2.setup(flat);
    return true;
}

int main()
{

    bool success = true;
    success = success && testAkimaSpline();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
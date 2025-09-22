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

Copyright 2023 Erlend Andersen
*/

#include "xraymc/beams/filters/ctaecfilter.hpp"

#include <iostream>

bool testAECCTFilter()
{

    xraymc::CTAECFilter f;

    std::vector<double> w = { 5, 6, 7, 6, 3, 5, 6, 2, 6, 5, 6, 7, 8, 9, 9, 9 };

    double start = -5;
    double stop = 4;

    f.setData({ 0, 0, start }, { 0, 0, stop }, w);

    auto area = f.integrate();
    auto kk = (stop - start) / area;
    return kk > 0.95 && kk < 1.05;
}

bool testAECCTFilterBetween()
{

    xraymc::CTAECFilter f;

    std::vector<double> w(100, 2.0); //{ 5, 6, 7, 6, 3, 5, 6, 2, 6, 5, 6, 7, 8, 9, 9, 9 };

    double start = -5;
    double stop = 4;

    f.setData({ 0, 0, start }, { 0, 0, stop }, w);

    f.normalizeBetween({ 0, 0, -1 }, { 0, 0, 1 });

    auto area = f.integrate({ 0, 0, -1 }, { 0, 0, 1 });
    auto kk = 2.0 / area;
    bool success = kk > 0.95 && kk < 1.05;

    double mean = 0;
    constexpr int N = 100;
    for (int i = 0; i < N; ++i) {
        std::array pos = { 0., 0., -1.0 + (2.0 * i) / N };
        mean += f(pos);
    }
    mean /= N;
    success = success && mean > 0.95 && mean < 1.05;
    return success;
}

int main()
{

    bool success = true;
    success = success && testAECCTFilter();
    success = success && testAECCTFilterBetween();

    std::cout << "Test CT AEC Filter ";
    if (success)
        std::cout << "SUCCESS\n";
    else
        std::cout << "FAILURE\n";

    if (success) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}
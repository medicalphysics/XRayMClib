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

#include "dxmc/beams/utilities/spheresamplingcircularfield.hpp"
#include "dxmc/beams/utilities/spheresamplingrectangularfield.hpp"

#include <fstream>
#include <iostream>

bool testBeamSampling()
{

    dxmc::SphereSamplingRectangularField sampler;
    constexpr auto alpha = 15 * dxmc::DEG_TO_RAD();
    const auto h = 180 * std::tan(alpha);
    const auto y_ang_min = std::atan((h - 39.0 / 2) / 180.0);
    const auto y_ang_max = std::atan((h + 39.0 / 2) / 180.0);
    const auto x_ang = std::atan(39.0 / (2 * 180.0));
    //beam.setCollimationHalfAngles(-x_ang, y_ang_min, x_ang, y_ang_max);
    sampler.setData(-x_ang, y_ang_min, x_ang, y_ang_max);

    constexpr std::size_t N = 1e5;

    dxmc::RandomState state;
    std::vector<std::array<double, 3>> samps(N);
    std::vector<double> x(N);
    std::vector<double> y(N);
    std::vector<double> z(N);

    const std::array xdir = { 1.0, 0.0, 0.0 };
    const std::array ydir = { 0.0, 1.0, 0.0 };
    const auto dir_0 = dxmc::vectormath::cross(xdir, ydir);

    std::ofstream out("out.txt");
    // out << "x, y\n";
    for (std::size_t i = 0; i < N; ++i) {
        const auto dir = sampler(state);

        x[i] = dxmc::vectormath::dot(xdir, dir);
        y[i] = dxmc::vectormath::dot(ydir, dir);
        z[i] = dxmc::vectormath::dot(dir_0, dir);
        out << x[i] << "," << y[i] << "," << z[i] << '\n';
    }

    return false;
}

int main()
{
    bool success = true;
    success = success && testBeamSampling();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
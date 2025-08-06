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

#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"

#include <iostream>
#include <string>
#include <vector>

bool testMaterials()
{
    for (std::size_t i = 1; i < 101; ++i) {
        auto mat_opt = dxmc::Material<>::byZ(i);
        if (mat_opt) {
            auto& mat = mat_opt.value();
        } else
            return false;
    }
    return true;
}

bool testNistMaterials()
{
    for (const auto& name : dxmc::NISTMaterials::listNames()) {
        auto mat_opt = dxmc::Material<>::byNistName(name);
        if (mat_opt) {
            auto& mat = mat_opt.value();
        } else
            return false;
    }
    return true;
}
std::pair<double, std::map<std::size_t, double>> TG195_breast_tissue()
{
    std::map<std::size_t, double> adipose_w;
    adipose_w[1] = 11.2;
    adipose_w[6] = 61.9;
    adipose_w[7] = 1.7;
    adipose_w[8] = 25.1;
    adipose_w[15] = 0.025;
    adipose_w[16] = 0.025;
    adipose_w[19] = 0.025;
    adipose_w[20] = 0.025;

    const double adipose_d = 0.93;

    std::map<std::size_t, double> gland_w;
    gland_w[1] = 10.2;
    gland_w[6] = 18.4;
    gland_w[7] = 3.2;
    gland_w[8] = 67.7;
    gland_w[15] = 0.125;
    gland_w[16] = 0.125;
    gland_w[19] = 0.125;
    gland_w[20] = 0.125;

    const double gland_d = 1.04;

    // weighted 20% gland 80% adipose
    std::map<std::size_t, double> w;
    for (const auto [Z, n] : adipose_w) {
        if (!w.contains(Z))
            w[Z] = 0;
        w[Z] += n * 0.8;
    }
    for (const auto [Z, n] : gland_w) {
        if (!w.contains(Z))
            w[Z] = 0;
        w[Z] += n * 0.2;
    }
    const auto d = adipose_d * 0.8 + gland_d * 0.2;

    return std::make_pair(d, w);
}

void breastAtt()
{
    const auto [dens, weights] = TG195_breast_tissue();
    double e = 1;
    const double step = .5;
    const double end = 120;

    const auto breast = dxmc::Material<12>::byWeight(weights).value();
    while (e <= end) {
        const auto att = breast.attenuationValues(e);
        std::cout << e << ',' << att.photoelectric << ',' << att.incoherent << ',' << att.coherent << '\n';
        e = e + step;
    }
}

int main(int argc, char* argv[])
{
    bool success = true;
    // auto lead = dxmc::Material<12>::byZ(82);
    breastAtt();
    /*std::cout << "Basic tests of materials, please run without --fast_math flags: ";

    auto stxt = [](bool v) -> std::string { return v ? " SUCCSESS " : " FAILED "; };


    success = success && testMaterials();
    success = success && testNistMaterials();
    std::cout << stxt(success) << std::endl;
*/
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

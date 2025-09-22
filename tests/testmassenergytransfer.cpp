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

#include "xraymc/material/material.hpp"

#include <iostream>
#include <string>
#include <vector>

struct NISTData {
    NISTData(double e, double a, double aen)
        : energy(e * 1000.0)
        , att(a)
        , atten(aen)
    {
    }
    double energy, att, atten;
};

std::vector<NISTData> nistAir()
{
    std::vector<NISTData> data = {
        { 1.00000E-03, 3.606E+03, 3.599E+03 },
        { 1.50000E-03, 1.191E+03, 1.188E+03 },
        { 2.00000E-03, 5.279E+02, 5.262E+02 },
        { 3.00000E-03, 1.625E+02, 1.614E+02 },
        { 3.20290E-03, 1.340E+02, 1.330E+02 },
        { 4.00000E-03, 7.788E+01, 7.636E+01 },
        { 5.00000E-03, 4.027E+01, 3.931E+01 },
        { 6.00000E-03, 2.341E+01, 2.270E+01 },
        { 8.00000E-03, 9.921E+00, 9.446E+00 },
        { 1.00000E-02, 5.120E+00, 4.742E+00 },
        { 1.50000E-02, 1.614E+00, 1.334E+00 },
        { 2.00000E-02, 7.779E-01, 5.389E-01 },
        { 3.00000E-02, 3.538E-01, 1.537E-01 },
        { 4.00000E-02, 2.485E-01, 6.833E-02 },
        { 5.00000E-02, 2.080E-01, 4.098E-02 },
        { 6.00000E-02, 1.875E-01, 3.041E-02 },
        { 8.00000E-02, 1.662E-01, 2.407E-02 },
        { 1.00000E-01, 1.541E-01, 2.325E-02 },
        { 1.50000E-01, 1.356E-01, 2.496E-02 },
        { 2.00000E-01, 1.233E-01, 2.672E-02 },
        { 3.00000E-01, 1.067E-01, 2.872E-02 },
        { 4.00000E-01, 9.549E-02, 2.949E-02 },
        { 5.00000E-01, 8.712E-02, 2.966E-02 }
    };
    return data;
}

bool testMassEnergyTransferAir()
{
    const std::string air_str = "Air, Dry (near sea level)";
    auto air = xraymc::Material<5>::byNistName(air_str).value();
    bool success = true;

    // std::cout << "Energy, NIST, xraymc\n";
    for (const auto& d : nistAir()) {
        //    std::cout << d.energy << ", " << d.atten << ", " << air.massEnergyTransferAttenuation(d.energy) << std::endl;
        const auto xraymc_atten = air.massEnergyTransferAttenuation(d.energy);
        success = success && (1 - d.atten / xraymc_atten) * 100 < 1.0;
    }
    return success;
}

int main(int argc, char* argv[])
{
    std::cout << "Test mass attenuation transfer against NIST data" << std::endl;

    auto stxt = [](bool v) -> std::string { return v ? " SUCCSESS " : " FAILED "; };

    bool success = true;
    success = success && testMassEnergyTransferAir();
    std::cout << "Air material: " << stxt(success) << std::endl;

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

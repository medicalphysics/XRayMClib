

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

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/personaldosimeter.hpp"

#include "xraymc/xraymc.hpp"

#include "xraymc/beams/flatmonoenergyfield.hpp"

#include <iostream>

bool test_dose()
{
    xraymc::PersonalDosimeter<false> dos;
    auto dir = dos.normalVector();

    auto s_vec = xraymc::vectormath::scale(dir, -1.0);

    xraymc::PencilBeam<> beam;
    beam.setDirection(s_vec);
    beam.setPosition(xraymc::vectormath::add(dos.center(), dir));
    beam.setEnergy(60);
    auto kerma = 1;
    beam.setAirKerma(kerma);
    beam.setNumberOfExposures(1);
    beam.setNumberOfParticlesPerExposure(1E6);

    for (double e = 0; e < 500; e = e + 5) {

        std::cout << e << ", " << dos.hp10ConvertionFactor(e) << std::endl;
    }

    return false;
}

bool test_angles()
{
    xraymc::PersonalDosimeter dos;
    xraymc::Particle p;
    for (std::size_t i = 1; i < 100; ++i) {

        double ang = -std::numbers::pi_v<double> / 2 + i * std::numbers::pi_v<double> / 100.0;
        p.dir = xraymc::vectormath::rotate({ 0, 0, -1 }, { 1, 0, 0 }, ang);

        auto w = dos.angularResponseWeight(p);
        std::cout << ang << "," << w << std::endl;
    }
    return true;
}

bool geo_test_angles()
{

    using D = xraymc::PersonalDosimeter<false>;
    using World = xraymc::World<D>;

    World world(1);
    auto& dos = world.addItem<D>("Dosimeter");

    
    auto dir = dos.normalVector();
    auto s_vec = xraymc::vectormath::scale(dir, -1.0);
    beam.setDirection(s_vec);
    beam.setPosition(xraymc::vectormath::add(dos.center(), dir));
    beam.setEnergy(60);
    auto kerma = 1;
    beam.setAirKerma(kerma);
    beam.setNumberOfExposures(30);
    beam.setNumberOfParticlesPerExposure(1E6);


}

void testAtt()
{

    auto mat_cand = xraymc::Material<5>::byNistName("Air, Dry (near sea level)");
    if (mat_cand) {
        auto& mat = mat_cand.value();
        std::vector<double> energy {
            1.00,
            1.50,
            2.00,
            3.00,
            3.20,
            3.20,
            4.00,
            5.00,
            6.00,
            8.00,
            10.00,
            15.00,
            20.00,
            30.00,
            40.00,
            50.00,
            60.00,
            80.00,
            100.00,
            150.00
        };
        std::cout << "Air kerma [keV/g]\n";
        for (auto e : energy)
            std::cout << e << ", " << e * mat.massEnergyTransferAttenuation(e) << std::endl;
    }
}

int main()

{
    // testAtt();

    bool success = true;
    success = success && test_angles();
    // success = success && test_dose();

    if (success)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}
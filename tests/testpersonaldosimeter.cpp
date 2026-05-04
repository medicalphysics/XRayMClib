

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

template <bool DS=false>
bool geo_test_angles()
{

    using D = xraymc::PersonalDosimeter<DS>;
    using World = xraymc::World<D>;

    World world(1);
    auto& dos = world.addItem<D>("Dosimeter");
    dos.setDirectionCosines({ 0, 1, 0 }, { 0, 0, 1 });

    xraymc::FlatMonoEnergyField<false> beam;
    beam.setEnergy(60);
    beam.setAirKerma(1.0);
    beam.setNumberOfExposures(30);
    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setLenght(10, 10);
    constexpr std::array<double, 3> dirX = { 0, -1, 0 };
    constexpr std::array<double, 3> dirY = { 0, 0, 1 };
    constexpr std::array<double, 3> pos = { 10, 0, 0 };
    beam.setPosition(pos);
    beam.setDirectionCosines(dirX, dirY);
    
    for (std::size_t i = 0; i < 90; i = i + 5) {
        const double angle = i * xraymc::DEG_TO_RAD();
        dos.setDirectionCosines({ 0, 1, 0 }, { 0, 0, 1 });
        dos.rotate(angle, { 0, 0, 1 });
        world.build();
        dos.clearDoseScored();
        xraymc::Transport::runConsole(world, beam, 0);

        std::cout << "Angle," << i;
        std::cout << ",AirKerma," << dos.doseScored(1).dose() << "," << dos.doseScored(1).numberOfEvents();
        std::cout << ",Hp10," << dos.doseScored(0).dose() << std::endl;
        dos.clearDoseScored();

        xraymc::VisualizeWorld viz(world);
        auto buffer = viz.createBuffer(512, 512);

        auto corner1 = xraymc::vectormath::add(beam.position(), xraymc::vectormath::add(xraymc::vectormath::scale(beam.directionCosines()[1], -0.5 * beam.lenghtY()), xraymc::vectormath::scale(beam.directionCosines()[0], 0.5 * beam.lenghtX())));
        auto corner2 = xraymc::vectormath::add(beam.position(), xraymc::vectormath::add(xraymc::vectormath::scale(beam.directionCosines()[1], -0.5 * beam.lenghtY()), xraymc::vectormath::scale(beam.directionCosines()[0], -0.5 * beam.lenghtX())));
        auto corner3 = xraymc::vectormath::add(beam.position(), xraymc::vectormath::add(xraymc::vectormath::scale(beam.directionCosines()[1], 0.5 * beam.lenghtY()), xraymc::vectormath::scale(beam.directionCosines()[0], 0.5 * beam.lenghtX())));
        auto corner4 = xraymc::vectormath::add(beam.position(), xraymc::vectormath::add(xraymc::vectormath::scale(beam.directionCosines()[1], 0.5 * beam.lenghtY()), xraymc::vectormath::scale(beam.directionCosines()[0], -0.5 * beam.lenghtX())));

        viz.addLineSegment(corner1, xraymc::vectormath::add(corner1, xraymc::vectormath::scale(beam.direction(), 9.0)), .1);
        viz.addLineSegment(corner2, xraymc::vectormath::add(corner2, xraymc::vectormath::scale(beam.direction(), 9.0)), .1);
        viz.addLineSegment(corner3, xraymc::vectormath::add(corner3, xraymc::vectormath::scale(beam.direction(), 9.0)), .1);
        viz.addLineSegment(corner4, xraymc::vectormath::add(corner4, xraymc::vectormath::scale(beam.direction(), 9.0)), .1);

        viz.addAABBoutline(dos.AABB(), 0.1);

        viz.setAzimuthalAngleDeg(60);
        viz.setPolarAngleDeg(60);
        viz.suggestFOV(1);
        viz.generate(world, buffer);
        std::string name = std::to_string(i) + ".png";
        viz.savePNG(name, buffer);
    }
    /* xraymc::VisualizeWorld viz(world);
     auto buffer = viz.createBuffer(512, 512);
     viz.setAzimuthalAngleDeg(45);
     viz.setPolarAngleDeg(45);
     viz.suggestFOV(1);
     viz.generate(world, buffer);
     viz.savePNG("test.png", buffer);
 */
    return true;
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
    // success = success && test_angles();
    success = success && geo_test_angles<true>();

    if (success)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}
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

Copyright 2023 Erlend Andersen
*/

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

bool testForcedinteractions()
{
    using Sphere = dxmc::WorldSphere<5, 1, false>;
    using SphereF = dxmc::WorldSphere<5, 1, true>;

    using World = dxmc::World<Sphere>;
    using WorldF = dxmc::World<SphereF>;

    World w(2);
    WorldF wf(2);

    constexpr double radii = 0.3;
    auto& sphere1 = w.addItem<Sphere>({ radii });
    auto& sphere2 = w.addItem<Sphere>({ radii, { radii * 2 + radii * 100, 0, 0 } });
    auto& spheref1 = wf.addItem<SphereF>({ radii });
    auto& spheref2 = wf.addItem<SphereF>({ radii, { radii * 2 + radii * 100, 0, 0 } });

    auto material_water = dxmc::Material<5>::byNistName("Water, Liquid").value();
    sphere1.setMaterial(material_water, 1.2);
    sphere2.setMaterial(material_water, 1.2);
    spheref1.setMaterial(material_water, 1.2);
    spheref2.setMaterial(material_water, 1.2);

    w.build();
    wf.build();

    dxmc::PencilBeam<true> beam({ -100, 0, 0 }, { 1, 0, 0 }, 20);
    beam.setNumberOfExposures(48);
    beam.setNumberOfParticlesPerExposure(1e4);

    dxmc::Transport transport;
    transport.runConsole(wf, beam);
    transport.runConsole(w, beam);

    auto dose1 = sphere1.doseScored().dose();
    auto dose2 = sphere2.doseScored().dose();
    auto dosef1 = spheref1.doseScored().dose();
    auto dosef2 = spheref2.doseScored().dose();

    auto diff1 = std::abs(dose1 - dosef1);
    auto diff2 = std::abs(dose2 - dosef2);

    auto sttd1 = sphere1.doseScored().standardDeviation();
    auto sttdf1 = spheref1.doseScored().standardDeviation();

    auto sttd2 = sphere2.doseScored().standardDeviation();
    auto sttdf2 = spheref2.doseScored().standardDeviation();

    std::cout << "Testing forced interactions in worldsphere" << std::endl;
    std::cout << "Random " << dose1 << ", " << dose2 << " stddev: " << sttd1 << ", " << sttd2 << std::endl;
    std::cout << "Forced " << dosef1 << ", " << dosef2 << " stddev: " << sttdf1 << ", " << sttdf2 << std::endl;
    std::cout << "Diff: " << diff1 << ", " << diff2 << std::endl;

    constexpr auto test_coeff = 2.57;

    const auto test1 = diff1 / std::sqrt(sttd1 * sttd1 + sttdf1 * sttdf1);
    const auto test2 = diff2 / std::sqrt(sttd2 * sttd2 + sttdf2 * sttdf2);

    std::cout << "Test statistics: " << test1 << ", " << test2 << std::endl;
    auto success = (test1 < test_coeff) && (test2 < test_coeff);
    if (success)
        std::cout << "SUCCESS: differences within statistical limits" << std::endl;
    else
        std::cout << "FAILURE: differences outside statistical limits" << std::endl;
    return success;
}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && testForcedinteractions();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

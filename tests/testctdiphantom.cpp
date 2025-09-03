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
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

bool testTracker()
{

    using CTDI = dxmc::CTDIPhantom<5, 1, false>;
    using Sphere = dxmc::WorldSphere<5, 1>;
    using World = dxmc::World<CTDI, Sphere>;

    World w;
    w.reserveNumberOfItems(2);
    auto& ctdi = w.addItem<CTDI>({ 16, 15 });
    auto& sphere = w.addItem<Sphere>({ 3, { 0, 0, 15+3.01 } });
    w.build(100);

    dxmc::PencilBeam<true> beam({ -25, -25, -15 }, { 1, 1, 0.5 }, 60);
    beam.setNumberOfExposures(16);
    beam.setNumberOfParticlesPerExposure(1e4);

    dxmc::Transport transport;
    //transport.setNumberOfThreads(1);
    transport.runConsole(w, beam);

    dxmc::VisualizeWorld viz(w);

    const auto& tracker = sphere.particleTracker();
    viz.addParticleTracks(tracker, 0.02);

    auto buffer = viz.createBuffer(2048, 2048);
    viz.setAzimuthalAngleDeg(60);
    viz.setPolarAngleDeg(0);
    viz.setDistance(100);
    viz.suggestFOV(2.0);
    viz.generate(w, buffer);
    viz.savePNG("test.png", buffer);

    return true;
}

bool testForcedinteractions()
{
    using CTDIf = dxmc::CTDIPhantom<5, 1, true>;
    using CTDI = dxmc::CTDIPhantom<5, 1, false>;
    using Worldf = dxmc::World<CTDIf>;
    using World = dxmc::World<CTDI>;

    World w;
    auto& ctdi = w.addItem<CTDI>({ 8 });
    w.build();

    Worldf wf;
    auto& ctdif = wf.addItem<CTDIf>({ 8 });
    wf.build();

    dxmc::PencilBeam beam({ -100, -100, 0 }, { 1, 1, 0 }, 60);
    beam.setNumberOfExposures(480);
    beam.setNumberOfParticlesPerExposure(1e6);

    dxmc::Transport transport;

    transport(wf, beam);
    transport(w, beam);

    auto dose_centerf = ctdif.doseScored(1);
    auto dose_center = ctdi.doseScored(1);
    auto df = dose_centerf.dose();
    auto d = dose_center.dose();
    auto sd = dose_center.standardDeviation();
    auto sdf = dose_centerf.standardDeviation();

    std::cout << "Forced: " << dose_centerf.dose() << " " << dose_centerf.standardDeviation() << std::endl;
    std::cout << "Random: " << dose_center.dose() << " " << dose_center.standardDeviation() << std::endl;

    constexpr auto test_coeff = 2.57;
    auto test = std::abs(d - df) / std::sqrt(sd * sd + sdf * sdf);

    std::cout << "Statistic: " << test << ", threshold: " << test_coeff << std::endl;

    auto success = test < test_coeff;
    if (success)
        std::cout << "SUCCESS: differences within statistical limits" << std::endl;
    else
        std::cout << "FAILURE: differences outside statistical limits" << std::endl;
    return success;
}

int main(int argc, char* argv[])
{
    auto success = true;
    //success = success && testTracker();
    success = success && testForcedinteractions();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

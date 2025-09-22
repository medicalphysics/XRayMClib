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

#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"

#include <iostream>
#include <string>
#include <vector>

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    xraymc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <std::size_t N = 5, int L = 2>
void testfluencescore()
{
    using Box = xraymc::WorldBox<N, L>;
    using FluenceScore = xraymc::FluenceScore;
    using World = xraymc::World<Box, FluenceScore>;
    using Material = xraymc::Material<N>;

    World world;

    auto& box = world.template addItem<Box>({ 10 });
    auto water = Material::byNistName("Water, Liquid").value();
    const double density = 1;
    box.setMaterial(water, density);
    auto& fluence = world.template addItem<FluenceScore>({ 20, { 0, 0, 15 }, { 0, 0, 1 } });

    world.build();

    xraymc::VisualizeWorld viz(world);
    viz.addLineProp({ -1000, 0, 0 }, { 1, 0, 0 }, 1000, 0.1);
    int height = 528;
    int width = 528;
    std::vector<double> buffer(height * width * 4, 1);
    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(500);
        viz.setPolarAngle(std::numbers::pi_v<double> / 3.0);
        viz.setAzimuthalAngle(std::numbers::pi_v<double> * i / 6.0);
        viz.suggestFOV();
        viz.generate(world, buffer, width, height);
        std::string name = "fluence_score" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer, width, height);
    }

    xraymc::Tube tube;
    tube.setVoltage(140);
    tube.setAlFiltration(8);
    tube.setEnergyResolution(1);
    const auto ts = tube.getSpecter();

    xraymc::IsotropicBeam<> beam({ -1000, 0, 0 }, { 0, 1, 0, 0, 0, 1 });
    beam.setEnergySpecter(ts);
    beam.setCollimationHalfAngles(0, 0, 0, 0);

    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setNumberOfExposures(100);

    xraymc::Transport transport;
    runDispatcher(transport, world, beam);

    const auto specter = fluence.getSpecter();
    const auto N = beam.numberOfExposures() * beam.numberOfParticles();

    std::cout << "Energy, IntensityIn, Energy, IntensityOut\n";
    for (std::size_t i = 0; i < ts.size(); ++i) {
        std::cout << ts[i].first << ", " << N * ts[i].second << ", ";
        std::cout << specter[i + 1].first << ", " << specter[i + 1].second << std::endl;
    }
}

int main()
{
    std::cout << "Testing fluence scoring\n";
    bool success = true;

    testfluencescore<5, 2>();
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
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

#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>
#include <string>
#include <vector>

template <std::size_t N = 15, int L = 2>
void testfluencescore()
{
    using Sphere = xraymc::WorldSphere<N, L, false>;
    using FluenceScore = xraymc::FluenceScore;
    using World = xraymc::World<Sphere, FluenceScore>;
    using Material = xraymc::Material<N>;

    World world;
    world.reserveNumberOfItems(2);
    auto& sphere = world.template addItem<Sphere>();
    sphere.setRadius(1);

    auto& fscore = world.template addItem<FluenceScore>();
    fscore.setRadius(1);
    fscore.setCenter({ 0, 2, 0 });
    fscore.setPlaneNormal({ 0, -1, 0 });
    fscore.setEnergyStep(1.0);

    auto material = Material::byChemicalFormula("PbBiSi").value();

    const double density = 11.0;
    sphere.setMaterial(material, density);

    world.build();

    xraymc::PencilBeam<false> beam({ -2, 0, 0 }, { 1, 0, 0 }, 135);
    beam.setNumberOfExposures(500);
    beam.setNumberOfParticlesPerExposure(1E6);

    xraymc::Transport transport;
    // transport.runConsole(world, beam);

    xraymc::VisualizeWorld viz(world);

    viz.addLineProp({ -2, 0, 0 }, { 1, 0, 0 }, 1, 0.1);
    int height = 528 * 2;
    int width = 528 * 2;
    std::vector<double> buffer(height * width * 4, 1);
    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(50);
        viz.setPolarAngleDeg(30);
        viz.setAzimuthalAngleDeg(i * 360.0 / 12.0);
        viz.suggestFOV(3);
        viz.generate(world, buffer, width, height);
        std::string name = "fluence_score" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer, width, height);
    }

    const auto specter = fscore.getFluenceSpecter();
    std::cout << "Energy keV, I (" << N << ")" << std::endl;
    for (const auto [energy, I] : specter) {
        if (energy < beam.energy() + 1)
            std::cout << energy << "," << I << std::endl;
    }
}

int main()
{
    std::cout << "Testing high Z material scatter\n";
    bool success = true;

    testfluencescore<63, 2>();
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
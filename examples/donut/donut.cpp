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

Copyright 2025 Erlend Andersen
*/

#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

// Example of xraymc
//  This example uses a voxelized donut to illustrate some usecases of xraymc

void example1()
{
    // make and display a voxelized donut with an aluminium ball suspended in the middle
    // first we define some constants

    // Select low energy correction, possible values are
    // 0: No electron binding corrections
    // 1: Livermore correction for whole atoms
    // 2: Correction by individual atomic shells
    constexpr int LOWENERGYCORRECTION = 2;

    // lets start by defining a voxelized donut as an axis aligned grid
    using VoxelItem = xraymc::AAVoxelGrid<16, 2, 0>;
    using Sphere = xraymc::WorldSphere<16, 2, false>; // If we want to put an aluminum ball inside the donut

    // dimensions of the voxelgrid
    std::array<std::size_t, 3> dimensions = { 128, 128, 128 };
    const auto N_voxels = dimensions[0] * dimensions[1] * dimensions[2];

    // voxel spacing in cm
    std::array<double, 3> spacing = { 0.1, 0.1, 0.1 };

    // Then we crate to buffers, one with material indices and one with densities in each voxel
    std::vector<double> densities(N_voxels);
    std::vector<std::uint8_t> materialIndices(N_voxels);

    // We also need a vector of materials the voxel grid consists of
    using Material = xraymc::Material<>;
    std::vector<Material> materials;
    // air is the first material
    materials.push_back(Material::byNistName("Air, Dry (near sea level)").value());
    // water is the second
    materials.push_back(Material::byChemicalFormula("H2O").value());
    // fat is the third material (it's a donut after all)
    materials.push_back(Material::byChemicalFormula("H(CHOH)3H").value());

    // Lets also get material densities (in g/cm3) for convienence
    std::vector<double> material_densities(3);
    material_densities[0] = xraymc::NISTMaterials::density("Air, Dry (near sea level)"); // air
    material_densities[1] = 1.0; // water
    material_densities[2] = 0.92; // fat

    // making the donut by assigning values to the densities and materialIndices buffers
    for (std::size_t i = 0; i < dimensions[0]; ++i)
        for (std::size_t j = 0; j < dimensions[1]; ++j)
            for (std::size_t k = 0; k < dimensions[2]; ++k) {

                // finding donut center
                const std::array center = {
                    dimensions[0] * spacing[0] * 0.5,
                    dimensions[1] * spacing[1] * 0.5,
                    dimensions[2] * spacing[2] * 0.5
                };
                const auto flatBufferIndex = i + j * dimensions[0] + k * dimensions[0] * dimensions[1];

                // cordinates of each voxel with respect to center
                const auto x = i * spacing[0] - center[0];
                const auto y = j * spacing[1] - center[1];
                const auto z = k * spacing[2] - center[2];

                constexpr double R = 4.0; // cm
                constexpr double r = 2; // cm

                bool inside_donut = std::pow(std::sqrt(x * x + y * y) - R, 2) + z * z < r * r;
                if (inside_donut) {
                    if (z >= 0.0) {
                        // The top consists of fat
                        materialIndices[flatBufferIndex] = 2;
                        densities[flatBufferIndex] = material_densities[2];
                    } else {
                        // Bottom is water
                        materialIndices[flatBufferIndex] = 1;
                        densities[flatBufferIndex] = material_densities[1];
                    }
                } else {
                    // Outside donut is air
                    materialIndices[flatBufferIndex] = 0;
                    densities[flatBufferIndex] = material_densities[0];
                }
            }

    // We have a donut and can start defining a world
    
    xraymc::World<VoxelItem, Sphere> world;
    world.reserveNumberOfItems(2);
    auto& voxelItem = world.template addItem<VoxelItem>({ dimensions, spacing, densities, materialIndices, materials }, "Donut");

    // adding an aluminum ball inside the donut
    Sphere sphere(1, { 0, 0, 0 });
    auto aluminum = xraymc::Material<16>::byZ(13).value();
    sphere.setMaterial(aluminum, 2.7); // aluminum
    world.addItem(sphere, "Sphere");

    world.build();

    // lets show the world
    xraymc::VisualizeWorld viz(world);
    auto buffer = viz.template createBuffer<double>(1024, 1024); // image size
    viz.setDistance(60); // camera distance
    viz.setAzimuthalAngleDeg(40); // camera angles
    viz.setPolarAngleDeg(30);
    viz.suggestFOV(2); // zoom = 2
    viz.setColorOfItem<std::uint8_t>(world.getItemPointerFromName("Donut"), { 255, 192, 203 }); // making it pink
    viz.generate(world, buffer);
    viz.savePNG("donut.png", buffer);

    // lets irradiate the donut with a xraybeam
    xraymc::DXBeam beam;
    beam.setPosition({ 0, 0, 30 }); // 30 cm above the donut
    beam.setDirectionCosines({ 1, 0, 0, 0, -1, 0 }); // field directions, the cross product is the beam direction
    beam.setCollimationHalfAnglesDeg({ 5, 5 }); // Half angles of the beam collimation
    beam.setNumberOfExposures(100); // Number of jobs
    beam.setNumberOfParticlesPerExposure(1e6); // particles per job
    beam.setTubeVoltage(60); // tube voltage
    beam.addTubeFiltrationMaterial(13, 2.0); // 2.0 mm Al filtration

    // setting up simulation manager
    xraymc::Transport transport;
    const auto number_of_threads = std::thread::hardware_concurrency(); // Using all threads for your machine
    auto simulation_time = transport.runConsole(world, beam, number_of_threads);

    // finding max donut dose
    double max_dose = 0;
    for (const auto& doseScore : voxelItem.getDoseScores())
        max_dose = std::max(max_dose, doseScore.dose());
    // Show the result
    viz.addLineProp(beam, 30, 0.02);
    viz.addColorByValueItem(world.getItemPointerFromName("Donut"));
    viz.addColorByValueItem(world.getItemPointerFromName("Sphere"));
    viz.setColorByValueMinMax(0.0, max_dose);
    viz.generate(world, buffer);
    viz.savePNG("donut_dose.png", buffer);
}

int main()
{
    example1();
    return EXIT_SUCCESS;
}
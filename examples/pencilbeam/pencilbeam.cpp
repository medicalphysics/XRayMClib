

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"

#include <algorithm>
#include <iostream>
#include <thread>

void example()
{
    // Create a thin aluminium cylinder and print depth dose
    std::cout << "Example Pencilbeam\nTransport of monoenergetic photons in a thin long aluminium cylinder\n";

    // Start creating a depthscore cylinder
    constexpr int N_ATOMIC_SHELLS = 5; // Number of atomic shells to consider binding energies for (5 is more than sufficient).
    constexpr int LOW_ENERGY_CORRECTION = 1; /* 0: No binding energy correction, 1: Livermore correction, 2; impulse approx. correction*/

    using Cylinder = dxmc::DepthDose<N_ATOMIC_SHELLS, LOW_ENERGY_CORRECTION>;
    using Room = dxmc::EnclosedRoom<N_ATOMIC_SHELLS, LOW_ENERGY_CORRECTION>;

    // Create a world that can consist of one or more depthdose objects
    dxmc::World<Cylinder, Room> world;
    // Reserve number of items in world.
    world.reserveNumberOfItems(2);

    // Adding a depthdose object
    auto& cylinder = world.template addItem<Cylinder>({ 1 /* cm radius */, 10 /* cm lenght */ }, "Cylinder");

    // Set material and density

    auto aluminium = dxmc::Material<N_ATOMIC_SHELLS>::byZ(13).value();
    cylinder.setMaterial(aluminium, 2.27 /* g/cm3 */);

    // Optional set cylinder material to water
    // auto water = dxmc::Material<N_ATOMIC_SHELLS>::byChemicalFormula("H2O").value();
    // cylinder.setMaterial(water, 1.0);

    // Adding room with walls of concrete
    auto& room = world.template addItem<Room>({ 2 /*cm wall thickness*/, 200 /*cm inner walls sizes*/ }, "Room");
    auto concrete = dxmc::Material<N_ATOMIC_SHELLS>::byNistName("Concrete, Ordinary").value();
    auto concrete_density = dxmc::NISTMaterials::density("Concrete, Ordinary");
    room.setMaterial(concrete, concrete_density);
    //    Example for constructing other materials
    //    auto water = dxmc::Material<N_ATOMIC_SHELLS>::byChemicalFormula("H2O").value();

    std::cout << "with radius " << cylinder.radius() << " cm and height " << cylinder.length() << " cm\n";

    // Building world
    world.build();

    // Define a radiation source
    dxmc::PencilBeam<> beam({ 0, 0, -10 } /* position */, { 0, 0, 1 } /* direction */);
    beam.setNumberOfExposures(64); // number of jobs
    beam.setNumberOfParticlesPerExposure(1000000); // histories per job

    // Run simulation
    auto nThreads = std::max(std::thread::hardware_concurrency(), std::uint32_t { 1 });
    auto time_elapsed = dxmc::Transport::runConsole(world, beam, nThreads, true);

    // Get max dose and print some values
    std::cout << "Depth dose in cylinder for " << beam.numberOfParticles() << " photons of " << beam.energy() << " keV\n";
    std::cout << "Simulation time: " << time_elapsed << std::endl;
    std::cout << "Depth [cm], Dose per Air Kerma [mGy/mGy], Uncertainty [%], #Events\n";
    double max_dose = 0;
    for (std::size_t i = 0; i < cylinder.resolution(); ++i) {
        auto dose_val = cylinder.doseScored(i).dose();
        max_dose = std::max(max_dose, dose_val);
        std::cout << (cylinder.length() * (i + 0.5)) / cylinder.resolution() << ", ";
        std::cout << cylinder.doseScored(i).dose() << ", ";
        std::cout << cylinder.doseScored(i).relativeUncertainty() * 100 << ", ";
        std::cout << cylinder.doseScored(i).numberOfEvents() << std::endl;
    }

    // Generate some images
    dxmc::VisualizeWorld viz(world);
    auto buffer = viz.template createBuffer<double>(1024, 1024);
    viz.setDistance(60);
    viz.setAzimuthalAngleDeg(90);
    viz.setPolarAngleDeg(30);
    viz.suggestFOV(4); // zoom = 2
    viz.setColorOfItem<std::uint8_t>(world.getItemPointerFromName("Cylinder"), { 255, 192, 203 }); // making it pink
    viz.generate(world, buffer);
    viz.savePNG("cylinder.png", buffer);

    viz.setColorByValueMinMax(0.0, max_dose); // color by dose value
    viz.addColorByValueItem(world.getItemPointerFromName("Cylinder"));
    viz.generate(world, buffer);
    viz.savePNG("cylinder_dose.png", buffer);
}

int main(int argc, char* argv[])
{
    example();
    return EXIT_SUCCESS;
}
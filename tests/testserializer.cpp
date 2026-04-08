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

#include "xraymc/serializer.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include "xraymc/beams/cbctbeam.hpp"
#include "xraymc/beams/ctsequentialbeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/ctspiraldualenergybeam.hpp"
#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/isotropicbeamcircle.hpp"
#include "xraymc/beams/isotropiccircularbeam.hpp"
#include "xraymc/beams/isotropiccircularmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"

#include <iostream>

template <typename T>
bool testItem(const T& item)
{

    xraymc::Serializer s;

    auto buffer = xraymc::Serializer::getEmptyBuffer();
    xraymc::Serializer::serializeItem(item, buffer);

    xraymc::Serializer::write("testitem.xr", buffer);

    auto rbuffer = s.read("testitem.xr").value();

    auto name = T::magicID();
    std::vector<char> itemBuffer;
    s.deserializeItem(name, itemBuffer, rbuffer);
    auto item_opt = T::deserialize(itemBuffer);
    if (item_opt) {
        auto& item_read = item_opt.value();
        return item_read == item;
    }

    return false;
}

bool testString()
{
    xraymc::Serializer s;

    std::vector<std::string> ss;
    ss.push_back("Test 1");
    ss.push_back("Testing 2 ");

    auto buffer = s.getEmptyBuffer();
    s.serialize(ss, buffer);

    std::vector<std::string> ss_out;
    auto end = s.deserialize(ss_out, buffer);
    return ss_out == ss;
}

template <xraymc::BeamType B>
bool testBeams()
{
    B beam;

    xraymc::Serializer s;
    auto save_buffer = xraymc::Serializer::getEmptyBuffer();
    const auto beam_buffer = beam.serialize();
    xraymc::Serializer::serializeItem(beam.magicID(), beam_buffer, save_buffer);

    auto id = beam.magicID();
    auto debuffer = xraymc::Serializer::getEmptyBuffer();
    auto debuffer_end = xraymc::Serializer::deserializeItem(id, debuffer, save_buffer);
    auto recon_beam_opt = B::deserialize(debuffer);
    return recon_beam_opt ? true : false;
}

void example()
{
    // Create a thin aluminium cylinder and print depth dose
    std::cout << "Example Pencilbeam\nTransport of monoenergetic photons in a thin long aluminium cylinder\n";

    // Start creating a depthscore cylinder
    constexpr int N_ATOMIC_SHELLS = 5; // Number of atomic shells to consider binding energies for (5 is more than sufficient).
    constexpr int LOW_ENERGY_CORRECTION = 1; /* 0: No binding energy correction, 1: Livermore correction, 2; impulse approx. correction*/

    using Cylinder = xraymc::DepthDose<N_ATOMIC_SHELLS, LOW_ENERGY_CORRECTION>;
    using Room = xraymc::EnclosedRoom<N_ATOMIC_SHELLS, LOW_ENERGY_CORRECTION>;

    // Create a world that can consist of one or more depthdose objects
    xraymc::World<Cylinder, Room> world;
    // Reserve number of items in world.
    world.reserveNumberOfItems(2);

    // Adding a depthdose object
    auto& cylinder = world.template addItem<Cylinder>({ 1 /* cm radius */, 10 /* cm lenght */ }, "Cylinder");

    // Set material and density
    auto aluminium = xraymc::Material<N_ATOMIC_SHELLS>::byZ(13).value();
    cylinder.setMaterial(aluminium, 2.27 /* g/cm3 */);

    // Optional set cylinder material to water
    // auto water = xraymc::Material<N_ATOMIC_SHELLS>::byChemicalFormula("H2O").value();
    // cylinder.setMaterial(water, 1.0);

    // Adding room with walls of concrete
    auto& room = world.template addItem<Room>({ 2 /*cm wall thickness*/, 200 /*cm inner walls sizes*/ }, "Room");
    auto concrete = xraymc::Material<N_ATOMIC_SHELLS>::byNistName("Concrete, Ordinary").value();
    auto concrete_density = xraymc::NISTMaterials::density("Concrete, Ordinary");
    room.setMaterial(concrete, concrete_density);
    //    Example for constructing other materials
    //    auto water = xraymc::Material<N_ATOMIC_SHELLS>::byChemicalFormula("H2O").value();

    std::cout << "with radius " << cylinder.radius() << " cm and height " << cylinder.length() << " cm\n";

    // Building world
    world.build();

    auto save_buffer = xraymc::Serializer::getEmptyBuffer();
    const auto world_buffer = world.serialize();
    xraymc::Serializer::serializeItem(world.magicID(), world_buffer, save_buffer);

    auto dename = world.magicID();
    std::vector<char> debuffer;
    auto test = xraymc::Serializer::deserializeItem(dename, debuffer, save_buffer);

    auto newWorld = world.deserialize(debuffer);

    // Define a radiation source
    xraymc::PencilBeam<> beam({ 0, 0, -10 } /* position */, { 0, 0, 1 } /* direction */);
    beam.setNumberOfExposures(64); // number of jobs
    beam.setNumberOfParticlesPerExposure(1000000); // histories per job
    beam.setEnergy(60.0);

    // Run simulation
    auto nThreads = std::max(std::thread::hardware_concurrency(), std::uint32_t { 1 });
    auto time_elapsed = xraymc::Transport::runConsole(world, beam, nThreads, true);

    // Get max dose and print some values
    std::cout << "Depth dose in cylinder for " << beam.numberOfParticles() << " photons of " << beam.energy() << " keV\n";
    std::cout << "Material: mass. att. coeff: " << aluminium.attenuationValues(beam.energy()).sum() << ", density: " << cylinder.density() << std::endl;
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
    xraymc::VisualizeWorld viz(world);
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

int main()
{

    // example();

    bool success = testString();
    success = success && testBeams<xraymc::PencilBeam<>>();
    success = success && testBeams<xraymc::IsotropicMonoEnergyBeamCircle<>>();
    success = success && testBeams<xraymc::IsotropicMonoEnergyBeam<>>();
    success = success && testBeams<xraymc::IsotropicCircularMonoEnergyBeam<>>();
    success = success && testBeams<xraymc::IsotropicCircularBeam<>>();
    success = success && testBeams<xraymc::IsotropicBeamCircle<>>();
    success = success && testBeams<xraymc::IsotropicBeam<>>();
    success = success && testBeams<xraymc::CTSpiralDualEnergyBeam<>>();
    success = success && testBeams<xraymc::CTSpiralBeam<>>();
    success = success && testBeams<xraymc::CTSequentialBeam<>>();
    success = success && testBeams<xraymc::CBCTBeam<>>();
    success = success && testBeams<xraymc::DXBeam<>>();

    xraymc::WorldSphere<12, 2, true> sphere;
    sphere.setCenter({ 1, 2, 3 });
    sphere.setMaterialDensity(1.9);
    sphere.setRadius(5);

    success = success && testItem(sphere);

    xraymc::WorldCylinder<16, 2> cylinder;
    success = success && testItem(cylinder);

    xraymc::AAVoxelGrid<12, 2, 255> vgrid;
    success = success && testItem(vgrid);

    xraymc::WorldBoxGrid<15, 2> boxgrid;
    success = success && testItem(boxgrid);

    xraymc::WorldBox<12, 2> box;
    success = success && testItem(box);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
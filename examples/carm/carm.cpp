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

#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/particletracker.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/triangulatedopensurface.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include "phantomreader.hpp"

#include <iostream>
#include <vector>

dxmc::AAVoxelGrid<5, 1, 0> testPhantom()
{
    auto d = ICRP110PhantomReader::readFemalePhantom("AF.dat", "AF_media.dat", "AF_organs.dat");

    dxmc::AAVoxelGrid<5, 1, 0> phantom;
    using Material = dxmc::Material<5>;
    std::vector<Material> materials;
    for (auto& w : d.mediaComposition()) {
        auto mat_cand = Material::byWeight(w);
        if (mat_cand)
            materials.push_back(mat_cand.value());
        else
            throw std::runtime_error("error");
    }
    phantom.setData(d.dimensions(), d.densityData(), d.mediaData(), materials);
    phantom.setSpacing(d.spacing());
    return phantom;
}

template <bool TRACK>
dxmc::TetrahedalMesh<5, 1, !TRACK> readICRP145Phantom(std::array<int, 3> depth = { 8, 8, 8 }, bool female = true)
{
    const std::string name = female ? "MRCP_AF" : "MRCP_AM";
    const std::string elefile = name + ".ele";
    const std::string nodefile = name + ".node";
    const std::string mediafile = name + "_media.dat";
    const std::string organfile = "icrp145organs.csv";

    dxmc::TetrahedalmeshReader<5, 1, !TRACK> reader(nodefile, elefile, mediafile, organfile);
    reader.rotate({ 0, 0, 1 }, std::numbers::pi_v<double>);
    return reader.getMesh(depth);
}

template <bool TRACK>
dxmc::TetrahedalMesh<5, 1, !TRACK> readICRP145Phantom(int depth = 8, bool female = true)
{
    return readICRP145Phantom<TRACK>({ depth, depth, depth }, female);
}

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    dxmc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << std::flush << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <bool TRACK = false>
void vizualize()
{
    using CTDIPhantom = dxmc::CTDIPhantom<5, 1>;
    using Mesh = dxmc::TriangulatedMesh<5, 1>;
    using Surface = dxmc::TriangulatedOpenSurface<5, 1>;
    using Sphere = dxmc::WorldSphere<5, 1>;
    using VGrid = dxmc::AAVoxelGrid<5, 1, 0>;
    using Room = dxmc::EnclosedRoom<5, 1>;
    using TetMesh = dxmc::TetrahedalMesh<5, 1, !TRACK>;
    using Box = dxmc::WorldBox<5, 1>;
    using World = dxmc::World<Mesh, Sphere, VGrid, Room, Surface, TetMesh, Box>;

    World world {};
    world.reserveNumberOfItems(8);

    const auto carbon = dxmc::Material<5>::byZ(6).value();
    const auto carbon_atom = dxmc::AtomHandler::Atom(6);
    const auto carbon_dens = dxmc::AtomHandler::Atom(6).standardDensity;

    auto& carm = world.template addItem<Mesh>({ "carm.stl" });
    carm.setMaterial(carbon, carbon_dens);
    carm.translate({ 10, 0, 0 });

    auto& table = world.template addItem<Mesh>({ "table.stl" });
    table.translate({ -30, 0, -20 });
    table.setMaterial(carbon, carbon_dens);

    const auto lead = dxmc::Material<5>::byZ(82).value();
    const auto lead_atom = dxmc::AtomHandler::Atom(82);
    const auto lead_dens = dxmc::AtomHandler::Atom(82).standardDensity;

    dxmc::STLReader stlreader;

    stlreader.setFilePath("ceilingshield.stl");
    auto ceilingshield_tri = stlreader();
    std::for_each(std::execution::par_unseq, ceilingshield_tri.begin(), ceilingshield_tri.end(), [](auto& tri) {
        tri.rotate(std::numbers::pi_v<double> / 2, { 1, 0, 0 });
        tri.rotate(std::numbers::pi_v<double> / 2 + std::numbers::pi_v<double> / 4, { 0, 0, 1 });
    });
    auto& ceilingshield = world.template addItem<Surface>({ ceilingshield_tri });
    ceilingshield.translate({ -5, -35, 50 });
    ceilingshield.scale(0.5);
    ceilingshield.setMaterial(lead, lead_dens);
    ceilingshield.setSurfaceThickness(0.1);

    auto& tableBox = world.template addItem<Box>({ { -80, -28, -120, 0, -27.9, 0 } });
    tableBox.setMaterial(lead, lead_dens);

    auto& room = world.template addItem<Room>();
    room.setInnerRoomAABB({ -350, -300, -120.1, 350, 300, 150 });
    room.setWallThickness(2);
    room.setMaterial(lead, lead_dens * 0.2 / 2.0);

    auto& phantom = world.addItem(testPhantom());
    phantom.rollAxis(2, 0);
    phantom.rollAxis(2, 1);
    phantom.flipAxis(2);
    auto table_aabb = table.AABB();
    auto phantom_aabb = phantom.AABB();
    phantom.translate({ -40, 0, table_aabb[5] - phantom_aabb[2] });

    stlreader.setFilePath("blanket.stl");
    auto blanket_tri = stlreader();
    // std::for_each(std::execution::par_unseq, blanket_tri.begin(), blanket_tri.end(), [](auto& tri) { tri.rotate(std::numbers::pi_v<double>, { 0, 0, 1 }); });
    auto& blanket = world.template addItem<Surface>({ blanket_tri });
    auto bc = blanket.center();
    blanket.translate(dxmc::vectormath::scale(bc, -1.0));
    blanket.translate({ -20, 0, table_aabb[5] - phantom_aabb[2] + 7 });
    blanket.setMaterial(lead, lead_dens);
    blanket.setSurfaceThickness(0.1);

    auto& doctor = world.template addItem<TetMesh>(readICRP145Phantom<TRACK>({ 64, 64, 256 }, true));
    const auto doctor_aabb = doctor.AABB();
    doctor.translate({ -40, -40, -doctor_aabb[2] - 120 });

    world.build();

    // adding beam
    using Beam = dxmc::DXBeam<TRACK>;
    const std::array<double, 3> source_pos = { 10, 0, -64 };
    Beam beam(source_pos);
    beam.setBeamSize(6, 6, 114);
    if constexpr (TRACK) {
        beam.setNumberOfExposures(56);
        beam.setNumberOfParticlesPerExposure(100000);
    } else {
        beam.setNumberOfExposures(2048 * 4);
        beam.setNumberOfParticlesPerExposure(1000000);
    }
    beam.setDAPvalue(25);

    dxmc::Transport transport;
    runDispatcher(transport, world, beam);

    double max_doctor_dose = 0;
    for (const auto& tet : doctor.tetrahedrons()) {
        max_doctor_dose = std::max(max_doctor_dose, tet.doseScored().dose());
    }

    std::cout << "Max doctor dose " << max_doctor_dose << " mGy, dose norm: " << std::endl;

    dxmc::VisualizeWorld viz(world);
    if constexpr (TRACK) {
        viz.addParticleTracks(doctor.particleTracker(), 0.1);
    } else {
        for (const auto& item : world.items()) {
            if (std::holds_alternative<TetMesh>(item))
                viz.addColorByValueItem(&item);
            if (std::holds_alternative<VGrid>(item))
                viz.addColorByValueItem(&item);
        }
    }

    viz.setColorByValueMinMax(0, 0.00001);
    constexpr int res = 2;
    constexpr double zoom = 3;
    auto buffer = viz.template createBuffer<double>(1024 * res, 1024 * res);
    viz.addLineProp(beam, 114, .2);
    viz.setDistance(400);

    std::vector<double> angles;
    for (std::size_t i = 0; i < 12; ++i)
        angles.push_back(i * 30);

    viz.setAzimuthalAngleDeg(60);
    for (auto a : angles) {
        viz.setPolarAngleDeg(a);
        viz.suggestFOV(zoom);
        viz.generate(world, buffer);
        std::string prefix = TRACK ? "Track" : "Dose";
        prefix += "Upper";
        std::string name = prefix + std::to_string(int(a)) + ".png";
        viz.savePNG(name, buffer);
        std::cout << "Rendertime " << buffer.renderTime.count() << " ms"
                  << "(" << 1000.0 / buffer.renderTime.count() << " fps)" << std::endl;
    }

    viz.setAzimuthalAngleDeg(120);
    for (auto a : angles) {
        viz.setPolarAngleDeg(a);
        viz.suggestFOV(zoom);
        viz.generate(world, buffer);
        std::string prefix = TRACK ? "Track" : "Dose";
        prefix += "Lower";
        std::string name = prefix + std::to_string(int(a)) + ".png";
        viz.savePNG(name, buffer);
        std::cout << "Rendertime " << buffer.renderTime.count() << " ms"
                  << "(" << 1000.0 / buffer.renderTime.count() << " fps)" << std::endl;
    }
}

int main()
{
    vizualize<true>();
    vizualize<false>();
    return 0;
}
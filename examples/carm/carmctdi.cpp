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
#include "dxmc/transport.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

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

void carmScatter()
{
    using CTDIPhantom = dxmc::CTDIPhantom<5, 1>;
    using Mesh = dxmc::TriangulatedMesh<5, 1>;
    using Sphere = dxmc::WorldSphere<5, 1>;
    using Room = dxmc::EnclosedRoom<5, 1>;
    using Box = dxmc::WorldBox<5, 1>;
    using Sphere = dxmc::WorldSphere<5, 1>;
    using World = dxmc::World<Mesh, Sphere, CTDIPhantom, Room, Box, Sphere>;
    using Viz = dxmc::VisualizeWorld;

    World world {};
    world.reserveNumberOfItems(6);

    // Adding c-arm
    auto& carm = world.addItem<Mesh>({ "carm.stl" });

    // Adding table
    auto& table = world.addItem<Box>();
    table.translate({ -30, 0, 0 });

    // Adding room
    auto& room = world.addItem<Room>();
    room.setInnerRoomAABB({ -350, -300, -150, 350, 300, 150 });
    room.setWallThickness(2);
    const auto lead = dxmc::Material<double, 5>::byZ(82).value();
    const auto lead_atom = dxmc::AtomHandler<double>::Atom(82);
    const auto lead_dens = dxmc::AtomHandler<double>::Atom(82).standardDensity;
    room.setMaterial(lead, lead_dens * 0.2 / 2.0);

    // Adding phantom
    auto& phantom = world.addItem<CTDIPhantom>();
    auto table_aabb = table.AABB();
    auto phantom_aabb = phantom.AABB();
    phantom.translate({ -40, 0, table_aabb[5] - phantom_aabb[2] });

    // Adding mesurements spheres

    world.build();

    // adding beam
    using Beam = dxmc::DXBeam;
    const std::array<double, 3> source_pos = { 0, 0, -70 };
    Beam beam(source_pos);
    beam.setBeamSize(6, 6, 114);
    beam.setNumberOfExposures(2000);
    beam.setNumberOfParticlesPerExposure(1000000);
    beam.setDAPvalue(25);

    dxmc::Transport transport;
    runDispatcher(transport, world, beam);
}

int main()
{

    return EXIT_SUCCESS;
}
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

Copyright 2022 Erlend Andersen
*/

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"
#include "dxmc/world/worlditems/fluencescore.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/triangulatedopensurface.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldboxgrid.hpp"
#include "dxmc/world/worlditems/worldcylinder.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

bool testWorld()
{
    using Sphere = dxmc::WorldSphere<5, 2>;

    using World = dxmc::World<Sphere>;

    World world;
    world.reserveNumberOfItems(2);
    auto& sphere = world.addItem<Sphere>({ 16, { 0, 0, 0 } });
    auto& sphere2 = world.addItem<Sphere>({ 16, { 0, 0, 0 } });
    sphere.translate({ 0, 0, 32 });
    sphere.setRadius(16);

    world.build();

    dxmc::PencilBeam beam;
    beam.setPosition({ -100, 0, 0 });
    beam.setDirection({ 1, 0, 0 });
    beam.setNumberOfExposures(16);
    beam.setNumberOfParticlesPerExposure(10000);

    dxmc::Transport::run(world, beam);

    dxmc::Particle p = { .pos = { -100, 0, 0 }, .dir = { 1, 0, 0 } };
    auto intersect = world.intersect(p);

    dxmc::VisualizeWorld viz(world);
    viz.setDistance(150);
    viz.setPolarAngleDeg(180);
    viz.setAzimuthalAngleDeg(0);
    viz.suggestFOV();
    viz.setCameraPosition({ 100, 100, 100 });

    auto buffer = viz.createBuffer(2048, 2048);
    viz.generate(world, buffer);
    viz.savePNG("test.png", buffer);

    bool test = true;
    return test;
}

int main(int argc, char* argv[])
{
 
    auto success = true;

    success = success && testWorld();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

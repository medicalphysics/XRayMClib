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

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <iostream>

bool testWorld()
{
    using Sphere = xraymc::WorldSphere<5, 2>;

    using World = xraymc::World<Sphere>;

    World world;
    world.reserveNumberOfItems(2);
    auto& sphere = world.addItem<Sphere>({ 16, { 0, 0, 0 } });
    auto& sphere2 = world.addItem<Sphere>({ 16, { 0, 0, 0 } });
    sphere.translate({ 0, 0, 32 });
    sphere.setRadius(16);

    world.build();

    xraymc::PencilBeam beam;
    beam.setPosition({ -100, 0, 0 });
    beam.setDirection({ 1, 0, 0 });
    beam.setNumberOfExposures(16);
    beam.setNumberOfParticlesPerExposure(10000);

    xraymc::Transport::run(world, beam);

    xraymc::Particle p = { .pos = { -100, 0, 0 }, .dir = { 1, 0, 0 } };
    auto intersect = world.intersect(p);

    xraymc::VisualizeWorld viz(world);
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

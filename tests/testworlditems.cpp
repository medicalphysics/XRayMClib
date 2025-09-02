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

Copyright 2023 Erlend Andersen
*/

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
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

template <typename U>
bool testItem()
{
    xraymc::World<U> world;
    world.reserveNumberOfItems(1);
    auto& item = world.template addItem<U>({});

    if constexpr (std::is_same_v<U, xraymc::TriangulatedMesh<>> || std::is_same_v<U, xraymc::TriangulatedOpenSurface<>>) {
        std::vector<std::array<double, 3>> p;
        constexpr double d = 30;
        p.push_back({ 1, 1, 0 }); // 0
        p.push_back({ 1, -1, 0 }); // 1
        p.push_back({ -1, -1, 0 }); // 2
        p.push_back({ -1, 1, 0 }); // 3
        p.push_back({ 0, 0, 1 });
        for (auto& i : p)
            for (auto& j : i)
                j *= d;

        std::vector<xraymc::Triangle> t;
        t.push_back({ p[0], p[1], p[4] });
        t.push_back({ p[1], p[2], p[4] });
        t.push_back({ p[2], p[3], p[4] });
        t.push_back({ p[3], p[0], p[4] });

        if constexpr (std::is_same_v<U, xraymc::TriangulatedMesh<>>) {
            // underside
            t.push_back({ p[0], p[3], p[2] });
            t.push_back({ p[2], p[1], p[0] });
        }
        item.setData(t);

    } else if constexpr (std::is_same_v<U, xraymc::TetrahedalMesh<>>) {
        std::vector<std::array<double, 3>> v(8);
        v[0] = { -1, 1, 1 };
        v[1] = { 1, -1, 1 };
        v[2] = { 1, 1, 1 };
        v[3] = { -1, -1, -1 };
        v[4] = { 1, -1, -1 };
        v[5] = { -1, -1, 1 };
        v[6] = { -1, 1, -1 };
        v[7] = { 1, 1, -1 };

        for (auto& i : v)
            for (auto& n : i)
                n *= 10;

        std::vector<xraymc::Tetrahedron> t(6);
        t[0] = { v[1], v[7], v[0], v[2], 0, 0 }; //*
        t[1] = { v[7], v[3], v[0], v[6], 0, 0 }; //*
        t[2] = { v[1], v[3], v[0], v[4], 0, 0 }; //*
        t[3] = { v[1], v[7], v[4], v[0], 0, 0 }; //*
        t[4] = { v[7], v[3], v[4], v[0], 0, 0 }; //*
        t[5] = { v[1], v[3], v[5], v[0], 0, 0 }; //*

        std::vector<double> dens(1, 1);
        std::vector<xraymc::Material<>> mats;
        mats.push_back(xraymc::Material<>::byChemicalFormula("H2O").value());

        item.setData(t, dens, mats);
    }

    item.translate({ 1, 1, 1 });
    item.translate({ -1, -1, -1 });
    auto center = item.center();
    auto aabb = item.AABB();
    xraymc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 } };
    auto intersection = item.intersect(p);
    auto intersectionViz = item.template intersectVisualization<U>(p);

    world.build();
    xraymc::PencilBeam<> beam;
    beam.setNumberOfExposures(1);
    beam.setNumberOfParticlesPerExposure(8);
    xraymc::Transport transport;
    transport.setNumberOfThreads(1);
    transport(world, beam);

    auto energy = item.energyScored();
    item.addEnergyScoredToDoseScore();
    item.clearEnergyScored();
    auto dose = item.doseScored(0);
    item.clearDoseScored();
    return true;
}

bool basicTestAllItems()
{
    auto success = true;
    success = success && testItem<xraymc::AAVoxelGrid<>>();
    success = success && testItem<xraymc::CTDIPhantom<>>();
    success = success && testItem<xraymc::DepthDose<>>();
    success = success && testItem<xraymc::FluenceScore>();
    success = success && testItem<xraymc::TetrahedalMesh<>>();
    success = success && testItem<xraymc::TriangulatedMesh<>>();
    success = success && testItem<xraymc::TriangulatedOpenSurface<>>();
    success = success && testItem<xraymc::WorldBox<>>();
    success = success && testItem<xraymc::WorldBoxGrid<>>();
    success = success && testItem<xraymc::WorldCylinder<>>();
    success = success && testItem<xraymc::WorldSphere<>>();
    success = success && testItem<xraymc::EnclosedRoom<>>();

    return success;
}

int main(int argc, char* argv[])
{
    std::cout << "Test world items ";
    auto success = true;

    success = success && basicTestAllItems();

    if (success)
        std::cout << "SUCCESS\n";
    else
        std::cout << "FAILURE\n";

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

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

        auto data = xraymc::TetrahedalMeshData {};
        data.nodes.resize(8);

        data.nodes[0] = { -1, 1, 1 };
        data.nodes[1] = { 1, -1, 1 };
        data.nodes[2] = { 1, 1, 1 };
        data.nodes[3] = { -1, -1, -1 };
        data.nodes[4] = { 1, -1, -1 };
        data.nodes[5] = { -1, -1, 1 };
        data.nodes[6] = { -1, 1, -1 };
        data.nodes[7] = { 1, 1, -1 };

        auto& v = data.nodes;
        for (auto& i : v)
            for (auto& n : i)
                n *= 10;

        auto& t = data.elements;
        t.resize(6);
        t[0] = { 1, 7, 0, 2 }; //*
        t[1] = { 7, 3, 0, 6 }; //*
        t[2] = { 1, 3, 0, 4 }; //*
        t[3] = { 1, 7, 4, 0 }; //*
        t[4] = { 7, 3, 4, 0 }; //*
        t[5] = { 1, 3, 5, 0 }; //*

        data.collectionIndices.resize(6);
        for (auto& c : data.collectionIndices)
            c = 1;

        data.collectionDensities.resize(1);
        data.collectionDensities[0] = 1;
        data.collectionMaterialComposition.resize(1);
        data.collectionMaterialComposition[0] = xraymc::Material<>::parseCompoundStr("H2O");

        std::vector<xraymc::Material<>> mats;
        mats.push_back(xraymc::Material<>::byChemicalFormula("H2O").value());

        data.collectionIndices.resize(6);

        item.setData(data);
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

template <typename U>
bool longTest()
{
    if constexpr (std::is_base_of_v<xraymc::DepthDose<>, U> == true) {

        xraymc::World<U> world;
        auto& item = world.addItem<U>();
        item.setRadius(2);
        item.setDirection({ 0, 0, -1 });
        item.setLenght(10);
        item.setResolution(100);
        auto aluminum_material = xraymc::Material<>::byZ(13).value();
        auto aluminum_density = xraymc::AtomHandler::Atom(13).standardDensity;
        item.setMaterial(aluminum_material, aluminum_density);

        world.build(1.0);

        xraymc::PencilBeam<> beam;
        beam.setPosition({ 0, 0, 10 });
        beam.setDirection({ 0, 0, -1 });
        beam.setEnergy(70); // keV
        beam.setAirKerma(1.0);
        beam.setNumberOfExposures(100);
        beam.setNumberOfParticlesPerExposure(1E6);
        xraymc::Transport transport;
        transport.runConsole(world, beam);

        double max_dose = 0;
        std::cout << "Depth [cm], Dose [keV], NumberOfEvents, Relative uncertanty [%]\n";
        for (const auto [depth, dose] : item.depthDoseScored()) {
            std::cout << depth << ", " << dose.dose();
            std::cout << ", " << dose.numberOfEvents();
            std::cout << ", " << dose.relativeUncertainty() << std::endl;
            max_dose = std::max(max_dose, dose.dose());
        }

        xraymc::VisualizeWorld viz(world);
        viz.setDistance(300);
        viz.setAzimuthalAngleDeg(60);
        viz.suggestFOV(1);
        auto buffer = viz.createBuffer(1024, 1024);
        viz.addLineProp(beam.position(), beam.direction(), 10, 0.1);
        viz.generate(world, buffer);
        viz.savePNG("pencilbeam.png", buffer);

        viz.addColorByValueItem(world.getItemPointers()[0]);
        viz.setColorByValueMinMax(0, max_dose);
        viz.generate(world, buffer);
        viz.savePNG("pencilbeam_color.png", buffer);
        return true;
    }
    return false;
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

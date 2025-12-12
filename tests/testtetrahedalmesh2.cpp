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

#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"

#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"

#include <iostream>
#include <string>
#include <vector>

using index_t = std::uint32_t;

xraymc::TetrahedalMeshData tetCube()
{
    std::vector<std::array<double, 3>> v(16);
    v[0] = { -1, 1, 1 };
    v[1] = { 1, -1, 1 };
    v[2] = { 1, 1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, -1, 1 };
    v[6] = { -1, 1, -1 };
    v[7] = { 1, 1, -1 };

    v[8] = { -1, 1, 3 };
    v[9] = { 1, -1, 3 };
    v[10] = { 1, 1, 3 };
    v[11] = { -1, -1, 1 };
    v[12] = { 1, -1, 1 };
    v[13] = { -1, -1, 3 };
    v[14] = { -1, 1, 1 };
    v[15] = { 1, 1, 1 };

    for (auto& i : v)
        for (auto& n : i)
            n *= 1;

    std::vector<std::array<index_t, 4>> t(12);

    t[0] = { 1, 7, 0, 2 }; //*
    t[1] = { 7, 3, 0, 6 }; //*
    t[2] = { 1, 3, 0, 4 }; //*
    t[3] = { 1, 7, 4, 0 }; //*
    t[4] = { 7, 3, 4, 0 }; //*
    t[5] = { 1, 3, 5, 0 }; //*

    t[6] = { 1 + 8, 7 + 8, 0 + 8, 2 + 8 }; //*
    t[7] = { 7 + 8, 3 + 8, 0 + 8, 6 + 8 }; //*
    t[8] = { 1 + 8, 3 + 8, 0 + 8, 4 + 8 }; //*
    t[9] = { 1 + 8, 7 + 8, 4 + 8, 0 + 8 }; //*
    t[10] = { 7 + 8, 3 + 8, 4 + 8, 0 + 8 }; //*
    t[11] = { 1 + 8, 3 + 8, 5 + 8, 0 + 8 }; //*

    xraymc::TetrahedalMeshData data;
    data.nodes = v;
    data.elements = t;
    data.collectionIndices.resize(t.size(), 0);
    return data;
}

bool testTetrahedalMeshGrid()
{

    auto data = tetCube();

    // xraymc::TetrahedalMeshGrid2 grid(data);
    return false;
}

void testWalk()
{
    std::string material_file = "/home/erlend/mrcptest/MRCP_AF_media.dat";
    std::string organ_file = "/home/erlend/mrcptest/icrp145organs.csv";
    std::string node_file = "/home/erlend/mrcptest/MRCP_AF.node";
    std::string element_file = "/home/erlend/mrcptest/MRCP_AF.ele";

    // std::string material_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF_media.dat";
    // std::string organ_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\icrp145organs.csv";
    // std::string node_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.node";
    // std::string element_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.ele";

    xraymc::TetrahedalMeshReader testreader(node_file, element_file, material_file, organ_file);
    testreader.rotate({ 0, 0, 1 }, std::numbers::pi_v<double>);
    using Mesh = xraymc::TetrahedalMesh<5, 2, true>;
    xraymc::World<Mesh> world;
    auto& item = world.template addItem<Mesh>(testreader.data());
    const auto aabb = item.AABB();
    item.translate({ -100, -30 - aabb[4] - 4, -aabb[2] - 120 });

    xraymc::Particle p1;
    p1.pos = { -95.259739884039476, -46.999651240345457, -118.48122907673229 };
    p1.dir = { -0.023143632090657524, 0.99961050744105795, 0.015595053930109647 };

    xraymc::Particle p2;
    p1.pos = { -86.73504383202426, -48.192271830600511, -5.0936455819747941 };
    p1.dir = { 0.021716261745313829, -0.46651056408553976, 0.88424900201945422 };
    xraymc::RandomState state;
    item.transport(p1, state);
}

void showPhantom()
{
    // std::string node_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.node";
    // std::string element_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.ele";
    // xraymc::TetrahedalMeshReader2 testreader(node_file, element_file);

    std::string material_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF_media.dat";
    std::string organ_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\icrp145organs.csv";
    std::string node_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.node";
    std::string element_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.ele";

    // std::string material_file = "/home/erlend/mrcptest/MRCP_AF_media.dat";
    // std::string organ_file = "/home/erlend/mrcptest/icrp145organs.csv";
    // std::string node_file = "/home/erlend/mrcptest/MRCP_AF.node";
    // std::string element_file = "/home/erlend/mrcptest/MRCP_AF.ele";

    xraymc::TetrahedalMeshReader testreader(node_file, element_file, material_file, organ_file);

    using Mesh = xraymc::TetrahedalMesh<5, 2, true>;
    xraymc::World<Mesh> world;

    auto tetdata = testreader.data();
    
    // tetdata.collectionNameMustContainFilter("Skin legs sensitive 50-100um", true);
    //tetdata.collectionNameMustContainFilter("blood", false);
    
    auto& item = world.template addItem<Mesh>(tetdata);
    item.setDisplayCollectionIndexFilter(79);
    // item.translate({ 10, 10, 10 });
    // item.rotate({ 0, 0, 1 }, 3.14 / 4);
    world.build();

    xraymc::VisualizeWorld viz(world);

    viz.setAzimuthalAngleDeg(80);
    viz.setPolarAngleDeg(0);
    viz.setDistance(1000);
    viz.suggestFOV(1);
    auto buffer = viz.template createBuffer<double>(1024, 1024);

    std::cout << "Start generating images" << std::endl;
    viz.generate(world, buffer);
    viz.savePNG("test0.png", buffer);
    std::cout << "Done 1" << std::endl;
    viz.setPolarAngleDeg(90);
    viz.generate(world, buffer);
    viz.savePNG("test90.png", buffer);
    std::cout << "Done 2" << std::endl;

    viz.setPolarAngleDeg(180);
    viz.generate(world, buffer);
    viz.savePNG("test180.png", buffer);
    std::cout << "Done 3\n";
    viz.setPolarAngleDeg(270);
    viz.generate(world, buffer);
    viz.savePNG("test270.png", buffer);
    std::cout << "Done 4\n";
}

void testTiming()
{
    // std::string material_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF_media.dat";
    // std::string organ_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\icrp145organs.csv";
    // std::string node_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.node";
    // std::string element_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.ele";

    std::string material_file = "/home/erlend/mrcptest/MRCP_AF_media.dat";
    std::string organ_file = "/home/erlend/mrcptest/icrp145organs.csv";
    std::string node_file = "/home/erlend/mrcptest/MRCP_AF.node";
    std::string element_file = "/home/erlend/mrcptest/MRCP_AF.ele";

    xraymc::TetrahedalMeshReader testreader(node_file, element_file, material_file, organ_file);

    // using Mesh = xraymc::TetrahedalMesh2<5, 2, false>;
    using Mesh = xraymc::TetrahedalMesh<5, 2, true>;
    xraymc::World<Mesh> world;
    auto& item = world.template addItem<Mesh>(testreader.data());
    world.build(1000);

    xraymc::DXBeam beam;
    beam.setPosition({ 0, -1800, 40 });
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 0, -1 });
    beam.setTubeVoltage(80);
    beam.setNumberOfExposures(10);
    beam.setNumberOfParticlesPerExposure(1000000);
    beam.setCollimationHalfAnglesDeg(0.3, 0.3);

    xraymc::Transport transport;
    const auto time = transport.runConsole(world, beam).count();

    std::cout << "Total time " << time << std::endl;

    xraymc::VisualizeWorld viz(world);

    std::vector<double> doses;
    doses.reserve(item.outerContourTetrahedronIndices().size());
    for (auto i : item.outerContourTetrahedronIndices()) {
        doses.push_back(item.doseScored(i).dose());
    }
    const auto max_dose = *std::max_element(doses.cbegin(), doses.cend());

    viz.setAzimuthalAngleDeg(80);
    viz.setDistance(1000);
    viz.suggestFOV(10);
    auto buffer = viz.template createBuffer<double>(1024, 1024);
    viz.addLineProp(beam, 180, .2);

    const int degstep = 45;
    for (int i = 0; i < 360; i = i + degstep) {
        viz.setPolarAngleDeg(i);
        viz.generate(world, buffer);
        std::string name = "test";
        name += std::to_string(i);
        name += ".png";
        viz.savePNG(name, buffer);
    }
    viz.addColorByValueItem(world.getItemPointers()[0]);
    viz.setColorByValueMinMax(0.0, max_dose);
    for (int i = 0; i < 360; i = i + degstep) {
        viz.setPolarAngleDeg(i);
        viz.generate(world, buffer);
        std::string name = "dose";
        name += std::to_string(i);
        name += ".png";
        viz.savePNG(name, buffer);
    }

    auto data = item.collectionData();
    std::cout << "Name, volume, density, dose, std, nevents\n";
    for (auto& d : data) {
        std::cout << d.name << ", ";
        std::cout << d.volume << ", ";
        std::cout << d.density << ", ";
        std::cout << d.dose << ", ";
        std::cout << std::sqrt(d.doseVariance) << ", ";
        std::cout << d.numberOfEvents << std::endl;
    }
}

int main()
{
    // testWalk();
    showPhantom();
    // testTiming();

    /*
        std::string material_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF_media.dat";
        std::string organ_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\icrp145organs.csv";
        std::string node_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.node";
        std::string element_file = "C:\\Users\\ander\\OneDrive\\phantomsMNCP\\adult\\MRCP_AF\\MRCP_AF.ele";

        // xraymc::TetrahedalMeshReader2 testreader(node_file, element_file, material_file, organ_file);



        node_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.node";
        element_file = "C:\\Users\\ander\\OneDrive\\tetgentest\\torus.1.ele";
        xraymc::TetrahedalMeshReader2 testreader(node_file, element_file);
        xraymc::TetrahedalMeshGrid2 grid(testreader.data(), { 32, 32, 16 });

        auto test = grid.pointInside({ -1, 0, 0 });

        xraymc::Particle p { .pos = { 0, -100, 0 }, .dir = { 0, 1, 0 } };
        auto res = grid.intersect(p);
    */
    bool success = true;

    // auto data = tetCube();
    // success = success && testTetrahedalMeshGrid();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
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

#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh2.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"

#include <iostream>
#include <string>
#include <vector>

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

    std::vector<std::array<std::uint32_t, 4>> t(12);

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
    data.collectionMaterialComposition.resize(1);
    data.collectionMaterialComposition[0][1] = 0.111894;
    data.collectionMaterialComposition[0][8] = 0.888106;
    data.collectionDensities.resize(1, 1.0);
    data.makeGenericCollectionNames();
    return data;
}
template <std::size_t N = 5, int L = 2, bool Fluence = true>
xraymc::TetrahedalMesh<N, L, Fluence> simpletetrahedron()
{
    auto tets = tetCube();

    xraymc::TetrahedalMesh<N, L, Fluence> mesh(tets, { 8, 8, 8 });

    return mesh;
}

template <std::size_t N = 5, int L = 2, bool F = false>
xraymc::WorldBox<N, L, F> simplebox()
{
    xraymc::WorldBox<N, L, F> box({ -1, -1, -1, 1, 1, 3 });
    auto water = xraymc::Material<N>::byNistName("Water, Liquid").value();
    box.setMaterial(water, 1);
    return box;
}

template <std::size_t N = 5, int L = 2>
xraymc::AAVoxelGrid<N, L, 255> simplegrid()
{

    xraymc::AAVoxelGrid<N, L, 255> grid;
    std::vector<std::uint8_t> m(1000, 0);
    std::vector<double> d(1000, 1);
    std::vector<xraymc::Material<N>> mats;
    mats.push_back(xraymc::Material<N>::byNistName("Water, Liquid").value());
    grid.setData({ 10, 10, 10 }, d, m, mats);
    grid.setSpacing({ .2, .2, .4 });
    grid.translate({ 0, 0, 1 });
    return grid;
}

bool testTransport()
{

    using M1 = xraymc::TetrahedalMesh<5, 1, true>;
    using M2 = xraymc::TetrahedalMesh<5, 1, false>;
    using T1 = xraymc::TetrahedalMesh2<5, 1, true>;
    using T2 = xraymc::TetrahedalMesh2<5, 1, false>;

    using B = xraymc::WorldBox<5, 1, false>;
    using BF = xraymc::WorldBox<5, 1, true>;
    using G = xraymc::AAVoxelGrid<5, 1, 255>;

    constexpr std::size_t N_HIST = 200000;
    constexpr std::size_t N_EXP = 32;

    xraymc::PencilBeam<> beam({ .050, .050, -1000 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(N_EXP);
    beam.setNumberOfParticlesPerExposure(N_HIST);

    xraymc::TransportProgress progress;

    {
        xraymc::World<M1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M1 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh forced " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<M2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M2 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh random " << dose << " " << progress.humanTotalTime() << std::endl;
    }
    {
        xraymc::World<T1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(T1 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh2 forced " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<T2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(T2 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh2 random " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<G> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(simplegrid<5, 1>());
        w.build();

        xraymc::Transport transport;
        transport(w, beam, &progress);
        auto s = mesh.size();
        double dose = 0;
        for (std::size_t i = 0; i < s; ++i)
            dose += mesh.doseScored(i).dose();
        std::cout << "Grid random " << dose / s << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<B> w;
        w.reserveNumberOfItems(1);
        auto& box = w.addItem(simplebox<5, 1, false>());
        w.build();

        xraymc::Transport transport;
        transport(w, beam, &progress);
        auto dose = box.doseScored().dose();
        std::cout << "Box random " << dose << " " << progress.humanTotalTime() << std::endl;
    }
    {
        xraymc::World<BF> w;
        w.reserveNumberOfItems(1);
        auto& box = w.addItem(simplebox<5, 1, true>());
        w.build();

        xraymc::Transport transport;
        transport(w, beam, &progress);
        auto dose = box.doseScored().dose();
        std::cout << "Box forced " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    return false;
}

xraymc::TetrahedalMeshReader readICRPPregnantPhantom(bool female = false)
{
    const std::string name = female ? "38wF" : "38wM";
    const std::string elefile = name + ".ele";
    const std::string nodefile = name + ".node";
    const std::string materialfile = name + ".material";

    const auto elepath = elefile;
    const auto nodepath = nodefile;
    const auto materialpath = materialfile;

    xraymc::TetrahedalMeshReader reader(nodepath, elepath, materialpath);
    return reader;
}

xraymc::TetrahedalMeshData getICRP()
{
    const std::string n = "C:/Users/ander/OneDrive/phantomsMNCP/adult/MRCP_AF/MRCP_AF.node";
    const std::string e = "C:/Users/ander/OneDrive/phantomsMNCP/adult/MRCP_AF/MRCP_AF.ele";

    xraymc::TetrahedalMeshReader reader(n, e);
    auto data = reader.data();
    data.collectionIndices.resize(data.elements.size(), 0);
    data.collectionMaterialComposition.resize(1);
    data.collectionMaterialComposition[0][1] = 0.111894;
    data.collectionMaterialComposition[0][8] = 0.888106;
    data.collectionDensities.resize(1, 1.0);
    data.makeGenericCollectionNames();
    return data;
}

bool testTransportICRP()
{

    using M1 = xraymc::TetrahedalMesh<5, 1, true>;
    using M2 = xraymc::TetrahedalMesh<5, 1, false>;
    using T1 = xraymc::TetrahedalMesh2<5, 1, true>;
    using T2 = xraymc::TetrahedalMesh2<5, 1, false>;

    constexpr std::size_t N_HIST = 10000;
    constexpr std::size_t N_EXP = 32;

    xraymc::PencilBeam<> beam({ .050, .050, -1000 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(N_EXP);
    beam.setNumberOfParticlesPerExposure(N_HIST);

    xraymc::TransportProgress progress;

    const auto data = getICRP();

    {
        xraymc::World<M1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M1 { data });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh forced " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<M2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M2 { data });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh random " << dose << " " << progress.humanTotalTime() << std::endl;
    }
    {
        xraymc::World<T1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(T1 { data });
        w.build();

        double dose = 0;
        xraymc::Transport transport;

        transport(w, beam, &progress);
        for (std::size_t i = 0; i < mesh.numberOfTetrahedra(); ++i)
            dose += mesh.doseScored(i).dose();
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh2 forced " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        xraymc::World<T2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(T2 { data });
        w.build();

        double dose = 0;
        xraymc::Transport transport;
        // transport.setNumberOfThreads(1);
        transport(w, beam, &progress);
        for (std::size_t i = 0; i < mesh.numberOfTetrahedra(); ++i)
            dose += mesh.doseScored(i).dose();
        std::cout << "Mesh2 random " << dose << " " << progress.humanTotalTime() << std::endl;
    }

    return false;
}

int main()
{

    auto reader = readICRPPregnantPhantom();
    auto& data = reader.data();
    data.collectionNameMustContainFilter("Fetal");

    bool success = true;

    success = success && testTransport();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
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
#include "xraymc/world/worlditems/worldbox.hpp"

#include <iostream>
#include <string>
#include <vector>

xraymc::TetrahedalMeshData tetCube()
{
    std::vector<std::array<double, 3>> v(10);

    v[0] = { 1, -1, 1 };
    v[1] = { -1, 1, -1 };
    v[2] = { -1, -1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, 1, 1 };
    v[6] = { 1, 1, 1 };
    v[7] = { 1, 1, -1 };
    v[8] = { 0, 0, 1 };
    v[9] = { 0, 0, -1 };

    std::vector<std::array<std::uint32_t, 4>> t(12);

    t[0] = { 0, 7, 4, 9 };
    t[1] = { 1, 6, 5, 8 };
    t[2] = { 1, 6, 8, 9 };
    t[3] = { 6, 1, 7, 9 };
    t[4] = { 1, 5, 2, 8 };
    t[5] = { 3, 0, 4, 9 };
    t[6] = { 3, 1, 2, 9 };
    t[7] = { 2, 1, 8, 9 };
    t[8] = { 7, 0, 6, 9 };
    t[9] = { 6, 0, 8, 9 };
    t[10] = { 0, 3, 2, 9 };
    t[11] = { 0, 2, 8, 9 };

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

template <std::size_t N = 5, int L = 2, bool F = false>
xraymc::WorldBox<N, L, F> simplebox()
{
    xraymc::WorldBox<N, L, F> box({ -1, -1, -1, 1, 1, 1 });
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
    grid.setSpacing({ .2, .2, .2 });

    return grid;
}

bool testTransport()
{
    constexpr int INTER = 1;

    using M1 = xraymc::TetrahedalMesh<12, INTER, true>;
    using M2 = xraymc::TetrahedalMesh<12, INTER, false>;

    using B = xraymc::WorldBox<12, INTER, false>;
    using BF = xraymc::WorldBox<12, INTER, true>;
    using G = xraymc::AAVoxelGrid<12, INTER, 255>;

    constexpr std::size_t N_THREADS = 0;

    constexpr std::size_t N_HIST = 1000000;
    constexpr std::size_t N_EXP = 32;

    xraymc::PencilBeam<> beam({ -10, .0001, .0001 }, { 1, 0, 0 }, 60);

    beam.setNumberOfExposures(N_EXP);
    beam.setNumberOfParticlesPerExposure(N_HIST);

    xraymc::TransportProgress progress;

    bool success = true;

    {
        xraymc::World<M1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M1 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;
        transport.setNumberOfThreads(N_THREADS);
        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh forced " << dose << " " << progress.humanTotalTime() << std::endl;
        success = success && dose > 16.0 && dose < 16.1;
    }

    {
        xraymc::World<M2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(M2 { tetCube() });
        w.build();

        double dose = 0;
        xraymc::Transport transport;
        transport.setNumberOfThreads(N_THREADS);
        transport(w, beam, &progress);
        auto coll = mesh.collectionData();
        dose = coll[0].dose;
        std::cout << "Mesh random " << dose << " " << progress.humanTotalTime() << std::endl;
        success = success && dose > 16.0 && dose < 16.1;
    }

    {
        xraymc::World<G> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(simplegrid<12, INTER>());
        w.build();

        xraymc::Transport transport;
        transport.setNumberOfThreads(N_THREADS);
        transport(w, beam, &progress);
        auto s = mesh.size();
        double dose = 0;
        for (std::size_t i = 0; i < s; ++i)
            dose += mesh.doseScored(i).dose();
        std::cout << "Grid random " << dose / s << " " << progress.humanTotalTime() << std::endl;
        success = success && dose / s > 16.0 && dose / s < 16.1;
    }

    {
        xraymc::World<B> w;
        w.reserveNumberOfItems(1);
        auto& box = w.addItem(simplebox<12, INTER, false>());
        w.build();

        xraymc::Transport transport;
        transport.setNumberOfThreads(N_THREADS);
        transport(w, beam, &progress);
        auto dose = box.doseScored().dose();
        std::cout << "Box random " << dose << " " << progress.humanTotalTime() << std::endl;
        success = success && dose > 16.0 && dose < 16.1;
    }
    {
        xraymc::World<BF> w;
        w.reserveNumberOfItems(1);
        auto& box = w.addItem(simplebox<12, INTER, true>());
        w.build();

        xraymc::Transport transport;
        transport.setNumberOfThreads(N_THREADS);
        transport(w, beam, &progress);
        auto dose = box.doseScored().dose();
        std::cout << "Box forced " << dose << " " << progress.humanTotalTime() << std::endl;
        success = success && dose > 16.0 && dose < 16.1;
    }

    return success;
}

int main()
{

    bool success = true;

    success = success && testTransport();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
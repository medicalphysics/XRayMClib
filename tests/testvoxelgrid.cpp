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

#include "testutils.hpp"

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using xraymctest::checkStatisticallyEqual;

using Grid = AAVoxelGrid<5, 1, 255>;

Grid makeGrid(const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing)
{
    const auto size = dim[0] * dim[1] * dim[2];
    std::vector<double> dens(size, 1.0);
    std::vector<std::uint8_t> mat(size, 0);
    std::vector<Material<5>> materials = { Material<5>::byNistName("Water, Liquid").value() };
    return Grid(dim, spacing, dens, mat, materials);
}

void testIndexing()
{
    const std::array<std::size_t, 3> dim = { 7, 5, 3 };
    const std::array<double, 3> spacing = { 0.2, 0.5, 1.0 };
    auto grid = makeGrid(dim, spacing);
    check(grid.size() == 7 * 5 * 3, "VoxelGrid: size");
    bool roundTrip = true;
    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const std::array<std::size_t, 3> ind = { x, y, z };
                const auto f = grid.flatIndex(ind);
                roundTrip = roundTrip && grid.index(f) == ind && f == x + dim[0] * (y + dim[1] * z);
            }
    check(roundTrip, "VoxelGrid: flat index round trip and x fastest ordering");

    // grid is centered at the origin
    const auto aabb = grid.AABB();
    checkNear(aabb[0], -0.7, 1E-12, "VoxelGrid: AABB x min");
    checkNear(aabb[4], 1.25, 1E-12, "VoxelGrid: AABB y max");
    checkNear(aabb[5], 1.5, 1E-12, "VoxelGrid: AABB z max");
    // positions map to voxel indices
    check(grid.index(std::array<double, 3> { -0.69, -1.24, -1.49 }) == std::array<std::size_t, 3> { 0, 0, 0 }, "VoxelGrid: position in first voxel");
    check(grid.index(std::array<double, 3> { 0.69, 1.24, 1.49 }) == std::array<std::size_t, 3> { 6, 4, 2 }, "VoxelGrid: position in last voxel");
    check(grid.index(std::array<double, 3> { 0.01, 0.01, 0.01 }) == std::array<std::size_t, 3> { 3, 2, 1 }, "VoxelGrid: position in center voxel");

    grid.translate({ 1, 2, 3 });
    check(grid.index(std::array<double, 3> { 1.01, 2.01, 3.01 }) == std::array<std::size_t, 3> { 3, 2, 1 }, "VoxelGrid: index follows translation");

    // setData rejects inconsistent input
    Grid g2;
    check(!g2.setData({ 2, 2, 2 }, std::vector<double>(7, 1.0), std::vector<std::uint8_t>(8, 0), { Material<5>::byZ(1).value() }), "VoxelGrid: setData rejects wrong density size");
    check(!g2.setData({ 2, 2, 2 }, std::vector<double>(8, 1.0), std::vector<std::uint8_t>(8, 1), { Material<5>::byZ(1).value() }), "VoxelGrid: setData rejects material index out of range");
}

void testIntersection()
{
    auto grid = makeGrid({ 64, 64, 64 }, { 0.2, 0.2, 0.2 });
    Particle p { .pos = { -100, 0, 0 }, .dir = { 1, 0, 0 }, .energy = 60, .weight = 1 };
    auto res = grid.intersect(p);
    check(res.valid() && !res.rayOriginIsInsideItem, "VoxelGrid: ray from outside hits");
    checkNear(res.intersection, 100 - 6.4, 1E-9, "VoxelGrid: distance to grid along x");
    p.pos = { 100, 0, 0 };
    p.dir = { -1, 0, 0 };
    checkNear(grid.intersect(p).intersection, 100 - 6.4, 1E-9, "VoxelGrid: distance to grid along -x");
    p.pos = { -100, -100, -100 };
    p.dir = vectormath::normalized(std::array<double, 3> { 1, 1, 1 });
    checkNear(grid.intersect(p).intersection, std::sqrt(3.0) * (100 - 6.4), 1E-9, "VoxelGrid: distance to grid along diagonal");
    p.pos = { 10, 0, -100 };
    p.dir = { 0, 0, 1 };
    check(!grid.intersect(p).valid(), "VoxelGrid: ray beside grid misses");
}

void testCTNumbers()
{
    std::vector<double> hu = { -1000, -700, -50, 0, 40, 300, 1000, 2000 };
    AAVoxelGrid<5, 1, 255> grid({ 8, 1, 1 }, { 1, 1, 1 }, hu);
    const auto dens = grid.getDensity();
    const auto mats = grid.getMaterialIndex();
    bool densIncreasing = true;
    for (std::size_t i = 1; i < dens.size(); ++i)
        densIncreasing = densIncreasing && dens[i] > dens[i - 1];
    check(densIncreasing, "VoxelGrid: density increases with CT number");
    checkNear(dens[3], 1.018, 1E-9, "VoxelGrid: density of water at 0 HU");
    checkNear(dens[1], 0.3, 1E-9, "VoxelGrid: density of lung at -700 HU");
    check(mats[0] == 0 && mats[1] == 1, "VoxelGrid: air and lung material classes");
    bool matIncreasing = true;
    for (std::size_t i = 1; i < mats.size(); ++i)
        matIncreasing = matIncreasing && mats[i] >= mats[i - 1];
    check(matIncreasing, "VoxelGrid: material class increases with CT number");

    // Schneider et al. 2000 tissue compositions
    const auto w = Grid::shcneiderMaterialWeights();
    const auto& air = w[0];
    check(air.contains(7) && air.contains(8) && air.contains(18) && !air.contains(1), "VoxelGrid: Schneider air consists of N, O and Ar");
    bool noArgon = true, organic = true, normalized = true;
    for (std::size_t i = 1; i < w.size(); ++i) {
        noArgon = noArgon && !w[i].contains(18);
        organic = organic && w[i].contains(1) && w[i].contains(6) && w[i].contains(7) && w[i].contains(8);
        const auto sum = std::transform_reduce(w[i].cbegin(), w[i].cend(), 0.0, std::plus<>(), [](const auto& p) { return p.second; });
        normalized = normalized && std::abs(sum - 100) < 1.0;
    }
    check(noArgon, "VoxelGrid: Schneider tissues contain no argon");
    check(organic, "VoxelGrid: Schneider tissues contain H, C, N and O");
    check(normalized, "VoxelGrid: Schneider tissue weights sum to 100 %");
    for (std::size_t i = 9; i < w.size(); ++i)
        check(w[i].contains(20) && (i == 9 || w[i].at(20) > w[i - 1].at(20)), "VoxelGrid: Skeletal tissue " + std::to_string(i - 8) + " calcium content increases");
}

// A voxel grid where half of the voxels are transparent equals a box covering the other half
void testTransparentVoxels()
{
    using TGrid = AAVoxelGrid<5, 1, 1>;
    using Box = WorldBox<5, 1>;
    const auto water = Material<5>::byNistName("Water, Liquid").value();
    const std::array<std::size_t, 3> dim = { 4, 4, 4 };
    std::vector<double> dens(64, 1.0);
    std::vector<std::uint8_t> mat(64, 1);
    for (std::size_t i = 0; i < 64; ++i) {
        const auto z = i / 16;
        if (z >= 2)
            mat[i] = 0; // upper half in z is water, lower half transparent
    }
    PencilBeam<> beam({ 0.1, 0.2, -10 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(40);
    beam.setNumberOfParticlesPerExposure(20000);
    Transport transport;

    World<TGrid> wg;
    auto& grid = wg.addItem<TGrid>({ dim, { 0.5, 0.5, 0.5 }, dens, mat, { water, water } });
    wg.setMaterialDensity(1E-9);
    wg.build();
    transport(wg, beam, nullptr, false);
    double eGrid = 0, vGrid = 0, eTransparent = 0;
    for (std::size_t i = 0; i < 64; ++i) {
        eGrid += grid.energyScored(i).energyImparted();
        vGrid += grid.energyScored(i).variance();
        if (mat[i] == 1)
            eTransparent += grid.energyScored(i).energyImparted();
    }
    check(eTransparent == 0, "VoxelGrid: no energy in transparent voxels");

    World<Box> wb;
    auto& box = wb.addItem<Box>({ { -1, -1, 0, 1, 1, 1 }, water, 1.0 });
    wb.setMaterialDensity(1E-9);
    wb.build();
    transport(wb, beam, nullptr, false);
    // voxels of the same history are correlated, increase the variance estimate
    checkStatisticallyEqual(eGrid, 2 * std::sqrt(vGrid), box.energyScored().energyImparted(), box.energyScored().standardDeviation(), 4, "VoxelGrid: transparent voxels do not attenuate");
}

void testDoseScoring()
{
    auto grid = makeGrid({ 3, 3, 3 }, { 1, 1, 1 });
    World<Grid> world;
    auto& g = world.addItem(std::move(grid));
    world.setMaterialDensity(1E-9);
    world.build();
    PencilBeam<> beam({ 0, 0, -10 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(10);
    beam.setNumberOfParticlesPerExposure(10000);
    Transport transport;
    transport(world, beam);
    // the beam passes through the central column of voxels
    const auto center0 = g.doseScored(g.flatIndex(std::array<std::size_t, 3> { 1, 1, 0 })).dose();
    const auto center2 = g.doseScored(g.flatIndex(std::array<std::size_t, 3> { 1, 1, 2 })).dose();
    const auto corner = g.doseScored(g.flatIndex(std::array<std::size_t, 3> { 0, 0, 0 })).dose();
    check(center0 > center2, "VoxelGrid: dose decreases with depth along beam");
    check(center2 > 10 * corner, "VoxelGrid: dose is highest along the beam");
    const auto scores = g.getDoseScores();
    check(scores.size() == 27 && scores[13].dose() > 0, "VoxelGrid: dose scores for all voxels");
}

int main()
{
    std::cout << "Testing voxel grid\n";
    testIndexing();
    testIntersection();
    testCTNumbers();
    testTransparentVoxels();
    testDoseScoring();
    return xraymctest::summary("testvoxelgrid");
}

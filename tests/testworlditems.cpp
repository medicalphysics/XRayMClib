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

// Generic tests that every world item must pass: geometry, intersections,
// transport, scoring and serialization.

#include "testgeometry.hpp"
#include "testutils.hpp"

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/personaldosimeter.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>
#include <string>

using xraymctest::check;
using Vec = std::array<double, 3>;
namespace vm = xraymc::vectormath;

// Item specific behaviour that differs from the generic contract
struct ItemOptions {
    std::size_t scoringIndices = 1; // number of energy/dose scoring indices to test
    bool scoresEnergy = true; // a beam aimed at the center deposits energy
    bool hitsInsideAABB = true; // intersection points lie within the item AABB
    bool clearsEnergy = true; // clearEnergyScored resets the energy score
};

template <typename U>
double totalEnergy(const U& item, std::size_t n)
{
    double e = 0;
    for (std::size_t i = 0; i < n; ++i)
        e += item.energyScored(i).energyImparted();
    return e;
}

template <typename U>
void testItem(U item, const std::string& name, ItemOptions opt = {})
{
    // Geometry
    const auto aabb0 = item.AABB();
    const auto c0 = item.center();
    check(xraymc::basicshape::AABB::pointInside(c0, aabb0), name + ": center inside AABB");
    bool aabbValid = true;
    for (std::size_t i = 0; i < 3; ++i)
        aabbValid = aabbValid && aabb0[i] <= aabb0[i + 3];
    check(aabbValid, name + ": AABB is well formed");

    const Vec shift = { 1.5, -2.25, 3.0 };
    item.translate(shift);
    const auto aabb1 = item.AABB();
    const auto c1 = item.center();
    bool shifted = true;
    for (std::size_t i = 0; i < 3; ++i) {
        shifted = shifted && std::abs(aabb1[i] - aabb0[i] - shift[i]) < 1E-9 && std::abs(aabb1[i + 3] - aabb0[i + 3] - shift[i]) < 1E-9;
        shifted = shifted && std::abs(c1[i] - c0[i] - shift[i]) < 1E-9;
    }
    check(shifted, name + ": translate moves AABB and center");
    item.translate(vm::scale(shift, -1.0));
    bool restored = true;
    for (std::size_t i = 0; i < 6; ++i)
        restored = restored && std::abs(item.AABB()[i] - aabb0[i]) < 1E-9;
    check(restored, name + ": translating back restores AABB");

    // Intersections: rays from outside the AABB aimed at the center hit inside the AABB,
    // rays pointing away from the AABB miss.
    auto state = xraymctest::seededState();
    const auto [lo, hi] = vm::splice(aabb0);
    const auto diag = vm::length(vm::subtract(hi, lo));
    bool hitsOk = true, missOk = true;
    std::size_t nHits = 0;
    for (std::size_t i = 0; i < 500; ++i) {
        const auto cost = state.randomUniform(-1.0, 1.0);
        const auto phi = state.randomUniform(2 * std::numbers::pi);
        const auto sint = std::sqrt(1 - cost * cost);
        const Vec u = { sint * std::cos(phi), sint * std::sin(phi), cost };
        const auto start = vm::add(c0, vm::scale(u, diag * 2));
        const xraymc::Particle p { .pos = start, .dir = vm::scale(u, -1.0), .energy = 60, .weight = 1 };
        const auto res = item.intersect(p);
        if (res.valid()) {
            ++nHits;
            const auto hit = vm::add(p.pos, vm::scale(p.dir, res.intersection));
            const std::array<double, 6> padded = { lo[0] - 1E-6, lo[1] - 1E-6, lo[2] - 1E-6, hi[0] + 1E-6, hi[1] + 1E-6, hi[2] + 1E-6 };
            hitsOk = hitsOk && res.intersection > 0 && !res.rayOriginIsInsideItem && xraymc::basicshape::AABB::pointInside(hit, padded);
        }
        const xraymc::Particle q { .pos = start, .dir = u, .energy = 60, .weight = 1 };
        missOk = missOk && !item.intersect(q).valid();
    }
    check(nHits > 0, name + ": rays aimed at center hit item (" + std::to_string(nHits) + " of 500)");
    if (opt.hitsInsideAABB)
        check(hitsOk, name + ": intersection points are on the item bounding box or inside it");
    check(missOk, name + ": rays pointing away from item miss");

    const xraymc::Particle pv { .pos = vm::add(c0, Vec { 0, 0, -diag * 2 }), .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto viz = item.template intersectVisualization<U>(pv);
    if (viz.valid())
        check(std::abs(vm::length(viz.normal) - 1) < 1E-6, name + ": visualization normal is unit vector");

    // Transport and scoring in a world
    xraymc::World<U> world;
    world.reserveNumberOfItems(1);
    auto& witem = world.addItem(std::move(item), "item");
    world.build();
    check(world.getItemPointerFromName("item") != nullptr && world.getItemPointerFromName("other") == nullptr, name + ": world item lookup by name");

    xraymc::PencilBeam<> beam(vm::add(c0, Vec { 0.01, 0.02, -diag * 2 }), { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(2000);
    xraymc::Transport transport;
    transport.setNumberOfThreads(2);
    transport(world, beam, nullptr, false);
    const auto energy = totalEnergy(witem, opt.scoringIndices);
    check(energy >= 0 && std::isfinite(energy), name + ": energy scored is finite and non-negative");
    if (opt.scoresEnergy)
        check(energy > 0, name + ": energy is scored for beam through item");
    check(energy <= 60.0 * beam.numberOfParticles(), name + ": energy scored does not exceed beam energy");

    witem.addEnergyScoredToDoseScore(1.0);
    bool doseOk = true;
    for (std::size_t i = 0; i < opt.scoringIndices; ++i)
        doseOk = doseOk && witem.doseScored(i).dose() >= 0 && std::isfinite(witem.doseScored(i).dose());
    check(doseOk, name + ": dose is finite and non-negative");

    // Serialization round trip preserves geometry and scored dose
    const auto buffer = witem.serialize();
    const auto copy = U::deserialize(buffer);
    check(copy.has_value(), name + ": deserialize succeeds");
    if (copy) {
        const auto& c = copy.value();
        bool aabbEqual = true;
        for (std::size_t i = 0; i < 6; ++i)
            aabbEqual = aabbEqual && std::abs(c.AABB()[i] - witem.AABB()[i]) < 1E-9;
        check(aabbEqual, name + ": deserialized AABB equals original");
        bool doseEqual = true;
        for (std::size_t i = 0; i < opt.scoringIndices; ++i)
            doseEqual = doseEqual && std::abs(c.doseScored(i).dose() - witem.doseScored(i).dose()) <= 1E-12 * std::abs(witem.doseScored(i).dose());
        check(doseEqual, name + ": deserialized dose equals original");
        bool interEqual = true;
        for (std::size_t i = 0; i < 200; ++i) {
            const Vec pos = { state.randomUniform(lo[0] - 1, hi[0] + 1), state.randomUniform(lo[1] - 1, hi[1] + 1), state.randomUniform(lo[2] - 1, hi[2] + 1) };
            auto dir = Vec { state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0) };
            vm::normalize(dir);
            const xraymc::Particle p { .pos = pos, .dir = dir, .energy = 60, .weight = 1 };
            const auto a = witem.intersect(p);
            const auto b = c.intersect(p);
            interEqual = interEqual && a.valid() == b.valid() && (!a.valid() || std::abs(a.intersection - b.intersection) < 1E-9);
        }
        check(interEqual, name + ": deserialized item has identical intersections");
        check(U::validMagicID(std::span<const char>(U::magicID().data(), 32)), name + ": magic ID is valid");
    }

    witem.clearEnergyScored();
    if (opt.clearsEnergy)
        check(totalEnergy(witem, opt.scoringIndices) == 0, name + ": clearEnergyScored resets energy");
    witem.clearDoseScored();
    bool cleared = true;
    for (std::size_t i = 0; i < opt.scoringIndices; ++i)
        cleared = cleared && witem.doseScored(i).dose() == 0;
    check(cleared, name + ": clearDoseScored resets dose");
}

int main()
{
    std::cout << "Testing world items\n";
    const auto water5 = xraymc::Material<5>::byNistName("Water, Liquid").value();
    {
        xraymc::WorldBox<5, 1> item({ -2, -1, -3, 2, 1, 3 });
        item.setMaterial(water5, 1.0);
        testItem(item, "WorldBox");
    }
    {
        xraymc::WorldBox<5, 1, true> item({ -2, -1, -3, 2, 1, 3 });
        item.setMaterial(water5, 1.0);
        testItem(item, "WorldBox (forced)");
    }
    {
        xraymc::WorldSphere<5, 1> item(3, { 1, 2, 3 });
        item.setMaterial(water5, 1.0);
        testItem(item, "WorldSphere");
    }
    {
        xraymc::WorldSphere<5, 1, true> item(3, { 1, 2, 3 });
        item.setMaterial(water5, 1.0);
        testItem(item, "WorldSphere (forced)");
    }
    {
        xraymc::WorldCylinder<5, 1> item(2, 6, { 0, 0, 0 }, { 0, 1, 1 });
        item.setMaterial(water5, 1.0);
        testItem(item, "WorldCylinder");
    }
    {
        xraymc::WorldBoxGrid<5, 1> item({ -2, -2, -2, 2, 2, 2 });
        item.setMaterial(water5, 1.0);
        item.setVoxelDimensions({ 4, 3, 2 });
        testItem(item, "WorldBoxGrid", { .scoringIndices = 24 });
    }
    {
        std::vector<double> dens(4 * 5 * 6, 1.0);
        std::vector<std::uint8_t> mat(dens.size(), 0);
        std::vector<xraymc::Material<5>> mats = { water5 };
        xraymc::AAVoxelGrid<5, 1, 255> item({ 4, 5, 6 }, { 0.5, 0.5, 0.5 }, dens, mat, mats);
        testItem(item, "AAVoxelGrid", { .scoringIndices = dens.size() });
    }
    {
        xraymc::CTDIPhantom<5, 1> item(8, 15);
        testItem(item, "CTDIPhantom", { .scoringIndices = 5 });
    }
    {
        xraymc::DepthDose<5, 1> item(2, 10, 20);
        item.setMaterial(water5, 1.0);
        testItem(item, "DepthDose", { .scoringIndices = 20 });
    }
    {
        xraymc::EnclosedRoom<5, 1> item(5, 50);
        testItem(item, "EnclosedRoom");
    }
    {
        xraymc::TriangulatedMesh<5, 1> item(xraymctest::boxTriangles(2));
        item.setMaterial(water5, 1.0);
        testItem(item, "TriangulatedMesh");
    }
    {
        xraymc::TriangulatedOpenSurface<5, 1> item(xraymctest::planeTriangles(2), 0.5);
        item.setMaterial(water5, 1.0);
        // The surface is a slab of finite thickness, the entry point of grazing rays is
        // moved back along the ray by thickness / (2 cos(theta)) and may lie outside the AABB.
        testItem(item, "TriangulatedOpenSurface", { .hitsInsideAABB = false });
    }
    {
        xraymc::TetrahedalMesh<5, 1, false> item(xraymctest::tetCube(2));
        testItem(item, "TetrahedalMesh", { .scoringIndices = 12 });
    }
    {
        xraymc::TetrahedalMesh<5, 1, true> item(xraymctest::tetCube(2));
        testItem(item, "TetrahedalMesh (forced)", { .scoringIndices = 12 });
    }
    {
        xraymc::FlatDetector item;
        item.setDetectorDimensions(8, 8);
        item.setPixelSpacing(0.5, 0.5);
        testItem(item, "FlatDetector", { .scoringIndices = 64 });
    }
    {
        xraymc::FluenceScore item(3);
        // The fluence counter accumulates until clearDoseScored is called
        testItem(item, "FluenceScore", { .clearsEnergy = false });
    }
    {
        xraymc::PersonalDosimeter<true> item;
        // The dosimeter faces +z and is shielded from behind, see testpersonaldosimeter for response tests
        testItem(item, "PersonalDosimeter", { .scoringIndices = 2, .scoresEnergy = false });
    }
    return xraymctest::summary("testworlditems");
}

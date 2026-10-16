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
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>
#include <limits>
#include <string>
#include <variant>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using Vec = std::array<double, 3>;
namespace vm = xraymc::vectormath;

using Sphere = WorldSphere<5, 1>;
using Box = WorldBox<5, 1>;
using Cylinder = WorldCylinder<5, 1>;
using W = World<Sphere, Box, Cylinder>;

void testItems()
{
    W world;
    world.reserveNumberOfItems(3);
    auto& s = world.addItem<Sphere>({ 2, { 0, 0, 0 } });
    auto& b = world.addItem<Box>({ { 5, -1, -1, 7, 1, 1 } }, "Box");
    world.addItem<Cylinder>("Cyl");
    check(world.items().size() == 3, "World: number of items");
    check(world.itemNames()[0] == "Item 1" && world.itemNames()[1] == "Box" && world.itemNames()[2] == "Cyl", "World: item names");
    check(std::holds_alternative<Box>(*world.getItemPointerFromName("Box")), "World: item by name has correct type");
    check(&std::get<Box>(*world.getItemPointerFromName("Box")) == &b, "World: item by name is the added item");
    check(world.getItemPointerFromName("nothing") == nullptr, "World: unknown item name gives nullptr");
    check(world.getItemPointers().size() == 3 && world.getItemPointers()[0] == &world.items()[0], "World: item pointers");
    (void)s;

    W w2;
    w2.reserveNumberOfItems(2);
    w2.addItem<Sphere>({ 1, { 0, 0, 0 } });
    w2.addItem<Box>({ { 5, -1, -1, 7, 1, 1 } });
    w2.build(3);
    const auto aabb = w2.AABB();
    const std::array<double, 6> expected = { -4, -4, -4, 10, 4, 4 };
    bool aabbOk = true;
    for (std::size_t i = 0; i < 6; ++i)
        aabbOk = aabbOk && std::abs(aabb[i] - expected[i]) < 1E-6;
    check(aabbOk, "World: AABB is union of items plus padding");
    check(xraymctest::checkVectorNear(w2.center(), { 3, 0, 0 }, 1E-6, "World: center"), "World: center of AABB");
    w2.build({ -50, -1, -1, 1, 1, 60 });
    // build(aabb) uses the default padding of 10 cm around the items
    check(std::abs(w2.AABB()[0] + 50) < 1E-9 && std::abs(w2.AABB()[5] - 60) < 1E-9 && std::abs(w2.AABB()[3] - 17) < 1E-6, "World: build with AABB extends world");

    w2.translate({ 0, 0, 100 });
    const auto& sp = std::get<Sphere>(w2.items()[0]);
    checkNear(sp.center()[2], 100, 1E-12, "World: translate moves items");
    const Particle p { .pos = { 0, 0, 50 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto r = w2.intersect(p);
    check(r.valid() && std::abs(r.intersection - 49) < 1E-9, "World: translate moves kd-tree");

    w2.setMaterial(Material<16>::byNistName("Water, Liquid").value(), 0.5);
    checkNear(w2.fillMaterialDensity(), 0.5, 0, "World: fill material density");
    check(w2.fillMaterial().composition().contains(8), "World: fill material");
    w2.setMaterialByWeight({ { 82, 1.0 } }, 2.0);
    check(w2.fillMaterial().composition().size() == 1 && w2.fillMaterialDensity() == 2.0, "World: fill material by weight");
}

// The kd-tree must return the nearest intersection among all items
void testKDTree()
{
    W world;
    world.reserveNumberOfItems(60);
    auto state = xraymctest::seededState();
    for (int i = 0; i < 60; ++i) {
        const Vec c = { state.randomUniform(-30.0, 30.0), state.randomUniform(-30.0, 30.0), state.randomUniform(-30.0, 30.0) };
        const auto r = state.randomUniform(0.5, 3.0);
        if (i % 3 == 0)
            world.addItem<Sphere>({ r, c });
        else if (i % 3 == 1)
            world.addItem<Box>({ { c[0] - r, c[1] - r, c[2] - r, c[0] + r, c[1] + 0.5 * r, c[2] + 2 * r } });
        else
            world.addItem<Cylinder>({ r, 2 * r, c, { state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0), 1 } });
    }
    world.build();

    std::size_t wrongItem = 0, wrongDistance = 0, hits = 0, falseMiss = 0;
    for (int i = 0; i < 20000; ++i) {
        const auto cost = state.randomUniform(-1.0, 1.0);
        const auto phi = state.randomUniform(2 * std::numbers::pi);
        const Vec u = { std::sqrt(1 - cost * cost) * std::cos(phi), std::sqrt(1 - cost * cost) * std::sin(phi), cost };
        const Vec target = { state.randomUniform(-30.0, 30.0), state.randomUniform(-30.0, 30.0), state.randomUniform(-30.0, 30.0) };
        const Particle p { .pos = vm::add(target, vm::scale(u, 80.0)), .dir = vm::scale(u, -1.0), .energy = 60, .weight = 1 };

        // brute force nearest intersection
        double tmin = std::numeric_limits<double>::max();
        const std::variant<Sphere, Box, Cylinder>* nearest = nullptr;
        for (const auto& item : world.items()) {
            const auto r = std::visit([&](const auto& it) { return it.intersect(p); }, item);
            if (r.valid() && r.intersection < tmin) {
                tmin = r.intersection;
                nearest = &item;
            }
        }
        const auto res = world.intersect(p);
        if (nearest) {
            ++hits;
            if (!res.valid())
                ++falseMiss;
            else {
                if (res.item != nearest && std::abs(res.intersection - tmin) > 1E-6)
                    ++wrongItem;
                if (std::abs(res.intersection - tmin) > 1E-6)
                    ++wrongDistance;
            }
        } else if (res.valid()) {
            ++wrongItem;
        }
    }
    check(hits > 2000, "World kd-tree: random rays hit items (" + std::to_string(hits) + ")");
    check(falseMiss == 0, "World kd-tree: no missed intersections (" + std::to_string(falseMiss) + ")");
    check(wrongItem == 0, "World kd-tree: nearest item found (" + std::to_string(wrongItem) + " wrong)");
    check(wrongDistance == 0, "World kd-tree: nearest distance found (" + std::to_string(wrongDistance) + " wrong)");
}

void testTransportToWorld()
{
    W world;
    world.addItem<Sphere>({ 1, { 0, 0, 0 } });
    world.build(4);
    Particle p { .pos = { 0, 0, -100 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    check(world.transportParticleToWorld(p) && std::abs(p.pos[2] + 5) < 1E-5, "World: particle is moved to world boundary");
    Particle q { .pos = { 0, 0, -100 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    check(!world.transportParticleToWorld(q), "World: particle moving away from world is not transported");
    Particle inside { .pos = { 0, 0, 3 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    check(world.transportParticleToWorld(inside) && inside.pos[2] == 3, "World: particle inside world is not moved");

    // World fill material interacts: a dense water world attenuates the beam before an item
    W dense;
    auto& sphere = dense.addItem<Sphere>({ 1, { 0, 0, 0 } });
    dense.build(30);
    PencilBeam<> beam({ 0, 0, -30 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(10);
    beam.setNumberOfParticlesPerExposure(10000);
    Transport transport;
    transport(dense, beam, nullptr, false);
    const auto eAir = sphere.energyScored().energyImparted();
    dense.setMaterial(Material<16>::byNistName("Water, Liquid").value(), 1.0);
    transport(dense, beam, nullptr, false);
    const auto eWater = sphere.energyScored().energyImparted();
    check(eWater < 0.5 * eAir, "World: dense fill material attenuates beam before item");
}

int main()
{
    std::cout << "Testing world\n";
    testItems();
    testKDTree();
    testTransportToWorld();
    return xraymctest::summary("testworld");
}

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
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using xraymctest::checkStatisticallyEqual;

// Two small spheres on a line, the second one only receives scattered photons from the first
template <bool FORCED>
std::array<DoseScore, 2> twoSpheres()
{
    using Sphere = WorldSphere<5, 1, FORCED>;
    World<Sphere> world(2);
    constexpr double r = 0.3;
    auto& s1 = world.template addItem<Sphere>({ r });
    auto& s2 = world.template addItem<Sphere>({ r, { 0, 3, 0 } });
    const auto water = Material<5>::byNistName("Water, Liquid").value();
    s1.setMaterial(water, 1.2);
    s2.setMaterial(water, 1.2);
    world.build();
    PencilBeam<> beam({ -100, 0, 0 }, { 1, 0, 0 }, 40);
    beam.setNumberOfExposures(40);
    beam.setNumberOfParticlesPerExposure(50000);
    Transport transport;
    transport(world, beam);
    return { s1.doseScored(), s2.doseScored() };
}

void testForced()
{
    const auto random = twoSpheres<false>();
    const auto forced = twoSpheres<true>();
    checkStatisticallyEqual(forced[0].dose(), forced[0].standardDeviation(), random[0].dose(), random[0].standardDeviation(), 4, "WorldSphere: forced dose equals random dose in primary sphere");
    checkStatisticallyEqual(forced[1].dose(), forced[1].standardDeviation(), random[1].dose(), random[1].standardDeviation(), 4, "WorldSphere: forced dose equals random dose in scatter sphere");
    check(forced[0].relativeUncertainty() < random[0].relativeUncertainty(), "WorldSphere: forced interactions reduce uncertainty in primary sphere");
    check(random[1].dose() > 0 && random[1].dose() < random[0].dose() / 10, "WorldSphere: scatter sphere receives a small dose");
}

void testDose()
{
    // dose is energy per mass: compare with energy scored for known volume and density
    using Sphere = WorldSphere<5, 1>;
    World<Sphere> world;
    auto& s = world.addItem<Sphere>({ 2 });
    s.setMaterial(Material<5>::byNistName("Water, Liquid").value(), 1.5);
    world.build();
    PencilBeam<> beam({ -10, 0, 0 }, { 1, 0, 0 }, 60);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(10000);
    Transport transport;
    transport(world, beam, nullptr, false);
    const auto energy = s.energyScored().energyImparted();
    s.addEnergyScoredToDoseScore(2.0);
    const auto mass = 4.0 / 3.0 * std::numbers::pi * 8 * 1.5;
    checkRelative(s.doseScored().dose(), energy / mass * 2.0, 1E-12, "WorldSphere: dose is energy per mass times calibration");

    s.setRadius(5);
    checkNear(s.radius(), 5, 0, "WorldSphere: setRadius");
    checkNear(s.AABB()[3], 5, 1E-12, "WorldSphere: AABB follows radius");
    s.setCenter({ 1, 2, 3 });
    checkNear(s.AABB()[0], -4, 1E-12, "WorldSphere: AABB follows center");
    check(s.setNistMaterial("Water, Liquid") && !s.setNistMaterial("Unobtainium"), "WorldSphere: setNistMaterial");
}

int main()
{
    std::cout << "Testing world sphere\n";
    testForced();
    testDose();
    return xraymctest::summary("testworldsphere");
}

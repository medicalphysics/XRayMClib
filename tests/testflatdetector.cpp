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

Copyright 2026 Erlend Andersen
*/

#include "testutils.hpp"

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkStatisticallyEqual;

FlatDetector makeDetector()
{
    FlatDetector d;
    d.setCenter({ 0, 0, 20 });
    d.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    d.setPixelSpacing({ 0.1, 0.2 });
    d.setDetectorDimensions({ 40, 30 });
    return d;
}

// Pixel index for a point in the detector plane, x index major as in FlatDetector
std::size_t pixelIndex(const FlatDetector& d, double x, double y)
{
    const auto ix = static_cast<std::size_t>((x + 0.1 * 40 * 0.5) / 0.1);
    const auto iy = static_cast<std::size_t>((y + 0.2 * 30 * 0.5) / 0.2);
    return ix * d.detectorDimensions()[1] + iy;
}

void testGeometry()
{
    auto d = makeDetector();
    check(d.detectorDimensions() == std::array<std::size_t, 2> { 40, 30 }, "FlatDetector: dimensions");
    check(xraymctest::checkVectorNear(d.normal(), { 0, 0, 1 }, 1E-12, "FlatDetector: normal"), "FlatDetector: normal from direction cosines");
    const auto aabb = d.AABB();
    check(aabb[0] <= -2 && aabb[3] >= 2 && aabb[1] <= -3 && aabb[4] >= 3, "FlatDetector: AABB covers detector area");

    Particle p { .pos = { 0.55, -1.15, 0 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto r = d.intersect(p);
    check(r.valid() && std::abs(r.intersection - 20) < 1E-12, "FlatDetector: ray hits detector plane");
    p.pos = { 2.05, 0, 0 };
    check(!d.intersect(p).valid(), "FlatDetector: ray just outside detector misses");
    p.pos = { 0, 0, 30 };
    check(!d.intersect(p).valid(), "FlatDetector: ray behind detector misses");
    p.pos = { 0, 0, 0 };
    p.dir = { 1, 0, 0 };
    check(!d.intersect(p).valid(), "FlatDetector: ray parallel to detector misses");
}

void testPixelScoring()
{
    // In vacuum every photon is absorbed in the pixel it hits
    World<FlatDetector> world;
    auto& d = world.addItem(makeDetector());
    world.setMaterialDensity(1E-12);
    world.build();
    const double x = 0.55, y = -1.15;
    PencilBeam<> beam({ x, y, 0 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(2);
    beam.setNumberOfParticlesPerExposure(1000);
    Transport transport;
    transport(world, beam, nullptr, false);
    const auto idx = pixelIndex(d, x, y);
    checkNear(d.energyScored(idx).energyImparted(), 60.0 * 2000, 1E-6, "FlatDetector: all energy scored in hit pixel");
    double others = 0;
    for (std::size_t i = 0; i < 40 * 30; ++i)
        if (i != idx)
            others += d.energyScored(i).energyImparted();
    check(others == 0, "FlatDetector: no energy in other pixels");

    // the energy score is the KERMA for unit area per pixel when dose is computed
    d.addEnergyScoredToDoseScore(1.0);
    checkNear(d.doseScored(idx).dose(), 60.0 * 2000 / (0.1 * 0.2), 1E-6, "FlatDetector: dose is energy per pixel area");

    // Shadow of a sphere: the primary photons through the center of a water sphere are attenuated
    World<WorldSphere<5, 1>, FlatDetector> w2;
    w2.reserveNumberOfItems(2);
    auto& sphere = w2.addItem<WorldSphere<5, 1>>({ 2, { 0, 0, 10 } });
    sphere.setMaterial(Material<5>::byNistName("Water, Liquid").value(), 1.0);
    auto& d2 = w2.addItem(makeDetector());
    w2.setMaterialDensity(1E-12);
    w2.build();
    PencilBeam<> b2({ 0.05, 0.1, -10 }, { 0, 0, 1 }, 60);
    b2.setNumberOfExposures(20);
    b2.setNumberOfParticlesPerExposure(10000);
    transport(w2, b2, nullptr, false);
    const auto mu = Material<5>::byNistName("Water, Liquid").value().attenuationValues(60).sum();
    const auto expected = 60.0 * b2.numberOfParticles() * std::exp(-mu * 4);
    const auto& e = d2.energyScored(pixelIndex(d2, 0.05, 0.1));
    // primary photons plus a small amount of scatter into the pixel
    check(e.energyImparted() > expected * 0.99 && e.energyImparted() < expected * 1.05, "FlatDetector: pixel behind sphere receives attenuated primary beam (" + std::to_string(e.energyImparted() / expected) + " of exp(-mu t))");
}

int main()
{
    std::cout << "Testing flat detector\n";
    testGeometry();
    testPixelScoring();
    return xraymctest::summary("testflatdetector");
}

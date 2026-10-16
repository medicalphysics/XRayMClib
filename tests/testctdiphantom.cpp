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

#include "testutils.hpp"

#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkStatisticallyEqual;

// Mean and standard error of the energy in each hole, estimated from independent batches.
// Energy deposition in the small air holes is a mix of many small and few large events,
// the batch estimate of the uncertainty is more robust than the per event estimate.
struct HoleEnergy {
    std::array<double, 5> mean = { 0, 0, 0, 0, 0 };
    std::array<double, 5> sem = { 0, 0, 0, 0, 0 };
};

// Irradiates the phantom with fan beams from four directions
template <bool FORCED>
HoleEnergy irradiate(const std::array<double, 3>& translation, bool pmmaHoles = false)
{
    using Phantom = CTDIPhantom<5, 1, FORCED>;
    World<Phantom> world;
    auto& phantom = world.template addItem<Phantom>({ 16, 15 });
    phantom.translate(translation);
    if (pmmaHoles) {
        const std::string pmma = "Polymethyl Methacralate (Lucite, Perspex)";
        phantom.setHoleMaterial(pmma, NISTMaterials::density(pmma));
    }
    world.build();
    IsotropicMonoEnergyBeam<> beam;
    beam.setEnergy(60);
    const auto a = std::atan(17.0 / 60);
    beam.setCollimationHalfAngles(-a, -0.02, a, 0.02);
    beam.setNumberOfExposures(64);
    beam.setNumberOfParticlesPerExposure(10000);
    Transport transport;
    constexpr std::size_t B = 12;
    std::array<double, 5> sum = { 0, 0, 0, 0, 0 }, sum2 = { 0, 0, 0, 0, 0 };
    for (std::size_t b = 0; b < B; ++b) {
        std::array<double, 5> batch = { 0, 0, 0, 0, 0 };
        for (int k = 0; k < 4; ++k) {
            const auto ang = k * std::numbers::pi / 2;
            const auto src = vectormath::rotate(std::array<double, 3> { -60, 0, 0 }, { 0, 0, 1 }, ang);
            beam.setPosition(vectormath::add(src, translation));
            beam.setDirectionCosines(vectormath::rotate(std::array<double, 3> { 0, 1, 0 }, { 0, 0, 1 }, ang), { 0, 0, 1 });
            transport(world, beam, nullptr, false);
            for (std::size_t i = 0; i < 5; ++i)
                batch[i] += phantom.energyScored(i).energyImparted();
        }
        for (std::size_t i = 0; i < 5; ++i) {
            sum[i] += batch[i];
            sum2[i] += batch[i] * batch[i];
        }
    }
    HoleEnergy res;
    for (std::size_t i = 0; i < 5; ++i) {
        res.mean[i] = sum[i] / B;
        res.sem[i] = std::sqrt(std::max(sum2[i] / B - res.mean[i] * res.mean[i], 0.0) / (B - 1));
    }
    return res;
}

void testGeometry()
{
    CTDIPhantom<5, 1> phantom(16, 15);
    const auto aabb = phantom.AABB();
    checkNear(aabb[3] - aabb[0], 32, 1E-9, "CTDIPhantom: diameter");
    checkNear(aabb[5] - aabb[2], 15, 1E-9, "CTDIPhantom: length");
    Particle p { .pos = { -50, 0, 0 }, .dir = { 1, 0, 0 }, .energy = 60, .weight = 1 };
    const auto r = phantom.intersect(p);
    check(r.valid() && std::abs(r.intersection - 34) < 1E-9, "CTDIPhantom: distance to phantom surface");
    p.pos = { -50, 0, 10 };
    check(!phantom.intersect(p).valid(), "CTDIPhantom: ray beside phantom misses");
}

void testHoles()
{
    const auto forced = irradiate<true>({ 0, 0, 0 });
    // With beams from four directions all peripheral holes receive the same dose
    for (std::size_t i = 2; i < 5; ++i)
        checkStatisticallyEqual(forced.mean[i], forced.sem[i], forced.mean[1], forced.sem[1], 5, "CTDIPhantom: peripheral hole " + std::to_string(i) + " equals peripheral hole 1");
    check(forced.mean[0] < forced.mean[1], "CTDIPhantom: center hole receives less than peripheral holes");

    // Interactions in the small air holes are too rare for analog sampling to be compared reliably,
    // with PMMA filled holes both estimators have many events
    const auto randomPMMA = irradiate<false>({ 0, 0, 0 }, true);
    const auto forcedPMMA = irradiate<true>({ 0, 0, 0 }, true);
    for (std::size_t i = 0; i < 5; ++i)
        checkStatisticallyEqual(randomPMMA.mean[i], randomPMMA.sem[i], forcedPMMA.mean[i], forcedPMMA.sem[i], 5, "CTDIPhantom: forced interactions equal random sampling in PMMA hole " + std::to_string(i));
    check(forcedPMMA.mean[0] > forced.mean[0], "CTDIPhantom: PMMA holes absorb more energy than air holes");

    // A translated phantom gives the same hole doses
    const auto moved = irradiate<true>({ 5, -3, 20 });
    for (std::size_t i = 0; i < 5; ++i)
        checkStatisticallyEqual(moved.mean[i], moved.sem[i], forced.mean[i], forced.sem[i], 5, "CTDIPhantom: translated phantom scores the same in hole " + std::to_string(i));
}

int main()
{
    std::cout << "Testing CTDI phantom\n";
    testGeometry();
    testHoles();
    return xraymctest::summary("testctdiphantom");
}

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

// Monte Carlo transport tests against analytical results and between equivalent geometries.

#include "testgeometry.hpp"
#include "testutils.hpp"

#include "xraymc/beams/flatmonoenergyfield.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using xraymctest::checkStatisticallyEqual;

constexpr int NS = 12; // number of material shells
constexpr int LE = 2; // full low energy correction

// Narrow beam geometry: the fraction of photons reaching a small, distant detector without
// losing energy equals the uncollided transmission exp(-mu * rho * t).
void testNarrowBeamTransmission()
{
    using Box = WorldBox<NS, LE>;
    const auto water = Material<NS>::byNistName("Water, Liquid").value();
    for (double energy : { 30.0, 60.0, 100.0 }) {
        World<Box, FluenceScore> world;
        world.reserveNumberOfItems(2);
        auto& slab = world.addItem<Box>({ { -5, -5, 0, 5, 5, 4 }, water, 1.0 });
        auto& detector = world.addItem<FluenceScore>({ 0.25, { 0, 0, 104 }, { 0, 0, -1 } });
        world.setMaterialDensity(1E-9); // vacuum
        world.build();

        PencilBeam<> beam({ 0, 0, -10 }, { 0, 0, 1 }, energy);
        beam.setNumberOfExposures(50);
        beam.setNumberOfParticlesPerExposure(20000);
        Transport transport;
        transport(world, beam, nullptr, false);

        const auto spec = detector.getSpecter();
        const auto counted = static_cast<double>(spec[static_cast<std::size_t>(energy)].second);
        const double N = static_cast<double>(beam.numberOfParticles());
        const auto T = std::exp(-water.attenuationValues(energy).sum() * 1.0 * 4.0);
        const auto t = counted / N;
        const auto label = "Narrow beam (water 4 cm, " + std::to_string(int(energy)) + " keV)";
        checkStatisticallyEqual(t, std::sqrt(t * (1 - t) / N), T, T * 0.003, 4, label + ": transmission equals exp(-mu rho t)");
        (void)slab;
    }
}

// In a thin slab energy is imparted at the first collision: E_imp = N E (mu_tr / mu) (1 - exp(-mu rho t))
void testThinSlabEnergyImparted()
{
    using Box = WorldBox<NS, LE>;
    const auto water = Material<NS>::byNistName("Water, Liquid").value();
    for (double energy : { 20.0, 60.0, 120.0 }) {
        World<Box> world;
        const double t = 0.02;
        auto& slab = world.addItem<Box>({ { -50, -50, 0, 50, 50, t }, water, 1.0 });
        world.setMaterialDensity(1E-9);
        world.build();
        PencilBeam<> beam({ 0, 0, -1 }, { 0, 0, 1 }, energy);
        beam.setNumberOfExposures(40);
        beam.setNumberOfParticlesPerExposure(100000);
        Transport transport;
        transport(world, beam, nullptr, false);

        const auto att = water.attenuationValues(energy);
        const auto mu = att.sum();
        const auto expected = beam.numberOfParticles() * energy * water.massEnergyTransferAttenuation(energy) / mu * (1 - std::exp(-mu * t));
        const auto& e = slab.energyScored();
        const auto label = "Thin slab (water, " + std::to_string(int(energy)) + " keV)";
        // allow 1 % for scattered photons interacting a second time in the slab
        checkStatisticallyEqual(e.energyImparted(), e.standardDeviation(), expected, 0.01 * expected, 4, label + ": energy imparted equals first collision KERMA");
    }
}

// Calibrated dose in a thin air slab equals the prescribed air KERMA
void testAirKermaCalibration()
{
    using Box = WorldBox<NS, LE>;
    const auto air = Material<NS>::byNistName("Air, Dry (near sea level)").value();
    const double density = 1.0; // a dense air equivalent gives enough interactions in a thin slab
    const double t = 0.01;
    const double energy = 50;
    {
        // pencil beam: the calibration is the KERMA for a fluence of one particle per cm2, a box of 1 cm2 cross section
        World<Box> world;
        auto& box = world.addItem<Box>({ { -0.5, -0.5, 0, 0.5, 0.5, t }, air, density });
        world.setMaterialDensity(1E-9);
        world.build();
        PencilBeam<> beam({ 0, 0, -1 }, { 0, 0, 1 }, energy);
        beam.setAirKerma(2.0);
        beam.setNumberOfExposures(40);
        beam.setNumberOfParticlesPerExposure(100000);
        Transport transport;
        transport(world, beam);
        const auto mu = air.attenuationValues(energy).sum() * density * t;
        const auto expected = 2.0 * (1 - std::exp(-mu)) / mu;
        checkStatisticallyEqual(box.doseScored().dose(), box.doseScored().standardDeviation(), expected, 0.01 * expected, 4, "PencilBeam: calibrated dose in thin air slab equals air KERMA");
    }
    {
        // flat field: box with the same cross section as the field
        World<Box> world;
        auto& box = world.addItem<Box>({ { -2, -1, 0, 2, 1, t }, air, density });
        world.setMaterialDensity(1E-9);
        world.build();
        FlatMonoEnergyField<> beam({ 0, 0, -1 }, energy);
        beam.setLenght(2, 1); // half widths
        beam.setAirKerma(3.0);
        beam.setNumberOfExposures(40);
        beam.setNumberOfParticlesPerExposure(100000);
        Transport transport;
        transport(world, beam);
        const auto mu = air.attenuationValues(energy).sum() * density * t;
        const auto expected = 3.0 * (1 - std::exp(-mu)) / mu;
        checkStatisticallyEqual(box.doseScored().dose(), box.doseScored().standardDeviation(), expected, 0.01 * expected, 4, "FlatMonoEnergyField: calibrated dose in thin air slab equals air KERMA");
    }
}

// A point source in the center of a thick sphere: all energy is absorbed in the sphere
template <std::size_t N>
void testEnergyConservation(const Material<N>& material, double density, double radius, double energy, const std::string& name)
{
    using Sphere = WorldSphere<N, LE>;
    World<Sphere> world;
    auto& sphere = world.addItem<Sphere>({ radius, { 0, 0, 0 } });
    sphere.setMaterial(material, density);
    world.setMaterialDensity(1E-9);
    world.build();
    IsotropicMonoEnergyBeamCircle<> beam({ 0, 0, 0 }, { 0, 0, 1 });
    beam.setCollimationHalfAngle(std::numbers::pi);
    beam.setEnergy(energy);
    beam.setNumberOfExposures(20);
    beam.setNumberOfParticlesPerExposure(20000);
    Transport transport;
    transport(world, beam, nullptr, false);
    const auto total = beam.numberOfParticles() * energy;
    const auto& e = sphere.energyScored();
    check(e.energyImparted() > 0.5 * total, name + ": source inside sphere deposits energy in sphere");
    checkStatisticallyEqual(e.energyImparted(), e.standardDeviation(), total, 0.002 * total, 4, name + ": all energy from a source inside a thick sphere is absorbed");
}

// Sum of energy over scoring indices of an item
template <typename U>
double itemEnergy(const U& item, std::size_t n)
{
    double e = 0;
    for (std::size_t i = 0; i < n; ++i)
        e += item.energyScored(i).energyImparted();
    return e;
}

// Energy imparted in a water cube [-1, 1]^3 by a diverging beam, estimated in batches so that
// correlations between scoring elements of the same history do not bias the uncertainty.
template <typename U>
std::pair<double, double> cubeEnergy(U item, std::size_t nIndices)
{
    World<U> world;
    auto& it = world.addItem(std::move(item));
    world.setMaterialDensity(1E-9);
    world.build();
    IsotropicMonoEnergyBeam<> beam({ 0, 0, -10 }, { { { 1, 0, 0 }, { 0, 1, 0 } } }, 60);
    const auto a = std::atan(1.2 / 9.0);
    beam.setCollimationHalfAngles(-a, -a, a, a);
    beam.setNumberOfExposures(20);
    beam.setNumberOfParticlesPerExposure(5000);
    constexpr std::size_t B = 16;
    double sum = 0, sum2 = 0;
    Transport transport;
    for (std::size_t b = 0; b < B; ++b) {
        transport(world, beam, nullptr, false);
        const auto e = itemEnergy(it, nIndices);
        sum += e;
        sum2 += e * e;
    }
    const auto mean = sum / B;
    const auto sem = std::sqrt((sum2 / B - mean * mean) / (B - 1));
    return { mean, sem };
}

void testGeometryEquivalence()
{
    const auto water = Material<NS>::byNistName("Water, Liquid").value();
    WorldBox<NS, LE> box({ -1, -1, -1, 1, 1, 1 }, water, 1.0);
    const auto ref = cubeEnergy(box, 1);
    check(ref.first > 0, "Geometry equivalence: reference box scores energy");

    auto compare = [&](const std::pair<double, double>& r, const std::string& name) {
        // the batch standard error has 15 degrees of freedom, use a wide limit
        checkStatisticallyEqual(r.first, r.second, ref.first, ref.second, 5, "Geometry equivalence: " + name + " equals WorldBox");
    };
    {
        WorldBox<NS, LE, true> item({ -1, -1, -1, 1, 1, 1 }, water, 1.0);
        compare(cubeEnergy(item, 1), "WorldBox (forced)");
    }
    {
        TriangulatedMesh<NS, LE> item(xraymctest::boxTriangles(1));
        item.setMaterial(water, 1.0);
        compare(cubeEnergy(item, 1), "TriangulatedMesh");
    }
    {
        TetrahedalMesh<NS, LE, false> item(xraymctest::tetCube(1));
        compare(cubeEnergy(item, 12), "TetrahedalMesh");
    }
    {
        TetrahedalMesh<NS, LE, true> item(xraymctest::tetCube(1));
        compare(cubeEnergy(item, 12), "TetrahedalMesh (forced)");
    }
    {
        std::vector<double> dens(4 * 4 * 4, 1.0);
        std::vector<std::uint8_t> mat(dens.size(), 0);
        AAVoxelGrid<NS, LE, 255> item({ 4, 4, 4 }, { 0.5, 0.5, 0.5 }, dens, mat, { water });
        compare(cubeEnergy(item, 64), "AAVoxelGrid");
    }
    {
        // voxels of material index 1 are transparent for TRANSPARENTVOXELS = 1, use material 0 only
        std::vector<double> dens(4 * 4 * 4, 1.0);
        std::vector<std::uint8_t> mat(dens.size(), 0);
        AAVoxelGrid<NS, LE, 1> item({ 4, 4, 4 }, { 0.5, 0.5, 0.5 }, dens, mat, { water, water });
        compare(cubeEnergy(item, 64), "AAVoxelGrid (with transparent material index)");
    }
    {
        WorldBoxGrid<NS, LE> item({ -1, -1, -1, 1, 1, 1 }, water, 1.0);
        item.setVoxelDimensions({ 3, 4, 5 });
        compare(cubeEnergy(item, 60), "WorldBoxGrid");
    }
}

void testProgress()
{
    using Box = WorldBox<5, 1>;
    World<Box> world;
    world.addItem<Box>({ 2 });
    world.build();
    PencilBeam<> beam({ 0, 0, -10 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(13);
    beam.setNumberOfParticlesPerExposure(1000);
    TransportProgress progress;
    Transport transport;
    transport.setNumberOfThreads(3);
    check(transport.numberOfThreads() == 3, "Transport: number of threads");
    transport(world, beam, &progress);
    const auto [done, total] = progress.progress();
    check(done == beam.numberOfParticles() && total == beam.numberOfParticles(), "TransportProgress: all particles are reported as completed");
    check(!progress.message().empty(), "TransportProgress: progress message");
}

int main()
{
    std::cout << "Testing Monte Carlo transport\n";
    testNarrowBeamTransmission();
    testThinSlabEnergyImparted();
    testAirKermaCalibration();
    testEnergyConservation(Material<NS>::byNistName("Water, Liquid").value(), 1.0, 50, 60, "Energy conservation (water, 60 keV)");
    // Lead K-shell fluorescence is emitted with the total photon energy per vacancy and the
    // weight multiplied by the number of photons per vacancy, see testinteractions.
    testEnergyConservation(Material<NS>::byZ(82).value(), 11.35, 10, 100, "Energy conservation (lead, 100 keV)");
    testGeometryEquivalence();
    testProgress();
    return xraymctest::summary("testtransport");
}

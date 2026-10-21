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

#include "xraymc/beams/flatmonoenergyfield.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aafluencegrid.hpp"
#include "xraymc/world/worlditems/worlditemtype.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

static_assert(WorldItemType<AAFluenceGrid<1>>);
static_assert(WorldItemType<AAFluenceGrid<-3>>);

// Axis bookkeeping for every allowed normal
static_assert(AAFluenceGrid<1>::normalAxisIndex() == 0 && AAFluenceGrid<1>::firstAxisIndex() == 1 && AAFluenceGrid<1>::secondAxisIndex() == 2);
static_assert(AAFluenceGrid<2>::normalAxisIndex() == 1 && AAFluenceGrid<2>::firstAxisIndex() == 0 && AAFluenceGrid<2>::secondAxisIndex() == 2);
static_assert(AAFluenceGrid<3>::normalAxisIndex() == 2 && AAFluenceGrid<3>::firstAxisIndex() == 0 && AAFluenceGrid<3>::secondAxisIndex() == 1);
static_assert(AAFluenceGrid<-2>::normalAxisIndex() == 1 && AAFluenceGrid<-2>::normalAxisDirection() == -1);
static_assert(AAFluenceGrid<1>::planeNormal() == std::array<double, 3> { 1, 0, 0 });
static_assert(AAFluenceGrid<-1>::planeNormal() == std::array<double, 3> { -1, 0, 0 });
static_assert(AAFluenceGrid<2>::planeNormal() == std::array<double, 3> { 0, 1, 0 });
static_assert(AAFluenceGrid<-3>::planeNormal() == std::array<double, 3> { 0, 0, -1 });

template <int N>
std::uint64_t totalCount(const AAFluenceGrid<N>& g)
{
    std::uint64_t sum = 0;
    const auto& res = g.resolution();
    for (std::uint64_t i = 0; i < res[0] * res[1]; ++i) {
        const auto s = g.getSpecter(i);
        sum = std::accumulate(s.cbegin(), s.cend(), sum);
    }
    return sum;
}

void testGeometry()
{
    // normal +y, first axis x, second axis z
    AAFluenceGrid<2> g({ 1, 2, 3 }, { 0.5, 2 }, { 4, 3 });
    const auto& aabb = g.AABB();
    checkNear(aabb[0], 0, 1E-12, "AAFluenceGrid: AABB xmin");
    checkNear(aabb[3], 2, 1E-12, "AAFluenceGrid: AABB xmax");
    checkNear(aabb[2], 0, 1E-12, "AAFluenceGrid: AABB zmin");
    checkNear(aabb[5], 6, 1E-12, "AAFluenceGrid: AABB zmax");
    check(aabb[1] < 2 && aabb[4] > 2 && aabb[4] - aabb[1] < 1E-5, "AAFluenceGrid: AABB is thin along normal");
    checkNear(g.pixelArea(), 1, 1E-12, "AAFluenceGrid: pixel area");

    Particle p { .pos = { 1.5, 7, 1 }, .dir = { 0, -1, 0 }, .energy = 60, .weight = 1 };
    const auto r = g.intersect(p);
    check(r.valid() && std::abs(r.intersection - 5) < 1E-12, "AAFluenceGrid: ray against normal hits");

    p.dir = { 0, 1, 0 };
    p.pos[1] = -3;
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray along normal is not scored");

    p.dir = { 0, -1, 0 };
    p.pos[1] = 0;
    check(!g.intersect(p).valid(), "AAFluenceGrid: grid behind ray misses");

    p.pos = { 2.5, 7, 1 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray outside first axis misses");
    p.pos = { 1.5, 7, 6.5 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray outside second axis misses");

    p.pos = { -1, 2, 1 };
    p.dir = { 1, 0, 0 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray parallel to grid misses");

    p.pos = { 1.5, 7, 1 };
    p.dir = { 0, -1, 0 };
    const auto v = g.intersectVisualization<double>(p);
    check(v.valid() && v.normal == std::array<double, 3> { 0, 1, 0 }, "AAFluenceGrid: visualization normal");

    g.translate({ 10, 0, 0 });
    check(!g.intersect(p).valid() && std::abs(g.center()[0] - 11) < 1E-12, "AAFluenceGrid: translate moves grid");
    p.pos[0] = 11.5;
    check(g.intersect(p).valid(), "AAFluenceGrid: ray hits translated grid");

    g.setResolution(8, 3);
    checkNear(g.AABB()[0], 9, 1E-12, "AAFluenceGrid: setResolution updates AABB");
    g.setSpacing(1);
    checkNear(g.AABB()[5], 4.5, 1E-12, "AAFluenceGrid: setSpacing updates AABB");

    // negative normal flips which direction is scored
    AAFluenceGrid<-3> gn({ 0, 0, 0 }, { 1, 1 }, { 2, 2 });
    Particle q { .pos = { 0.5, 0.5, -1 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    check(gn.intersect(q).valid(), "AAFluenceGrid<-3>: ray along +z hits");
    q.pos[2] = 1;
    q.dir = { 0, 0, -1 };
    check(!gn.intersect(q).valid(), "AAFluenceGrid<-3>: ray along -z is not scored");
}

void testIntersectionEdgeCases()
{
    // normal +z, x and y in [-1, 1]
    AAFluenceGrid<3> g({ 0, 0, 0 }, { 0.5, 0.5 }, { 4, 4 });

    Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray starting on plane misses");

    // oblique ray hitting at (0.5, -0.5, 0)
    p.pos = { -0.5, 0.5, 2 };
    p.dir = { 1, -1, -2 };
    vectormath::normalize(p.dir);
    const auto r = g.intersect(p);
    check(r.valid() && std::abs(r.intersection - std::sqrt(6.0)) < 1E-12, "AAFluenceGrid: oblique ray hits at correct distance");

    // grid edges are excluded
    p.pos = { 1, 0, 1 };
    p.dir = { 0, 0, -1 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray on first axis edge misses");
    p.pos = { 0, -1, 1 };
    check(!g.intersect(p).valid(), "AAFluenceGrid: ray on second axis edge misses");
    p.pos = { 0.999, -0.999, 1 };
    check(g.intersect(p).valid(), "AAFluenceGrid: ray just inside corner hits");

    // nearly parallel rays are rejected
    p.pos = { 0, 0, 1E-8 };
    p.dir = { 1, 0, -1E-7 };
    vectormath::normalize(p.dir);
    check(!g.intersect(p).valid(), "AAFluenceGrid: nearly parallel ray misses");

    // missed rays give an empty visualization result
    p.pos = { 5, 5, 1 };
    p.dir = { 0, 0, -1 };
    check(!g.intersectVisualization<double>(p).valid(), "AAFluenceGrid: visualization miss is invalid");
}

// Shoots a ray from the scored side through a known pixel center and checks hit distance,
// AABB and that transport bins the photon in that pixel, for any normal axis.
template <int N>
void testAxisVariant()
{
    using Grid = AAFluenceGrid<N>;
    constexpr auto F = Grid::firstAxisIndex();
    constexpr auto S = Grid::secondAxisIndex();
    constexpr auto A = Grid::normalAxisIndex();
    const std::string name = "AAFluenceGrid<" + std::to_string(N) + ">: ";

    const std::array<double, 3> center = { 1, 2, 3 };
    // first axis in [c - 1, c + 1] with 0.5 cm pixels, second axis in [c - 2, c + 2] with 2 cm pixels
    Grid g(center, { 0.5, 2 }, { 4, 2 });
    checkNear(g.AABB()[F], center[F] - 1, 1E-12, name + "AABB first axis min");
    checkNear(g.AABB()[F + 3], center[F] + 1, 1E-12, name + "AABB first axis max");
    checkNear(g.AABB()[S], center[S] - 2, 1E-12, name + "AABB second axis min");
    checkNear(g.AABB()[S + 3], center[S] + 2, 1E-12, name + "AABB second axis max");

    // center of pixel (2, 1)
    Particle p { .pos = center, .dir = { 0, 0, 0 }, .energy = 30.2, .weight = 1 };
    p.pos[F] += 0.25;
    p.pos[S] += 1;
    p.pos[A] += 4 * Grid::normalAxisDirection();
    p.dir[A] = -Grid::normalAxisDirection();
    const auto r = g.intersect(p);
    check(r.valid() && std::abs(r.intersection - 4) < 1E-12, name + "ray against normal hits");
    if (!r.valid())
        return;
    p.translate(r.intersection);
    auto state = xraymctest::seededState();
    g.transport(p, state);
    check(g.getSpecter(2, 1)[30] == 1 && totalCount(g) == 1, name + "photon binned in the pixel hit");
}

void testSetters()
{
    AAFluenceGrid<2> g;
    check(g.resolution() == std::array<std::uint64_t, 2> { 100, 100 }, "AAFluenceGrid: default resolution");
    check(g.spacing() == std::array<double, 2> { 1, 1 }, "AAFluenceGrid: default spacing");
    checkNear(g.energySpacing(), 1, 1E-12, "AAFluenceGrid: default energy spacing");
    check(g.getSpectrumEnergy().size() == static_cast<std::size_t>(MAX_ENERGY()) + 1, "AAFluenceGrid: default number of energy bins");

    g.setSpacing(0);
    check(g.spacing() == std::array<double, 2> { 0.001, 0.001 }, "AAFluenceGrid: spacing clamped to minimum");
    g.setSpacing({ 0.5, 3 });
    checkNear(g.pixelArea(), 1.5, 1E-12, "AAFluenceGrid: setSpacing from array");
    g.setResolution(0, 0);
    check(g.resolution() == std::array<std::uint64_t, 2> { 1, 1 }, "AAFluenceGrid: zero resolution replaced by one");
    g.setResolution({ 3, 5 });
    check(g.resolution() == std::array<std::uint64_t, 2> { 3, 5 }, "AAFluenceGrid: setResolution from array");
    g.setEnergySpacing(0.01);
    checkNear(g.energySpacing(), 0.1, 1E-12, "AAFluenceGrid: energy spacing clamped to minimum");
    g.setEnergySpacing(7);
    const auto e = g.getSpectrumEnergy();
    check(e.size() == static_cast<std::size_t>(MAX_ENERGY() / 7) + 1 && e.back() + 3.5 >= MAX_ENERGY(), "AAFluenceGrid: energy bins cover MAX_ENERGY");

    // x in [-0.75, 0.75], z in [-7.5, 7.5]
    g.setCenter({ 1, 2, 3 });
    check(g.center() == std::array<double, 3> { 1, 2, 3 }, "AAFluenceGrid: setCenter");
    checkNear(g.AABB()[0], 0.25, 1E-12, "AAFluenceGrid: setCenter updates AABB first axis");
    checkNear(g.AABB()[5], 10.5, 1E-12, "AAFluenceGrid: setCenter updates AABB second axis");
    check(g.AABB()[1] < 2 && g.AABB()[4] > 2, "AAFluenceGrid: setCenter updates AABB normal axis");

    // changing spacing keeps counts
    auto state = xraymctest::seededState();
    Particle p { .pos = { 1, 2, 3 }, .dir = { 0, -1, 0 }, .energy = 60, .weight = 1 };
    g.transport(p, state);
    g.setSpacing(1);
    check(totalCount(g) == 1, "AAFluenceGrid: setSpacing keeps counts");
}

void testDose()
{
    // normal +z, x and y in [-1, 1] with 1 x 0.5 cm pixels
    AAFluenceGrid<3> g({ 0, 0, 0 }, { 1, 0.5 }, { 2, 4 });
    auto state = xraymctest::seededState();
    Particle p;
    for (int i = 0; i < 3; ++i) {
        p = { .pos = { 0.5, 0.25, 0 }, .dir = { 0, 0, -1 }, .energy = 60.2, .weight = 1 };
        g.transport(p, state);
    }
    p = { .pos = { 0.5, 0.25, 0 }, .dir = { 0, 0, -1 }, .energy = 30.7, .weight = 1 };
    g.transport(p, state);

    for (std::size_t i = 0; i < 8; ++i)
        check(g.doseScored(i).dose() == 0, "AAFluenceGrid: dose is zero before conversion");

    const auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
    // dose is evaluated at bin centers
    const double kerma = (3 * 60.5 * air.massEnergyTransferAttenuation(60.5) + 30.5 * air.massEnergyTransferAttenuation(30.5)) / g.pixelArea();
    const std::size_t hit = 1 + 2 * 2;

    g.addEnergyScoredToDoseScore(2);
    checkRelative(g.doseScored(hit).dose(), 2 * kerma, 1E-9, "AAFluenceGrid: pixel air kerma from spectrum");
    bool othersZero = true;
    for (std::size_t i = 0; i < 8; ++i)
        if (i != hit)
            othersZero = othersZero && g.doseScored(i).dose() == 0;
    check(othersZero, "AAFluenceGrid: no kerma in pixels without counts");

    g.addEnergyScoredToDoseScore(2);
    checkRelative(g.doseScored(hit).dose(), 2 * kerma, 1E-9, "AAFluenceGrid: repeated conversion does not accumulate");

    p = { .pos = { 0.5, 0.25, 1 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    const auto v = g.intersectVisualization<double>(p);
    check(v.valid(), "AAFluenceGrid: visualization ray hits");
    checkRelative(v.value, 2 * kerma, 1E-9, "AAFluenceGrid: visualization value is pixel kerma");
    p.pos = { -0.5, 0.25, 1 };
    check(g.intersectVisualization<double>(p).value == 0, "AAFluenceGrid: visualization value of empty pixel");

    g.clearEnergyScored();
    check(g.energyScored().numberOfEvents() == 0 && totalCount(g) == 4, "AAFluenceGrid: clearEnergyScored keeps spectrum counts");

    g.clearDoseScored();
    check(g.doseScored(hit).dose() == 0 && totalCount(g) == 0, "AAFluenceGrid: clearDoseScored resets dose and counts");

    g.transport(p, state);
    g.addEnergyScoredToDoseScore();
    g.setResolution(2, 4);
    check(g.doseScored(0).dose() == 0, "AAFluenceGrid: setResolution resets dose");
}

void testCalibratedAirKerma()
{
    // A calibrated flat field covering the grid gives the prescribed air kerma in every pixel
    World<AAFluenceGrid<-3>> world;
    auto& g = world.addItem<AAFluenceGrid<-3>>({ { 0, 0, 10 }, { 1, 1 }, { 4, 4 } });
    world.setMaterialDensity(1E-12);
    world.build();
    FlatMonoEnergyField<> beam({ 0, 0, 0 }, 60.5);
    beam.setLenght(3);
    beam.setAirKerma(5);
    beam.setNumberOfExposures(16);
    beam.setNumberOfParticlesPerExposure(50000);
    Transport transport;
    transport(world, beam);

    double mean = 0;
    double maxDev = 0;
    for (std::size_t i = 0; i < 16; ++i) {
        const auto d = g.doseScored(i).dose();
        mean += d / 16;
        maxDev = std::max(maxDev, std::abs(d / 5 - 1));
    }
    checkRelative(mean, 5, 0.01, "AAFluenceGrid: mean pixel kerma equals beam air kerma");
    check(maxDev < 0.05, "AAFluenceGrid: every pixel kerma equals beam air kerma (max relative deviation " + std::to_string(maxDev) + ")");
}

void testBinning()
{
    // normal -x, first axis y, second axis z
    AAFluenceGrid<-1> g({ 0, 0, 0 }, { 1, 2 }, { 4, 2 });
    g.setEnergySpacing(2);
    const auto energies = g.getSpectrumEnergy();
    checkNear(energies[0], 1, 1E-12, "AAFluenceGrid: first energy bin center");
    checkNear(energies[30], 61, 1E-12, "AAFluenceGrid: energy bin center");
    checkNear(g.energySpacing(), 2, 1E-12, "AAFluenceGrid: energy spacing");

    auto state = xraymctest::seededState();
    // y in [-2, 2] with 1 cm pixels, z in [-2, 2] with 2 cm pixels
    Particle p { .pos = { 0, 1.5, -1 }, .dir = { 1, 0, 0 }, .energy = 61, .weight = 1 };
    g.transport(p, state);
    check(p.pos[0] > 0, "AAFluenceGrid: transport moves particle past grid");
    p = { .pos = { 0, -1.5, 1.5 }, .dir = { 1, 0, 0 }, .energy = 10.5, .weight = 1 };
    g.transport(p, state);
    g.transport(p, state);

    const auto s30 = g.getSpecter(3, 0);
    check(s30[30] == 1, "AAFluenceGrid: photon counted in its pixel and energy bin");
    const auto s01 = g.getSpecter(0, 1);
    check(s01[5] == 2, "AAFluenceGrid: second pixel counted");
    check(g.getSpecter(0 + 4 * 1) == s01, "AAFluenceGrid: flat pixel index matches 2D index");
    check(totalCount(g) == 3, "AAFluenceGrid: no counts in other bins");
    check(g.energyScored().numberOfEvents() == 3, "AAFluenceGrid: energy score event count");
    checkNear(g.energyScored().energyImparted(), 61 + 2 * 10.5, 1E-9, "AAFluenceGrid: energy score sum");

    // particles on the outer edge are clamped into the edge pixels
    p = { .pos = { 0, 2, 2 }, .dir = { 1, 0, 0 }, .energy = MAX_ENERGY(), .weight = 1 };
    g.transport(p, state);
    const auto sEdge = g.getSpecter(3, 1);
    check(sEdge.back() == 1, "AAFluenceGrid: edge photon clamped into last pixel and energy bin");

    g.clearDoseScored();
    check(totalCount(g) == 0 && g.energyScored().numberOfEvents() == 0, "AAFluenceGrid: clearDoseScored resets counts");

    p = { .pos = { 0, 0.5, 0.5 }, .dir = { 1, 0, 0 }, .energy = 60, .weight = 1 };
    g.transport(p, state);
    g.setEnergySpacing(1);
    check(totalCount(g) == 0, "AAFluenceGrid: setEnergySpacing resets counts");
    g.transport(p, state);
    g.setResolution(2, 2);
    check(totalCount(g) == 0, "AAFluenceGrid: setResolution resets counts");
}

void testSerialization()
{
    AAFluenceGrid<3> g({ 1, -2, 3 }, { 0.5, 0.25 }, { 6, 4 });
    g.setEnergySpacing(5);
    auto state = xraymctest::seededState();
    Particle p { .pos = { 1.1, -1.9, 3 }, .dir = { 0, 0, -1 }, .energy = 42, .weight = 1 };
    g.transport(p, state);
    g.addEnergyScoredToDoseScore(3);

    const auto buffer = g.serialize();
    auto d = AAFluenceGrid<3>::deserialize(buffer);
    check(d.has_value(), "AAFluenceGrid: deserialize succeeds");
    if (!d)
        return;
    check(d->center() == g.center(), "AAFluenceGrid: deserialize center");
    check(d->spacing() == g.spacing(), "AAFluenceGrid: deserialize spacing");
    check(d->resolution() == g.resolution(), "AAFluenceGrid: deserialize resolution");
    check(d->AABB() == g.AABB(), "AAFluenceGrid: deserialize AABB");
    check(d->getSpectrumEnergy() == g.getSpectrumEnergy(), "AAFluenceGrid: deserialize energy bins");
    check(d->getSpecter(3, 2) == g.getSpecter(3, 2) && totalCount(*d) == 1, "AAFluenceGrid: deserialize counts");
    check(d->energyScored().numberOfEvents() == 1, "AAFluenceGrid: deserialize energy score");
    const auto hit = 3 + 6 * 2;
    check(g.doseScored(hit).dose() > 0 && d->doseScored(hit) == g.doseScored(hit), "AAFluenceGrid: deserialize pixel dose");

    p = { .pos = { 1.1, -1.9, 4 }, .dir = { 0, 0, -1 }, .energy = 42, .weight = 1 };
    const auto v = d->intersectVisualization<double>(p);
    check(v.valid() && v.value == g.doseScored(hit).dose(), "AAFluenceGrid: visualization of deserialized grid");

    const auto id = AAFluenceGrid<3>::magicID();
    std::vector<char> tagged(id.cbegin(), id.cend());
    check(AAFluenceGrid<3>::validMagicID(tagged), "AAFluenceGrid: valid magic ID");
    tagged[0] = 'X';
    check(!AAFluenceGrid<3>::validMagicID(tagged), "AAFluenceGrid: invalid magic ID rejected");
}

void testMonoenergeticFluence()
{
    World<AAFluenceGrid<-3>> world;
    auto& g = world.addItem<AAFluenceGrid<-3>>({ { 0, 0, 10 }, { 1, 1 }, { 4, 4 } });
    world.setMaterialDensity(1E-12);
    world.build();
    PencilBeam<> beam({ 0.5, -1.5, 0 }, { 0, 0, 1 }, 60.3);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(2500);
    Transport transport;
    transport(world, beam, nullptr, false);

    // x in [-2, 2] -> pixel 2, y in [-2, 2] -> pixel 0
    const auto spec = g.getSpecter(2, 0);
    check(spec[60] == 10000, "AAFluenceGrid: every photon counted in its pixel and energy bin");
    check(totalCount(g) == 10000, "AAFluenceGrid: no photons in other pixels or bins");
}

void testUniformFluence()
{
    // A narrow isotropic beam far away gives nearly uniform fluence over a small grid
    World<AAFluenceGrid<-3>> world;
    auto& g = world.addItem<AAFluenceGrid<-3>>({ { 0, 0, 100 }, { 1, 1 }, { 4, 4 } });
    world.setMaterialDensity(1E-12);
    world.build();
    IsotropicBeam<> beam({ 0, 0, 0 });
    beam.setEnergySpecter(std::vector<std::pair<double, double>> { { 59.0, 1.0 }, { 61.0, 1.0 } });
    beam.setCollimationHalfAngles(0.03);
    beam.setNumberOfExposures(16);
    beam.setNumberOfParticlesPerExposure(50000);
    Transport transport;
    transport(world, beam, nullptr, false);

    const auto total = static_cast<double>(totalCount(g));
    check(total > 0, "AAFluenceGrid: isotropic beam reaches grid");
    const auto expected = total / 16;
    double chi2 = 0;
    for (std::uint64_t i = 0; i < 16; ++i) {
        const auto s = g.getSpecter(i);
        const auto obs = static_cast<double>(std::accumulate(s.cbegin(), s.cend(), std::uint64_t { 0 }));
        chi2 += (obs - expected) * (obs - expected) / expected;
    }
    const double dof = 15;
    const auto z = (chi2 - dof) / std::sqrt(2.0 * dof);
    check(z < 5, "AAFluenceGrid: uniform fluence over pixels (chi-square z = " + std::to_string(z) + ")");
}

int main()
{
    std::cout << "Testing axis aligned fluence grid\n";
    testGeometry();
    testIntersectionEdgeCases();
    testAxisVariant<1>();
    testAxisVariant<-1>();
    testAxisVariant<2>();
    testAxisVariant<-2>();
    testAxisVariant<3>();
    testAxisVariant<-3>();
    testSetters();
    testBinning();
    testDose();
    testSerialization();
    testMonoenergeticFluence();
    testUniformFluence();
    testCalibratedAirKerma();
    return xraymctest::summary("testaafluencegrid");
}

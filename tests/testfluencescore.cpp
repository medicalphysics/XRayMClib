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

#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

void testGeometry()
{
    FluenceScore f(2, { 1, 1, 5 }, { 0, 0, -1 });
    checkNear(f.area(), std::numbers::pi * 4, 1E-12, "FluenceScore: disc area");
    Particle p { .pos = { 2, 2, 0 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto r = f.intersect(p);
    check(r.valid() && std::abs(r.intersection - 5) < 1E-12, "FluenceScore: ray through disc hits");
    p.pos = { 3, 3, 0 };
    check(!f.intersect(p).valid(), "FluenceScore: ray outside disc radius misses");
    p.pos = { 1, 1, 0 };
    p.dir = { 1, 0, 0 };
    check(!f.intersect(p).valid(), "FluenceScore: ray parallel to disc misses");
    f.setRadius(1);
    checkNear(f.area(), std::numbers::pi, 1E-12, "FluenceScore: setRadius");
}

void testMonoenergeticFluence()
{
    World<FluenceScore> world;
    auto& f = world.addItem<FluenceScore>({ 1, { 0, 0, 10 }, { 0, 0, -1 } });
    world.setMaterialDensity(1E-12);
    world.build();
    PencilBeam<> beam({ 0.1, 0, 0 }, { 0, 0, 1 }, 60.3);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(2500);
    Transport transport;
    transport(world, beam, nullptr, false);
    const auto spec = f.getSpecter();
    check(spec[60].second == 10000, "FluenceScore: every photon counted in its energy bin");
    checkNear(spec[60].first, 60.5, 1E-12, "FluenceScore: bin energy is the bin center");
    std::uint64_t total = 0;
    for (const auto& [e, n] : spec)
        total += n;
    check(total == 10000, "FluenceScore: no photons in other bins");
    const auto fl = f.getFluenceSpecter(10000);
    checkRelative(fl[60].second, 1.0 / f.area(), 1E-12, "FluenceScore: fluence is count per primary per area");

    f.clearDoseScored();
    std::uint64_t cleared = 0;
    for (const auto& [e, n] : f.getSpecter())
        cleared += n;
    check(cleared == 0, "FluenceScore: clearDoseScored resets spectrum");
}

void testSpectrum()
{
    // In vacuum the recorded spectrum equals the beam spectrum
    World<FluenceScore> world;
    auto& f = world.addItem<FluenceScore>({ 50, { 0, 0, 10 }, { 0, 0, -1 } });
    world.setMaterialDensity(1E-12);
    world.build();
    Tube tube(100);
    tube.setAlFiltration(2.5);
    tube.setEnergyResolution(1.0);
    IsotropicBeam<> beam({ 0, 0, 0 });
    beam.setEnergySpecter(tube.getSpecter());
    beam.setCollimationHalfAngles(0.1);
    beam.setNumberOfExposures(20);
    beam.setNumberOfParticlesPerExposure(50000);
    Transport transport;
    transport(world, beam, nullptr, false);
    const auto spec = f.getSpecter();
    const auto ts = tube.getSpecter();
    const double N = static_cast<double>(beam.numberOfParticles());
    double chi2 = 0;
    std::size_t dof = 0;
    for (const auto& [e, w] : ts) {
        // tube bin [e, e+1) is sampled uniformly, fluence bins are [i, i+1)
        const auto idx = static_cast<std::size_t>(e);
        const auto expected = w * N;
        if (expected > 50 && idx < spec.size()) {
            const auto obs = static_cast<double>(spec[idx].second);
            chi2 += (obs - expected) * (obs - expected) / expected;
            ++dof;
        }
    }
    const auto z = (chi2 - dof) / std::sqrt(2.0 * dof);
    check(z < 5, "FluenceScore: recorded spectrum equals beam spectrum (chi-square z = " + std::to_string(z) + ")");
}

int main()
{
    std::cout << "Testing fluence scoring\n";
    testGeometry();
    testMonoenergeticFluence();
    testSpectrum();
    return xraymctest::summary("testfluencescore");
}

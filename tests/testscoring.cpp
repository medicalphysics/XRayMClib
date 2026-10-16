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

#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"

#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

void testEnergyScore()
{
    xraymc::EnergyScore e;
    check(e.energyImparted() == 0 && e.numberOfEvents() == 0 && e.variance() == 0, "EnergyScore: default is empty");
    e.scoreEnergy(2.0);
    e.scoreEnergy(3.0);
    e.scoreEnergy(0.0);
    e.scoreEnergy(-1.0);
    check(e.numberOfEvents() == 2, "EnergyScore: zero and negative energies are not scored");
    checkNear(e.energyImparted(), 5.0, 0, "EnergyScore: sum of energies");
    checkNear(e.energyImpartedSquared(), 13.0, 0, "EnergyScore: sum of squared energies");

    // For many events the variance estimate of the sum approaches the sum of squares
    xraymc::EnergyScore big;
    auto state = xraymctest::seededState();
    double sum2 = 0;
    for (int i = 0; i < 100000; ++i) {
        const auto v = state.randomUniform(0.5, 1.5);
        big.scoreEnergy(v);
        sum2 += v * v;
    }
    checkRelative(big.variance(), sum2, 1E-4, "EnergyScore: variance of sum approximates sum of squares");
    checkNear(big.standardDeviation(), std::sqrt(big.variance()), 1E-12, "EnergyScore: standard deviation");
    checkNear(big.relativeUncertainty(), 1.96 * big.standardDeviation() / big.energyImparted(), 1E-15, "EnergyScore: relative uncertainty is 95% confidence");

    big.clear();
    check(big.energyImparted() == 0 && big.numberOfEvents() == 0 && big.energyImpartedSquared() == 0, "EnergyScore: clear");
    big.set(4, 10, 3);
    check(big.energyImparted() == 4 && big.energyImpartedSquared() == 10 && big.numberOfEvents() == 3, "EnergyScore: set");

    // Scoring is thread safe
    xraymc::EnergyScore shared;
    constexpr int nThreads = 8;
    constexpr int perThread = 100000;
    {
        std::vector<std::jthread> threads;
        for (int t = 0; t < nThreads; ++t)
            threads.emplace_back([&]() {
                for (int i = 0; i < perThread; ++i)
                    shared.scoreEnergy(1.0);
            });
    }
    check(shared.numberOfEvents() == nThreads * perThread, "EnergyScore: concurrent scoring counts all events");
    checkNear(shared.energyImparted(), nThreads * perThread, 0, "EnergyScore: concurrent scoring sums all energy");
}

void testDoseScore()
{
    xraymc::EnergyScore e;
    for (int i = 0; i < 1000; ++i)
        e.scoreEnergy(2.0);
    xraymc::DoseScore d;
    check(d.dose() == 0 && d.numberOfEvents() == 0 && d.relativeUncertainty() == 0, "DoseScore: default is empty");
    const double volume = 2.0, density = 0.5, calibration = 3.0;
    d.addScoredEnergy(e, volume, density, calibration);
    // dose = energy / mass * calibration
    checkRelative(d.dose(), 2000.0 / (volume * density) * calibration, 1E-12, "DoseScore: dose is energy per mass times calibration");
    checkRelative(d.variance(), e.variance() * std::pow(calibration / (volume * density), 2), 1E-12, "DoseScore: variance scales with squared factor");
    check(d.numberOfEvents() == 1000, "DoseScore: number of events");
    checkNear(d.relativeUncertainty(), 1.96 * d.standardDeviation() / d.dose(), 1E-15, "DoseScore: relative uncertainty");

    d.addScoredEnergy(e, volume, density, calibration);
    checkRelative(d.dose(), 2 * 2000.0 / (volume * density) * calibration, 1E-12, "DoseScore: repeated scoring accumulates dose");
    check(d.numberOfEvents() == 2000, "DoseScore: repeated scoring accumulates events");

    d.clear();
    check(d.dose() == 0 && d.variance() == 0 && d.numberOfEvents() == 0, "DoseScore: clear");
    d.set(1, 2, 3);
    check(d.dose() == 1 && d.variance() == 2 && d.numberOfEvents() == 3, "DoseScore: set");
    xraymc::DoseScore d2;
    d2.set(1, 2, 3);
    check(d == d2, "DoseScore: equality");
}

int main()
{
    std::cout << "Testing energy and dose scoring\n";
    testEnergyScore();
    testDoseScore();
    return xraymctest::summary("testscoring");
}

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

#include "xraymc/interactions.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/material/material.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using xraymctest::checkStatisticallyEqual;

using Material = xraymc::Material<16>;

xraymc::Particle photon(double energy)
{
    return xraymc::Particle { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
}

// Integral of f over [-1, 1] using composite Gauss-Legendre quadrature
double integrateCos(const std::function<double(double)>& f, int intervals = 40)
{
    double sum = 0;
    for (int i = 0; i < intervals; ++i) {
        const double a = -1 + 2.0 * i / intervals;
        const double b = a + 2.0 / intervals;
        sum += xraymc::gaussIntegration(a, b, [&](double c) { return f(c); });
    }
    return sum;
}

// Running mean and standard error of the mean
struct Stat {
    double sum = 0, sum2 = 0;
    std::size_t n = 0;
    void add(double v)
    {
        sum += v;
        sum2 += v * v;
        ++n;
    }
    double mean() const { return sum / n; }
    double sem() const { return std::sqrt((sum2 / n - mean() * mean()) / n); }
};

// Klein-Nishina differential cross section per unit cos(theta), up to a constant
double kleinNishina(double k, double c)
{
    const auto e = 1 / (1 + k * (1 - c));
    return e * e * (e + 1 / e - (1 - c * c));
}

template <int L>
void testComptonKinematics(const Material& mat, double energy, const std::string& name)
{
    auto state = xraymctest::seededState();
    const auto k = energy / xraymc::ELECTRON_REST_MASS();
    constexpr std::size_t N = 1000000;
    Stat cosStat, eStat;
    bool kinematicsOk = true;
    bool unitOk = true;
    bool energyImpartedOk = true;
    for (std::size_t i = 0; i < N; ++i) {
        auto p = photon(energy);
        const auto Ei = xraymc::interactions::comptonScatter<16, L>(p, mat, state);
        const auto c = p.dir[2];
        const auto eExpected = 1 / (1 + k * (1 - c));
        kinematicsOk = kinematicsOk && std::abs(p.energy / energy - eExpected) < 1E-6;
        unitOk = unitOk && std::abs(xraymc::vectormath::length(p.dir) - 1) < 1E-9;
        energyImpartedOk = energyImpartedOk && std::abs(Ei - (energy - p.energy)) < 1E-9;
        cosStat.add(c);
        eStat.add(p.energy / energy);
    }
    check(kinematicsOk, name + ": scattered energy follows Compton formula for sampled angle");
    check(unitOk, name + ": scattered direction is unit vector");
    check(energyImpartedOk, name + ": energy imparted equals energy lost by photon");

    // Analytical expectations
    std::function<double(double)> pdf;
    if constexpr (L == 0) {
        pdf = [=](double c) { return kleinNishina(k, c); };
    } else {
        pdf = [&, k](double c) {
            const auto x = Material::momentumTransferCosAngle(energy, c);
            return kleinNishina(k, c) * mat.scatterFactor(std::max(x, 1E-6));
        };
    }
    const auto norm = integrateCos(pdf);
    const auto meanCos = integrateCos([&](double c) { return c * pdf(c); }) / norm;
    const auto meanE = integrateCos([&](double c) { return pdf(c) / (1 + k * (1 - c)); }) / norm;
    checkStatisticallyEqual(cosStat.mean(), cosStat.sem(), meanCos, 0, 4, name + ": mean cos(theta) agrees with analytical distribution");
    checkStatisticallyEqual(eStat.mean(), eStat.sem(), meanE, 0, 4, name + ": mean energy ratio agrees with analytical distribution");
}

void testComptonIA()
{
    const auto water = Material::byNistName("Water, Liquid").value();
    const auto lead = Material::byZ(82).value();
    auto state = xraymctest::seededState();
    // The material table is derived from the EPDL incoherent scatter function while the
    // impulse approximation includes Doppler broadening, the models differ more for high Z.
    for (const auto& [mat, name, tolerance] : { std::make_tuple(&water, std::string("water"), 0.01), std::make_tuple(&lead, std::string("lead"), 0.03) }) {
        for (double energy : { 20.0, 60.0, 140.0 }) {
            constexpr std::size_t N = 500000;
            Stat eStat;
            bool bounds = true;
            bool balance = true;
            for (std::size_t i = 0; i < N; ++i) {
                auto p = photon(energy);
                const auto Ei = xraymc::interactions::comptonScatter<16, 2>(p, *mat, state);
                bounds = bounds && p.energy >= 0 && p.energy <= energy;
                balance = balance && std::abs(Ei - (energy - p.energy)) < 1E-9;
                eStat.add(p.energy);
            }
            const auto label = "Compton IA (" + name + ", " + std::to_string(int(energy)) + " keV)";
            check(bounds, label + ": scattered energy within [0, E]");
            check(balance, label + ": energy imparted equals energy lost by photon");
            // The material table holds the mean scattered energy from the same physics data
            checkRelative(eStat.mean(), mat->meanIncoherentScatterEnergy(energy), tolerance, label + ": mean scattered energy agrees with material table");
        }
    }
}

template <int L>
void testRayleigh(const Material& mat, double energy, const std::string& name)
{
    auto state = xraymctest::seededState();
    constexpr std::size_t N = 1000000;
    Stat cosStat, cos2Stat;
    bool energyOk = true;
    bool unitOk = true;
    for (std::size_t i = 0; i < N; ++i) {
        auto p = photon(energy);
        xraymc::interactions::rayleightScatter<16, L>(p, mat, state);
        energyOk = energyOk && p.energy == energy;
        unitOk = unitOk && std::abs(xraymc::vectormath::length(p.dir) - 1) < 1E-9;
        cosStat.add(p.dir[2]);
        cos2Stat.add(p.dir[2] * p.dir[2]);
    }
    check(energyOk, name + ": coherent scatter preserves energy");
    check(unitOk, name + ": scattered direction is unit vector");
    std::function<double(double)> pdf;
    if constexpr (L == 0) {
        pdf = [](double c) { return 1 + c * c; };
    } else {
        pdf = [&](double c) {
            const auto x = Material::momentumTransferCosAngle(energy, c);
            const auto ff = mat.formFactor(std::max(x, 1E-6));
            return (1 + c * c) * ff * ff;
        };
    }
    const auto norm = integrateCos(pdf);
    const auto meanCos = integrateCos([&](double c) { return c * pdf(c); }) / norm;
    const auto meanCos2 = integrateCos([&](double c) { return c * c * pdf(c); }) / norm;
    // The form factor is sampled from a tabulated inverse CDF, allow a small systematic error
    const double sys = L == 0 ? 0 : 0.005;
    checkStatisticallyEqual(cosStat.mean(), cosStat.sem(), meanCos, sys, 4, name + ": mean cos(theta) agrees with analytical distribution");
    checkStatisticallyEqual(cos2Stat.mean(), cos2Stat.sem(), meanCos2, sys, 4, name + ": mean cos^2(theta) agrees with analytical distribution");
}

void testPhotoelectric()
{
    const auto water = Material::byNistName("Water, Liquid").value();
    auto state = xraymctest::seededState();

    // Without binding effects the photon is absorbed
    {
        auto p = photon(60);
        p.weight = 0.5;
        const auto pe = water.attenuationValues(60).photoelectric;
        const auto Ei = xraymc::interactions::photoelectricEffect<16, 0>(pe, p, water, state);
        check(p.energy == 0, "Photoelectric (no binding): photon is absorbed");
        checkNear(Ei, 30, 1E-12, "Photoelectric (no binding): energy imparted is weighted photon energy");
    }

    // Water has no fluorescence above the energy cut off, everything is absorbed
    {
        bool absorbed = true;
        const auto pe = water.attenuationValues(60).photoelectric;
        for (int i = 0; i < 100000; ++i) {
            auto p = photon(60);
            const auto Ei = xraymc::interactions::photoelectricEffectIA(pe, p, water, state);
            absorbed = absorbed && p.energy == 0 && std::abs(Ei - 60) < 1E-12;
        }
        check(absorbed, "Photoelectric IA (water): photon fully absorbed");
    }

    // Lead above the K-edge: compare emission probability with shell cross sections and check energy balance
    {
        const auto lead = Material::byZ(82).value();
        const double E = 100;
        const auto pe = lead.attenuationValues(E).photoelectric;
        const auto kProb = lead.attenuationPhotoelectricShell(0, E) / pe;
        const auto& kShell = lead.shell(0);

        constexpr std::size_t N = 1000000;
        std::size_t nK = 0;
        bool belowIncident = true;
        Stat balance;
        for (std::size_t i = 0; i < N; ++i) {
            auto p = photon(E);
            const auto Ei = xraymc::interactions::photoelectricEffectIA(pe, p, lead, state);
            if (p.energy > 0 && std::abs(p.energy - kShell.energyOfPhotonsPerInitVacancy) < 1E-9)
                ++nK;
            belowIncident = belowIncident && p.energy < E && Ei >= 0;
            // energy imparted + energy carried by the (weighted) fluorescence photon
            balance.add(Ei + p.energy * p.weight);
        }
        check(belowIncident, "Photoelectric IA (lead): fluorescence energy below incident energy");
        const double fK = static_cast<double>(nK) / N;
        checkStatisticallyEqual(fK, std::sqrt(fK * (1 - fK) / N), kProb, 0, 4, "Photoelectric IA (lead): K-shell selection probability equals shell cross section fraction");
        // Energy conservation: on average the energy leaving the interaction cannot exceed the photon energy.
        // EADL gives the total energy of fluorescence photons per initial vacancy (77.4 keV for the lead
        // K-shell, with 1.32 photons per vacancy). Emitting a photon with the total energy and weight
        // multiplied by the number of photons per vacancy overestimates the fluorescence energy.
        check(balance.mean() <= E * (1 + 1E-6), "Photoelectric IA (lead): energy imparted + fluorescence energy <= incident energy (mean " + std::to_string(balance.mean()) + " keV for " + std::to_string(E) + " keV photons)");
    }
}

void testInteract()
{
    const auto water = Material::byNistName("Water, Liquid").value();
    auto state = xraymctest::seededState();
    for (double energy : { 15.0, 40.0, 100.0 }) {
        const auto att = water.attenuationValues(energy);
        constexpr std::size_t N = 1000000;
        std::size_t npe = 0, ninc = 0, ncoh = 0;
        bool consistent = true;
        for (std::size_t i = 0; i < N; ++i) {
            auto p = photon(energy);
            const auto res = xraymc::interactions::interact<16, 2>(att, p, water, state);
            const int types = res.interactionWasPhotoelectric + res.interactionWasIncoherent + res.interactionWasCoherent;
            consistent = consistent && types == 1 && res.particleDirectionChanged;
            consistent = consistent && res.energyImparted >= 0 && p.energy <= energy;
            consistent = consistent && (res.particleAlive ? p.energy >= xraymc::MIN_ENERGY() : p.energy == 0);
            if (res.interactionWasCoherent)
                consistent = consistent && p.energy == energy && res.energyImparted == 0 && !res.particleEnergyChanged;
            npe += res.interactionWasPhotoelectric;
            ninc += res.interactionWasIncoherent;
            ncoh += res.interactionWasCoherent;
        }
        const auto label = "interact (water, " + std::to_string(int(energy)) + " keV)";
        check(consistent, label + ": result flags consistent with particle state");
        for (const auto& [n, mu, t] : { std::make_tuple(npe, att.photoelectric, "photoelectric"), std::make_tuple(ninc, att.incoherent, "incoherent"), std::make_tuple(ncoh, att.coherent, "coherent") }) {
            const double f = static_cast<double>(n) / N;
            const double fexp = mu / att.sum();
            checkStatisticallyEqual(f, std::sqrt(fexp * (1 - fexp) / N), fexp, 0, 4, label + ": " + t + " fraction equals attenuation fraction");
        }
    }
}

void testRussianRoulette()
{
    // Low weight particles are killed with probability 0.9, survivors get 10x weight: the expected weight is unchanged
    const auto water = Material::byNistName("Water, Liquid").value();
    auto state = xraymctest::seededState();
    // Coherent only attenuation makes the interaction energy independent of the outcome
    xraymc::AttenuationValues att { .photoelectric = 0, .incoherent = 0, .coherent = 1 };
    constexpr std::size_t N = 1000000;
    const double w0 = 0.05;
    Stat weight;
    std::size_t survived = 0;
    for (std::size_t i = 0; i < N; ++i) {
        auto p = photon(60);
        p.weight = w0;
        const auto res = xraymc::interactions::interact<16, 2>(att, p, water, state);
        weight.add(res.particleAlive ? p.weight : 0);
        survived += res.particleAlive;
    }
    const double fs = static_cast<double>(survived) / N;
    const double fsExp = 1 - xraymc::interactions::russianRuletteProbability();
    checkStatisticallyEqual(fs, std::sqrt(fsExp * (1 - fsExp) / N), fsExp, 0, 4, "Russian roulette: survival probability");
    checkStatisticallyEqual(weight.mean(), weight.sem(), w0, 0, 4, "Russian roulette: expected weight is preserved");

    // Particles above the threshold are not affected
    std::size_t killed = 0;
    for (std::size_t i = 0; i < 10000; ++i) {
        auto p = photon(60);
        p.weight = 0.5;
        const auto res = xraymc::interactions::interact<16, 2>(att, p, water, state);
        killed += !res.particleAlive;
    }
    check(killed == 0, "Russian roulette: not applied above weight threshold");
}

void testForced()
{
    // The forced interaction estimator must have the same expectation value as analog sampling
    // of the first interaction in a slab of thickness L
    const auto water = Material::byNistName("Water, Liquid").value();
    auto state = xraymctest::seededState();
    for (double energy : { 20.0, 60.0 }) {
        const auto att = water.attenuationValues(energy);
        const double density = 1;
        const double L = 0.5;
        constexpr std::size_t N = 2000000;
        Stat analog, forced;
        bool forcedOk = true;
        for (std::size_t i = 0; i < N; ++i) {
            {
                auto p = photon(energy);
                const auto s = -std::log(state.randomUniform()) / (att.sum() * density);
                double Ei = 0;
                if (s < L) {
                    p.translate(s);
                    Ei = xraymc::interactions::interact<16, 2>(att, p, water, state).energyImparted;
                }
                analog.add(Ei);
            }
            {
                auto p = photon(energy);
                const auto res = xraymc::interactions::interactForced<16, 2>(L, density, att, p, water, state);
                forced.add(res.energyImparted);
                // particle is either still within the slab after a scatter event or moved beyond it
                forcedOk = forcedOk && p.pos[2] >= 0 && (res.particleDirectionChanged ? p.pos[2] <= L : p.pos[2] > L);
            }
        }
        const auto label = "interactForced (water, " + std::to_string(int(energy)) + " keV)";
        check(forcedOk, label + ": particle position consistent with outcome");
        checkStatisticallyEqual(forced.mean(), forced.sem(), analog.mean(), analog.sem(), 4, label + ": mean energy imparted equals analog sampling");
        check(forced.sem() < analog.sem(), label + ": forced estimator has lower variance");
    }
}

void testTracking()
{
    const auto water = Material::byNistName("Water, Liquid").value();
    auto state = xraymctest::seededState();
    xraymc::ParticleTrack p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    for (int i = 0; i < 3; ++i) {
        p.translate(1);
        xraymc::interactions::rayleightScatter<16, 2>(p, water, state);
    }
    check(p.getSize() == 3, "ParticleTrack: scatter events register positions");
    const auto h = p.getHistory();
    check(std::abs(h[0][2] - 1) < 1E-12, "ParticleTrack: first registered position is first interaction point");
    for (int i = 0; i < 10; ++i)
        p.registerPosition();
    check(p.getSize() == xraymc::ParticleTrack::N, "ParticleTrack: history is capped at N positions");
}

int main()
{
    std::cout << "Testing photon interactions\n";
    const auto water = Material::byNistName("Water, Liquid").value();
    const auto lead = Material::byZ(82).value();
    testComptonKinematics<0>(water, 60, "Compton Klein-Nishina (60 keV)");
    testComptonKinematics<0>(water, 140, "Compton Klein-Nishina (140 keV)");
    testComptonKinematics<1>(water, 30, "Compton Livermore (water, 30 keV)");
    testComptonKinematics<1>(lead, 100, "Compton Livermore (lead, 100 keV)");
    testComptonIA();
    testRayleigh<0>(water, 60, "Rayleigh dipole (60 keV)");
    testRayleigh<2>(water, 20, "Rayleigh form factor (water, 20 keV)");
    testRayleigh<2>(lead, 80, "Rayleigh form factor (lead, 80 keV)");
    testPhotoelectric();
    testInteract();
    testRussianRoulette();
    testForced();
    testTracking();
    return xraymctest::summary("testinteractions");
}

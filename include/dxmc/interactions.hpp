/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

namespace dxmc {
namespace interactions {

    constexpr static double russianRuletteProbability()
    {
        return 0.9;
    }

    constexpr static double russianRuletteWeightThreshold()
    {
        return 0.1;
    }

    template <std::size_t Nshells, int LOWENERGYCORRECTION = 2, ParticleType P = Particle>
    void rayleightScatter(P& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            particle.registerPosition();
        }

        if constexpr (LOWENERGYCORRECTION == 0) {
            bool reject;
            double theta;
            do {
                constexpr auto extreme = (4 * std::numbers::sqrt2_v<double>) / (3 * std::numbers::sqrt3_v<double>);
                const auto r1 = state.randomUniform<double>(0, extreme);
                theta = state.randomUniform<double>(0, PI_VAL<double>());
                const auto sinang = std::sin(theta);
                reject = r1 > ((2 - sinang * sinang) * sinang);
            } while (reject);
            // calc angle and add randomly 90 degrees since dist i symetrical
            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            particle.dir = vectormath::peturb(particle.dir, theta, phi);

        } else {
            // theta is scattering angle
            // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

            // finding qmax
            const auto qmax = material.momentumTransferMax(particle.energy);
            const auto qmax_squared = qmax * qmax;

            double cosAngle;
            do {
                const auto q_squared = material.sampleSquaredMomentumTransferFromFormFactorSquared(qmax_squared, state);
                cosAngle = 1 - 2 * q_squared / qmax_squared;
            } while ((1 + cosAngle * cosAngle) * 0.5 < state.randomUniform());

            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            const auto theta = std::acos(cosAngle);
            particle.dir = vectormath::peturb(particle.dir, theta, phi);
        }
    }

    template <std::size_t Nshells>
    int comptonScatterIA_NRC_sample_shell(const ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // Binding energy og last shell should be set to zero
        const auto& shells = material.shells();
        int shell_idx;
        do {
            shell_idx = 0;
            const auto r = state.randomUniform();
            double acc = shells[0].numberOfElectronsFraction;
            while (acc < r) {
                ++shell_idx;
                acc += shells[shell_idx].numberOfElectronsFraction;
            }
        } while (shells[shell_idx].bindingEnergy > particle.energy && shell_idx != material.numberOfShells());
        return shell_idx;
    }

    template <std::size_t Nshells>
    double comptonScatterIA(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        const auto k = particle.energy / ELECTRON_REST_MASS();
        const auto E = particle.energy;
        double S;
        double cosTheta;
        double e;
        double alpha;
        double p;
        double expb;
        double pi;
        double kc;
        double U;
        int shell_idx;
        do {
            // sample cosTheta
            shell_idx = comptonScatterIA_NRC_sample_shell(particle, material, state);
            const auto e_min = 1 / (1 + 2 * k);
            const auto g_max = 1 / e_min + e_min;
            double g;

            do {
                const double r1 = state.randomUniform();
                e = r1 + (1 - r1) * e_min;
                cosTheta = 1 - std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with arg > 2 (better way?)
                const double sinThetaSqr = 1 - cosTheta * cosTheta;
                g = 1 / e + e - sinThetaSqr;
            } while (state.randomUniform(g_max) > g);

            // Binding energy in units of mec2
            U = shell_idx < material.numberOfShells() - 1 ? material.shells()[shell_idx].bindingEnergy / ELECTRON_REST_MASS() : 0.0;

            // calculate pz max value; pi
            pi = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);

            // Calculate S
            const auto J0 = material.shells()[shell_idx].HartreeFockOrbital_0;
            kc = k * e;
            const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
            alpha = qc / k + kc * (kc - k * cosTheta) / (k * qc);
            const auto b_part = (1 + 2 * J0 * std::abs(pi));
            expb = std::exp(-(0.5 * b_part * b_part - 0.5));

            p = std::sqrt((U + 1) * (U + 1) - 1);

            constexpr bool USE_SECOND_ORDER_APPROXIMATION = false;
            if (pi < -p) {
                S = (1 - alpha * p) * expb * 0.5;
            } else if (pi <= 0) {
                const auto S1 = (1 - alpha * pi) * expb * 0.5;
                if constexpr (USE_SECOND_ORDER_APPROXIMATION) {
                    const double pik = std::exp(0.5) / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                    const auto S2 = -alpha / (4 * J0) * pik;
                    const auto errf1 = std::erf((1 + 2 * J0 * p) / std::numbers::sqrt2);
                    const auto errf2 = std::erf((1 + 2 * J0 * std::abs(pi)) / std::numbers::sqrt2);
                    S = S1 + S2 * (errf1 - errf2);
                } else {
                    S = S1;
                }
            } else if (pi < p) {
                const auto S1 = (1 - alpha * pi) * expb * 0.5;
                if constexpr (USE_SECOND_ORDER_APPROXIMATION) {
                    const double pik = std::exp(0.5) / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                    const auto S2 = -alpha / (4 * J0) * pik;
                    const auto errf1 = std::erf((1 + 2 * J0 * p) / std::numbers::sqrt2);
                    const auto errf2 = std::erf((1 + 2 * J0 * std::abs(pi)) / std::numbers::sqrt2);
                    S = 1 - S1 + S2 * (errf1 - errf2);
                } else {
                    S = 1 - S1;
                }
            } else {
                S = 1 - (1 + alpha * p) * expb * 0.5;
            }
        } while (state.randomUniform() > S);

        // sample pz
        double pz;

        constexpr bool F_CORRECTION = false;
        if constexpr (F_CORRECTION) {
            const auto J0 = material.shells()[shell_idx].HartreeFockOrbital_0;
            double Fmax, F;
            if (pi < -p) {
                Fmax = 1 - alpha * p;
            } else if (pi < p) {
                Fmax = 1 + alpha * pi;
            } else {
                Fmax = 1 + alpha * p;
            }
            do {
                const auto r = state.randomUniform(S);
                if (r < 0.5) {
                    const auto logarg = std::max(2 * r, std::numeric_limits<double>::epsilon());
                    const auto sqrtarg = 1 - 2 * std::log(logarg);
                    pz = (1 - std::sqrt(sqrtarg)) / (2 * J0);
                } else {
                    const auto logarg = std::max(2 * (1 - r), std::numeric_limits<double>::epsilon());
                    const auto sqrtarg = 1 - 2 * std::log(logarg);
                    pz = (std::sqrt(sqrtarg) - 1) / (2 * J0);
                }
                if (pz < -p) {
                    F = 1 - alpha * p;
                } else if (pz < p) {
                    F = 1 + alpha * pz;
                } else {
                    F = 1 + alpha * p;
                }
            } while (state.randomUniform(Fmax) > F);
        } else {
            const auto J0 = material.shells()[shell_idx].HartreeFockOrbital_0;
            const auto r = state.randomUniform(S);
            if (r < 0.5) {
                const auto logarg = std::max(2 * r, std::numeric_limits<double>::epsilon());
                const auto sqrtarg = 1 - 2 * std::log(logarg);
                pz = (1 - std::sqrt(sqrtarg)) / (2 * J0);
            } else {
                const auto logarg = std::max(2 * (1 - r), std::numeric_limits<double>::epsilon());
                const auto sqrtarg = 1 - 2 * std::log(logarg);
                pz = (std::sqrt(sqrtarg) - 1) / (2 * J0);
            }
        }

        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);

        // calculating kbar
        const auto bar_part = 1 - 2 * e * cosTheta + e * e * (1 - pz * pz * (1 - cosTheta * cosTheta));
        const auto kbar = kc * (1 - pz * pz * e * cosTheta + pz * std::sqrt(bar_part)) / (1 - pz * pz * e * e);
        particle.energy = kbar * ELECTRON_REST_MASS();
        return (E - particle.energy) * particle.weight;
    }

    template <std::size_t Nshells, int LOWENERGYCORRECTION = 2, ParticleType P = Particle>
    auto comptonScatter(P& particle, const Material<Nshells>& material, RandomState& state) noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            particle.registerPosition();
        }

        if constexpr (LOWENERGYCORRECTION == 2) {
            return comptonScatterIA(particle, material, state);
        } else {
            const auto k = particle.energy / ELECTRON_REST_MASS();
            const auto emin = 1 / (1 + 2 * k);
            const auto gmaxInv = emin / (1 + emin * emin);

            double e, cosTheta;
            bool rejected;
            do {
                const auto r1 = state.randomUniform();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                if constexpr (LOWENERGYCORRECTION == 0) {
                    rejected = state.randomUniform() > g;
                } else { // Livermore
                    const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                    const auto scatterFactor = material.scatterFactor(q);
                    rejected = state.randomUniform(material.effectiveZ()) > (g * scatterFactor);
                }
            } while (rejected);

            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            const auto theta = std::acos(cosTheta);
            particle.dir = vectormath::peturb(particle.dir, theta, phi);

            const auto E = particle.energy;
            particle.energy *= e;
            return (E - particle.energy) * particle.weight;
        }
    }

    template <std::size_t Nshells>
    auto photoelectricEffectIA(const double totalPhotoCrossSection, ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // finding shell based on photoelectric cross section
        const std::uint_fast8_t max_shell = material.numberOfShells();
        std::uint_fast8_t shell = 0;
        auto prob = state.randomUniform();
        bool next;
        do {
            const auto& sh = material.shell(shell);
            if (sh.bindingEnergy < MIN_ENERGY()) {
                shell = max_shell;
            }
            next = shell != max_shell;
            if (next && sh.bindingEnergy < particle.energy) {
                const auto shellCS = material.attenuationPhotoelectricShell(shell, particle.energy);
                const auto shellProb = shellCS / totalPhotoCrossSection;
                prob -= shellProb;
                next = prob > 0;
            }
            if (next)
                ++shell;
        } while (next);

        auto E = particle.energy * particle.weight;
        particle.energy = 0;
        if (shell != max_shell) {
            const auto& s = material.shell(shell);
            if (s.energyOfPhotonsPerInitVacancy > MIN_ENERGY()) {
                particle.energy = s.energyOfPhotonsPerInitVacancy;
                E -= particle.energy * particle.weight;
                particle.weight *= s.numberOfPhotonsPerInitVacancy;
                const auto costheta = 1 - 2 * state.randomUniform();
                const auto theta = std::acos(costheta);
                const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
                particle.dir = vectormath::peturb(particle.dir, theta, phi);
            }
        }
        return E;
    }

    template <int Nshells, int LOWENERGYCORRECTION = 2, ParticleType P = Particle>
    auto photoelectricEffect(const double totalPhotoCrossSection, P& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            particle.registerPosition();
        }

        if constexpr (LOWENERGYCORRECTION == 2) {
            return photoelectricEffectIA(totalPhotoCrossSection, particle, material, state);
        } else {
            const auto E = particle.energy * particle.weight;
            particle.energy = 0;
            return E;
        }
    }

    struct InteractionResult {
        double energyImparted = 0;
        bool particleAlive = true;
        bool particleEnergyChanged = false;
        bool particleDirectionChanged = false;
    };

    template <std::size_t Nshells, int LOWENERGYCORRECTION = 2>
    InteractionResult interact(const AttenuationValues& attenuation, ParticleType auto& particle, const Material<Nshells>& material, RandomState& state)
    {
        InteractionResult res;
        const auto r2 = state.randomUniform(attenuation.sum());
        if (r2 < attenuation.photoelectric) {
            const auto Ei = interactions::photoelectricEffect<Nshells, LOWENERGYCORRECTION>(attenuation.photoelectric, particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
            res.particleDirectionChanged = true;
        } else if (r2 < (attenuation.photoelectric + attenuation.incoherent)) {
            const auto Ei = interactions::comptonScatter<Nshells, LOWENERGYCORRECTION>(particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
            res.particleDirectionChanged = true;
        } else {
            interactions::rayleightScatter<Nshells, LOWENERGYCORRECTION>(particle, material, state);
            res.particleDirectionChanged = true;
        }
        if (particle.energy < MIN_ENERGY()) {
            res.particleAlive = false;
            res.energyImparted += particle.energy * particle.weight;
        } else {
            if (particle.weight < interactions::russianRuletteWeightThreshold() && res.particleAlive) {
                if (state.randomUniform() < interactions::russianRuletteProbability()) {
                    res.particleAlive = false;
                    particle.energy = 0;
                } else {
                    constexpr auto factor = 1 / (1 - interactions::russianRuletteProbability());
                    particle.weight *= factor;
                    res.particleAlive = true;
                }
            }
        }
        return res;
    }

    template <std::size_t NMaterialShells, int LOWENERGYCORRECTION = 2>
    InteractionResult interactForced(double maxStepLen, double materialDensity, const AttenuationValues& attenuation, ParticleType auto& particle, const Material<NMaterialShells>& material, RandomState& state)
    {
        InteractionResult intRes;
        const auto relativePeProbability = attenuation.photoelectric / attenuation.sum();
        const auto attSum = attenuation.sum() * materialDensity;
        const auto probNotInteraction = std::exp(-attSum * maxStepLen);
        // Forced photoelectric effect
        intRes.energyImparted = particle.energy * particle.weight * (1 - probNotInteraction) * relativePeProbability;

        // Remainder probability
        const auto p1 = state.randomUniform();
        if (p1 > probNotInteraction) {
            if (state.randomUniform() > relativePeProbability) { // scatter interaction happens
                // Translate particle to interaction point
                const auto stepLen = -std::log(p1) / attSum;
                particle.translate(stepLen);

                // Decide what scatter interaction
                const auto r2 = state.randomUniform(attenuation.incoherent + attenuation.coherent);
                if (r2 < attenuation.incoherent) {
                    const auto Ei = interactions::comptonScatter<NMaterialShells, LOWENERGYCORRECTION>(particle, material, state);
                    intRes.energyImparted += Ei;
                    intRes.particleEnergyChanged = true;
                    intRes.particleDirectionChanged = true;
                } else {
                    interactions::rayleightScatter<NMaterialShells, LOWENERGYCORRECTION>(particle, material, state);
                    intRes.particleDirectionChanged = true;
                }

                // Handle cutoff values
                if (particle.energy < MIN_ENERGY()) {
                    intRes.particleAlive = false;
                    intRes.energyImparted += particle.energy;
                } else {
                    if (particle.weight < interactions::russianRuletteWeightThreshold() && intRes.particleAlive) {
                        if (state.randomUniform() < interactions::russianRuletteProbability()) {
                            intRes.particleAlive = false;
                        } else {
                            constexpr auto factor = 1 / (1 - interactions::russianRuletteProbability());
                            particle.weight *= factor;
                            intRes.particleAlive = true;
                        }
                    }
                }
            } else {
                // real photoelectric effect event
                // we don't score energy since it's already done. But terminates particle to prevent bias.
                intRes.particleAlive = false;
                particle.energy = 0;
            }
        } else {
            // No interaction event, transport particle to border
            particle.border_translate(maxStepLen);
        }
        return intRes;
    }

    template <std::size_t Nshells, int LOWENERGYCORRECTION = 2>
    InteractionResult interactForced(double maxStepLenght, double density, ParticleType auto& particle, const Material<Nshells>& material, RandomState& state)
    {
        const auto attenuation = material.attenuationValues(particle.energy);
        return interactForced<Nshells, LOWENERGYCORRECTION>(maxStepLenght, density, attenuation, particle, material, state);
    }
}
}

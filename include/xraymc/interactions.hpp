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

/**
 * @file interactions.hpp
 * @brief Photon interaction samplers for analog and forced Monte Carlo transport.
 *
 * Implements the three photon interaction processes needed for diagnostic X-ray
 * Monte Carlo simulation:
 * - **Rayleigh (coherent) scattering** — `rayleightScatter()`
 * - **Compton (incoherent) scattering** — `comptonScatter()` / `comptonScatterIA()`
 * - **Photoelectric absorption** — `photoelectricEffect()` / `photoelectricEffectIA()`
 *
 * Two transport modes are provided:
 * - **Analog** — `interact()`: samples one interaction per call from the full cross section.
 * - **Forced** — `interactForced()`: scores photoelectric energy analytically and samples
 *   scatter events stochastically for variance reduction.
 *
 * All samplers share the `LOWENERGYCORRECTION` compile-time flag (see namespace doc).
 */

#pragma once

#include "xraymc/constants.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

namespace xraymc {

/**
 * @namespace xraymc::interactions
 * @brief Free functions implementing photon interaction physics for Monte Carlo transport.
 *
 * All public samplers are function templates parameterised on:
 * - `Nshells` — number of electron shells in the `Material` model (compile-time constant).
 * - `LOWENERGYCORRECTION` — selects the physics model complexity (see table below).
 * - `P` — particle type satisfying `ParticleType` (default `Particle`; use `ParticleTrack`
 *   to record interaction positions automatically).
 *
 * ### `LOWENERGYCORRECTION` flag
 *
 * | Value | Compton                                | Rayleigh                        | Photoelectric              |
 * |-------|----------------------------------------|---------------------------------|----------------------------|
 * | `0`   | Klein–Nishina only                     | Dipole distribution             | Full absorption            |
 * | `1`   | Klein–Nishina + incoherent scatter factor (Livermore) | Form-factor corrected | Full absorption |
 * | `2`   | Full Relativistic Impulse Approximation | Form-factor corrected | Shell-resolved + fluorescence |
 *
 * The default throughout the library is `2` (highest fidelity).
 *
 * ### Transport flow
 *
 * **Analog** (`interact`):
 * ```
 * interact()
 * ├── photoelectricEffect()  →  photoelectricEffectIA()
 * ├── comptonScatter()       →  comptonScatterIA()  or  Klein-Nishina
 * └── rayleightScatter()
 * ```
 * After each interaction, Russian roulette is applied if the particle weight
 * drops below `russianRuletteWeightThreshold()`.
 *
 * **Forced** (`interactForced`):
 * ```
 * interactForced()
 * ├── Photoelectric energy scored analytically (Beer-Lambert expectation)
 * ├── comptonScatter()   (stochastic, at sampled step position)
 * └── rayleightScatter() (stochastic, at sampled step position)
 * ```
 * The forced variant eliminates photoelectric variance at the cost of always
 * translating the particle to a scatter point within the step.
 *
 * ### Return values
 * All interaction functions return the energy imparted to the medium in **keV**,
 * already multiplied by the particle weight. The particle is updated in place.
 * High-level callers receive an `InteractionResult` from `interact()` /
 * `interactForced()` which also carries flags for each interaction type.
 */
namespace interactions {

    /**
     * @brief Survival probability used in Russian roulette weight reduction.
     *
     * When a particle's weight falls below `russianRuletteWeightThreshold()`, a
     * Russian roulette game is played: the particle survives with this probability
     * and its weight is boosted by 1/(1 − probability), or it is killed otherwise.
     *
     * @return The survival probability (0.9).
     */
    constexpr static double russianRuletteProbability()
    {
        return 0.9;
    }

    /**
     * @brief Weight threshold below which Russian roulette is applied.
     *
     * If a particle's statistical weight drops below this value after an interaction,
     * Russian roulette is invoked to either terminate the particle or restore its
     * weight to an unbiased level.
     *
     * @return The weight threshold (0.1).
     */
    constexpr static double russianRuletteWeightThreshold()
    {
        return 0.1;
    }

    /**
     * @brief Samples a Rayleigh (coherent) scattering deflection and updates the particle direction.
     *
     * Two sampling modes are selected at compile time via `LOWENERGYCORRECTION`:
     * - `0`: Classical dipole distribution — rejection samples θ from the (2 − sin²θ)sinθ
     *   distribution without form-factor correction.
     * - `1` or `2` (default): Form-factor–corrected sampling following the EGSnrc method
     *   (SLAC-730). The squared momentum transfer q² is sampled from the squared atomic
     *   form factor, and the scattering angle is accepted/rejected against the dipole term.
     *
     * If `P` is `ParticleTrack`, the current position is recorded before scattering.
     *
     * @tparam Nshells           Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag selecting the sampling algorithm (default 2).
     * @tparam P                 Particle type satisfying `ParticleType` (default `Particle`).
     * @param particle  The photon to scatter; its direction is updated in place.
     * @param material  The material providing form-factor and momentum-transfer data.
     * @param state     The PRNG state used for random sampling.
     */
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

    /**
     * @brief Samples the electron shell index for an Impulse Approximation Compton scatter event.
     *
     * Selects a shell stochastically, weighted by each shell's fractional electron occupancy,
     * skipping shells whose binding energy exceeds the photon energy (no interaction possible).
     *
     * @tparam Nshells Number of electron shells in the material model.
     * @param particle  The photon; only its energy is read.
     * @param material  The material providing shell binding energies and electron fractions.
     * @param state     The PRNG state used for random sampling.
     * @return Zero-based index of the sampled shell.
     */
    template <std::size_t Nshells>
    int comptonScatterIA_NRC_sample_shell(const ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        const auto& shells = material.shells();
        int shell = 0;
        double min_acc = 0;
        while (particle.energy < shells[shell].bindingEnergy) {
            min_acc += shells[shell].numberOfElectronsFraction;
            ++shell;
        }
        const auto r = state.randomUniform(min_acc, 1.0);
        do {
            min_acc += shells[shell].numberOfElectronsFraction;
            ++shell;
        } while (r > min_acc);

        return shell - 1;
    }

    /**
     * @brief Samples an Impulse Approximation (IA) incoherent (Compton) scatter event.
     *
     * Implements the NRC/EGSnrc Compton IA model (see PIRS-701). For each candidate
     * scatter angle:
     *   1. A shell is selected via `comptonScatterIA_NRC_sample_shell`.
     *   2. The Klein–Nishina scattered photon energy ratio `e` and cos(θ) are sampled.
     *   3. The Compton profile S(p_z) for the selected shell is evaluated and used as
     *      an acceptance probability (first-order approximation; second-order is
     *      compiled out by default).
     *   4. The longitudinal momentum component p_z is sampled from the shell's
     *      Hartree–Fock orbital, and the corrected scattered energy k̄ is computed.
     *
     * The particle's energy and direction are updated in place.  For high-Z / low-energy
     * photons where the non-relativistic p_z treatment breaks down, the scattered energy
     * is clamped to [0, E_initial].
     *
     * @tparam Nshells Number of electron shells in the material model.
     * @param particle  The photon to scatter; energy and direction are updated in place.
     * @param material  The material providing shell data and Hartree–Fock orbitals.
     * @param state     The PRNG state used for random sampling.
     * @return Energy imparted to the medium in keV (weighted: (E_i − E_f) × weight).
     */
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
            alpha = std::max(0.0, qc / k + kc * (kc - k * cosTheta) / (k * qc)); // alpha must be >= 0
            const auto b_part = (1 + 2 * J0 * std::abs(pi));
            expb = std::exp(-(0.5 * b_part * b_part - 0.5));

            p = std::sqrt((U + 1) * (U + 1) - 1);

            constexpr bool USE_SECOND_ORDER_APPROXIMATION = false;
            if constexpr (USE_SECOND_ORDER_APPROXIMATION) {
                auto Erf = [](const auto expb_arg, const auto J0_arg, const auto pi_arg) -> double {
                    // approximation of the error function
                    constexpr double sqrt_e = 1.648721270;
                    const auto t = 1.0 / (1.0 + 0332673 * (1 + 2 * J0_arg * std::abs(pi_arg)));
                    constexpr auto a1 = 0.34802;
                    constexpr auto a2 = -0.0958798;
                    constexpr auto a3 = 0.7478556;
                    return sqrt_e - expb_arg * t * (a1 + a2 * t + a3 * t * t);
                };
                if (pi < -p) {
                    S = (1 - alpha * p) * expb * 0.5;
                } else if (pi <= 0) {
                    const auto S1 = (1 + alpha * pi) * expb * 0.5;
                    constexpr double pik = 1.0 / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                    const auto S2 = -alpha / (4 * J0) * pik;
                    const auto errf1 = Erf(expb, J0, p);
                    const auto errf2 = Erf(expb, J0, pi);
                    S = S1 + S2 * (errf1 - errf2);
                } else if (pi < p) {
                    const auto S1 = (1 + alpha * pi) * expb * 0.5;
                    constexpr double pik = 1.0 / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                    const auto S2 = -alpha / (4 * J0) * pik;
                    const auto errf1 = Erf(expb, J0, p);
                    const auto errf2 = Erf(expb, J0, pi);
                    S = 1 - S1 + S2 * (errf1 - errf2);
                } else {
                    S = 1 - (1 + alpha * p) * expb * 0.5;
                }
            } else {
                if (pi < -p) {
                    S = (1 - alpha * p) * expb * 0.5;
                } else if (pi <= 0) {
                    S = (1 + alpha * pi) * expb * 0.5;
                } else if (pi < p) {
                    S = 1 - (1 + alpha * pi) * expb * 0.5;
                } else {
                    S = 1 - (1 + alpha * p) * expb * 0.5;
                }
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
            } while (state.randomUniform() > F / Fmax);
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

        // calculating kbar
        const auto bar_part = 1 - 2 * e * cosTheta + e * e * (1 - pz * pz * (1 - cosTheta * cosTheta));
        const auto kbar = kc * (1 - pz * pz * e * cosTheta + pz * std::sqrt(bar_part)) / (1 - pz * pz * e * e);

        // For high Z materials and low energy photons, we can have |pz| > 1. This is not
        // valid for our non relativistic treatment of pz and tehre is a small chance of
        // ketting a negative kbar. Hence we limit the photon energy
        particle.energy = std::clamp(kbar * ELECTRON_REST_MASS(), 0.0, E);

        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);

        return (E - particle.energy) * particle.weight;
    }

    /**
     * @brief Samples an incoherent (Compton) scatter event and updates the particle.
     *
     * Dispatches to one of three sampling algorithms based on `LOWENERGYCORRECTION`:
     * - `2` (default): Full Impulse Approximation via `comptonScatterIA` (NRC/EGSnrc model).
     * - `1`: Livermore model — Klein–Nishina sampling with incoherent scatter-factor
     *   rejection using `material.scatterFactor()`.
     * - `0`: Simple Klein–Nishina sampling without any scatter-factor correction.
     *
     * If `P` is `ParticleTrack`, the current position is recorded before scattering.
     *
     * @tparam Nshells           Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag selecting the sampling algorithm (default 2).
     * @tparam P                 Particle type satisfying `ParticleType` (default `Particle`).
     * @param particle  The photon to scatter; energy and direction are updated in place.
     * @param material  The material providing attenuation and scatter-factor data.
     * @param state     The PRNG state used for random sampling.
     * @return Energy imparted to the medium in keV (weighted: (E_i − E_f) × weight).
     */
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

    /**
     * @brief Samples a photoelectric absorption event with shell-resolved fluorescence (IA model).
     *
     * Selects the absorbing shell stochastically, weighted by each shell's partial
     * photoelectric cross section relative to @p totalPhotoCrossSection. If the selected
     * shell has a mean fluorescence photon energy above `MIN_ENERGY()`, a characteristic
     * X-ray is emitted isotropically with the shell's fluorescence yield and energy;
     * otherwise the photon is fully absorbed. The particle energy is set to the
     * fluorescence photon energy (or zero on full absorption).
     *
     * @tparam Nshells Number of electron shells in the material model.
     * @param totalPhotoCrossSection Total photoelectric cross section at the photon energy (cm²/g).
     * @param particle  The photon; energy is set to the fluorescence photon energy or zero.
     * @param material  The material providing shell cross sections and fluorescence data.
     * @param state     The PRNG state used for random sampling.
     * @return Energy imparted to the medium in keV (weighted by particle weight).
     */
    template <std::size_t Nshells>
    auto photoelectricEffectIA(const double totalPhotoCrossSection, ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // finding shell based on photoelectric cross section

        std::uint_fast8_t shell = 0;
        auto prob = state.randomUniform();
        do {
            const auto& sh = material.shell(shell);
            if (sh.bindingEnergy > particle.energy) {
                shell++;
            } else {
                const auto shellCS = material.attenuationPhotoelectricShell(shell, particle.energy);
                const auto shellProb = shellCS / totalPhotoCrossSection;
                prob -= shellProb;
                if (prob > 0.0)
                    ++shell;
            }
        } while (prob > 0.0 && shell != material.numberOfShells() - 1);
        if (prob > 0)
            ++shell;

        auto E = particle.energy * particle.weight;
        particle.energy = 0;
        if (shell != material.numberOfShells()) {
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

    /**
     * @brief Samples a photoelectric absorption event and updates the particle.
     *
     * Dispatches based on `LOWENERGYCORRECTION`:
     * - `2` (default): Full IA model with shell selection and fluorescence via
     *   `photoelectricEffectIA`.
     * - `0` or `1`: Simple full absorption — the photon energy is deposited entirely
     *   and the particle energy is set to zero.
     *
     * If `P` is `ParticleTrack`, the current position is recorded before the interaction.
     *
     * @tparam Nshells           Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag selecting the sampling algorithm (default 2).
     * @tparam P                 Particle type satisfying `ParticleType` (default `Particle`).
     * @param totalPhotoCrossSection Total photoelectric cross section at the photon energy (cm²/g).
     * @param particle  The photon; energy is set to fluorescence photon energy or zero.
     * @param material  The material providing shell and fluorescence data.
     * @param state     The PRNG state used for random sampling.
     * @return Energy imparted to the medium in keV (weighted by particle weight).
     */
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

    /**
     * @brief Outcome of a single photon interaction step.
     *
     * Returned by `interact()` and `interactForced()`. Contains the energy deposited
     * in the current step and flags describing what happened to the particle.
     */
    struct InteractionResult {
        // the size of InteractionResult is 16 bytes anyway, why not throw in
        // type of interaction, even it's almost not used
        double energyImparted = 0; ///< Energy deposited in the medium this step (keV, weighted).
        bool particleAlive = true; ///< False if the particle was absorbed or killed by Russian roulette.
        bool particleEnergyChanged = false; ///< True if the photon energy changed (photoelectric or Compton).
        bool particleDirectionChanged = false; ///< True if the photon direction changed (any interaction).
        bool interactionWasPhotoelectric = false; ///< True if the interaction was photoelectric absorption.
        bool interactionWasCoherent = false; ///< True if the interaction was Rayleigh (coherent) scatter.
        bool interactionWasIncoherent = false; ///< True if the interaction was Compton (incoherent) scatter.
    };

    /**
     * @brief Samples one analog photon interaction and updates the particle state.
     *
     * Selects the interaction type — photoelectric, incoherent (Compton), or coherent
     * (Rayleigh) — by comparing a uniform random number against the partial cross
     * sections in @p attenuation. After the interaction, Russian roulette is applied
     * if the particle weight falls below `russianRuletteWeightThreshold()`.
     *
     * @tparam Nshells           Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag forwarded to the individual interaction
     *                            samplers (default 2 = full IA model).
     * @param attenuation Pre-computed attenuation values (photoelectric, incoherent, coherent)
     *                    at the particle energy (cm²/g).
     * @param particle    The photon to interact; energy and direction updated in place.
     * @param material    The material providing cross-section and shell data.
     * @param state       The PRNG state used for random sampling.
     * @return An `InteractionResult` describing the deposited energy and interaction type.
     */
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
            res.interactionWasPhotoelectric = true;
        } else if (r2 < (attenuation.photoelectric + attenuation.incoherent)) {
            const auto Ei = interactions::comptonScatter<Nshells, LOWENERGYCORRECTION>(particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
            res.particleDirectionChanged = true;
            res.interactionWasIncoherent = true;
        } else {
            interactions::rayleightScatter<Nshells, LOWENERGYCORRECTION>(particle, material, state);
            res.particleDirectionChanged = true;
            res.interactionWasCoherent = true;
        }
        if (particle.energy < MIN_ENERGY()) {
            res.particleAlive = false;
            res.energyImparted += particle.energy * particle.weight;
            particle.energy = 0;
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

    /**
     * @brief Samples a forced-interaction step for variance reduction.
     *
     * Implements the interaction forcing technique: the photoelectric energy deposition
     * is scored analytically as the expected contribution over the step, while a scatter
     * event (Compton or Rayleigh) is sampled stochastically at a random point within the
     * step. If no stochastic event occurs the particle is transported to the far boundary.
     *
     * The forced photoelectric contribution is:
     *   E_imparted += E × weight × (1 − exp(−μ_total × ρ × L)) × (μ_pe / μ_total)
     *
     * Russian roulette is applied after any stochastic scatter if the weight is low.
     *
     * @tparam NMaterialShells   Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag forwarded to interaction samplers (default 2).
     * @param maxStepLen      Length of the transport step (cm).
     * @param materialDensity Density of the current material (g/cm³).
     * @param attenuation     Pre-computed attenuation values at the particle energy (cm²/g).
     * @param particle        The photon; position, energy, and direction updated in place.
     * @param material        The material providing cross-section and shell data.
     * @param state           The PRNG state used for random sampling.
     * @return An `InteractionResult` describing the deposited energy and interaction type.
     */
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
                    intRes.interactionWasIncoherent = true;
                } else {
                    interactions::rayleightScatter<NMaterialShells, LOWENERGYCORRECTION>(particle, material, state);
                    intRes.particleDirectionChanged = true;
                    intRes.interactionWasCoherent = true;
                }

                // Handle cutoff values
                if (particle.energy < MIN_ENERGY()) {
                    intRes.particleAlive = false;
                    intRes.energyImparted += particle.energy * particle.weight;
                    particle.energy = 0;
                } else {
                    if (particle.weight < interactions::russianRuletteWeightThreshold()) {
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
                // we don't score energy since it's already done, else do as usual to prevent bias.
                const auto Ei = interactions::photoelectricEffect<NMaterialShells, LOWENERGYCORRECTION>(attenuation.photoelectric, particle, material, state);
                // In case of characteristic radiation
                intRes.particleEnergyChanged = true;
                intRes.particleDirectionChanged = true;
                intRes.interactionWasPhotoelectric = true;

                if (particle.energy < MIN_ENERGY()) {
                    intRes.particleAlive = false;
                    intRes.energyImparted += particle.energy * particle.weight;
                    particle.energy = 0;
                } else {
                    if (particle.weight < interactions::russianRuletteWeightThreshold()) {
                        if (state.randomUniform() < interactions::russianRuletteProbability()) {
                            intRes.particleAlive = false;
                        } else {
                            constexpr auto factor = 1 / (1 - interactions::russianRuletteProbability());
                            particle.weight *= factor;
                            intRes.particleAlive = true;
                        }
                    }
                }
            }
        } else {
            // No interaction event, transport particle to border
            particle.border_translate(maxStepLen);
        }
        return intRes;
    }

    /**
     * @brief Convenience overload of `interactForced` that computes attenuation internally.
     *
     * Looks up the attenuation values for the particle energy from @p material and
     * delegates to the primary `interactForced` overload.
     *
     * @tparam Nshells           Number of electron shells in the material model.
     * @tparam LOWENERGYCORRECTION Compile-time flag forwarded to interaction samplers (default 2).
     * @param maxStepLenght   Length of the transport step (cm).
     * @param density         Density of the current material (g/cm³).
     * @param particle        The photon; position, energy, and direction updated in place.
     * @param material        The material providing cross-section and shell data.
     * @param state           The PRNG state used for random sampling.
     * @return An `InteractionResult` describing the deposited energy and interaction type.
     */
    template <std::size_t Nshells, int LOWENERGYCORRECTION = 2>
    InteractionResult interactForced(double maxStepLenght, double density, ParticleType auto& particle, const Material<Nshells>& material, RandomState& state)
    {
        const auto attenuation = material.attenuationValues(particle.energy);
        return interactForced<Nshells, LOWENERGYCORRECTION>(maxStepLenght, density, attenuation, particle, material, state);
    }
}
}

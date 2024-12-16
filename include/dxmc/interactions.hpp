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
    double comptonScatterIA_NRC2(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {

        // sample shell_idx
        int shell_idx = 0;
        {
            const double rand = state.randomUniform();
            double acc = material.shells()[0].numberOfElectronsFraction;

            while (acc < rand && shell_idx < Nshells) {

                acc += material.shells()[shell_idx].numberOfElectronsFraction;
                ++shell_idx;
            }
        }

        // calculate scatter angle
        const auto k = particle.energy / ELECTRON_REST_MASS();
        double cosTheta;
        {

            const auto emin = 1 / (1 + 2 * k);
            const auto gmaxInv = emin / (1 + emin * emin);
            bool rejected;
            do {
                const auto r1 = state.randomUniform();
                const auto e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;

                rejected = state.randomUniform() > g;

            } while (rejected);
        }

        // calculate upper limit for electron momentum p_i
        const auto U = material.shells()[shell_idx].bindingEnergy / ELECTRON_REST_MASS();
        const double p_i = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);

        // calculate S_i, scatter function for shell i
        const auto p = std::sqrt((U + 1) * (U + 1) - 1);
        const auto kc = k / (1 + k * (1 - cosTheta));
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        const auto alpha = (qc / k) * (1 + (kc * kc - kc * k * cosTheta) / (qc * qc));

        const auto J = material.shells()[shell_idx].HartreeFockOrbital_0;
        const auto b_part = 1 + 2 * J * std::abs(p_i);
        const auto b = 0.5 * (b_part * b_part) - 0.5;

        double S_i;
        {

            if (p_i <= -p) {
                S_i = (1 - alpha * p) * std::exp(-b) * 0.5;
            } else if (p_i <= 0) {
                const auto t = 1 / (1 + 0.332673 * b_part);
                constexpr auto a1 = 0.34802;
                constexpr auto a2 = -0.0958798;
                constexpr auto a3 = 0.7478556;
                const auto E = std::exp(-0.5) - std::exp(-b) * t * (a1 + a2 * t + a3 * t * t);
                const auto tp = 1 / (1 + 0.332673 * (1 + 2 * J * p));
                const auto Ep = std::exp(-0.5) - std::exp(-(1 + 2 * J * p)) * tp * (a1 + a2 * tp + a3 * tp * tp);
                S_i = (1 - alpha * p) * std::exp(-b) * 0.5 - alpha * std::sqrt(PI_VAL() / 2) * (Ep - E) / (4 * J);
            } else if (p_i < p) {
                const auto t = 1 / (1 + 0.332673 * b_part);
                constexpr auto a1 = 0.34802;
                constexpr auto a2 = -0.0958798;
                constexpr auto a3 = 0.7478556;
                const auto E = std::exp(-0.5) - std::exp(-b) * t * (a1 + a2 * t + a3 * t * t);
                const auto tp = 1 / (1 + 0.332673 * (1 + 2 * J * p));
                const auto Ep = std::exp(-0.5) - std::exp(-(1 + 2 * J * p)) * tp * (a1 + a2 * tp + a3 * tp * tp);
                S_i = 1 - (1 - alpha * p) * std::exp(-b) * 0.5 - alpha * std::sqrt(PI_VAL() / 2) * (Ep - E) / (4 * J);
            } else {
                S_i = 1 - (1 - alpha * p) * std::exp(-b) * 0.5;
            }
        }

        // sampling pz
        double pz;
        {
            double ratio, F, Fmax;
            const auto expnb = std::exp(-b);
            do {
                const auto r = expnb * state.randomUniform();
                if (r < 0.5) {
                    pz = (1 - std::sqrt(1 - 2 * std::log(2 * r))) / (2 * J);
                } else {
                    pz = (std::sqrt(1 - 2 * std::log(2 - 2 * r)) - 1) / (2 * J);
                }

                if (pz < -p) {
                    F = 1 - alpha * p;
                } else if (pz < p) {
                    F = 1 + alpha * pz;
                } else {
                    F = 1 + alpha * p;
                }

                if (p_i < -p) {
                    Fmax = 1 - alpha * p;
                } else if (p_i < p) {
                    Fmax = 1 + alpha * p_i;
                } else {
                    Fmax = 1 + alpha * p;
                }

                ratio = F / Fmax;

            } while (state.randomUniform() > ratio);
        }

        // calculate k_post
        const auto eps = kc / k;
        const auto k_post = kc * (1 - pz * pz * eps * cosTheta + pz * std::sqrt(1 - 2 * eps * cosTheta + eps * eps * (1 - pz * pz * (1 - cosTheta * cosTheta)))) / (1 - pz * pz * eps * eps);
        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);

        const auto E = particle.energy;
        particle.energy = k_post * ELECTRON_REST_MASS();
        return (E - particle.energy) * particle.weight;
    }

    template <std::size_t Nshells>
    double comptonScatterIA_NRC3(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {

        // analytical compton profiles
        const auto nia_func = [](double pz, double J0) -> double {
            constexpr auto d2 = std::numbers::sqrt2;
            constexpr auto d1 = 1 / d2;
            double part = d1 - d2 * J0 * pz;
            double n_neg = 0.5 * std::exp(d1 * d1 - part * part);
            return pz < 0.0 ? n_neg : 1 - n_neg;
        };

        // sample pz from analytical profiles
        const auto pz_func = [](double A, double J0) -> double {
            constexpr auto d2 = std::numbers::sqrt2;
            constexpr auto d1 = 1 / d2;
            if (A < 0.5) {
                return (d1 - std::sqrt(d1 * d1 - std::numbers::ln2 * A)) / (d2 * J0);
            } else {
                return (std::sqrt(d1 * d1 - std::numbers::ln2 * (1 - A)) - d1) / (d2 * J0);
            }
        };

        // maximum momentum transfer
        const auto pi_max_func = [](double k, double U, double cosTheta) {
            double part = k * (k - U) * (1 - cosTheta);
            return (part - U) / std::sqrt(2 * part + U * U);
        };

        // maximum momentum transfer at cosTheta = -1
        const auto pi_max_func_reverse = [](double k, double U) {
            double part = k * (k - U) * 2;
            return (part - U) / std::sqrt(2 * part + U * U);
        };

        const auto k = particle.energy / ELECTRON_REST_MASS();

        // precompute binding energiens in units of electron rest mass
        std::array<double, Nshells + 1> U, S_reverse, S;
        for (std::size_t i = 0; i < U.size(); ++i) {
            const auto& shell = material.shells()[i];
            U[i] = shell.bindingEnergy / ELECTRON_REST_MASS();
            if (U[i] <= k) {
                double pirev = pi_max_func(k, U[i], -1.0);
                S_reverse[i] = shell.numberOfElectrons * nia_func(pirev, shell.HartreeFockOrbital_0);
            } else {
                S_reverse[i] = 0;
            }
        }
        const double S_reverse_sum = std::reduce(std::execution::unseq, S_reverse.cbegin(), S_reverse.cend(), 0.0);

        // sample cosTheta
        double cosTheta, S_sum, kc, e;
        {
            const auto emin = 1 / (1 + 2 * k);
            const auto gmaxInv = emin / (1 + emin * emin);
            bool rejected;
            do {
                do {
                    const auto r1 = state.randomUniform();
                    e = r1 + (1 - r1) * emin;
                    const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                    cosTheta = 1 - t;
                    const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                    const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                    rejected = state.randomUniform() > g;
                } while (rejected);

                // we have an angle, lets sample p

                for (std::size_t i = 0; i < U.size(); ++i) {
                    const auto& shell = material.shells()[i];
                    if (U[i] <= k) {
                        double pi_max = pi_max_func(k, U[i], cosTheta);
                        S[i] = shell.numberOfElectrons * nia_func(pi_max, shell.HartreeFockOrbital_0);
                    } else {
                        S[i] = 0;
                    }
                }
                S_sum = std::reduce(std::execution::unseq, S.cbegin(), S.cend(), 0.0);

                const auto T = (S_sum / S_reverse_sum) * (1 - ((1 - e) * ((2 * k + 1) * e - 1) / (k * k * e * (1 + e * e))));
                // const auto T = (S_sum / S_reverse_sum);
                rejected = state.randomUniform() > T;

            } while (rejected);
            kc = k * e;
        }

        // sample shell
        double pz;
        {
            const double qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
            const double alpha = qc * (1 + kc * (kc - k * cosTheta) / (qc * qc)) / k;
            double Fmax = std::max(1 + alpha, 1 - alpha); // can be better, see egsnrc

            bool rejected;
            do {
                do {
                    // sample shell
                    const auto thres = state.randomUniform(S_sum);
                    int shell_idx = 0;
                    double accumulator = S[0];
                    while (accumulator < thres) {
                        shell_idx++;
                        accumulator += S[shell_idx];
                    }
                    // sample pz
                    const auto n_electrons = material.shells()[shell_idx].numberOfElectrons;
                    const auto HF_orb = material.shells()[shell_idx].HartreeFockOrbital_0;
                    const auto nia_pz = state.randomUniform() * S[shell_idx] / n_electrons;
                    pz = pz_func(nia_pz, HF_orb);
                    rejected = pz < -1;
                } while (rejected);
                const double F = 1 + alpha * pz;
                rejected = state.randomUniform() * Fmax > F;
            } while (rejected);
        }

        const double s = pz >= 0 ? 1 : -1;
        const double kbar_part = 1 - pz * pz * e * cosTheta;
        const double ebar = e * (kbar_part - s * std::sqrt(kbar_part * kbar_part - (1 - pz * pz * e * e) * (1 - pz * pz))) / (1 - pz * pz * e * e);

        const auto E = particle.energy;
        particle.energy *= ebar;
        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);
        auto Ei = (E - particle.energy) * particle.weight;
        return Ei;
    }

    template <std::size_t Nshells>
    int comptonScatterIA_NRC_sample_shell(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
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
        } while (shells[shell_idx].bindingEnergy > particle.energy && shell_idx != Nshells + 1);
        return shell_idx;
    }

    template <std::size_t Nshells>
    double comptonScatterIA_NRC(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        constexpr double fine_structure = 0.007297;
        const auto k = particle.energy / ELECTRON_REST_MASS();
        double S;
        double cosTheta;
        double e;
        double alpha;
        double p;
        double expb;
        double pi;
        double J0;
        do {
            // sample cosTheta
            const auto shell_idx = comptonScatterIA_NRC_sample_shell(particle, material, state);
            const auto e_min = 1 / (1 + 2 * k);
            const auto g_max = 1 / e_min + e_min;
            double g;
            do {
                const double r1 = state.randomUniform();
                e = r1 + (1 - r1) * e_min;
                cosTheta = 1 - (1 - e) / (e * k);
                const double sinThetaSqr = 1 - cosTheta * cosTheta;
                g = 1 / e + e - sinThetaSqr;
            } while (state.randomUniform(g_max) > g);

            // calculate pz max value; pi
            const auto U = material.shells()[shell_idx].bindingEnergy / ELECTRON_REST_MASS();
            pi = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);

            // Calculate S
            J0 = material.shells()[shell_idx].HartreeFockOrbital_0 / fine_structure;
            const auto kc = k * e;
            const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
            alpha = qc / k + kc * (kc - k * cosTheta) / (k * qc);
            const auto b_part = (1 + 2 * J0 * std::abs(pi));
            expb = std::exp(-(0.5 * b_part * b_part - 0.5));

            p = std::sqrt((U + 1) * (U + 1) - 1);

            if (pi < -p) {
                S = (1 - 2 * alpha * p) * expb * 0.5;
            } else if (pi <= 0) {
                const double pik = std::exp(0.5) / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                const auto S1 = (1 - 2 * alpha * pi) * expb * 0.5;
                const auto S2 = -alpha / (4 * J0) * pik;
                const auto errf1 = std::erf((1 + 2 * J0 * p) / std::numbers::sqrt2);
                const auto errf2 = std::erf((1 + 2 * J0 * std::abs(pi)) / std::numbers::sqrt2);
                S = S1 + S2 * (errf1 - errf2);
            } else if (pi < p) {
                const double pik = std::exp(0.5) / (std::numbers::sqrt2 * std::numbers::inv_sqrtpi);
                const auto S1 = (1 - 2 * alpha * pi) * expb * 0.5;
                const auto S2 = -alpha / (4 * J0) * pik;
                const auto errf1 = std::erf((1 + 2 * J0 * p) / std::numbers::sqrt2);
                const auto errf2 = std::erf((1 + 2 * J0 * std::abs(pi)) / std::numbers::sqrt2);
                S = 1 - S1 + S2 * (errf1 - errf2);
            } else {
                S = 1 - (1 + alpha * p) * expb * 0.5;
            }
        } while (state.randomUniform() > S);

        // sample pz
        double pz;

        double Fmax, F;
        if (pi < -p) {
            Fmax = 1 - alpha * p;
        } else if (pi < p) {
            Fmax = 1 + alpha * pi;
        } else {
            Fmax = 1 + alpha * p;
        }
        do {
            const auto r = state.randomUniform(expb);
            if (r < 0.5) {
                pz = fine_structure * (1 - std::sqrt(1 - 2 * std::log(2 * r))) / (2 * J0);
            } else {
                pz = fine_structure * (std::sqrt(1 - 2 * std::log(2 * (1 - r))) - 1) / (2 * J0);
            }
            if (pz < -p) {
                F = 1 - alpha * p;
            } else if (pz < p) {
                F = 1 + alpha * pz;
            } else {
                F = 1 + alpha * p;
            }

        } while (state.randomUniform(Fmax) > F);
        // pz = 0;
        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);

        // calkulating kbar
        const auto ebar_part = 1 - 2 * e * cosTheta + e * e * (1 - pz * pz * (1 - cosTheta * cosTheta));
        const auto ebar = e * (1 - pz * pz * e * cosTheta + pz * std::sqrt(ebar_part)) / (1 - pz * pz * e * e);

        const auto E = particle.energy;
        particle.energy *= ebar;
        return (E - particle.energy) * particle.weight;
    }
    template <std::size_t Nshells>
    auto comptonScatterIA(ParticleType auto& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // Penelope model for compton scattering. Note spelling error in manual for sampling of pz.
        // In addition we use hartree scatter factors instead of integrating all compton profiles when sampling cosTheta
        const auto k = particle.energy / ELECTRON_REST_MASS();

        // sample cosangle and e
        double e, cosTheta;
        {
            const auto emin = 1 / (1 + 2 * k);
            bool rejected;
            const auto gmaxInv = emin / (1 + emin * emin);
            do {
                const auto r1 = state.randomUniform();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                const auto scatterFactor = material.scatterFactor(q);
                rejected = state.randomUniform(material.effectiveZ()) > (g * scatterFactor);
            } while (rejected);
        }
        const auto kc = k * e;
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        // sample shell and pz
        const auto nia = [](const auto pz, const auto J0) -> double {
            constexpr auto d1 = 1 / std::numbers::sqrt2_v<double>;
            constexpr auto d2 = std::numbers::sqrt2_v<double>;
            const auto p1 = d1 + d2 * J0 * std::abs(pz);
            const auto p2 = 0.5 * std::exp(std::max(0.5 - p1 * p1, 0.0));
            return pz > 0 ? 1 - p2 : p2;
        };
        std::uint_fast16_t shellIdx;
        double pz;
        {
            std::array<double, Nshells + 1> shell_nia;
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), shell_nia.begin(), [&nia, cosTheta, k](const auto& shell) -> double {
                const auto U = shell.bindingEnergy / ELECTRON_REST_MASS();
                if (U > k)
                    return 0.0;
                const auto pi_max = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);
                const auto ni = nia(pi_max, shell.HartreeFockOrbital_0);
                return ni;
            });
            std::array<double, Nshells + 1> shell_probs;
            std::transform(shell_nia.cbegin(), shell_nia.cend(), material.shells().cbegin(), shell_probs.begin(), [](const auto n, const auto& shell) -> double {
                return n * shell.numberOfElectrons;
            });
            const auto shell_probs_sum = std::reduce(std::execution::unseq, shell_probs.cbegin(), shell_probs.cend(), 0.0);
            bool rejected;
            do {
                do {
                    const auto shell_probs_r = state.randomUniform(shell_probs_sum);
                    shellIdx = 0;
                    auto accum = shell_probs[shellIdx];
                    while (shell_probs_r > accum) {
                        ++shellIdx;
                        accum += shell_probs[shellIdx];
                    }

                    const auto nia_shell = shell_nia[shellIdx];
                    const auto A = state.randomUniform(nia_shell);
                    constexpr auto d1 = 1 / std::numbers::sqrt2_v<double>;
                    constexpr auto d2 = std::numbers::sqrt2_v<double>;
                    const auto J0 = material.shell(shellIdx).HartreeFockOrbital_0;
                    if (A < 0.5) {
                        pz = (d1 - std::sqrt(0.5 - std::log(2 * A))) / (d2 * J0);
                    } else {
                        pz = (std::sqrt(0.5 - std::log(2 * (1 - A))) - d1) / (d2 * J0);
                    }
                } while (pz < -1);

                const auto F_eval = [k, kc, qc, cosTheta](const auto p) {
                    const auto pn = std::clamp(p, -0.2, 0.2);
                    const auto f = 1 + qc * pn * (1 + kc * (kc - k * cosTheta) / (qc * qc)) / k;
                    return std::max(f, 0.0);
                };
                const auto Fmax = std::max(F_eval(0.2), F_eval(-0.2));
                const auto F = F_eval(pz);
                rejected = state.randomUniform(Fmax) > F;
            } while (rejected);
        }
        const auto t = pz * pz;
        const auto eb1 = 1 - t * e * e;
        const auto eb2 = 1 - t * e * cosTheta;
        const auto sign = pz > 0 ? 1 : -1;

        const auto sqrt_factor = std::max(eb2 * eb2 - eb1 * (1 - t), 0.0); // may have small negative facors, fix

        const auto eb = e / eb1 * (eb2 + sign * std::sqrt(sqrt_factor));

        const auto E = particle.energy;
        particle.energy *= eb;
        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);
        auto Ei = (E - particle.energy) * particle.weight;
        return Ei;
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
            return comptonScatterIA_NRC(particle, material, state);
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
                const auto theta = state.randomUniform(PI_VAL());
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

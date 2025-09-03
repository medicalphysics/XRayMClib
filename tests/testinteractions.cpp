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

Copyright 2022 Erlend Andersen
*/

#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"

#include <iostream>
#include <thread>
#include <chrono>
#include <utility>

using TimeType = std::chrono::high_resolution_clock::time_point;
auto timeNow = []() { return std::chrono::high_resolution_clock::now(); };
auto duration = [](const auto a) { return std::chrono::duration_cast<std::chrono::nanoseconds>(a).count(); };




std::pair<double, std::map<std::size_t, double>> TG195_air()
{
    std::map<std::size_t, double> w;
    w[6] = 0.0124;
    w[7] = 75.5268;
    w[8] = 23.1781;
    w[18] = 1.2827;

    return std::make_pair(0.001205, w);
}

double testS()
{
   const dxmc::Particle particle = { .pos = { 0, 0, 0 }, .dir = { 0, 1, 0 }, .energy = 56.4, .weight = 1 };
    const auto [air_dens, air_comp] = TG195_air();
    const auto material = dxmc::Material<5>::byWeight(air_comp).value();
    dxmc::RandomState state;

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
    double pz;
    do {
        const auto k = particle.energy / dxmc::ELECTRON_REST_MASS();
        const auto E = particle.energy;
        do {
            // sample cosTheta
            shell_idx = dxmc::interactions::comptonScatterIA_NRC_sample_shell(particle, material, state);
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
            U = shell_idx < material.numberOfShells() - 1 ? material.shells()[shell_idx].bindingEnergy / dxmc::ELECTRON_REST_MASS() : 0.0;

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
                const double r = 0; // std::nextafter( state.randomUniform(S), 0.5);
                if (r < 0.5) {
                    pz = (1 - std::sqrt(1 - 2 * std::log(2 * r))) / (2 * J0);
                } else {
                    pz = (std::sqrt(1 - 2 * std::log(2 * (1 - r))) - 1) / (2 * J0);
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
            const double r = std::nextafter(0.0, 0.5);            //state.randomUniform(S);
            if (r < 0.5) {
                pz = (1 - std::sqrt(1 - 2 * std::log(2 * r))) / (2 * J0);
            } else {
                pz = (std::sqrt(1 - 2 * std::log(2 * (1 - r))) - 1) / (2 * J0);
            }
        }
    } while (pz >= -p && pz <= p);
    return S;
}

bool longTestComptonIA()
{

    constexpr std::size_t N = 1E7;
    const auto [air_dens, air_comp] = TG195_air();
    const auto m = dxmc::Material<5>::byWeight(air_comp).value();
    dxmc::RandomState state;

    bool res = true;
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p = { .pos = { 0, 0, 0 }, .dir = { 0, 1, 0 }, .energy = 56.4, .weight = 1 };
        int teller = 50;
        while (p.energy > 1.0 && teller > 0 && res) {
            auto E = dxmc::interactions::comptonScatterIA<5>(p, m, state);
            if (p.energy < 0.0 || p.energy > 56.4) {
                std::cout << p.energy << " " << E << std::endl;
                res = false;
            }
        }
    }
    return res;
}

bool longTestInteractions()
{
    constexpr std::size_t N = 1E8;
    const auto [air_dens, air_comp] = TG195_air();
    const auto m = dxmc::Material<5>::byWeight(air_comp).value();
    dxmc::RandomState state;
    dxmc::Particle p;
    dxmc::interactions::InteractionResult res;
    for (std::size_t i = 0; i < N; ++i) {
        p = { .pos = { 0, 0, 0 }, .dir = { 0, 1, 0 }, .energy = 56.4, .weight = 1 };
        bool alive = true;
        while (alive) {
            const auto att = m.attenuationValues(p.energy);
            res = dxmc::interactions::interact<5, 2>(att, p, m, state);
            alive = res.particleAlive;
            if (p.energy < 0 || p.energy > 56.4) {
                std::cout << p.energy << std::endl;
                alive = false;
            }
        }
    }
    return false;
}

bool longTestInteractionsThreaded()
{
    for (int ii = 0; ii < 6; ++ii) {
        auto n = std::thread::hardware_concurrency();
        std::vector<std::jthread> jobs;
        jobs.reserve(n - 1);
        auto t0 = timeNow();
        for (int i = 0; i < n - 1; ++i)
            jobs.emplace_back(longTestInteractions);
        longTestInteractions();

        for (auto& j : jobs)
            j.join();
        auto t1 = timeNow();
        std::cout << "Time " << duration(t1 - t0) << std::endl;
    }
    return false;
}

int main()
{
    bool success = true;
    success = success && longTestInteractionsThreaded();
    success = success && longTestComptonIA();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
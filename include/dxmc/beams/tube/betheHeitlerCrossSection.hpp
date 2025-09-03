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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"

#include <algorithm>
#include <cmath>
#include <execution>

namespace dxmc {
namespace BetheHeitlerCrossSection {

    // ELECTRON DIFFUSION CALCULATIONS IN TUNGSTEN
    static constexpr double SIMULATED_ENERGY()
    {
        return 100.0; // Elektron energy [keV] for MonteCarlo simulation of electron penetration in tungsten
    }

    static constexpr std::array<double, 5> THOMSONWIDDINGTONCONSTANT_T()
    {
        std::array<double, 5> res = { 50, 80, 100, 120, 150 };
        return res;
    }

    static constexpr std::array<double, 5> THOMSONWIDDINGTONCONSTANT_C()
    {
        std::array<double, 5> res = { 565, 710, 792, 865, 964 };
        return res;
    }
    // Electron density distributions:

    static constexpr std::array<double, 12> CP100_F_x()
    {
        std::array<double, 12> res = { 0, 0.965, 1.93, 2.895, 3.86, 4.825, 5.79, 7.72, 9.65, 11.58, 13.51, 15.44 }; // [mg/cm2]
        return res;
    }

    static constexpr std::array<double, 45> CP100_u()
    {
        std::array<double, 45> res = { 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.73, 0.75, 0.77, 0.79, 0.81, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.97, 0.99 };
        return res;
    }

    static constexpr std::array<double, 12> CP100_M_x()
    {
        std::array<double, 12> res = { 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8 }; // [mg/cm2]
        return res;
    }

    static constexpr auto CP100_F()
    {
        std::array<double, CP100_F_x().size() * CP100_u().size()> res = {
            0, 0.029, 0.04, 0.049, 0.06, 0.077, 0.096, 0.137, 0.186, 0.242, 0.381, 0.5,
            0, 0.025, 0.035, 0.046, 0.06, 0.079, 0.102, 0.155, 0.22, 0.31, 0.412, 0.61,
            0, 0.022, 0.034, 0.045, 0.064, 0.086, 0.114, 0.179, 0.248, 0.357, 0.501, 0.725,
            0, 0.021, 0.031, 0.045, 0.065, 0.096, 0.126, 0.201, 0.296, 0.431, 0.566, 0.821,
            0, 0.018, 0.03, 0.045, 0.066, 0.099, 0.14, 0.236, 0.354, 0.479, 0.645, 1.073,
            0, 0.017, 0.029, 0.045, 0.075, 0.112, 0.158, 0.264, 0.406, 0.573, 0.783, 1.101,
            0, 0.016, 0.027, 0.046, 0.079, 0.125, 0.174, 0.296, 0.45, 0.64, 0.936, 1.293,
            0, 0.015, 0.026, 0.049, 0.086, 0.136, 0.198, 0.34, 0.518, 0.745, 1.09, 1.413,
            0, 0.014, 0.027, 0.05, 0.094, 0.148, 0.223, 0.386, 0.594, 0.848, 1.109, 1.656,
            0, 0.013, 0.026, 0.052, 0.103, 0.168, 0.243, 0.429, 0.657, 0.917, 1.294, 1.582,
            0, 0.012, 0.026, 0.054, 0.111, 0.187, 0.273, 0.489, 0.732, 1.03, 1.398, 1.816,
            0, 0.011, 0.026, 0.058, 0.12, 0.206, 0.303, 0.542, 0.812, 1.15, 1.627, 1.853,
            0, 0.01, 0.027, 0.063, 0.133, 0.227, 0.336, 0.604, 0.906, 1.273, 1.672, 2.261,
            0, 0.011, 0.028, 0.067, 0.147, 0.25, 0.368, 0.64, 1.009, 1.381, 1.825, 2.302,
            0, 0.01, 0.027, 0.076, 0.159, 0.275, 0.413, 0.728, 1.071, 1.512, 1.989, 2.38,
            0, 0.01, 0.029, 0.081, 0.176, 0.308, 0.451, 0.799, 1.181, 1.646, 2.172, 2.481,
            0, 0.01, 0.03, 0.089, 0.193, 0.33, 0.499, 0.859, 1.268, 1.717, 2.197, 2.788,
            0, 0.01, 0.032, 0.096, 0.216, 0.373, 0.539, 0.929, 1.4, 1.894, 2.352, 2.628,
            0, 0.011, 0.035, 0.109, 0.234, 0.41, 0.607, 1.012, 1.497, 1.977, 2.379, 2.412,
            0, 0.011, 0.039, 0.118, 0.267, 0.45, 0.66, 1.116, 1.598, 2.051, 2.429, 2.793,
            0, 0.011, 0.04, 0.132, 0.297, 0.501, 0.73, 1.2, 1.722, 2.235, 2.575, 2.477,
            0, 0.012, 0.046, 0.149, 0.326, 0.553, 0.797, 1.308, 1.824, 2.304, 2.455, 2.307,
            0, 0.012, 0.052, 0.168, 0.368, 0.617, 0.872, 1.421, 1.925, 2.371, 2.612, 2.344,
            0, 0.014, 0.057, 0.193, 0.414, 0.672, 0.953, 1.531, 2.088, 2.451, 2.492, 1.922,
            0, 0.016, 0.067, 0.219, 0.46, 0.75, 1.042, 1.662, 2.194, 2.454, 2.337, 1.857,
            0, 0.017, 0.078, 0.245, 0.523, 0.834, 1.16, 1.774, 2.274, 2.483, 2.149, 1.394,
            0, 0.021, 0.09, 0.287, 0.59, 0.927, 1.283, 1.924, 2.407, 2.457, 2.013, 1.197,
            0, 0.022, 0.106, 0.328, 0.672, 1.04, 1.402, 2.065, 2.452, 2.375, 1.73, 0.894,
            0, 0.025, 0.125, 0.386, 0.762, 1.159, 1.565, 2.227, 2.552, 2.229, 1.372, 0.523,
            0, 0.029, 0.148, 0.454, 0.868, 1.297, 1.72, 2.379, 2.543, 2.068, 1.073, 0.298,
            0, 0.036, 0.18, 0.532, 0.998, 1.469, 1.906, 2.51, 2.514, 1.786, 0.732, 0.161,
            0, 0.042, 0.221, 0.638, 1.143, 1.661, 2.11, 2.677, 2.426, 1.452, 0.424, 0.05,
            0, 0.052, 0.275, 0.758, 1.336, 1.893, 2.37, 2.783, 2.336, 1.074, 0.197, 0.014,
            0, 0.062, 0.338, 0.914, 1.574, 2.159, 2.624, 2.905, 2.02, 0.678, 0.067, 0.009,
            0, 0.082, 0.43, 1.11, 1.837, 2.483, 2.944, 2.942, 1.6, 0.32, 0.008, 0.009,
            0, 0.104, 0.556, 1.367, 2.194, 2.881, 3.274, 2.842, 1.089, 0.084, 0.005, 0.028,
            0, 0.139, 0.737, 1.71, 2.668, 3.385, 3.658, 2.55, 0.529, 0.003, 0.001, 0.014,
            0, 0.195, 0.98, 2.174, 3.282, 3.996, 4.011, 1.912, 0.103, 0.002, 0.001, 0,
            0, 0.283, 1.356, 2.844, 4.149, 4.694, 4.161, 0.956, 0.001, 0, 0.002, 0.005,
            0, 0.43, 1.946, 3.889, 5.343, 5.38, 3.662, 0.092, 0, 0, 0, 0.009,
            0, 0.711, 2.97, 5.677, 6.98, 5.312, 1.734, 0, 0, 0, 0, 0,
            0, 1.291, 5.064, 8.943, 8.291, 2.196, 0, 0, 0, 0, 0, 0,
            0, 2.829, 10.642, 15.399, 2.419, 0, 0, 0, 0, 0, 0, 0,
            0, 8.298, 22.892, 0.201, 0, 0, 0, 0, 0, 0, 0, 0,
            50, 34.973, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        return res;
    }

    static constexpr auto CP100_M()
    {
        std::array<double, CP100_M_x().size() * CP100_u().size()> res = {
            0.109, 0.44, 0.528, 0.617, 0.706, 0.808, 0.9, 1.1, 1.335, 1.595, 1.916, 2.177,
            0.113, 0.425, 0.525, 0.616, 0.718, 0.816, 0.925, 1.152, 1.361, 1.642, 2.046, 2.571,
            0.121, 0.417, 0.521, 0.626, 0.733, 0.843, 0.958, 1.182, 1.44, 1.715, 1.995, 2.503,
            0.125, 0.416, 0.522, 0.64, 0.753, 0.871, 0.996, 1.215, 1.503, 1.743, 2.132, 2.535,
            0.138, 0.414, 0.536, 0.657, 0.777, 0.896, 1.021, 1.286, 1.534, 1.833, 2.132, 2.535,
            0.15, 0.422, 0.552, 0.674, 0.808, 0.926, 1.057, 1.329, 1.577, 1.898, 2.111, 2.6,
            0.162, 0.425, 0.565, 0.702, 0.832, 0.967, 1.099, 1.37, 1.626, 1.933, 2.305, 2.75,
            0.175, 0.439, 0.58, 0.727, 0.862, 0.993, 1.124, 1.412, 1.692, 2.003, 2.383, 2.65,
            0.191, 0.451, 0.599, 0.75, 0.898, 1.027, 1.159, 1.423, 1.708, 2.011, 2.291, 2.557,
            0.207, 0.465, 0.621, 0.776, 0.921, 1.066, 1.198, 1.48, 1.75, 2.055, 2.322, 2.406,
            0.224, 0.477, 0.634, 0.802, 0.953, 1.099, 1.245, 1.507, 1.792, 2.074, 2.331, 2.585,
            0.255, 0.49, 0.657, 0.828, 0.987, 1.128, 1.277, 1.534, 1.824, 2.083, 2.209, 2.424,
            0.273, 0.51, 0.687, 0.859, 1.023, 1.159, 1.297, 1.574, 1.836, 2.072, 2.318, 2.331,
            0.3, 0.532, 0.704, 0.887, 1.042, 1.202, 1.338, 1.607, 1.872, 2.072, 2.217, 2.288,
            0.327, 0.556, 0.73, 0.916, 1.083, 1.231, 1.374, 1.637, 1.899, 2.045, 2.205, 2.091,
            0.356, 0.574, 0.758, 0.948, 1.108, 1.259, 1.394, 1.645, 1.867, 2.007, 2.046, 2.012,
            0.392, 0.602, 0.787, 0.974, 1.139, 1.289, 1.425, 1.649, 1.866, 1.936, 1.916, 1.762,
            0.426, 0.63, 0.814, 1.003, 1.172, 1.322, 1.454, 1.659, 1.818, 1.877, 1.879, 1.59,
            0.457, 0.657, 0.851, 1.044, 1.208, 1.333, 1.47, 1.666, 1.812, 1.834, 1.699, 1.482,
            0.498, 0.692, 0.884, 1.08, 1.239, 1.379, 1.486, 1.689, 1.762, 1.78, 1.618, 1.25,
            0.54, 0.725, 0.924, 1.107, 1.27, 1.408, 1.509, 1.66, 1.729, 1.693, 1.512, 1.096,
            0.588, 0.759, 0.965, 1.157, 1.301, 1.433, 1.528, 1.661, 1.681, 1.568, 1.32, 0.888,
            0.64, 0.803, 1.007, 1.192, 1.343, 1.456, 1.546, 1.648, 1.617, 1.462, 1.15, 0.834,
            0.692, 0.839, 1.049, 1.232, 1.375, 1.475, 1.561, 1.63, 1.537, 1.302, 0.956, 0.695,
            0.748, 0.893, 1.095, 1.284, 1.415, 1.501, 1.571, 1.595, 1.476, 1.201, 0.819, 0.458,
            0.81, 0.935, 1.14, 1.321, 1.439, 1.533, 1.576, 1.541, 1.394, 1.044, 0.723, 0.308,
            0.889, 0.992, 1.205, 1.371, 1.481, 1.55, 1.577, 1.518, 1.297, 0.936, 0.508, 0.229,
            0.975, 1.061, 1.262, 1.426, 1.509, 1.566, 1.565, 1.449, 1.154, 0.771, 0.406, 0.132,
            1.052, 1.126, 1.322, 1.469, 1.55, 1.581, 1.555, 1.373, 1.003, 0.621, 0.223, 0.097,
            1.144, 1.195, 1.392, 1.519, 1.575, 1.58, 1.537, 1.296, 0.899, 0.45, 0.143, 0.036,
            1.273, 1.289, 1.477, 1.582, 1.605, 1.586, 1.496, 1.187, 0.762, 0.336, 0.106, 0.025,
            1.39, 1.379, 1.546, 1.627, 1.63, 1.582, 1.469, 1.091, 0.603, 0.208, 0.04, 0.014,
            1.517, 1.482, 1.642, 1.685, 1.656, 1.551, 1.401, 0.944, 0.434, 0.126, 0.009, 0.007,
            1.683, 1.601, 1.723, 1.743, 1.656, 1.532, 1.324, 0.796, 0.29, 0.051, 0.003, 0.029,
            1.842, 1.726, 1.828, 1.781, 1.66, 1.471, 1.227, 0.633, 0.162, 0.021, 0.004, 0.025,
            2.046, 1.873, 1.926, 1.824, 1.655, 1.397, 1.08, 0.457, 0.071, 0.001, 0.002, 0.011,
            2.262, 2.058, 2.032, 1.859, 1.616, 1.285, 0.927, 0.265, 0.018, 0.001, 0.002, 0.004,
            2.515, 2.24, 2.136, 1.881, 1.553, 1.125, 0.702, 0.115, 0.001, 0, 0.000, 0.000,
            2.786, 2.451, 2.239, 1.866, 1.42, 0.902, 0.442, 0.026, 0, 0, 0, 0.000,
            3.076, 2.659, 2.306, 1.79, 1.181, 0.603, 0.188, 0.001, 0, 0, 0, 0,
            3.395, 2.872, 2.332, 1.592, 0.813, 0.247, 0.025, 0, 0, 0, 0, 0,
            3.604, 3.035, 2.217, 1.168, 0.324, 0.021, 0, 0, 0, 0, 0, 0,
            3.625, 3.051, 1.73, 0.397, 0.009, 0, 0, 0, 0, 0, 0.000, 0.000,
            3.329, 2.534, 0.45, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2.578, 0.388, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        return res;
    }

    static double ThomsonWiddingtonRange(double T0)
    // T0: electron energy, returns Thomson widdington range [mg/cm2] for a electron with initial energy T0 in tungsten. Valid for range [50, 150] keV
    {
        return 0.0119 * std::pow(T0, 1.513);
    }

    static double ThomsonWiddingtonLaw(double x, double tubeVoltage)
    {
        const auto& twct_array = THOMSONWIDDINGTONCONSTANT_T();
        auto twct_begin = twct_array.cbegin();
        auto twct_end = twct_array.cend();

        auto idx1 = std::lower_bound(twct_begin, twct_end, tubeVoltage);
        if (auto n = std::distance(idx1, twct_end); n < 2) {
            idx1 = twct_end;
            std::advance(idx1, -3 + n);
        }
        auto idx2 = idx1;
        ++idx2;

        const auto& twcc_array = THOMSONWIDDINGTONCONSTANT_C();
        auto idy1 = twcc_array.cbegin();
        std::advance(idy1, std::distance(twct_begin, idx1));
        auto idy2 = idy1;
        ++idy2;

        const auto C = *idy1 + ((*idy2 - *idy1) / (*idx2 - *idx1)) * (tubeVoltage - *idx1);

        const auto twl = (tubeVoltage * tubeVoltage - C * x) / (tubeVoltage * tubeVoltage);
        return twl < 0 ? 0 : twl;
    }

    static double numberFractionF(double x, double tubeVoltage)
    {
        constexpr auto L = 1.753;
        return std::pow(ThomsonWiddingtonLaw(x, tubeVoltage), L);
    }

    static double numberFractionM(double x, double tubeVoltage)
    {
        constexpr auto K = 18.0;
        constexpr auto Bd = 0.584;
        constexpr auto B0 = 0.5;
        const auto exp = 1 - std::exp(-K * x / ThomsonWiddingtonRange(tubeVoltage));
        const auto F = Bd * exp;
        const auto B = B0 + (Bd - B0) * exp;
        return numberFractionF(x, tubeVoltage) * B * (F + 1) / (1 - B * F);
    }

    static double bilinearInterpolate(double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2, double x, double y)
    {
        const auto xf1 = ((x2 - x) / (x2 - x1));
        const auto xf2 = ((x - x1) / (x2 - x1));

        const auto r1 = xf1 * q11 + xf2 * q21;
        const auto r2 = xf1 * q12 + xf2 * q22;

        const auto invdiv = 1 / (y2 - y1);

        return ((y2 - y) * invdiv) * r1 + ((y - y1) * invdiv) * r2;
    }

    static double electronDensity_F(double uval, double xval)
    {
        const auto& cp_f_x = CP100_F_x();
        const auto x = std::clamp(xval, cp_f_x.front(), cp_f_x.back());

        const auto& cp_u = CP100_u();
        const auto u = std::clamp(uval, cp_u.front(), cp_u.back());

        auto idx1 = std::lower_bound(cp_f_x.cbegin(), cp_f_x.cend(), x);
        if (auto n = std::distance(idx1, cp_f_x.cend()); n < 2) {
            idx1 = cp_f_x.cend();
            std::advance(idx1, -3 + n);
        }
        auto idx2 = idx1;
        ++idx2;

        auto idu1 = std::lower_bound(cp_u.cbegin(), cp_u.cend(), u);
        if (auto n = std::distance(idu1, cp_u.cend()); n < 2) {
            idu1 = cp_u.end();
            std::advance(idu1, -3 + n);
        }
        auto idu2 = idu1;
        ++idu2;

        const auto& cp_f = CP100_F();
        auto idq = cp_f.cbegin();
        std::advance(idq, std::distance(cp_u.cbegin(), idu1) * cp_f_x.size() + std::distance(cp_f_x.cbegin(), idx1)); // moving iterator to first part of 2d function
        const auto q11 = *idq;
        std::advance(idq, 1);
        const auto q21 = *idq;
        std::advance(idq, cp_f_x.size());
        const auto q22 = *idq;
        std::advance(idq, -1);
        const auto q12 = *idq;
        return bilinearInterpolate(q11, q12, q21, q22, *idx1, *idx2, *idu1, *idu2, x, u);
    }

    static double electronDensity_M(double uval, double xval)
    {
        const auto& cp_m_x = CP100_M_x();
        const auto x = std::clamp(xval, cp_m_x.front(), cp_m_x.back());

        const auto& cp_u = CP100_u();
        const auto u = std::clamp(uval, cp_u.front(), cp_u.back());

        auto idx1 = std::lower_bound(cp_m_x.cbegin(), cp_m_x.cend(), x);
        if (auto n = std::distance(idx1, cp_m_x.cend()); n < 2) {
            idx1 = cp_m_x.cend();
            std::advance(idx1, -3 + n);
        }
        auto idx2 = idx1;
        ++idx2;

        auto idu1 = std::lower_bound(cp_u.cbegin(), cp_u.cend(), u);
        if (auto n = std::distance(idu1, cp_u.cend()); n < 2) {
            idu1 = cp_u.cend();
            std::advance(idu1, -3 + n);
        }
        auto idu2 = idu1;
        ++idu2;

        const auto& cp_m = CP100_M();
        auto idq = cp_m.cbegin();
        std::advance(idq, std::distance(cp_u.cbegin(), idu1) * cp_m_x.size() + std::distance(cp_m_x.cbegin(), idx1)); // moving iterator to first part of 2d function
        const auto q11 = *idq;
        std::advance(idq, 1);
        const auto q21 = *idq;
        std::advance(idq, cp_m_x.size());
        const auto q22 = *idq;
        std::advance(idq, -1);
        const auto q12 = *idq;
        return bilinearInterpolate(q11, q12, q21, q22, *idx1, *idx2, *idu1, *idu2, x, u);
    }

    static double electronDensity(double u, double x, double tubeVoltage)
    {
        const auto f = ThomsonWiddingtonRange(SIMULATED_ENERGY()) / ThomsonWiddingtonRange(tubeVoltage);
        return numberFractionF(x, tubeVoltage) * electronDensity_F(u, x * f) + numberFractionM(x, tubeVoltage) * electronDensity_M(u, x * f);
    }

    // END ELECTRON DIFFUSION CALCULATIONS IN TUNGSTEN END

    static constexpr std::size_t TUNGSTEN_ATOMIC_NUMBER = 74;

    static constexpr double FINE_STRUCTURE_CONSTANT()
    {
        return 7.29735308E-03;
    }

    static constexpr double CLASSIC_ELECTRON_RADIUS()
    {
        return 2.81794092E-15; // [m]
    }

    static constexpr double PHI_BAR()
    {
        return (TUNGSTEN_ATOMIC_NUMBER * TUNGSTEN_ATOMIC_NUMBER) * CLASSIC_ELECTRON_RADIUS() * CLASSIC_ELECTRON_RADIUS() * FINE_STRUCTURE_CONSTANT();
    }

    // BEGIN SEMIRELATIVISTIC BETHE HEITLER CROSSECTION CALCULATION

    static double tungstenFiltration(double tungstenAtt, double x, double sintakeoffAngle)
    {
        // filtrates a photon in the tungsten anode. hv :photon energy [keV], x depth in tungsten (density*distance) [mg/cm2], takeoffAngle: scatter angle in radians
        return std::exp(-tungstenAtt * x * 0.001 / sintakeoffAngle); // 0.001 is for mg/cm2 -> g/cm2
    }

    static double betheHeitlerCrossSection(double hv, double Ti)
    {
        constexpr auto phi_bar_const = (PHI_BAR() * 2) / 3;
        const auto Ei = ELECTRON_REST_MASS() + Ti;
        const auto Ef = Ei - hv;
        const auto pic_sqr = Ei * Ei - ELECTRON_REST_MASS() * ELECTRON_REST_MASS();
        const auto pic = std::sqrt(pic_sqr);
        const auto pfc_sqr = Ef * Ef - ELECTRON_REST_MASS() * ELECTRON_REST_MASS();
        if (pfc_sqr <= 0)
            return 0;
        const auto pfc = std::sqrt(pfc_sqr);

        const auto L = 2 * std::log((Ei * Ef + pic * pfc - ELECTRON_REST_MASS() * ELECTRON_REST_MASS()) / (ELECTRON_REST_MASS() * hv));

        const auto columb_correction = pic / pfc;
        return phi_bar_const * (4 * Ei * Ef * L - 7 * pic * pfc) / (hv * pic * pic) * columb_correction;
    }

    static double betheHeitlerSpectra(double T0, double hv, double takeoffAngle)
    {

        const auto& tungsten = AtomHandler::Atom(TUNGSTEN_ATOMIC_NUMBER);
        const auto photo = interpolate(tungsten.photoel, hv);
        const auto coherent = interpolate(tungsten.coherent, hv);
        const auto incoherent = interpolate(tungsten.incoherent, hv);

        const auto tungstenTotAtt = photo + coherent + incoherent;
        constexpr auto xmax = 14.0;
        constexpr auto umax = 1.0;
        constexpr auto xstep = 0.1;
        constexpr auto ustep = 0.005;

        const auto sinTakeoffAngle = std::sin(takeoffAngle);

        double x = 0;
        double I_obs = 0;
        if (hv <= 0)
            return 0;
        while (x <= xmax) {
            double I_step = 0;
            double u = ustep;
            while (u <= umax) {
                I_step = I_step + betheHeitlerCrossSection(hv, T0 * u) * electronDensity(u, x, T0) * ustep;
                u = u + ustep;
            }
            I_obs = I_obs + I_step * tungstenFiltration(tungstenTotAtt, x, sinTakeoffAngle) * xstep;
            x = x + xstep;
        }
        return I_obs;
    }

    static std::array<std::pair<double, double>, 5> characteristicTungstenKedge(double T0, double takeoffAngle)
    {
        constexpr std::array<double, 4> k_edge_energies { 59.3, 58.0, 67.2, 69.1 };
        constexpr std::array<double, 4> k_edge_fractions { 0.505, 0.291, 0.162, 0.042 };
        constexpr auto P = 0.33;
        constexpr auto omega_k = 0.94;
        constexpr auto rk = 4.4;

        std::array<std::pair<double, double>, 5> char_rad;
        std::transform(std::execution::par_unseq, k_edge_energies.cbegin(), k_edge_energies.cend(), k_edge_fractions.cbegin(), char_rad.begin(),
            [=](const auto energy, const auto fraction) -> std::pair<double, double> { return std::make_pair(energy, (1 + P) * fraction * rk * omega_k * betheHeitlerSpectra(T0, energy, takeoffAngle)); });
        return char_rad;
    }

    // END SEMIRELATIVISTIC BETHE HEITLER CROSS SECTION CALCULATION

}
}

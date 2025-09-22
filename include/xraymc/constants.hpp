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

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "xraymc/floating.hpp"
#include <charconv>
#include <numbers>
#include <optional>
#include <string_view>

namespace xraymc {

template <Floating T = double>
constexpr T GEOMETRIC_ERROR()
{
    return T { 1E-6 };
}

template <Floating T = double>
constexpr T MAX_ENERGY()
{
    return T { XRAYMCLIB_MAXENERGY };
}

template <Floating T = double>
constexpr T MIN_ENERGY()
{
    return T { XRAYMCLIB_MINENERGY };
}

template <Floating T = double>
consteval T KEV_TO_ANGSTROM()
{
    /* consteval T hc_si = T { 1.239841193E-6 }; // ev*m
    consteval T m2A = T { 1E10 }; // meters to Angstrom
    consteval T eV2keV = T { 1E-3 }; // eV to keV
    consteval T hc = hc_si * m2A * eV2keV; // kev*Angstrom
    consteval T hc_inv = T { 1.0 } / hc;
    */
    return T { 12.398520 };
}

template <Floating T = double>
consteval T PI_VAL()
{
    return std::numbers::pi_v<T>;
}

template <Floating T = double>
consteval T DEG_TO_RAD()
{
    return PI_VAL<T>() / T { 180 };
}

template <Floating T = double>
consteval T RAD_TO_DEG()
{
    return T { 180 } / PI_VAL<T>();
}

template <Floating T = double>
consteval T KEV_TO_MJ()
{
    return T { 1.6021773e-13 }; // milli Joules}
}

template <Floating T = double>
consteval T MJ_TO_KEV()
{
    return T { 1 } / KEV_TO_MJ<T>();
}

template <Floating T = double>
consteval T ELECTRON_REST_MASS()
{
    return T { 510.9989461 }; // kev/c^2
}

}
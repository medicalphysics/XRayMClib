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

/**
 * @brief Small offset added to ray origins to avoid self-intersection artifacts.
 *
 * When a ray is advanced to a geometry boundary, its origin is nudged by this
 * amount (1×10⁻⁶ cm) to ensure it clears the surface and does not re-intersect
 * the same face in the next transport step.
 *
 * @tparam T Floating-point type (default: double).
 * @return The geometric error tolerance in cm.
 */
template <Floating T = double>
constexpr T GEOMETRIC_ERROR()
{
    return T { 1E-6 };
}

/**
 * @brief Maximum photon energy supported by the cross-section tables (in keV).
 *
 * Set at compile time via the `XRAYMCLIB_MAXENERGY` preprocessor macro.
 * Particles with energies above this value cannot be transported.
 *
 * @tparam T Floating-point type (default: double).
 * @return The maximum supported energy in keV.
 */
template <Floating T = double>
constexpr T MAX_ENERGY()
{
    return T { XRAYMCLIB_MAXENERGY };
}

/**
 * @brief Minimum photon energy supported by the cross-section tables (in keV).
 *
 * Set at compile time via the `XRAYMCLIB_MINENERGY` preprocessor macro.
 * Particles with energies below this value are considered absorbed.
 *
 * @tparam T Floating-point type (default: double).
 * @return The minimum supported energy in keV.
 */
template <Floating T = double>
constexpr T MIN_ENERGY()
{
    return T { XRAYMCLIB_MINENERGY };
}

/**
 * @brief Conversion factor from photon energy (keV) to wavelength (Å).
 *
 * Derived from the relation λ = hc / E, where hc = 12.398520 keV·Å.
 * Multiply an energy in keV by this value to obtain the corresponding
 * wavelength in Ångström.
 *
 * @tparam T Floating-point type (default: double).
 * @return hc in units of keV·Å (12.398520).
 */
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

/**
 * @brief The mathematical constant π.
 *
 * @tparam T Floating-point type (default: double).
 * @return π to the precision of `T`, via `std::numbers::pi_v<T>`.
 */
template <Floating T = double>
consteval T PI_VAL()
{
    return std::numbers::pi_v<T>;
}

/**
 * @brief Conversion factor from degrees to radians (π / 180).
 *
 * @tparam T Floating-point type (default: double).
 * @return The factor to multiply degrees by to obtain radians.
 */
template <Floating T = double>
consteval T DEG_TO_RAD()
{
    return PI_VAL<T>() / T { 180 };
}

/**
 * @brief Conversion factor from radians to degrees (180 / π).
 *
 * @tparam T Floating-point type (default: double).
 * @return The factor to multiply radians by to obtain degrees.
 */
template <Floating T = double>
consteval T RAD_TO_DEG()
{
    return T { 180 } / PI_VAL<T>();
}

/**
 * @brief Conversion factor from keV to millijoules (mJ).
 *
 * 1 keV = 1.6021773×10⁻¹³ mJ. Used when converting deposited energy to
 * SI-compatible dose units.
 *
 * @tparam T Floating-point type (default: double).
 * @return The factor to multiply keV by to obtain mJ.
 */
template <Floating T = double>
consteval T KEV_TO_MJ()
{
    return T { 1.6021773e-13 }; // milli Joules}
}

/**
 * @brief Conversion factor from millijoules (mJ) to keV.
 *
 * Reciprocal of `KEV_TO_MJ<T>()`.
 *
 * @tparam T Floating-point type (default: double).
 * @return The factor to multiply mJ by to obtain keV.
 */
template <Floating T = double>
consteval T MJ_TO_KEV()
{
    return T { 1 } / KEV_TO_MJ<T>();
}

/**
 * @brief Rest mass energy of the electron in keV.
 *
 * mₑc² = 510.9989461 keV. Used in Compton and pair-production kinematics.
 *
 * @tparam T Floating-point type (default: double).
 * @return The electron rest mass energy in keV/c².
 */
template <Floating T = double>
consteval T ELECTRON_REST_MASS()
{
    return T { 510.9989461 }; // kev/c^2
}

}
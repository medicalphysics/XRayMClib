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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace xraymc {

/**
 * @brief Concept that constrains a type to standard floating-point types.
 *
 * Satisfied by `float`, `double`, and `long double`. Used throughout XRayMClib
 * as a template constraint to ensure numeric precision compatibility for
 * physical quantities such as energy, position, and cross-section values.
 *
 * @tparam T The type to check.
 */
template <typename T>
concept Floating = std::is_floating_point<T>::value;

/**
 * @brief Returns true if @p value is neither NaN nor ±infinity.
 *
 * Unlike `std::isfinite`, this inspects the IEEE-754 exponent bits directly
 * through integer operations, so the test cannot be optimized away by flags
 * such as `-ffast-math`, `-ffinite-math-only` or `/fp:fast` that let the
 * compiler assume NaN and infinity never occur.
 *
 * @tparam T `float` or `double`.
 * @param value Value to test.
 * @return True if the exponent field is not all ones.
 */
template <typename T>
    requires(std::is_same_v<T, float> || std::is_same_v<T, double>)
constexpr bool isFinite(const T value) noexcept
{
    static_assert(std::numeric_limits<T>::is_iec559, "isFinite requires IEEE-754 floating point");
    using U = std::conditional_t<std::is_same_v<T, float>, std::uint32_t, std::uint64_t>;
    constexpr int mantissaBits = std::numeric_limits<T>::digits - 1;
    constexpr int exponentBits = static_cast<int>(sizeof(T) * 8) - 1 - mantissaBits;
    constexpr U exponentMask = ((U { 1 } << exponentBits) - 1) << mantissaBits;
    return (std::bit_cast<U>(value) & exponentMask) != exponentMask;
}

}
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

#include "xraymc/floating.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

// Header library for simple 3D vector math
namespace xraymc {
namespace vectormath {

    template <typename T>
    concept Index = std::is_integral_v<T> && std::is_same<bool, T>::value == false;

    template <Floating T>
    [[nodiscard]] inline constexpr T length_sqr(const std::array<T, 3>& vec) noexcept
    {
        return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    }

    template <Floating T>
    [[nodiscard]] inline T length(const std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return std::sqrt(lsqr);
    }

    template <Floating T>
    [[nodiscard]] inline std::pair<std::array<T, 3>, std::array<T, 3>> splice(const std::array<T, 6>& a)
    {
        std::array a1 { a[0], a[1], a[2] };
        std::array a2 { a[3], a[4], a[5] };
        return std::make_pair(a1, a2);
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 6> join(const std::array<T, 3>& a, const std::array<T, 3>& b)
    {
        return std::array { a[0], a[1], a[2], b[0], b[1], b[2] };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> add(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
    {
        return std::array<T, 3> { v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
    }

    template <Floating T, typename... Args>
    [[nodiscard]] inline constexpr std::array<T, 3> add(const std::array<T, 3>& v1, const std::array<T, 3>& v2, Args... args)
    {
        auto v = add(v1, v2);
        return add(v, args...);
    }

    template <Floating T>
    [[nodiscard]] inline constexpr auto add(const std::array<T, 3>& v1, T v2) noexcept
    {
        return std::array { v1[0] + v2, v1[1] + v2, v1[2] + v2 };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr auto add(T v2, const std::array<T, 3>& v1) noexcept
    {
        return add(v1, v2);
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> subtract(const std::array<T, 3>& v1, const std::array<T, 3>& v2) noexcept
    {
        return std::array<T, 3> { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(const std::array<T, 3>& v, T s)
    {
        return std::array { v[0] * s, v[1] * s, v[2] * s };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(T s, const std::array<T, 3>& v)
    {
        return std::array { v[0] * s, v[1] * s, v[2] * s };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(const std::array<T, 3>& v, const std::array<T, 3>& s)
    {
        return std::array { v[0] * s[0], v[1] * s[1], v[2] * s[2] };
    }

    template <Floating T>
    inline void normalize(std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        constexpr T one { 1 };
        const T norm = one / std::sqrt(lsqr);
        vec[0] *= norm;
        vec[1] *= norm;
        vec[2] *= norm;
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> normalized(const std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        const T norm = T { 1 } / std::sqrt(lsqr);
        return scale(vec, norm);
    }

    template <Floating T>
    [[nodiscard]] inline constexpr T dot(const std::array<T, 3>& v1, const std::array<T, 3>& v2) noexcept
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<T, 3>& v1, const std::array<T, 3> v2) noexcept
    {
        return std::array<T, 3> {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<T, 6>& v1) noexcept
    {
        return std::array<T, 3> {
            v1[1] * v1[5] - v1[2] * v1[4],
            v1[2] * v1[3] - v1[0] * v1[5],
            v1[0] * v1[4] - v1[1] * v1[3]
        };
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<std::array<T, 3>, 2>& v) noexcept
    {
        return cross(v[0], v[1]);
    }

    template <Floating T>
    [[nodiscard]] inline constexpr T tripleProduct(const std::array<T, 3>& v1, const std::array<T, 3>& v2, const std::array<T, 3>& v3) noexcept
    {
        return dot(v1, cross(v2, v3));
    }

    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> rotate(const std::array<T, 3>& vec, const std::array<T, 3>& axis, const T sinAngle, const T cosAngle) noexcept
    {
        const auto ax = cross(axis, vec);
        const auto va = dot(vec, axis);
        const auto v2 = cross(ax, axis);
        return std::array<T, 3> {
            va * axis[0] + cosAngle * v2[0] + sinAngle * ax[0],
            va * axis[1] + cosAngle * v2[1] + sinAngle * ax[1],
            va * axis[2] + cosAngle * v2[2] + sinAngle * ax[2]
        };
    }

    template <Floating T>
    [[nodiscard]] inline std::array<T, 3> rotate(const std::array<T, 3>& vec, const std::array<T, 3>& axis, const T angle) noexcept
    {
        const T sang = std::sin(angle);
        const T cang = std::cos(angle);
        return rotate(vec, axis, sang, cang);
    }

    template <Floating T>
    [[nodiscard]] inline T angleBetween(const std::array<T, 3>& vec1, const std::array<T, 3>& vec2) noexcept
    {
        // Herons formula for numeric stable angle computation
        // Do not edit parenthesis and such

        const auto vec3 = subtract(vec1, vec2);
        const T a = length(vec1);
        const T b = length(vec2);
        const T c = length(vec3);

        const T u = b >= c ? c - (a - b) : b - (a - c);

        const T nom = ((a - b) + c) * u;
        const T den = (a + (b + c)) * ((a - c) + b);
        return T { 2 } * std::atan(std::sqrt(nom / den));
    }

    template <Index U = std::size_t, Floating T>
    [[nodiscard]] inline U argmin3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x <= y ? x <= z ? 0 : 2 : y <= z ? 1
                                                : 2;
    }

    template <Index U = std::size_t, Floating T>
    [[nodiscard]] inline U argmax3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x >= y ? x >= z ? 0 : 2 : y >= z ? 1
                                                : 2;
    }
    template <Index U = std::size_t, Index T>
    [[nodiscard]] inline U argmax3(const std::array<T, 3>& vec) noexcept
    {
        const T x = vec[0];
        const T y = vec[1];
        const T z = vec[2];
        return x >= y ? x >= z ? 0 : 2 : y >= z ? 1
                                                : 2;
    }

    template <Floating T>
    [[nodiscard]] inline std::array<T, 3> changeBasis(const std::array<T, 3>& b1, const std::array<T, 3>& b2, const std::array<T, 3>& b3, const std::array<T, 3>& vector) noexcept
    {
        std::array<T, 3> res {
            b1[0] * vector[0] + b2[0] * vector[1] + b3[0] * vector[2],
            b1[1] * vector[0] + b2[1] * vector[1] + b3[1] * vector[2],
            b1[2] * vector[0] + b2[2] * vector[1] + b3[2] * vector[2]
        };
        return res;
    }

    template <Floating T>
    [[nodiscard]] inline std::array<T, 3> changeBasisInverse(const std::array<T, 3>& b1, const std::array<T, 3>& b2, const std::array<T, 3>& b3, const std::array<T, 3>& vector) noexcept
    {
        std::array<T, 3> res {
            b1[0] * vector[0] + b1[1] * vector[1] + b1[2] * vector[2],
            b2[0] * vector[0] + b2[1] * vector[1] + b2[2] * vector[2],
            b3[0] * vector[0] + b3[1] * vector[1] + b3[2] * vector[2]
        };
        return res;
    }

    template <Floating T>
    [[nodiscard]] inline std::array<T, 3> peturb(const std::array<T, 3>& vec, const T theta, const T phi) noexcept
    {
        // rotates a unit vector theta degrees from its current direction
        // phi degrees about a arbitrary axis orthogonal to the direction vector

        const auto minInd = argmin3<std::uint_fast32_t, T>(vec);
        std::array<T, 3> k { 0, 0, 0 };
        k[minInd] = T { 1 };

        auto vec_xy_raw = vectormath::cross(vec, k);
        normalize(vec_xy_raw);

        // rotating the arbitrary orthogonal axis about vector direction
        const auto vec_xy = rotate(vec_xy_raw, vec, phi);

        auto res = rotate(vec, vec_xy, theta);
        // We normalize result in case of multiple calls on same vector
        normalize(res);
        return res;
    }

}
}
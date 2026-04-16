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

/**
 * @file vectormath.hpp
 * @brief Header-only 3D vector math utilities used throughout XRayMClib.
 *
 * All functions operate on `std::array<T, 3>` (or `std::array<T, 6>` for packed pairs)
 * and are templated on a `Floating` type (typically `double`). The design keeps
 * vector values on the stack and avoids heap allocation.
 *
 * ## Functional groups
 *
 * | Group              | Functions                                                          |
 * |--------------------|--------------------------------------------------------------------|
 * | Length / norm      | `length_sqr`, `length`, `normalize`, `normalized`                 |
 * | Arithmetic         | `add`, `subtract`, `scale`                                        |
 * | Products           | `dot`, `cross`, `tripleProduct`                                   |
 * | Rotation           | `rotate` (sin/cos or angle overload), `peturb`                    |
 * | Angle              | `angleBetween`                                                     |
 * | Arg-extrema        | `argmin3`, `argmax3`                                              |
 * | Basis change       | `changeBasis`, `changeBasisInverse`                               |
 * | Pack / unpack      | `splice`, `join`                                                  |
 *
 * ## Design notes
 * - All functions are `[[nodiscard]]` where they return a value.
 * - `normalize` modifies its argument in place; `normalized` returns a copy.
 * - `angleBetween` uses Heron's numerically stable formula and should be preferred
 *   over `std::acos(dot(a, b))` for nearly parallel or anti-parallel vectors.
 * - `peturb` is the primary entry point for scattering deflections: it maps a
 *   polar angle @p theta and azimuthal angle @p phi onto a new direction without
 *   requiring explicit coordinate-frame bookkeeping.
 */
#pragma once

#include "xraymc/floating.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace xraymc {

/**
 * @namespace xraymc::vectormath
 * @brief Header-only 3D vector math utilities used by the transport core.
 *
 * Vectors are represented as `std::array<T, 3>` where `T` satisfies the
 * `Floating` concept (see floating.hpp). Packed AABB-style values use
 * `std::array<T, 6>` with the layout `[xmin, ymin, zmin, xmax, ymax, zmax]`.
 *
 * The namespace intentionally avoids a dedicated vector class so that all
 * geometry data can be stored in plain aggregate types, simplifying
 * serialization and SIMD-friendly memory layouts.
 */
namespace vectormath {

    /**
     * @brief Concept satisfied by non-boolean integral types used as array indices.
     * @tparam T The candidate type.
     */
    template <typename T>
    concept Index = std::is_integral_v<T> && std::is_same<bool, T>::value == false;

    /**
     * @brief Returns the squared Euclidean length of a 3D vector.
     * @tparam T Floating-point type.
     * @param vec Input vector.
     * @return vec·vec (no square root taken).
     */
    template <Floating T>
    [[nodiscard]] inline constexpr T length_sqr(const std::array<T, 3>& vec) noexcept
    {
        return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    }

    /**
     * @brief Returns the Euclidean length of a 3D vector.
     * @tparam T Floating-point type.
     * @param vec Input vector.
     * @return √(vec·vec).
     */
    template <Floating T>
    [[nodiscard]] inline T length(const std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return std::sqrt(lsqr);
    }

    /**
     * @brief Splits a 6-element array into two 3-element arrays.
     * @tparam T Floating-point type.
     * @param a Input array [x0, y0, z0, x1, y1, z1].
     * @return Pair of {first three elements, last three elements}.
     */
    template <Floating T>
    [[nodiscard]] inline std::pair<std::array<T, 3>, std::array<T, 3>> splice(const std::array<T, 6>& a)
    {
        std::array a1 { a[0], a[1], a[2] };
        std::array a2 { a[3], a[4], a[5] };
        return std::make_pair(a1, a2);
    }

    /**
     * @brief Joins two 3-element arrays into a single 6-element array.
     * @tparam T Floating-point type.
     * @param a First vector [x0, y0, z0].
     * @param b Second vector [x1, y1, z1].
     * @return Array [x0, y0, z0, x1, y1, z1].
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 6> join(const std::array<T, 3>& a, const std::array<T, 3>& b)
    {
        return std::array { a[0], a[1], a[2], b[0], b[1], b[2] };
    }

    /**
     * @brief Returns the element-wise sum of two 3D vectors.
     * @tparam T Floating-point type.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @return v1 + v2.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> add(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
    {
        return std::array<T, 3> { v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
    }

    /**
     * @brief Returns the element-wise sum of three or more 3D vectors (variadic).
     * @tparam T    Floating-point type.
     * @tparam Args Additional `std::array<T, 3>` arguments.
     * @param v1   First vector.
     * @param v2   Second vector.
     * @param args Further vectors to accumulate.
     * @return Element-wise sum of all arguments.
     */
    template <Floating T, typename... Args>
    [[nodiscard]] inline constexpr std::array<T, 3> add(const std::array<T, 3>& v1, const std::array<T, 3>& v2, Args... args)
    {
        auto v = add(v1, v2);
        return add(v, args...);
    }

    /**
     * @brief Adds a scalar to every element of a 3D vector.
     * @tparam T Floating-point type.
     * @param v1 Input vector.
     * @param v2 Scalar to add to each component.
     * @return {v1[0]+v2, v1[1]+v2, v1[2]+v2}.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr auto add(const std::array<T, 3>& v1, T v2) noexcept
    {
        return std::array { v1[0] + v2, v1[1] + v2, v1[2] + v2 };
    }

    /**
     * @brief Adds a scalar to every element of a 3D vector (scalar-first overload).
     * @tparam T Floating-point type.
     * @param v2 Scalar to add to each component.
     * @param v1 Input vector.
     * @return {v1[0]+v2, v1[1]+v2, v1[2]+v2}.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr auto add(T v2, const std::array<T, 3>& v1) noexcept
    {
        return add(v1, v2);
    }

    /**
     * @brief Returns the element-wise difference of two 3D vectors.
     * @tparam T Floating-point type.
     * @param v1 Minuend vector.
     * @param v2 Subtrahend vector.
     * @return v1 − v2.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> subtract(const std::array<T, 3>& v1, const std::array<T, 3>& v2) noexcept
    {
        return std::array<T, 3> { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
    }

    /**
     * @brief Scales a 3D vector by a scalar (vector-first).
     * @tparam T Floating-point type.
     * @param v Input vector.
     * @param s Scale factor.
     * @return v * s.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(const std::array<T, 3>& v, T s)
    {
        return std::array { v[0] * s, v[1] * s, v[2] * s };
    }

    /**
     * @brief Scales a 3D vector by a scalar (scalar-first).
     * @tparam T Floating-point type.
     * @param s Scale factor.
     * @param v Input vector.
     * @return v * s.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(T s, const std::array<T, 3>& v)
    {
        return std::array { v[0] * s, v[1] * s, v[2] * s };
    }

    /**
     * @brief Element-wise multiplies two 3D vectors.
     * @tparam T Floating-point type.
     * @param v Input vector.
     * @param s Per-component scale factors.
     * @return {v[0]*s[0], v[1]*s[1], v[2]*s[2]}.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> scale(const std::array<T, 3>& v, const std::array<T, 3>& s)
    {
        return std::array { v[0] * s[0], v[1] * s[1], v[2] * s[2] };
    }

    /**
     * @brief Normalizes a 3D vector to unit length in place.
     *
     * Multiplies each component by 1/‖vec‖. Behaviour is undefined if @p vec
     * is the zero vector.
     *
     * @tparam T Floating-point type.
     * @param vec Vector to normalize; modified in place.
     */
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

    /**
     * @brief Returns a normalized copy of a 3D vector.
     *
     * Does not modify the input. Behaviour is undefined if @p vec is the zero vector.
     *
     * @tparam T Floating-point type.
     * @param vec Input vector.
     * @return Unit vector in the direction of @p vec.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> normalized(const std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        const T norm = T { 1 } / std::sqrt(lsqr);
        return scale(vec, norm);
    }

    /**
     * @brief Returns the dot (inner) product of two 3D vectors.
     * @tparam T Floating-point type.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @return v1 · v2.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr T dot(const std::array<T, 3>& v1, const std::array<T, 3>& v2) noexcept
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    /**
     * @brief Returns the cross product of two 3D vectors.
     * @tparam T Floating-point type.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @return v1 × v2.
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<T, 3>& v1, const std::array<T, 3> v2) noexcept
    {
        return std::array<T, 3> {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
    }

    /**
     * @brief Returns the cross product of two vectors packed into a 6-element array.
     *
     * Interprets elements [0,1,2] as the first vector and [3,4,5] as the second.
     *
     * @tparam T Floating-point type.
     * @param v1 Packed pair of 3D vectors [x0,y0,z0, x1,y1,z1].
     * @return v1[0:3] × v1[3:6].
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<T, 6>& v1) noexcept
    {
        return std::array<T, 3> {
            v1[1] * v1[5] - v1[2] * v1[4],
            v1[2] * v1[3] - v1[0] * v1[5],
            v1[0] * v1[4] - v1[1] * v1[3]
        };
    }

    /**
     * @brief Returns the cross product of two vectors stored in a 2-element array of arrays.
     * @tparam T Floating-point type.
     * @param v Array of two 3D vectors {v[0], v[1]}.
     * @return v[0] × v[1].
     */
    template <Floating T>
    [[nodiscard]] inline constexpr std::array<T, 3> cross(const std::array<std::array<T, 3>, 2>& v) noexcept
    {
        return cross(v[0], v[1]);
    }

    /**
     * @brief Returns the scalar triple product v1 · (v2 × v3).
     *
     * Geometrically equals the signed volume of the parallelepiped spanned by the
     * three vectors.
     *
     * @tparam T Floating-point type.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param v3 Third vector.
     * @return v1 · (v2 × v3).
     */
    template <Floating T>
    [[nodiscard]] inline constexpr T tripleProduct(const std::array<T, 3>& v1, const std::array<T, 3>& v2, const std::array<T, 3>& v3) noexcept
    {
        return dot(v1, cross(v2, v3));
    }

    /**
     * @brief Rotates a 3D vector about an axis by an angle given as sin/cos components.
     *
     * Uses Rodrigues' rotation formula decomposed as:
     *   result = (vec · axis) · axis + cosAngle · (axis × vec) × axis + sinAngle · (axis × vec)
     *
     * @tparam T Floating-point type.
     * @param vec      Vector to rotate.
     * @param axis     Unit rotation axis.
     * @param sinAngle Sine of the rotation angle.
     * @param cosAngle Cosine of the rotation angle.
     * @return Rotated vector.
     */
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

    /**
     * @brief Rotates a 3D vector about an axis by a given angle in radians.
     * @tparam T Floating-point type.
     * @param vec   Vector to rotate.
     * @param axis  Unit rotation axis.
     * @param angle Rotation angle in radians.
     * @return Rotated vector.
     */
    template <Floating T>
    [[nodiscard]] inline std::array<T, 3> rotate(const std::array<T, 3>& vec, const std::array<T, 3>& axis, const T angle) noexcept
    {
        const T sang = std::sin(angle);
        const T cang = std::cos(angle);
        return rotate(vec, axis, sang, cang);
    }

    /**
     * @brief Computes the angle between two 3D vectors in radians.
     *
     * Uses the numerically stable Heron's formula variant to avoid cancellation
     * errors for nearly parallel or anti-parallel vectors.
     *
     * @tparam T Floating-point type.
     * @param vec1 First vector (need not be unit length).
     * @param vec2 Second vector (need not be unit length).
     * @return Angle in radians in [0, π].
     */
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

    /**
     * @brief Returns the index of the component with the smallest absolute value.
     *
     * Useful for constructing an arbitrary orthogonal vector without cancellation.
     *
     * @tparam U Return type, must satisfy `Index` (default `std::size_t`).
     * @tparam T Floating-point element type.
     * @param vec Input vector.
     * @return Index in {0, 1, 2} of the component with the smallest |value|.
     */
    template <Index U = std::size_t, Floating T>
    [[nodiscard]] inline U argmin3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x <= y ? x <= z ? 0 : 2 : y <= z ? 1
                                                : 2;
    }

    /**
     * @brief Returns the index of the component with the largest absolute value (floating-point).
     * @tparam U Return type, must satisfy `Index` (default `std::size_t`).
     * @tparam T Floating-point element type.
     * @param vec Input vector.
     * @return Index in {0, 1, 2} of the component with the largest |value|.
     */
    template <Index U = std::size_t, Floating T>
    [[nodiscard]] inline U argmax3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x >= y ? x >= z ? 0 : 2 : y >= z ? 1
                                                : 2;
    }

    /**
     * @brief Returns the index of the component with the largest value (integer overload).
     * @tparam U Return type, must satisfy `Index` (default `std::size_t`).
     * @tparam T Integer element type satisfying `Index`.
     * @param vec Input vector.
     * @return Index in {0, 1, 2} of the component with the largest value.
     */
    template <Index U = std::size_t, Index T>
    [[nodiscard]] inline U argmax3(const std::array<T, 3>& vec) noexcept
    {
        const T x = vec[0];
        const T y = vec[1];
        const T z = vec[2];
        return x >= y ? x >= z ? 0 : 2 : y >= z ? 1
                                                : 2;
    }

    /**
     * @brief Transforms a vector from local to world coordinates.
     *
     * Applies the 3×3 matrix whose columns are the basis vectors @p b1, @p b2, @p b3:
     *   result = b1·v[0] + b2·v[1] + b3·v[2]  (column-major multiply)
     *
     * @tparam T Floating-point type.
     * @param b1     First column basis vector (local x-axis in world coords).
     * @param b2     Second column basis vector (local y-axis in world coords).
     * @param b3     Third column basis vector (local z-axis in world coords).
     * @param vector Vector expressed in the local basis.
     * @return Vector expressed in world coordinates.
     */
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

    /**
     * @brief Transforms a vector from world to local coordinates (inverse basis transform).
     *
     * Applies the transpose of the 3×3 matrix whose columns are @p b1, @p b2, @p b3
     * (which equals the inverse when the basis is orthonormal):
     *   result[i] = bi · vector
     *
     * @tparam T Floating-point type.
     * @param b1     First basis vector (local x-axis in world coords).
     * @param b2     Second basis vector (local y-axis in world coords).
     * @param b3     Third basis vector (local z-axis in world coords).
     * @param vector Vector expressed in world coordinates.
     * @return Vector expressed in the local basis.
     */
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

    /**
     * @brief Deflects a unit direction vector by polar angle @p theta and azimuthal angle @p phi.
     *
     * Constructs an arbitrary orthogonal axis by crossing @p vec with the world axis
     * whose component in @p vec is smallest (to avoid numerical cancellation), then:
     * 1. Rotates that axis by @p phi around @p vec to obtain the scattering plane normal.
     * 2. Rotates @p vec by @p theta around the chosen axis.
     * The result is re-normalized to correct accumulated floating-point drift.
     *
     * Used by the interaction samplers to apply a scattering deflection in the particle's
     * local reference frame without requiring explicit Euler-angle bookkeeping.
     *
     * @tparam T Floating-point type.
     * @param vec   Incident unit direction vector.
     * @param theta Polar deflection angle from the incident direction (radians).
     * @param phi   Azimuthal rotation angle around the incident direction (radians).
     * @return New unit direction vector after the deflection.
     */
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
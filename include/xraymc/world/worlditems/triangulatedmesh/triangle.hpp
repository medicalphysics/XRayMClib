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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <optional>

namespace xraymc {

/**
 * @brief A single triangle in 3-D space, used as the primitive element of a triangulated mesh.
 *
 * Stores three vertex positions and provides geometric queries (normal, centroid, AABB,
 * area) as well as affine transformations (translate, scale, rotate, mirror) and a
 * Möller–Trumbore ray–triangle intersection test.
 */
class Triangle {
public:
    /**
     * @brief Constructs a triangle from three individual vertex positions.
     * @param first  First vertex in cm.
     * @param second Second vertex in cm.
     * @param third  Third vertex in cm.
     */
    Triangle(const std::array<double, 3>& first, const std::array<double, 3>& second, const std::array<double, 3>& third)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
    }

    /**
     * @brief Constructs a triangle from a packed array of three vertex positions.
     * @param vertices Array of three vertices, each as {x, y, z} in cm.
     */
    Triangle(const std::array<std::array<double, 3>, 3>& vertices)
        : m_vertices(vertices)
    {
    }

    /**
     * @brief Constructs a triangle by reading nine consecutive doubles from a raw pointer.
     *
     * Vertices are read in row-major order: v0={[0],[1],[2]}, v1={[3],[4],[5]}, v2={[6],[7],[8]}.
     * @param first_element Pointer to the first of nine contiguous double values.
     */
    Triangle(const double* first_element)
    {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto flatIdx = i * 3 + j;
                m_vertices[i][j] = *(first_element + flatIdx);
            }
    }

    /// @brief Defaulted three-way comparison (lexicographic over vertex coordinates).
    auto operator<=>(const Triangle& other) const = default;

    /**
     * @brief Translates all three vertices by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }

    /**
     * @brief Mirrors all three vertices through a world-space point.
     *
     * Each vertex v is mapped to v + (-2)(v - point) = 2·point - v.
     * @param point The point of reflection in cm.
     */
    void mirror(const std::array<double, 3>& point)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            const auto d = vectormath::subtract(m_vertices[i], point);
            m_vertices[i] = vectormath::add(m_vertices[i], vectormath::scale(d, -2.0));
        }
    }

    /**
     * @brief Mirrors all three vertices about a plane perpendicular to axis @p dim at @p value.
     * @param value Coordinate of the mirror plane along @p dim.
     * @param dim   Axis index: 0 = x, 1 = y, 2 = z.
     */
    void mirror(const double value, const std::uint_fast32_t dim)
    {
        for (std::uint_fast32_t i = 0; i < 3; ++i) {
            const auto d = m_vertices[i][dim] - value;
            m_vertices[i][dim] += d * -2.0;
        }
    }

    /**
     * @brief Uniformly scales all vertex coordinates by @p scale.
     * @param scale Scale factor applied to every coordinate of every vertex.
     */
    void scale(double scale)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] *= scale;
            }
        });
    }

    /**
     * @brief Rotates all vertices by @p radians around @p axis (through the origin).
     * @param radians Rotation angle in radians.
     * @param axis    Rotation axis (need not be normalized).
     */
    void rotate(double radians, const std::array<double, 3>& axis)
    {
        std::transform(std::execution::unseq, m_vertices.cbegin(), m_vertices.cend(), m_vertices.begin(), [&](auto& vert) {
            return vectormath::rotate(vert, axis, radians);
        });
    }

    /**
     * @brief Returns the normalized surface normal computed from the cross product of two edges.
     * @return Unit normal vector (v1-v0) × (v2-v0), normalized.
     */
    std::array<double, 3> planeVector() const noexcept
    {
        const auto a = vectormath::subtract(m_vertices[1], m_vertices[0]);
        const auto b = vectormath::subtract(m_vertices[2], m_vertices[0]);
        const auto n = vectormath::cross(a, b);
        return vectormath::normalized(n);
    }

    /// @brief Returns the three vertex positions as a packed 3×3 array.
    const std::array<std::array<double, 3>, 3>& vertices() const
    {
        return m_vertices;
    }

    /// @brief Returns the centroid of the triangle (mean of the three vertices) in cm.
    std::array<double, 3> center() const
    {
        std::array<double, 3> cent { 0, 0, 0 };
        for (const auto& vert : m_vertices) {
            for (std::size_t i = 0; i < 3; i++)
                cent[i] += vert[i];
        }
        constexpr double factor { 1 / 3.0 };
        for (std::size_t i = 0; i < 3; i++)
            cent[i] *= factor;
        return cent;
    }

    /// @brief Returns the axis-aligned bounding box of the triangle as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
    std::array<double, 6> AABB() const
    {
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (std::size_t j = 0; j < 3; j++) {
            for (std::size_t i = 0; i < 3; i++) {
                aabb[i] = std::min(aabb[i], m_vertices[j][i]);
            }
            for (std::size_t i = 0; i < 3; i++) {
                const auto idx = i + 3;
                aabb[idx] = std::max(aabb[idx], m_vertices[j][i]);
            }
        }
        return aabb;
    }

    /**
     * @brief Returns the surface area of the triangle in cm².
     *
     * Computed as |AB × AC| / 2.
     * @return Triangle area in cm².
     */
    double area() const
    { // Triangle of ABC, area = |AB x AC|/2
        const auto AB = vectormath::subtract(m_vertices[1], m_vertices[0]);
        const auto AC = vectormath::subtract(m_vertices[2], m_vertices[0]);
        return vectormath::length(vectormath::cross(AB, AC)) / 2;
    }

    /**
     * @brief Tests a particle ray against the triangle using the Möller–Trumbore algorithm.
     * @param p Particle whose position and direction define the ray.
     * @return The ray parameter t ≥ 0 at the intersection point, or std::nullopt if the
     *         ray misses the triangle or hits from behind.
     */
    std::optional<double> intersect(const ParticleType auto& p) const
    {
        // from moller trombore paper
        const auto& v0 = m_vertices[0];
        const auto& v1 = m_vertices[1];
        const auto& v2 = m_vertices[2];

        const auto E1 = vectormath::subtract(v1, v0);
        const auto TT = vectormath::subtract(p.pos, v0);
        const auto Q = vectormath::cross(TT, E1);

        const auto E2 = vectormath::subtract(v2, v0);
        const auto P = vectormath::cross(p.dir, E2);

        const auto det_inv = 1 / vectormath::dot(P, E1);

        const auto v = vectormath::dot(Q, p.dir) * det_inv;
        const auto u = vectormath::dot(P, TT) * det_inv;

        return v >= 0 && u >= 0 && (u + v) <= 1 ? std::make_optional(vectormath::dot(Q, E2) * det_inv)
                                                : std::nullopt;
    }

private:
    std::array<std::array<double, 3>, 3> m_vertices;
};
}
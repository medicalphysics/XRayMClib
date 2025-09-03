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

#pragma once

#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <optional>

namespace dxmc {

class Triangle {
public:
    Triangle(const std::array<double, 3>& first, const std::array<double, 3>& second, const std::array<double, 3>& third)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
    }

    Triangle(const std::array<std::array<double, 3>, 3>& vertices)
        : m_vertices(vertices)
    {
    }

    Triangle(const double* first_element)
    {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto flatIdx = i * 3 + j;
                m_vertices[i][j] = *(first_element + flatIdx);
            }
    }

    auto operator<=>(const Triangle& other) const = default;

    void translate(const std::array<double, 3>& dist)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }

    void mirror(const std::array<double, 3>& point)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            const auto d = vectormath::subtract(m_vertices[i], point);
            m_vertices[i] = vectormath::add(m_vertices[i], vectormath::scale(d, -2.0));
        }
    }

    void mirror(const double value, const std::uint_fast32_t dim)
    {
        for (std::uint_fast32_t i = 0; i < 3; ++i) {
            const auto d = m_vertices[i][dim] - value;
            m_vertices[i][dim] += d * -2.0;
        }
    }

    void scale(double scale)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] *= scale;
            }
        });
    }

    void rotate(double radians, const std::array<double, 3>& axis)
    {
        std::transform(std::execution::unseq, m_vertices.cbegin(), m_vertices.cend(), m_vertices.begin(), [&](auto& vert) {
            return vectormath::rotate(vert, axis, radians);
        });
    }

    std::array<double, 3> planeVector() const noexcept
    {
        const auto a = vectormath::subtract(m_vertices[1], m_vertices[0]);
        const auto b = vectormath::subtract(m_vertices[2], m_vertices[0]);
        const auto n = vectormath::cross(a, b);
        return vectormath::normalized(n);
    }

    const std::array<std::array<double, 3>, 3>& vertices() const
    {
        return m_vertices;
    }

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

    double area() const
    { // Triangle of ABC, area = |AB x AC|/2
        const auto AB = vectormath::subtract(m_vertices[1], m_vertices[0]);
        const auto AC = vectormath::subtract(m_vertices[2], m_vertices[0]);
        return vectormath::length(vectormath::cross(AB, AC)) / 2;
    }

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
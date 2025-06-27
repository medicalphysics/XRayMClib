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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {

struct SphereSamplingRectangularField {
    // Sampling of uniform points on a sphere constrained by a rectangular field
    // Based on: An Area-Preserving Parametrization for Spherical Rectangles by
    // Carlos Ureña1 and Marcos Fajardo2 and Alan King
    SphereSamplingRectangularField(
        const std::array<double, 3>& plane_lower_left,
        const std::array<double, 3>& plane_cosinex,
        const std::array<double, 3>& plane_cosiney,
        const std::array<double, 3>& source_pos)
    {
        setData(plane_lower_left, plane_cosinex, plane_cosiney, source_pos);
    }
    SphereSamplingRectangularField(const std::array<double, 2>& angles)
    {
        setData(angles[0], angles[1]);
    }
    SphereSamplingRectangularField(double collimationHalfAngle_x = 0, double collimationHalfAngle_y = 0)
    {
        setData(collimationHalfAngle_x, collimationHalfAngle_y);
    }
    SphereSamplingRectangularField(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }
    SphereSamplingRectangularField(double x_min, double y_min, double x_max, double y_max)
    {
        setData(x_min, y_min, x_max, y_max);
    }
    void setData(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }
    void setData(double x_min_ang, double y_min_ang, double x_max_ang, double y_max_ang)
    {
        const double x_min = std::tan(x_min_ang);
        const double x_max = std::tan(x_max_ang);
        const double y_min = std::tan(y_min_ang);
        const double y_max = std::tan(y_max_ang);
        std::array<double, 3> cos_x = { x_max - x_min, 0, 0 };
        std::array<double, 3> cos_y = { 0, y_max - y_min, 0 };
        std::array<double, 3> p = { x_min, y_min, 1 };
        std::array<double, 3> o = { 0, 0, 0 };
        setData(p, cos_x, cos_y, o);
    }

    void setData(double collimationHalfAngle_x, double collimationHalfAngle_y)
    {
        const double x = std::tan(collimationHalfAngle_x);
        const double y = std::tan(collimationHalfAngle_y);
        std::array<double, 3> cos_x = { 2 * x, 0, 0 };
        std::array<double, 3> cos_y = { 0, 2 * y, 0 };
        std::array<double, 3> p = { -x, -y, 1 };
        std::array<double, 3> o = { 0, 0, 0 };
        setData(p, cos_x, cos_y, o);
    }

    void setData(
        const std::array<double, 3>& plane_lower_left,
        const std::array<double, 3>& plane_cosinex,
        const std::array<double, 3>& plane_cosiney,
        const std::array<double, 3>& source_pos)
    {
        const double exl = vectormath::length(plane_cosinex);
        const double eyl = vectormath::length(plane_cosiney);
        m_o = source_pos;
        m_x = vectormath::normalized(plane_cosinex);
        m_y = vectormath::normalized(plane_cosiney);
        m_z = vectormath::cross(m_x, m_y);

        const auto d = vectormath::subtract(plane_lower_left, source_pos);
        m_z0 = vectormath::dot(d, m_z);
        if (m_z0 > 0) {
            m_z = vectormath::scale(m_z, -1.0);
            m_z0 *= -1;
        }

        m_z0sq = m_z0 * m_z0;
        m_x0 = vectormath::dot(d, m_x);
        m_y0 = vectormath::dot(d, m_y);
        m_x1 = m_x0 + exl;
        m_y1 = m_y0 + eyl;
        m_y0sq = m_y0 * m_y0;
        m_y1sq = m_y1 * m_y1;

        std::array<double, 3> v00 = { m_x0, m_y0, m_z0 };
        std::array<double, 3> v01 = { m_x0, m_y1, m_z0 };
        std::array<double, 3> v10 = { m_x1, m_y0, m_z0 };
        std::array<double, 3> v11 = { m_x1, m_y1, m_z0 };

        auto n0 = vectormath::normalized(vectormath::cross(v00, v10));
        auto n1 = vectormath::normalized(vectormath::cross(v10, v11));
        auto n2 = vectormath::normalized(vectormath::cross(v11, v01));
        auto n3 = vectormath::normalized(vectormath::cross(v01, v00));

        double g0 = std::acos(-vectormath::dot(n0, n1));
        double g1 = std::acos(-vectormath::dot(n1, n2));
        double g2 = std::acos(-vectormath::dot(n2, n3));
        double g3 = std::acos(-vectormath::dot(n3, n0));
        m_b0 = n0[2];
        m_b1 = n2[2];
        m_b0sq = m_b0 * m_b0;

        constexpr auto pi2 = std::numbers::pi_v<double> * 2;
        m_k = pi2 - g2 - g3;
        m_S = g0 + g1 - m_k;
        return;
    }

    std::array<double, 3> operator()(RandomState& state) const
    {
        const auto u = state.randomUniform();
        const auto au = u * m_S + m_k;
        const auto fu = (std::cos(au) * m_b0 - m_b1) / std::sin(au);
        const auto cu = std::clamp(1.0 / std::sqrt(fu * fu + m_b0sq), 0.0, 1.0) * (fu > 0 ? 1 : -1);
        const auto xu = std::clamp(-(cu * m_z0) / std::sqrt(1.0 - cu * cu), m_x0, m_x1);

        const auto v = state.randomUniform();
        const auto d = std::sqrt(xu * xu + m_z0sq);
        const auto h0 = m_y0 / std::sqrt(d * d + m_y0sq);
        const auto h1 = m_y1 / std::sqrt(d * d + m_y1sq);
        const auto hv = h0 + v * (h1 - h0);
        const auto hv2 = hv * hv;
        constexpr auto hv2lim = 1.0 - std::numeric_limits<double>::epsilon();
        const auto yv = hv2 < hv2lim ? (hv * d) / std::sqrt(1 - hv2) : m_y1;

        std::array<double, 3> dir = {
            m_o[0] + xu * m_x[0] + yv * m_y[0] + m_z0 * m_z[0],
            m_o[1] + xu * m_x[1] + yv * m_y[1] + m_z0 * m_z[1],
            m_o[2] + xu * m_x[2] + yv * m_y[2] + m_z0 * m_z[2]
        };
        vectormath::normalize(dir);
        return dir;
    }
    std::array<double, 3> m_o;
    std::array<double, 3> m_x;
    std::array<double, 3> m_y;
    std::array<double, 3> m_z;
    double m_z0;
    double m_z0sq;
    double m_x0;
    double m_y0;
    double m_y0sq;
    double m_x1;
    double m_y1;
    double m_y1sq;
    double m_b0;
    double m_b1;
    double m_b0sq;
    double m_k;
    double m_S;
};
struct SphereSamplingSquareField {
    // Sampling of uniform points on a sphere constrained by a rectangular field
    // Based on random sampling of points on a whole sphere by sampling x, y, z in [-1, 1] and normalizing the vector
    // Uses simple random sampling of x, y, z inside a constrained box and rejecting points outside field
    // perhaps https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone is better?

    SphereSamplingSquareField(double half_angle = 0)
    {
        setData(half_angle);
    }

    void setData(double half_angle)
    {
        border = std::max(std::sin(half_angle), std::numeric_limits<double>::min());
        const auto sinz = std::sqrt(2 * border * border);
        cosz = std::sqrt(1 - sinz * sinz);
    }

    std::array<double, 3> operator()(RandomState& state) const
    {
        std::array<double, 3> dir;
        double x, y, z;

        do {
            const auto costheta = state.randomUniform(cosz, 1.0);
            const auto sintheta = std::sqrt(1 - costheta * costheta);
            constexpr auto pi2 = std::numbers::pi_v<double> * 2;
            double phi = state.randomUniform(0.0, pi2);
            double cosphi = std::cos(phi);
            double sinphi = std::sin(phi);

            x = sintheta * cosphi;
            y = sintheta * sinphi;
            z = costheta;

        } while (std::abs(x) > border || std::abs(y) > border);
        dir = { x, y, z };
        return dir;
    }
    double border = 0;
    double cosz = 1;
};

}
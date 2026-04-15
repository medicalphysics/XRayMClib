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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <numeric>

namespace xraymc {

/**
 * @brief Uniform direction sampler constrained to an axis-aligned spherical rectangle.
 *
 * Samples unit directions uniformly over the solid angle subtended by a rectangular
 * field whose edges are aligned with the x- and y-axes and whose centre is on the +z
 * axis. The field is specified by its angular half-extents or by explicit min/max angles
 * in x and y.
 *
 * Implements the area-preserving parametrisation described in:
 *   "An Area-Preserving Parametrization for Spherical Rectangles"
 *   C. Ureña, M. Fajardo, A. King (2013).
 *
 * The algorithm maps two independent uniform variates (u, v) ∈ [0,1)² to a point on
 * the spherical cap via an analytic CDF inversion — no rejection sampling is needed.
 *
 * For a field that is not axis-aligned, use `SphereSamplingRectangularGeneralField`.
 */
struct SphereSamplingRectangularField {
    /**
     * @brief Constructs a symmetric field from a pair of half-angles {x, y} [rad].
     * @param angles  `{collimationHalfAngle_x, collimationHalfAngle_y}` in radians.
     */
    SphereSamplingRectangularField(const std::array<double, 2>& angles)
    {
        setData(angles[0], angles[1]);
    }

    /**
     * @brief Constructs a symmetric field from separate x and y half-angles [rad].
     * @param collimationHalfAngle_x  Half-angle in the x direction [rad]. Default: 0.
     * @param collimationHalfAngle_y  Half-angle in the y direction [rad]. Default: 0.
     */
    SphereSamplingRectangularField(double collimationHalfAngle_x = 0, double collimationHalfAngle_y = 0)
    {
        setData(collimationHalfAngle_x, collimationHalfAngle_y);
    }

    /**
     * @brief Constructs an asymmetric field from a four-element angle array [rad].
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds in radians.
     */
    SphereSamplingRectangularField(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }

    /**
     * @brief Constructs an asymmetric field from explicit min/max angles [rad].
     * @param x_min  Minimum x angle [rad].
     * @param y_min  Minimum y angle [rad].
     * @param x_max  Maximum x angle [rad].
     * @param y_max  Maximum y angle [rad].
     */
    SphereSamplingRectangularField(double x_min, double y_min, double x_max, double y_max)
    {
        setData(x_min, y_min, x_max, y_max);
    }

    /**
     * @brief Sets an asymmetric field from a four-element angle array [rad].
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds in radians.
     */
    void setData(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }

    /**
     * @brief Sets a symmetric field from x and y half-angles [rad].
     *
     * Expands to `setData(-hx, -hy, +hx, +hy)`.
     *
     * @param collimationHalfAngle_x  Half-angle in the x direction [rad].
     * @param collimationHalfAngle_y  Half-angle in the y direction [rad].
     */
    void setData(double collimationHalfAngle_x, double collimationHalfAngle_y)
    {
        setData(-collimationHalfAngle_x, -collimationHalfAngle_y, collimationHalfAngle_x, collimationHalfAngle_y);
    }

    /**
     * @brief Samples a uniformly distributed unit direction within the rectangular field.
     *
     * Draws two independent uniform variates u and v, inverts the area-preserving CDF
     * analytically to obtain the x- and y-direction tangent components, then normalises
     * to produce a unit vector pointing into the +z hemisphere.
     *
     * @param state  Per-thread PRNG state.
     * @return Unit direction vector {x, y, z} with z > 0.
     */
    std::array<double, 3> operator()(RandomState& state) const
    {
        const auto u = state.randomUniform();
        const auto au = u * m_S + m_k;
        const auto fu = (std::cos(au) * m_b0 - m_b1) / std::sin(au);
        const auto cu = 1 / std::sqrt(fu * fu + m_b0sq) * (fu > 0 ? 1 : -1);
        const auto xu = cu / std::sqrt(1 - cu * cu);

        const auto v = state.randomUniform();
        const auto d = std::sqrt(xu * xu + 1);
        const auto h0 = m_y0 / std::sqrt(d * d + m_y0sq);
        const auto h1 = m_y1 / std::sqrt(d * d + m_y1sq);
        const auto hv = h0 + v * (h1 - h0);
        const auto hv2 = hv * hv;
        constexpr auto hv2lim = 1 - std::numeric_limits<double>::epsilon();
        const auto yv = hv2 < hv2lim ? (hv * d) / std::sqrt(1 - hv2) : m_y1;

        const auto norm = 1 / std::sqrt(xu * xu + yv * yv + 1);
        std::array<double, 3> dir = { xu * norm, yv * norm, norm };
        return dir;
    }

    /**
     * @brief Sets an asymmetric field from explicit angular bounds and precomputes sampling constants.
     *
     * Converts angular bounds to tangent-plane coordinates, computes the four bounding-plane
     * normals of the spherical rectangle, then derives the area-preserving parametrisation
     * constants m_S and m_k used by `operator()`.
     *
     * Very small fields (|max − min| < 1e-6 rad in either dimension) are expanded slightly
     * to avoid numerical degenerate cases.
     *
     * @param x_min_ang  Minimum x angular bound [rad].
     * @param y_min_ang  Minimum y angular bound [rad].
     * @param x_max_ang  Maximum x angular bound [rad].
     * @param y_max_ang  Maximum y angular bound [rad].
     */
    void setData(double x_min_ang, double y_min_ang, double x_max_ang, double y_max_ang)
    {
        auto away_zero = [](double& min, double& max) -> void {
            constexpr double minang = 0.000001; // 100 * std::numeric_limits<double>::epsilon();
            if (std::abs(max - min) < minang) {
                const auto mid = std::midpoint(min, max);
                min = mid - minang;
                max = mid + minang;
            }
        };
        away_zero(x_min_ang, x_max_ang);
        away_zero(y_min_ang, y_max_ang);

        const auto exl = std::tan(x_max_ang) - std::tan(x_min_ang);
        const auto eyl = std::tan(y_max_ang) - std::tan(y_min_ang);

        // const std::array<double, 3> d = { std::tan(x_min_ang), std::tan(y_min_ang), 1 };

        m_x0 = std::tan(x_min_ang);
        m_y0 = std::tan(y_min_ang);
        m_x1 = m_x0 + exl;
        m_y1 = m_y0 + eyl;
        m_y0sq = m_y0 * m_y0;
        m_y1sq = m_y1 * m_y1;

        std::array<double, 3> v00 = { m_x0, m_y0, -1 };
        std::array<double, 3> v01 = { m_x0, m_y1, -1 };
        std::array<double, 3> v10 = { m_x1, m_y0, -1 };
        std::array<double, 3> v11 = { m_x1, m_y1, -1 };

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

    double m_x0;    ///< Tangent-plane x coordinate of the left edge.
    double m_y0;    ///< Tangent-plane y coordinate of the bottom edge.
    double m_y0sq;  ///< m_y0².
    double m_x1;    ///< Tangent-plane x coordinate of the right edge.
    double m_y1;    ///< Tangent-plane y coordinate of the top edge.
    double m_y1sq;  ///< m_y1².
    double m_b0;    ///< z-component of the bottom bounding-plane normal n0.
    double m_b1;    ///< z-component of the top bounding-plane normal n2.
    double m_b0sq;  ///< m_b0².
    double m_k;     ///< CDF offset constant (2π − g2 − g3).
    double m_S;     ///< Solid angle of the spherical rectangle [sr].
};

/**
 * @brief Uniform direction sampler constrained to an arbitrarily oriented spherical rectangle.
 *
 * Generalises `SphereSamplingRectangularField` to support rectangular fields that are
 * not aligned with the coordinate axes. The field is defined by a lower-left corner in
 * world space, two orthogonal edge vectors (`plane_cosinex`, `plane_cosiney`), and the
 * source position. Sampled world-space directions are returned as unit vectors.
 *
 * Implements the same area-preserving CDF-inversion algorithm as
 * `SphereSamplingRectangularField`, extended to arbitrary orientation by transforming
 * into a local frame aligned with the plane normal.
 *
 * Reference: "An Area-Preserving Parametrization for Spherical Rectangles",
 * C. Ureña, M. Fajardo, A. King (2013).
 */
struct SphereSamplingRectangularGeneralField {
    /**
     * @brief Constructs the sampler from an explicit world-space rectangle description.
     *
     * @param plane_lower_left  World-space position of the rectangle's lower-left corner.
     * @param plane_cosinex     Edge vector along the local x axis (length = field width).
     * @param plane_cosiney     Edge vector along the local y axis (length = field height).
     *                          Must be orthogonal to @p plane_cosinex.
     * @param source_pos        World-space source (focal spot) position.
     */
    SphereSamplingRectangularGeneralField(
        const std::array<double, 3>& plane_lower_left,
        const std::array<double, 3>& plane_cosinex,
        const std::array<double, 3>& plane_cosiney,
        const std::array<double, 3>& source_pos)
    {
        setData(plane_lower_left, plane_cosinex, plane_cosiney, source_pos);
    }

    /**
     * @brief Constructs a symmetric axis-aligned field from a pair of half-angles {x, y} [rad].
     * @param angles  `{collimationHalfAngle_x, collimationHalfAngle_y}` in radians.
     */
    SphereSamplingRectangularGeneralField(const std::array<double, 2>& angles)
    {
        setData(angles[0], angles[1]);
    }

    /**
     * @brief Constructs a symmetric axis-aligned field from x and y half-angles [rad].
     * @param collimationHalfAngle_x  Half-angle in the x direction [rad]. Default: 0.01.
     * @param collimationHalfAngle_y  Half-angle in the y direction [rad]. Default: 0.01.
     */
    SphereSamplingRectangularGeneralField(double collimationHalfAngle_x = .01, double collimationHalfAngle_y = .01)
    {
        setData(collimationHalfAngle_x, collimationHalfAngle_y);
    }

    /**
     * @brief Constructs an asymmetric axis-aligned field from a four-element angle array [rad].
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds in radians.
     */
    SphereSamplingRectangularGeneralField(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }

    /**
     * @brief Constructs an asymmetric axis-aligned field from explicit min/max angles [rad].
     * @param x_min  Minimum x angle [rad].
     * @param y_min  Minimum y angle [rad].
     * @param x_max  Maximum x angle [rad].
     * @param y_max  Maximum y angle [rad].
     */
    SphereSamplingRectangularGeneralField(double x_min, double y_min, double x_max, double y_max)
    {
        setData(x_min, y_min, x_max, y_max);
    }

    /**
     * @brief Sets an asymmetric axis-aligned field from a four-element angle array [rad].
     * @param angles  `{x_min, y_min, x_max, y_max}` angular bounds in radians.
     */
    void setData(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }

    /**
     * @brief Sets an asymmetric axis-aligned field from explicit angular bounds [rad].
     *
     * Converts angular bounds to tangent-plane coordinates, constructs a world-space
     * rectangle on the z = 1 plane centred on the origin, and delegates to the
     * world-space `setData` overload.
     *
     * @param x_min_ang  Minimum x angular bound [rad].
     * @param y_min_ang  Minimum y angular bound [rad].
     * @param x_max_ang  Maximum x angular bound [rad].
     * @param y_max_ang  Maximum y angular bound [rad].
     */
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

    /**
     * @brief Sets a symmetric axis-aligned field from x and y half-angles [rad].
     *
     * Constructs a world-space rectangle on the z = 1 plane of full size
     * 2·tan(hx) × 2·tan(hy) centred on the origin and delegates to the
     * world-space `setData` overload.
     *
     * @param collimationHalfAngle_x  Half-angle in the x direction [rad].
     * @param collimationHalfAngle_y  Half-angle in the y direction [rad].
     */
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

    /**
     * @brief Sets the sampler from a world-space rectangle and source position.
     *
     * Computes the local frame (x̂, ŷ, ẑ) aligned with the plane, projects the
     * rectangle corners into that frame, builds the four bounding-plane normals,
     * and derives the area-preserving parametrisation constants m_S and m_k.
     * The plane normal direction is flipped if necessary so that m_z0 < 0
     * (i.e. the source lies on the positive-normal side of the plane).
     *
     * @param plane_lower_left  World-space lower-left corner of the rectangle.
     * @param plane_cosinex     Edge vector along local x (length = field width).
     * @param plane_cosiney     Edge vector along local y (length = field height).
     * @param source_pos        World-space source position.
     */
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

    /**
     * @brief Samples a uniformly distributed unit direction within the rectangular field.
     *
     * Draws two independent uniform variates u and v, inverts the area-preserving CDF
     * analytically in the local frame, then transforms the result back to world space
     * and normalises. Clamps intermediate values to avoid numerical overflow at the
     * field boundary.
     *
     * @param state  Per-thread PRNG state.
     * @return Unit direction vector in world space pointing into the rectangular field.
     */
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

    std::array<double, 3> m_o; ///< Source (origin) position in world space.
    std::array<double, 3> m_x; ///< Local x-axis unit vector (along plane_cosinex).
    std::array<double, 3> m_y; ///< Local y-axis unit vector (along plane_cosiney).
    std::array<double, 3> m_z; ///< Local z-axis unit vector (plane normal, pointing toward source).
    double m_z0;    ///< Signed distance from source to the plane along m_z (always ≤ 0).
    double m_z0sq;  ///< m_z0².
    double m_x0;    ///< Local-frame x coordinate of the left edge.
    double m_y0;    ///< Local-frame y coordinate of the bottom edge.
    double m_y0sq;  ///< m_y0².
    double m_x1;    ///< Local-frame x coordinate of the right edge.
    double m_y1;    ///< Local-frame y coordinate of the top edge.
    double m_y1sq;  ///< m_y1².
    double m_b0;    ///< z-component of the bottom bounding-plane normal n0.
    double m_b1;    ///< z-component of the top bounding-plane normal n2.
    double m_b0sq;  ///< m_b0².
    double m_k;     ///< CDF offset constant (2π − g2 − g3).
    double m_S;     ///< Solid angle of the spherical rectangle [sr].
};

/**
 * @brief Uniform direction sampler constrained to a square field via rejection sampling.
 *
 * Samples directions uniformly within the square cap defined by |x| ≤ sin(half_angle)
 * and |y| ≤ sin(half_angle) on the unit sphere. The algorithm draws directions
 * uniformly from the enclosing circular cone (radius sin(half_angle)·√2) and rejects
 * those whose projected x or y component exceeds the square boundary.
 *
 * @note This is a rejection sampler. For small or very large half-angles the acceptance
 * rate approaches π/4 ≈ 0.785 (the ratio of a square to its circumscribed circle).
 * For fields that are axis-aligned rectangular (not necessarily square) consider using
 * `SphereSamplingRectangularField` which is rejection-free.
 */
struct SphereSamplingSquareField {
    /**
     * @brief Constructs the sampler for the given square half-angle.
     * @param half_angle  Half-width of the square field [rad]. Default: 0.01.
     */
    SphereSamplingSquareField(double half_angle = .01)
    {
        setData(half_angle);
    }

    /**
     * @brief Sets the square half-angle and precomputes the enclosing-cone cosine.
     *
     * Computes `border = sin(half_angle)` (clamped away from zero) and
     * `cosz = cos(asin(border·√2))` — the cosine of the enclosing circular cone
     * used as the outer bound for the rejection loop.
     *
     * @param half_angle  Half-width of the square field [rad].
     */
    void setData(double half_angle)
    {
        border = std::max(std::sin(half_angle), std::numeric_limits<double>::min());
        const auto sinz = std::sqrt(2 * border * border);
        cosz = std::sqrt(1 - sinz * sinz);
    }

    /**
     * @brief Samples a uniformly distributed unit direction within the square field.
     *
     * Repeatedly draws directions from the enclosing circular cone until both |x| and
     * |y| are within `border`. Thread-safe when @p state is not shared.
     *
     * @param state  Per-thread PRNG state.
     * @return Unit direction vector {x, y, z} with |x| ≤ border and |y| ≤ border.
     */
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
    double border = 0; ///< sin(half_angle) — maximum allowed |x| and |y| component.
    double cosz = 1;   ///< cos of the enclosing cone angle (cos(asin(border·√2))).
};

}
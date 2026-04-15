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

#include "xraymc/constants.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {

/**
 * @brief Uniform direction sampler constrained to a circular cone around the +z axis.
 *
 * Samples directions uniformly over the spherical cap defined by a half-opening
 * angle θ from the +z axis (i.e. the set of unit vectors with cos θ ≥ cos(angle)).
 *
 * The method exploits the fact that for a uniform distribution on a sphere the
 * z-coordinate is uniformly distributed on [−1, 1]. Constraining the cap is
 * therefore equivalent to drawing z uniformly from [cos(angle), 1] and sampling
 * the azimuthal angle φ uniformly from [0, 2π):
 *
 *   z   ~ Uniform(cos(angle), 1)
 *   φ   ~ Uniform(0, 2π)
 *   dir = (√(1−z²)·cos φ,  √(1−z²)·sin φ,  z)
 *
 * This produces a rejection-free, O(1) sampler with no approximation.
 */
struct SphereSamplingCircularField {
    /**
     * @brief Constructs the sampler for the given cone half-angle.
     * @param angle  Half-opening angle of the cone [rad]. A value of 0 produces
     *               directions exactly along +z; π produces the full sphere.
     *               Default: 0.
     */
    SphereSamplingCircularField(double angle = 0)
    {
        setData(angle);
    }

    /**
     * @brief Sets the cone half-angle.
     *
     * Precomputes and stores cos(angle) to avoid repeated trigonometric evaluation
     * inside `operator()`.
     *
     * @param angle  Half-opening angle of the cone [rad].
     */
    void setData(double angle)
    {
        m_cosz = std::cos(angle);
    }

    /**
     * @brief Samples a uniformly distributed unit direction within the cone.
     *
     * @param state  Per-thread PRNG state.
     * @return Unit direction vector {x, y, z} with z ∈ [cos(angle), 1].
     */
    std::array<double, 3> operator()(RandomState& state) const
    {
        const auto z = state.randomUniform(m_cosz, 1.0);
        const auto phi = state.randomUniform(2 * PI_VAL());
        const auto cosphi = std::cos(phi);
        const auto sinphi = std::sin(phi);
        const auto sinz = std::sqrt(1 - z * z);
        std::array<double, 3> dir = { sinz * cosphi, sinz * sinphi, z };
        return dir;
    }

    double m_cosz = 0; ///< Cosine of the cone half-angle; lower bound for the z-coordinate draw.
};
}
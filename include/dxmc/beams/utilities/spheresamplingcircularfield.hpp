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

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {
struct SphereSamplingCircularField {
    // Sampling of uniform points on a sphere constrained by a rectangular field
    // Based on random sampling of points on a whole sphere by sampling x, y, z in [-1, 1] and normalizing the vector
    // Uses simple random sampling of x, y, z inside a constrained box and rejecting points outside field
    // perhaps https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone is better?

    SphereSamplingCircularField(double angle = 0)
    {
        setData(angle);
    }

    void setData(double angle)
    {
        m_cosz = std::cos(angle);
    }
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

    double m_cosz = 0;
};
}
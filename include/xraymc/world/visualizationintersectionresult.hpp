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

#include <array>

namespace xraymc {

/**
 * @brief Result of a KD-tree ray intersection query for visualization rendering.
 *
 * Extends the information returned by a transport intersection with a surface normal
 * and a scalar value (e.g. absorbed dose) used for shading in VisualizeWorld. A
 * result is considered valid when `intersectionValid` is true.
 *
 * @tparam U The world item type pointed to (e.g. `std::variant<Us...>`).
 */
template <typename U>
struct VisualizationIntersectionResult {
    std::array<double, 3> normal = { 0, 0, 0 }; ///< Unit outward surface normal at the intersection point.
    double intersection = 0;                      ///< Ray parameter t at the intersection point in cm.
    const U* item = nullptr;                      ///< Pointer to the intersected item, or nullptr on miss.
    bool rayOriginIsInsideItem = false;           ///< True if the ray origin is inside the intersected item.
    bool intersectionValid = false;               ///< True when a valid intersection was found.
    double value = 0;                             ///< Scalar item value (e.g. dose in keV/g) used for color mapping.

    /// @brief Returns true if a valid intersection was found.
    inline bool valid() const
    {
        return intersectionValid;
    }
};
}
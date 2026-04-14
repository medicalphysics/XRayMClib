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

namespace xraymc {

/**
 * @brief Result of a KD-tree ray intersection query for transport simulation.
 *
 * Returned by KDTree, KDTreeFlat, and the mesh KD-tree variants after a ray–world
 * intersection test. A result is considered valid when `item` is non-null.
 *
 * @tparam U The world item type pointed to (e.g. `std::variant<Us...>` or a triangle type).
 */
template <typename U>
struct KDTreeIntersectionResult {
    U* item = nullptr;             ///< Pointer to the closest intersected item, or nullptr on miss.
    double intersection = 0;       ///< Ray parameter t at the intersection point in cm.
    bool rayOriginIsInsideItem = false; ///< True if the ray origin is inside the intersected item.

    /// @brief Returns true if an intersection was found (item is non-null).
    bool valid() const
    {
        return item != nullptr;
    }
};
}
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
 * @brief Result of a ray–world-item intersection test used during particle transport.
 *
 * Returned by each world item's `intersect()` method. A result is considered valid
 * when `intersectionValid` is true. The `rayOriginIsInsideItem` flag is used by the
 * transport loop to determine whether the particle is entering or exiting the item.
 */
struct WorldIntersectionResult {
    double intersection = 0;          ///< Ray parameter t at the intersection point in cm.
    bool rayOriginIsInsideItem = false; ///< True if the ray origin is inside the item at the time of the query.
    bool intersectionValid = false;   ///< True when a valid intersection was found.

    /// @brief Returns true if a valid intersection was found.
    inline bool valid() const
    {
        return intersectionValid;
    }
};
}
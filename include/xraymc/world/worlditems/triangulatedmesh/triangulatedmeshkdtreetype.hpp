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

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "xraymc/particle.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"

#include <array>

namespace xraymc {

/**
 * @brief Concept constraining the element type used by MeshKDTree and MeshKDTreeFlat.
 *
 * A type @p U satisfies MeshKDTreeType if it provides:
 * - Three-way comparison (`<=>`) for ordering.
 * - `translate(vec)` — displaces the element by a 3-D vector in cm.
 * - `scale(s)` — uniformly scales all coordinates by a scalar factor.
 * - `intersect(p)` — returns `std::optional<double>` (the ray parameter t) for a
 *   Möller–Trumbore or equivalent ray–element intersection test.
 * - `center()` — returns the element centroid as `std::array<double,3>` in cm.
 * - `AABB()` — returns the axis-aligned bounding box as `std::array<double,6>`
 *   ({xmin,ymin,zmin,xmax,ymax,zmax}) in cm.
 * - `planeVector()` — returns the unit surface normal as `std::array<double,3>`.
 */
template <typename U>
concept MeshKDTreeType = requires(U u, Particle p, std::array<double, 3> vec, double scale) {
    u <=> u;
    u.translate(vec);
    u.scale(scale);
    {
        u.intersect(p)
    } -> std::same_as<std::optional<double>>;
    {
        u.center()
    } -> std::convertible_to<std::array<double, 3>>;
    {
        u.AABB()
    } -> std::convertible_to<std::array<double, 6>>;
    {
        u.planeVector()
    } -> std::convertible_to<std::array<double, 3>>;
};
}
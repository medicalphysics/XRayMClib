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

#include "xraymc/world/collisions/meshcollision.hpp"
#include "xraymc/world/world.hpp"

namespace xraymc {

class CollisionTests {
public:
    // CollisionTest() = default;
    template <WorldItemType F, WorldItemType U>
    static bool test(const U& item1, const F& item2)
    {
        const auto aabb1 = item1.AABB();
        const auto aabb2 = item2.AABB();

        return basicshape::AABB::collide(aabb1, aabb2) ? collision::collide(item1, item2) : false;
    }

protected:
private:
};
}

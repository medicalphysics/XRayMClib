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

#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/basicshapes/sphere.hpp"

#include "xraymc/xraymcrandom.hpp"

#include <iostream>

bool testSphereIntersection()
{
    bool success = true;
    xraymc::Particle p;

    double radii = 2;
    std::array<double, 3> center = { 0, 0, 0 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::sphere::intersect(p, center, radii);
    success = success && !t.valid();

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray sphere intersections." << std::endl;

    return success;
}

bool testAABBIntersection()
{
    bool success = true;
    xraymc::Particle p;

    std::array<double, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = xraymc::basicshape::AABB::intersect(p, aabb);
    success = success && !t.valid();

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray AABB intersections. " << std::endl;

    return success;
}

bool testCylindarIntersection()
{
    bool success = true;
    xraymc::Particle p;

    auto cylinder = xraymc::basicshape::cylinder::Cylinder({ 0, 0, 0 }, { 0, 0, 1 }, 2, 2);

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -8000, 0, 0 };
    p.dir = { 1, 0, 0 };
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 8000, 0 };
    p.dir = { 0, -1, 0 };
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height - 1 };
    p.dir = { 1, 0, 1 };
    xraymc::vectormath::normalize(p.dir);
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && std::abs(t.intersection - std::sqrt(2.0)) <= std::numeric_limits<double>::epsilon() * 10;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height + 1 };
    p.dir = { 1, 0, 1 };
    xraymc::vectormath::normalize(p.dir);
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && std::abs(t.intersection - 3 * std::sqrt(2.0)) <= std::numeric_limits<double>::epsilon() * 10;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height - 1 };
    p.dir = { -1, 0, -1 };
    xraymc::vectormath::normalize(p.dir);
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && !t.valid();

    p.pos = { -1, -1, 0 };
    p.dir = { 0, 0, -1 };
    xraymc::vectormath::normalize(p.dir);
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == cylinder.half_height;

    p.pos = { -1, -1, -cylinder.half_height * 2 };
    p.dir = { 0, 0, -1 };
    xraymc::vectormath::normalize(p.dir);
    t = xraymc::basicshape::cylinder::intersect(p, cylinder);
    success = success && !t.valid();

    xraymc::RandomState state;
    for (std::size_t i = 0; i < 1E6; ++i) {
        std::array<double, 3> pos_init = {
            state.randomUniform(-cylinder.radius, cylinder.radius),
            state.randomUniform(-cylinder.radius, cylinder.radius),
            state.randomUniform(-cylinder.half_height - 1, cylinder.half_height + 1),
        };

        pos_init = {
            -0.216281652,
            1.74022770,
            -1.04047060
        };

        const auto pos = xraymc::vectormath::add(pos_init, cylinder.center);
        const bool inside = xraymc::basicshape::cylinder::pointInside(pos, cylinder);

        bool inside_t = pos_init[0] * pos_init[0] + pos_init[1] * pos_init[1] <= cylinder.radius * cylinder.radius;
        inside_t = inside_t && -cylinder.half_height <= pos_init[2] && pos_init[2] <= cylinder.half_height;

        success = success && inside == inside_t;
        if (!success)
            auto test = false;
    }

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray cylinder intersections." << std::endl;

    return success;
}

bool testCylinderForwardIntersection()
{

    bool success = true;

    const double radii = 4;
    const std::array<double, 3> center = { 8, 0, 0 };
    const double half_height = 6;

    const auto cylinder = xraymc::basicshape::cylinder::Cylinder(center, { 0, 0, 1 }, radii, half_height);

    xraymc::Particle p;

    std::optional<std::array<double, 2>> res;

    std::array dummy = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
    std::array<double, 2> t;

    p.pos = { -10, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = xraymc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 14;
    success = success && t[1] == 22;

    p.pos = { 6, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = xraymc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    p.pos = { 8, 0, 0 };
    p.dir = { 0, 0, 1 };
    res = xraymc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    xraymc::RandomState state;
    for (std::size_t i = 0; i < 1e6; ++i) {
        for (std::size_t i = 0; i < 3; ++i) {
            p.pos[i] = state.randomUniform(center[i] - radii * 2, center[i] + radii * 2);
            p.dir[i] = state.randomUniform(-1.0, 1.0);
        }
        xraymc::vectormath::normalize(p.dir);
        auto r = xraymc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
        if (r) {
            auto t2 = r.value();
            if (t2[0] < 0 || t2[1] < 0 || t2[1] <= t2[0])
                success = false;
        }
    }

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray cylinder forward intersections." << std::endl;

    return success;
}

int main()
{
    std::cout << "Testing ray intersection on basic shapes\n";

    bool success = true;

    success = success && testSphereIntersection();

    success = success && testAABBIntersection();

    success = success && testCylindarIntersection();

    success = success && testCylinderForwardIntersection();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
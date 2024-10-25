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

Copyright 2022 Erlend Andersen
*/

#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/basicshapes/sphere.hpp"

#include "dxmc/dxmcrandom.hpp"

#include <iostream>

bool testSphereIntersection()
{
    bool success = true;
    dxmc::Particle p;

    double radii = 2;
    std::array<double, 3> center = { 0, 0, 0 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
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
    dxmc::Particle p;

    std::array<double, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
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
    dxmc::Particle p;

    auto cylinder = dxmc::basicshape::cylinder::Cylinder({ 0, 0, 0 }, { 0, 0, 1 }, 2, 2);

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -8000, 0, 0 };
    p.dir = { 1, 0, 0 };
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 8000, 0 };
    p.dir = { 0, -1, 0 };
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height - 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && std::abs(t.intersection - std::sqrt(2.0)) <= std::numeric_limits<double>::epsilon() * 10;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height + 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && std::abs(t.intersection - 3 * std::sqrt(2.0)) <= std::numeric_limits<double>::epsilon() * 10;

    p.pos = { -1, 0, cylinder.center[0] - cylinder.half_height - 1 };
    p.dir = { -1, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && !t.valid();

    p.pos = { -1, -1, 0 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && t.valid() && t.intersection == cylinder.half_height;

    p.pos = { -1, -1, -cylinder.half_height * 2 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, cylinder);
    success = success && !t.valid();

    dxmc::RandomState state;
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

        const auto pos = dxmc::vectormath::add(pos_init, cylinder.center);
        const bool inside = dxmc::basicshape::cylinder::pointInside(pos, cylinder);

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

    const auto cylinder = dxmc::basicshape::cylinder::Cylinder(center, { 0, 0, 1 }, radii, half_height);

    dxmc::Particle p;

    std::optional<std::array<double, 2>> res;

    std::array dummy = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
    std::array<double, 2> t;

    p.pos = { -10, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 14;
    success = success && t[1] == 22;

    p.pos = { 6, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    p.pos = { 8, 0, 0 };
    p.dir = { 0, 0, 1 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    p.pos = { 4.0839869976043701, -0.55314152687788010, 7.5147741995751858 };
    p.dir = { 0.54546342835878225, -0.13294172771255339, -0.82752410560411160 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    dxmc::RandomState state;
    for (std::size_t i = 0; i < 1e6; ++i) {
        for (std::size_t i = 0; i < 3; ++i) {
            p.pos[i] = state.randomUniform(center[i] - radii * 2, center[i] + radii * 2);
            p.dir[i] = state.randomUniform(-1.0, 1.0);
        }
        dxmc::vectormath::normalize(p.dir);
        auto r = dxmc::basicshape::cylinder::intersectForwardInterval(p, cylinder);
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
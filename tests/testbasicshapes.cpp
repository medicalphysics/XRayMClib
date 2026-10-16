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

#include "testutils.hpp"

#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/basicshapes/cylinderz.hpp"
#include "xraymc/world/basicshapes/sphere.hpp"
#include "xraymc/world/basicshapes/tetrahedron.hpp"
#include "xraymc/world/collision.hpp"

#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

using xraymctest::check;
using xraymctest::checkNear;
using Vec = std::array<double, 3>;

namespace bs = xraymc::basicshape;

constexpr double EPS = 1E-9;

xraymc::Particle ray(const Vec& pos, Vec dir)
{
    xraymc::vectormath::normalize(dir);
    return xraymc::Particle { .pos = pos, .dir = dir, .energy = 60, .weight = 1 };
}

Vec randomDirection(xraymc::RandomState& state)
{
    const auto cost = state.randomUniform(-1.0, 1.0);
    const auto phi = state.randomUniform(2 * std::numbers::pi);
    const auto sint = std::sqrt(1 - cost * cost);
    return { sint * std::cos(phi), sint * std::sin(phi), cost };
}

// Generic consistency test between a shape's intersect() and pointInside().
// For random rays we march along the ray and find the first point where the
// inside/outside status changes. If such a crossing exists the shape must report
// a valid intersection at that distance, and the rayOriginIsInsideItem flag must
// agree with pointInside at the ray origin.
struct RayConsistency {
    std::size_t rays = 0;
    std::size_t hits = 0;
    std::size_t missedHits = 0;
    std::size_t falseHits = 0;
    std::size_t wrongDistance = 0;
    std::size_t wrongInsideFlag = 0;
};

RayConsistency rayConsistency(
    const std::function<xraymc::WorldIntersectionResult(const xraymc::Particle&)>& intersect,
    const std::function<bool(const Vec&)>& inside,
    const std::array<double, 6>& samplingBox, double maxLength, std::size_t nRays, std::uint64_t seed)
{
    auto state = xraymctest::seededState(seed);
    RayConsistency res;
    constexpr double step = 1E-3;
    for (std::size_t i = 0; i < nRays; ++i) {
        const Vec pos = {
            state.randomUniform(samplingBox[0], samplingBox[3]),
            state.randomUniform(samplingBox[1], samplingBox[4]),
            state.randomUniform(samplingBox[2], samplingBox[5])
        };
        const auto p = ray(pos, randomDirection(state));
        const bool startInside = inside(p.pos);
        double tcross = -1;
        for (double t = step; t < maxLength; t += step) {
            const auto pt = xraymc::vectormath::add(p.pos, xraymc::vectormath::scale(p.dir, t));
            if (inside(pt) != startInside) {
                tcross = t;
                break;
            }
        }
        const auto r = intersect(p);
        ++res.rays;
        if (tcross > 0) {
            ++res.hits;
            if (!r.valid()) {
                ++res.missedHits;
            } else {
                if (std::abs(r.intersection - tcross) > 2 * step)
                    ++res.wrongDistance;
                if (r.rayOriginIsInsideItem != startInside)
                    ++res.wrongInsideFlag;
            }
        } else if (r.valid() && !startInside) {
            // a hit was reported but the marching found no crossing; allow grazing hits
            // where the chord is shorter than the marching step
            ++res.falseHits;
        }
    }
    return res;
}

void reportConsistency(const RayConsistency& r, const std::string& name, std::size_t allowedFalseHits = 2)
{
    check(r.hits > r.rays / 20, name + ": ray sample includes hits (" + std::to_string(r.hits) + " of " + std::to_string(r.rays) + ")");
    check(r.missedHits == 0, name + ": no missed intersections (" + std::to_string(r.missedHits) + ")");
    check(r.wrongDistance == 0, name + ": intersection distance agrees with ray marching (" + std::to_string(r.wrongDistance) + " wrong)");
    check(r.wrongInsideFlag == 0, name + ": rayOriginIsInsideItem agrees with pointInside (" + std::to_string(r.wrongInsideFlag) + " wrong)");
    check(r.falseHits <= allowedFalseHits, name + ": no false intersections (" + std::to_string(r.falseHits) + ")");
}

void testSphere()
{
    const Vec center = { 0, 0, 0 };
    const double radius = 2;

    // Large distances trigger catastrophic cancellation in the naive quadratic solution
    auto t = bs::sphere::intersect(ray({ 0, 0, -8000 }, { 0, 0, 1 }), center, radius);
    check(t.valid() && !t.rayOriginIsInsideItem, "Sphere: hit from far outside");
    checkNear(t.intersection, 7998, 1E-9, "Sphere: distance from far outside is accurate");

    t = bs::sphere::intersect(ray({ 0, 0, 8000 }, { 0, 0, -1 }), center, radius);
    checkNear(t.intersection, 7998, 1E-9, "Sphere: distance from far outside, opposite side");

    t = bs::sphere::intersect(ray({ 0, 0, 8000 }, { 0, 0, 1 }), center, radius);
    check(!t.valid(), "Sphere: ray pointing away misses");

    t = bs::sphere::intersect(ray({ 0, 3, -10 }, { 0, 0, 1 }), center, radius);
    check(!t.valid(), "Sphere: ray passing beside misses");

    t = bs::sphere::intersect(ray({ 0, 0, -1 }, { 0, 0, 1 }), center, radius);
    check(t.valid() && t.rayOriginIsInsideItem, "Sphere: origin inside, moving towards center is inside");
    checkNear(t.intersection, 3, EPS, "Sphere: exit distance moving towards center");

    t = bs::sphere::intersect(ray({ 0, 0, 1 }, { 0, 0, 1 }), center, radius);
    check(t.valid() && t.rayOriginIsInsideItem, "Sphere: origin inside, moving away from center is inside");
    checkNear(t.intersection, 1, EPS, "Sphere: exit distance moving away from center");

    t = bs::sphere::intersect(ray({ 1, 0, 0 }, { 0, 1, 0 }), center, radius);
    check(t.valid() && t.rayOriginIsInsideItem, "Sphere: origin inside, moving perpendicular to center");
    checkNear(t.intersection, std::sqrt(3.0), EPS, "Sphere: exit distance moving perpendicular to center");

    // Offset center and oblique ray
    const Vec c2 = { 1, -2, 3 };
    const auto p2 = ray({ 1 - 10, -2, 3 + 10 }, { 1, 0, -1 });
    t = bs::sphere::intersect(p2, c2, 1.5);
    checkNear(t.intersection, std::sqrt(200.0) - 1.5, EPS, "Sphere: offset center, oblique ray");

    check(bs::sphere::pointInside({ 0.5, 0.5, 0.5 }, center, radius), "Sphere: pointInside for inner point");
    check(!bs::sphere::pointInside({ 1.5, 1.5, 0.5 }, center, radius), "Sphere: pointInside false for outer point");

    const auto interval = bs::sphere::intersectForwardInterval(ray({ 0, 0, -5 }, { 0, 0, 1 }), center, radius);
    check(interval.has_value(), "Sphere: forward interval for hitting ray");
    if (interval) {
        checkNear(interval.value()[0], 3, EPS, "Sphere: forward interval entry");
        checkNear(interval.value()[1], 7, EPS, "Sphere: forward interval exit");
    }

    const auto viz = bs::sphere::intersectVisualization<int>(ray({ 0, 0, -5 }, { 0, 0, 1 }), center, radius);
    check(viz.valid(), "Sphere: visualization intersection valid");
    checkNear(xraymc::vectormath::length(viz.normal), 1, 1E-12, "Sphere: visualization normal is unit length");
    checkNear(std::abs(viz.normal[2]), 1, 1E-12, "Sphere: visualization normal is radial");

    const auto r = rayConsistency([&](const auto& p) { return bs::sphere::intersect(p, c2, 1.5); },
        [&](const auto& pos) { return bs::sphere::pointInside(pos, c2, 1.5); },
        { c2[0] - 3, c2[1] - 3, c2[2] - 3, c2[0] + 3, c2[1] + 3, c2[2] + 3 }, 12, 3000, 1);
    reportConsistency(r, "Sphere");
}

void testAABB()
{
    const std::array<double, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    auto t = bs::AABB::intersect(ray({ 0, 0, -8000 }, { 0, 0, 1 }), aabb);
    check(t.valid() && !t.rayOriginIsInsideItem, "AABB: hit from far outside");
    checkNear(t.intersection, 7998, EPS, "AABB: distance from far outside");

    t = bs::AABB::intersect(ray({ 0, 0, 8000 }, { 0, 0, 1 }), aabb);
    check(!t.valid(), "AABB: ray pointing away misses");

    t = bs::AABB::intersect(ray({ 0, 0, 1 }, { 0, 0, 1 }), aabb);
    check(t.valid() && t.rayOriginIsInsideItem, "AABB: origin inside is flagged");
    checkNear(t.intersection, 1, EPS, "AABB: exit distance from inside");

    t = bs::AABB::intersect(ray({ 0, 0, 1 }, { 0, 0, -1 }), aabb);
    checkNear(t.intersection, 3, EPS, "AABB: exit distance from inside, reversed");

    t = bs::AABB::intersect(ray({ -10, -10, -10 }, { 1, 1, 1 }), aabb);
    checkNear(t.intersection, std::sqrt(3.0) * 8, EPS, "AABB: diagonal ray");

    // Axis aligned rays passing beside the box. These have zero direction
    // components and must still be tested against the corresponding slab.
    t = bs::AABB::intersect(ray({ 5, 0, -10 }, { 0, 0, 1 }), aabb);
    check(!t.valid(), "AABB: axis aligned ray beside box in x misses");
    t = bs::AABB::intersect(ray({ 0, -5, -10 }, { 0, 0, 1 }), aabb);
    check(!t.valid(), "AABB: axis aligned ray beside box in y misses");
    t = bs::AABB::intersect(ray({ 3, 3, -10 }, { 0, 0, 1 }), aabb);
    check(!t.valid(), "AABB: axis aligned ray beside box in x and y misses");
    t = bs::AABB::intersect(ray({ -10, 0, 5 }, { 1, 0, 0 }), aabb);
    check(!t.valid(), "AABB: axis aligned ray along x beside box misses");
    t = bs::AABB::intersect(ray({ -10, 5, 0 }, { 1, 1, 0 }), aabb);
    check(!t.valid(), "AABB: ray in xy-plane passing beside box misses");
    t = bs::AABB::intersect(ray({ -10, 1, 1 }, { 1, 0, 0 }), aabb);
    check(t.valid(), "AABB: axis aligned ray through box hits");

    check(bs::AABB::pointInside({ 0, 1, -1 }, aabb), "AABB: pointInside inner point");
    check(bs::AABB::pointInside({ 2, 2, 2 }, aabb), "AABB: pointInside includes boundary");
    check(!bs::AABB::pointInside({ 0, 2.1, 0 }, aabb), "AABB: pointInside outer point");

    const std::array<double, 6> b2 = { 1, 1, 1, 3, 3, 3 };
    const std::array<double, 6> b3 = { 2.5, -1, -1, 4, 1, 1 };
    check(bs::AABB::overlap(aabb, b2) && bs::AABB::collide(aabb, b2), "AABB: overlapping boxes");
    check(!bs::AABB::overlap(aabb, b3) && !bs::AABB::collide(aabb, b3), "AABB: separated boxes");
    check(bs::AABB::overlap(b2, b3) == bs::AABB::collide(b2, b3), "AABB: overlap and collide agree");

    const auto interval = bs::AABB::intersectForwardInterval(ray({ -5, 0, 0 }, { 1, 0, 0 }), aabb);
    check(interval.has_value() && std::abs(interval.value()[0] - 3) < EPS && std::abs(interval.value()[1] - 7) < EPS, "AABB: forward interval");

    const auto viz = bs::AABB::intersectVisualization<int>(ray({ -5, 0.5, 0.3 }, { 1, 0, 0 }), aabb);
    check(viz.valid() && std::abs(std::abs(viz.normal[0]) - 1) < EPS, "AABB: visualization normal is face normal");

    const auto r = rayConsistency([&](const auto& p) { return bs::AABB::intersect(p, b3); },
        [&](const auto& pos) { return bs::AABB::pointInside(pos, b3); },
        { 0, -3, -3, 6, 3, 3 }, 12, 3000, 2);
    reportConsistency(r, "AABB");
}

void testCylinder()
{
    const auto cylinder = bs::cylinder::Cylinder({ 0, 0, 0 }, { 0, 0, 1 }, 2, 2);

    auto t = bs::cylinder::intersect(ray({ 0, 0, -8000 }, { 0, 0, 1 }), cylinder);
    check(t.valid() && !t.rayOriginIsInsideItem, "Cylinder: hit end cap from far along axis");
    checkNear(t.intersection, 7998, 1E-8, "Cylinder: distance to end cap");

    t = bs::cylinder::intersect(ray({ -8000, 0, 0 }, { 1, 0, 0 }), cylinder);
    checkNear(t.intersection, 7998, 1E-8, "Cylinder: distance to wall along x");

    t = bs::cylinder::intersect(ray({ 0, 8000, 0 }, { 0, -1, 0 }), cylinder);
    checkNear(t.intersection, 7998, 1E-8, "Cylinder: distance to wall along -y");

    t = bs::cylinder::intersect(ray({ -1, 0, -3 }, { 1, 0, 1 }), cylinder);
    check(t.valid(), "Cylinder: oblique hit on end cap");
    checkNear(t.intersection, std::sqrt(2.0), 1E-12, "Cylinder: oblique distance to end cap");

    t = bs::cylinder::intersect(ray({ -1, 0, -1 }, { 1, 0, 1 }), cylinder);
    check(t.valid() && t.rayOriginIsInsideItem, "Cylinder: origin inside is flagged");
    checkNear(t.intersection, 3 * std::sqrt(2.0), 1E-12, "Cylinder: oblique exit through wall");

    t = bs::cylinder::intersect(ray({ -1, 0, -3 }, { -1, 0, -1 }), cylinder);
    check(!t.valid(), "Cylinder: ray pointing away misses");

    t = bs::cylinder::intersect(ray({ -1, -1, 0 }, { 0, 0, -1 }), cylinder);
    check(t.valid() && t.rayOriginIsInsideItem, "Cylinder: inside, parallel to axis");
    checkNear(t.intersection, cylinder.half_height, 1E-12, "Cylinder: exit through end cap parallel to axis");

    t = bs::cylinder::intersect(ray({ -1, -1, -4 }, { 0, 0, -1 }), cylinder);
    check(!t.valid(), "Cylinder: outside, parallel to axis pointing away misses");

    t = bs::cylinder::intersect(ray({ 3, 0, -10 }, { 0, 0, 1 }), cylinder);
    check(!t.valid(), "Cylinder: parallel to axis outside radius misses");

    checkNear(cylinder.volume(), std::numbers::pi * 4 * 4, 1E-12, "Cylinder: volume");

    // pointInside compared with analytical test for random points
    auto state = xraymctest::seededState();
    bool insideOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        const Vec pos = {
            state.randomUniform(-cylinder.radius, cylinder.radius),
            state.randomUniform(-cylinder.radius, cylinder.radius),
            state.randomUniform(-cylinder.half_height - 1, cylinder.half_height + 1),
        };
        bool expected = pos[0] * pos[0] + pos[1] * pos[1] <= cylinder.radius * cylinder.radius;
        expected = expected && -cylinder.half_height <= pos[2] && pos[2] <= cylinder.half_height;
        insideOk = insideOk && bs::cylinder::pointInside(pos, cylinder) == expected;
    }
    check(insideOk, "Cylinder: pointInside agrees with analytical test");

    // Rotated and translated cylinder
    const auto tilted = bs::cylinder::Cylinder({ 1, 2, 3 }, { 1, 1, 0.3 }, 1.5, 2.5);
    const auto aabb = bs::cylinder::cylinderAABB(tilted);
    bool aabbOk = true;
    for (std::size_t i = 0; i < 20000; ++i) {
        const Vec pos = { state.randomUniform(-4.0, 6.0), state.randomUniform(-3.0, 7.0), state.randomUniform(-2.0, 8.0) };
        if (bs::cylinder::pointInside(pos, tilted))
            aabbOk = aabbOk && bs::AABB::pointInside(pos, aabb);
    }
    check(aabbOk, "Cylinder: AABB contains all inside points of tilted cylinder");

    auto rc = rayConsistency([&](const auto& p) { return bs::cylinder::intersect(p, tilted); },
        [&](const auto& pos) { return bs::cylinder::pointInside(pos, tilted); },
        { -3, -2, -1, 5, 6, 7 }, 14, 3000, 3);
    reportConsistency(rc, "Cylinder (tilted)");

    // forward interval is always ordered and non-negative
    bool forwardOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        const Vec pos = { state.randomUniform(-8.0, 8.0), state.randomUniform(-8.0, 8.0), state.randomUniform(-8.0, 8.0) };
        const auto r = bs::cylinder::intersectForwardInterval(ray(pos, randomDirection(state)), cylinder);
        if (r) {
            const auto v = r.value();
            forwardOk = forwardOk && v[0] >= 0 && v[1] > 0 && v[1] >= v[0];
        }
    }
    check(forwardOk, "Cylinder: forward intervals are non-negative and ordered");

    const auto viz = bs::cylinder::intersectVisualization<int>(ray({ -5, 0, 0.5 }, { 1, 0, 0 }), cylinder);
    check(viz.valid() && std::abs(std::abs(viz.normal[0]) - 1) < 1E-6, "Cylinder: visualization normal on wall is radial");
}

void testCylinderZ()
{
    const Vec center = { 1, -1, 2 };
    const double radius = 2;
    const double half = 3;

    auto t = bs::cylinderZ::intersect(ray({ 1, -1, -20 }, { 0, 0, 1 }), center, radius, half);
    check(t.valid() && !t.rayOriginIsInsideItem, "CylinderZ: hit end cap along axis");
    checkNear(t.intersection, 19, EPS, "CylinderZ: distance to end cap");

    t = bs::cylinderZ::intersect(ray({ -20, -1, 2 }, { 1, 0, 0 }), center, radius, half);
    checkNear(t.intersection, 19, EPS, "CylinderZ: distance to wall");

    t = bs::cylinderZ::intersect(ray({ 1, -1, 2 }, { 1, 0, 0 }), center, radius, half);
    check(t.valid() && t.rayOriginIsInsideItem, "CylinderZ: origin inside is flagged");
    checkNear(t.intersection, 2, EPS, "CylinderZ: exit through wall");

    t = bs::cylinderZ::intersect(ray({ 5, -1, -20 }, { 0, 0, 1 }), center, radius, half);
    check(!t.valid(), "CylinderZ: parallel to axis outside radius misses");

    const auto r = rayConsistency([&](const auto& p) { return bs::cylinderZ::intersect(p, center, radius, half); },
        [&](const auto& pos) { return bs::cylinderZ::pointInside(pos, center, radius, half); },
        { -3, -5, -3, 5, 3, 7 }, 16, 3000, 4);
    reportConsistency(r, "CylinderZ");

    const auto fi = bs::cylinderZ::intersectForwardInterval(ray({ -20, -1, 2 }, { 1, 0, 0 }), center, radius, half);
    check(fi.has_value() && std::abs(fi.value()[0] - 19) < EPS && std::abs(fi.value()[1] - 23) < EPS, "CylinderZ: forward interval");
}

void testTetrahedron()
{
    const Vec v0 = { 0, 0, 0 };
    const Vec v1 = { 1, 0, 0 };
    const Vec v2 = { 0, 1, 0 };
    const Vec v3 = { 0, 0, 1 };

    check(bs::tetrahedron::tetrahedronNormalsOutsideTest(v0, v1, v2, v3), "Tetrahedron: reference tetrahedron has outward normals");
    check(!bs::tetrahedron::tetrahedronNormalsOutsideTest(v0, v2, v1, v3), "Tetrahedron: swapped vertices give inward normals");
    checkNear(bs::tetrahedron::volume(v0, v1, v2, v3), 1.0 / 6.0, 1E-15, "Tetrahedron: volume of unit corner tetrahedron");
    checkNear(bs::tetrahedron::volume(v0, v2, v1, v3), 1.0 / 6.0, 1E-15, "Tetrahedron: volume is independent of orientation");

    const auto aabb = bs::tetrahedron::AABB(v0, v1, v2, v3);
    check(aabb == std::array<double, 6> { 0, 0, 0, 1, 1, 1 }, "Tetrahedron: AABB");

    check(bs::tetrahedron::pointInside(v0, v1, v2, v3, { 0.1, 0.1, 0.1 }), "Tetrahedron: pointInside inner point");
    check(!bs::tetrahedron::pointInside(v0, v1, v2, v3, { 0.5, 0.5, 0.5 }), "Tetrahedron: pointInside outer point");
    check(!bs::tetrahedron::pointInside(v0, v1, v2, v3, { -0.1, 0.1, 0.1 }), "Tetrahedron: pointInside point behind face");

    auto t = bs::tetrahedron::intersect(v0, v1, v2, v3, ray({ 0.2, 0.2, -5 }, { 0, 0, 1 }));
    check(t.valid() && !t.rayOriginIsInsideItem, "Tetrahedron: hit from outside");
    checkNear(t.intersection, 5, 1E-12, "Tetrahedron: distance to face");

    t = bs::tetrahedron::intersect(v0, v1, v2, v3, ray({ 0.2, 0.2, 0.2 }, { 0, 0, 1 }));
    check(t.valid() && t.rayOriginIsInsideItem, "Tetrahedron: origin inside is flagged");
    checkNear(t.intersection, 0.4, 1E-12, "Tetrahedron: exit distance through slanted face");

    t = bs::tetrahedron::intersect(v0, v1, v2, v3, ray({ 2, 2, -5 }, { 0, 0, 1 }));
    check(!t.valid(), "Tetrahedron: ray beside tetrahedron misses");

    // Volume from Monte Carlo sampling of pointInside
    auto state = xraymctest::seededState();
    const Vec w0 = { -1, -1, 0 }, w1 = { 2, -0.5, 0.2 }, w2 = { 0.3, 2, -0.4 }, w3 = { 0.1, 0.2, 2.5 };
    const auto [a, b, c, d] = bs::tetrahedron::tetrahedronNormalsOutsideTest(w0, w1, w2, w3) ? std::array { w0, w1, w2, w3 } : std::array { w0, w2, w1, w3 };
    const auto box = bs::tetrahedron::AABB(a, b, c, d);
    const auto boxVolume = (box[3] - box[0]) * (box[4] - box[1]) * (box[5] - box[2]);
    constexpr std::size_t N = 1000000;
    std::size_t nInside = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const Vec pos = { state.randomUniform(box[0], box[3]), state.randomUniform(box[1], box[4]), state.randomUniform(box[2], box[5]) };
        if (bs::tetrahedron::pointInside(a, b, c, d, pos))
            ++nInside;
    }
    const auto vol = bs::tetrahedron::volume(a, b, c, d);
    const auto frac = vol / boxVolume;
    checkNear(boxVolume * nInside / N, vol, 5 * boxVolume * std::sqrt(frac * (1 - frac) / N), "Tetrahedron: pointInside volume agrees with analytical volume");

    const auto r = rayConsistency([&](const auto& p) { return bs::tetrahedron::intersect(a, b, c, d, p); },
        [&](const auto& pos) { return bs::tetrahedron::pointInside(a, b, c, d, pos); },
        { box[0] - 1, box[1] - 1, box[2] - 1, box[3] + 1, box[4] + 1, box[5] + 1 }, 10, 3000, 5);
    reportConsistency(r, "Tetrahedron");

    // Triangle and AABB overlap
    check(bs::tetrahedron::triangleAABBoverlap(v0, v1, v2, { -0.5, -0.5, -0.5, 0.5, 0.5, 0.5 }), "Tetrahedron: triangle overlaps box");
    check(!bs::tetrahedron::triangleAABBoverlap(v0, v1, v2, { 0.6, 0.6, -0.5, 1, 1, 0.5 }), "Tetrahedron: triangle hypotenuse misses box corner");
    check(bs::tetrahedron::insideAABB(v0, v1, v2, v3, { -1, -1, -1, 2, 2, 2 }), "Tetrahedron: insideAABB for enclosing box");
    check(!bs::tetrahedron::insideAABB(v0, v1, v2, v3, { 0.5, 0.5, 0.5, 2, 2, 2 }), "Tetrahedron: insideAABB false for box beyond slanted face");

    // Cross validate the triangle box test against the separating axis test in collision.hpp
    std::size_t disagree = 0;
    std::size_t overlaps = 0;
    constexpr std::size_t NT = 200000;
    for (std::size_t i = 0; i < NT; ++i) {
        auto rv = [&]() { return Vec { state.randomUniform(-2.0, 2.0), state.randomUniform(-2.0, 2.0), state.randomUniform(-2.0, 2.0) }; };
        const auto t0 = rv(), t1 = rv(), t2 = rv();
        const Vec bc = rv();
        const Vec bh = { state.randomUniform(0.05, 1.0), state.randomUniform(0.05, 1.0), state.randomUniform(0.05, 1.0) };
        const std::array<double, 6> box = { bc[0] - bh[0], bc[1] - bh[1], bc[2] - bh[2], bc[0] + bh[0], bc[1] + bh[1], bc[2] + bh[2] };
        const bool ref = xraymc::collision::testCollision(xraymc::Triangle(t0, t1, t2), box);
        const bool res = bs::tetrahedron::triangleAABBoverlap(t0, t1, t2, box);
        if (ref)
            ++overlaps;
        if (ref != res)
            ++disagree;
    }
    check(overlaps > NT / 10 && overlaps < NT - NT / 10, "Tetrahedron: random triangle/box sample has both outcomes");
    check(disagree == 0, "Tetrahedron: triangleAABBoverlap agrees with collision::testCollision (" + std::to_string(disagree) + " of " + std::to_string(NT) + " disagree)");
}

int main()
{
    std::cout << "Testing ray intersection on basic shapes\n";
    testSphere();
    testAABB();
    testCylinder();
    testCylinderZ();
    testTetrahedron();
    return xraymctest::summary("testbasicshapes");
}

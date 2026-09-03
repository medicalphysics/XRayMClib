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

#include "xraymc/world/collision.hpp"

#include <array>
#include <iostream>
#include <string>
#include <vector>

using Vec = std::array<double, 3>;

static int g_failures = 0;

static void check(bool value, bool expected, const std::string& name)
{
    const bool ok = value == expected;
    if (!ok)
        ++g_failures;
    std::cout << (ok ? "[ PASS ] " : "[ FAIL ] ") << name
              << " (got " << std::boolalpha << value << ", expected " << expected << ")\n";
}

static xraymc::Triangle tri(const Vec& a, const Vec& b, const Vec& c)
{
    return xraymc::Triangle(a, b, c);
}

// ---------------------------------------------------------------------------
// AABB overlap tests
// ---------------------------------------------------------------------------
static void testAABB()
{
    // layout: {min_x, min_y, min_z, max_x, max_y, max_z}
    const std::array<double, 6> unit { 0, 0, 0, 1, 1, 1 };

    check(xraymc::collision::testCollision(unit, unit), true, "AABB: identical boxes overlap");
    check(xraymc::collision::testCollision(unit, { 0.5, 0.5, 0.5, 2, 2, 2 }), true, "AABB: partial overlap");
    check(xraymc::collision::testCollision(unit, { 0.25, 0.25, 0.25, 0.75, 0.75, 0.75 }), true, "AABB: fully contained");
    check(xraymc::collision::testCollision(unit, { 1, 1, 1, 2, 2, 2 }), true, "AABB: face/corner touching counts as overlap");
    check(xraymc::collision::testCollision(unit, { 2, 0, 0, 3, 1, 1 }), false, "AABB: separated along x");
    check(xraymc::collision::testCollision(unit, { 0, 0, 5, 1, 1, 6 }), false, "AABB: separated along z");
    check(xraymc::collision::testCollision(unit, { -3, -3, -3, -1, -1, -1 }), false, "AABB: separated in all axes");
    // overlaps on two axes but not the third -> no collision
    check(xraymc::collision::testCollision(unit, { 0.5, 0.5, 2, 1.5, 1.5, 3 }), false, "AABB: overlap x/y but gap in z");
}

// ---------------------------------------------------------------------------
// Triangle-triangle tests: non-coplanar
// ---------------------------------------------------------------------------
static void testTriangleNonCoplanar()
{
    // Two parallel triangles in z=0 and z=5 planes -> all vertices on one side.
    {
        auto t1 = tri({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 });
        auto t2 = tri({ 0, 0, 5 }, { 1, 0, 5 }, { 0, 1, 5 });
        check(xraymc::collision::testCollision(t1, t2), false, "Tri: parallel, separated in z");
        check(xraymc::collision::testCollision(t2, t1), false, "Tri: parallel, separated in z (swapped)");
    }

    // A vertical triangle piercing through the interior of a horizontal one.
    {
        auto floor = tri({ -1, -1, 0 }, { 1, -1, 0 }, { 0, 1, 0 });
        auto blade = tri({ 0, 0, -1 }, { 0, 0, 1 }, { 0, 1, 0 });
        check(xraymc::collision::testCollision(floor, blade), true, "Tri: blade pierces floor interior");
        check(xraymc::collision::testCollision(blade, floor), true, "Tri: blade pierces floor interior (swapped)");
    }

    // Each triangle straddles the other's *plane*, but their intervals on the
    // plane-plane intersection line are disjoint -> NOT a collision. This guards
    // the interval-overlap step of the Moeller test.
    {
        auto near = tri({ -10, -1, 0 }, { -10, 1, 0 }, { 10, 0, 0 }); // thin, along x, in z=0
        auto far = tri({ 0, 50, -1 }, { 0, 50, 1 }, { 0, 52, 0 }); // in x=0 plane, far in +y
        check(xraymc::collision::testCollision(near, far), false, "Tri: cross planes but disjoint intervals (skew, far apart)");
        check(xraymc::collision::testCollision(far, near), false, "Tri: cross planes but disjoint intervals (swapped)");
    }

    // Triangle whose plane the other only touches at a single vertex.
    {
        auto floor = tri({ -1, -1, 0 }, { 3, -1, 0 }, { -1, 3, 0 });
        auto tip = tri({ 0.2, 0.2, 0 }, { 0.2, 0.2, 2 }, { 1, 1, 2 }); // one vertex on z=0 inside floor
        check(xraymc::collision::testCollision(floor, tip), true, "Tri: vertex touches interior of other triangle");
    }

    // Triangle above the other's plane, tilted, no contact.
    {
        auto floor = tri({ -1, -1, 0 }, { 1, -1, 0 }, { 0, 1, 0 });
        auto roof = tri({ -1, -1, 1 }, { 1, -1, 2 }, { 0, 1, 1.5 });
        check(xraymc::collision::testCollision(floor, roof), false, "Tri: tilted triangle entirely above");
    }
}

// ---------------------------------------------------------------------------
// Triangle-triangle tests: coplanar
// ---------------------------------------------------------------------------
static void testTriangleCoplanar()
{
    // Overlapping coplanar triangles in z=0.
    {
        auto t1 = tri({ 0, 0, 0 }, { 4, 0, 0 }, { 0, 4, 0 });
        auto t2 = tri({ 1, 1, 0 }, { 5, 1, 0 }, { 1, 5, 0 });
        check(xraymc::collision::testCollision(t1, t2), true, "Tri: coplanar overlapping");
        check(xraymc::collision::testCollision(t2, t1), true, "Tri: coplanar overlapping (swapped)");
    }

    // Disjoint coplanar triangles in z=0. Guards against the coplanar branch
    // always reporting a hit.
    {
        auto t1 = tri({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 });
        auto t2 = tri({ 10, 10, 0 }, { 11, 10, 0 }, { 10, 11, 0 });
        check(xraymc::collision::testCollision(t1, t2), false, "Tri: coplanar disjoint");
        check(xraymc::collision::testCollision(t2, t1), false, "Tri: coplanar disjoint (swapped)");
    }

    // One coplanar triangle fully inside another.
    {
        auto big = tri({ -10, -10, 2 }, { 10, -10, 2 }, { 0, 10, 2 });
        auto small = tri({ -1, -1, 2 }, { 1, -1, 2 }, { 0, 1, 2 });
        check(xraymc::collision::testCollision(big, small), true, "Tri: coplanar, one inside the other");
        check(xraymc::collision::testCollision(small, big), true, "Tri: coplanar, one inside the other (swapped)");
    }

    // Identical triangle.
    {
        auto t = tri({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 });
        check(xraymc::collision::testCollision(t, t), true, "Tri: identical triangles");
    }

    // Coplanar triangles sharing only an edge (touching, not overlapping).
    {
        auto t1 = tri({ 0, 0, 0 }, { 2, 0, 0 }, { 0, 2, 0 });
        auto t2 = tri({ 2, 0, 0 }, { 0, 2, 0 }, { 2, 2, 0 });
        check(xraymc::collision::testCollision(t1, t2), true, "Tri: coplanar sharing an edge");
    }

    // Parallel but offset coplanar-plane triangles (same z, no xy overlap already
    // covered; here different z so definitely no contact).
    {
        auto t1 = tri({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 });
        auto t2 = tri({ 0, 0, 3 }, { 1, 0, 3 }, { 0, 1, 3 });
        check(xraymc::collision::testCollision(t1, t2), false, "Tri: parallel planes, no contact");
    }
}

// ---------------------------------------------------------------------------
// Mesh-mesh tests
// ---------------------------------------------------------------------------
std::vector<xraymc::Triangle> getPyramid()
{
    std::vector<Vec> p;
    constexpr auto d = 30.0;
    p.push_back({ 1, 1, 0 }); // 0
    p.push_back({ 1, -1, 0 }); // 1
    p.push_back({ -1, -1, 0 }); // 2
    p.push_back({ -1, 1, 0 }); // 3
    p.push_back({ 0, 0, 1 });
    for (auto& i : p)
        for (auto& j : i)
            j *= d;

    std::vector<xraymc::Triangle> t;
    t.push_back({ p[0], p[1], p[4] });
    t.push_back({ p[1], p[2], p[4] });
    t.push_back({ p[2], p[3], p[4] });
    t.push_back({ p[3], p[0], p[4] });

    // underside
    t.push_back({ p[0], p[3], p[2] });
    t.push_back({ p[2], p[1], p[0] });

    return t;
}

std::vector<xraymc::Triangle> getBox(double scale = 1)
{
    std::vector<Vec> p;
    p.push_back({ 1, 1, 1 }); // 0
    p.push_back({ 1, 1, -1 }); // 1
    p.push_back({ 1, -1, 1 }); // 2
    p.push_back({ -1, 1, 1 }); // 3
    p.push_back({ -1, -1, 1 }); // 4
    p.push_back({ -1, 1, -1 }); // 5
    p.push_back({ 1, -1, -1 }); // 6
    p.push_back({ -1, -1, -1 }); // 7
    for (auto& i : p)
        for (auto& j : i)
            j *= scale;

    std::vector<xraymc::Triangle> t;
    t.push_back({ p[0], p[3], p[4] });
    t.push_back({ p[0], p[4], p[2] });
    t.push_back({ p[6], p[2], p[4] });
    t.push_back({ p[6], p[4], p[7] });
    t.push_back({ p[7], p[4], p[3] });
    t.push_back({ p[7], p[3], p[5] });
    t.push_back({ p[5], p[1], p[6] });
    t.push_back({ p[5], p[6], p[7] });
    t.push_back({ p[1], p[0], p[2] });
    t.push_back({ p[1], p[2], p[6] });
    t.push_back({ p[5], p[3], p[0] });
    t.push_back({ p[5], p[0], p[1] });
    return t;
}

// A unit cube [-1, 1]^3 split into 12 tetrahedra (two per face, meeting at the
// two face centers on z = +/-1). Same helper shape as tetCube() in
// testtetrahedalmesh.cpp.
xraymc::TetrahedalMeshData tetCube()
{
    std::vector<Vec> v(10);
    v[0] = { 1, -1, 1 };
    v[1] = { -1, 1, -1 };
    v[2] = { -1, -1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, 1, 1 };
    v[6] = { 1, 1, 1 };
    v[7] = { 1, 1, -1 };
    v[8] = { 0, 0, 1 };
    v[9] = { 0, 0, -1 };

    std::vector<std::array<std::uint32_t, 4>> t(12);
    t[0] = { 0, 7, 4, 9 };
    t[1] = { 1, 6, 5, 8 };
    t[2] = { 1, 6, 8, 9 };
    t[3] = { 6, 1, 7, 9 };
    t[4] = { 1, 5, 2, 8 };
    t[5] = { 3, 0, 4, 9 };
    t[6] = { 3, 1, 2, 9 };
    t[7] = { 2, 1, 8, 9 };
    t[8] = { 7, 0, 6, 9 };
    t[9] = { 6, 0, 8, 9 };
    t[10] = { 0, 3, 2, 9 };
    t[11] = { 0, 2, 8, 9 };

    xraymc::TetrahedalMeshData data;
    data.nodes = v;
    data.elements = t;
    data.collectionIndices.resize(t.size(), 0);
    data.collectionMaterialComposition.resize(1);
    data.collectionMaterialComposition[0][1] = 0.111894;
    data.collectionMaterialComposition[0][8] = 0.888106;
    data.collectionDensities.resize(1, 1.0);
    data.makeGenericCollectionNames();
    return data;
}

// ---------------------------------------------------------------------------
// Tetrahedral-mesh tests (collision uses the outer-contour triangles)
// ---------------------------------------------------------------------------
static void testTetrahedral()
{
    using TetMesh = xraymc::TetrahedalMesh<5, 1, false>;
    using TriMesh = xraymc::TriangulatedMesh<5, 1>;

    // Two identical unit cubes -> overlap everywhere.
    {
        TetMesh a(tetCube());
        TetMesh b(tetCube());
        check(xraymc::collision::testCollision(a, b), true, "Tet: identical cubes");
    }

    // Cubes offset so they interpenetrate.
    {
        TetMesh a(tetCube());
        TetMesh b(tetCube());
        b.translate({ 1.0, 0.3, 0.2 });
        check(xraymc::collision::testCollision(a, b), true, "Tet: two cubes interpenetrating");
        check(xraymc::collision::testCollision(b, a), true, "Tet: two cubes interpenetrating (swapped)");
    }

    // Cubes side by side with a clear gap between their surfaces.
    {
        TetMesh a(tetCube());
        TetMesh b(tetCube());
        b.translate({ 2.5, 0, 0 }); // half-size 1 each -> 0.5 gap
        check(xraymc::collision::testCollision(a, b), false, "Tet: two cubes with a surface gap");
    }

    // Far apart -> rejected by the AABB pre-test.
    {
        TetMesh a(tetCube());
        TetMesh b(tetCube());
        b.translate({ 1000, 0, 0 });
        check(xraymc::collision::testCollision(a, b), false, "Tet: two cubes far away, AABB reject");
    }

    // Tetrahedral mesh vs triangulated mesh, both orderings.
    {
        TetMesh tet(tetCube());
        TriMesh box(getBox());
        box.translate({ 1.0, 0.3, 0.2 });
        check(xraymc::collision::testCollision(tet, box), true, "Tet/Tri: overlapping");
        check(xraymc::collision::testCollision(box, tet), true, "Tri/Tet: overlapping (swapped)");
    }
    {
        TetMesh tet(tetCube());
        TriMesh box(getBox());
        box.translate({ 1000, 0, 0 });
        check(xraymc::collision::testCollision(tet, box), false, "Tet/Tri: far apart, AABB reject");
        check(xraymc::collision::testCollision(box, tet), false, "Tri/Tet: far apart, AABB reject (swapped)");
    }
}

static void testMesh()
{
    using Mesh = xraymc::TriangulatedMesh<5, 1>;

    // Box overlapping the flank of the pyramid.
    {
        Mesh a(getPyramid());
        Mesh b(getBox());
        b.translate({ 3, 0, -1 });
        check(xraymc::collision::testCollision(a, b), true, "Mesh: box overlaps pyramid flank");
        check(xraymc::collision::testCollision(b, a), true, "Mesh: box overlaps pyramid flank (swapped)");
    }

    // Box far away -> rejected by the AABB pre-test.
    {
        Mesh a(getPyramid());
        Mesh b(getBox());
        b.translate({ 1000, 0, 0 });
        check(xraymc::collision::testCollision(a, b), false, "Mesh: box far away, AABB reject");
        check(xraymc::collision::testCollision(b, a), false, "Mesh: box far away, AABB reject (swapped)");
    }

    // Two boxes whose AABBs overlap but whose surfaces do not touch
    // (small box sitting in the gap between the pyramid apex and base is awkward,
    // so use two boxes side by side with a clear gap).
    {
        Mesh a(getBox());
        Mesh b(getBox());
        b.translate({ 2.5, 0, 0 }); // half-size 1 each -> 0.5 gap between surfaces
        check(xraymc::collision::testCollision(a, b), false, "Mesh: two boxes with a surface gap");
    }

    // Two boxes clearly interpenetrating.
    {
        Mesh a(getBox());
        Mesh b(getBox());
        b.translate({ 1.0, 0.3, 0.2 });
        check(xraymc::collision::testCollision(a, b), true, "Mesh: two boxes interpenetrating");
    }

    // Identical meshes.
    {
        Mesh a(getBox());
        Mesh b(getBox());
        check(xraymc::collision::testCollision(a, b), true, "Mesh: identical boxes");
    }

    // A box scaled small and moved so it straddles one face of the big box.
    {
        Mesh a(getBox(4.0));
        Mesh b(getBox(1.0));
        b.translate({ 4.0, 0, 0 }); // centered on the +x face of the big box
        check(xraymc::collision::testCollision(a, b), true, "Mesh: small box straddling a face of a big box");
    }
}

int main()
{
    testAABB();
    testTriangleNonCoplanar();
    testTriangleCoplanar();
    testMesh();
    testTetrahedral();

    if (g_failures == 0) {
        std::cout << "\nAll collision tests passed.\n";
        return EXIT_SUCCESS;
    }
    std::cout << "\n"
              << g_failures << " collision test(s) failed.\n";
    return EXIT_FAILURE;
}

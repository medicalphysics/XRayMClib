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

#include "testgeometry.hpp"
#include "testutils.hpp"

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkStatisticallyEqual;
using Vec = std::array<double, 3>;
namespace vm = xraymc::vectormath;

using Mesh = TriangulatedMesh<5, 1>;

Particle randomRay(RandomState& state, double extent)
{
    const Vec pos = { state.randomUniform(-extent, extent), state.randomUniform(-extent, extent), state.randomUniform(-extent, extent) };
    Vec dir = { state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0) };
    vm::normalize(dir);
    return { .pos = pos, .dir = dir, .energy = 60, .weight = 1 };
}

void testIntersections()
{
    // A triangulated box must intersect exactly as the corresponding AABB
    Mesh mesh(xraymctest::boxTriangles(2));
    mesh.translate({ 1, -1, 0.5 });
    const std::array<double, 6> aabb = { -1, -3, -1.5, 3, 1, 2.5 };
    bool aabbOk = true;
    for (std::size_t i = 0; i < 6; ++i)
        aabbOk = aabbOk && std::abs(mesh.AABB()[i] - aabb[i]) < 1E-5; // the mesh AABB is padded slightly
    check(aabbOk, "TriangulatedMesh: AABB of translated box");
    check(xraymctest::checkVectorNear(mesh.center(), { 1, -1, 0.5 }, 1E-9, "TriangulatedMesh: center"), "TriangulatedMesh: center of translated box");

    auto state = xraymctest::seededState();
    std::size_t mismatches = 0, hits = 0;
    for (std::size_t i = 0; i < 20000; ++i) {
        const auto p = randomRay(state, 5);
        const auto a = mesh.intersect(p);
        const auto b = basicshape::AABB::intersect(p, aabb);
        if (b.valid())
            ++hits;
        if (a.valid() != b.valid() || (a.valid() && (std::abs(a.intersection - b.intersection) > 1E-9 || a.rayOriginIsInsideItem != b.rayOriginIsInsideItem)))
            ++mismatches;
    }
    check(hits > 1000, "TriangulatedMesh: random rays hit the box");
    check(mismatches == 0, "TriangulatedMesh: intersections equal AABB intersections (" + std::to_string(mismatches) + " mismatches)");

    // Intersections do not depend on the kd-tree depth
    Mesh shallow(xraymctest::pyramidTriangles(3), 1);
    Mesh deep(xraymctest::pyramidTriangles(3), 8);
    std::size_t depthMismatch = 0;
    for (std::size_t i = 0; i < 20000; ++i) {
        const auto p = randomRay(state, 6);
        const auto a = shallow.intersect(p);
        const auto b = deep.intersect(p);
        if (a.valid() != b.valid() || (a.valid() && std::abs(a.intersection - b.intersection) > 1E-9))
            ++depthMismatch;
    }
    check(depthMismatch == 0, "TriangulatedMesh: intersections independent of kd-tree depth (" + std::to_string(depthMismatch) + " mismatches)");
}

void testTransforms()
{
    Mesh mesh(xraymctest::boxTriangles(1));
    mesh.scale(3);
    checkNear(mesh.AABB()[3], 3, 1E-5, "TriangulatedMesh: scale");
    Mesh tall(xraymctest::pyramidTriangles(1));
    tall.rotate(std::numbers::pi / 2, { 1, 0, 0 });
    // the apex at z = 1 is rotated to y = -1
    checkNear(tall.AABB()[1], -1, 1E-6, "TriangulatedMesh: rotate moves apex");
    checkNear(tall.AABB()[5], 1, 1E-6, "TriangulatedMesh: rotate moves base");
    Mesh mir(xraymctest::pyramidTriangles(1));
    mir.mirror(0.0, 2);
    checkNear(mir.AABB()[2], -1, 1E-6, "TriangulatedMesh: mirror in z plane");
    // a mirrored mesh must still be a closed surface with outward normals
    Particle p { .pos = { 0, 0, -5 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto r = mir.intersect(p);
    check(r.valid() && !r.rayOriginIsInsideItem && std::abs(r.intersection - 4) < 1E-6, "TriangulatedMesh: mirrored mesh intersects from outside");
    p.pos = { 0, 0, -0.5 };
    check(mir.intersect(p).rayOriginIsInsideItem, "TriangulatedMesh: mirrored mesh has outward normals");
}

void writeAsciiSTL(const std::filesystem::path& path, const std::vector<Triangle>& tris, bool crlf)
{
    std::ofstream f(path, std::ios::binary);
    const std::string nl = crlf ? "\r\n" : "\n";
    f << "solid test" << nl;
    for (const auto& t : tris) {
        const auto n = t.planeVector();
        f << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << nl;
        f << "    outer loop" << nl;
        for (const auto& v : t.vertices())
            f << "      vertex " << v[0] << " " << v[1] << " " << v[2] << nl;
        f << "    endloop" << nl;
        f << "  endfacet" << nl;
    }
    f << "endsolid test" << nl;
}

void writeBinarySTL(const std::filesystem::path& path, const std::vector<Triangle>& tris, const std::string& header)
{
    std::ofstream f(path, std::ios::binary);
    std::array<char, 80> h {};
    std::memcpy(h.data(), header.data(), std::min(header.size(), h.size()));
    f.write(h.data(), h.size());
    const std::uint32_t n = static_cast<std::uint32_t>(tris.size());
    f.write(reinterpret_cast<const char*>(&n), 4);
    for (const auto& t : tris) {
        std::array<float, 12> data;
        const auto nv = t.planeVector();
        for (std::size_t i = 0; i < 3; ++i)
            data[i] = static_cast<float>(nv[i]);
        std::size_t k = 3;
        for (const auto& v : t.vertices())
            for (auto c : v)
                data[k++] = static_cast<float>(c);
        f.write(reinterpret_cast<const char*>(data.data()), sizeof(float) * 12);
        const std::uint16_t attr = 0;
        f.write(reinterpret_cast<const char*>(&attr), 2);
    }
}

bool sameTriangles(const std::vector<Triangle>& a, const std::vector<Triangle>& b)
{
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                if (std::abs(a[i].vertices()[j][k] - b[i].vertices()[j][k]) > 1E-5)
                    return false;
    return true;
}

void testSTL()
{
    const auto tris = xraymctest::pyramidTriangles(2.5);
    {
        const auto path = xraymctest::tempPath("ascii.stl");
        writeAsciiSTL(path, tris, false);
        STLReader reader;
        check(sameTriangles(reader(path.string()), tris), "STLReader: ASCII file");
        writeAsciiSTL(path, tris, true);
        check(sameTriangles(reader(path.string()), tris), "STLReader: ASCII file with CRLF line endings");
        Mesh m(path.string());
        check(m.triangles().size() == tris.size(), "TriangulatedMesh: constructed from STL file");
        std::filesystem::remove(path);
    }
    {
        const auto path = xraymctest::tempPath("binary.stl");
        writeBinarySTL(path, tris, "binary stl header");
        STLReader reader;
        check(sameTriangles(reader(path.string()), tris), "STLReader: binary file");
        // Many exporters write "solid" in the header of binary files
        writeBinarySTL(path, tris, "solid exported by a CAD program");
        check(sameTriangles(reader(path.string()), tris), "STLReader: binary file with header starting with solid");
        std::filesystem::remove(path);
    }
    STLReader reader;
    check(reader(xraymctest::tempPath("missing.stl").string()).empty(), "STLReader: missing file gives no triangles");
}

void testOpenSurface()
{
    // A thin lead surface absorbs the same energy as a lead box of the same thickness
    using Surface = TriangulatedOpenSurface<5, 1>;
    using Box = WorldBox<5, 1>;
    const auto lead = Material<5>::byZ(82).value();
    const double t = 0.01;
    PencilBeam<> beam({ 0.1, 0.2, -10 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(20);
    beam.setNumberOfParticlesPerExposure(20000);
    Transport transport;

    World<Surface> ws;
    auto& surface = ws.addItem<Surface>({ xraymctest::planeTriangles(1), t });
    surface.setMaterial(lead, 11.35);
    ws.setMaterialDensity(1E-9);
    ws.build();
    transport(ws, beam, nullptr, false);

    World<Box> wb;
    auto& box = wb.addItem<Box>({ { -1, -1, -t / 2, 1, 1, t / 2 }, lead, 11.35 });
    wb.setMaterialDensity(1E-9);
    wb.build();
    transport(wb, beam, nullptr, false);
    const auto& es = surface.energyScored();
    const auto& eb = box.energyScored();
    checkStatisticallyEqual(es.energyImparted(), es.standardDeviation(), eb.energyImparted(), eb.standardDeviation(), 4, "TriangulatedOpenSurface: energy equals thin box of same thickness");

    // oblique incidence through the surface also equals the box
    PencilBeam<> oblique({ -5, 0.2, -5 }, { 1, 0, 1 }, 60);
    oblique.setNumberOfExposures(20);
    oblique.setNumberOfParticlesPerExposure(20000);
    transport(ws, oblique, nullptr, false);
    transport(wb, oblique, nullptr, false);
    checkStatisticallyEqual(es.energyImparted(), es.standardDeviation(), eb.energyImparted(), eb.standardDeviation(), 4, "TriangulatedOpenSurface: energy equals thin box at oblique incidence");
}

int main()
{
    std::cout << "Testing triangulated mesh\n";
    testIntersections();
    testTransforms();
    testSTL();
    testOpenSurface();
    return xraymctest::summary("testtriangulatedmesh");
}

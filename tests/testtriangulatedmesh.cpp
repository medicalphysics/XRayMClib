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
#include <sstream>
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
    check(Mesh::surfaceValidationError(mir.triangles()).empty(), "TriangulatedMesh: mirrored mesh passes surface check");
    Mesh pointMirror(xraymctest::pyramidTriangles(1));
    pointMirror.mirror({ 0, 0, 0 });
    check(Mesh::surfaceValidationError(pointMirror.triangles()).empty(), "TriangulatedMesh: mesh mirrored through a point passes surface check");
    Mesh negative(xraymctest::pyramidTriangles(1));
    negative.scale(-2);
    check(Mesh::surfaceValidationError(negative.triangles()).empty(), "TriangulatedMesh: negatively scaled mesh passes surface check");
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

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream f(path, std::ios::binary);
    f << text;
}

std::string asciiFacet(const Triangle& t, const std::string& indent = "  ")
{
    std::ostringstream s;
    s << indent << "facet normal 0 0 0\n"
      << indent << "outer loop\n";
    for (const auto& v : t.vertices())
        s << indent << indent << "vertex " << v[0] << " " << v[1] << " " << v[2] << "\n";
    s << indent << "endloop\n"
      << indent << "endfacet\n";
    return s.str();
}

void testSTLMissingElements()
{
    const auto tris = xraymctest::boxTriangles(2);
    const auto path = xraymctest::tempPath("broken.stl");
    {
        // tab indented files must not lose vertices
        std::string text = "solid tabs\n";
        for (const auto& t : tris)
            text += asciiFacet(t, "\t");
        text += "endsolid tabs\n";
        writeText(path, text);
        STLReader reader;
        check(sameTriangles(reader(path.string()), tris), "STLReader: ASCII file with tab indentation");
    }
    {
        // a facet missing a vertex is skipped without shifting the following triangles
        std::string text = "solid missing\n";
        for (std::size_t i = 0; i < tris.size(); ++i) {
            auto facet = asciiFacet(tris[i]);
            if (i == 3)
                facet.erase(facet.find("    vertex"), facet.find("\n", facet.find("    vertex")) + 1 - facet.find("    vertex"));
            text += facet;
        }
        text += "endsolid missing\n";
        writeText(path, text);
        STLReader reader;
        const auto read = reader(path.string());
        auto expected = tris;
        expected.erase(expected.begin() + 3);
        check(sameTriangles(read, expected), "STLReader: facet with missing vertex is skipped, other triangles intact");
        check(reader.numberOfDroppedTriangles() == 1 && !reader.message().empty(), "STLReader: skipped facet is reported (" + reader.message() + ")");
        Mesh mesh(path.string());
        check(!mesh.valid() && mesh.validationError().find("not closed") != std::string::npos, "TriangulatedMesh: mesh with missing triangle is not closed (" + mesh.validationError() + ")");
        check(mesh.validationError().find("Skipped 1 facet") != std::string::npos, "TriangulatedMesh: validation error includes STL reader message");
    }
    {
        // vertex count not a multiple of three at end of file
        std::string text = "solid truncated\n" + asciiFacet(tris[0]) + "facet normal 0 0 1\nouter loop\nvertex 1 2 3\n";
        writeText(path, text);
        STLReader reader;
        const auto read = reader(path.string());
        check(read.size() == 1 && reader.numberOfDroppedTriangles() == 1, "STLReader: truncated ASCII file");
    }
    {
        // degenerate triangle is dropped
        std::string text = "solid degenerate\n" + asciiFacet(tris[0]) + asciiFacet(Triangle({ 0, 0, 0 }, { 1, 1, 1 }, { 2, 2, 2 })) + "endsolid degenerate\n";
        writeText(path, text);
        STLReader reader;
        const auto read = reader(path.string());
        check(read.size() == 1 && reader.numberOfDroppedTriangles() == 1, "STLReader: degenerate triangle is skipped");
    }
    {
        // binary file shorter than its triangle count
        writeBinarySTL(path, tris, "binary");
        std::filesystem::resize_file(path, std::filesystem::file_size(path) - 25);
        STLReader reader;
        check(reader(path.string()).empty() && reader.message().find("truncated") != std::string::npos, "STLReader: truncated binary file is rejected (" + reader.message() + ")");
        // binary file with padding after the last triangle
        writeBinarySTL(path, tris, "binary");
        {
            std::ofstream f(path, std::ios::binary | std::ios::app);
            f << "padding";
        }
        check(sameTriangles(reader(path.string()), tris), "STLReader: binary file with trailing bytes");
    }
    std::filesystem::remove(path);
}

void testSurfaceValidation()
{
    const auto tris = xraymctest::boxTriangles(1);
    check(Mesh(tris).valid(), "TriangulatedMesh: closed box is valid (" + Mesh(tris).validationError() + ")");
    check(Mesh(xraymctest::pyramidTriangles(1)).valid(), "TriangulatedMesh: closed pyramid is valid");
    check(Mesh::surfaceValidationError(xraymctest::pyramidTriangles(1)).empty(), "TriangulatedMesh: pyramid fixture has outward normals");
    check(!Mesh().valid(), "TriangulatedMesh: empty mesh is invalid");
    checkNear(vm::length(Mesh().center()), 0, 1E-12, "TriangulatedMesh: center of empty mesh");
    {
        auto t = tris;
        t.pop_back();
        check(Mesh(t).validationError().find("not closed") != std::string::npos, "TriangulatedMesh: missing triangle is detected");
    }
    {
        auto t = tris;
        const auto& v = t[4].vertices();
        t[4] = Triangle(v[0], v[2], v[1]);
        check(Mesh(t).validationError().find("consistently oriented") != std::string::npos, "TriangulatedMesh: flipped triangle is detected");
    }
    {
        auto t = tris;
        for (auto& tri : t) {
            const auto& v = tri.vertices();
            tri = Triangle(v[0], v[2], v[1]);
        }
        check(Mesh::surfaceValidationError(t).find("inwards") != std::string::npos, "TriangulatedMesh: inward normals are detected");
        // the mesh flips inward normals automatically
        Mesh inward(t);
        check(inward.valid() && Mesh::surfaceValidationError(inward.triangles()).empty(), "TriangulatedMesh: inward normals are flipped on construction");
        const Particle p { .pos = { 0, 0, 0.5 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
        check(inward.intersect(p).rayOriginIsInsideItem, "TriangulatedMesh: flipped mesh tracks particles inside");
        Mesh viaSetData;
        viaSetData.setData(t);
        check(viaSetData.valid() && Mesh::surfaceValidationError(viaSetData.triangles()).empty(), "TriangulatedMesh: setData flips inward normals");
    }
    {
        auto t = tris;
        t.push_back(t[0]);
        check(Mesh(t).validationError().find("not manifold") != std::string::npos, "TriangulatedMesh: duplicated triangle is detected");
    }
    check(Mesh(xraymctest::planeTriangles(1)).validationError().find("not closed") != std::string::npos, "TriangulatedMesh: open surface is not closed");
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
    testSTLMissingElements();
    testSurfaceValidation();
    testOpenSurface();
    return xraymctest::summary("testtriangulatedmesh");
}

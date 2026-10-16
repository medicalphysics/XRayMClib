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

Copyright 2023 Erlend Andersen
*/

#include "testgeometry.hpp"
#include "testutils.hpp"

#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using Vec = std::array<double, 3>;
namespace vm = xraymc::vectormath;

using Mesh = TetrahedalMesh<5, 1, false>;

void testStructure()
{
    const auto data = xraymctest::tetCube(1);
    Mesh mesh(data);
    check(mesh.numberOfThetrahedrons() == 12, "TetrahedalMesh: number of tetrahedra");
    double volume = 0;
    for (std::uint32_t i = 0; i < mesh.numberOfThetrahedrons(); ++i)
        volume += mesh.tetrahedalVolume(i);
    checkNear(volume, 8, 1E-12, "TetrahedalMesh: sum of tetrahedron volumes equals cube volume");
    const auto coll = mesh.collectionData();
    check(coll.size() == 1, "TetrahedalMesh: one collection");
    checkNear(coll[0].volume, 8, 1E-12, "TetrahedalMesh: collection volume");
    checkNear(coll[0].density, 1, 1E-12, "TetrahedalMesh: collection density");
    check(mesh.collectionNames().size() == 1 && mesh.tetrahedronCollectionIndex(mesh.collectionNames()[0]).has_value(), "TetrahedalMesh: collection names");
    check(!mesh.tetrahedronCollectionIndex("no such collection").has_value(), "TetrahedalMesh: unknown collection name");
    const auto aabb = mesh.AABB();
    bool aabbOk = true;
    for (std::size_t i = 0; i < 3; ++i)
        aabbOk = aabbOk && std::abs(aabb[i] + 1) < 1E-5 && std::abs(aabb[i + 3] - 1) < 1E-5;
    check(aabbOk, "TetrahedalMesh: AABB of cube");
    check(!mesh.outerContourTetrahedronIndices().empty(), "TetrahedalMesh: outer contour tetrahedra");
    double area = 0;
    for (const auto& t : mesh.constructOuterContourTriangles())
        area += 0.5 * vm::length(vm::cross(vm::subtract(t[1], t[0]), vm::subtract(t[2], t[0])));
    checkNear(area, 24, 1E-9, "TetrahedalMesh: outer contour triangles cover the cube surface");
    const auto copy = mesh.copyData();
    check(copy.nodes == data.nodes && copy.elements.size() == data.elements.size(), "TetrahedalMesh: copyData returns mesh data");
}

void testIntersections()
{
    Mesh mesh(xraymctest::tetCube(1));
    const std::array<double, 6> aabb = { -1, -1, -1, 1, 1, 1 };
    auto state = xraymctest::seededState();
    std::size_t mismatches = 0, hits = 0;
    for (std::size_t i = 0; i < 20000; ++i) {
        const Vec pos = { state.randomUniform(-3.0, 3.0), state.randomUniform(-3.0, 3.0), state.randomUniform(-3.0, 3.0) };
        Vec dir = { state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0), state.randomUniform(-1.0, 1.0) };
        vm::normalize(dir);
        const Particle p { .pos = pos, .dir = dir, .energy = 60, .weight = 1 };
        const auto a = mesh.intersect(p);
        const auto b = basicshape::AABB::intersect(p, aabb);
        if (b.valid())
            ++hits;
        if (a.valid() != b.valid() || (a.valid() && (std::abs(a.intersection - b.intersection) > 1E-6 || a.rayOriginIsInsideItem != b.rayOriginIsInsideItem)))
            ++mismatches;
    }
    check(hits > 1000, "TetrahedalMesh: random rays hit the cube");
    check(mismatches == 0, "TetrahedalMesh: intersections equal AABB intersections (" + std::to_string(mismatches) + " mismatches)");

    mesh.translate({ 10, 0, 0 });
    const Particle p { .pos = { 10, 0, -5 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    check(mesh.intersect(p).valid() && std::abs(mesh.intersect(p).intersection - 4) < 1E-6, "TetrahedalMesh: intersection after translation");
    mesh.rotate({ 0, 0, 1 }, std::numbers::pi / 4);
    // rotated about its center, the cube diagonal now spans sqrt(2) in x
    checkNear(mesh.AABB()[3] - mesh.AABB()[0], 2 * std::sqrt(2.0), 1E-5, "TetrahedalMesh: rotate updates AABB");
}

void testDoseInTetrahedra()
{
    // A pencil beam through the cube only deposits energy in tetrahedra along the beam
    World<Mesh> world;
    auto& mesh = world.addItem<Mesh>({ xraymctest::tetCube(1) });
    world.setMaterialDensity(1E-9);
    world.build();
    PencilBeam<> beam({ 0.5, -0.4, -5 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(10);
    beam.setNumberOfParticlesPerExposure(20000);
    Transport transport;
    transport(world, beam, nullptr, false);
    double total = 0, onPath = 0;
    for (std::uint32_t i = 0; i < mesh.numberOfThetrahedrons(); ++i) {
        const auto e = mesh.energyScored(i).energyImparted();
        total += e;
        const auto t = mesh.tetrahedron(i);
        // tetrahedra that contain a point on the beam line
        bool contains = false;
        for (double z = -0.99; z < 1; z += 0.01)
            contains = contains || basicshape::tetrahedron::pointInside(t[0], t[1], t[2], t[3], { 0.5, -0.4, z });
        if (contains)
            onPath += e;
    }
    check(total > 0, "TetrahedalMesh: energy deposited in mesh");
    check(onPath / total > 0.8, "TetrahedalMesh: most energy deposited in tetrahedra along the beam (" + std::to_string(onPath / total) + ")");
    mesh.addEnergyScoredToDoseScore(1.0);
    const auto coll = mesh.collectionData();
    checkNear(coll[0].dose, total / 8.0, 1E-9 * total, "TetrahedalMesh: collection dose is energy per collection mass");
}

void writeTetgenFiles(const std::filesystem::path& node, const std::filesystem::path& ele, const TetrahedalMeshData& data, bool trailingNewline)
{
    {
        std::ofstream f(node);
        f << data.nodes.size() << " 3 0 0\n";
        for (std::size_t i = 0; i < data.nodes.size(); ++i) {
            f << i << " " << data.nodes[i][0] << " " << data.nodes[i][1] << " " << data.nodes[i][2];
            if (trailingNewline || i + 1 < data.nodes.size())
                f << "\n";
        }
    }
    {
        std::ofstream f(ele);
        f << data.elements.size() << " 4 1\n";
        for (std::size_t i = 0; i < data.elements.size(); ++i) {
            const auto& e = data.elements[i];
            f << i << " " << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << " " << (i < 6 ? 0 : 1);
            if (trailingNewline || i + 1 < data.elements.size())
                f << "\n";
        }
    }
}

void testReader()
{
    const auto data = xraymctest::tetCube(1);
    const auto node = xraymctest::tempPath("cube.node");
    const auto ele = xraymctest::tempPath("cube.ele");
    for (bool newline : { true, false }) {
        writeTetgenFiles(node, ele, data, newline);
        TetrahedalMeshReader reader(node.string(), ele.string());
        const std::string label = newline ? "" : " (no trailing newline)";
        check(reader.valid(), "TetrahedalMeshReader: valid" + label);
        check(reader.data().nodes == data.nodes, "TetrahedalMeshReader: nodes read" + label);
        check(reader.data().elements == data.elements, "TetrahedalMeshReader: elements read" + label);
        // region attributes become collections
        check(reader.data().collectionDensities.size() == 2, "TetrahedalMeshReader: region attributes become collections" + label);
        if (reader.valid() && reader.data().elements.size() == 12) {
            Mesh mesh(reader.data());
            const auto coll = mesh.collectionData();
            check(coll.size() == 2 && std::abs(coll[0].volume + coll[1].volume - 8) < 1E-9, "TetrahedalMeshReader: mesh from file has cube volume" + label);
        }
    }
    std::filesystem::remove(node);
    std::filesystem::remove(ele);
    TetrahedalMeshReader missing(xraymctest::tempPath("missing.node").string(), xraymctest::tempPath("missing.ele").string());
    check(!missing.valid(), "TetrahedalMeshReader: missing files are not valid");
}

int main()
{
    std::cout << "Testing tetrahedal mesh\n";
    testStructure();
    testIntersections();
    testDoseInTetrahedra();
    testReader();
    return xraymctest::summary("testtetrahedalmesh");
}

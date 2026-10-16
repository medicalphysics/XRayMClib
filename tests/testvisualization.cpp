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

#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <filesystem>
#include <iostream>
#include <string>

using namespace xraymc;
using xraymctest::check;

using Sphere = WorldSphere<5, 1>;

template <typename U>
std::array<U, 3> pixel(const VisualizationBuffer<U>& b, std::size_t x, std::size_t y)
{
    const auto i = (x + y * b.width) * 4;
    return { b.buffer[i], b.buffer[i + 1], b.buffer[i + 2] };
}

void testRender()
{
    World<Sphere> world;
    auto& sphere = world.addItem<Sphere>({ 5, { 0, 0, 0 } }, "Sphere");
    world.build();
    VisualizeWorld viz(world);
    viz.setDistance(100);
    viz.setPolarAngleDeg(60);
    viz.setAzimuthalAngleDeg(30);
    viz.suggestFOV(1);
    viz.setBackgroundColor({ 10, 20, 30 });
    viz.setColorOfItem<std::uint8_t>(world.getItemPointerFromName("Sphere"), { 255, 0, 0 });
    auto buffer = viz.createBuffer<std::uint8_t>(65, 65);
    viz.generate(world, buffer);
    check(buffer.buffer.size() == 65 * 65 * 4, "Visualization: buffer size");
    const auto corner = pixel(buffer, 0, 0);
    check(corner == std::array<std::uint8_t, 3> { 10, 20, 30 }, "Visualization: background color in corner");
    const auto center = pixel(buffer, 32, 32);
    check(center[0] > 50 && center[1] == 0 && center[2] == 0, "Visualization: item color at image center");
    bool opaque = true;
    for (std::size_t i = 3; i < buffer.buffer.size(); i += 4)
        opaque = opaque && buffer.buffer[i] == 255;
    check(opaque, "Visualization: alpha channel is opaque");
    // the sphere is round: the center row has item pixels in a symmetric span
    std::size_t first = 65, last = 0;
    for (std::size_t x = 0; x < 65; ++x)
        if (pixel(buffer, x, 32) != std::array<std::uint8_t, 3> { 10, 20, 30 }) {
            first = std::min(first, x);
            last = std::max(last, x);
        }
    check(first < 32 && last > 32 && std::abs(static_cast<int>(first + last) - 64) <= 2, "Visualization: sphere is centered in image");

    // double buffers hold values in [0, 1]
    auto dbuffer = viz.createBuffer<double>(32, 24);
    viz.generate(world, dbuffer);
    bool range = dbuffer.buffer.size() == 32 * 24 * 4;
    for (auto v : dbuffer.buffer)
        range = range && v >= 0 && v <= 1;
    check(range, "Visualization: double buffer values in [0, 1]");

    // a line prop in front of the sphere is drawn
    viz.addLineProp({ -200, 0, 0 }, { 1, 0, 0 }, 400, 0.5);
    viz.setPolarAngleDeg(90);
    viz.setAzimuthalAngleDeg(90);
    viz.setPropColor({ 0, 255, 0 });
    viz.generate(world, buffer);
    std::size_t green = 0;
    for (std::size_t x = 0; x < 65; ++x)
        for (std::size_t y = 0; y < 65; ++y) {
            const auto c = pixel(buffer, x, y);
            if (c[1] > 0 && c[0] == 0)
                ++green;
        }
    check(green > 0, "Visualization: line prop is rendered");
    viz.clearLineProps();

    // color by value adds a color bar
    viz.addColorByValueItem(world.getItemPointerFromName("Sphere"));
    viz.setColorByValueMinMax(0, 1);
    viz.generate(world, buffer);
    check(pixel(buffer, 32, 32) != center, "Visualization: color by value changes item color");

    const auto path = xraymctest::tempPath("viz.png");
    const bool saved = VisualizeWorld<Sphere>::savePNG(path.string(), buffer);
#ifdef XRAYMCLIB_USE_LOADPNG
    check(saved && std::filesystem::exists(path) && std::filesystem::file_size(path) > 0, "Visualization: PNG written");
#else
    check(!saved, "Visualization: PNG not written without lodepng");
#endif
    std::filesystem::remove(path);
    (void)sphere;
}

// All world item types can be rendered
void testAllItems()
{
    using Box = WorldBox<5, 1>;
    using Cyl = WorldCylinder<5, 1>;
    using CTDI = CTDIPhantom<5, 1>;
    using Mesh = TriangulatedMesh<5, 1>;
    using Tet = TetrahedalMesh<5, 1, false>;
    using Grid = AAVoxelGrid<5, 1, 0>;
    World<Sphere, Box, Cyl, CTDI, Mesh, Tet, Grid> world;
    world.reserveNumberOfItems(7);
    world.addItem<Sphere>({ 2, { -20, 0, 0 } });
    world.addItem<Box>({ 2, { -10, 0, 0 } });
    world.addItem<Cyl>({ 2, 4, { 0, 0, 0 } });
    world.addItem<CTDI>({ 8, 15, { 20, 0, 0 } });
    auto& mesh = world.addItem<Mesh>({ xraymctest::boxTriangles(2) });
    mesh.translate({ 0, 10, 0 });
    auto& tet = world.addItem<Tet>({ xraymctest::tetCube(2) });
    tet.translate({ 0, -10, 0 });
    std::vector<double> dens(27, 1.0);
    std::vector<std::uint8_t> mat(27, 1);
    mat[13] = 0; // a transparent voxel in the center
    auto& grid = world.addItem<Grid>({ { 3, 3, 3 }, { 1, 1, 1 }, dens, mat, { Material<5>::byZ(1).value(), Material<5>::byNistName("Water, Liquid").value() } });
    grid.translate({ 0, 0, 10 });
    world.build();
    VisualizeWorld viz(world);
    viz.setDistance(150);
    viz.setPolarAngleDeg(70);
    viz.setAzimuthalAngleDeg(20);
    viz.suggestFOV(1);
    auto buffer = viz.createBuffer<std::uint8_t>(128, 128);
    viz.generate(world, buffer);
    std::size_t nonBackground = 0;
    const auto bg = pixel(buffer, 0, 0);
    for (std::size_t x = 0; x < 128; ++x)
        for (std::size_t y = 0; y < 128; ++y)
            if (pixel(buffer, x, y) != bg)
                ++nonBackground;
    check(nonBackground > 500, "Visualization: all item types are rendered (" + std::to_string(nonBackground) + " item pixels)");
}

int main()
{
    std::cout << "Testing visualization\n";
    testRender();
    testAllItems();
    return xraymctest::summary("testvisualization");
}

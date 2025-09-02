

#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

std::vector<xraymc::Triangle> getPyramid()
{
    std::vector<std::array<double, 3>> p;
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
    std::vector<std::array<double, 3>> p;
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
std::vector<xraymc::Triangle> getPlane(double scale = 1)
{
    std::vector<std::array<double, 3>> p;
    p.push_back({ -1, -1, 0 }); // 0
    p.push_back({ 1, -1, 0 }); // 1
    p.push_back({ 1, 1, 0 }); // 2
    p.push_back({ -1, -1, 0 }); // 3
    p.push_back({ 1, 1, 0 }); // 4
    p.push_back({ -1, 1, 0 }); // 5

    for (auto& i : p)
        for (auto& j : i)
            j *= scale;

    std::vector<xraymc::Triangle> t;
    t.push_back({ p[0], p[1], p[2] });
    t.push_back({ p[3], p[4], p[5] });
    return t;
}

template <std::size_t N = 5, int L = 2>
void testMeshVisualization()
{

    using Mesh = xraymc::TriangulatedMesh<N, L>;
    using World = xraymc::World<Mesh>;
    using Material = xraymc::Material<N>;

    World world;

    const auto waterComp = xraymc::NISTMaterials::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    int option;
    option = 0; // Box
    // option = 1; // Triangle
    // option = 2; // Bunny
    // option = 3; // bunny_low
    // option = 4; // duck

    if (option == 0) {
        const auto triangles = getBox();
        Mesh mesh(triangles);
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 1) {
        const auto triangles = getPyramid();
        Mesh mesh(triangles);
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 2) {
        Mesh mesh("bunny.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 3) {
        Mesh mesh("bunny_low.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 4) {
        Mesh mesh("duck.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    }

    world.build();

    xraymc::VisualizeWorld viz(world);

    auto buffer = viz.createBuffer();

    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(500);
        viz.setPolarAngle(std::numbers::pi_v<double> / 3.0);
        viz.setAzimuthalAngle(std::numbers::pi_v<double> * i / 6.0);
        // viz.setCameraPosition({ -60, -30, -10 });
        viz.suggestFOV();
        viz.generate(world, buffer);
        std::string name = "color_" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer);
        // writeImage(buffer, name);
    }
}

double testScoring()
{
    using Mesh = xraymc::TriangulatedMesh<5, 2>;
    using Box = xraymc::WorldBox<5, 2>;
    using World = xraymc::World<Mesh, Box>;
    using Material = xraymc::Material<5>;

    using Beam = xraymc::IsotropicMonoEnergyBeam<>;
    Beam beam({ 0, 0, -1000 }, { 1, 0, 0, 0, 1, 0 }, 60);
    beam.setNumberOfExposures(12);
    beam.setNumberOfParticlesPerExposure(1E5);

    const auto waterComp = xraymc::NISTMaterials::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    {
        World world;
        auto tri = getBox();
        auto& item = world.addItem<Mesh>({ tri });
        item.setMaterial(water, 1);
        world.build();
        xraymc::Transport transport;
        auto milli = transport.runConsole(world, beam, 1);
        std::cout << "mesh: " << item.doseScored().dose() << " time: " << milli << std::endl;
    }
    {
        World world;
        auto& item = world.addItem<Box>({ 1 });
        item.setMaterial(water, 1);
        world.build();
        xraymc::Transport transport;
        auto milli = transport.runConsole(world, beam, 1);
        std::cout << "box: " << item.doseScored().dose() << " time: " << milli << std::endl;
    }

    return 0;
}

bool testOpenSurface()
{
    /*std::vector<xraymc::Triangle> tris;
    tris.push_back({ { -0.3333333, -1, 0.5743564 },
        { 0.33333337, -1, 0 },
        { 0.33333337, 0, -0.5099079 } });
    tris.push_back({ { -0.3333333, -1, 0.5743564 },
        { 0.33333337, 0, -0.5099079 },
        { -0.3333333, 0, 0 } });

    tris.push_back({ { 0.33333337, -1, 0 },
        { 1, -1, 0.65487117 },
        { 1, 0, 0 } });
    tris.push_back({ { 0.33333337, -1, 0 },
        { 1, 0, 0 },
        { 0.33333337, 0, -0.5099079 } });
    tris.push_back({ { -0.3333333, 0, 0 },
        { 0.33333337, 0, -0.5099079 },
        { 0.33333337, 1, 0 } });
    tris.push_back({ { -0.3333333, 0, 0 },
        { 0.33333337, 1, 0 },
        { -0.3333333, 1, 0.6040865 } });
    tris.push_back({ { 0.33333337, 0, -0.5099079 },
        { 1, 0, 0 },
        { 1, 1, 0.63710684 } });
    tris.push_back({ { 0.33333337, 0, -0.5099079 },
        { 1, 1, 0.63710684 },
        { 0.33333337, 1, 0 } });
*/

    constexpr int NSHELL = 12;
    constexpr int MODEL = 2;

    using Sphere = xraymc::WorldSphere<NSHELL, MODEL>;

    std::cout << "Test surface mesh\n";

    xraymc::World<xraymc::TriangulatedOpenSurface<NSHELL, MODEL>, Sphere> w1;
    xraymc::World<xraymc::WorldBox<NSHELL, MODEL>, Sphere> w2;

    w1.reserveNumberOfItems(2);
    w2.reserveNumberOfItems(2);

    auto lead = xraymc::Material<NSHELL>::byZ(82).value();
    double lead_dens = 11.34;

    auto& mesh = w1.template addItem<xraymc::TriangulatedOpenSurface<NSHELL, MODEL>>(getPlane());
    auto& box = w2.template addItem<xraymc::WorldBox<NSHELL, MODEL>>();

    auto& sphere_mesh = w1.template addItem<Sphere>({ 1, { 2, 2, -2 } });
    auto& sphere_box = w2.template addItem<Sphere>({ 1, { 2, 2, -2 } });

    mesh.setMaterial(lead, lead_dens);
    mesh.setSurfaceThickness(0.05);

    box.setAABB({ -1, -1, 0, 1, 1, 0.05 });
    box.setMaterial(lead, lead_dens);

    w1.build();
    w2.build();

    xraymc::PencilBeam<false> beam;
    beam.setPosition({ 0, 0, -10 });
    beam.setDirection({ 0, 0, 1 });
    beam.setNumberOfExposures(160);
    beam.setNumberOfParticlesPerExposure(1000000);

    xraymc::Transport transport;
    transport.runConsole(w1, beam, 1);
    transport.runConsole(w2, beam, 1);

    std::cout << "Mesh: " << mesh.doseScored(0).dose() << ", Scatter: " << sphere_mesh.doseScored(0).dose() << ", Events: " << sphere_mesh.doseScored(0).numberOfEvents() << std::endl;
    std::cout << "Box: " << box.doseScored(0).dose() << ", Scatter: " << sphere_box.doseScored(0).dose() << ", Events: " << sphere_box.doseScored(0).numberOfEvents() << std::endl;

    if (std::abs(mesh.doseScored(0).dose() - box.doseScored(0).dose()) < 0.01) {
        std::cout << "SUCCESS\n";
        return true;
    } else {
        std::cout << "FAILURE\n";
        return false;
    }

    /*xraymc::VisualizeWorld viz(w1);
     auto buffer = viz.createBuffer();

     for (std::size_t i = 0; i < 180; i += 30) {
         viz.setDistance(500);
         viz.setPolarAngleDeg(60);
         viz.setAzimuthalAngleDeg(i);

         viz.suggestFOV();
         viz.generate(w1, buffer);
         std::string name = "surf_" + std::to_string(i) + ".png";
         viz.savePNG(name, buffer);
         // writeImage(buffer, name);
     }*/
}

int main(int argc, char* argv[])
{
    std::cout << "Testing triangulated mesh\n";
    bool success = true;
    success = success && testOpenSurface();
    // testScoring();
    // testMeshPlaneVisualization();
    // testMeshVisualization();

    /*
        std::cout << "Testing dose scoring of mesh\n";
        bool success = true;
        const auto dmesh = testScoring<true>();
        const auto dbox = testScoring<false>();
        const auto ddiff = (dmesh / dbox - 1) * 100;
        success = success && std::abs(ddiff) < 0.1;
        if (success)
            std::cout << "SUCCESS: ";
        else
            std::cout << "FAILURE: ";
        std::cout << "Dose mesh: " << dmesh << ", dose box: ";
        std::cout << dbox << ", difference[%] " << ddiff << std::endl;

    if (success)
        return EXIT_SUCCESS;
        */
    return EXIT_FAILURE;
}

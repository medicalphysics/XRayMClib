

#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/triangulatedopensurface.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    dxmc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

std::vector<dxmc::Triangle> getPyramid()
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

    std::vector<dxmc::Triangle> t;
    t.push_back({ p[0], p[1], p[4] });
    t.push_back({ p[1], p[2], p[4] });
    t.push_back({ p[2], p[3], p[4] });
    t.push_back({ p[3], p[0], p[4] });

    // underside
    t.push_back({ p[0], p[3], p[2] });
    t.push_back({ p[2], p[1], p[0] });

    return t;
}

std::vector<dxmc::Triangle> getBox(double scale = 1)
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

    std::vector<dxmc::Triangle> t;
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
std::vector<dxmc::Triangle> getPlane(double scale = 1)
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

    std::vector<dxmc::Triangle> t;
    t.push_back({ p[0], p[1], p[2] });
    t.push_back({ p[3], p[4], p[5] });
    return t;
}

template <std::size_t N = 5, int L = 2>
void testMeshVisualization()
{

    using Mesh = dxmc::TriangulatedMesh<N, L>;
    using World = dxmc::World<Mesh>;
    using Material = dxmc::Material<N>;

    World world;

    const auto waterComp = dxmc::NISTMaterials::Composition("Water, Liquid");
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

    dxmc::VisualizeWorld viz(world);

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

template <std::size_t N = 5, int L = 2>
void testMeshPlaneVisualization()
{

    using Plane = dxmc::TriangulatedOpenSurface<N, L>;
    using World = dxmc::World<Plane>;
    using Material = dxmc::Material<N>;

    World world;

    const auto waterComp = dxmc::NISTMaterials::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    const auto tri = getPlane(5);
    world.reserveNumberOfItems(1);
    auto& plane = world.template addItem<Plane>({ tri });
    plane.setMaterial(water);

    world.build();

    dxmc::VisualizeWorld viz(world);

    dxmc::IsotropicMonoEnergyBeam<true> beam;
    beam.setPosition({ 0, 0, 10 });
    beam.setDirectionCosines({ -1, 0, 0, 0, 1, 0 });
    beam.setNumberOfExposures(24);
    beam.setNumberOfParticlesPerExposure(100);
    beam.setCollimationHalfAngles({ 0.1, 0.1 });

    dxmc::Transport transport;
    transport.setNumberOfThreads(1);
    transport(world, beam);
    viz.addLineProp(beam, 100, 1);

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
    using Mesh = dxmc::TriangulatedMesh<5, 2>;
    using Box = dxmc::WorldBox<5, 2>;
    using World = dxmc::World<Mesh, Box>;
    using Material = dxmc::Material<5>;

    using Beam = dxmc::IsotropicMonoEnergyBeam<>;
    Beam beam({ 0, 0, -1000 }, { 1, 0, 0, 0, 1, 0 }, 60);
    beam.setNumberOfExposures(12);
    beam.setNumberOfParticlesPerExposure(1E5);

    const auto waterComp = dxmc::NISTMaterials::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    {
        World world;
        auto tri = getBox();
        /*std::string path = R"(C:\Users\ander\Downloads\carm.stl)";
        const std::size_t max_tree_dept = 4;
        auto& item = world.addItem<Mesh>(Mesh(path, max_tree_dept));
        */
        auto& item = world.addItem<Mesh>({ tri });
        item.setMaterial(water, 1);
        world.build();
        dxmc::Transport transport;
        auto milli = runDispatcher(transport, world, beam);
        std::cout << "mesh: " << item.doseScored().dose() << " time: " << milli << std::endl;
    }
    {
        World world;
        auto& item = world.addItem<Box>({ 1 });
        item.setMaterial(water, 1);
        world.build();
        dxmc::Transport transport;
        auto milli = runDispatcher(transport, world, beam);
        std::cout << "box: " << item.doseScored().dose() << " time: " << milli << std::endl;
    }

    return 0;
}

bool testOpenSurface()
{
    std::vector<dxmc::Triangle> tris;
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

    dxmc::World<dxmc::TriangulatedOpenSurface<>> world;

    world.template addItem<dxmc::TriangulatedOpenSurface<>>(tris);
    world.build(0);

    dxmc::VisualizeWorld viz(world);
    auto buffer = viz.createBuffer();

    for (std::size_t i = 0; i < 180; i += 30) {
        viz.setDistance(500);
        viz.setPolarAngleDeg(60);
        viz.setAzimuthalAngleDeg(i);

        viz.suggestFOV();
        viz.generate(world, buffer);
        std::string name = "surf_" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer);
        // writeImage(buffer, name);
    }
    return false;
}

int main(int argc, char* argv[])
{
    std::cout << "Testing triangulated mesh\n";
    // testOpenSurface();
    testScoring();
    testMeshPlaneVisualization();
    testMeshVisualization();

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

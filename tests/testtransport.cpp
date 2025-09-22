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

#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    xraymc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

void saveBinaryArray(const std::vector<double>& data, const std::string& name)
{
    auto myfile = std::fstream(name, std::ios::out | std::ios::binary);
    const auto bytes = data.size() * sizeof(double);
    myfile.write((char*)&data[0], bytes);
    myfile.close();
}

bool testTriangularMesh()
{
    using Mesh = xraymc::TriangulatedMesh<5, 1>;
    using World = xraymc::World<Mesh>;
    using Triangle = xraymc::Triangle;
    using Box = xraymc::WorldBox<5, 1>;

    std::vector<Triangle> triangles;
    /* triangles.push_back({ { 1, 1, 1 }, { -1, 1, 1 }, { -1, -1, 1 } });
    triangles.push_back({ { 1, 1, 1 }, { -1, -1, 1 }, { 1, -1, 1 } });
    triangles.push_back({ { 1, -1, -1 }, { 1, -1, 1 }, { -1, -1, 1 } });
    triangles.push_back({ { 1, -1, -1 }, { -1, -1, 1 }, { -1, -1, -1 } });
    triangles.push_back({ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, 1 } });
    triangles.push_back({ { -1, -1, -1 }, { -1, 1, 1 }, { -1, 1, -1 } });
    triangles.push_back({ { -1, 1, -1 }, { 1, 1, -1 }, { 1, -1, -1 } });
    triangles.push_back({ { -1, 1, -1 }, { 1, -1, -1 }, { -1, -1, -1 } });
    triangles.push_back({ { 1, 1, -1 }, { 1, 1, 1 }, { 1, -1, 1 } });
    triangles.push_back({ { 1, 1, -1 }, { 1, -1, 1 }, { 1, -1, -1 } });
    triangles.push_back({ { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, 1 } });
    triangles.push_back({ { -1, 1, -1 }, { 1, 1, 1 }, { 1, 1, -1 } });
    */

    triangles.push_back({ { 0, 0, 0 }, { 5, 0, 0 }, { 0, 5, 0 } });
    triangles.push_back({ { 0, 0, 0 }, { 5, 0, 0 }, { 0, 0, 5 } });
    triangles.push_back({ { 5, 0, 0 }, { 0, 5, 0 }, { 0, 0, 5 } });
    triangles.push_back({ { 0, 0, 0 }, { 0, 0, 5 }, { 0, 5, 0 } });

    World world;
    // Mesh mesh(triangles);
    auto& mesh = world.template addItem<Mesh>({ triangles });
    mesh.setNistMaterial("Water, Liquid");
    world.build();

    xraymc::World<Box> worldBox;
    auto& box = worldBox.template addItem<Box>({ 1 });
    box.setNistMaterial("Water, Liquid");
    worldBox.build();

    using Beam = xraymc::IsotropicMonoEnergyBeam<>;

    Beam beam;
    beam.setPosition({ 0, 0, -1 });
    beam.setDirectionCosines({ 1, 0, 0, 0, 1, 0 });
    beam.setEnergy(60);
    beam.setNumberOfExposures(32);
    beam.setNumberOfParticlesPerExposure(1E5);

    xraymc::Transport transport;
    // transport(world, beam);
    // transport(worldBox, beam);

    const auto& dose = mesh.doseScored();
    const auto& doseBox = box.doseScored();

    std::cout << "Mesh: " << dose.dose() << ", [" << dose.standardDeviation() << "]\n";
    std::cout << "Box: " << doseBox.dose() << ", [" << doseBox.standardDeviation() << "]\n";
    auto test = (dose.dose() - doseBox.dose()) / std::sqrt((dose.variance() + doseBox.variance()) / 2);
    std::cout << "t-test: " << test << "\n";

    xraymc::VisualizeWorld viz(world);
    viz.setPolarAngle(std::numbers::pi_v<double> * 1.0f / 3.0f);
    viz.setAzimuthalAngle(std::numbers::pi_v<double> * 1.0f / 3.0f);
    viz.setDistance(60);
    viz.setCameraPosition({ 30, 30, -10 });
    viz.suggestFOV(5);
    int height = 1024;
    int width = 1024;
    std::vector<double> buffer(height * width * 4, 1);
    viz.generate(world, buffer, width, height);
    saveBinaryArray(buffer, "color.bin");

    return false;
}

bool testCylinder()
{
    using Cylinder = xraymc::WorldCylinder<5, 1>;

    xraymc::World<Cylinder> world;

    auto& cylinder = world.template addItem<Cylinder>({ 4, 60 });

    world.build(120);

    xraymc::Transport transport;

    std::uint64_t nPart = 0;
    std::chrono::milliseconds time;

    constexpr double dist = 60;

    xraymc::IsotropicMonoEnergyBeam<> beam;
    const auto collAngleY = std::atan(2.0f / dist);
    const auto collAngleZ = std::atan(2.0f / dist);
    beam.setCollimationHalfAngles({ -collAngleY, -collAngleZ, collAngleY, collAngleZ });

    beam.setNumberOfParticlesPerExposure(1e6);
    beam.setNumberOfExposures(32);

    std::cout << "Angle, EnergyImparted, nEvents, StdDev*2" << std::endl;

    std::vector<xraymc::EnergyScore> res(8 * 3);
    const double angStep = 360 / res.size();
    for (int ang = 0; ang < res.size(); ++ang) {
        const auto a = xraymc::DEG_TO_RAD() * ang * angStep;
        constexpr std::array<double, 3> pos = { -dist, 0, 0 };
        constexpr std::array<double, 3> cosy = { 0, 1, 0 };
        constexpr std::array<double, 3> cosz = { 0, 0, 1 };

        auto rpos = xraymc::vectormath::rotate(pos, { 0, 0, 1 }, a);
        auto rcosy = xraymc::vectormath::rotate(cosy, { 0, 0, 1 }, a);
        auto rcosz = xraymc::vectormath::rotate(cosz, { 0, 0, 1 }, a);
        beam.setPosition(rpos);
        beam.setDirectionCosines(rcosy, rcosz);

        time = runDispatcher(transport, world, beam);
        res[ang] = cylinder.energyScored();

        cylinder.clearEnergyScored();

        std::cout << ang * angStep << ", " << res[ang].energyImparted() << ", " << res[ang].numberOfEvents();
        std::cout << ", " << res[ang].standardDeviation() << std::endl;
    }

    return true;
}

bool testCTDI()
{
    using CTDI = xraymc::CTDIPhantom<5, 1>;
    using Box = xraymc::WorldBox<4, 1>;

    xraymc::World<CTDI, Box> world;
    auto& phantom = world.template addItem<CTDI>({});

    phantom.setHoleMaterial("Polymethyl Methacralate (Lucite, Perspex)", xraymc::NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)"));

    world.build(30);

    xraymc::Transport transport;

    std::uint64_t nPart = 0;
    std::chrono::milliseconds time;

    xraymc::IsotropicMonoEnergyBeam<> beam;

    const auto collAngleY = std::atan(16.0 / 60);
    const auto collAngleZ = 0; //    std::atan(4.0f / 60);
    beam.setCollimationHalfAngles({ -collAngleY, -collAngleZ, collAngleY, collAngleZ });
    beam.setNumberOfParticlesPerExposure(1e6);
    beam.setNumberOfExposures(32);

    std::vector<std::array<double, 6>> res(8);
    constexpr double angStep = 45;
    for (int ang = 0; ang < res.size(); ++ang) {
        const auto a = xraymc::DEG_TO_RAD() * ang * angStep;
        constexpr std::array<double, 3> pos = { -60, 0, 0 };
        constexpr std::array<double, 3> cosy = { 0, 1, 0 };
        constexpr std::array<double, 3> cosz = { 0, 0, 1 };

        auto rpos = xraymc::vectormath::rotate(pos, { 0, 0, 1 }, a);
        auto rcosy = xraymc::vectormath::rotate(cosy, { 0, 0, 1 }, a);
        auto rcosz = xraymc::vectormath::rotate(cosz, { 0, 0, 1 }, a);
        beam.setPosition(rpos);
        beam.setDirectionCosines(rcosy, rcosz);

        time = runDispatcher(transport, world, beam);
        for (int i = 0; i < 6; ++i) {
            res[ang][i] = (phantom.energyScored(i).energyImparted() * 1000) / beam.numberOfParticles();
        }

        std::cout << ang * angStep << ", ";
        for (auto a : res[ang]) {
            std::cout << a << ", ";
        }
        std::cout << std::endl;
        phantom.clearEnergyScored();
    }

    for (int i = 0; i < res.size(); ++i) {
        std::cout << i * angStep << ", ";
        for (auto a : res[i]) {
            std::cout << a << ", ";
        }
        std::cout << std::endl;
    }

    return true;
}

bool testDepth(bool print = false)
{
    using Cylinder = xraymc::DepthDose<5, 2>;
    using World = xraymc::World<Cylinder>;
    using Beam = xraymc::PencilBeam<>;

    World world;
    auto& cylinder = world.template addItem<Cylinder>({ 0.1, 20, 20 });
    auto material = xraymc::Material<5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)");
    if (material) {
        cylinder.setMaterial(material.value());
        cylinder.setMaterialDensity(1.19);
    }
    world.build();

    const double energy = 60;

    Beam beam;
    beam.setEnergy(energy);
    beam.setPosition({ 0, 0, -100 });
    beam.setDirection({ 0, 0, 1 });

    beam.setNumberOfExposures(40);
    beam.setNumberOfParticlesPerExposure(1e4);

    // beam.setParticleWeight(T { 0.5 });

    xraymc::Transport transport;
    auto time = runDispatcher(transport, world, beam);

    if (print)
        std::cout << "pos, energy, std, nevents\n";

    double dose0 = -1;
    double pos0 = -1;
    const auto att = material.value().attenuationValues(beam.energy()).sum();

    bool success = true;

    for (const auto& [pos, d] : cylinder.depthEnergyScored()) {
        if (dose0 < 0) {
            dose0 = d.energyImparted();
            pos0 = pos;
        }
        if (print) {
            std::cout << pos << ", ";
            std::cout << d.energyImparted() << ", ";
            std::cout << d.standardDeviation() << ", ";
            std::cout << d.numberOfEvents() << "\n";
        }
        success = success && d.energyImparted() / dose0 - std::exp(-(pos - pos0) * att) < 0.01;
    }
    return success;
}

template <typename T>
void writeImage(const std::vector<T>& buffer, const std::string& name)
{
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

template <typename T = int, typename U>
std::vector<T> generateDonut(const std::array<std::size_t, 3>& dim, const std::array<U, 3>& spacing = { 1, 1, 1 })
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });

    const U R = U { 0.25 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;
    const U r = U { 0.1 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const auto flat_ind = x + y * dim[0] + z * dim[0] * dim[1];
                const auto xc = x * spacing[0] + spacing[0] / 2 - (dim[0] * spacing[0]) / 2;
                const auto yc = y * spacing[1] + spacing[0] / 2 - (dim[1] * spacing[1]) / 2;
                const auto zc = z * spacing[2] + spacing[0] / 2 - (dim[2] * spacing[2]) / 2;

                const auto p1 = R - std::sqrt(xc * xc + yc * yc);
                if (p1 * p1 + zc * zc < r * r)
                    d[flat_ind] = 1;
            }
    return d;
}

template <typename T = int>
std::vector<T> generateEdges(const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing = { 1, 1, 1 })
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });

    const auto dim_max = std::max(dim[0], std::max(dim[1], dim[2]));

    const std::size_t c0 = 3; // dim_max * 1 / 4;
    const std::size_t c1 = dim_max - 4; // dim_max * 3 / 4;

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const auto flat_ind = x + y * dim[0] + z * dim[0] * dim[1];
                if (c0 <= x && x <= c1 && (y == c0 || y == c1) && (z == c0 || z == c1)) {
                    d[flat_ind] = 1;
                }
                if (c0 <= y && y <= c1 && (x == c0 || x == c1) && (z == c0 || z == c1)) {
                    d[flat_ind] = 1;
                }
                if (c0 <= z && z <= c1 && (y == c0 || y == c1) && (x == c0 || x == c1)) {
                    d[flat_ind] = 1;
                }
            }
    return d;
}

template <std::uint_fast8_t TRANSPARENT = 0>
bool testAAVoxelGridTransport()
{
    bool success = true;

    using AAVoxelGrid = xraymc::AAVoxelGrid<5, 1, TRANSPARENT>;
    using Cylinder = xraymc::WorldCylinder<5, 2>;
    using Sphere = xraymc::WorldSphere<5, 2>;
    using World = xraymc::World<AAVoxelGrid, Sphere>;

    World world;
    auto& grid = world.addItem(AAVoxelGrid());
    auto& sphere = world.addItem(Sphere(3));

    auto air = xraymc::Material<5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = xraymc::Material<5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = xraymc::NISTMaterials::density("Air, Dry (near sea level)");
    const auto pmma_dens = xraymc::NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");

    sphere.setMaterial(pmma);
    sphere.setMaterialDensity(pmma_dens);

    const std::array<std::size_t, 3> dim = { 64, 64, 64 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<double, 3> spacing = { .5f, .5f, .5f };

    // setting up grid
    auto matInd = generateDonut<std::uint8_t>(dim, spacing);
    // auto matInd = generateEdges<std::uint8_t>(dim, spacing);
    std::vector<double> dens(size, air_dens);
    std::transform(matInd.cbegin(), matInd.cend(), dens.begin(), [=](const auto i) { return i == 0 ? air_dens : pmma_dens; });
    std::vector<xraymc::Material<5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);
    grid.setData(dim, dens, matInd, materials);
    grid.setSpacing(spacing);

    xraymc::PencilBeam<> beam({ 0, 0, -100 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(64);
    beam.setNumberOfParticlesPerExposure(100000);
    xraymc::Transport transport;
    world.build();
    auto time = runDispatcher(transport, world, beam);
    std::cout << std::format("Total time: {}", time) << std::endl;

    std::vector<double> doseArray(dens.size());
    for (std::size_t i = 0; i < size; ++i) {
        const auto& dose = grid.energyScored(i);
        doseArray[i] = dose.energyImparted();
    }
    writeImage(doseArray, "dose.bin");

    return success;
}

bool testAAVoxelGrid()
{
    xraymc::AAVoxelGrid<5> item;

    auto air = xraymc::Material<5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = xraymc::Material<5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = xraymc::NISTMaterials::density("Air, Dry (near sea level)");
    const auto pmma_dens = xraymc::NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");

    const std::array<std::size_t, 3> dim = { 64, 64, 64 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<double, 3> spacing = { .20f, .20f, .20f };

    // material arrays
    std::vector<double> dens(size, air_dens);
    std::vector<std::uint8_t> materialIdx(size, 0);
    std::vector<xraymc::Material<5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);

    // test indices and assign material
    bool success = true;
    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const std::array tind = { x, y, z };
                const auto find = item.flatIndex(tind);
                const auto ind = item.index(find);
                success = success && ind == tind;
                /* if (x > dim[0] / 2 && y > dim[1] / 2 && z > dim[2] / 2) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                } else if (x < dim[0] / 2 && y < dim[1] / 2 && z < dim[2] / 2) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                }*/
                if (x == 1 && y == 1 && z == 1) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                }
            }

    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    xraymc::Particle p;

    p.pos = { -100, 0, 0 };
    p.dir = { 1, 0, 0 };
    auto res = item.intersect(p);
    success = success && res.valid() && res.intersection == 99;

    p.pos = { 100, 0, 0 };
    p.dir = { -1, 0, 0 };
    res = item.intersect(p);
    success = success && res.valid() && res.intersection == 99;

    p.pos = { -100, -100, -100 };
    p.dir = { 1, 1, 1 };
    xraymc::vectormath::normalize(p.dir);
    res = item.intersect(p);
    auto val = std::sqrt(3 * 100 * 100.0) - std::sqrt(3 * spacing[0] * spacing[0] / 4);
    success = success && res.valid() && val - res.intersection < 1E-6;

    return success;
}

int main()
{
    bool success = true;

    success = success && testTriangularMesh();

    success = success && testCylinder();
    success = success && testCTDI();

    success = success && testAAVoxelGridTransport<0>();
    success = success && testAAVoxelGridTransport<255>();

    success = success && testAAVoxelGrid();

    success = success && testDepth();

    if (success)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}
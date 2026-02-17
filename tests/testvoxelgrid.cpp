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
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>
#include <string>
#include <vector>

template <typename T>
void writeImage(const std::vector<T>& buffer, const std::string& name)
{
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

template <typename T = int, typename U>
std::vector<T> generateSphere(const std::array<std::size_t, 3>& dim, const std::array<U, 3>& spacing = { 1, 1, 1 })
{
    const double R = dim[0] / 2;

    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });
    for (int x = 0; x < dim[0]; ++x)
        for (int y = 0; y < dim[1] - 1; ++y)
            for (int z = 0; z < dim[2]; ++z) {
                const auto xc = x * spacing[0] - (dim[0] * spacing[0]) / 2 + spacing[0] / 2;
                const auto yc = y * spacing[1] - (dim[1] * spacing[1]) / 2 + spacing[1] / 2;
                const auto zc = z * spacing[2] - (dim[2] * spacing[2]) / 2 + spacing[2] / 2;
                bool k = xc * xc + yc * yc + zc * zc < R * R;
                d[x + y * dim[0] + z * dim[0] * dim[1]] = k ? 1 : 0;
            }
    return d;
}

template <typename T = int, typename U>
std::vector<T> generateBox(const std::array<std::size_t, 3>& dim, const std::array<U, 3>& spacing = { 1, 1, 1 }, bool empty = false)
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    if (empty) {
        std::vector<T> d(s, T { 0 });
        for (int x = 0; x < dim[0]; ++x)
            for (int y = 0; y < dim[1]; ++y)
                for (int z = 0; z < dim[2]; ++z)
                    if (y % 2 == 1)
                        d[x + y * dim[0] + z * dim[0] * dim[1]] = 1;
        return d;
    } else {
        std::vector<T> d(s, T { 1 });
        for (int x = 1; x < dim[0] - 1; ++x)
            for (int y = 1; y < dim[1] - 1; ++y)
                for (int z = 1; z < dim[2] - 1; ++z)
                    d[x + y * dim[0] + z * dim[0] * dim[1]] = 0;

        return d;
    }
}

template <typename T = int, typename U>
std::vector<T> generateDonut(const std::array<std::size_t, 3>& dim, const std::array<U, 3>& spacing = { 1, 1, 1 })
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });

    const U R = U { 0.3 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;
    const U r = U { 0.2 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const auto flat_ind = x + y * dim[0] + z * dim[0] * dim[1];
                const auto xc = x * spacing[0] - (dim[0] * spacing[0]) / 2 + spacing[0] / 2;
                const auto yc = y * spacing[1] - (dim[1] * spacing[1]) / 2 + spacing[1] / 2;
                const auto zc = z * spacing[2] - (dim[2] * spacing[2]) / 2 + spacing[2] / 2;

                const auto p1 = R - std::sqrt(xc * xc + yc * yc);
                if (p1 * p1 + zc * zc < r * r)
                    d[flat_ind] = 1;
            }
    return d;
}

template <std::uint8_t TRANSPARENTVOXELS = 0>
bool dose_test()
{
    constexpr std::size_t N = 16;

    const std::string suff = std::to_string(TRANSPARENTVOXELS) + ".png";

    xraymc::World<xraymc::AAVoxelGrid<N, 2, TRANSPARENTVOXELS>, xraymc::FlatDetector, xraymc::WorldSphere<N, 2, true>> world;
    world.reserveNumberOfItems(3);

    auto mat = generateDonut<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 });
    // auto mat = generateSphere<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 });
    //  auto mat = generateBox<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 });

    std::vector<double> dens(mat.size());
    std::vector<xraymc::Material<N>> materials;
    materials.push_back(xraymc::Material<N>::byNistName("Air, Dry (near sea level)").value());
    materials.push_back(xraymc::Material<N>::byNistName("Water, Liquid").value());
    std::transform(mat.cbegin(), mat.cend(), dens.begin(), [](const auto v) { return v == 0 ? 0.001 : 1.0; });
    auto& donut = world.addItem<xraymc::AAVoxelGrid<N, 2, TRANSPARENTVOXELS>>({ { 64, 64, 64 }, { 1.0, 1.0, 1.0 }, dens, mat, materials });

    auto& sphere = world.addItem<xraymc::WorldSphere<N, 2, true>>("Sphere");
    sphere.setCenter({ 0, 0, 40 });
    sphere.setRadius(7.5);
    sphere.setMaterial(materials[1], 1.0);

    auto& detector = world.addItem<xraymc::FlatDetector>("Detector");
    detector.setCenter({ 0, 0, 64.5 });
    detector.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    detector.setDetectorDimensions({ 512, 512 });
    const double spacing = 64.0 / 512.0;
    detector.setPixelSpacing({ spacing, spacing });
    world.setMaterialDensity(0.001);
    world.build();

    xraymc::IsotropicMonoEnergyBeam<> beam;
    beam.setEnergy(78);
    beam.setPosition({ 0, 0, -200 });
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    const auto SID = beam.position()[2] + 64 * 1.0 / 2.0 + 1.0;
    const auto angle = std::tan(32.0 / SID);

    beam.setCollimationHalfAngles(angle, angle);
    beam.setNumberOfExposures(100);
    beam.setNumberOfParticlesPerExposure(1000000);

    xraymc::VisualizeWorld viz(world);
    viz.setPolarAngleDeg(60);
    viz.setAzimuthalAngleDeg(60);
    viz.setDistance(200);
    viz.suggestFOV(1);
    auto buffer = viz.createBuffer(2048, 2048);
    viz.addLineProp(beam, 200, 1);
    viz.generate(world, buffer);

    viz.savePNG("test" + suff, buffer);

    xraymc::Transport transport;
    std::cout << "Running transport...";
    auto time = transport.runConsole(world, beam, 0, true, 2000);
    std::cout << " finished in " << time;

    std::vector<std::uint8_t> image_buffer(detector.detectorDimensions()[0] * detector.detectorDimensions()[1] * 4);
    std::vector<double> dose_buffer(detector.detectorDimensions()[0] * detector.detectorDimensions()[1]);
    for (std::size_t i = 0; i < dose_buffer.size(); ++i) {
        dose_buffer[i] = detector.doseScored(i).dose();
    }

    std::cout << " Total detector dose: " << std::accumulate(dose_buffer.cbegin(), dose_buffer.cend(), 0.0) << " mGy";
    std::cout << " Sphere dose: " << sphere.doseScored().dose() << " mGy\n";

    const auto max_dose = *std::max_element(dose_buffer.cbegin(), dose_buffer.cend());
    const auto min_dose = *std::min_element(dose_buffer.cbegin(), dose_buffer.cend());
    for (std::size_t i = 0; i < dose_buffer.size(); ++i) {
        const auto d = static_cast<std::uint8_t>(255.0 * (dose_buffer[i] - min_dose) / (max_dose - min_dose));
        image_buffer[i * 4 + 0] = d;
        image_buffer[i * 4 + 1] = d;
        image_buffer[i * 4 + 2] = d;
        image_buffer[i * 4 + 3] = 255;
    }
    xraymc::xraymclodepng::savePNG("image" + suff, image_buffer, detector.detectorDimensions()[0], detector.detectorDimensions()[1]);

    viz.addColorByValueItem(world.getItemPointerFromName("Detector"));
    viz.setColorByValueMinMax(min_dose, max_dose);
    viz.generate(world, buffer);
    viz.savePNG("dose" + suff, buffer);

    return true;
}

template <std::uint8_t TRANSPARENTVOXELS = 0>
bool pencil_test()
{
    constexpr std::size_t N = 16;

    const std::string suff = std::to_string(TRANSPARENTVOXELS) + ".png";

    xraymc::World<xraymc::AAVoxelGrid<N, 2, TRANSPARENTVOXELS>, xraymc::FlatDetector, xraymc::WorldSphere<N, 2, true>> world;
    world.reserveNumberOfItems(3);

    // auto mat = generateDonut<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 });
    // auto mat = generateSphere<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 });
    auto mat = generateBox<std::uint8_t, double>({ 64, 64, 64 }, { 1.0, 1.0, 1.0 }, true);

    std::vector<double> dens(mat.size());
    std::vector<xraymc::Material<N>> materials;
    materials.push_back(xraymc::Material<N>::byNistName("Air, Dry (near sea level)").value());
    materials.push_back(xraymc::Material<N>::byNistName("Water, Liquid").value());
    std::transform(mat.cbegin(), mat.cend(), dens.begin(), [](const auto v) { return v == 0 ? 0.001 : 1.0; });
    auto& donut = world.addItem<xraymc::AAVoxelGrid<N, 2, TRANSPARENTVOXELS>>({ { 64, 64, 64 }, { 1.0, 1.0, 1.0 }, dens, mat, materials });

    auto& sphere = world.addItem<xraymc::WorldSphere<N, 2, true>>("Sphere");
    sphere.setCenter({ 0, 40, 0 });
    sphere.setRadius(2);
    sphere.setMaterial(materials[0], 0.001);

    auto& detector = world.addItem<xraymc::FlatDetector>("Detector");
    detector.setCenter({ 0, 64.5, 0 });
    detector.setDirectionCosines({ 1, 0, 0 }, { 0, 0, 1 });
    detector.setDetectorDimensions({ 512, 512 });
    const double spacing = 64.0 / 512.0;
    detector.setPixelSpacing({ spacing, spacing });
    world.setMaterialDensity(0.001);

    xraymc::PencilBeam beam;
    beam.setEnergy(100);
    beam.setPosition({ -50, -100, -30 });
    const auto beam_dir = xraymc::vectormath::scale(beam.position(), -1.0);
    beam.setDirection(beam_dir);
    beam.setNumberOfExposures(100);
    beam.setNumberOfParticlesPerExposure(1000000);

    detector.setCenter({ 50, 100, 30 });
    detector.setDirectionCosines(xraymc::vectormath::cross(beam.direction(), { 0, 0, 1 }), { 0, 0, 1 });
    sphere.setCenter(xraymc::vectormath::add(beam.position(), xraymc::vectormath::scale(beam.direction(), 200.0)));

    world.build();

    xraymc::VisualizeWorld viz(world);
    viz.setPolarAngleDeg(60);
    viz.setAzimuthalAngleDeg(40);
    viz.setDistance(200);
    viz.suggestFOV(1);
    auto buffer = viz.createBuffer(2048, 2048);
    viz.addLineProp(beam.position(), beam.direction(), 200, 1);

    viz.generate(world, buffer);

    viz.savePNG("test" + suff, buffer);

    xraymc::Transport transport;
    std::cout << "Running transport...";
    auto time = transport.runConsole(world, beam, 0, true, 2000);
    std::cout << " finished in " << time;

    std::vector<std::uint8_t> image_buffer(detector.detectorDimensions()[0] * detector.detectorDimensions()[1] * 4);
    std::vector<double> dose_buffer(detector.detectorDimensions()[0] * detector.detectorDimensions()[1]);
    for (std::size_t i = 0; i < dose_buffer.size(); ++i) {
        dose_buffer[i] = detector.doseScored(i).dose();
    }

    std::cout << " Total detector dose: " << std::accumulate(dose_buffer.cbegin(), dose_buffer.cend(), 0.0) << " mGy";
    const auto sphere_dose = sphere.doseScored().dose();
    std::cout << " Total sphere dose: " << sphere_dose << " mGy" << std::endl;

    const auto max_dose = *std::max_element(dose_buffer.cbegin(), dose_buffer.cend());
    const auto min_dose = *std::min_element(dose_buffer.cbegin(), dose_buffer.cend());
    for (std::size_t i = 0; i < dose_buffer.size(); ++i) {
        const auto d = static_cast<std::uint8_t>(255.0 * (dose_buffer[i] - min_dose) / (max_dose - min_dose));
        image_buffer[i * 4 + 0] = d;
        image_buffer[i * 4 + 1] = d;
        image_buffer[i * 4 + 2] = d;
        image_buffer[i * 4 + 3] = 255;
    }
    xraymc::xraymclodepng::savePNG("image" + suff, image_buffer, detector.detectorDimensions()[0], detector.detectorDimensions()[1]);

    viz.addColorByValueItem(world.getItemPointerFromName("Detector"));
    viz.setColorByValueMinMax(min_dose, max_dose);
    viz.generate(world, buffer);
    viz.savePNG("dose" + suff, buffer);

    return true;
}

int main()
{

    bool success = true;
    success = success && dose_test<0>();
    success = success && dose_test<255>();
    // success = success && pencil_test<0>();
    // success = success && pencil_test<255>();
    //  testGeometryColor();
    //  testGeometryDistance();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
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

Copyright 2026 Erlend Andersen
*/

#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#ifdef XRAYMCLIB_USE_LOADPNG
#include "lodepng/xraymclodepngwrapper.hpp"
#endif

template <typename U>
    requires(std::same_as<U, double> || std::same_as<U, std::uint8_t>)
static bool savePNG(const std::string& filename, const std::vector<U>& buffer, std::size_t width, std::size_t height)
{
#ifdef XRAYMCLIB_USE_LOADPNG
    return xraymc::xraymclodepng::savePNG(filename, buffer, width, height);
#else
    return false;
#endif
}

void testDetector()
{

    using Sphere = xraymc::WorldSphere<16, 2, false>;
    using Detector = xraymc::FlatDetector;
    using Material = xraymc::Material<16>;

    auto water = Material::byChemicalFormula("H2O").value();

    xraymc::World<Sphere, Detector> world;
    world.reserveNumberOfItems(2);

    auto& sphere = world.addItem<Sphere>();
    sphere.setCenter({ 0, 0, 10 });
    sphere.setRadius(2);
    sphere.setMaterial(water, 1.0);

    auto& detector = world.addItem<Detector>();
    detector.setCenter({ 0, 0, 20 });
    detector.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    detector.setPixelSpacing({ .1, .1 });
    detector.setDetectorDimensions({ 100, 100 });

    world.build(100);

    xraymc::DXBeam<false> beam({ 0, 0, -30 });
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    const auto angle = std::tan(3 / 50.0);
    beam.setCollimationHalfAngles({ angle, angle });
    beam.setTubeVoltage(120);
    beam.setNumberOfExposures(36);
    beam.setNumberOfParticlesPerExposure(1e3);

    xraymc::Transport transport;
    transport.runConsole(world, beam);

    std::vector<double> buffer, dose;
    const auto nPixels = detector.detectorDimensions()[0] * detector.detectorDimensions()[1];
    dose.resize(nPixels);
    for (std::size_t i = 0; i < nPixels; ++i) {
        dose[i] = detector.doseScored(i).dose();
    }
    auto max = std::max_element(dose.cbegin(), dose.cend());

    std::cout << "Max dose on detector: " << *max << " mGy" << std::endl;
    buffer.resize(nPixels * 4);
    for (std::size_t i = 0; i < nPixels; ++i) {
        const auto energy = 1 - dose[i] / (*max);
        buffer[i * 4] = energy;
        buffer[i * 4 + 1] = energy;
        buffer[i * 4 + 2] = energy;
        buffer[i * 4 + 3] = 1.0;
    }
    bool test = savePNG("detector_output.png", buffer, detector.detectorDimensions()[0], detector.detectorDimensions()[1]);
    return;
}

int main()
{
    testDetector();
    return 0;
}
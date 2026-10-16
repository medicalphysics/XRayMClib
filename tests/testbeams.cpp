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

#include "testutils.hpp"

#include "xraymc/beams/beamtype.hpp"
#include "xraymc/beams/cbctbeam.hpp"
#include "xraymc/beams/ctsequentialbeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/ctspiraldualenergybeam.hpp"
#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/flatmonoenergyfield.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/isotropicbeamcircle.hpp"
#include "xraymc/beams/isotropiccircularbeam.hpp"
#include "xraymc/beams/isotropiccircularmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/vectormath.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using Vec = std::array<double, 3>;

namespace vm = xraymc::vectormath;

// Generic properties of particles sampled from every exposure of a beam
template <BeamType B>
void checkBeamParticles(const B& beam, const std::string& name, double emin, double emax,
    const std::function<bool(std::size_t, const Particle&)>& geometry = {}, std::size_t particlesPerExposure = 2000)
{
    auto state = xraymctest::seededState();
    bool unit = true, energy = true, weight = true, position = true, geom = true;
    std::uint64_t totalParticles = 0;
    const auto nExp = beam.numberOfExposures();
    const std::size_t stride = std::max<std::size_t>(1, nExp / 20);
    for (std::size_t i = 0; i < nExp; ++i) {
        const auto exp = beam.exposure(i);
        totalParticles += exp.numberOfParticles();
        if (i % stride != 0)
            continue;
        for (std::size_t j = 0; j < particlesPerExposure; ++j) {
            const auto p = exp.sampleParticle(state);
            unit = unit && std::abs(vm::length(p.dir) - 1) < 1E-9;
            energy = energy && p.energy >= emin && p.energy <= emax;
            weight = weight && p.weight > 0 && std::isfinite(p.weight);
            if (geometry) {
                Particle q { .pos = p.pos, .dir = p.dir, .energy = p.energy, .weight = p.weight };
                geom = geom && geometry(i, q);
            }
        }
    }
    check(totalParticles == beam.numberOfParticles(), name + ": numberOfParticles equals sum over exposures");
    check(unit, name + ": particle directions are unit vectors");
    check(energy, name + ": particle energies within [" + std::to_string(emin) + ", " + std::to_string(emax) + "] keV");
    check(weight, name + ": particle weights positive and finite");
    if (geometry)
        check(geom, name + ": particle geometry within beam definition");
    const auto cal = beam.calibrationFactor();
    check(cal > 0 && std::isfinite(cal), name + ": calibration factor positive and finite");
}

// Angles of direction d projected on the two direction cosines relative to the central direction
std::array<double, 2> fieldAngles(const Vec& d, const Vec& cx, const Vec& cy, const Vec& dir)
{
    const auto z = vm::dot(d, dir);
    return { std::atan(vm::dot(d, cx) / z), std::atan(vm::dot(d, cy) / z) };
}

void testPencilBeam()
{
    PencilBeam<> beam({ 1, 2, 3 }, { 0, 3, 4 }, 55);
    beam.setNumberOfExposures(10);
    beam.setNumberOfParticlesPerExposure(100);
    check(beam.numberOfParticles() == 1000, "PencilBeam: number of particles");
    checkNear(vm::length(beam.direction()), 1, 1E-12, "PencilBeam: direction is normalized");
    const auto cos = beam.directionCosines();
    checkNear(vm::dot(cos[0], cos[1]), 0, 1E-12, "PencilBeam: direction cosines are orthogonal");
    check(xraymctest::checkVectorNear(vm::cross(cos[0], cos[1]), beam.direction(), 1E-12, "PencilBeam: direction cosines span beam direction"), "PencilBeam: direction cosines consistent");

    checkBeamParticles(beam, "PencilBeam", 55, 55, [&](std::size_t, const Particle& p) {
        return p.pos == beam.position() && std::abs(vm::dot(p.dir, beam.direction()) - 1) < 1E-12;
    });

    // Calibration converts to air KERMA for unit area: K = E * N * mu_tr(E)
    beam.setAirKerma(2.0);
    const auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
    checkRelative(beam.calibrationFactor(), 2.0 / (55 * 1000 * air.massEnergyTransferAttenuation(55)), 1E-9, "PencilBeam: calibration factor from air KERMA");

    PencilBeam<true> tracked({ 0, 0, 0 }, { 0, 0, 1 }, 60);
    auto state = xraymctest::seededState();
    const auto p = tracked.exposure(0).sampleParticle(state);
    check(p.getSize() == 1, "PencilBeam: tracking particle has registered start position");
    beam.setEnergy(-70);
    checkNear(beam.energy(), 70, 0, "PencilBeam: negative energy is stored as absolute value");
}

void testFlatField()
{
    FlatMonoEnergyField<> beam({ 0, 0, -10 }, 45);
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    beam.setLenght(3, 2);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(50000);
    checkBeamParticles(beam, "FlatMonoEnergyField", 45, 45, [&](std::size_t, const Particle& p) {
        const auto d = vm::subtract(p.pos, beam.position());
        return std::abs(d[0]) <= 3 && std::abs(d[1]) <= 2 && std::abs(d[2]) < 1E-12 && p.dir == Vec { 0, 0, 1 };
    });
    // positions are uniform over the field: the variance of a uniform distribution on [-l, l] is l^2 / 3
    auto state = xraymctest::seededState();
    const auto exp = beam.exposure(0);
    double sx = 0, sxx = 0, syy = 0;
    constexpr std::size_t N = 200000;
    for (std::size_t i = 0; i < N; ++i) {
        const auto p = exp.sampleParticle(state);
        sx += p.pos[0];
        sxx += p.pos[0] * p.pos[0];
        syy += p.pos[1] * p.pos[1];
    }
    checkNear(sx / N, 0, 0.02, "FlatMonoEnergyField: field is centered on position");
    checkRelative(sxx / N, 9.0 / 3.0, 0.01, "FlatMonoEnergyField: field half width in x");
    checkRelative(syy / N, 4.0 / 3.0, 0.01, "FlatMonoEnergyField: field half width in y");

    // Calibration: air KERMA for the fluence over the full field area (2 lx) * (2 ly)
    beam.setAirKerma(1.0);
    const auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
    const auto fluence = beam.numberOfParticles() / (2 * 3.0 * 2 * 2.0);
    checkRelative(beam.calibrationFactor(), 1.0 / (45 * fluence * air.massEnergyTransferAttenuation(45)), 1E-9, "FlatMonoEnergyField: calibration uses full field area");
}

template <typename B>
void testIsotropicRectangular(B& beam, const std::string& name, double emin, double emax)
{
    const std::array<double, 4> angles = { -0.1, -0.2, 0.3, 0.05 };
    beam.setPosition({ 5, -3, 2 });
    beam.setDirectionCosines({ 0, 1, 0 }, { 0, 0, 1 }); // beam along +x
    beam.setCollimationHalfAngles(angles);
    beam.setNumberOfExposures(3);
    beam.setNumberOfParticlesPerExposure(10000);
    const auto& cos = beam.directionCosines();
    const auto dir = vm::cross(cos[0], cos[1]);
    check(xraymctest::checkVectorNear(dir, { 1, 0, 0 }, 1E-12, name + ": beam direction from direction cosines"), name + ": direction");
    checkBeamParticles(beam, name, emin, emax, [&](std::size_t, const Particle& p) {
        const auto a = fieldAngles(p.dir, cos[0], cos[1], dir);
        constexpr double eps = 1E-9;
        return p.pos == Vec { 5, -3, 2 } && a[0] >= angles[0] - eps && a[0] <= angles[2] + eps && a[1] >= angles[1] - eps && a[1] <= angles[3] + eps;
    });
}

void testIsotropicBeams()
{
    {
        IsotropicMonoEnergyBeam<> beam;
        beam.setEnergy(80);
        testIsotropicRectangular(beam, "IsotropicMonoEnergyBeam", 80, 80);
        beam.setEnergy(1000);
        checkNear(beam.energy(), MAX_ENERGY(), 0, "IsotropicMonoEnergyBeam: energy is clamped to MAX_ENERGY");
    }
    {
        Tube tube(90);
        tube.setAlFiltration(3);
        IsotropicBeam<> beam;
        beam.setEnergySpecter(tube.getSpecter());
        testIsotropicRectangular(beam, "IsotropicBeam", 0, 90);

        // mean sampled energy equals mean of the spectrum
        auto state = xraymctest::seededState();
        const auto exp = beam.exposure(0);
        double sum = 0;
        constexpr std::size_t N = 500000;
        for (std::size_t i = 0; i < N; ++i)
            sum += exp.sampleParticle(state).energy;
        // energies are sampled uniformly within each 1 keV bin starting at the bin energy
        checkNear(sum / N, tube.meanSpecterEnergy() + 0.5, 0.2, "IsotropicBeam: mean sampled energy equals spectrum mean");

        IsotropicBeam<true> tracked(beam);
        const auto p = tracked.exposure(0).sampleParticle(state);
        check(p.getSize() == 1 && p.pos == beam.position(), "IsotropicBeam: copy to tracking beam keeps geometry");
        check(tracked.numberOfParticles() == beam.numberOfParticles(), "IsotropicBeam: copy to tracking beam keeps number of particles");
    }
    for (bool mono : { true, false }) {
        const double half = 0.25;
        const Vec dir = vm::normalized(Vec { 1, 1, 0 });
        const std::string name = mono ? "IsotropicMonoEnergyBeamCircle" : "IsotropicBeamCircle";
        auto geometry = [&](std::size_t, const Particle& p) {
            return vm::angleBetween(p.dir, dir) <= half + 1E-9 && p.pos == Vec { 0, 0, 0 };
        };
        if (mono) {
            IsotropicMonoEnergyBeamCircle<> beam({ 0, 0, 0 }, dir);
            beam.setCollimationHalfAngle(half);
            beam.setEnergy(70);
            beam.setNumberOfExposures(2);
            beam.setNumberOfParticlesPerExposure(10000);
            checkBeamParticles(beam, name, 70, 70, geometry);
        } else {
            IsotropicBeamCircle<> beam({ 0, 0, 0 }, dir);
            beam.setCollimationHalfAngle(half);
            beam.setEnergySpecter(Tube(100).getSpecter());
            beam.setNumberOfExposures(2);
            beam.setNumberOfParticlesPerExposure(10000);
            checkBeamParticles(beam, name, 0, 100, geometry);
        }
    }
    for (bool mono : { true, false }) {
        // Particles start on a circle around the z axis through position and point inwards
        const Vec center = { 1, 2, 3 };
        const double radius = 50;
        const std::array<double, 4> angles = { -0.05, -0.1, 0.05, 0.1 };
        const std::string name = mono ? "IsotropicCircularMonoEnergyBeam" : "IsotropicCircularBeam";
        auto geometry = [&](std::size_t, const Particle& p) {
            const auto d = vm::subtract(p.pos, center);
            const bool onCircle = std::abs(std::sqrt(d[0] * d[0] + d[1] * d[1]) - radius) < 1E-9 && std::abs(d[2]) < 1E-9;
            const auto inward = vm::normalized(Vec { -d[0], -d[1], 0 });
            const auto tangent = vm::cross(Vec { 0, 0, 1 }, inward);
            const auto a = fieldAngles(p.dir, tangent, { 0, 0, 1 }, inward);
            const bool within = vm::dot(p.dir, inward) > 0 && std::abs(a[0]) <= 0.05 + 1E-9 && std::abs(a[1]) <= 0.1 + 1E-9;
            return onCircle && within;
        };
        if (mono) {
            IsotropicCircularMonoEnergyBeam<> beam(center);
            beam.setRadius(radius);
            beam.setEnergy(50);
            beam.setCollimationHalfAngles(angles);
            beam.setNumberOfExposures(2);
            beam.setNumberOfParticlesPerExposure(10000);
            checkBeamParticles(beam, name, 50, 50, geometry);
        } else {
            IsotropicCircularBeam<> beam(center);
            beam.setRadius(radius);
            beam.setEnergySpecter(Tube(80).getSpecter());
            beam.setCollimationHalfAngles(angles);
            beam.setNumberOfExposures(2);
            beam.setNumberOfParticlesPerExposure(10000);
            checkBeamParticles(beam, name, 0, 80, geometry);
        }
    }
}

void testDXBeam()
{
    DXBeam<> beam({ 0, 0, -100 });
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 });
    beam.setTubeVoltage(85);
    beam.setBeamSize(20, 10, 100);
    checkNear(beam.collimationHalfAngles()[0], std::atan(10.0 / 100.0), 1E-12, "DXBeam: setBeamSize sets x collimation");
    checkNear(beam.collimationHalfAngles()[1], std::atan(5.0 / 100.0), 1E-12, "DXBeam: setBeamSize sets y collimation");
    beam.setCollimationHalfAnglesDeg(5, 3);
    checkNear(beam.collimationHalfAnglesDeg()[0], 5, 1E-9, "DXBeam: collimation in degrees");
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(10000);
    const auto ax = beam.collimationHalfAngles()[0];
    const auto ay = beam.collimationHalfAngles()[1];
    checkBeamParticles(beam, "DXBeam", 0, 85, [&](std::size_t, const Particle& p) {
        const auto a = fieldAngles(p.dir, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 });
        return p.pos == Vec { 0, 0, -100 } && std::abs(a[0]) <= ax + 1E-9 && std::abs(a[1]) <= ay + 1E-9;
    });
    checkNear(beam.tube().voltage(), 85, 0, "DXBeam: tube voltage");
    beam.addTubeFiltrationMaterial(13, 2.5);
    checkNear(beam.tubeFiltration(13), 2.5, 1E-12, "DXBeam: tube filtration");
    check(beam.tubeAlHalfValueLayer() > 1 && beam.tubeMeanSpecterEnergy() > 20, "DXBeam: tube HVL and mean energy");

    // the DAP calibration scales with the DAP value and inversely with the number of particles
    beam.setDAPvalue(1.0);
    const auto c1 = beam.calibrationFactor();
    beam.setDAPvalue(3.0);
    checkRelative(beam.calibrationFactor(), 3 * c1, 1E-12, "DXBeam: calibration proportional to DAP");
    beam.setNumberOfParticlesPerExposure(20000);
    checkRelative(beam.calibrationFactor(), 1.5 * c1, 1E-12, "DXBeam: calibration inversely proportional to number of particles");
}

void testCBCTBeam()
{
    CBCTBeam<> beam;
    const Vec iso = { 1, 2, 3 };
    beam.setIsocenter(iso);
    beam.setRotationAxis({ 0, 0, 1 });
    beam.setTubeVoltage(100);
    beam.setStartAngleDeg(-90);
    beam.setStopAngleDeg(90);
    beam.setStepAngleDeg(1);
    beam.setSourceDetectorDistance(120);
    beam.setCollimationHalfAnglesDeg(10, 8);
    beam.setNumberOfParticlesPerExposure(1000);
    check(beam.numberOfExposures() == 180, "CBCTBeam: number of exposures for 180 deg arc with 1 deg step (" + std::to_string(beam.numberOfExposures()) + ")");
    const auto a = beam.collimationHalfAngles();
    checkBeamParticles(beam, "CBCTBeam", 0, 100, [&](std::size_t, const Particle& p) {
        const auto d = vm::subtract(iso, p.pos);
        const bool distance = std::abs(vm::length(d) - 60) < 1E-9 && std::abs(d[2]) < 1E-9;
        const auto inward = vm::normalized(d);
        const auto tangent = vm::cross(Vec { 0, 0, 1 }, inward);
        const auto ang = fieldAngles(p.dir, tangent, { 0, 0, 1 }, inward);
        return distance && std::abs(ang[0]) <= a[0] + 1E-9 && std::abs(ang[1]) <= a[1] + 1E-9;
    },
        500);
    // consecutive exposures are separated by the step angle around the isocenter
    const auto d0 = vm::subtract(beam.exposure(10).position(), iso);
    const auto d1 = vm::subtract(beam.exposure(11).position(), iso);
    checkNear(vm::angleBetween(d0, d1), DEG_TO_RAD(), 1E-9, "CBCTBeam: step angle between exposures");
    const auto dFirst = vm::subtract(beam.exposure(0).position(), iso);
    const auto dLast = vm::subtract(beam.exposure(beam.numberOfExposures() - 1).position(), iso);
    checkNear(vm::angleBetween(dFirst, dLast), DEG_TO_RAD() * 179, 1E-9, "CBCTBeam: arc between first and last exposure");
}

void testCTSpiralBeam()
{
    CTSpiralBeam<> beam;
    beam.setStartStopPosition({ 0, 0, -5 }, { 0, 0, 5 });
    beam.setCollimation(2);
    beam.setPitch(1.0);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);
    beam.setScanFieldOfView(50);
    beam.setTubeVoltage(120);
    beam.setNumberOfParticlesPerExposure(500);
    // 10 cm scan length with 2 cm table feed per rotation is 5 rotations with 72 exposures each
    check(beam.numberOfExposures() == 360, "CTSpiralBeam: number of exposures (" + std::to_string(beam.numberOfExposures()) + ")");
    const auto angles = beam.collimationHalfAngles();
    checkNear(angles[0], std::atan(50.0 / 115.0), 1E-12, "CTSpiralBeam: fan half angle from field of view");
    checkNear(angles[1], std::atan(0.5 * 2.0 / 115.0), 1E-12, "CTSpiralBeam: slice half angle from collimation");

    auto geometry = [&](std::size_t i, const Particle& p) {
        const auto z = p.pos[2];
        const Vec axisPoint = { 0, 0, z };
        const auto d = vm::subtract(axisPoint, p.pos);
        const bool distance = std::abs(vm::length(d) - 115.0 / 2) < 1E-9;
        const bool zpos = std::abs(z - (-5 + 2.0 * i * 5.0 / 360.0)) < 1E-9;
        const auto inward = vm::normalized(d);
        const auto tangent = vm::cross(Vec { 0, 0, 1 }, inward);
        const auto a = fieldAngles(p.dir, tangent, { 0, 0, 1 }, inward);
        return distance && zpos && std::abs(a[0]) <= angles[0] + 1E-9 && std::abs(a[1]) <= angles[1] + 1E-9;
    };
    checkBeamParticles(beam, "CTSpiralBeam", 0, 120, geometry, 200);

    const auto d0 = beam.exposure(0).position();
    const auto d72 = beam.exposure(72).position();
    checkNear(d72[2] - d0[2], 2.0, 1E-9, "CTSpiralBeam: table feed per rotation is pitch times collimation");
    checkNear(std::atan2(d72[1], d72[0]), std::atan2(d0[1], d0[0]), 1E-9, "CTSpiralBeam: same gantry angle after one rotation");

    // Exposure weight follows the AEC profile
    std::vector<double> aec(11);
    for (std::size_t i = 0; i < aec.size(); ++i)
        aec[i] = 1.0 + i;
    beam.setAECFilterData({ 0, 0, -5 }, { 0, 0, 5 }, aec);
    const auto wStart = beam.exposure(0).weight();
    const auto wEnd = beam.exposure(359).weight();
    check(wEnd > 5 * wStart, "CTSpiralBeam: exposure weight follows AEC profile");
}

void testCTSequentialBeam()
{
    CTSequentialBeam<> beam;
    beam.setPosition({ 0, 0, 0 });
    beam.setScanNormal({ 0, 0, 1 });
    beam.setNumberOfSlices(3);
    beam.setSliceSpacing(4);
    beam.setCollimation(1);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(100);
    beam.setTubeVoltage(100);
    beam.setNumberOfParticlesPerExposure(500);
    check(beam.numberOfExposures() == 3 * 72, "CTSequentialBeam: number of exposures for 3 slices with 5 deg step (" + std::to_string(beam.numberOfExposures()) + ")");

    // every exposure k * 72 + j belongs to slice k
    bool slicesOk = true;
    for (std::size_t i = 0; i < beam.numberOfExposures(); ++i) {
        const auto slice = i / 72;
        slicesOk = slicesOk && std::abs(beam.exposure(i).position()[2] - 4.0 * slice) < 1E-9;
    }
    check(slicesOk, "CTSequentialBeam: exposures are assigned to the correct slice");
    const auto angles = beam.collimationHalfAngles();
    checkBeamParticles(beam, "CTSequentialBeam", 0, 100, [&](std::size_t, const Particle& p) {
        const Vec axisPoint = { 0, 0, p.pos[2] };
        const auto d = vm::subtract(axisPoint, p.pos);
        const auto inward = vm::normalized(d);
        const auto tangent = vm::cross(Vec { 0, 0, 1 }, inward);
        const auto a = fieldAngles(p.dir, tangent, { 0, 0, 1 }, inward);
        return std::abs(vm::length(d) - 50) < 1E-9 && std::abs(a[0]) <= angles[0] + 1E-9 && std::abs(a[1]) <= angles[1] + 1E-9;
    },
        200);
}

void testCTDualEnergyBeam()
{
    CTSpiralDualEnergyBeam<> beam;
    beam.setStartStopPosition({ 0, 0, 0 }, { 0, 0, 4 });
    beam.setCollimation(2);
    beam.setPitch(1);
    beam.setStepAngleDeg(10);
    beam.setTubeAVoltage(140);
    beam.setTubeBVoltage(80);
    beam.setTubeBoffsetAngleDeg(90);
    beam.setScanFieldOfViewB(32);
    beam.setRelativeMasTubeA(2);
    beam.setRelativeMasTubeB(1);
    beam.setNumberOfParticlesPerExposure(500);
    check(beam.numberOfExposures() == 2 * 72, "CTSpiralDualEnergyBeam: exposures for both tubes (" + std::to_string(beam.numberOfExposures()) + ")");
    checkBeamParticles(beam, "CTSpiralDualEnergyBeam", 0, 140, {}, 200);

    // even exposures are tube A, odd exposures are tube B
    auto state = xraymctest::seededState();
    double maxB = 0;
    for (std::size_t i = 1; i < beam.numberOfExposures(); i += 2) {
        const auto e = beam.exposure(i);
        for (int j = 0; j < 200; ++j)
            maxB = std::max(maxB, e.sampleParticle(state).energy);
    }
    check(maxB <= 80, "CTSpiralDualEnergyBeam: tube B energies below tube B voltage");
    const auto pa = beam.exposure(0).position();
    const auto pb = beam.exposure(1).position();
    checkNear(vm::angleBetween(Vec { pa[0], pa[1], 0 }, Vec { pb[0], pb[1], 0 }), std::numbers::pi / 2, 1E-9, "CTSpiralDualEnergyBeam: tube B offset angle");
    // relative weights include tube output and field of view and are normalized to sum to two
    checkNear(beam.tubeRelativeWeightA() + beam.tubeRelativeWeightB(), 2, 1E-12, "CTSpiralDualEnergyBeam: relative tube weights sum to two");
    const auto ratio = beam.tubeRelativeWeightA() / beam.tubeRelativeWeightB();
    beam.setRelativeMasTubeB(2);
    checkRelative(beam.tubeRelativeWeightA() / beam.tubeRelativeWeightB(), ratio / 2, 1E-9, "CTSpiralDualEnergyBeam: relative weight proportional to relative mAs");
}

int main()
{
    std::cout << "Testing beams\n";
    testPencilBeam();
    testFlatField();
    testIsotropicBeams();
    testDXBeam();
    testCBCTBeam();
    testCTSpiralBeam();
    testCTSequentialBeam();
    testCTDualEnergyBeam();
    return xraymctest::summary("testbeams");
}

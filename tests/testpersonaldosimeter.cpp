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

#include "testutils.hpp"

#include "xraymc/beams/flatmonoenergyfield.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/personaldosimeter.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;
using xraymctest::checkStatisticallyEqual;
namespace vm = xraymc::vectormath;

void testGeometry()
{
    PersonalDosimeter<true> d;
    const auto n = d.normalVector();
    const auto& c = d.directionCosines();
    checkNear(vm::length(n), 1, 1E-12, "Dosimeter: normal is unit vector");
    checkNear(vm::dot(n, c[0]), 0, 1E-12, "Dosimeter: normal is orthogonal to x cosine");
    checkNear(vm::dot(n, c[1]), 0, 1E-12, "Dosimeter: normal is orthogonal to y cosine");
    const auto aabb = d.AABB();
    checkNear(aabb[3] - aabb[0], PersonalDosimeter<true>::width(), 1E-12, "Dosimeter: AABB width");
    checkNear(aabb[4] - aabb[1], PersonalDosimeter<true>::height(), 1E-12, "Dosimeter: AABB height");

    d.rotate(std::numbers::pi / 3, { 1, 0, 0 });
    const auto n2 = d.normalVector();
    checkNear(vm::angleBetween(n, n2), std::numbers::pi / 3, 1E-9, "Dosimeter: rotate turns normal");
    checkNear(vm::dot(d.directionCosines()[0], d.directionCosines()[1]), 0, 1E-12, "Dosimeter: rotated cosines stay orthogonal");

    d.setCenter({ 1, 2, 3 });
    check(xraymctest::checkVectorNear(d.center(), { 1, 2, 3 }, 1E-12, "Dosimeter: setCenter"), "Dosimeter: center");

    // Only the front face is hit
    PersonalDosimeter<true> f;
    const Particle front { .pos = { 0, 0, 5 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    const Particle side { .pos = { 3, 0, 5 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
    check(f.intersect(front).valid(), "Dosimeter: ray towards front face hits");
    checkNear(f.intersect(front).intersection, 5, 1E-9, "Dosimeter: distance to face");
    check(!f.intersect(side).valid(), "Dosimeter: ray beside dosimeter misses");
}

void testResponse()
{
    PersonalDosimeter<true> d;
    const auto n = d.normalVector();
    const auto cx = d.directionCosines()[0];
    auto incidence = [&](double angle) {
        // direction towards the face, rotated by angle about the y cosine
        const auto dir = vm::rotate(vm::scale(n, -1.0), d.directionCosines()[1], angle);
        return Particle { .pos = { 0, 0, 0 }, .dir = dir, .energy = 60, .weight = 1 };
    };
    const auto w0 = d.angularResponseWeight(incidence(0));
    checkNear(w0, std::pow(std::cos(15.0 * std::numbers::pi / 180), 0.2), 1E-9, "Dosimeter: response at normal incidence");
    check(d.angularResponseWeight(Particle { .pos = { 0, 0, 0 }, .dir = n, .energy = 60, .weight = 1 }) == 0, "Dosimeter: no response from behind");
    bool decreasing = true, symmetric = true, bounded = true;
    double last = w0;
    for (double a = 0.1; a < 1.5; a += 0.1) {
        const auto w = d.angularResponseWeight(incidence(a));
        decreasing = decreasing && w < last;
        symmetric = symmetric && std::abs(w - d.angularResponseWeight(incidence(-a))) < 1E-12;
        bounded = bounded && w > 0 && w <= 1;
        last = w;
    }
    check(decreasing, "Dosimeter: response decreases with incidence angle");
    check(symmetric, "Dosimeter: response symmetric in incidence angle");
    check(bounded, "Dosimeter: response in (0, 1]");
    (void)cx;

    // Hp(10) conversion coefficients (Solc et al. 2025), interpolated by a compact spline
    checkRelative(d.hp10ConvertionFactor(20), 0.6087834, 0.01, "Dosimeter: Hp(10)/Ka at 20 keV");
    checkRelative(d.hp10ConvertionFactor(30), 1.1195263, 0.01, "Dosimeter: Hp(10)/Ka at 30 keV");
    checkRelative(d.hp10ConvertionFactor(60), 1.9001095, 0.01, "Dosimeter: Hp(10)/Ka at 60 keV");
    checkRelative(d.hp10ConvertionFactor(80), 1.9075632, 0.01, "Dosimeter: Hp(10)/Ka at 80 keV");
    checkNear(d.hp10ConvertionFactor(5), 0, 0, "Dosimeter: Hp(10)/Ka is zero below 7.5 keV");
    bool nonNegative = true;
    for (double e = 7.5; e < 150; e += 0.25)
        nonNegative = nonNegative && d.hp10ConvertionFactor(e) >= 0;
    check(nonNegative, "Dosimeter: Hp(10)/Ka is never negative");
}

// Air KERMA and Hp(10) scored in a calibrated flat field
template <bool DIRECTIONAL>
void testCalibratedField(double angle)
{
    using D = PersonalDosimeter<DIRECTIONAL>;
    World<D> world;
    auto& d = world.template addItem<D>();
    d.rotate(angle, d.directionCosines()[1]);
    world.setMaterialDensity(1E-9);
    world.build();

    // Field perpendicular to the beam direction, larger than the dosimeter
    const auto n = PersonalDosimeter<DIRECTIONAL>().normalVector();
    FlatMonoEnergyField<> beam(vm::scale(n, 10.0), 60);
    beam.setDirectionCosines({ 0, 1, 0 }, { 1, 0, 0 }); // beam travels along -z
    beam.setLenght(8, 8);
    beam.setAirKerma(1.0);
    beam.setNumberOfExposures(20);
    beam.setNumberOfParticlesPerExposure(50000);
    Transport transport;
    transport(world, beam);

    const auto label = std::string("Dosimeter") + (DIRECTIONAL ? " (directional)" : "") + " at " + std::to_string(int(std::round(angle * 180 / std::numbers::pi))) + " deg";
    const auto& ka = d.airKermaScored();
    const auto& hp = d.hp10DoseScored();
    // Without directional dependence every crossing photon counts fully, a tilted dosimeter
    // intercepts cos(angle) of the fluence. The directional model cancels this area effect
    // and applies the angular response instead.
    double expected = std::cos(angle);
    if constexpr (DIRECTIONAL) {
        const Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, -1 }, .energy = 60, .weight = 1 };
        expected = d.angularResponseWeight(p);
    }
    checkStatisticallyEqual(ka.dose(), ka.standardDeviation(), expected, 0, 4, label + ": scored air KERMA equals field air KERMA times response");
    checkRelative(hp.dose(), ka.dose() * d.hp10ConvertionFactor(60), 1E-9, label + ": Hp(10) equals air KERMA times conversion coefficient");
}

int main()
{
    std::cout << "Testing personal dosimeter\n";
    testGeometry();
    testResponse();
    testCalibratedField<false>(0);
    testCalibratedField<false>(std::numbers::pi / 4);
    testCalibratedField<true>(0);
    testCalibratedField<true>(std::numbers::pi / 4);
    return xraymctest::summary("testpersonaldosimeter");
}

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

#include "testutils.hpp"

#include "xraymc/beams/filters/bowtiefilter.hpp"
#include "xraymc/beams/filters/ctaecfilter.hpp"
#include "xraymc/beams/filters/ctorganaecfilter.hpp"
#include "xraymc/constants.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

constexpr double pi = std::numbers::pi_v<double>;

// Mean value of f over [a, b] by the midpoint rule
template <typename F>
double meanValue(const F& f, double a, double b, std::size_t N = 200000)
{
    double sum = 0;
    const double h = (b - a) / N;
    for (std::size_t i = 0; i < N; ++i)
        sum += f(a + (i + 0.5) * h);
    return sum / N;
}

void testBowtie()
{
    // Measured bowtie profile, angle [rad] and relative intensity
    const std::vector<double> weight = { 9084.336977, 9056.683016, 9010.555199, 8976.621955, 8937.858092, 8847.304732, 8754.532826, 8656.860432, 8527.312593, 8422.551773, 8306.212116, 8125.62487, 7931.430865, 7786.041916, 7639.243664, 7424.480585, 7197.561479, 7024.472266, 6851.343852, 6621.742936, 6388.017234, 6152.921949, 5974.012763, 5794.802125, 5568.133743, 5353.997039, 5134.814809, 4912.471357, 4695.75262, 4498.066906, 4310.238528, 4119.691564, 3935.618728, 3777.102292, 3629.564305, 3476.787516, 3289.33132, 3128.526529, 3013.761207, 2897.170244, 2762.216686, 2654.233847, 2548.605629, 2441.276609, 2178.949877, 2116.337359, 2015.555302, 1964.260525, 1900.254578, 1857.235929, 1790.316289, 1756.816535, 1691.007325, 1654.787224, 1621.951184, 1579.044149, 1547.725401, 1520.275338, 1492.313526, 1473.686923, 1446.729063, 1427.908945, 1399.920027, 1248.556091 };
    std::vector<double> angle;
    for (std::size_t i = 0; i < weight.size(); ++i)
        angle.push_back(0.002678189 + i * 0.005356377);

    xraymc::BowtieFilter filter(angle, weight);
    const auto a0 = angle.front(), a1 = angle.back();
    checkNear(meanValue([&](double a) { return filter(a); }, a0, a1), 1.0, 1E-3, "Bowtie: mean weight over measured angles is one");
    check(filter(0.1) == filter(-0.1), "Bowtie: one sided filter is symmetric");
    check(filter(a0) > filter(0.15) && filter(0.15) > filter(0.3), "Bowtie: weight decreases away from center");
    // The filter is a compact 6 segment spline through the measured profile, end points are exact
    checkRelative(filter(angle[10]) / filter(angle[40]), weight[10] / weight[40], 0.02, "Bowtie: relative weights approximate measured relative intensity");
    checkRelative(filter(a1) / filter(a0), weight.back() / weight.front(), 1E-9, "Bowtie: relative weight at last knot equals measured relative intensity");

    xraymc::BowtieFilter def;
    check(def(0) > def(0.2) && def(0.2) > 0, "Bowtie: default filter is peaked at center");
    checkNear(meanValue([&](double a) { return def(a); }, 0, 0.390607044), 1.0, 1E-3, "Bowtie: default filter mean weight is one");

    const auto buffer = filter.serialize();
    const auto copy = xraymc::BowtieFilter::deserialize(buffer);
    check(copy.has_value() && copy.value()(0.123) == filter(0.123), "Bowtie: serialize round trip");
}

void testOrganAEC()
{
    xraymc::CTOrganAECFilter f;
    const auto deg45 = pi / 4;
    const auto deg10 = xraymc::DEG_TO_RAD() * 10;
    f.setLowWeightFactor(0.25);
    f.setStartStopAngles(-deg45, deg45);
    f.setRampAngle(deg10);

    checkNear(f(0), 0.25, 1E-12, "Organ AEC: low weight inside zone");
    checkNear(f(pi), f.maxWeight(), 1E-12, "Organ AEC: max weight opposite to zone");
    check(f.maxWeight() > 1, "Organ AEC: compensating filter has max weight above one");
    checkNear(f(deg45 + deg10 / 2), (0.25 + f.maxWeight()) / 2, 1E-9, "Organ AEC: linear ramp midpoint");
    checkNear(f(deg45 + 1E-9), 0.25, 1E-6, "Organ AEC: ramp is continuous at low zone edge");
    checkNear(f(deg45 + deg10 - 1E-9), f.maxWeight(), 1E-6, "Organ AEC: ramp is continuous at high zone edge");
    checkNear(f(-deg45 - deg10 / 2), (0.25 + f.maxWeight()) / 2, 1E-9, "Organ AEC: ramp on both sides");
    checkNear(f(2 * pi), f(0), 1E-12, "Organ AEC: angles are periodic");
    checkNear(meanValue([&](double a) { return f(a); }, -pi, pi), 1.0, 1E-4, "Organ AEC: mean weight over a rotation is one");

    f.setCompensateOutside(false);
    checkNear(f.maxWeight(), 1.0, 0, "Organ AEC: non compensating filter has max weight one");
    checkNear(f(pi), 1.0, 1E-12, "Organ AEC: non compensating filter outside zone");

    f.setLowWeightFactor(5);
    checkNear(f.lowWeight(), 1.0, 0, "Organ AEC: low weight factor is clamped to one");
    f.setRampAngleDeg(80);
    checkNear(f.rampAngle(), pi / 4, 1E-12, "Organ AEC: ramp angle is clamped to 45 deg");
    f.setStartAngleDeg(-30);
    checkNear(f.startAngleDeg(), -30, 1E-9, "Organ AEC: start angle in degrees");

    xraymc::CTOrganAECFilter g(0.4, -1.0, 0.5, 0.2);
    g.setUseFilter(true);
    const auto buffer = g.serialize();
    const auto copy = xraymc::CTOrganAECFilter::deserialize(buffer);
    // compare stored state exactly; evaluated weights may differ by an ulp under -ffast-math
    check(copy.has_value()
            && copy->startAngle() == g.startAngle()
            && copy->stopAngle() == g.stopAngle()
            && copy->rampAngle() == g.rampAngle()
            && copy->lowWeight() == g.lowWeight()
            && copy->compensateOutside() == g.compensateOutside()
            && copy->useFilter(),
        "Organ AEC: serialize round trip");
    checkNear(copy.value()(2.0), g(2.0), 1E-12, "Organ AEC: deserialized filter gives same weight");
}

void testCTAEC()
{
    const std::array<double, 3> start = { 0, 0, -5 };
    const std::array<double, 3> stop = { 0, 0, 4 };

    xraymc::CTAECFilter def;
    check(def.isEmpty(), "CT AEC: default filter is empty");
    checkNear(def({ 0, 0, 0.3 }), 1.0, 1E-12, "CT AEC: default filter weight is one");

    // A constant profile must give weight one everywhere
    xraymc::CTAECFilter c(start, stop, std::vector<double>(16, 2.0));
    bool flat = true;
    for (double z = -5; z <= 4; z += 0.1)
        flat = flat && std::abs(c({ 0, 0, z }) - 1) < 1E-9;
    check(flat, "CT AEC: constant profile gives unit weight along the whole axis");
    check(!c.isEmpty() && c.size() == 16, "CT AEC: filter size");
    checkNear(c.length(), 9, 1E-12, "CT AEC: axis length");

    std::vector<double> w = { 5, 6, 7, 6, 3, 5, 6, 2, 6, 5, 6, 7, 8, 9, 9, 9 };
    xraymc::CTAECFilter f(start, stop, w);
    checkNear(f.integrate(), f.length(), 1E-9, "CT AEC: integral equals axis length");
    checkNear(meanValue([&](double z) { return f({ 0, 0, z }); }, -5, 4), 1.0, 1E-6, "CT AEC: mean weight along the axis is one");
    // relative weights at the sample positions are preserved
    const double step = 9.0 / 15.0;
    checkRelative(f({ 0, 0, -5 + step * 4 }) / f({ 0, 0, -5 + step * 12 }), w[4] / w[12], 1E-9, "CT AEC: relative weights at sample positions preserved");
    checkRelative(f({ 0, 0, 4 }) / f({ 0, 0, -5 }), w.back() / w.front(), 1E-9, "CT AEC: last sample is at the stop position");
    checkNear(f({ 0, 0, -50 }), f({ 0, 0, -5 }), 1E-12, "CT AEC: positions before start are clamped");
    checkNear(f({ 0, 0, 50 }), f({ 0, 0, 4 }), 1E-12, "CT AEC: positions after stop are clamped");
    checkNear(f({ 3, -2, 0 }), f({ 0, 0, 0 }), 1E-12, "CT AEC: weight depends on axial position only");

    xraymc::CTAECFilter g(start, stop, std::vector<double>(100, 2.0));
    g.normalizeBetween({ 0, 0, -1 }, { 0, 0, 1 });
    checkNear(g.integrate({ 0, 0, -1 }, { 0, 0, 1 }), 2.0, 1E-9, "CT AEC: normalizeBetween gives integral equal to sub-range length");
    checkNear(meanValue([&](double z) { return g({ 0, 0, z }); }, -1, 1), 1.0, 1E-6, "CT AEC: normalizeBetween gives unit mean in sub-range");

    const auto buffer = f.serialize();
    const auto copy = xraymc::CTAECFilter::deserialize(buffer);
    check(copy.has_value() && copy.value()({ 0, 0, 1.234 }) == f({ 0, 0, 1.234 }) && copy.value().length() == f.length(), "CT AEC: serialize round trip");
}

int main()
{
    std::cout << "Testing beam filters\n";
    testBowtie();
    testOrganAEC();
    testCTAEC();
    return xraymctest::summary("testbeamfilters");
}

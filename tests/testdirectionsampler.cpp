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

#include "xraymc/beams/utilities/spheresamplingcircularfield.hpp"
#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using xraymctest::check;
using xraymctest::checkNear;

// Solid angle of the rectangle [x1, x2] x [y1, y2] in the plane z = 1 seen from the origin
double rectangleSolidAngle(double x1, double y1, double x2, double y2)
{
    auto F = [](double x, double y) { return std::atan(x * y / std::sqrt(1 + x * x + y * y)); };
    return F(x2, y2) - F(x1, y2) - F(x2, y1) + F(x1, y1);
}

// Samples directions and tests that they fall within the angular field and that the
// fraction of directions within a grid of sub-fields equals the sub-field solid angle.
template <typename Sampler>
void testRectangular(const Sampler& sampler, double xmin, double ymin, double xmax, double ymax, const std::string& name)
{
    auto state = xraymctest::seededState();
    constexpr std::size_t N = 2000000;
    constexpr std::size_t G = 6;
    const double tx0 = std::tan(xmin), tx1 = std::tan(xmax), ty0 = std::tan(ymin), ty1 = std::tan(ymax);
    std::vector<double> counts(G * G, 0);
    bool unit = true;
    bool inside = true;
    for (std::size_t i = 0; i < N; ++i) {
        const auto d = sampler(state);
        unit = unit && std::abs(xraymc::vectormath::length(d) - 1) < 1E-9;
        const auto px = d[0] / d[2];
        const auto py = d[1] / d[2];
        const double eps = 1E-9;
        inside = inside && d[2] > 0 && px >= tx0 - eps && px <= tx1 + eps && py >= ty0 - eps && py <= ty1 + eps;
        const auto ix = std::min(static_cast<std::size_t>(std::max((px - tx0) / (tx1 - tx0), 0.0) * G), G - 1);
        const auto iy = std::min(static_cast<std::size_t>(std::max((py - ty0) / (ty1 - ty0), 0.0) * G), G - 1);
        counts[ix + iy * G] += 1;
    }
    check(unit, name + ": sampled directions are unit vectors");
    check(inside, name + ": sampled directions within collimation angles");

    const auto total = rectangleSolidAngle(tx0, ty0, tx1, ty1);
    double chi2 = 0;
    for (std::size_t iy = 0; iy < G; ++iy)
        for (std::size_t ix = 0; ix < G; ++ix) {
            const double x1 = tx0 + (tx1 - tx0) * ix / G;
            const double x2 = tx0 + (tx1 - tx0) * (ix + 1) / G;
            const double y1 = ty0 + (ty1 - ty0) * iy / G;
            const double y2 = ty0 + (ty1 - ty0) * (iy + 1) / G;
            const auto e = N * rectangleSolidAngle(x1, y1, x2, y2) / total;
            chi2 += (counts[ix + iy * G] - e) * (counts[ix + iy * G] - e) / e;
        }
    const double dof = G * G - 1;
    const auto z = (chi2 - dof) / std::sqrt(2 * dof);
    check(z < 5, name + ": directions are uniform in solid angle (chi-square z = " + std::to_string(z) + ")");
}

// The square field sampler bounds the direction components |x|, |y| <= sin(half_angle).
// Uniform in solid angle means the density in (x, y) is proportional to 1 / z.
void testSquare(double halfAngle)
{
    xraymc::SphereSamplingSquareField sampler(halfAngle);
    auto state = xraymctest::seededState();
    const auto b = std::sin(halfAngle);
    constexpr std::size_t N = 1000000;
    constexpr std::size_t G = 5;
    std::vector<double> counts(G * G, 0);
    bool inside = true;
    for (std::size_t i = 0; i < N; ++i) {
        const auto d = sampler(state);
        inside = inside && std::abs(d[0]) <= b && std::abs(d[1]) <= b && d[2] > 0 && std::abs(xraymc::vectormath::length(d) - 1) < 1E-9;
        const auto ix = std::min(static_cast<std::size_t>((d[0] + b) / (2 * b) * G), G - 1);
        const auto iy = std::min(static_cast<std::size_t>((d[1] + b) / (2 * b) * G), G - 1);
        counts[ix + iy * G] += 1;
    }
    const auto name = "SphereSamplingSquareField (" + std::to_string(halfAngle) + " rad)";
    check(inside, name + ": directions within |x|, |y| <= sin(half angle)");

    // cell probabilities by midpoint integration of 1 / sqrt(1 - x^2 - y^2)
    std::vector<double> probs(G * G, 0);
    constexpr std::size_t M = 60;
    const double h = 2 * b / (G * M);
    for (std::size_t iy = 0; iy < G * M; ++iy)
        for (std::size_t ix = 0; ix < G * M; ++ix) {
            const double x = -b + (ix + 0.5) * h;
            const double y = -b + (iy + 0.5) * h;
            probs[ix / M + (iy / M) * G] += 1 / std::sqrt(1 - x * x - y * y);
        }
    const auto total = std::reduce(probs.cbegin(), probs.cend(), 0.0);
    double chi2 = 0;
    for (std::size_t i = 0; i < G * G; ++i) {
        const auto e = N * probs[i] / total;
        chi2 += (counts[i] - e) * (counts[i] - e) / e;
    }
    const double dof = G * G - 1;
    const auto z = (chi2 - dof) / std::sqrt(2 * dof);
    check(z < 5, name + ": directions are uniform in solid angle (chi-square z = " + std::to_string(z) + ")");
}

void testCircular()
{
    for (double angle : { 0.05, 0.5, 1.2 }) {
        xraymc::SphereSamplingCircularField sampler(angle);
        auto state = xraymctest::seededState();
        constexpr std::size_t N = 1000000;
        bool inside = true;
        double sumCos = 0;
        double sumX = 0, sumY = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const auto d = sampler(state);
            inside = inside && d[2] >= std::cos(angle) - 1E-12 && std::abs(xraymc::vectormath::length(d) - 1) < 1E-9;
            sumCos += d[2];
            sumX += d[0];
            sumY += d[1];
        }
        const auto name = "SphereSamplingCircularField (" + std::to_string(angle) + " rad)";
        check(inside, name + ": directions within cone");
        // uniform in solid angle means cos(theta) is uniform in [cos(angle), 1]
        const auto c0 = std::cos(angle);
        const auto sd = (1 - c0) / std::sqrt(12.0 * N);
        checkNear(sumCos / N, (1 + c0) / 2, 5 * sd, name + ": mean cos(theta) of uniform cone");
        const auto sx = std::sin(angle) / std::sqrt(2.0 * N);
        checkNear(sumX / N, 0, 5 * sx, name + ": symmetric in x");
        checkNear(sumY / N, 0, 5 * sx, name + ": symmetric in y");
    }
}

int main()
{
    std::cout << "Testing direction samplers\n";
    {
        const double ax = 0.3, ay = 0.1;
        xraymc::SphereSamplingRectangularField s(ax, ay);
        testRectangular(s, -ax, -ay, ax, ay, "SphereSamplingRectangularField symmetric");
    }
    {
        xraymc::SphereSamplingRectangularField s(-0.1, -0.4, 0.5, 0.2);
        testRectangular(s, -0.1, -0.4, 0.5, 0.2, "SphereSamplingRectangularField asymmetric");
    }
    {
        // field not containing the central axis
        const double alpha = 15 * xraymc::DEG_TO_RAD();
        const auto h = 180 * std::tan(alpha);
        const auto ymin = std::atan((h - 39.0 / 2) / 180.0);
        const auto ymax = std::atan((h + 39.0 / 2) / 180.0);
        const auto ax = std::atan(39.0 / (2 * 180.0));
        xraymc::SphereSamplingRectangularField s;
        s.setData(-ax, ymin, ax, ymax);
        testRectangular(s, -ax, ymin, ax, ymax, "SphereSamplingRectangularField off-axis field");
    }
    {
        xraymc::SphereSamplingRectangularField s(std::array<double, 2> { 0.7, 0.7 });
        testRectangular(s, -0.7, -0.7, 0.7, 0.7, "SphereSamplingRectangularField wide field");
    }
    {
        xraymc::SphereSamplingRectangularGeneralField s(-0.2, -0.3, 0.4, 0.1);
        testRectangular(s, -0.2, -0.3, 0.4, 0.1, "SphereSamplingRectangularGeneralField asymmetric");
    }
    testSquare(0.05);
    testSquare(0.6);
    testCircular();
    return xraymctest::summary("testdirectionsampler");
}

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

#include "xraymc/interpolation.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

void testLinearInterp()
{
    checkNear(xraymc::interp(1.0, 3.0, 10.0, 20.0, 2.0), 15.0, 1E-12, "interp scalar midpoint");
    checkNear(xraymc::interp(1.0, 3.0, 10.0, 20.0, 4.0), 25.0, 1E-12, "interp scalar extrapolates");
    double x[2] = { 0, 2 };
    double y[2] = { -1, 1 };
    checkNear(xraymc::interp(x, y, 0.5), -0.5, 1E-12, "interp array overload");
    checkNear(xraymc::interp(std::make_pair(0.0, 1.0), std::make_pair(1.0, 3.0), 0.25), 1.5, 1E-12, "interp pair overload");
}

void testLogLogInterp()
{
    // A power law is a straight line in log-log space and must be reproduced exactly
    auto f = [](double x) { return 3.0 * std::pow(x, -2.5); };
    const double x0 = 2, x1 = 50;
    bool exact = true;
    for (double x = x0; x <= x1; x += 1.7)
        exact = exact && std::abs(xraymc::logloginterp(x0, x1, f(x0), f(x1), x) / f(x) - 1) < 1E-10;
    check(exact, "logloginterp reproduces power law");
    double xa[2] = { x0, x1 };
    double ya[2] = { f(x0), f(x1) };
    checkRelative(xraymc::logloginterp(xa, ya, 7.0), f(7.0), 1E-10, "logloginterp array overload");
    check(std::isnan(xraymc::logloginterp(1.0, 2.0, -1.0, 2.0, 1.5)), "logloginterp with negative value gives NaN");
}

void testInterpolateTables()
{
    const std::vector<std::pair<double, double>> data = { { 0, 0 }, { 1, 2 }, { 2, 3 }, { 4, 1 } };
    checkNear(xraymc::interpolate(data, 0.5), 1.0, 1E-12, "interpolate table first segment");
    checkNear(xraymc::interpolate(data, 1.5), 2.5, 1E-12, "interpolate table middle segment");
    checkNear(xraymc::interpolate(data, 3.0), 2.0, 1E-12, "interpolate table last segment");
    checkNear(xraymc::interpolate(data, 1.0), 2.0, 1E-12, "interpolate table at knot");
    checkNear(xraymc::interpolate(data, -1.0), -2.0, 1E-12, "interpolate table extrapolates below");
    checkNear(xraymc::interpolate(data, 5.0), 0.0, 1E-12, "interpolate table extrapolates above");
    checkNear(xraymc::interpolate<double, false, true>(data, -1.0), 0.0, 0, "interpolate without downward extrapolation returns 0");
    checkNear(xraymc::interpolate<double, true, false>(data, 5.0), 0.0, 0, "interpolate without upward extrapolation returns 0");

    const std::vector<double> xs = { 0.1, 0.5, 1.2, 2.2, 3.9 };
    const auto ys = xraymc::interpolate(data, xs);
    bool same = ys.size() == xs.size();
    for (std::size_t i = 0; i < xs.size() && same; ++i)
        same = std::abs(ys[i] - xraymc::interpolate(data, xs[i])) < 1E-12;
    check(same, "vector interpolate equals scalar interpolate");

    const std::vector<double> xv = { 0, 1, 2, 4 };
    const std::vector<double> yv = { 0, 2, 3, 1 };
    checkNear(xraymc::interpolate(xv.cbegin(), xv.cend(), yv.cbegin(), yv.cend(), 3.0), 2.0, 1E-12, "iterator interpolate");
    checkNear(xraymc::interpolate(xv.cbegin(), xv.cend(), yv.cbegin(), yv.cend(), -3.0), 0.0, 1E-12, "iterator interpolate clamps below");
    checkNear(xraymc::interpolate(xv.cbegin(), xv.cend(), yv.cbegin(), yv.cend(), 9.0), 1.0, 1E-12, "iterator interpolate clamps above");
}

void testRemovePoints()
{
    std::vector<std::pair<double, double>> line;
    for (int i = 0; i <= 10; ++i)
        line.push_back({ i, 2.0 * i + 1 });
    xraymc::removeUnneededInterpolationPoints(line);
    check(line.size() == 2 && line.front().first == 0 && line.back().first == 10, "collinear points reduced to end points");

    std::vector<std::pair<double, double>> kink = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 1 } };
    xraymc::removeUnneededInterpolationPoints(kink);
    check(kink.size() == 3 && kink[1].first == 2, "kink point is kept");
}

void testVectorOps()
{
    const std::vector<double> a = { 1, 2, 3 };
    const std::vector<double> b = { 10, 20, 30 };
    const std::vector<double> c = { 100, 200, 300 };
    check(xraymc::addVectors(a, b) == std::vector<double> { 11, 22, 33 }, "addVectors of two vectors");
    check(xraymc::addVectors(a, b, c) == std::vector<double> { 111, 222, 333 }, "addVectors of three vectors");
}

void testTrapz()
{
    // trapezoidal rule is exact for piecewise linear functions
    std::vector<double> x, f;
    std::vector<std::pair<double, double>> fp;
    for (int i = 0; i <= 20; ++i) {
        const double xi = i * 0.5;
        x.push_back(xi);
        f.push_back(3 * xi + 1);
        fp.push_back({ xi, 3 * xi + 1 });
    }
    const auto analytic = [](double a, double b) { return 1.5 * (b * b - a * a) + (b - a); };
    checkNear(xraymc::trapz(x, f), analytic(0, 10), 1E-10, "trapz of linear function is exact");
    checkNear(xraymc::trapz(fp), analytic(0, 10), 1E-10, "trapz of pairs is exact");
    const auto cum = xraymc::trapz_cum(x, f);
    check(cum.size() == x.size() && cum.front() == 0, "trapz_cum starts at zero");
    checkNear(cum.back(), analytic(0, 10), 1E-10, "trapz_cum ends at total integral");
    checkNear(cum[6], analytic(0, 3), 1E-10, "trapz_cum intermediate value");
    checkNear(xraymc::trapz(fp, 1.25, 7.75), analytic(1.25, 7.75), 1E-10, "trapz over sub interval between knots");
    checkNear(xraymc::trapz(fp, -5.0, 50.0), analytic(0, 10), 1E-10, "trapz sub interval is clamped to data");

    // quadratic function converges with second order accuracy
    std::vector<double> xq, fq;
    for (int i = 0; i <= 1000; ++i) {
        xq.push_back(i / 1000.0);
        fq.push_back(xq.back() * xq.back());
    }
    checkNear(xraymc::trapz(xq, fq), 1.0 / 3.0, 1E-6, "trapz of x^2 on [0, 1]");
}

void testGauss()
{
    // 20 point Gauss-Legendre integrates polynomials up to degree 39 exactly
    auto poly = [](double x) { return 5 * std::pow(x, 9) - 2 * std::pow(x, 4) + x - 7; };
    auto polyInt = [](double x) { return 0.5 * std::pow(x, 10) - 0.4 * std::pow(x, 5) + 0.5 * x * x - 7 * x; };
    checkRelative(xraymc::gaussIntegration(-1.5, 2.0, poly), polyInt(2.0) - polyInt(-1.5), 1E-12, "gaussIntegration exact for polynomial");
    checkNear(xraymc::gaussIntegration(0.0, std::numbers::pi, [](double x) { return std::sin(x); }), 2.0, 1E-12, "gaussIntegration of sin on [0, pi]");
    checkNear(xraymc::gaussIntegration(0.0, 1.0, [](double x) { return std::exp(x); }), std::numbers::e - 1, 1E-12, "gaussIntegration of exp on [0, 1]");

    const auto pts = xraymc::gaussIntegrationPoints(0.0, 2.0);
    bool inside = true;
    for (auto p : pts)
        inside = inside && p > 0 && p < 2;
    check(inside, "gauss nodes are inside the integration interval");
}

void testAkima()
{
    // Akima splines reproduce straight lines exactly and pass through every knot
    std::vector<std::pair<double, double>> line;
    for (int i = 0; i < 10; ++i)
        line.push_back({ i * 1.5, -2.0 * i * 1.5 + 4 });
    xraymc::AkimaSpline<double> sl(line);
    bool linear = true;
    for (double x = 0; x <= 13.5; x += 0.37)
        linear = linear && std::abs(sl(x) - (-2.0 * x + 4)) < 1E-10;
    check(linear, "AkimaSpline reproduces linear data");
    checkNear(sl.integral(1.0, 9.0), (-(81.0 - 1.0) + 4 * 8.0), 1E-9, "AkimaSpline integral of linear data");

    std::vector<std::pair<double, double>> data;
    for (int i = 0; i < 15; ++i) {
        const double x = i * 0.3;
        data.push_back({ x, std::sin(x) });
    }
    // knots given out of order must be sorted by the spline
    std::swap(data[3], data[10]);
    xraymc::AkimaSpline<double> sp(data);
    bool knots = true;
    for (const auto& [x, y] : data)
        knots = knots && std::abs(sp(x) - y) < 1E-10;
    check(knots, "AkimaSpline passes through all knots");
    // Akima uses one sided slopes at the end points, the first and last interval is less accurate
    double maxErr = 0, maxErrInterior = 0;
    for (double x = 0; x <= 4.2; x += 0.01) {
        const auto err = std::abs(sp(x) - std::sin(x));
        maxErr = std::max(maxErr, err);
        if (x > 0.3 && x < 3.9)
            maxErrInterior = std::max(maxErrInterior, err);
    }
    check(maxErrInterior < 2E-3, "AkimaSpline approximates sin in interior intervals (max err " + std::to_string(maxErrInterior) + ")");
    check(maxErr < 1E-2, "AkimaSpline approximates sin in end intervals (max err " + std::to_string(maxErr) + ")");
    checkNear(sp(4.2), std::sin(4.2), 1E-10, "AkimaSpline evaluates last knot from last segment");
    checkNear(sp(4.05), sp.integral(4.05 - 1E-4, 4.05 + 1E-4) / 2E-4, 1E-6, "AkimaSpline value in last interval consistent with integral");
    checkNear(sp.integral(0.0, 4.2), 1 - std::cos(4.2), 2E-3, "AkimaSpline integral approximates analytical integral");

    auto scaled = sp;
    scaled.scale(3.0);
    checkNear(scaled(1.234), 3 * sp(1.234), 1E-12, "AkimaSpline scale multiplies values");

    std::vector<double> xs, ys;
    for (const auto& [x, y] : data) {
        xs.push_back(x);
        ys.push_back(y);
    }
    xraymc::AkimaSpline<double> spxy(xs, ys);
    checkNear(spxy(2.05), sp(2.05), 1E-12, "AkimaSpline x/y constructor equals pair constructor");

    // fewer than five points are padded internally
    xraymc::AkimaSpline<double> few({ { 0, 1 }, { 1, -1 } });
    checkNear(few(0.25), 0.5, 1E-10, "AkimaSpline with two points is linear");
    checkNear(few(0.0), 1.0, 1E-10, "AkimaSpline with two points passes through first knot");
    checkNear(few(1.0), -1.0, 1E-10, "AkimaSpline with two points passes through last knot");
}

void testAkimaStatic()
{
    std::vector<std::pair<double, double>> line;
    for (int i = 0; i < 30; ++i)
        line.push_back({ i * 0.1, 0.5 * i * 0.1 + 2 });
    xraymc::AkimaSplineStatic<double, 8> sl(line);
    bool linear = true;
    for (double x = 0; x <= 2.9; x += 0.05)
        linear = linear && std::abs(sl(x) - (0.5 * x + 2)) < 1E-10;
    check(linear, "AkimaSplineStatic reproduces linear data");
    checkNear(sl.integral(0.0, 2.0), 0.25 * 4 + 4, 1E-9, "AkimaSplineStatic integral of linear data");

    std::vector<std::pair<double, double>> curve;
    for (int i = 0; i <= 100; ++i)
        curve.push_back({ i * 0.03, std::exp(-i * 0.03) });
    xraymc::AkimaSplineStatic<double, 12> sc(curve);
    double maxErr = 0;
    for (double x = 0; x <= 3; x += 0.01)
        maxErr = std::max(maxErr, std::abs(sc(x) - std::exp(-x)));
    check(maxErr < 1E-2, "AkimaSplineStatic approximates exp(-x) (max err " + std::to_string(maxErr) + ")");
    checkNear(sc(3.0), std::exp(-3.0), 1E-10, "AkimaSplineStatic evaluates last knot from last segment");
    checkNear(sc(2.9), sc.integral(2.9 - 1E-4, 2.9 + 1E-4) / 2E-4, 1E-6, "AkimaSplineStatic value in last interval consistent with integral");

    const auto data = sc.copyInternalData();
    const auto copy = xraymc::AkimaSplineStatic<double, 12>::fromInternalData(data);
    check(copy.has_value() && copy.value()(1.2345) == sc(1.2345), "AkimaSplineStatic internal data round trip");
    check(!xraymc::AkimaSplineStatic<double, 12>::fromInternalData(std::vector<double>(7)).has_value(), "AkimaSplineStatic rejects internal data of wrong size");
}

void testCubicSplines()
{
    std::vector<std::pair<double, double>> data;
    for (int i = 0; i <= 20; ++i) {
        const double x = i * 0.25;
        data.push_back({ x, std::cos(x) });
    }
    xraymc::CubicSplineInterpolator<double> cs(data);
    bool knots = true;
    for (const auto& [x, y] : data)
        knots = knots && std::abs(cs(x) - y) < 1E-9;
    check(knots, "CubicSplineInterpolator passes through all knots");
    double maxErr = 0;
    for (double x = 0.25; x <= 4.75; x += 0.01)
        maxErr = std::max(maxErr, std::abs(cs(x) - std::cos(x)));
    check(maxErr < 5E-3, "CubicSplineInterpolator approximates cos (max err " + std::to_string(maxErr) + ")");
    checkNear(cs.meanValue(), std::sin(5.0) / 5.0, 2E-3, "CubicSplineInterpolator mean value");

    std::vector<std::pair<double, double>> line;
    for (int i = 0; i <= 10; ++i)
        line.push_back({ i, 3.0 - 0.5 * i });
    xraymc::CubicSplineInterpolator<double> cl(line);
    bool linear = true;
    for (double x = 0; x <= 10; x += 0.13)
        linear = linear && std::abs(cl(x) - (3.0 - 0.5 * x)) < 1E-9;
    check(linear, "natural CubicSplineInterpolator reproduces linear data");

    xraymc::CubicSplineInterpolatorStatic<double, 30> css(0.0, std::numbers::pi, [](double x) { return std::sin(x); });
    double maxErrS = 0;
    for (double x = 0; x <= std::numbers::pi; x += 0.01)
        maxErrS = std::max(maxErrS, std::abs(css(x) - std::sin(x)));
    check(maxErrS < 1E-3, "CubicSplineInterpolatorStatic approximates sin (max err " + std::to_string(maxErrS) + ")");
}

void testLSInterpolator()
{
    // Least squares cubic spline through noisy-free samples of a smooth function
    std::vector<double> x, y;
    for (int i = 0; i <= 400; ++i) {
        x.push_back(i * 0.01);
        y.push_back(std::log(1 + x.back()) + 0.3 * x.back());
    }
    xraymc::CubicLSInterpolator<double> ls(x, y, 15, false);
    double maxErr = 0;
    for (std::size_t i = 0; i < x.size(); ++i)
        maxErr = std::max(maxErr, std::abs(ls(x[i]) - y[i]));
    check(maxErr < 1E-3, "CubicLSInterpolator fits smooth function (max err " + std::to_string(maxErr) + ")");
    check(!ls.getDataTable().empty(), "CubicLSInterpolator has a knot table");
}

void testMatrix()
{
    // A x = b with known solution
    xraymc::Matrix<double> A(3, 3, { 2, 1, -1, -3, -1, 2, -2, 1, 2 });
    const auto x = A.solve({ 8, -11, -3 });
    check(x.has_value(), "Matrix solve of regular system has a solution");
    if (x) {
        checkNear(x.value()[0], 2, 1E-12, "Matrix solve x0");
        checkNear(x.value()[1], 3, 1E-12, "Matrix solve x1");
        checkNear(x.value()[2], -1, 1E-12, "Matrix solve x2");
        const auto b = A * x.value();
        check(std::abs(b[0] - 8) < 1E-12 && std::abs(b[1] + 11) < 1E-12 && std::abs(b[2] + 3) < 1E-12, "Matrix multiply recovers right hand side");
    }
    xraymc::Matrix<double> S(2, 2, { 1, 2, 2, 4 });
    check(!S.solve({ 1, 2 }).has_value(), "Matrix solve of singular system returns nullopt");

    xraymc::Matrix<double> M(2, 3, { 1, 2, 3, 4, 5, 6 });
    const auto T = M.transpose();
    check(T(0, 1) == 4 && T(2, 0) == 3 && T(1, 1) == 5, "Matrix transpose");
}

int main()
{
    std::cout << "Testing interpolation and integration\n";
    testLinearInterp();
    testLogLogInterp();
    testInterpolateTables();
    testRemovePoints();
    testVectorOps();
    testTrapz();
    testGauss();
    testAkima();
    testAkimaStatic();
    testCubicSplines();
    testLSInterpolator();
    testMatrix();
    return xraymctest::summary("testinterpolation");
}

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

#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <string>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;

template <Floating T>
constexpr T tol()
{
    return std::is_same_v<T, float> ? T { 1E-5 } : T { 1E-12 };
}

template <Floating T>
bool near(const std::array<T, 3>& a, const std::array<T, 3>& b, T eps = tol<T>())
{
    for (std::size_t i = 0; i < 3; ++i)
        if (std::abs(a[i] - b[i]) > eps)
            return false;
    return true;
}

template <Floating T>
std::string tname()
{
    return std::is_same_v<T, float> ? " [float]" : " [double]";
}

template <Floating T>
void testBasicArithmetic()
{
    const std::array<T, 3> a = { 1, 2, 3 };
    const std::array<T, 3> b = { -4, 5, 0.5 };

    check(vectormath::add(a, b) == std::array<T, 3> { -3, 7, 3.5 }, "add two vectors" + tname<T>());
    check(vectormath::add(a, b, a) == std::array<T, 3> { -2, 9, 6.5 }, "add three vectors" + tname<T>());
    check(vectormath::add(a, T { 1 }) == std::array<T, 3> { 2, 3, 4 }, "add scalar" + tname<T>());
    check(vectormath::add(T { 1 }, a) == std::array<T, 3> { 2, 3, 4 }, "add scalar (commuted)" + tname<T>());
    check(vectormath::subtract(a, b) == std::array<T, 3> { 5, -3, 2.5 }, "subtract" + tname<T>());
    check(vectormath::scale(a, T { 2 }) == std::array<T, 3> { 2, 4, 6 }, "scale by scalar" + tname<T>());
    check(vectormath::scale(T { 2 }, a) == std::array<T, 3> { 2, 4, 6 }, "scale by scalar (commuted)" + tname<T>());
    check(vectormath::scale(a, b) == std::array<T, 3> { -4, 10, 1.5 }, "element-wise scale" + tname<T>());
    check(vectormath::dot(a, b) == T { -4 + 10 + 1.5 }, "dot product" + tname<T>());
    check(vectormath::length_sqr(a) == T { 14 }, "squared length" + tname<T>());
    check(std::abs(vectormath::length(a) - std::sqrt(T { 14 })) < tol<T>(), "length" + tname<T>());

    const auto [l, r] = vectormath::splice(std::array<T, 6> { 1, 2, 3, 4, 5, 6 });
    check(l == std::array<T, 3> { 1, 2, 3 } && r == std::array<T, 3> { 4, 5, 6 }, "splice of 6-array" + tname<T>());
    check(vectormath::join(l, r) == std::array<T, 6> { 1, 2, 3, 4, 5, 6 }, "join is inverse of splice" + tname<T>());
}

template <Floating T>
void testCrossAndTriple()
{
    const std::array<T, 3> x = { 1, 0, 0 };
    const std::array<T, 3> y = { 0, 1, 0 };
    const std::array<T, 3> z = { 0, 0, 1 };
    check(vectormath::cross(x, y) == z, "x cross y = z (right-handed)" + tname<T>());
    check(vectormath::cross(y, z) == x, "y cross z = x" + tname<T>());
    check(vectormath::cross(z, x) == y, "z cross x = y" + tname<T>());
    check(vectormath::cross(y, x) == vectormath::scale(z, T { -1 }), "cross product is anti-commutative" + tname<T>());
    check(vectormath::cross(vectormath::join(x, y)) == z, "cross of joined 6-array" + tname<T>());
    check(vectormath::cross(std::array<std::array<T, 3>, 2> { x, y }) == z, "cross of array pair" + tname<T>());
    check(vectormath::tripleProduct(x, y, z) == T { 1 }, "triple product of unit basis is 1" + tname<T>());
    check(vectormath::tripleProduct(x, x, z) == T { 0 }, "triple product with repeated vector is 0" + tname<T>());

    auto state = xraymctest::seededState();
    bool orthogonal = true;
    for (int i = 0; i < 1000; ++i) {
        const std::array<T, 3> a = { state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1) };
        const std::array<T, 3> b = { state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1) };
        const auto c = vectormath::cross(a, b);
        orthogonal = orthogonal && std::abs(vectormath::dot(a, c)) < 100 * tol<T>() && std::abs(vectormath::dot(b, c)) < 100 * tol<T>();
    }
    check(orthogonal, "cross product is orthogonal to both arguments" + tname<T>());
}

template <Floating T>
void testNormalize()
{
    std::array<T, 3> v = { 3, 0, 4 };
    const auto n = vectormath::normalized(v);
    check(near(n, std::array<T, 3> { T { 0.6 }, 0, T { 0.8 } }), "normalized returns unit vector" + tname<T>());
    check(v == std::array<T, 3> { 3, 0, 4 }, "normalized does not modify argument" + tname<T>());
    vectormath::normalize(v);
    check(near(v, n), "normalize in place equals normalized" + tname<T>());
    check(std::abs(vectormath::length(v) - 1) < tol<T>(), "normalized vector has unit length" + tname<T>());
}

template <Floating T>
void testArgMinMax()
{
    check(vectormath::argmin3(std::array<T, 3> { 3, -1, 2 }) == 1, "argmin3 uses absolute value" + tname<T>());
    check(vectormath::argmin3(std::array<T, 3> { -5, 4, 0.5 }) == 2, "argmin3 picks last component" + tname<T>());
    check(vectormath::argmax3(std::array<T, 3> { 3, -7, 2 }) == 1, "argmax3 uses absolute value" + tname<T>());
    check(vectormath::argmax3(std::array<T, 3> { 9, 4, 0.5 }) == 0, "argmax3 picks first component" + tname<T>());
    check(vectormath::argmax3(std::array<int, 3> { 1, 2, 9 }) == 2, "argmax3 for integers" + tname<T>());
}

template <Floating T>
void testRotate()
{
    constexpr auto pi = std::numbers::pi_v<T>;
    const std::array<T, 3> x = { 1, 0, 0 };
    const std::array<T, 3> y = { 0, 1, 0 };
    const std::array<T, 3> z = { 0, 0, 1 };

    check(near(vectormath::rotate(x, z, pi / 2), y), "rotate x 90 deg about z gives y" + tname<T>());
    check(near(vectormath::rotate(y, x, pi / 2), z), "rotate y 90 deg about x gives z" + tname<T>());
    check(near(vectormath::rotate(z, y, pi / 2), x), "rotate z 90 deg about y gives x" + tname<T>());
    check(near(vectormath::rotate(x, z, pi), vectormath::scale(x, T { -1 })), "rotate 180 deg flips vector" + tname<T>());
    check(near(vectormath::rotate(z, z, T { 1.234 }), z), "rotation about own axis is identity" + tname<T>());
    check(near(vectormath::rotate(x, z, std::sin(pi / 2), std::cos(pi / 2)), y), "rotate with precomputed sin/cos" + tname<T>());

    // many small rotations add up to a full revolution
    auto t = vectormath::normalized(std::array<T, 3> { 1, 2, 3 });
    const auto t0 = t;
    const auto axis = vectormath::normalized(std::array<T, 3> { -1, 1, 0.5 });
    constexpr int N = 60;
    for (int i = 0; i < N; ++i)
        t = vectormath::rotate(t, axis, 2 * pi / N);
    check(near(t, t0, 100 * tol<T>()), "60 rotations of 6 deg returns to start" + tname<T>());

    auto state = xraymctest::seededState();
    bool preservesLength = true;
    bool preservesAxisComponent = true;
    for (int i = 0; i < 1000; ++i) {
        const std::array<T, 3> v = { state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1) };
        auto a = std::array<T, 3> { state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1) };
        vectormath::normalize(a);
        const auto ang = state.randomUniform<T>(-pi, pi);
        const auto r = vectormath::rotate(v, a, ang);
        preservesLength = preservesLength && std::abs(vectormath::length(r) - vectormath::length(v)) < 100 * tol<T>();
        preservesAxisComponent = preservesAxisComponent && std::abs(vectormath::dot(r, a) - vectormath::dot(v, a)) < 100 * tol<T>();
    }
    check(preservesLength, "rotation preserves length" + tname<T>());
    check(preservesAxisComponent, "rotation preserves component along axis" + tname<T>());
}

template <Floating T>
void testAngleBetween()
{
    constexpr auto pi = std::numbers::pi_v<T>;
    const std::array<T, 3> x = { 1, 0, 0 };
    const std::array<T, 3> y = { 0, 1, 0 };
    check(std::abs(vectormath::angleBetween(x, y) - pi / 2) < tol<T>(), "angle between orthogonal vectors is pi/2" + tname<T>());
    check(std::abs(vectormath::angleBetween(x, x)) < tol<T>(), "angle between equal vectors is 0" + tname<T>());
    check(std::abs(vectormath::angleBetween(x, vectormath::scale(x, T { -1 })) - pi) < tol<T>(), "angle between opposite vectors is pi" + tname<T>());
    check(std::abs(vectormath::angleBetween(x, std::array<T, 3> { 1, 1, 0 }) - pi / 4) < 10 * tol<T>(), "angle is independent of length" + tname<T>());
    // Very small angles is where the naive acos(dot) formula fails
    const T small = std::is_same_v<T, float> ? T { 1E-3 } : T { 1E-7 };
    const auto v = vectormath::rotate(x, std::array<T, 3> { 0, 0, 1 }, small);
    check(std::abs(vectormath::angleBetween(x, v) / small - 1) < T { 1E-2 }, "angleBetween is accurate for small angles" + tname<T>());
}

template <Floating T>
void testChangeBasis()
{
    // An orthonormal basis rotated from the standard basis
    const auto axis = vectormath::normalized(std::array<T, 3> { 1, 1, 1 });
    const auto b1 = vectormath::rotate(std::array<T, 3> { 1, 0, 0 }, axis, T { 0.3 });
    const auto b2 = vectormath::rotate(std::array<T, 3> { 0, 1, 0 }, axis, T { 0.3 });
    const auto b3 = vectormath::rotate(std::array<T, 3> { 0, 0, 1 }, axis, T { 0.3 });
    const std::array<T, 3> v = { T { 0.2 }, T { -1.5 }, T { 3 } };
    const auto w = vectormath::changeBasis(b1, b2, b3, v);
    check(near(w, vectormath::add(vectormath::scale(b1, v[0]), vectormath::scale(b2, v[1]), vectormath::scale(b3, v[2])), 10 * tol<T>()), "changeBasis is linear combination of basis" + tname<T>());
    check(near(vectormath::changeBasisInverse(b1, b2, b3, w), v, 10 * tol<T>()), "changeBasisInverse undoes changeBasis" + tname<T>());
}

template <Floating T>
void testPeturb()
{
    constexpr auto pi = std::numbers::pi_v<T>;
    auto state = xraymctest::seededState();
    bool angleOk = true;
    bool unitOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        std::array<T, 3> vec = { state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1), state.randomUniform<T>(-1, 1) };
        vectormath::normalize(vec);
        const auto theta = state.randomUniform<T>(T { 0.01 }, pi - T { 0.01 });
        const auto phi = state.randomUniform<T>(2 * pi);
        const auto vt = vectormath::peturb(vec, theta, phi);
        angleOk = angleOk && std::abs(vectormath::angleBetween(vec, vt) - theta) < T { 1E-3 };
        unitOk = unitOk && std::abs(vectormath::length(vt) - 1) < 10 * tol<T>();
    }
    check(angleOk, "peturb deflects by exactly theta" + tname<T>());
    check(unitOk, "peturb returns unit vectors" + tname<T>());

    // Repeated peturbation must not drift away from unit length
    std::array<T, 3> vec = { 0, 0, 1 };
    for (std::size_t i = 0; i < 10000; ++i)
        vec = vectormath::peturb(vec, state.randomUniform<T>(pi), state.randomUniform<T>(2 * pi));
    check(std::abs(vectormath::length(vec) - 1) < 10 * tol<T>(), "repeated peturb stays normalized" + tname<T>());

    // Azimuthal angle phi must be uniformly distributed around the original direction
    const std::array<T, 3> z = { 0, 0, 1 };
    constexpr std::size_t N = 100000;
    T meanX = 0, meanY = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto v = vectormath::peturb(z, pi / 2, state.randomUniform<T>(2 * pi));
        meanX += v[0];
        meanY += v[1];
    }
    check(std::abs(meanX / N) < T { 0.01 } && std::abs(meanY / N) < T { 0.01 }, "peturb with uniform phi is symmetric about direction" + tname<T>());
}

template <Floating T>
void tests()
{
    testBasicArithmetic<T>();
    testCrossAndTriple<T>();
    testNormalize<T>();
    testArgMinMax<T>();
    testRotate<T>();
    testAngleBetween<T>();
    testChangeBasis<T>();
    testPeturb<T>();
}

int main()
{
    std::cout << "Testing vectormath\n";
    tests<float>();
    tests<double>();
    return xraymctest::summary("testvectormath");
}

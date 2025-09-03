

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"

#include <array>
#include <iostream>
#include <numbers>

using namespace dxmc;

template <Floating T>
bool equal(T lh, T rh, const T thres = 1E-5)
{
    return std::abs(lh - rh) < thres;
}

template <Floating T>
bool equal(const std::array<T, 3>& lh, const std::array<T, 3>& rh)
{
    bool res = true;
    for (int i = 0; i < 3; ++i)
        res = res && equal(lh[i], rh[i]);
    return res;
}

template <typename T>
bool testRotate()
{
    const std::array<T, 3> vec = { 0, 0, 1 };
    const std::array<T, 3> z = { 0, 0, 1 };
    const std::array<T, 3> y = { 0, 1, 0 };
    const std::array<T, 3> x = { 1, 0, 0 };

    bool success = true;
    success = success && equal(vec, vectormath::rotate(vec, x, 2 * std::numbers::pi_v<T>));
    success = success && equal(vec, vectormath::rotate(vec, y, 2 * std::numbers::pi_v<T>));
    success = success && equal(vec, vectormath::rotate(vec, z, 2 * std::numbers::pi_v<T>));

    std::array<T, 60> ang;
    for (int i = 0; i < ang.size(); ++i) {
        ang[i] = (2 * std::numbers::pi_v<T>) / ang.size();
    }
    auto t = vec;
    for (const auto a : ang) {
        t = vectormath::rotate(t, x, a);
    }
    for (const auto a : ang) {
        t = vectormath::rotate(t, y, a);
    }
    success = success && equal(vec, t);

    return success;
}

template <typename T>
bool testPeturb()
{

    dxmc::RandomState state;
    bool success = true;
    for (std::size_t i = 0; i < 1E6; ++i) {
        std::array<T, 3> vec = {
            state.randomUniform<T>(-1, 1),
            state.randomUniform<T>(-1, 1),
            state.randomUniform<T>(-1, 1)
        };
        dxmc::vectormath::normalize(vec);

        const auto angle = state.randomUniform<T>(2 * std::numbers::pi_v<T>) - std::numbers::pi_v<T>;
        const auto cosang = std::cos(angle);
        if (std::abs(cosang) < T { 0.99 }) {
            const auto vt = vectormath::peturb(vec, angle, state.randomUniform<T>(2 * std::numbers::pi_v<T>));
            const auto res = dxmc::vectormath::angleBetween(vec, vt);
            success = success && equal(res, std::abs(angle), T { 1E-2 });
            const auto vt_lenght = vectormath::length(vt);
            success = success && equal(vt_lenght, T { 1 });
            if (!success) {
                auto test = false;
            }
        }
    }

    std::array<T, 3> vec = { 0, 0, 1 };
    for (std::size_t i = 0; i < 1E2; ++i) {
        const auto angle = state.randomUniform<T>(2 * std::numbers::pi_v<T>) - std::numbers::pi_v<T>;

        const auto phi = state.randomUniform<T>(2 * std::numbers::pi_v<T>);
        vec = vectormath::peturb(vec, angle, phi);
        const auto lenght = vectormath::length(vec);
        success = success && equal(lenght, T { 1 }, T { 1E-3 });
        if (!success) {
            auto test = false;
        }
        vectormath::normalize(vec);
    }

    return success;
}

template <Floating T>
bool tests()
{
    bool success = true;
    success = success && testRotate<T>();
    success = success && testPeturb<T>();

    std::cout << "Test vectormat floting point size of " << sizeof(T);
    if (success)
        std::cout << " SUCCESS\n";
    else
        std::cout << " FAILURE\n";

    return success;
}

int main(int argc, char* argv[])
{
    auto succ = tests<float>() && tests<double>();
    if (succ)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

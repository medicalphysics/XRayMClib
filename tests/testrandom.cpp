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

#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <vector>

using namespace xraymc;

constexpr double ERRF = 1e-4;
inline bool isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

template <typename T>
bool testDistribution()
{
    std::vector<T> x;
    for (std::size_t i = 0; i < 200; ++i) {
        const T val = T { -3 } + (i * T { 6 }) / (200 - 1);
        const auto w = std::exp(-val * val);
        x.push_back(w);
    }

    xraymc::RandomDistribution dist(x);

    constexpr std::size_t N = 1e6;

    std::vector<T> y(x.size(), T { 0 });
    RandomState state;
    for (std::size_t i = 0; i < N; i++) {
        const auto ind = dist.sampleIndex(state);
        y[ind] += T { 1 };
    }
    for (auto& yy : y) {
        yy /= N;
    }

    const auto x_sum = std::reduce(x.cbegin(), x.cend(), T { 0 });

    bool success = true;

    for (std::size_t i = 0; i < 200; ++i) {
        const auto xx = x[i] / x_sum;
        const auto diff = std::abs(y[i] - xx);
        success = success && diff < T { 0.001 };
        // std::cout << "Sample: " << y[i] << ", Analytical: " << xx << ", Diff: " << std::abs(y[i] - xx)  << std::endl;
    }
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test random discreet distribution for sizeof(T) = " << sizeof(T) << std::endl;
    return success;
}

template <typename T>
bool testSpecterDistribution()
{
    std::vector<std::pair<T, T>> x;
    std::vector<T> pos;

    for (std::size_t i = 0; i < 200; ++i) {
        const T val = T { -3 } + (i * T { 6 }) / (200 - 1);
        const auto w = std::exp(-val * val);
        x.push_back(std::make_pair(val, w));
        pos.push_back(val);
    }

    xraymc::SpecterDistribution dist(x);

    constexpr std::size_t N = 1e6;

    std::vector<T> y(x.size(), T { 0 });
    RandomState state;
    for (std::size_t i = 0; i < N; i++) {
        const auto e = dist.sampleValue(state);
        const auto pIdx = std::upper_bound(pos.cbegin(), pos.cend() - 1, e);
        const auto ind = std::distance(pos.cbegin(), pIdx);
        y[ind] += T { 1 };
    }
    for (auto& yy : y) {
        yy /= N;
    }

    T x_sum = 0;
    for (const auto& xx : x)
        x_sum += xx.second;

    bool success = true;

    for (std::size_t i = 0; i < 200; ++i) {
        const auto xx = x[i].second / x_sum;
        const auto diff = std::abs(y[i] - xx);
        success = success && diff < T { 0.001 };
        // std::cout << "Sample: " << y[i] << ", Analytical: " << xx << ", Diff: " << std::abs(y[i] - xx) << std::endl;
    }
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test random discreet specter distribution for sizeof(T) = " << sizeof(T) << std::endl;
    return success;
}

template <typename T>
bool testUniform()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto r = state.randomUniform<T>();
        const auto t = static_cast<std::size_t>(std::nextafter(r * N, T { 0 }));
        pcg[t]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    // check total number
    const auto pcg_sum = std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 });
    bool success = pcg_sum == S;

    T max_deviation = 0;
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        // deviation i percent
        const auto dev = (static_cast<T>(t) * N) / static_cast<T>(S) - T { 1 };
        max_deviation = std::max(std::abs(dev), max_deviation);
    }
    // max deviation is less than 0.5 %
    success = success && max_deviation * T { 100 } < T { 0.5 };

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test random uniform, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

template <typename T>
bool testUniformRange()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    constexpr T r1 = 5;
    constexpr T r2 = 10;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = state.randomUniform(r1, r2);
        const auto idx = static_cast<std::size_t>(std::nextafter(N * ((t - r1) / (r2 - r1)), T { 0 }));
        pcg[idx]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    // check total number
    const auto pcg_sum = std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 });
    bool success = pcg_sum == S;

    T max_deviation = 0;
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        // deviation i percent
        const auto dev = (static_cast<T>(t) * N) / static_cast<T>(S) - T { 1 };
        max_deviation = std::max(std::abs(dev), max_deviation);
    }
    // max deviation is less than 0.5 %
    success = success && max_deviation * T { 100 } < T { 0.5 };

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test random uniform range, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

bool testrandomInteger(std::uint8_t max = 128)
{
    std::vector<std::size_t> vals(max, 0);
    RandomState state;
    const std::size_t S = 1E6 * max;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto ind = state.randomInteger(max);
        vals[ind]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();
    auto sum = std::accumulate(vals.cbegin(), vals.cend(), std::size_t { 0 });

    bool success = sum == S;

    double max_diff = 0;
    for (auto n : vals) {
        const auto per_n = static_cast<double>(S) / max;
        const auto diff = static_cast<double>(n) / per_n - 1;
        max_diff = std::max(std::abs(diff), max_diff);
    }
    success = success && max_diff * 100 < 1.0;

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test random uniform index, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

template <typename T>
bool testCPDFdist(bool print = false)
{
    auto f = [](T x) { return std::exp(-(x * x) / 2); };
    const T xmin = -3;
    const T xmax = 3;
    const std::size_t N = 50;
    xraymc::CPDFSampling samp(xmin, xmax, f);
    xraymc::RandomState state;
    std::vector<T> hist(N, T { 0 });
    const auto step = (xmax - xmin) / (N);
    constexpr std::size_t Nsamp = 1E7;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < Nsamp; ++i) {
        const T x = samp(state);
        const auto idx = static_cast<std::size_t>(std::nextafter((x - xmin) / step, T { 0 }));
        hist[idx] += 1;
    }
    const auto end = std::chrono::high_resolution_clock::now();
    std::for_each(hist.begin(), hist.end(), [](auto& h) { h /= Nsamp; });

    auto func = hist;
    T func_sum = 0;
    for (std::size_t i = 0; i < hist.size(); ++i) {
        const auto x = xmin + step * i + step / 2;
        const auto p = f(x);
        func_sum += p;
        func[i] = p;
    }
    std::for_each(func.begin(), func.end(), [&](auto& p) { p /= func_sum; });

    std::vector<T> error(hist.size(), 0);

    if (print)
        std::cout << "x, sampled PDF, analytical pdf" << std::endl;

    T maxError = 0;
    for (std::size_t i = 0; i < hist.size(); ++i) {
        const auto x = xmin + step * i + step / 2;
        error[i] = std::abs(hist[i] - func[i]);
        if (print)
            std::cout << x << ", " << hist[i] << ", " << func[i] << std::endl;
        maxError = std::max(error[i], maxError);
    }
    auto success = maxError < T { 0.01 };

    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCESS" : "FAILURE";
    std::cout << msg << " Test Gaussian analytical pdf, time: " << time.count() << "ms" << std::endl;
    return success;
}

int main(int argc, char* argv[])
{
    std::cout << "Testing random number generator:" << std::endl;

    bool success = true;
    success = success && testDistribution<float>();
    success = success && testDistribution<double>();
    success = success && testSpecterDistribution<float>();
    success = success && testSpecterDistribution<double>();
    success = success && testrandomInteger();
    success = success && testUniform<double>();
    success = success && testUniform<float>();
    success = success && testUniformRange<double>();
    success = success && testUniformRange<float>();
    success = success && testCPDFdist<double>();
    success = success && testCPDFdist<float>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

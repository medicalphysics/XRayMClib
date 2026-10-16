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

#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;

template <Floating T>
std::string tname()
{
    return std::is_same_v<T, float> ? " [float]" : " [double]";
}

// Pearson chi-square statistic of observed counts against expected probabilities.
// Returns the statistic normalized as z = (chi2 - dof) / sqrt(2 dof), which is
// approximately standard normal for many bins.
double chiSquareZ(const std::vector<double>& counts, const std::vector<double>& probs)
{
    const auto N = std::reduce(counts.cbegin(), counts.cend(), 0.0);
    double chi2 = 0;
    std::size_t dof = 0;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        const auto e = probs[i] * N;
        if (e > 0) {
            chi2 += (counts[i] - e) * (counts[i] - e) / e;
            ++dof;
        }
    }
    dof -= 1;
    return (chi2 - dof) / std::sqrt(2.0 * dof);
}

void testSeeding()
{
    auto s1 = xraymctest::seededState(42);
    auto s2 = xraymctest::seededState(42);
    auto s3 = xraymctest::seededState(43);
    bool same = true;
    bool different = false;
    for (int i = 0; i < 1000; ++i) {
        const auto a = s1.pcg32();
        const auto b = s2.pcg32();
        const auto c = s3.pcg32();
        same = same && a == b;
        different = different || a != c;
    }
    check(same, "equal seeds give equal sequences");
    check(different, "different seeds give different sequences");

    RandomState r1, r2;
    check(r1.pcg32() != r2.pcg32() || r1.pcg32() != r2.pcg32(), "default constructed states are independently seeded");
}

template <Floating T>
void testUniform()
{
    auto state = xraymctest::seededState();
    constexpr std::size_t S = 10000000;
    constexpr std::size_t N = 100;
    std::vector<double> hist(N, 0);
    bool inRange = true;
    double sum = 0, sum2 = 0;
    for (std::size_t i = 0; i < S; ++i) {
        const auto r = state.randomUniform<T>();
        inRange = inRange && r >= T { 0 } && r < T { 1 };
        const auto idx = std::min(static_cast<std::size_t>(r * N), N - 1);
        hist[idx] += 1;
        sum += r;
        sum2 += r * r;
    }
    check(inRange, "randomUniform() in [0, 1)" + tname<T>());
    const auto mean = sum / S;
    const auto var = sum2 / S - mean * mean;
    checkNear(mean, 0.5, 5 * std::sqrt(1.0 / 12 / S), "randomUniform() mean is 1/2" + tname<T>());
    checkNear(var, 1.0 / 12, 1E-4, "randomUniform() variance is 1/12" + tname<T>());
    const auto z = chiSquareZ(hist, std::vector<double>(N, 1.0 / N));
    check(std::abs(z) < 5, "randomUniform() histogram passes chi-square (z = " + std::to_string(z) + ")" + tname<T>());
}

// The largest possible 32 bit output must still map to a value below 1, in
// particular for float where 2^32 - 1 is not representable.
template <Floating T>
void testUniformUpperEdge()
{
    // Search for a generator state whose next output is within 64 of 2^32 - 1.
    auto search = xraymctest::seededState(7);
    std::uint64_t found[2] = { 0, 0 };
    bool haveState = false;
    for (std::size_t i = 0; i < 400000000 && !haveState; ++i) {
        const std::uint64_t s0 = search.m_state[0];
        const std::uint64_t s1 = search.m_state[1];
        if (search.pcg32() >= 0xFFFFFFC0u) {
            found[0] = s0;
            found[1] = s1;
            haveState = true;
        }
    }
    check(haveState, "found generator state producing maximal output" + tname<T>());
    if (haveState) {
        xraymc::RandomState state(found);
        const auto r = state.randomUniform<T>();
        check(r < T { 1 }, "randomUniform() < 1 for maximal generator output" + tname<T>());
    }
}

template <Floating T>
void testUniformRange()
{
    auto state = xraymctest::seededState();
    constexpr std::size_t S = 1000000;
    constexpr T r1 = 5;
    constexpr T r2 = 10;
    bool inRange = true;
    bool maxInRange = true;
    double sum = 0;
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = state.randomUniform(r1, r2);
        inRange = inRange && t >= r1 && t <= r2;
        sum += t;
        const auto m = state.randomUniform(T { 3 });
        maxInRange = maxInRange && m >= 0 && m <= 3;
    }
    check(inRange, "randomUniform(min, max) within bounds" + tname<T>());
    check(maxInRange, "randomUniform(max) within bounds" + tname<T>());
    checkNear(sum / S, 7.5, 5 * 5 / std::sqrt(12.0 * S), "randomUniform(min, max) mean is midpoint" + tname<T>());
}

void testIntegers()
{
    auto state = xraymctest::seededState();
    constexpr std::uint32_t max = 37; // not a power of two, exercises rejection
    constexpr std::size_t S = 3700000;
    std::vector<double> hist(max, 0);
    bool inRange = true;
    for (std::size_t i = 0; i < S; ++i) {
        const auto v = state.randomInteger(max);
        inRange = inRange && v < max;
        if (v < max)
            hist[v] += 1;
    }
    check(inRange, "randomInteger(max) < max");
    const auto z = chiSquareZ(hist, std::vector<double>(max, 1.0 / max));
    check(std::abs(z) < 5, "randomInteger is uniform by chi-square (z = " + std::to_string(z) + ")");

    std::vector<double> hist2(10, 0);
    bool inRange2 = true;
    for (std::size_t i = 0; i < 1000000; ++i) {
        const auto v = state.randomUniform<std::size_t>(10);
        inRange2 = inRange2 && v < 10;
        if (v < 10)
            hist2[v] += 1;
    }
    check(inRange2, "randomUniform<size_t>(max) < max");
    const auto z2 = chiSquareZ(hist2, std::vector<double>(10, 0.1));
    check(std::abs(z2) < 5, "randomUniform<size_t> is uniform by chi-square (z = " + std::to_string(z2) + ")");

    bool rangeOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        const auto v = state.randomUniform<int>(-3, 4);
        rangeOk = rangeOk && v >= -3 && v < 4;
    }
    check(rangeOk, "randomUniform<int>(min, max) in [min, max)");

    bool cappedOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        const auto v = state.randomInteger32BitCapped(std::uint64_t { 1000 });
        cappedOk = cappedOk && v < 1000;
    }
    check(cappedOk, "randomInteger32BitCapped(max) < max");
}

void testNormal()
{
    auto state = xraymctest::seededState();
    constexpr std::size_t S = 1000000;
    double sum = 0, sum2 = 0, sumxy = 0;
    for (std::size_t i = 0; i < S; ++i) {
        const auto [a, b] = state.randomNormal();
        sum += a + b;
        sum2 += a * a + b * b;
        sumxy += a * b;
    }
    const auto n = 2.0 * S;
    checkNear(sum / n, 0, 5 / std::sqrt(n), "randomNormal mean is 0");
    checkNear(sum2 / n, 1, 0.01, "randomNormal variance is 1");
    checkNear(sumxy / S, 0, 5 / std::sqrt(1.0 * S), "randomNormal pair is uncorrelated");
}

template <Floating T>
void testRandomDistribution()
{
    std::vector<T> w;
    for (std::size_t i = 0; i < 200; ++i) {
        const T val = T { -3 } + (i * T { 6 }) / (200 - 1);
        w.push_back(std::exp(-val * val));
    }
    // include a zero weight bin, which must never be sampled
    w[0] = 0;
    xraymc::RandomDistribution<T> dist(w);
    constexpr std::size_t N = 2000000;
    std::vector<double> hist(w.size(), 0);
    auto state = xraymctest::seededState();
    for (std::size_t i = 0; i < N; ++i)
        hist[dist.sampleIndex(state)] += 1;
    const auto sum = std::reduce(w.cbegin(), w.cend(), T { 0 });
    std::vector<double> probs(w.size());
    std::transform(w.cbegin(), w.cend(), probs.begin(), [=](auto v) { return v / sum; });
    const auto z = chiSquareZ(hist, probs);
    check(std::abs(z) < 5, "RandomDistribution matches weights by chi-square (z = " + std::to_string(z) + ")" + tname<T>());
    check(hist[0] == 0, "RandomDistribution never samples zero weight" + tname<T>());

    xraymc::RandomDistribution<T> single;
    check(single.sampleIndex(state) == 0, "default RandomDistribution always returns index 0" + tname<T>());
}

template <Floating T>
void testSpecterDistribution()
{
    std::vector<std::pair<T, T>> spec;
    std::vector<T> energies;
    for (std::size_t i = 0; i < 100; ++i) {
        const T e = T { 10 } + i;
        spec.push_back(std::make_pair(e, std::exp(-(e - 50) * (e - 50) / 200)));
        energies.push_back(e);
    }
    xraymc::SpecterDistribution<T> dist(spec);
    constexpr std::size_t N = 2000000;
    auto state = xraymctest::seededState();
    std::vector<double> hist(spec.size(), 0);
    bool inRange = true;
    for (std::size_t i = 0; i < N; ++i) {
        const auto idx = dist.sampleIndex(state);
        hist[idx] += 1;
        const auto e = dist.sampleValue(state);
        inRange = inRange && e >= energies.front() && e <= energies.back();
    }
    T wsum = 0;
    for (const auto& [e, w] : spec)
        wsum += w;
    std::vector<double> probs(spec.size());
    for (std::size_t i = 0; i < spec.size(); ++i)
        probs[i] = spec[i].second / wsum;
    const auto z = chiSquareZ(hist, probs);
    check(std::abs(z) < 5, "SpecterDistribution index matches weights by chi-square (z = " + std::to_string(z) + ")" + tname<T>());
    check(inRange, "SpecterDistribution values within energy range" + tname<T>());

    xraymc::SpecterDistribution<T> dist2(energies, std::vector<T>(energies.size(), T { 1 }));
    double mean = 0;
    for (std::size_t i = 0; i < N; ++i)
        mean += dist2.sampleValue(state);
    mean /= N;
    // Each bin i is sampled uniformly in [e_i, e_i+1], except the last bin that returns its energy
    double expected = 0;
    for (std::size_t i = 0; i < energies.size() - 1; ++i)
        expected += energies[i] + 0.5;
    expected += energies.back();
    expected /= energies.size();
    checkNear(mean, expected, 0.05, "SpecterDistribution flat spectrum mean value" + tname<T>());

    // round trip through the internal data representation
    const auto data = dist.copyInteralData();
    const auto copy = xraymc::SpecterDistribution<T>::fromInternalData(data);
    check(copy.has_value() && copy.value().copyInteralData() == data, "SpecterDistribution internal data round trip" + tname<T>());
}

template <Floating T>
void testCPDFdist()
{
    auto f = [](T x) { return std::exp(-(x * x) / 2); };
    const T xmin = -3;
    const T xmax = 3;
    constexpr std::size_t N = 50;
    xraymc::CPDFSampling<T> samp(xmin, xmax, f);
    auto state = xraymctest::seededState();
    std::vector<double> hist(N, 0);
    const auto step = (xmax - xmin) / N;
    constexpr std::size_t Nsamp = 2000000;
    bool inRange = true;
    for (std::size_t i = 0; i < Nsamp; ++i) {
        const T x = samp(state);
        inRange = inRange && x >= xmin && x <= xmax;
        const auto idx = std::min(static_cast<std::size_t>(std::max((x - xmin) / step, T { 0 })), N - 1);
        hist[idx] += 1;
    }
    check(inRange, "CPDFSampling samples within [xmin, xmax]" + tname<T>());

    // exact bin probabilities from the error function
    std::vector<double> probs(N);
    const auto norm = std::erf(xmax / std::numbers::sqrt2) - std::erf(xmin / std::numbers::sqrt2);
    for (std::size_t i = 0; i < N; ++i) {
        const double a = xmin + step * i;
        const double b = a + step;
        probs[i] = (std::erf(b / std::numbers::sqrt2) - std::erf(a / std::numbers::sqrt2)) / norm;
    }
    double maxDev = 0;
    for (std::size_t i = 0; i < N; ++i)
        maxDev = std::max(maxDev, std::abs(hist[i] / Nsamp - probs[i]));
    check(maxDev < 1E-3, "CPDFSampling of Gaussian matches analytical bin probabilities (max dev = " + std::to_string(maxDev) + ")" + tname<T>());

    // truncated sampling
    const T xcut = 0.5;
    bool truncOk = true;
    for (std::size_t i = 0; i < 100000; ++i) {
        const auto x = samp(xcut, state);
        truncOk = truncOk && x <= xcut && x >= xmin;
    }
    check(truncOk, "CPDFSampling with max_x never exceeds max_x" + tname<T>());
}

int main()
{
    std::cout << "Testing random number generation\n";
    testSeeding();
    testUniform<float>();
    testUniform<double>();
    testUniformUpperEdge<float>();
    testUniformUpperEdge<double>();
    testUniformRange<float>();
    testUniformRange<double>();
    testIntegers();
    testNormal();
    testRandomDistribution<float>();
    testRandomDistribution<double>();
    testSpecterDistribution<float>();
    testSpecterDistribution<double>();
    testCPDFdist<float>();
    testCPDFdist<double>();
    return xraymctest::summary("testrandom");
}

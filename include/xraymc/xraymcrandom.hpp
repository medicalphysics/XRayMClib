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

Copyright 2019 Erlend Andersen
*/

#pragma once // include guard

#include "xraymc/floating.hpp"

#include <assert.h>
#include <concepts>
#include <execution>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief Per-thread PCG32 pseudo-random number generator state.
 *
 * Wraps the "Really Minimal PCG32" generator (M.E. O'Neill, 2014, Apache-2.0).
 * Each transport worker thread owns one `RandomState` instance; the object is
 * intentionally non-copyable to prevent accidental sharing of PRNG state between
 * threads.
 *
 * Provides:
 * - `randomUniform()` — floating-point value in [0, 1).
 * - `randomUniform(max)` / `randomUniform(min, max)` — bounded uniform samples.
 * - `randomNormal()` — a pair of standard-normal values via Box-Muller.
 * - `randomInteger(max)` / `randomInteger32BitCapped(max)` — unbiased integer draws.
 * - `pcg32()` — direct access to the underlying 32-bit generator output.
 */
class RandomState {
public:
    /**
     * @brief Seeds the generator from the platform's non-deterministic random device.
     *
     * Uses `std::random_device` to produce two independent 64-bit seed values,
     * giving a different sequence on each construction.
     */
    RandomState()
    {
        std::random_device d;
        std::uniform_int_distribution<std::uint64_t> dist(0);
        m_state[0] = static_cast<std::uint64_t>(dist(d));
        m_state[1] = static_cast<std::uint64_t>(dist(d));
    }
    /**
     * @brief Seeds the generator from a caller-supplied 2-element array.
     *
     * @param state Pointer to an array of two `uint64_t` seed values.
     *              Neither value should be zero; the sequence is fully determined
     *              by these two values.
     */
    RandomState(std::uint64_t state[2])
    {
        m_state[0] = state[0];
        m_state[1] = state[1];
    }
    RandomState(const RandomState&) = delete; ///< Non-copyable: each thread must own its own state.
    RandomState& operator=(const RandomState&) = delete; ///< Non-copyable: each thread must own its own state.

    /**
     * @brief Generate a random floating point number in range [0, 1.0)
     * @tparam T Must satisfy std::is_floating_point<T>::value == True
     * @return Random floating point in range [0, 1)
     */
    template <typename T = double>
    inline T randomUniform() noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
        const auto ui = pcg32();
        constexpr T uiMaxInv = T { 2.32830643653869628906e-010 };
        return ui * uiMaxInv;
    }

    /**
     * @brief Generate random uniform number in interval from 0 to max, exclusive
     * @tparam T Type of number, either integral or floating point
     * @param max Max of range
     * @return Random number in range [0, max).
     */
    template <typename T>
    inline T randomUniform(const T max) noexcept
    {
        if constexpr (std::is_floating_point<T>::value) {
            return randomUniform<T>() * max;
        } else if constexpr (std::is_integral<T>::value) {
            const std::uint32_t threshold = static_cast<std::uint32_t>(-max % max);
            for (;;) {
                const auto r = pcg32();
                if (r >= threshold)
                    return static_cast<T>(r % static_cast<std::uint32_t>(max));
            }
        } else {
            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "Must be integral or floating point value.");
        }
    }

    /**
     * @brief Generate random uniform number in interval from 0 to max, exclusive
     * @tparam T Type of number, either integral or floating point.
     * @param min Min of range.
     * @param max Max of range.
     * @return Random number in range [min, max).
     */
    template <typename T>
    inline T randomUniform(const T min, const T max) noexcept
    {
        if constexpr (std::is_floating_point<T>::value) {
            const T r = randomUniform<T>();
            const T range = max - min;
            return min + r * range;
        } else if constexpr (std::is_integral<T>::value) {
            const T range = max - min;
            const T r = randomUniform<T>(range);
            return min + r;
        } else
            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "Must be integral or floating point value.");
    }

    /**
     * @brief Generates a pair of independent standard-normal samples via Box-Muller.
     *
     * Draws two uniform values u1, u2 ∈ (0, 1) and returns
     * {√(−2 ln u1) · sin(2π u2), √(−2 ln u1) · cos(2π u2)},
     * both distributed as N(0, 1).
     *
     * @return Pair of independent N(0, 1) random variates.
     */
    inline std::pair<double, double> randomNormal() noexcept
    {
        const auto u1 = randomUniform();
        const auto u2 = randomUniform();
        const auto R = std::sqrt(-2.0 * std::log(u1));
        constexpr double pi2 = std::numbers::pi_v<double> * 2;
        return std::make_pair(R * std::sin(pi2 * u2), R * std::cos(pi2 * u2));
    }

    /**
     * @brief Returns a uniformly distributed integer in [0, max).
     *
     * Uses rejection sampling to eliminate modulo bias. Constrained to types
     * of at most 32 bits; use `randomInteger32BitCapped` for wider types.
     *
     * @tparam T Unsigned integral type, sizeof(T) ≤ 4.
     * @param max Exclusive upper bound.
     * @return Unbiased random integer in [0, max).
     */
    template <std::unsigned_integral T>
    inline T randomInteger(const T max) noexcept
    {
        static_assert(sizeof(max) <= 4, "This prng only supports up to 32 bit random integers, for a capped to 32 bit random integer use randomInteger32BitCapped instead");
        const T threshold = (static_cast<T>(-max)) % max;
        for (;;) {
            const auto r = pcg32();
            if (r >= threshold)
                return r % max;
        }
    }

    /**
     * @brief Returns a uniformly distributed integer in [0, max) for 64-bit or wider types.
     *
     * Identical rejection-sampling logic to `randomInteger`, but the result is capped
     * to the 32-bit generator output range. Use this overload when @p max fits in a
     * `uint32_t` but the variable type is `uint64_t` or larger.
     *
     * @tparam T Unsigned integral type, sizeof(T) > 4.
     * @param max Exclusive upper bound (must fit within uint32_t).
     * @return Unbiased random integer in [0, max).
     */
    template <std::unsigned_integral T>
    inline T randomInteger32BitCapped(const T max) noexcept
    {
        static_assert(sizeof(max) > 4, "This function is intended for 64 bit values or greater, use randomInteger method instead");
        const T threshold = (static_cast<T>(-max)) % max;
        for (;;) {
            const auto r = pcg32();
            if (r >= threshold)
                return r % max;
        }
    }

    /**
     * @brief Advances the PCG32 state and returns a 32-bit pseudo-random integer.
     *
     * Implementation of the "Really Minimal PCG32" generator
     * (© 2014 M.E. O'Neill / pcg-random.org, Apache-2.0).
     * Uses the LCG multiplier 6364136223846793005 and the XSH-RR output
     * function for good statistical quality at low cost.
     *
     * @return A pseudo-random 32-bit unsigned integer.
     */
    inline std::uint32_t pcg32() noexcept
    {
        const std::uint64_t oldstate = m_state[0];
        // Advance internal state
        m_state[0] = oldstate * 6364136223846793005ULL + (m_state[1] | 1);
        // Calculate output function (XSH RR), uses old state for max ILP
        const std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
        const std::uint32_t rot = oldstate >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }
    std::uint64_t m_state[2]; ///< PCG32 state: [0] is the LCG accumulator, [1] is the stream selector.
};

/**
 * @brief O(1) sampler for a discrete probability distribution using Vose's alias method.
 *
 * Builds a two-column alias table from a weight vector during construction.
 * Each `sampleIndex()` call draws exactly two uniform values — one to select a
 * column and one to accept or redirect to the alias — giving O(1) sampling
 * regardless of the number of bins.
 *
 * @tparam T Floating-point type (default: double).
 */
template <Floating T = double>
class RandomDistribution {
public:
    /**
     * @brief Constructs the alias table from a weight vector.
     * @param weights Per-bin probability weights; automatically normalized to sum to 1.
     *                Must contain at least two elements.
     */
    RandomDistribution(const std::vector<T>& weights)
    {
        m_data = generateTable(weights);
    }

    /// @brief Default constructor — creates a trivial single-bin distribution.
    RandomDistribution()
    {
        m_data.push_back(std::make_pair(T { 1 }, 0));
    }

    /**
     * @brief Sample an index from 0 to size() - 1 according to weights. This function is thread safe.
     * @param state Random state for sampling of an index.
     * @return index with probability according to weights[index].
     */
    std::size_t sampleIndex(RandomState& state) const
    {
        const auto r = state.randomUniform<T>();
        const auto k = state.randomUniform<std::size_t>(m_data.size());
        // return r < m_probs[k] ? k : m_alias[k];
        return r < m_data[k].first ? k : m_data[k].second;
    }

    /**
     * @brief Builds the alias table for a given weight vector.
     *
     * Implements the "squaring the histogram" (Vose) alias method:
     * normalizes weights to sum to N, then iteratively pairs under-full and
     * over-full buckets until all probabilities are in [0, 1].
     *
     * @param weights Per-bin weights (need not be normalized).
     * @return Vector of {acceptance probability, alias index} pairs.
     */
    static std::vector<std::pair<T, std::uint64_t>> generateTable(const std::vector<T>& weights)
    {
        // Squaring the histogram method
        const auto m_size = weights.size();

        std::vector<std::pair<T, std::uint64_t>> data(m_size);

        std::vector<T> norm_probs(m_size);
        std::vector<std::int64_t> large_block(m_size);
        std::vector<std::int64_t> small_block(m_size);

        std::int64_t num_small_block = 0;
        std::int64_t num_large_block = 0;

        const auto sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end(), T { 0.0 });
        const auto scale_factor = m_size / sum;
        std::transform(std::execution::par_unseq, weights.cbegin(), weights.cend(), norm_probs.begin(), [=](const auto w) -> T { return w * scale_factor; });

        for (std::int64_t i = m_size - 1; i >= 0; i--) {
            if (norm_probs[i] < T { 1.0 }) {
                small_block[num_small_block++] = i;
            } else {
                large_block[num_large_block++] = i;
            }
        }

        while (num_small_block && num_large_block) {
            const auto cur_small_block = small_block[--num_small_block];
            const auto cur_large_block = large_block[--num_large_block];

            data[cur_small_block] = std::make_pair(norm_probs[cur_small_block], cur_large_block);

            norm_probs[cur_large_block] = norm_probs[cur_large_block] + norm_probs[cur_small_block] - 1;
            if (norm_probs[cur_large_block] < 1) {
                small_block[num_small_block++] = cur_large_block;
            } else {
                large_block[num_large_block++] = cur_large_block;
            }
        }

        while (num_large_block) {
            data[large_block[--num_large_block]].first = 1;
        }

        while (num_small_block) {
            data[small_block[--num_small_block]].first = 1;
        }
        return data;
    }

private:
    std::vector<std::pair<T, std::uint64_t>> m_data;
};

/**
 * @brief O(1) sampler for a continuous energy spectrum represented as a histogram.
 *
 * Stores a beam spectrum as an alias table over energy bins. `sampleValue()` first
 * draws a bin index via the alias method, then returns a uniform random energy
 * within that bin's energy interval — approximating sampling from the piecewise-
 * constant PDF defined by the input weights.
 *
 * Provides `copyInternalData()` / `fromInternalData()` for serialization.
 *
 * @tparam T Floating-point type (default: double).
 */
template <Floating T = double>
class SpecterDistribution {
public:
    /// @brief Default constructor — creates an empty distribution.
    SpecterDistribution()
    {
    }

    /**
     * @brief Constructs a spectrum distribution from separate energy and weight vectors.
     *
     * Builds an alias table from the supplied weights so that `sampleValue()` runs in O(1).
     * Weights are normalized internally; they need not sum to unity on input.
     *
     * @param energy  Energy values in keV for each bin, in monotonically increasing order.
     *                Must have the same size as @p weights.
     * @param weights Relative probability weight for each energy bin. Must be non-negative.
     *                Size must match @p energy and be at least 1.
     */
    SpecterDistribution(const std::vector<T>& energy, const std::vector<T>& weights)
    {
        if (energy.size() != weights.size() || energy.size() == 0)
            return;
        const auto table = RandomDistribution<T>::generateTable(weights);

        m_data.resize(table.size());
        std::transform(std::execution::par_unseq, table.cbegin(), table.cend(), energy.cbegin(), m_data.begin(), [](const auto& tab, const T energy) {
            return DataElement { .alias = tab.second, .prob = tab.first, .energy = energy };
        });
    }

    /**
     * @brief Constructs a spectrum distribution from a vector of (energy, weight) pairs.
     *
     * Convenience overload that accepts a combined vector of `{energy_keV, weight}` pairs.
     * Energies are extracted into a monotonically increasing sequence and weights are
     * normalized internally via the alias method.
     *
     * @param specter Vector of `{energy [keV], relative weight}` pairs, one per bin.
     *                Must have at least one element.
     */
    SpecterDistribution(const std::vector<std::pair<T, T>>& specter)
    {
        std::vector<T> weights(specter.size());
        std::vector<T> energies(specter.size());

        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), energies.begin(), [](const auto& p) { return p.first; });
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), weights.begin(), [](const auto& p) { return p.second; });
        const auto table = RandomDistribution<T>::generateTable(weights);

        m_data.resize(table.size());
        std::transform(std::execution::par_unseq, table.cbegin(), table.cend(), energies.cbegin(), m_data.begin(), [](const auto& tab, const T energy) {
            return DataElement { .alias = tab.second, .prob = tab.first, .energy = energy };
        });
    }

    /**
     * @brief Sample an index from 0 to size() - 1 according to weights. This function is thread safe.
     * @param state Random state for sampling of an index.
     * @return index with probability according to weights[index].
     */
    std::size_t sampleIndex(RandomState& state) const
    {
        const auto r = state.randomUniform<T>();
        const auto k = state.randomUniform<std::size_t>(m_data.size());
        return r < m_data[k].prob ? k : m_data[k].alias;
    }

    /**
     * @brief Sample an energy value according to weights probability. This function is thread safe.
     * The sampling is done by first randomly sample an index into the energy vector. A random uniform energy in the interval energies[sampleIndex] and energies[sampleIndex+1] is returned. If sampleIndex is the last index in weights, the last energy value is returned.
     * @return a random energy according to weights probability
     */
    T sampleValue(RandomState& state) const
    {
        const std::size_t ind = sampleIndex(state);
        return ind < m_data.size() - 1 ? state.randomUniform(m_data[ind].energy, m_data[ind + 1].energy) : m_data[ind].energy;
    }

    /**
     * @brief Copy internal data.
     * Copy internal data to reconstruct the object from serialized data.
     * @return a vector of doubles representing the object state
     */
    std::vector<double> copyInteralData() const
    {
        std::vector<double> data;
        data.reserve(m_data.size() * 3);
        for (const auto& d : m_data) {
            data.push_back(static_cast<double>(d.alias));
            data.push_back(static_cast<double>(d.prob));
            data.push_back(static_cast<double>(d.energy));
        }
        return data;
    }

    /**
     * @brief Reconstruct SpecterDistribution from internal data.
     * Reconstruct SpecterDistribution from internal data.
     * @return an optional<SpecterDistribution<T>>
     */
    static std::optional<SpecterDistribution<T>> fromInternalData(const std::vector<double>& data)
    {
        SpecterDistribution<T> item;
        item.m_data.resize(data.size() / 3);
        std::size_t idx = 0;
        for (auto& d : item.m_data) {
            d.alias = static_cast<std::uint64_t>(data.at(idx++));
            d.prob = static_cast<T>(data.at(idx++));
            d.energy = static_cast<T>(data.at(idx++));
        }
        return item;
    }

protected:
private:
    struct DataElement {
        std::uint64_t alias = 0;
        T prob = 0;
        T energy = 0;
    };
    std::vector<DataElement> m_data;
};

/**
 * @brief Adaptive-grid CDF-inversion sampler for an arbitrary continuous PDF.
 *
 * At construction, an N-point grid is built over [xmin, xmax] by iteratively
 * refining the interval with the largest Simpson-integration error until exactly
 * N grid points are placed. Each grid element stores the local CDF value and two
 * shape coefficients (a, b) that enable analytic inversion of the piecewise-
 * rational CDF approximation — giving O(log N) sampling via a binary search on the
 * CDF array followed by a closed-form root solve.
 *
 * The two `operator()` overloads allow sampling from the full domain or from a
 * truncated upper bound.
 *
 * @tparam T Floating-point type.
 * @tparam N Number of grid points (default: 20). Higher values improve accuracy at
 *           the cost of construction time and memory.
 */
template <Floating T, std::size_t N = 20>
class CPDFSampling {
public:
    /// @brief Default constructor — grid is zero-initialised; must be assigned before use.
    CPDFSampling() = default;
    bool operator==(const CPDFSampling<T, N>&) const = default;

    /**
     * @brief Constructs the sampler by building an adaptive N-point grid for @p pdf.
     *
     * Starts with N/2 uniformly spaced points and iteratively inserts a new midpoint
     * into the interval whose integrated-PDF estimate has the largest absolute error
     * (measured against a 51-point Simpson quadrature) until exactly N points remain.
     * The final grid stores per-element CDF values and shape coefficients for O(1)
     * analytic CDF inversion within each interval.
     *
     * @tparam F  Callable type; must accept and return `T` (enforced by concept).
     * @param xmin  Lower bound of the sampling domain.
     * @param xmax  Upper bound of the sampling domain.
     * @param pdf   Probability density function to sample from. Need not be normalised.
     */
    template <std::regular_invocable<T> F>
        requires std::is_same<std::invoke_result_t<F, T>, T>::value
    CPDFSampling(T xmin, T xmax, const F& pdf)
    {
        std::size_t n_start = N / 2;

        std::vector<T> points(n_start);
        for (std::size_t i = 0; i < n_start; ++i) {
            points[i] = xmin + i * (xmax - xmin) / (n_start - 1);
        }

        while (points.size() != N) {
            auto grid = buildGrid(points, pdf);
            std::vector<T> error(grid.size() - 1);
            for (std::size_t i = 0; i < grid.size() - 1; ++i) {
                const T est = integrateGrid(grid[i], grid[i + 1]);
                const auto inter = simpson_integral(points[i], points[i + 1], pdf);
                error[i] = std::abs(est - inter);
            }
            auto max = std::max_element(error.cbegin(), error.cend());
            const auto idx = std::distance(error.cbegin(), max);
            auto g0 = grid.cbegin() + idx;
            auto g1 = g0 + 1;
            const T new_point = g0->x + (g1->x - g0->x) * T { 0.5 };
            points.insert(points.cbegin() + idx + 1, new_point);
        }

        auto grid = buildGrid(points, pdf);
        std::copy(grid.cbegin(), grid.cend(), m_grid.begin());
    }
    /**
     * @brief Samples a value from the full distribution [xmin, xmax].
     *
     * Draws a uniform CDF variate, locates the enclosing grid interval via binary
     * search, then analytically inverts the piecewise-rational CDF to obtain the
     * corresponding x value. Thread-safe when @p state is not shared.
     *
     * @param state Per-thread PRNG state.
     * @return Random variate distributed according to the PDF supplied at construction.
     */
    T operator()(RandomState& state) const
    {
        const T r1 = state.randomUniform<T>(m_grid[0].e, m_grid[N - 1].e);
        auto pos1 = std::upper_bound(m_grid.cbegin() + 1, m_grid.cend() - 1, r1, [](const auto num, const auto& element) -> bool { return num < element.e; });
        auto pos0 = pos1 - 1;

        const auto v = r1 - pos0->e;
        const auto delta = pos1->e - pos0->e;

        const auto nom = (1 + pos0->a + pos0->b) * delta * v;
        const auto den = delta * delta + pos0->a * delta * v + pos0->b * v * v;
        const auto res = pos0->x + nom * (pos1->x - pos0->x) / den;
        return res;
    }
    /**
     * @brief Samples a value from the distribution truncated to [xmin, max_x].
     *
     * Restricts the CDF range to the grid interval containing @p max_x, samples a
     * CDF variate uniformly from that restricted range, and inverts analytically.
     * If the inversion overshoots @p max_x (possible due to the rational approximation),
     * the draw is rejected and retried — making this a rejection sampler within the
     * truncated domain. Thread-safe when @p state is not shared.
     *
     * @param max_x  Exclusive upper bound; the returned value is guaranteed ≤ @p max_x.
     * @param state  Per-thread PRNG state.
     * @return Random variate in [xmin, max_x] distributed proportionally to the PDF.
     */
    T operator()(T max_x, RandomState& state) const
    {
        const auto max_xpos = std::upper_bound(m_grid.cbegin() + 1, m_grid.cend() - 1, max_x, [](const auto num, const auto& element) -> bool { return num < element.x; });
        T res;
        do {
            const T r1 = state.randomUniform<T>(m_grid[0].e, max_xpos->e);
            auto pos1 = std::upper_bound(m_grid.cbegin() + 1, max_xpos, r1, [](const auto num, const auto& element) -> bool { return num < element.e; });
            auto pos0 = pos1 - 1;

            const auto v = r1 - pos0->e;
            const auto delta = pos1->e - pos0->e;

            const auto nom = (1 + pos0->a + pos0->b) * delta * v;
            const auto den = delta * delta + pos0->a * delta * v + pos0->b * v * v;
            res = pos0->x + nom * (pos1->x - pos0->x) / den;
        } while (res > max_x);
        return res;
    }

protected:
    /**
     * @brief One node of the adaptive CDF grid.
     *
     * Stores the PDF-domain coordinate @p x, the cumulative probability @p e at that
     * point, and two rational-approximation shape coefficients @p a and @p b used for
     * analytic CDF inversion within the interval [x, x_next].
     */
    struct GridEl {
        T x = 0; ///< PDF-domain coordinate of this grid node.
        T e = 0; ///< Cumulative probability (CDF value) at this node.
        T a = 0; ///< First rational-approximation shape coefficient.
        T b = 0; ///< Second rational-approximation shape coefficient.
        bool operator==(const GridEl&) const = default;
    };

    /**
     * @brief Estimates the integrated PDF over the interval [@p g0.x, @p g1.x].
     *
     * Uses a numerical quadrature over the rational CDF approximation defined by
     * the shape coefficients in @p g0, providing an estimate that can be compared
     * against a direct Simpson integration for error-driven grid refinement.
     *
     * @param g0      Left grid node.
     * @param g1      Right grid node.
     * @param Nsteps  Number of quadrature steps (default: 50).
     * @return Estimated integral of the PDF over the interval.
     */
    static T integrateGrid(const GridEl& g0, const GridEl& g1, std::size_t Nsteps = 50)
    {
        T sum = 0;
        const T step = (g1.x - g0.x) / (Nsteps - 1);
        for (std::size_t i = 1; i < Nsteps; ++i) {
            const T x = g0.x + step * i;
            const auto t = (x - g0.x) / (g1.x - g0.x);
            const auto a = g0.a;
            const auto b = g0.b;
            const auto p1 = 1 + a + b - a * t;
            const auto n = p1 / (2 * b * t) * (1 - std::sqrt(1 - (4 * b * t * t) / (p1 * p1)));

            const auto r1 = 1 + a * n + b * n * n;
            const auto r2 = (1 + a + b) * (1 - b * n * n);

            sum += r1 * r1 * (g1.e - g0.e) / (r2 * (g1.x - g0.x));
        }
        return sum * step;
    }
    /**
     * @brief Builds a CDF grid from a set of x-domain sample points and a PDF functor.
     *
     * Evaluates the PDF at each point and accumulates a running Simpson integral to
     * produce unnormalized CDF values. For each consecutive pair of points, computes
     * the rational-approximation shape coefficients @p a and @p b that allow analytic
     * CDF inversion within that interval.
     *
     * @tparam F   Callable accepting and returning `T`.
     * @param points  Ordered x-domain knot positions.
     * @param pdf     PDF functor; called once per knot.
     * @return Vector of @p GridEl nodes with x, cumulative-e, a, and b filled in.
     */
    template <typename F>
    static std::vector<GridEl> buildGrid(const std::vector<T> points, const F& pdf)
    {
        // cumulative prob (unnormalized)
        std::vector<T> e(points.size());
        e[0] = pdf(points[0]);
        for (std::size_t i = 1; i < e.size(); ++i) {
            const auto x0 = points[i - 1];
            const auto x1 = points[i];
            e[i] = e[i - 1] + simpson_integral(x0, x1, pdf);
        }

        std::vector<GridEl> grid;

        for (std::size_t i = 0; i < points.size() - 1; ++i) {

            const T ex = (e[i + 1] - e[i]) / (points[i + 1] - points[i]);
            const T b = T { 1 } - ex * ex / (pdf(points[i]) * pdf(points[i + 1]));
            const T a = ex / pdf(points[i]) - b - T { 1 };
            grid.push_back({ points[i], e[i], a, b });
        }
        grid.push_back({ points.back(), e.back(), T { 0 }, T { 0 } });

        return grid;
    }

    /**
     * @brief Computes the definite integral of @p pdf over [x0, x1] via composite Simpson's rule.
     *
     * Uses @p Nsimp uniformly spaced points (must be odd). Serves as the reference
     * integrator for error-driven grid refinement during construction.
     *
     * @tparam F      Callable accepting and returning `T`.
     * @tparam Nsimp  Number of quadrature points (default: 51, must be odd).
     * @param x0   Lower integration bound.
     * @param x1   Upper integration bound.
     * @param pdf  PDF functor to integrate.
     * @return Approximate integral of @p pdf over [x0, x1].
     */
    template <typename F, int Nsimp = 51>
    static T simpson_integral(T x0, T x1, const F& pdf)
    {
        const T h = (x1 - x0) / (Nsimp - 1);
        const T p1 = pdf(x0) + pdf(x1);
        T p2 = 0;
        for (std::size_t i = 1; i < Nsimp - 1; i += 2)
            p2 += pdf(x0 + i * h);
        T p3 = 0;
        for (std::size_t i = 2; i < Nsimp - 2; i += 2)
            p3 += pdf(x0 + i * h);
        return h * (p1 + 4 * p2 + 2 * p3) / T { 3 };
    }

private:
    std::array<GridEl, N> m_grid;
};

}
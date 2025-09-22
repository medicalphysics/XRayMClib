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

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <iterator>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace xraymc {

template <Floating T>
inline T interp(const T x0, const T x1, const T y0, const T y1, const T x)
{
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}
template <Floating T, Floating U>
inline U interp(T x[2], T y[2], U xi)
{
    return y[0] + (y[1] - y[0]) * (xi - x[0]) / (x[1] - x[0]);
}

template <Floating T>
inline T interp(const std::pair<T, T>& v1, const std::pair<T, T>& v2, T x)
{
    return v1.second + (v2.second - v1.second) * (x - v1.first) / (v2.first - v1.first);
}

template <Floating T>
inline T logloginterp(T x0, T x1, T y0, T y1, T x)
{
    // we do not test for negative values, instead we rely on std::log10 to return NaN values
    const double value = std::log10(y0) + (std::log10(y1 / y0) / std::log10(x1 / x0) * std::log10(x / x0)); // std::log 10 always promotes to double
    return std::pow(10., value); // std::pow always promotes to doubles
}

template <Floating T, Floating U>
inline T logloginterp(T x[2], T y[2], U xi)
{
    const T x0 = x[0];
    const T x1 = x[1];
    const T y0 = y[0];
    const T y1 = y[1];
    // we do not test for negative values, instead we rely on std::log10 to return NaN values
    const double value = std::log10(y0) + (std::log10(y1 / y0) / std::log10(x1 / x0) * std::log10(xi / x0)); // std::log 10 always promotes to double
    return std::pow(10., value); // std::pow always promotes to doubles
}

template <typename It, Floating T>
    requires std::is_same_v<typename std::iterator_traits<It>::value_type, T>
T interpolate(It xbegin, It xend, It ybegin, It yend, T xvalue)
{
    auto upper = std::upper_bound(xbegin, xend, xvalue);
    if (upper == xbegin)
        return *ybegin;
    if (upper == xend)
        return *(yend - 1);
    auto lower = upper;
    std::advance(lower, -1);

    auto lowery = ybegin;
    std::advance(lowery, std::distance(xbegin, lower));
    auto uppery = ybegin;
    std::advance(uppery, std::distance(xbegin, upper));

    return interp(*lower, *upper, *lowery, *uppery, xvalue);
}
template <Floating T, bool EXTRAPOLATE_DOWN = true, bool EXTRAPOLATE_UP = true>
inline T interpolate(const std::vector<std::pair<T, T>>& data, T x)
{

    if constexpr (!EXTRAPOLATE_DOWN) {
        if (x < data[0].first)
            return 0;
    }
    if constexpr (!EXTRAPOLATE_UP) {
        if (x > data.back().first)
            return 0;
    }

    const auto val = std::make_pair(x, T { 0 });
    auto upper = std::upper_bound(data.begin() + 1, data.end() - 1, val, [](const auto& lh, const auto& rh) -> bool { return lh.first < rh.first; });
    auto lower = upper - 1;

    return interp(*lower, *upper, x);
}

template <Floating T>
inline std::vector<T> interpolate(const std::vector<std::pair<T, T>>& data, const std::vector<T>& x)
{
    std::vector<T> y(x.size());
    if (data.size() < 4) {
        for (std::size_t i = 0; i < x.size(); ++i)
            y[i] = interpolate(data, x[i]);
        return y;
    }
    auto upper = data.cbegin() + 1;
    auto dStop = data.cend() - 1;
    for (std::size_t i = 0; i < x.size(); ++i) {
        while (x[i] > upper->first && upper != dStop)
            ++upper;
        const auto lower = upper - 1;
        y[i] = interp(lower->first, upper->first, lower->second, upper->second, x[i]);
    }
    return y;
}

template <Floating T>
void removeUnneededInterpolationPoints(std::vector<std::pair<T, T>>& data, T epsilon = 1E-6)
{
    std::size_t right = 2;
    while (right < data.size()) {
        const auto& v1 = data[right - 2];
        const auto& v12 = data[right - 1];
        const auto& v2 = data[right];
        const auto guess = interp(v1, v2, v12.first);
        if (std::abs(guess - v12.second) < epsilon * v12.second) {
            // remove point
            data.erase(data.cbegin() + (right - 1));
        } else {
            right++;
        }
    }
    data.shrink_to_fit();
}

// Variadic template to add vectors of equal size
template <Floating T>
std::vector<T> addVectors(const std::vector<T>& f, const std::vector<T>& l)
{
    auto res = f;
    std::transform(std::execution::par_unseq, l.cbegin(), l.end(), res.cbegin(), res.begin(), std::plus<T>());
    return res;
}
template <Floating T>
std::vector<T> addVectors(std::vector<T>& f, const std::vector<T>& l)
{
    std::transform(std::execution::par_unseq, l.cbegin(), l.end(), f.cbegin(), f.begin(), std::plus<T>());
    return f;
}
template <Floating T, typename... args>
std::vector<T> addVectors(const std::vector<T>& f, const std::vector<T>& l, args... other)
{
    auto res = addVectors(f, l);
    return addVectors(res, other...);
}

template <Floating T>
std::vector<T> trapz_cum(const std::vector<T>& x, const std::vector<T>& f)
{
    std::vector<T> integ(f.size());
    integ[0] = T { 0.0 };
    for (std::size_t i = 1; i < f.size(); ++i) {
        integ[i] = integ[i - 1] + (f[i - 1] + f[i]) * T { 0.5 } * (x[i] - x[i - 1]);
    }
    return integ;
}
template <Floating T>
T trapz(const std::vector<T>& x, const std::vector<T>& f)
{
    T integ = 0;
    for (std::size_t i = 1; i < f.size(); ++i) {
        integ += (f[i - 1] + f[i]) * T { 0.5 } * (x[i] - x[i - 1]);
    }
    return integ;
}

template <Floating T>
T trapz(const std::vector<std::pair<T, T>>& f)
{
    T integ = 0;
    for (std::size_t i = 1; i < f.size(); ++i) {
        integ += (f[i - 1].second + f[i].second) * T { 0.5 } * (f[i].first - f[i - 1].first);
    }
    return integ;
}

template <Floating T>
T trapz(const std::vector<std::pair<T, T>>& f, T start, T stop)
{
    T integ = 0;
    const auto start_x = std::clamp(start, f.front().first, f.back().first);
    const auto stop_x = std::clamp(stop, f.front().first, f.back().first);
    std::size_t startIdx = f.size() - 1;
    std::size_t stopIdx = 0;

    for (std::size_t i = 1; i < f.size(); ++i) {
        if (f[i - 1].first > start_x && f[i].first < stop_x) {
            integ += (f[i - 1].second + f[i].second) * T { 0.5 } * (f[i].first - f[i - 1].first);
            stopIdx = std::max(stopIdx, i);
            startIdx = std::min(startIdx, i - 1);
        }
    }
    // remainder
    const auto start_y = interpolate(f, start_x);
    const auto stop_y = interpolate(f, stop_x);
    if (f[startIdx].first - start_x > T { 0 })
        integ += (start_y + f[startIdx].second) * T { 0.5 } * (f[startIdx].first - start_x);
    if (stop_x - f[stopIdx].first > T { 0 })
        integ += (f[stopIdx].second + stop_y) * T { 0.5 } * (stop_x - f[stopIdx].first);
    return integ;
}

template <Floating T>
constexpr T gaussIntegration(const T start, const T stop, const std::array<T, 20>& gaussPoints)
{
    constexpr std::array<T, 20> weights = {
        1.5275338713072585E-01,
        1.4917298647260375E-01,
        1.4209610931838205E-01,
        1.3168863844917663E-01,
        1.1819453196151842E-01,
        1.0193011981724044E-01,
        8.3276741576704749E-02,
        6.2672048334109064E-02,
        4.0601429800386941E-02,
        1.7614007139152118E-02,
        1.5275338713072585E-01,
        1.4917298647260375E-01,
        1.4209610931838205E-01,
        1.3168863844917663E-01,
        1.1819453196151842E-01,
        1.0193011981724044E-01,
        8.3276741576704749E-02,
        6.2672048334109064E-02,
        4.0601429800386941E-02,
        1.7614007139152118E-02
    };
    const T interval_half = (stop - start) * T { 0.5 };
    const T value = std::transform_reduce(std::execution::unseq, weights.cbegin(), weights.cend(), gaussPoints.cbegin(), T { 0 }, std::plus<>(), std::multiplies<>());
    return value * interval_half;
}
template <Floating T>
constexpr std::array<T, 20> gaussIntegrationPoints(const T start, const T stop)
{
    constexpr std::array<T, 20> x_val = {
        7.6526521133497334E-02,
        2.2778585114164508E-01,
        3.7370608871541956E-01,
        5.1086700195082710E-01,
        6.3605368072651503E-01,
        7.4633190646015079E-01,
        8.3911697182221882E-01,
        9.1223442825132591E-01,
        9.6397192727791379E-01,
        9.9312859918509492E-01,
        -7.6526521133497334E-02,
        -2.2778585114164508E-01,
        -3.7370608871541956E-01,
        -5.1086700195082710E-01,
        -6.3605368072651503E-01,
        -7.4633190646015079E-01,
        -8.3911697182221882E-01,
        -9.1223442825132591E-01,
        -9.6397192727791379E-01,
        -9.9312859918509492E-01
    };
    std::array<T, 20> function_points;
    const T xi = (stop - start) * T { 0.5 };
    const T xii = (stop + start) * T { 0.5 };
    std::transform(std::execution::unseq, x_val.cbegin(), x_val.cend(), function_points.begin(), [=](const auto x) { return x * xi + xii; });
    return function_points;
}

template <Floating T, std::regular_invocable<T> F>
    requires std::is_same<std::invoke_result_t<F, T>, T>::value
constexpr T gaussIntegration(const T start, const T stop, const F function)
{
    // F is a class or function that can be called with a Floating type argument and will return a Floating type value
    auto function_points = gaussIntegrationPoints(start, stop);
    std::array<T, 20> function_values;
    std::transform(std::execution::unseq, function_points.cbegin(), function_points.cend(), function_values.begin(), [&](const auto x) { return function(x); });
    return gaussIntegration(start, stop, function_values);
}

template <Floating T>
class AkimaSpline {
public:
    AkimaSpline()
    {
        m_data.resize(2);
        m_data[1].x = 1;
    }
    AkimaSpline(const std::vector<std::pair<T, T>>& data)
    {
        setup(data);
    }

    AkimaSpline(const std::vector<T>& x, const std::vector<T>& y)
    {
        std::vector<std::pair<T, T>> data(std::min(x.size(), y.size()));
        for (std::size_t i = 0; i < std::min(x.size(), y.size()); ++i) {
            data[i].first = x[i];
            data[i].second = y[i];
        }

        setup(data);
    }

    T operator()(T x) const
    {
        const Interval v { .x = x };
        const auto i1 = std::upper_bound(m_data.cbegin() + 1, m_data.cend() - 1, v, [](const auto& lh, const auto& rh) {
            return lh.x < rh.x;
        });
        const auto i0 = i1 - 1;
        const auto& seg = *i0;
        const auto xd = x - seg.x;

        return seg.a + seg.b * xd + seg.c * xd * xd + seg.d * xd * xd * xd;
    }

    void scale(T s)
    {
        for (auto& v : m_data) {
            v.a *= s;
            v.b *= s;
            v.c *= s;
            v.d *= s;
        }
    }

    T integral(T start, T stop) const
    {
        start = std::max(start, m_data.front().x);

        T integrand = 0;
        for (std::size_t i = 0; i < m_data.size() - 1; ++i) {
            if (m_data[i].x <= start && start <= m_data[i + 1].x && stop > start) {

                const auto& x0 = m_data[i].x;
                const auto& a = m_data[i].a;
                const auto& b = m_data[i].b;
                const auto& c = m_data[i].c;
                const auto& d = m_data[i].d;

                const auto xa = start - x0;
                auto xb = m_data[i + 1].x - x0;
                if (stop < m_data[i + 1].x)
                    xb = stop - x0;

                const auto fa = a * xa + b * xa * xa / 2 + c * xa * xa * xa / 3 + d * xa * xa * xa * xa / 4;
                const auto fb = a * xb + b * xb * xb / 2 + c * xb * xb * xb / 3 + d * xb * xb * xb * xb / 4;

                integrand += (fb - fa);

                start = m_data[i + 1].x;
            }
        }
        if (start < stop) {
            // we have a rest
            const auto& x0 = m_data.back().x;
            const auto& a = m_data.back().a;
            const auto& b = m_data.back().b;
            const auto& c = m_data.back().c;
            const auto& d = m_data.back().d;

            const auto xa = start - x0;
            auto xb = stop - x0;

            const auto fa = a * xa + b * xa * xa / 2 + c * xa * xa * xa / 3 + d * xa * xa * xa * xa / 4;
            const auto fb = a * xb + b * xb * xb / 2 + c * xb * xb * xb / 3 + d * xb * xb * xb * xb / 4;

            integrand += (fb - fa);
        }
        return integrand;
    }

    void setup(std::vector<std::pair<T, T>> data)
    {
        if (data.size() < 5) {
            m_data.resize(2);
            m_data[1].x = 1;
            return;
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        m_data.resize(data.size() - 1);
        std::vector<T> m(data.size() - 1);
        for (std::size_t i = 0; i < data.size() - 1; ++i) {
            m[i] = (data[i + 1].second - data[i].second) / (data[i + 1].first - data[i].first);
        }
        std::vector<T> s(data.size());
        s[0] = m[0];
        s[1] = (m[0] + m[1]) / 2;
        for (std::size_t i = 2; i < data.size() - 2; ++i) {
            const auto den = std::abs(m[i + 1] - m[i]) + std::abs(m[i - 1] - m[i - 2]);
            if (den > 0.000001)
                s[i] = (std::abs(m[i + 1] - m[i]) * m[i - 1] + std::abs(m[i - 1] - m[i - 2]) * m[i]) / den;
            else
                s[i] = (m[i - 1] + m[i]) / 2;
        }
        s[data.size() - 2] = (m[data.size() - 3] + m[data.size() - 2]) / 2;
        s[data.size() - 1] = m[data.size() - 2];

        for (std::size_t i = 0; i < s.size() - 1; ++i) {
            m_data[i].x = data[i].first;
            m_data[i].a = data[i].second;
            m_data[i].b = s[i];
            const auto xdiff = data[i + 1].first - data[i].first;
            m_data[i].c = (3 * m[i] - 2 * s[i] - s[i + 1]) / xdiff;
            m_data[i].d = (s[i] + s[i + 1] - 2 * m[i]) / (xdiff * xdiff);
        }
    }

private:
    struct Interval {
        T x = 0, a = 1, b = 0, c = 0, d = 0;
    };
    std::vector<Interval> m_data;
};

template <Floating T, std::size_t N_KNOTS = 5>
class AkimaSplineStatic {
public:
    AkimaSplineStatic()
    {
        static_assert(N_KNOTS >= 5);
        for (std::size_t i = 0; i < N_KNOTS; ++i)
            m_data[i].x = static_cast<T>(i);
    }
    AkimaSplineStatic(const std::vector<std::pair<T, T>>& data)
    {
        static_assert(N_KNOTS >= 5);
        setup(data);
    }

    AkimaSplineStatic(const std::vector<T>& x, const std::vector<T>& y)
    {
        static_assert(N_KNOTS >= 5);
        std::vector<std::pair<T, T>> data(std::min(x.size(), y.size()));
        for (std::size_t i = 0; i < std::min(x.size(), y.size()); ++i) {
            data[i].first = x[i];
            data[i].second = y[i];
        }
        setup(data);
    }

    void scale(T s)
    {
        for (auto& v : m_data) {
            v.a *= s;
            v.b *= s;
            v.c *= s;
            v.d *= s;
        }
    }

    T integral(T start, T stop) const
    {
        start = std::max(start, m_data.front().x);

        T integrand = 0;
        for (std::size_t i = 0; i < m_data.size() - 1; ++i) {
            if (m_data[i].x <= start && start <= m_data[i + 1].x && stop > start) {

                const auto& x0 = m_data[i].x;
                const auto& a = m_data[i].a;
                const auto& b = m_data[i].b;
                const auto& c = m_data[i].c;
                const auto& d = m_data[i].d;

                const auto xa = start - x0;
                auto xb = m_data[i + 1].x - x0;
                if (stop < m_data[i + 1].x)
                    xb = stop - x0;

                const auto fa = a * xa + b * xa * xa / 2 + c * xa * xa * xa / 3 + d * xa * xa * xa * xa / 4;
                const auto fb = a * xb + b * xb * xb / 2 + c * xb * xb * xb / 3 + d * xb * xb * xb * xb / 4;

                integrand += (fb - fa);

                start = m_data[i + 1].x;
            }
        }
        if (start < stop) {
            // we have a rest
            const auto& x0 = m_data.back().x;
            const auto& a = m_data.back().a;
            const auto& b = m_data.back().b;
            const auto& c = m_data.back().c;
            const auto& d = m_data.back().d;

            const auto xa = start - x0;
            auto xb = stop - x0;

            const auto fa = a * xa + b * xa * xa / 2 + c * xa * xa * xa / 3 + d * xa * xa * xa * xa / 4;
            const auto fb = a * xb + b * xb * xb / 2 + c * xb * xb * xb / 3 + d * xb * xb * xb * xb / 4;

            integrand += (fb - fa);
        }
        return integrand;
    }

    T operator()(T x) const
    {
        const Interval v { .x = x };
        const auto i1 = std::upper_bound(m_data.cbegin() + 1, m_data.cend() - 1, v, [](const auto& lh, const auto& rh) {
            return lh.x < rh.x;
        });
        const auto i0 = i1 - 1;
        const auto& seg = *i0;
        const auto xd = x - seg.x;

        return seg.a + seg.b * xd + seg.c * xd * xd + seg.d * xd * xd * xd;
    }

    void setup(std::vector<std::pair<T, T>> data_r)
    {
        std::sort(data_r.begin(), data_r.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        std::vector<std::pair<T, T>> data(N_KNOTS + 1);
        if (data_r.size() < N_KNOTS + 1) {
            std::copy(data_r.cbegin(), data_r.cend(), data.begin());
            AkimaSpline<T> sp(data_r);
            const auto N_rest = N_KNOTS + 1 - data_r.size();
            for (std::size_t i = 0; i < N_rest; ++i) {
                const auto idx = data_r.size() * (i + 1) / (N_rest + 1);
                const auto x = (data_r[idx].first + data_r[idx + 1].first) / 2;
                const auto y = sp(x);
                data[i + data_r.size()] = std::make_pair(x, y);
            }
            std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        } else if (data_r.size() > N_KNOTS + 1) {
            for (std::size_t i = 0; i < N_KNOTS + 1; ++i) {
                const auto frac = static_cast<T>(i) / N_KNOTS;
                const auto dIdx = static_cast<std::size_t>(frac * (data_r.size() - 1));
                data[i] = data_r[dIdx];
            }
            std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        } else {
            data = data_r;
        }

        std::array<T, N_KNOTS> m;
        for (std::size_t i = 0; i < N_KNOTS; ++i) {
            m[i] = (data[i + 1].second - data[i].second) / (data[i + 1].first - data[i].first);
        }
        std::array<T, N_KNOTS + 1> s;

        for (std::size_t i = 0; i < N_KNOTS + 1; ++i) {
            if (i == 0)
                s[i] = m[i];
            else if (i == 1)
                s[i] = (m[i - 1] + m[i]) / 2;
            else if (i == N_KNOTS - 1)
                s[i] = (m[i - 1] + m[i]) / 2;
            else if (i == N_KNOTS)
                s[i] = m[i - 1];
            else {
                const auto den = std::abs(m[i + 1] - m[i]) + std::abs(m[i - 1] - m[i - 2]);
                if (den > 0.0001)
                    s[i] = (std::abs(m[i + 1] - m[i]) * m[i - 1] + std::abs(m[i - 1] - m[i - 2]) * m[i]) / den;
                else
                    s[i] = (m[i - 1] + m[i]) / 2;
            }
        }

        for (std::size_t i = 0; i < N_KNOTS; ++i) {
            m_data[i].x = data[i].first;
            m_data[i].a = data[i].second;
            m_data[i].b = s[i];
            const auto xdiff = data[i + 1].first - data[i].first;
            m_data[i].c = (3 * m[i] - 2 * s[i] - s[i + 1]) / xdiff;
            m_data[i].d = (s[i] + s[i + 1] - 2 * m[i]) / (xdiff * xdiff);
        }
    }

private:
    struct Interval {
        T x = 0, a = 1, b = 0, c = 0, d = 0;
    };
    std::array<Interval, N_KNOTS> m_data;
};

template <Floating T>
class CubicSplineInterpolator {
public:
    CubicSplineInterpolator()
    {
        m_x.resize(1);
        m_x[0] = 0;
        m_coefficients.resize(3);
        for (auto& c : m_coefficients) {
            c[0] = 1;
            for (std::size_t i = 1; i < 4; ++i)
                c[i] = 0;
        }
    }

    CubicSplineInterpolator(const std::vector<T>& x, const std::vector<T>& y)
    {
        const auto N = std::min(x.size(), y.size());
        std::vector<std::pair<T, T>> data(N);
        for (std::size_t i = 0; i < N; ++i) {
            data[i].first = x[i];
            data[i].second = y[i];
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        setup(data);
    }

    CubicSplineInterpolator(std::vector<std::pair<T, T>> data)
    {
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        setup(data);
    }

    void setup(const std::vector<std::pair<T, T>>& data)
    {

        std::vector<T> y(data.size());
        m_x.resize(data.size());

        for (std::size_t i = 0; i < m_x.size(); ++i) {
            m_x[i] = data[i].first;
            y[i] = data[i].second;
        }

        std::vector<T> h(m_x.size());
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            h[i] = m_x[i + 1] - m_x[i];
        }

        std::vector<T> delta(m_x.size(), T { 1 });
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }
        std::vector<T> H(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            H[i] = 2 * (h[i - 1] + h[i]);
        }
        H[0] = H[1];

        std::vector<T> d(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            d[i] = 6 * (delta[i] - delta[i - 1]);
        }
        d[m_x.size() - 1] = 0;
        d[0] = 0;
        auto s = thomasPenSplineElimination(h, H, d);

        m_coefficients.resize(m_x.size() + 1);
        m_coefficients[0] = { y[0], 0, 0, 0 };
        m_coefficients[m_x.size()] = { y[m_x.size() - 1], 0, 0, 0 };
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            m_coefficients[i + 1][0] = (s[i] * m_x[i + 1] * m_x[i + 1] * m_x[i + 1] - s[i + 1] * m_x[i] * m_x[i] * m_x[i] + 6 * (y[i] * m_x[i + 1] - y[i + 1] * m_x[i])) / (6 * h[i]);
            m_coefficients[i + 1][0] += h[i] * (s[i + 1] * m_x[i] - s[i] * m_x[i + 1]) / 6;
            m_coefficients[i + 1][1] = (s[i + 1] * m_x[i] * m_x[i] - s[i] * m_x[i + 1] * m_x[i + 1] + 2 * (y[i + 1] - y[i])) / (2 * h[i]) + h[i] * (s[i] - s[i + 1]) / 6;
            m_coefficients[i + 1][2] = (s[i] * m_x[i + 1] - s[i + 1] * m_x[i]) / (2 * h[i]);
            m_coefficients[i + 1][3] = (s[i + 1] - s[i]) / (6 * h[i]);
        }
    }

    T operator()(const T x) const
    {
        auto it = std::upper_bound(m_x.cbegin(), m_x.cend(), x);
        const auto i = std::distance(m_x.cbegin(), it);
        return m_coefficients[i][0] + m_coefficients[i][1] * x + m_coefficients[i][2] * x * x + m_coefficients[i][3] * x * x * x;
    }

    T meanValue() const
    {
        T integrand = 0;
        // integrate part
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            const auto x0 = m_x[i - 1];
            const auto x1 = m_x[i];

            integrand += m_coefficients[i][0] * (x1 - x0);
            integrand += m_coefficients[i][1] * (x1 * x1 - x0 * x0) / 2;
            integrand += m_coefficients[i][2] * (x1 * x1 * x1 - x0 * x0 * x0) / 3;
            integrand += m_coefficients[i][3] * (x1 * x1 * x1 * x1 - x0 * x0 * x0 * x0) / 4;
        }
        return integrand / (m_x.back() - m_x.front());
    }

    void scale(T value)
    {
        for (auto& el : m_coefficients)
            for (auto& c : el)
                c *= value;
    }

protected:
    static std::vector<T> thomasPenSplineElimination(const std::vector<T>& h_p, const std::vector<T>& H_p, const std::vector<T>& d_p)
    {
        // Thomas algorithm for gaussian elimination for a trigonal system of equations
        /*
        |b0 c0  0 0  ..  | x0 |   |d0|
        |a1 b1 c1 0  ... | x1 | = |d1|
        |0  a2 b2 c2  ...| x2 |   |d2|
        |0  0  a3 b3 c3  | .. |   |..|

        */
        std::vector<T> h = h_p;
        std::vector<T> H = H_p;
        std::vector<T> d = d_p;

        std::vector<T> b(d.size(), T { 2 });
        for (std::size_t i = 1; i < d.size(); ++i) {
            const T w = h[i - 1] / H[i - 1];
            H[i] -= w * h[i - 1];
            d[i] -= w * d[i - 1];
        }
        std::vector<T> x(d.size());
        x[d.size() - 1] = d[d.size() - 1] / H[d.size() - 1];
        for (int i = d.size() - 2; i >= 0; --i) {
            x[i] = (d[i] - h[i] * x[i + 1]) / H[i];
        }
        return x;
    }

private:
    std::vector<std::array<T, 4>> m_coefficients;
    std::vector<T> m_x;
};

template <Floating T, int N = 30>
class CubicSplineInterpolatorStatic {
public:
    CubicSplineInterpolatorStatic(const std::vector<T>& x, const std::vector<T>& y)
    {
        const auto n_ele = std::min(x.size(), y.size());
        std::vector<std::pair<T, T>> data(n_ele);
        for (std::size_t i = 0; i < n_ele; ++i) {
            data[i].first = x[i];
            data[i].second = y[i];
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        setup(data);
    }

    CubicSplineInterpolatorStatic(std::vector<std::pair<T, T>> data)
    {
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
        setup(data);
    }

    template <std::regular_invocable<T> F>
        requires std::is_same<std::invoke_result_t<F, T>, T>::value
    CubicSplineInterpolatorStatic(const T start, const T stop, F function)
    {
        setup(start, stop, function);
    }

    T operator()(const T x_val) const
    {
        const T x = std::clamp(x_val, m_start, m_stop);

        const std::size_t index = x > m_start ? static_cast<std::size_t>((x - m_start) / m_step) : 0;
        const std::size_t offset = index < N - 1 ? index * 4 : (N - 2) * 4;
        return m_coefficients[offset] + m_coefficients[offset + 1] * x + m_coefficients[offset + 2] * x * x + m_coefficients[offset + 3] * x * x * x;
    }

protected:
    template <std::regular_invocable<T> F>
        requires std::is_same<std::invoke_result_t<F, T>, T>::value
    void setup(const T start, const T stop, F function)
    {
        std::vector<T> y(N);
        m_start = start;
        m_step = (stop - start) / (N - 2);
        for (std::size_t i = 0; i < N; ++i) {
            m_x[i] = m_start + m_step * i;
            y[i] = function(m_x[i]);
        }
        m_stop = m_x.back();

        std::vector<T> h(m_x.size());
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            h[i] = m_x[i + 1] - m_x[i];
        }

        std::vector<T> delta(m_x.size(), T { 1 });
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }
        std::vector<T> H(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            H[i] = 2 * (h[i - 1] + h[i]);
        }
        H[0] = H[1];

        std::vector<T> d(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            d[i] = 6 * (delta[i] - delta[i - 1]);
        }
        d[m_x.size() - 1] = 0;
        d[0] = 0;
        auto s = thomasPenSplineElimination(h, H, d);

        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            const auto offset = i * 4;

            m_coefficients[offset + 0] = (s[i] * m_x[i + 1] * m_x[i + 1] * m_x[i + 1] - s[i + 1] * m_x[i] * m_x[i] * m_x[i] + 6 * (y[i] * m_x[i + 1] - y[i + 1] * m_x[i])) / (6 * h[i]);
            m_coefficients[offset + 0] += h[i] * (s[i + 1] * m_x[i] - s[i] * m_x[i + 1]) / 6;
            m_coefficients[offset + 1] = (s[i + 1] * m_x[i] * m_x[i] - s[i] * m_x[i + 1] * m_x[i + 1] + 2 * (y[i + 1] - y[i])) / (2 * h[i]) + h[i] * (s[i] - s[i + 1]) / 6;
            m_coefficients[offset + 2] = (s[i] * m_x[i + 1] - s[i + 1] * m_x[i]) / (2 * h[i]);
            m_coefficients[offset + 3] = (s[i + 1] - s[i]) / (6 * h[i]);
        }
    }

    void setup(const std::vector<std::pair<T, T>>& data)
    {
        CubicSplineInterpolator f(data);
        const auto start = data.front().first;
        const auto stop = data.back().first;
        setup(start, stop, f);
    }

    static std::vector<T> thomasPenSplineElimination(const std::vector<T>& h_p, const std::vector<T>& H_p, const std::vector<T>& d_p)
    {
        // Thomas algorithm for gaussian elimination for a trigonal system of equations
        /*
        |b0 c0  0 0  ..  | x0 |   |d0|
        |a1 b1 c1 0  ... | x1 | = |d1|
        |0  a2 b2 c2  ...| x2 |   |d2|
        |0  0  a3 b3 c3  | .. |   |..|

        */
        std::vector<T> h = h_p;
        std::vector<T> H = H_p;
        std::vector<T> d = d_p;

        std::vector<T> b(d.size(), T { 2 });
        for (std::size_t i = 1; i < d.size(); ++i) {
            const T w = h[i - 1] / H[i - 1];
            H[i] -= w * h[i - 1];
            d[i] -= w * d[i - 1];
        }
        std::vector<T> x(d.size());
        x[d.size() - 1] = d[d.size() - 1] / H[d.size() - 1];
        for (int i = d.size() - 2; i >= 0; --i) {
            x[i] = (d[i] - h[i] * x[i + 1]) / H[i];
        }
        return x;
    }

private:
    std::array<T, (N - 1) * 4> m_coefficients;
    std::array<T, N> m_x;
    T m_step = 0;
    T m_start = 0;
    T m_stop = 0;
};

template <Floating T>
class Matrix {
public:
    Matrix() { };
    Matrix(std::size_t R, std::size_t C)
        : m_r(R)
        , m_c(C)
    {
        m_data.resize(m_r * m_c, T { 0 });
    }
    Matrix(std::size_t R, std::size_t C, const std::vector<T>& d)
        : m_r(R)
        , m_c(C)
    {
        m_data.resize(m_r * m_c);
        std::copy(d.cbegin(), d.cend(), m_data.begin());
    }

    T& operator()(std::size_t r, std::size_t c)
    {
        return m_data[index(r, c)];
    }
    const T& operator()(std::size_t r, std::size_t c) const
    {

        return m_data[index(r, c)];
    }
    std::vector<T> operator*(const std::vector<T>& x)
    {
        std::vector<T> y(x.size(), 0);
        for (std::size_t r = 0; r < m_r; ++r) {
            for (std::size_t c = 0; c < m_c; ++c) {
                y[r] += m_data[index(r, c)] * x[c];
            }
        }
        return y;
    }
    std::optional<std::vector<T>> solve(const std::vector<T>& bvalues) const
    {
        // gaussian ellim

        auto b = bvalues;
        Matrix<T> A(m_r, m_c, m_data);

        for (std::size_t k = 0; k < m_r; ++k) {
            auto imax = A.find_row_argmax(k);
            if (std::abs(A(imax, k)) <= std::numeric_limits<T>::epsilon()) {
                return std::nullopt;
            }
            A.swapRows(k, imax);
            std::swap(b[k], b[imax]);
            for (std::size_t i = k + 1; i < m_r; ++i) {
                const auto c = A(i, k) / A(k, k);
                A(i, k) = T { 0 };
                for (std::size_t j = k + 1; j < m_r; ++j) {
                    A(i, j) -= c * A(k, j);
                }
                b[i] -= c * b[k];
            }
        }

        // backsubstitution before return
        std::vector<T> x(b.size());
        for (std::size_t r = m_r - 1; r < m_r; --r) { // note we rely on unsigned wraparound
            T s = 0;
            for (std::size_t c = r + 1; c < m_r; ++c) {
                s += A(r, c) * b[c];
            }
            b[r] = (b[r] - s) / A(r, r);
        }
        return b;
    }

    void fill(const T value)
    {
        std::fill(m_data.begin(), m_data.end(), value);
    }
    void setZero()
    {
        fill(T { 0 });
    }
    void setZero(std::size_t r, std::size_t c)
    {
        m_r = r;
        m_c = c;
        m_data.resize(m_r * m_c);
        fill(T { 0 });
    }

    std::size_t rows() const { return m_r; }
    std::size_t cols() const { return m_c; }

    Matrix<T> transpose() const
    {
        Matrix<T> t(m_c, m_r);
        for (std::size_t i = 0; i < m_r; ++i) {
            for (std::size_t j = 0; j < m_c; ++j) {
                t(j, i) = m_data[index(i, j)];
            }
        }
        return t;
    }

protected:
    void swapRows(std::size_t from, std::size_t to)
    {
        if (from == to)
            return;
        for (std::size_t i = 0; i < m_c; ++i) {
            std::swap(m_data[index(from, i)], m_data[index(to, i)]);
        }
    }

    std::size_t index(std::size_t r, std::size_t c) const
    {
        return r * m_c + c;
    }

    std::size_t find_row_argmax(std::size_t r) const
    {
        auto imax = r;
        T max_pivot = std::abs(m_data[index(r, r)]);
        for (std::size_t i = r + 1; i < m_r; ++i) {
            const T a = std::abs(m_data[index(i, r)]);
            if (a > max_pivot) {
                max_pivot = a;
                imax = i;
            }
        }
        return imax;
    }

private:
    std::vector<T> m_data;
    std::size_t m_r = 0;
    std::size_t m_c = 0;
};

template <Floating T>
class CubicLSInterpolator {
public:
    CubicLSInterpolator(const std::vector<std::pair<T, T>>& data, std::size_t n_knots = 30, bool maybe_discont = true)
    {
        std::vector<T> x, y;
        x.reserve(data.size());
        y.reserve(data.size());
        for (const auto& p : data) {
            x.push_back(p.first);
            y.push_back(p.second);
        }
        setupSplines(x, y, n_knots, maybe_discont);
    }
    CubicLSInterpolator(const std::vector<T>& x, const std::vector<T>& y, std::size_t n_knots = 30, bool maybe_discont = true)
    {
        setupSplines(x, y, n_knots, maybe_discont);
    }

    CubicLSInterpolator(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& t, bool maybe_discont = true)
    {
        setupSplines(x, y, t, maybe_discont);
    }

    T operator()(const T x) const
    {
        return evaluateSpline(x, m_data);
    }

    const std::vector<std::array<T, 3>>& getDataTable() const { return m_data; }

    static inline T evaluateSpline(const T x, const std::vector<std::array<T, 3>>& data)
    {
        return evaluateSpline(x, data.cbegin(), data.cend());
    }

    template <std::random_access_iterator It>
        requires std::is_same_v<typename std::iterator_traits<It>::value_type, std::array<T, 3>>
    static inline T evaluateSpline(const T x, const It begin, const It end)
    {
        const std::array<T, 3> x_comp { x, 0, 0 };
        auto it = std::upper_bound(begin + 1, end - 1, x_comp, [](const auto& left, const auto& right) -> bool { return left[0] < right[0]; });
        const auto& d1 = *it;
        const auto& d0 = *(--it);

        const auto h = d1[0] - d0[0];
        const auto u = (x - d0[0]) / h;
        const auto v = u - 1;
        const auto uu = u * u;
        const auto vv = v * v;

        const auto f1 = (2 * u + 1) * vv * d0[1];
        const auto f2 = uu * (1 - 2 * v) * d1[1];
        const auto f3 = h * u * vv * d0[2];
        const auto f4 = h * uu * v * d1[2];
        return f1 + f2 + f3 + f4;
    }

protected:
    struct Spline {
        Spline(const std::vector<T>& t, const std::vector<T>& p, const std::vector<T>& pz)
        {
            m_data.reserve(t.size());
            for (std::size_t i = 0; i < t.size(); ++i) {
                std::array<T, 3> val { t[i], p[i], pz[i] };
                m_data.push_back(val);
            }
        }
        Spline() { };
        Spline(const std::vector<std::array<T, 3>>& data)
            : m_data(data)
        {
        }
        Spline operator+(const Spline& right) const
        {
            Spline res(m_data);
            res.m_data.insert(res.m_data.end(), right.m_data.cbegin(), right.m_data.cend());
            return res;
        }
        void operator+=(const Spline& right)
        {
            m_data.insert(m_data.end(), right.m_data.cbegin(), right.m_data.cend());
        }
        std::vector<std::array<T, 3>> m_data; // m_t, m_z, m_zp
    };

    template <typename U>
        requires(std::is_convertible<U, std::vector<T>>::value || std::is_convertible<U, std::size_t>::value)
    void setupSplines(const std::vector<T>& x, const std::vector<T>& y, const U& t, bool maybe_discont)
    {
        Spline spline;
        if (maybe_discont) {
            // all this to handle dicontinous functions designated by a equal x value
            std::size_t x_start = 0;
            for (std::size_t i = 0; i < x.size() - 1; ++i) {
                constexpr auto epsilon = std::numeric_limits<T>::epsilon();
                if (x[i + 1] - x[i] <= epsilon) {
                    // ups, we have a discontinous function
                    auto x_beg = x_start;
                    x_start = i + 1;

                    std::vector<T> x_red, y_red;
                    for (std::size_t j = x_beg; j <= i; ++j) {
                        x_red.push_back(x[j]);
                        y_red.push_back(y[j]);
                    }
                    if constexpr (std::is_convertible<U, std::size_t>::value) {
                        auto spline_opt = calculateLSSplinePart(x_red, y_red, t);
                        if (spline_opt)
                            spline += spline_opt.value();
                    } else {
                        std::vector<T> t_red;
                        t_red.push_back(x[x_beg]);
                        for (std::size_t j = 0; j < t.size(); ++j) {
                            if ((t[j] < x[i]) && (t[j] > x[x_beg])) {
                                t_red.push_back(t[j]);
                            }
                        }
                        t_red.push_back(x[i]);
                        spline += calculateLSSplinePart(x_red, y_red, t_red);
                    }
                }
            }
            // did we se any discontinuities?
            if (x_start == 0) {
                auto spline_opt = calculateLSSplinePart(x, y, t);
                if (spline_opt)
                    spline += spline_opt.value();
            } else {
                std::vector<T> x_red, y_red;
                for (std::size_t j = x_start; j < x.size(); ++j) {
                    x_red.push_back(x[j]);
                    y_red.push_back(y[j]);
                }
                if constexpr (std::is_convertible<U, std::size_t>::value) {
                    auto spline_opt = calculateLSSplinePart(x_red, y_red, t);
                    if (spline_opt)
                        spline += spline_opt.value();
                } else {
                    std::vector<T> t_red;
                    t_red.push_back(x[x_start]);
                    for (std::size_t j = 0; j < t.size(); ++j) {
                        if (t[j] > x[x_start]) {
                            t_red.push_back(t[j]);
                        }
                    }
                    spline += calculateLSSplinePart(x_red, y_red, t_red);
                }
            }
        } else {
            auto spline_opt = calculateLSSplinePart(x, y, t);
            if (spline_opt)
                spline += spline_opt.value();
        }

        m_data = spline.m_data;
    }

    std::optional<Spline> calculateLSSplinePart(const std::vector<T>& x, const std::vector<T>& y, std::size_t N = 0) const
    {
        const auto Nlim = (x.size() - 1) / 3;
        if (N <= 1)
            N = Nlim;

        N = std::min(N, Nlim);

        if (x.size() < 4 || N < 2) {
            return std::nullopt;
        }

        std::vector<T> t(N);
        t[0] = x.front();
        t[N - 1] = x.back();
        for (std::size_t i = 1; i < N - 1; i++) {
            t[i] = x.front() + (i * (x.back() - x.front())) / (N - 1);
        }

        return calculateLSSplinePart(x, y, t);
    }

    std::optional<Spline> calculateLSSplinePart(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& t) const
    { // assume x is sorted
        const std::size_t N = t.size() - 1;

        std::vector<T> h(N);
        for (std::size_t i = 0; i < N; ++i)
            h[i] = t[i + 1] - t[i];

        // deal with t intervals without data
        std::vector<std::size_t> pn(N, 0);
        std::vector<std::size_t> p(N, 0);
        std::size_t j = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i] <= t[j + 1]) {
                pn[j]++;
            } else {
                j++;
                p[j] = i;
                pn[j] = 0;
                while (x[i] > t[j + 1]) {
                    j++;
                    p[j] = i;
                    pn[j] = 0;
                }
                pn[j] = 1;
            }
        }

        // setting up matrices
        std::vector<T> taa(N, T { 0 });
        std::vector<T> tab(N, T { 0 });
        std::vector<T> tag(N, T { 0 });
        std::vector<T> tad(N, T { 0 });
        std::vector<T> tbb(N, T { 0 });
        std::vector<T> tbg(N, T { 0 });
        std::vector<T> tbd(N, T { 0 });
        std::vector<T> tgg(N, T { 0 });
        std::vector<T> tgd(N, T { 0 });
        std::vector<T> tdd(N, T { 0 });

        std::vector<T> D(N + 1, T { 0 });
        std::vector<T> G(N + 1, T { 0 });

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = p[i]; j < p[i] + pn[i]; ++j) {
                const auto u = (x[j] - t[i]) / h[i];
                const auto v = u - 1;
                const auto alpha = (2 * u + 1) * v * v;
                const auto beta = u * u * (1 - 2 * v);
                const auto gamma = h[i] * u * v * v;
                const auto delta = h[i] * u * u * v;
                taa[i] += alpha * alpha;
                tab[i] += alpha * beta;
                tag[i] += alpha * gamma;
                tad[i] += alpha * delta;
                tbb[i] += beta * beta;
                tbg[i] += beta * gamma;
                tbd[i] += beta * delta;
                tgg[i] += gamma * gamma;
                tgd[i] += gamma * delta;
                tdd[i] += delta * delta;
                D[i] += 2 * y[j] * alpha;
                D[i + 1] += 2 * y[j] * beta;
                G[i] += 2 * y[j] * gamma;
                G[i + 1] += 2 * y[j] * delta;
            }
        }

        using Matrix = Matrix<T>;

        Matrix A, B, E, C, F, Z;
        A.setZero(N + 1, N + 1);
        B.setZero(N + 1, N + 1);
        E.setZero(N + 1, N + 1);
        C.setZero(N - 1, N + 1);
        F.setZero(N - 1, N + 1);
        Z.setZero(N - 1, N + 1);

        for (std::size_t i = 0; i < N; ++i) {
            A(i, i) += 2 * taa[i];
            A(i + 1, i + 1) = 2 * tbb[i];
            A(i, i + 1) = 2 * tab[i];
            A(i + 1, i) = A(i, i + 1);

            B(i, i) += 2 * tag[i];
            B(i + 1, i + 1) = 2 * tbd[i];
            B(i, i + 1) = 2 * tbg[i];
            B(i + 1, i) = 2 * tad[i];

            E(i, i) += 2 * tgg[i];
            E(i + 1, i + 1) = 2 * tdd[i];
            E(i, i + 1) = 2 * tgd[i];
            E(i + 1, i) = E(i, i + 1);
        }
        for (std::size_t i = 0; i < N - 1; ++i) {
            C(i, i) = 3 * h[i + 1] / h[i];
            C(i, i + 2) = -3 * h[i] / h[i + 1];
            C(i, i + 1) = -C(i, i) - C(i, i + 2);

            F(i, i) = h[i + 1];
            F(i, i + 1) = 2 * (h[i] + h[i + 1]);
            F(i, i + 2) = h[i];
        }

        auto BT = B.transpose();
        auto CT = C.transpose();
        auto FT = F.transpose();

        Matrix sol(A.rows() + BT.rows() + C.rows(), A.cols() + BT.cols() + CT.cols());
        sol.setZero();

        for (std::size_t r = 0; r < sol.rows(); ++r) {
            for (std::size_t c = 0; c < sol.cols(); ++c) {
                if ((r < N + 1) && (c < N + 1)) {
                    sol(r, c) = A(r, c);
                } else if ((r < N + 1) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = BT(r, c - (N + 1));
                } else if ((r < N + 1) && (c >= 2 * (N + 1))) {
                    sol(r, c) = CT(r, c - (2 * (N + 1)));
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c < N + 1)) {
                    sol(r, c) = B(r - (N + 1), c);
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = E(r - (N + 1), c - (N + 1));
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c >= 2 * (N + 1))) {
                    sol(r, c) = FT(r - (N + 1), c - 2 * (N + 1));
                } else if ((r >= 2 * (N + 1)) && (c < N + 1)) {
                    sol(r, c) = C(r - 2 * (N + 1), c);
                } else if ((r >= 2 * (N + 1)) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = F(r - 2 * (N + 1), c - (N + 1));
                } else if ((r >= 2 * (N + 1)) && (c >= 2 * (N + 1))) {
                    sol(r, c) = 0;
                }
            }
        }

        std::vector<T> bval(sol.rows(), T { 0 });

        std::copy(D.begin(), D.end(), bval.begin());
        std::copy(G.begin(), G.end(), bval.begin() + N + 1);

        auto res_opt = sol.solve(bval);
        if (!res_opt)
            return std::nullopt;
        auto& res = res_opt.value();
        std::vector<T> m_p(N + 1);
        std::vector<T> m_pz(N + 1);

        std::copy(res.begin(), res.begin() + N + 1, m_p.begin());
        std::copy(res.begin() + N + 1, res.begin() + 2 * (N + 1), m_pz.begin());
        Spline spline(t, m_p, m_pz);
        return std::make_optional(spline);
    }

private:
    std::vector<std::array<T, 3>> m_data; // array of m_t, m_z, m_zp
};

}
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

#pragma once

#include "xraymc/constants.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/vectormath.hpp"

#include <array>
#include <vector>

namespace xraymc {

class CTAECFilter {
public:
    CTAECFilter()
    {
        m_data.resize(2);
        m_data[0] = std::make_pair(0.0, 1.0);
        m_data[1] = std::make_pair(1.0, 1.0);
    }

    CTAECFilter(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& values)
    {
        setData(start, stop, values);
    }

    void normalizeBetween(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        if (!isEmpty()) {
            normalize(start, stop);
        }
    }

    std::size_t size() const
    {
        return m_data.size();
    }

    bool isEmpty() const
    {
        return m_data.size() == 2;
    }

    std::vector<double> weights() const
    {
        std::vector<double> w(m_data.size());
        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), w.begin(), [](const auto& v) { return v.second; });
        return w;
    }

    const std::array<double, 3>& start() const
    {
        return m_start;
    }

    std::array<double, 3> stop() const
    {
        return xraymc::vectormath::add(m_start, vectormath::scale(m_dir, m_length));
    }

    double length() const
    {
        return m_length;
    }

    void setData(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& data)
    {
        m_start = start;
        const auto dir = vectormath::subtract(stop, start);
        if (vectormath::length_sqr(dir) < GEOMETRIC_ERROR()) {
            m_dir = { 0, 0, 1 };
            m_length = 1;
        } else {
            m_dir = vectormath::normalized(dir);
            m_length = vectormath::length(dir);
        }

        if (data.size() < 2) {
            m_data.resize(2);
            m_data[0] = std::make_pair(0.0, 1.0);
            m_data[1] = std::make_pair(1.0, 1.0);
            return;
        }

        m_data.resize(data.size());
        const auto step = m_length / data.size();
        for (std::size_t i = 0; i < data.size(); ++i) {
            m_data[i].first = step * i;
            m_data[i].second = data[i];
        }
        normalize();
    }

    double operator()(const std::array<double, 3>& pos) const
    {
        return this->operator()(positionToIndex(pos));
    }

    double operator()(double d) const
    {
        return interpolate(m_data, d);
    }

    double integrate() const
    {
        return trapz(m_data);
    }

    double integrate(const std::array<double, 3>& start_p, const std::array<double, 3>& stop_p) const
    {
        const auto start = positionToIndex(start_p);
        const auto stop = positionToIndex(stop_p);
        return trapz(m_data, start, stop);
    }

protected:
    double positionToIndex(const std::array<double, 3>& start) const
    {
        const auto dist_start = vectormath::subtract(start, m_start);
        return std::clamp(vectormath::dot(dist_start, m_dir), 0.0, m_length);
    }

    void normalize()
    {
        const auto area = integrate();
        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto k = m_length / area;
        for (auto& d : m_data)
            d.second *= k;
    }

    void normalize(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        const auto p_start = positionToIndex(start);
        const auto p_stop = positionToIndex(stop);

        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto area = integrate(start, stop);
        const auto k = std::abs(p_stop - p_start) / area;
        for (auto& d : m_data)
            d.second *= k;
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    double m_length = 1;
    std::vector<std::pair<double, double>> m_data;
};
}
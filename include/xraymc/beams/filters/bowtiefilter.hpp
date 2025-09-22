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

#include "xraymc/interpolation.hpp"

#include <array>

namespace xraymc {

class BowtieFilter {
public:
    template <bool ONESIDED = true>
    BowtieFilter(const std::vector<double>& angles_r, const std::vector<double>& intensity_r)
    {
        setData<ONESIDED>(angles_r, intensity_r);
    }
    template <bool ONESIDED = true>
    BowtieFilter(const std::vector<std::pair<double, double>>& data)
    {
        setData<ONESIDED>(data);
    }

    BowtieFilter()
    {
        // generic filter from a Siemens Definition Flash
        static std::vector<std::pair<double, double>> data = {
            { 0.166511074, 3.53208 },
            { 0.000000000, 13.9167 },
            { 0.041992107, 12.5868 },
            { 0.083836642, 9.41943 },
            { 0.246954945, 1.96665 },
            { 0.324269441, 1.27605 },
            { 0.390607044, 0.947716 }
        };
        setData<true>(data);
    }

    template <bool ONESIDED = true>
    double operator()(double angle) const
    {
        if constexpr (ONESIDED)
            return m_inter(std::abs(angle));
        else
            return m_inter(angle);
    }

    template <bool ONESIDED = true>
    void setData(std::vector<std::pair<double, double>> data)
    {
        for (auto& d : data) {
            if constexpr (ONESIDED)
                d.first = std::abs(d.first);
            d.second = std::abs(d.second);
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        m_inter.setup(data);

        const auto start = data.front().first;
        const auto stop = data.back().first;
        const auto area = m_inter.integral(start, stop);

        m_inter.scale((stop - start) / area);
    }

    template <bool ONESIDED = true>
    void setData(const std::vector<double>& angles_r, const std::vector<double>& intensity_r)
    {
        const auto N = std::min(angles_r.size(), intensity_r.size());
        std::vector<std::pair<double, double>> data(N);

        for (std::size_t i = 0; i < N; ++i) {
            data[i].first = angles_r[i];
            data[i].second = std::abs(intensity_r[i]);
        }

        setData<ONESIDED>(data);
    }

private:
    AkimaSplineStatic<double, 6> m_inter;
};
}
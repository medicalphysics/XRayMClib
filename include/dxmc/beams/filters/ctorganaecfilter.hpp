/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/vectormath.hpp"

#include <array>
#include <numbers>
#include <vector>

namespace dxmc {

class CTOrganAECFilter {
public:
    CTOrganAECFilter()
    {
    }

    CTOrganAECFilter(double low_weight, double start_angle, double stop_angle, double ramp_angle = 0)
    {
        setLowWeightFactor(low_weight);
        setStartStopAngles(start_angle, stop_angle);
        setRampAngle(ramp_angle);
    }

    double operator()(double angle) const
    {
        const auto nang = normalize_angle(angle);
        if (m_start_angle < nang && nang < m_stop_angle) {
            // low region
            return m_weight_factor;
        } else if (m_start_angle - m_ramp_angle < nang && nang < m_stop_angle + m_ramp_angle) {
            // ramp region
            const auto delta = (maxWeight() - m_weight_factor) / m_ramp_angle;
            const auto dang = nang > m_stop_angle ? nang - m_stop_angle : m_start_angle - nang;
            return dang * delta + m_weight_factor;
        } else {
            // high region
            return maxWeight();
        }
    }

    void setUseFilter(bool on)
    {
        m_useFilter = on;
    }

    bool useFilter() const { return m_useFilter; }

    void setCompensateOutside(bool on)
    {
        m_compensate_outside = on;
    }

    bool compensateOutside() const
    {
        return m_compensate_outside;
    }

    void setStartAngle(double ang)
    {
        m_start_angle = normalize_angle(ang);
    }
    void setStartAngleDeg(double ang)
    {
        m_start_angle = normalize_angle(ang * DEG_TO_RAD<double>());
    }

    void setStopAngle(double ang)
    {
        m_stop_angle = normalize_angle(ang);
    }
    void setStopAngleDeg(double ang)
    {
        m_stop_angle = normalize_angle(ang * DEG_TO_RAD<double>());
    }
    double stopAngle() const { return m_stop_angle; }
    double startAngle() const { return m_start_angle; }
    double stopAngleDeg() const { return m_stop_angle * RAD_TO_DEG<double>(); }
    double startAngleDeg() const { return m_start_angle * RAD_TO_DEG<double>(); }

    void setStartStopAngles(double min, double max)
    {
        setStartAngle(min);
        setStopAngle(max);
    }

    void setRampAngle(double ang)
    {
        m_ramp_angle = std::clamp(std::abs(ang), 0.0, std::numbers::pi_v<double> / 4.0);
    }
    void setRampAngleDeg(double ang) { setRampAngle(ang * DEG_TO_RAD<double>()); }
    double rampAngle() const { return m_ramp_angle; }
    double rampAngleDeg() const { return m_ramp_angle * RAD_TO_DEG<double>(); }

    void setLowWeightFactor(double w)
    {
        m_weight_factor = std::clamp(std::abs(w), 0.001, 1.0);
    }
    double lowWeight() const { return m_weight_factor; }

    double maxWeight() const
    {
        if (m_compensate_outside) {
            const auto al = std::abs(m_start_angle - m_stop_angle);
            const auto ar = m_ramp_angle;
            const auto ah = 2 * std::numbers::pi_v<double> - 2 * ar - al;
            return (std::numbers::pi_v<double> * 2 - m_weight_factor * (al + ar)) / (ah + ar);
        } else {
            return 1.0;
        }
    }

protected:
    static double normalize_angle(double angle)
    {
        while (angle < std::numbers::pi_v<double>) {
            angle += std::numbers::pi_v<double> * 2;
        }
        while (angle > std::numbers::pi_v<double>) {
            angle -= std::numbers::pi_v<double> * 2;
        }
        return angle;
    }

private:
    double m_start_angle = 0;
    double m_stop_angle = std::numbers::pi_v<double> / 2;
    double m_ramp_angle = 0;
    double m_weight_factor = 0.6;
    bool m_compensate_outside = true;
    bool m_useFilter = false;
};
}
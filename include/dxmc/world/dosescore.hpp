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

Copyright 2023 Erlend Andersen
*/

#pragma once
#include "dxmc/world/energyscore.hpp"
#include <atomic>

namespace dxmc {

class DoseScore {
public:
    DoseScore() { }

    void addScoredEnergy(const EnergyScore& energy, double volume, double density, double calibrationfactor = 1)
    {
        const auto mass = volume * density; // grams
        const auto factor = calibrationfactor / mass;
        m_dose += energy.energyImparted() * factor;
        m_doseVariance += energy.variance() * factor * factor;
        m_nEvents += energy.numberOfEvents();
    }

    auto dose() const
    {
        return m_dose;
    }

    auto variance() const
    {
        return m_doseVariance;
    }

    auto standardDeviation() const
    {
        return std::sqrt(variance());
    }

    auto relativeUncertainty() const
    {
        return m_nEvents > 0 ? 1.96 * standardDeviation() / dose() : 0.0;
    }

    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

    void clear()
    {
        m_nEvents = 0;
        m_dose = 0;
        m_doseVariance = 0;
    }

private:
    double m_dose = 0;
    double m_doseVariance = 0;
    std::uint64_t m_nEvents = 0;
};
}

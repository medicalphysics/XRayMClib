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
#include "xraymc/world/energyscore.hpp"
#include <atomic>

namespace xraymc {

/**
 * @brief Accumulates absorbed dose and its statistical uncertainty from EnergyScore data.
 *
 * Converts energy imparted (from an EnergyScore) to absorbed dose by dividing by the
 * mass of the scoring volume (volume × density). Variance is propagated correctly so
 * that relativeUncertainty() returns the 95 % confidence interval half-width as a
 * fraction of the mean dose.
 *
 * Units:
 * - Dose:     eV/g  (energy in eV divided by mass in g)
 * - Variance: (eV/g)²
 */
class DoseScore {
public:
    /// @brief Default constructor; initialises all accumulators to zero.
    DoseScore() { }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const DoseScore&) const = default;

    /**
     * @brief Converts an EnergyScore to dose and accumulates it.
     *
     * Dose increment = energyImparted × calibrationfactor / (volume × density).
     * Variance is scaled by the square of the same factor.
     * @param energy            EnergyScore holding the accumulated energy and its variance.
     * @param volume            Scoring volume in cm³.
     * @param density           Material density in g/cm³.
     * @param calibrationfactor Optional multiplicative calibration factor (default 1).
     */
    void addScoredEnergy(const EnergyScore& energy, double volume, double density, double calibrationfactor = 1)
    {
        const auto mass = volume * density; // grams
        const auto factor = calibrationfactor / mass;
        m_dose += energy.energyImparted() * factor;
        m_doseVariance += energy.variance() * factor * factor;
        m_nEvents += energy.numberOfEvents();
    }

    /// @brief Returns the accumulated mean dose in eV/g.
    auto dose() const
    {
        return m_dose;
    }

    /// @brief Returns the accumulated dose variance in (eV/g)².
    auto variance() const
    {
        return m_doseVariance;
    }

    /// @brief Returns the standard deviation of the dose in eV/g.
    auto standardDeviation() const
    {
        return std::sqrt(variance());
    }

    /**
     * @brief Returns the 95 % confidence interval half-width as a fraction of the mean dose.
     *
     * Computed as 1.96 × standardDeviation / dose. Returns 0 when no events have been scored.
     */
    auto relativeUncertainty() const
    {
        return m_nEvents > 0 ? 1.96 * standardDeviation() / dose() : 0.0;
    }

    /// @brief Returns the total number of interaction events contributing to the score.
    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

    /// @brief Resets all accumulators to zero.
    void clear()
    {
        m_nEvents = 0;
        m_dose = 0;
        m_doseVariance = 0;
    }

    /**
     * @brief Directly sets the dose, variance, and event count (used during deserialization).
     * @param dose      Mean dose in eV/g.
     * @param variance  Dose variance in (eV/g)².
     * @param n_events  Number of scored events.
     */
    void set(double dose, double variance, std::uint64_t n_events)
    {
        m_dose = dose;
        m_doseVariance = variance;
        m_nEvents = n_events;
    }

private:
    double m_dose = 0;
    double m_doseVariance = 0;
    std::uint64_t m_nEvents = 0;
};
}

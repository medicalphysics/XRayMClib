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

#pragma once
#include <atomic>
#include <cmath>

namespace xraymc {

/**
 * @brief Thread-safe accumulator for energy imparted and its statistical variance.
 *
 * Records the sum of deposited energies and the sum of their squares using
 * `std::atomic_ref` so that multiple transport threads can score concurrently
 * without external synchronisation. The stored sums are used to estimate the
 * variance of the total energy via the law of total variance, assuming that the
 * number of interaction events follows a Poisson distribution.
 *
 * Units: all energies in eV.
 */
class EnergyScore {
public:
    /// @brief Default constructor; initialises all accumulators to zero.
    EnergyScore() { }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const EnergyScore&) const = default;

    /**
     * @brief Records a single energy-deposition event in a thread-safe manner.
     *
     * Adds @p energy to the running sum and @p energy² to the sum-of-squares,
     * and increments the event counter. Values ≤ 0 are silently ignored.
     * All three updates are performed with relaxed memory order via `std::atomic_ref`.
     * @param energy Energy deposited in this event in eV; must be > 0 to be recorded.
     */
    void scoreEnergy(double energy)
    {
        if (energy <= 0.0)
            return;

        {
            auto aref = std::atomic_ref(m_energyImparted);
            aref.fetch_add(energy, std::memory_order_relaxed);
        }
        {
            auto aref = std::atomic_ref(m_energyImpartedSquared);
            aref.fetch_add(energy * energy, std::memory_order_relaxed);
        }
        {
            auto aref = std::atomic_ref(m_nEvents);
            ++aref;
        }
    }

    /// @brief Returns the total energy imparted (sum of all scored energies) in eV.
    double energyImparted() const
    {
        return m_energyImparted;
    }

    /// @brief Returns the sum of squared scored energies in eV².
    double energyImpartedSquared() const
    {
        return m_energyImpartedSquared;
    }

    /**
     * @brief Estimates the variance of the total energy imparted in eV².
     *
     * Applies the law of total variance for a random sum Sₙ = X₁ + … + Xₙ where
     * the number of events N is Poisson-distributed with mean and variance equal to
     * the observed event count Z:
     *
     *   Var(Sₙ) = Z · Var(X) + E[X]² · Z
     *
     * where E[X] = energyImparted / Z and
     * Var(X) = (E[X²] − E[X]²) / (Z − 1).
     *
     * Returns 0 when fewer than two events have been scored.
     */
    double variance() const
    {
        // We will calculate the the variance of a random sum Sn = X1+X2+...+Xn
        // where N is random integer larger than 0
        // X1, ..., Xn is a random variable that is normal distributed

        // The law of total variance states that
        // Var(Sn) = E[N] * Var(X) + E(X)*E(X) * Var(N)
        // let number of interaction  events be Z
        // Assuming N is poisson distributed with E[N] = Z, and
        // Var(N) = Z
        // Var(Sn) = Z*Var(X) + E(X)*E(X) * Z
        // E(X) is the expectation value of X, i.e (energy imparted) / (number of events)
        // Var(X) is the variance of X, Var(X) = ( E(X**2) - E(X)*E(X) ) / ( N-1 )
        if (numberOfEvents() > 1) {

            // expected energy per event
            const auto e = energyImparted() / numberOfEvents();

            // expected squared energy per event
            const auto e2 = energyImpartedSquared() / numberOfEvents();

            // variance of X per event
            const auto var_per_event = (e2 - e * e) / (numberOfEvents() - 1);

            // variance of sum of events
            auto var = numberOfEvents() * var_per_event + e * e * numberOfEvents();

            return var;

            // For a fixed number of N the variance is simply the Var(X)*N*N
            // Var(Sn) = N*(e2-e*e)/(N-1)
        }
        return 0;
    }

    /// @brief Returns the standard deviation of the total energy imparted in eV.
    double standardDeviation() const
    {
        return std::sqrt(variance());
    }

    /**
     * @brief Returns the 95 % confidence interval half-width as a fraction of the mean energy.
     *
     * Computed as 1.96 × standardDeviation / energyImparted.
     */
    double relativeUncertainty() const
    {
        return 1.96 * standardDeviation() / energyImparted();
    }

    /// @brief Returns the total number of scored interaction events.
    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

    /**
     * @brief Directly sets all accumulators (used during deserialization).
     * @param energy        Total energy imparted in eV.
     * @param energySquared Sum of squared energies in eV².
     * @param n_events      Number of scored events.
     */
    void set(double energy, double energySquared, std::uint64_t n_events)
    {
        m_nEvents = n_events;
        m_energyImparted = energy;
        m_energyImpartedSquared = energySquared;
    }

    /// @brief Resets all accumulators to zero.
    void clear()
    {
        m_nEvents = 0;
        m_energyImparted = 0;
        m_energyImpartedSquared = 0;
    }

private:
    std::uint64_t m_nEvents = 0;
    double m_energyImparted = 0;
    double m_energyImpartedSquared = 0;
};
}

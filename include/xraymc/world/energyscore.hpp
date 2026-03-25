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

class EnergyScore {
public:
    EnergyScore() { }
    bool operator==(const EnergyScore&) const = default;

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
    double energyImparted() const
    {
        return m_energyImparted;
    }

    double energyImpartedSquared() const
    {
        return m_energyImpartedSquared;
    }

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

    double standardDeviation() const
    {
        return std::sqrt(variance());
    }

    double relativeUncertainty() const
    {
        return 1.96 * standardDeviation() / energyImparted();
    }

    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

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

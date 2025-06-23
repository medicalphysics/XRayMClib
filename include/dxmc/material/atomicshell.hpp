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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include <concepts>
#include <cstdint>
#include <utility>
#include <vector>

namespace dxmc {

struct AtomicShellRadiativeEmission {    
    std::uint64_t vacancy = 0; // shell from where the vacancy is filled
    double energy = 0;
    double probability = 0;
};

struct AtomicShell {
    AtomicShell(std::uint64_t shell = 0)
        : shell(shell)
    {
    }
    std::uint64_t shell = 0;
    double numberOfElectrons = 0;
    double bindingEnergy = 0;
    double kineticEnergy = 0;
    double HartreeFockOrbital_0 = 0;
    double numberOfPhotonsPerInitVacancy = 0;
    double energyOfPhotonsPerInitVacancy = 0;
    std::vector<std::pair<double, double>> photoel;
    std::vector<AtomicShellRadiativeEmission> radiativeEmissions;
};

}
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

#include "atomicshell.hpp"

#include <concepts>
#include <cstdint>
#include <map>
#include <vector>

namespace dxmc {

struct AtomicElement {
    std::uint64_t Z = 0;
    double atomicWeight = 0;
    double standardDensity = 0;
    std::vector<std::pair<double, double>> coherent;
    std::vector<std::pair<double, double>> incoherent;
    std::vector<std::pair<double, double>> photoel;
    std::vector<std::pair<double, double>> formFactor;
    std::vector<std::pair<double, double>> incoherentSF;
    std::vector<std::pair<double, double>> incoherentMeanScatterEnergy;
    std::map<std::uint64_t, AtomicShell> shells;
};

}
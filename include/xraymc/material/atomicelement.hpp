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

#include "atomicshell.hpp"

#include <concepts>
#include <cstdint>
#include <map>
#include <vector>

namespace xraymc {

/**
 * @brief Tabulated photon-interaction data for a single chemical element.
 *
 * Stores the EPICS cross-section tables needed for analog photon transport:
 * coherent (Rayleigh) scattering, incoherent (Compton) scattering, photoelectric
 * absorption, form factors, and incoherent scattering functions. Each table is a
 * vector of `{energy [keV], cross-section or function value}` pairs sampled on a
 * log-spaced energy grid.
 *
 * The `shells` map provides per-subshell binding energies and fluorescence yields
 * for photoelectric edge sampling, keyed by subshell identifier.
 */
struct AtomicElement {
    std::uint8_t Z = 0;                                                         ///< Atomic number.
    double atomicWeight = 0;                                                     ///< Atomic weight [g/mol].
    double standardDensity = 0;                                                  ///< Standard elemental density [g/cm³].
    std::vector<std::pair<double, double>> coherent;                             ///< Coherent (Rayleigh) scattering cross-section table: {energy [keV], μ/ρ [cm²/g]}.
    std::vector<std::pair<double, double>> incoherent;                           ///< Incoherent (Compton) scattering cross-section table: {energy [keV], μ/ρ [cm²/g]}.
    std::vector<std::pair<double, double>> photoel;                              ///< Photoelectric absorption cross-section table: {energy [keV], μ/ρ [cm²/g]}.
    std::vector<std::pair<double, double>> formFactor;                           ///< Coherent scattering form factor table: {momentum transfer x, F(x)}.
    std::vector<std::pair<double, double>> incoherentSF;                        ///< Incoherent scattering function table: {momentum transfer x, S(x)}.
    std::vector<std::pair<double, double>> incoherentMeanScatterEnergy;         ///< Mean scattered photon energy for Compton interactions: {energy [keV], mean scatter energy [keV]}.
    std::map<std::uint64_t, AtomicShell> shells;                                ///< Per-subshell data (binding energy, fluorescence yield) keyed by subshell identifier.
};

}
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

#include <concepts>
#include <cstdint>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief A single radiative (fluorescence) transition that fills a vacancy in an atomic subshell.
 *
 * When a photoelectric vacancy is created in an inner shell, it may be filled by an
 * electron from a higher shell with emission of a characteristic X-ray photon. This
 * struct records the source shell, photon energy, and relative probability for one
 * such transition.
 */
struct AtomicShellRadiativeEmission {
    std::uint64_t vacancy = 0;   ///< Subshell identifier of the shell from which the vacancy is filled.
    double energy = 0;           ///< Energy of the emitted characteristic X-ray photon [keV].
    double probability = 0;      ///< Relative probability of this transition (unnormalised within the shell).
};

/**
 * @brief Tabulated data for a single atomic subshell used in photoelectric and fluorescence sampling.
 *
 * Stores the orbital parameters, per-subshell photoelectric cross-section table, and the
 * list of radiative (fluorescence) emission transitions that can fill a vacancy created
 * in this subshell. Used by `AtomicElement::shells` to model detailed photoelectric
 * interactions including characteristic X-ray emission.
 */
struct AtomicShell {
    /**
     * @brief Constructs an `AtomicShell` with the given subshell identifier.
     * @param shell  Subshell identifier (e.g. EADL subshell code). Default: 0.
     */
    AtomicShell(std::uint64_t shell = 0)
        : shell(shell)
    {
    }
    std::uint64_t shell = 0;                                        ///< Subshell identifier (EADL subshell code).
    double numberOfElectrons = 0;                                   ///< Number of electrons occupying this subshell.
    double bindingEnergy = 0;                                       ///< Electron binding energy [keV].
    double kineticEnergy = 0;                                       ///< Mean kinetic energy of an electron in this subshell [keV].
    double HartreeFockOrbital_0 = 0;                                ///< Hartree–Fock orbital momentum parameter p₀ [a.u.].
    double numberOfPhotonsPerInitVacancy = 0;                       ///< Average number of fluorescence photons emitted per initial vacancy (fluorescence yield).
    double energyOfPhotonsPerInitVacancy = 0;                       ///< Average total energy carried away by fluorescence photons per initial vacancy [keV].
    std::vector<std::pair<double, double>> photoel;                 ///< Per-subshell photoelectric cross-section table: {energy [keV], μ/ρ [cm²/g]}.
    std::vector<AtomicShellRadiativeEmission> radiativeEmissions;   ///< List of radiative transitions that can fill a vacancy in this subshell.
};

}
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

#include "xraymc/material/atomicelement.hpp"
#include "xraymc/material/atomserializer.hpp"

#include <concepts>
#include <execution>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace xraymc {

/**
 * @brief Singleton registry of photon-interaction data for all chemical elements.
 *
 * Loads the binary physics-data file (`physicslists.bin`) once at program startup
 * and provides static accessors to `AtomicElement` records keyed by atomic number Z.
 * The data file is located by first checking the current working directory, then
 * falling back to the compile-time path `XRAYMCLIB_PHYSICSLISTSPATH`.
 *
 * All public methods are static; the singleton is managed internally via `Instance()`.
 * Copy construction and assignment are deleted to enforce the singleton pattern.
 */
class AtomHandler {
public:
    /**
     * @brief Returns the `AtomicElement` for the given atomic number.
     *
     * If @p Z is not found in the loaded data (e.g. because the element is outside
     * the supported range or the data file failed to load), a reference to an empty
     * dummy element is returned.
     *
     * @param Z  Atomic number (integral type).
     * @return Const reference to the corresponding `AtomicElement`, or a dummy if absent.
     */
    static const AtomicElement& Atom(std::integral auto Z)
    {
        auto& instance = Instance();
        if (instance.m_elements.contains(Z)) {
            return instance.m_elements.at(Z);
        }
        return instance.m_dummyElement;
    }

    /**
     * @brief Returns true if data for atomic number @p Z has been loaded.
     * @param Z  Atomic number (integral type).
     * @return True if the element exists in the registry; false otherwise.
     */
    static bool atomExists(std::integral auto Z)
    {
        auto& instance = Instance();
        return instance.m_elements.contains(Z);
    }

    /**
     * @brief Returns the full map of all loaded elements, keyed by atomic number.
     * @return Const reference to the internal `{Z → AtomicElement}` map.
     */
    static const std::map<std::uint8_t, AtomicElement>& allAtoms()
    {
        const auto& instance = Instance();
        return instance.m_elements;
    }

    /**
     * @brief Converts an atomic number to its chemical symbol string.
     *
     * Supports Z in [1, 100] (H through Fm). Returns an empty string for
     * values outside this range.
     *
     * @param Z  Atomic number (integral type).
     * @return Chemical symbol (e.g. "H", "Fe", "Au"), or "" if Z is out of range.
     */
    static std::string toSymbol(std::integral auto Z)
    {
        std::string res;
        if (0 < Z && Z < 101) {
            const std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };
            res = S[Z - 1];
        }
        return res;
    }

    /// @brief Copy construction is deleted (singleton).
    AtomHandler(const AtomHandler&) = delete;
    /// @brief Copy assignment is deleted (singleton).
    void operator=(const AtomHandler&) = delete;

protected:
    /**
     * @brief Returns the singleton `AtomHandler` instance, constructing it on first call.
     * @return Const reference to the singleton.
     */
    static const AtomHandler& Instance()
    {
        static AtomHandler instance;
        return instance;
    }

    /**
     * @brief Constructs the singleton by locating and loading the physics-data file.
     *
     * Searches for `physicslists.bin` first in the current working directory, then
     * at the compile-time path `XRAYMCLIB_PHYSICSLISTSPATH`. If found, the file is
     * read into memory and deserialized via `AtomSerializer::deserializeAtoms`.
     * If the file cannot be found or read, `m_elements` remains empty and all
     * `Atom()` queries return the dummy element.
     */
    AtomHandler()
    {
        // reading data
        std::string datapath;
        if (std::filesystem::exists("physicslists.bin")) {
            datapath = "physicslists.bin";
        } else {
            const std::string EPICSdataBuildPath { XRAYMCLIB_PHYSICSLISTSPATH };
            if (std::filesystem::exists(EPICSdataBuildPath)) {
                datapath = EPICSdataBuildPath;
            }
        }

        // reading buffer
        std::ifstream buffer_file(datapath, std::ios::binary);
        if (buffer_file.good()) {
            std::vector<char> data(std::istreambuf_iterator<char>(buffer_file), {});
            m_elements = AtomSerializer::deserializeAtoms<std::uint8_t>(data);
        }
    }

private:
    AtomicElement m_dummyElement;                        ///< Returned for unknown atomic numbers.
    std::map<std::uint8_t, AtomicElement> m_elements;   ///< Loaded element data keyed by atomic number Z.
};

}
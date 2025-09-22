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

class AtomHandler {
public:
    static const AtomicElement& Atom(std::uint64_t Z)
    {
        auto& instance = Instance();
        if (instance.m_elements.contains(Z)) {
            return instance.m_elements.at(Z);
        }
        return instance.m_dummyElement;
    }

    static bool atomExists(std::uint64_t Z)
    {
        auto& instance = Instance();
        return instance.m_elements.contains(Z);
    }

    static const std::map<std::uint64_t, AtomicElement>& allAtoms()
    {
        const auto& instance = Instance();
        return instance.m_elements;
    }

    static std::string toSymbol(std::integral auto Z)
    {
        std::string res;
        if (0 < Z && Z < 101) {
            const std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };
            res = S[Z - 1];
        }
        return res;
    }

    AtomHandler(const AtomHandler&) = delete;
    void operator=(const AtomHandler&) = delete;

protected:
    static const AtomHandler& Instance()
    {
        static AtomHandler instance;
        return instance;
    }
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

            m_elements = AtomSerializer::deserializeAtoms(data);
        }
    }

private:
    AtomicElement m_dummyElement;
    std::map<std::uint64_t, AtomicElement> m_elements;
};

}
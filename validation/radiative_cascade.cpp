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

Copyright 2024 Erlend Andersen
*/

/**
 * @file radiative_cascade.cpp
 * @brief Validation tool that tabulates the radiative (fluorescence) cascade of an element.
 *
 * For every subshell of the element, an initial vacancy is assumed and all chains of
 * radiative transitions that move the vacancy outwards are followed recursively.
 * Each emitted characteristic X-ray line is written to `cascade.txt` in the current
 * working directory together with its cumulative probability per initial vacancy.
 *
 * Only radiative branches are followed; vacancies created by non-radiative (Auger and
 * Coster-Kronig) transitions are not included.
 */

#include "xraymc/material/atomhandler.hpp"
#include <fstream>
#include <iostream>

/**
 * @brief One characteristic X-ray line in a radiative cascade.
 */
struct Line {
    std::uint64_t from;   ///< Subshell holding the vacancy before the transition.
    std::uint64_t to;     ///< Subshell the filling electron comes from, i.e. where the vacancy moves to.
    std::uint64_t parent; ///< Subshell of the initial vacancy that started the cascade.
    double probability;   ///< Cumulative probability of this line per initial vacancy in @ref parent.
    double energy;        ///< Energy of the emitted photon [keV].
};

/**
 * @brief Recursively collects the radiative transitions that fill a vacancy.
 *
 * Appends one Line for each radiative transition from subshell @p from. Each line
 * is weighted by the probability of reaching @p from. The function then recurses
 * into the subshell that receives the new vacancy.
 *
 * @param atom   Atomic data for the element.
 * @param parent Subshell of the initial vacancy, recorded in every line.
 * @param from   Subshell currently holding the vacancy.
 * @param prob   Probability of a vacancy in @p from per initial vacancy in @p parent.
 * @param lines  Output vector that receives the emitted lines.
 */
void write_cascade_worker(const xraymc::AtomicElement& atom, std::uint64_t parent, std::uint64_t from, double prob, std::vector<Line>& lines)
{
    const auto& shell = atom.shells.at(from);
    for (const auto& em : shell.radiativeEmissions) {
        Line line;
        line.energy = em.energy;
        line.probability = em.probability * prob;
        line.from = from;
        line.to = em.vacancy;
        line.parent = parent;
        lines.push_back(line);
        write_cascade_worker(atom, parent, em.vacancy, prob * em.probability, lines);
    }
}

/**
 * @brief Writes the radiative cascades of all subshells of an element to `cascade.txt`.
 *
 * Starts one cascade per subshell and writes every resulting line as CSV with the
 * columns `Z,parent,from,to,probability,energy`. An existing file is overwritten.
 *
 * @param Z Atomic number of the element.
 */
void write_cascade(std::uint8_t Z)
{
    const auto a = xraymc::AtomHandler::Atom(Z);
    std::vector<Line> lines;
    for (const auto& [shell, shell_data] : a.shells) {
        write_cascade_worker(a, shell, shell, 1.0, lines);
    }

    std::ofstream f("cascade.txt");
    f << "Z,parent,from,to,probability,energy\n";
    for (const auto& l : lines) {
        f << static_cast<std::uint64_t>(Z) << ',' << l.parent << ',' << l.from << ',' << l.to << ',';
        f << l.probability << ',' << l.energy << '\n';
    }
}

/**
 * @brief Writes the radiative cascade for lead (Z = 82).
 */
int main()
{
    write_cascade(82);
    return EXIT_SUCCESS;
}

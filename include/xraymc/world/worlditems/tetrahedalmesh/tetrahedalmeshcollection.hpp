
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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include <cstdint>
#include <string>

namespace xraymc {
/**
 * @brief Aggregate dose and geometry data for a named collection of tetrahedra.
 *
 * A TetrahedalMesh partitions its tetrahedra into named collections, each representing
 * a distinct tissue or material region. This struct holds the per-collection summary
 * produced after dose conversion: the material density, total volume, mean dose,
 * dose variance, interaction event count, and the collection name.
 */
struct TetrahedalMeshCollection {
    double density = 0;          ///< Material mass density of the collection in g/cm³.
    double volume = 0;           ///< Total volume of all tetrahedra in the collection in cm³.
    double dose = 0;             ///< Mean absorbed dose in the collection in Gy (or simulation units).
    double doseVariance = 0;     ///< Variance of the absorbed dose estimate.
    std::uint64_t numberOfEvents = 0; ///< Total number of energy-deposition events scored in this collection.
    std::string name;            ///< Human-readable name identifying the collection.
};
}
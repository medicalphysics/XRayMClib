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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include <algorithm>
#include <execution>
#include <limits>
#include <map>
#include <ranges>
#include <vector>

namespace dxmc {

struct TetrahedalMeshData {
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::array<std::uint32_t, 4>> elements;
    std::vector<std::uint32_t> collectionIndices; // same size as elements

    // specify collection/organ properties
    std::vector<double> collectionDensities;
    std::vector<std::map<std::size_t, double>> collectionMaterialComposition;
    std::vector<std::string> collectionNames;

    auto maxCollectionNumber() const
    {
        auto iter = std::max_element(std::execution::par_unseq, collectionIndices.cbegin(), collectionIndices.cend());
        return *iter;
    }
    bool valid() const
    {
        // test if element indices can index nodes
        const auto n_nodes = nodes.size();
        auto pos = std::find_if(std::execution::par_unseq, elements.cbegin(), elements.cend(), [n_nodes](const auto& e) {
            bool v = true;
            for (const auto& p : e)
                v = v && p < n_nodes;
            return !v;
        });

        bool val = elements.cend() == pos;

        // test if collection indices is complete
        const auto max_idx = maxCollectionNumber() + 1;
        val = val && collectionDensities.size() == max_idx;
        // val = val && collectionMaterialComposition.size() == max_idx;
        val = val && collectionNames.size() == max_idx;
        return val;
    }
};
}
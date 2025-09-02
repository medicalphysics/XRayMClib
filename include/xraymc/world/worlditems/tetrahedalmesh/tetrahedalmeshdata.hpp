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

#include <algorithm>
#include <array>
#include <execution>
#include <limits>
#include <map>
#include <ranges>
#include <string>
#include <vector>

namespace xraymc {

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
        val = val && collectionMaterialComposition.size() == max_idx;
        val = val && collectionNames.size() == max_idx;
        return val;
    }

    void makeGenericCollectionNames()
    {
        const auto N = maxCollectionNumber() + 1;
        if (collectionNames.size() < N)
            collectionNames.resize(N);
        for (std::size_t i = 0; i < N; ++i) {
            if (collectionNames[i].empty()) {
                collectionNames[i] = "Generic name " + std::to_string(i + 1);
            }
        }
    }

    void collectionNameMustContainFilter(const std::string& filter)
    {
        std::vector<std::uint32_t> indicesToKeep;
        for (std::uint32_t i = 0; i < collectionNames.size(); ++i) {
            const auto& collname = collectionNames[i];
            if (collname.contains(filter)) {
                indicesToKeep.push_back(i);
            }
        }

        if (indicesToKeep.size() == 0) {
            return;
        }
        std::sort(indicesToKeep.begin(), indicesToKeep.end());

        struct FilterData {
            std::array<std::uint32_t, 4> element;
            std::uint32_t collectionIndex;
        };

        std::vector<FilterData> filtered(elements.size());
        std::transform(std::execution::par_unseq, elements.cbegin(), elements.cend(), collectionIndices.cbegin(), filtered.begin(), [](const auto& el, const auto i) {
            return FilterData { .element = el, .collectionIndex = i };
        });

        auto last = std::remove_if(std::execution::par_unseq, filtered.begin(), filtered.end(), [&indicesToKeep](const auto& f) {
            return !std::binary_search(indicesToKeep.cbegin(), indicesToKeep.cend(), f.collectionIndex);
        });
        filtered.erase(last, filtered.end());
        elements.clear();
        collectionIndices.clear();
        for (const auto& f : filtered) {
            elements.push_back(f.element);
            collectionIndices.push_back(f.collectionIndex);
        }

        // remove unneeded vertices
        std::vector<std::uint32_t> indices_to_keep;
        indices_to_keep.reserve(elements.size() * 4);
        for (const auto& el : elements)
            for (auto i : el)
                indices_to_keep.push_back(i);
        std::sort(indices_to_keep.begin(), indices_to_keep.end());
        indices_to_keep.erase(std::unique(indices_to_keep.begin(), indices_to_keep.end()), indices_to_keep.end());
        std::map<std::uint32_t, std::uint32_t> indices_map;
        for (std::uint32_t i = 0; i < indices_to_keep.size(); ++i)
            indices_map[indices_to_keep[i]] = i;
        // renaming indices
        for (auto& el : elements)
            for (std::uint32_t i = 0; i < 4; ++i)
                el[i] = indices_map[el[i]];
        // removing nodes
        std::vector<std::array<double, 3>> nodes_to_keep;
        nodes_to_keep.reserve(indices_map.size());
        for (std::uint32_t i = 0; i < nodes.size(); ++i)
            if (std::binary_search(indices_to_keep.cbegin(), indices_to_keep.cend(), i))
                nodes_to_keep.push_back(nodes[i]);
        nodes = nodes_to_keep;
    }
};
}
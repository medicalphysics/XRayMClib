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

#include "xraymc/world/basicshapes/tetrahedron.hpp"

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

    bool testTetrahedronNormals() const
    {
        return std::all_of(std::execution::par_unseq, elements.cbegin(), elements.cend(), [this](const auto& el) {
            return basicshape::tetrahedron::tetrahedronNormalsOutsideTest(
                nodes[el[0]],
                nodes[el[1]],
                nodes[el[2]],
                nodes[el[3]]);
        });
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

    void collectionNameMustContainFilter(const std::string& filter, bool caseSensitive = true)
    {
        std::vector<std::uint32_t> indicesToKeep;
        if (caseSensitive) {
            for (std::uint32_t i = 0; i < collectionNames.size(); ++i) {
                const auto& collname = collectionNames[i];
                if (collname.contains(filter)) {
                    indicesToKeep.push_back(i);
                }
            }
        } else {
            auto caseFilter = filter;
            for (auto& c : caseFilter)
                c = std::tolower(c);
            for (std::uint32_t i = 0; i < collectionNames.size(); ++i) {
                auto collname = collectionNames[i];
                std::transform(collname.cbegin(), collname.cend(), collname.begin(), [](auto& c) { return std::tolower(c); });
                if (collname.contains(caseFilter)) {
                    indicesToKeep.push_back(i);
                }
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

        // pruning nodes, elements and collection indices
        // bilding swap table
        std::map<std::uint32_t, std::uint32_t> element_map;
        std::map<std::uint32_t, std::uint32_t> collection_map;
        {
            std::size_t i = 0;
            std::size_t c = 0;
            for (const auto& f : filtered) {
                if (!collection_map.contains(f.collectionIndex))
                    collection_map[f.collectionIndex] = c++;

                for (const auto& e : f.element)
                    if (!element_map.contains(e))
                        element_map[e] = i++;
            }
        }

        // swapping:
        for (const auto [key, value] : element_map)
            std::swap(nodes[key], nodes[value]);
        nodes.resize(element_map.size());
        for (const auto [key, value] : collection_map) {
            std::swap(collectionDensities[key], collectionDensities[value]);
            std::swap(collectionMaterialComposition[key], collectionMaterialComposition[value]);
            std::swap(collectionNames[key], collectionNames[value]);
        }
        collectionDensities.resize(collection_map.size());
        collectionMaterialComposition.resize(collection_map.size());
        collectionNames.resize(collection_map.size());

        elements.clear();
        collectionIndices.clear();
        for (const auto& f : filtered) {
            std::array<std::uint32_t, 4> m;
            for (std::size_t i = 0; i < 4; ++i) {
                m[i] = element_map.at(f.element[i]);
            }
            elements.push_back(m);
            collectionIndices.push_back(collection_map.at(f.collectionIndex));
        }
    }
};
}
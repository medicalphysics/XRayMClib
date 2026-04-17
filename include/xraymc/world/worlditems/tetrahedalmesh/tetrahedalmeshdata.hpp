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

/**
 * @brief Raw geometry and per-collection material data for a tetrahedral mesh.
 *
 * This plain data struct holds all the inputs needed to construct a TetrahedalMesh:
 * the vertex node positions, the four-vertex element connectivity, a per-element
 * collection index, and three parallel arrays that describe each named collection
 * (density, material composition by elemental weight fractions, and display name).
 * Helper methods validate consistency, fill missing names, and filter the mesh down
 * to a named subset of collections.
 */
struct TetrahedalMeshData {
    std::vector<std::array<double, 3>> nodes; ///< Vertex positions in cm.
    std::vector<std::array<std::uint32_t, 4>> elements; ///< Tetrahedron vertex index quads.
    std::vector<std::uint32_t> collectionIndices; ///< Per-element collection index; same size as elements.

    std::vector<double> collectionDensities; ///< Mass density in g/cm³ for each collection.
    std::vector<std::map<std::uint8_t, double>> collectionMaterialComposition; ///< Elemental weight fractions for each collection.
    std::vector<std::string> collectionNames; ///< Human-readable name for each collection.

    /**
     * @brief Returns the largest collection index present in collectionIndices.
     * @return Maximum collection index value.
     */
    auto maxCollectionNumber() const
    {
        auto iter = std::max_element(std::execution::par_unseq, collectionIndices.cbegin(), collectionIndices.cend());
        return *iter;
    }

    /**
     * @brief Checks whether the mesh data is internally consistent.
     *
     * Verifies that all element vertex indices are within the bounds of the nodes
     * array, and that collectionDensities, collectionMaterialComposition, and
     * collectionNames each have exactly maxCollectionNumber() + 1 entries.
     * @return true if the data is valid and ready to use; false otherwise.
     */
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

    /**
     * @brief Replaces the elemental weight-fraction composition of a collection.
     * @param collectionIndex Index of the collection to modify.
     * @param newComposition  New elemental weight fractions keyed by atomic number.
     * @throws std::out_of_range if @p collectionIndex is out of bounds.
     */
    void changeMaterialComposition(std::uint8_t collectionIndex, const std::map<std::uint8_t, double>& newComposition)
    {
        auto& comp = collectionMaterialComposition.at(collectionIndex);
        comp.clear();
        for (const auto& [key, value] : newComposition)
            comp[key] = value;
    }

    /**
     * @brief Tests whether all tetrahedra have outward-pointing face normals.
     *
     * Uses the geometric normal-orientation test from basicshape::tetrahedron for
     * every element. A mesh that passes this test can be used for inside/outside
     * determination without further preprocessing.
     * @return true if every tetrahedron has consistently outward normals; false otherwise.
     */
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

    /**
     * @brief Fills any empty collection names with the placeholder "Generic name N".
     *
     * Resizes collectionNames to maxCollectionNumber() + 1 if needed, then assigns
     * a numbered placeholder to every entry that is currently empty.
     */
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

    /**
     * @brief Filters the mesh to retain only collections whose name contains @p filter.
     *
     * Elements belonging to collections that do not match are removed. The surviving
     * collections are remapped to a contiguous index range starting at 0, and the node
     * array is pruned to remove vertices no longer referenced by any remaining element.
     * If no collection matches, the mesh is left unchanged.
     * @param filter        Substring that a collection name must contain to be kept.
     * @param caseSensitive When false, both the names and the filter are lowercased
     *                      before comparison.
     */
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

        {
            std::uint32_t teller = 0;
            for (std::uint32_t i = 0; i < collectionIndices.size(); ++i) {
                if (std::binary_search(indicesToKeep.cbegin(), indicesToKeep.cend(), collectionIndices[i])) {
                    std::swap(collectionIndices[i], collectionIndices[teller]);
                    std::swap(elements[i], elements[teller]);
                    teller++;
                }
            }
            collectionIndices.resize(teller);
            elements.resize(teller);
        }

        // swapping:
        std::map<std::uint32_t, std::uint32_t> collection_map;
        for (std::uint32_t i = 0; i < indicesToKeep.size(); ++i)
            collection_map[indicesToKeep[i]] = i;
        for (const auto [key, value] : collection_map) {
            std::swap(collectionDensities[key], collectionDensities[value]);
            std::swap(collectionMaterialComposition[key], collectionMaterialComposition[value]);
            std::swap(collectionNames[key], collectionNames[value]);
        }
        collectionDensities.resize(collection_map.size());
        collectionMaterialComposition.resize(collection_map.size());
        collectionNames.resize(collection_map.size());
        std::transform(std::execution::par_unseq, collectionIndices.cbegin(), collectionIndices.cend(), collectionIndices.begin(), [&](auto c) {
            return collection_map.at(c);
        });

        // pruning nodes
        std::vector<std::uint32_t> indices;
        indices.reserve(elements.size() * 4);
        for (const auto& e : elements)
            for (const auto& n : e)
                indices.push_back(n);

        std::sort(std::execution::par_unseq, indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
        for (std::uint32_t i = 0; i < indices.size(); ++i)
            std::swap(nodes[i], nodes[indices[i]]);
        nodes.resize(indices.size());

        std::map<std::uint32_t, std::uint32_t> indices_map;
        for (std::uint32_t i = 0; i < indices.size(); ++i)
            indices_map[indices[i]] = i;

        std::transform(std::execution::par_unseq, elements.cbegin(), elements.cend(), elements.begin(), [&](auto el) {
            for (auto& n : el)
                n = indices_map.at(n);
            return el;
        });
    }
};
}
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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/particle.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/tetrahedron.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh2/tetrahedalmeshdata.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>
#include <vector>

namespace dxmc {

class TetrahedalMeshGrid2 {

public:
    TetrahedalMeshGrid2() { }
    TetrahedalMeshGrid2(const TetrahedalMeshData& data)
    {
        setData(data);
    }
    void setData(const TetrahedalMeshData& data)
    {
        m_nodes = data.nodes;
        m_elements = data.elements;
        calculateAABB();
        updateGrid();
    }

protected:
    void calculateAABB()
    {
        m_aabb = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (const auto& p : m_nodes) {
            for (std::size_t i = 0; i < 3; ++i) {
                m_aabb[i] = std::min(m_aabb[i], p[i]);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], p[i]);
            }
        }
        // padding AABB
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] -= GEOMETRIC_ERROR();
            m_aabb[i + 3] += GEOMETRIC_ERROR();
        }

        // update grid resolution
        for (std::size_t i = 0; i < 3; ++i)
            m_gridSpacing[i] = (m_aabb[i + 3] - m_aabb[i]) / m_gridDimensions[i];
    }

    std::array<std::uint32_t, 3> gridIndex(const std::array<double, 3>& p) const
    {
        std::array<std::uint32_t, 3> d = {
            static_cast<std::uint32_t>((p[0] - m_aabb[0]) / m_gridSpacing[0]),
            static_cast<std::uint32_t>((p[1] - m_aabb[1]) / m_gridSpacing[1]),
            static_cast<std::uint32_t>((p[2] - m_aabb[2]) / m_gridSpacing[2])
        };
        return d;
    }
    std::uint32_t gridIndex(const std::array<std::uint32_t, 3>& ind) const
    {
        return ind[0] + m_gridDimensions[0] * (ind[1] + m_gridDimensions[1] * ind[2]);
    }

    void updateGrid()
    {
        const auto N = std::reduce(m_gridDimensions.cbegin(), m_gridDimensions.cend(), std::uint32_t { 1 }, std::multiplies {});

        // finding aabbs
        std::vector<std::array<double, 6>> aabbs(m_elements.size());
        std::transform(std::execution::par_unseq, m_elements.cbegin(), m_elements.cend(), aabbs.begin(), [&](const auto& element) {
            const auto& v0 = m_nodes[element[0]];
            const auto& v1 = m_nodes[element[1]];
            const auto& v2 = m_nodes[element[2]];
            const auto& v3 = m_nodes[element[3]];
            return basicshape::tetrahedron::AABB(v0, v1, v2, v3);
        });

        // assigning voxels
        std::vector<std::vector<std::uint32_t>> indices(N);

        for (std::uint32_t i = 0; i < m_elements.size(); ++i) {
            const auto& [l, r] = vectormath::splice(aabbs[i]);
            const auto find = gridIndex(l);
            const auto rind = gridIndex(r);
            for (std::uint32_t z = find[2]; z <= rind[2]; ++z)
                for (std::uint32_t y = find[1]; y <= rind[1]; ++y)
                    for (std::uint32_t x = find[0]; x <= rind[0]; ++x) {
                        std::array<double, 6> vox_aabb = {
                            x * m_gridSpacing[0] + m_aabb[0],
                            x * m_gridSpacing[1] + m_aabb[1],
                            x * m_gridSpacing[2] + m_aabb[2],
                            (x + 1) * m_gridSpacing[0] + m_aabb[0],
                            (x + 1) * m_gridSpacing[1] + m_aabb[1],
                            (x + 1) * m_gridSpacing[2] + m_aabb[2],
                        };
                        const auto& v0 = m_nodes[m_elements[i][0]];
                        const auto& v1 = m_nodes[m_elements[i][1]];
                        const auto& v2 = m_nodes[m_elements[i][2]];
                        const auto& v3 = m_nodes[m_elements[i][3]];
                        if (basicshape::tetrahedron::insideAABB(v0, v1, v2, v3, vox_aabb)) {
                            const std::array ind = { x, y, z };
                            const auto ind_flat = gridIndex(ind);
                            indices[ind_flat].push_back(i);
                        }
                    }
        }
        // sorting and removing duplicates if any
        std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [](auto& i_vec) {
            std::sort(i_vec.begin(), i_vec.end());
            auto last = std::unique(i_vec.begin(), i_vec.end());
            i_vec.erase(last, i_vec.end());
        });

        // assigning grid Indexes
        m_gridIndices.resize(N + 1, 0);

        auto N_indices = std::transform_reduce(std::execution::par_unseq, indices.cbegin(), indices.cend(), std::size_t { 0 }, std::plus {}, [](const auto& i_vec) { return i_vec.size(); });
        m_gridElements.reserve(N_indices);

        std::uint32_t teller = 0;
        for (std::size_t i = 0; i < indices.size(); ++i) {
            m_gridIndices[i] = teller;
            for (auto ind : indices[i]) {
                m_gridElements.push_back(ind);
                teller++;
            }
        }
        m_gridIndices[N] = teller;
    }

private:
    std::vector<std::array<double, 3>> m_nodes;
    std::vector<std::array<std::uint32_t, 4>> m_elements;

    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<double, 3> m_gridSpacing = { 1, 1, 1 };
    std::array<std::uint32_t, 3> m_gridDimensions = { 2, 2, 2 };
    std::vector<std::uint32_t> m_gridIndices; // same size as grid;
    std::vector<std::uint32_t> m_gridElements; // most likely larger than elements
};
}
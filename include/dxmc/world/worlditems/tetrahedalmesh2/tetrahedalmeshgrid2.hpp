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
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>
#include <vector>

namespace dxmc {

class TetrahedalMeshGrid2 {

public:
    TetrahedalMeshGrid2(const std::array<std::uint32_t, 3>& dimensions = { 8, 8, 8 })
        : m_gridDimensions(dimensions)
    {
    }
    TetrahedalMeshGrid2(const TetrahedalMeshData& data, const std::array<std::uint32_t, 3>& dimensions)
        : m_gridDimensions(dimensions)
    {
        setData(data);
    }
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
    void setData(const TetrahedalMeshData& data, const std::array<std::uint32_t, 3>& dimensions)
    {
        m_gridDimensions = dimensions;
        m_nodes = data.nodes;
        m_elements = data.elements;
        calculateAABB();
        updateGrid();
    }
    void setDimensions(const std::array<std::uint32_t, 3>& dimensions)
    {
        m_gridDimensions = dimensions;
        updateGrid();
    }

    const auto& dimensions() const
    {
        return m_gridDimensions;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void translate(const std::array<double, 3>& vec)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += vec[i];
            m_aabb[i + 3] += vec[i];
        }
        std::for_each(std::execution::par_unseq, m_nodes.begin(), m_nodes.end(), [&vec](auto& n) {
            n = vectormath::add(n, vec);
        });
    }

    void rotate(const std::array<double, 3>& axis, double angle)
    {
        std::transform(std::execution::par_unseq, m_nodes.cbegin(), m_nodes.cend(), m_nodes.begin(), [angle, &axis](const auto& v) {
            return vectormath::rotate(v, axis, angle);
        });
        calculateAABB();
        updateGrid();
    }

    std::array<std::uint32_t, 3> gridIndex(const std::array<double, 3>& p) const
    {
        std::array<std::uint32_t, 3> d = {
            static_cast<std::uint32_t>(std::clamp((p[0] - m_aabb[0]) / m_gridSpacing[0], 0.0, static_cast<double>(m_gridDimensions[0] - 1))),
            static_cast<std::uint32_t>(std::clamp((p[1] - m_aabb[1]) / m_gridSpacing[1], 0.0, static_cast<double>(m_gridDimensions[1] - 1))),
            static_cast<std::uint32_t>(std::clamp((p[2] - m_aabb[2]) / m_gridSpacing[2], 0.0, static_cast<double>(m_gridDimensions[2] - 1)))
        };
        return d;
    }
    std::uint32_t gridIndex(const std::array<std::uint32_t, 3>& ind) const
    {
        return ind[0] + m_gridDimensions[0] * (ind[1] + m_gridDimensions[1] * ind[2]);
    }
    std::uint32_t numberOfTetrahedrons(std::uint32_t index) const
    {
        auto start = m_gridIndices[index];
        auto stop = m_gridIndices[index + 1];
        return stop - start;
    }

    const std::vector<std::array<double, 3>>& nodes() const
    {
        return m_nodes;
    }
    const std::vector<std::array<std::uint32_t, 4>>& elements() const
    {
        return m_elements;
    }

    std::vector<double> volumes() const
    {
        std::vector<double> v(m_elements.size());
        std::transform(std::execution::par_unseq, m_elements.cbegin(), m_elements.cend(), v.begin(), [&](const auto& e) {
            const auto& v0 = m_nodes[e[0]];
            const auto& v1 = m_nodes[e[1]];
            const auto& v2 = m_nodes[e[2]];
            const auto& v3 = m_nodes[e[3]];
            return basicshape::tetrahedron::volume(v0, v1, v2, v3);
        });
        return v;
    }

    std::optional<std::uint32_t> pointInside(const std::array<double, 3>& pos) const
    {
        const auto gridIdx = gridIndex(gridIndex(pos));
        const auto start = m_gridIndices[gridIdx];
        const auto stop = m_gridIndices[gridIdx + 1];

        for (auto i = start; i < stop; ++i) {
            auto elIdx = m_gridElements[i];
            const auto& v0 = m_nodes[m_elements[elIdx][0]];
            const auto& v1 = m_nodes[m_elements[elIdx][1]];
            const auto& v2 = m_nodes[m_elements[elIdx][2]];
            const auto& v3 = m_nodes[m_elements[elIdx][3]];
            if (basicshape::tetrahedron::pointInside(v0, v1, v2, v3, pos)) {
                return elIdx;
            }
        }
        return std::nullopt;
    }

    KDTreeIntersectionResult<const std::uint32_t> intersect(const ParticleType auto& particle) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<false>(particle, m_aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<const std::uint32_t> {};
    }

protected:
    static inline int argmin3(const std::array<double, 3>& a)
    {
        return a[0] < a[2] ? a[0] < a[1] ? 0 : 1 : a[2] < a[1] ? 2
                                                               : 1;
    }

    KDTreeIntersectionResult<const std::uint32_t> intersect(const ParticleType auto& p, const std::array<double, 2>& t) const
    {
        auto idx = gridIndex(vectormath::add(p.pos, vectormath::scale(p.dir, t[0])));
        const std::array<int, 3> step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };
        const std::array<double, 3> delta = {
            m_gridSpacing[0] / std::abs(p.dir[0]),
            m_gridSpacing[1] / std::abs(p.dir[1]),
            m_gridSpacing[2] / std::abs(p.dir[2])
        };

        std::array<double, 3> tmax;
        for (int i = 0; i < 3; ++i) {
            tmax[i] = step[i] > 0 ? (m_aabb[i] + (idx[i] + 1) * m_gridSpacing[i] - p.pos[i]) / p.dir[i] : (m_aabb[i] + idx[i] * m_gridSpacing[i] - p.pos[i]) / p.dir[i];
        };

        int dimension = argmin3(tmax);
        KDTreeIntersectionResult<const std::uint32_t> res;
        res.intersection = std::numeric_limits<double>::max();
        bool cont = true;
        while (cont) {
            // we have a valid voxel, check intersections
            const auto grid_ind = gridIndex(idx);
            for (auto flat_idx = m_gridIndices[grid_ind]; flat_idx < m_gridIndices[grid_ind + 1]; ++flat_idx) {
                const auto elementIdx = m_gridElements[flat_idx];
                const auto& v0 = m_nodes[m_elements[elementIdx][0]];
                const auto& v1 = m_nodes[m_elements[elementIdx][1]];
                const auto& v2 = m_nodes[m_elements[elementIdx][2]];
                const auto& v3 = m_nodes[m_elements[elementIdx][3]];

                const auto res_cand = basicshape::tetrahedron::intersect(v0, v1, v2, v3, p);
                if (res_cand.valid() && res_cand.intersection <= tmax[dimension] && res_cand.intersection < res.intersection) {
                    res.intersection = res_cand.intersection;
                    res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                    res.item = &m_gridElements[flat_idx];
                    if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                        return res;
                }
            }
            if (res.valid()) {
                cont = false;
            } else {
                idx[dimension] += step[dimension];
                if (idx[dimension] < m_gridDimensions[dimension]) { // we test for less than zero by overflow
                    tmax[dimension] += delta[dimension];
                    dimension = argmin3(tmax);
                } else {
                    cont = false;
                }
            }
        }
        return res;
    }

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
                            y * m_gridSpacing[1] + m_aabb[1],
                            z * m_gridSpacing[2] + m_aabb[2],
                            (x + 1) * m_gridSpacing[0] + m_aabb[0],
                            (y + 1) * m_gridSpacing[1] + m_aabb[1],
                            (z + 1) * m_gridSpacing[2] + m_aabb[2],
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
        m_gridIndices.clear();
        m_gridIndices.resize(N + 1, 0);

        auto N_indices = std::transform_reduce(std::execution::par_unseq, indices.cbegin(), indices.cend(), std::size_t { 0 }, std::plus {}, [](const auto& i_vec) { return i_vec.size(); });
        m_gridElements.clear();
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
    std::array<std::uint32_t, 3> m_gridDimensions = { 8, 8, 8 };
    std::vector<std::uint32_t> m_gridIndices; // same size as grid + 1;
    std::vector<std::uint32_t> m_gridElements; // most likely larger than elements
};
}
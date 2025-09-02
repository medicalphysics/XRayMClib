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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "xraymc/particle.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/kdtreeintersectionresult.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>
#include <vector>

namespace xraymc {

class TetrahedalMeshGrid {
public:
    TetrahedalMeshGrid() { }
    TetrahedalMeshGrid(const std::vector<Tetrahedron>& tets, int depth = 8)
    {
        setData(tets, depth);
    }
    TetrahedalMeshGrid(const std::vector<Tetrahedron>& tets, int dx, int dy, int dz)
    {
        setData(tets, dx, dy, dz);
    }
    TetrahedalMeshGrid(const std::vector<Tetrahedron>& tets, const std::array<int, 3> depth)
    {
        setData(tets, depth);
    }

    void setData(const std::vector<Tetrahedron>& tets, int depth = 8)
    {
        std::array<int, 3> d = { depth, depth, depth };
        setData(tets, d);
    }

    void setData(const std::vector<Tetrahedron>& tets, int dx, int dy, int dz)
    {
        std::array<int, 3> d = { dx, dy, dz };
        setData(tets, d);
    }

    void setData(const std::vector<Tetrahedron>& tets, const std::array<int, 3>& depth)
    {
        for (int i = 0; i < 3; ++i)
            m_N[i] = std::clamp(depth[i], 1, 1000);

        m_tets = tets;
        // sorting tets by diagonal
        std::sort(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [](const auto& lh, const auto& rh) {
            constexpr std::array<double, 3> n = { 1, 1, 1 };
            return vectormath::dot(lh.center(), n) < vectormath::dot(rh.center(), n);
        });

        calculateAABB();
        assignGrid();
    }

    std::vector<Tetrahedron>& tetrahedrons()
    {
        return m_tets;
    }

    const std::vector<Tetrahedron>& tetrahedrons() const
    {
        return m_tets;
    }

    void clearDoseScored()
    {
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
            tri.clearDoseScored();
        });
    }

    void clearEnergyScored()
    {
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
            tri.clearEnergyScored();
        });
    }

    void resample(const std::array<int, 3>& depth)
    {
        for (int i = 0; i < 3; ++i)
            m_N[i] = std::clamp(depth[i], 1, 1000);
        assignGrid();
    }

    std::size_t maxThetrahedronsVoxelCount() const
    {
        std::uint32_t c = 0;
        for (const auto v : m_gridIndices)
            c = std::max(c, v.end - v.begin);
        return c;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    template <std::regular_invocable<Tetrahedron&> F>
    void apply(const F& func, auto policy = std::execution::par_unseq)
    {
        std::for_each(policy, m_tets.begin(), m_tets.end(), func);
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (std::uint_fast32_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
            tri.translate(dist);
        });
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<const Tetrahedron> intersect(const ParticleType auto& particle) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<false>(particle, m_aabb);
        return inter ? intersect<COLLECTION>(particle, *inter) : KDTreeIntersectionResult<const Tetrahedron> {};
    }

    const Tetrahedron* pointInside(const std::array<double, 3>& pos) const
    {
        const auto idx_flat = m_gridIndices[getIndicesFlat<true>(pos)];
        for (auto gidx = idx_flat.begin; gidx < idx_flat.end; ++gidx) {
            const auto tet_idx = m_tetIndices[gidx];
            const auto& tet = m_tets[tet_idx];
            if (tet.pointInside(pos)) {
                return &tet;
            }
        }
        return nullptr;
    }

    Tetrahedron* pointInside(const std::array<double, 3>& pos)
    {
        const auto idx_flat = m_gridIndices[getIndicesFlat<true>(pos)];
        for (auto gidx = idx_flat.begin; gidx < idx_flat.end; ++gidx) {
            const auto tet_idx = m_tetIndices[gidx];
            auto& tet = m_tets[tet_idx];
            if (tet.pointInside(pos)) {
                return &tet;
            }
        }
        return nullptr;
    }

    template <bool CHECK_BOUNDS = false>
    std::array<int, 3> getIndices(const std::array<double, 3>& pos) const
    {
        if constexpr (CHECK_BOUNDS) {
            std::array<int, 3> res = {
                std::clamp(static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]), 0, m_N[0] - 1),
                std::clamp(static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]), 0, m_N[1] - 1),
                std::clamp(static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2]), 0, m_N[2] - 1)
            };
            return res;
        } else {
            std::array<int, 3> res = {
                static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]),
                static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]),
                static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2])
            };
            return res;
        }
    }

    template <bool CHECK_BOUNDS = false>
    int getIndicesFlat(const std::array<double, 3>& pos) const
    {
        return getIndicesFlat(getIndices<CHECK_BOUNDS>(pos));
    }

    int getIndicesFlat(const std::array<int, 3>& idx) const
    {
        return idx[0] + m_N[0] * (idx[1] + m_N[1] * idx[2]);
    }

    int getIndicesFlat(int x, int y, int z) const
    {
        return x + m_N[0] * (y + m_N[1] * z);
    }

protected:
    static inline int argmin3(const std::array<double, 3>& a)
    {
        return a[0] < a[2] ? a[0] < a[1] ? 0 : 1 : a[2] < a[1] ? 2
                                                               : 1;
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<const Tetrahedron> intersect(const ParticleType auto& p, const std::array<double, 2>& t) const
    {
        auto idx = getIndices<true>(vectormath::add(p.pos, vectormath::scale(p.dir, t[0])));
        const std::array<int, 3> step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };
        const std::array<double, 3> delta = {
            m_spacing[0] / std::abs(p.dir[0]),
            m_spacing[1] / std::abs(p.dir[1]),
            m_spacing[2] / std::abs(p.dir[2])
        };

        std::array<double, 3> tmax;
        for (int i = 0; i < 3; ++i) {
            tmax[i] = step[i] > 0 ? (m_aabb[i] + (idx[i] + 1) * m_spacing[i] - p.pos[i]) / p.dir[i] : (m_aabb[i] + idx[i] * m_spacing[i] - p.pos[i]) / p.dir[i];
        };

        int dimension = argmin3(tmax);
        KDTreeIntersectionResult<const Tetrahedron> res;
        res.intersection = std::numeric_limits<double>::max();
        bool cont = true;
        while (cont) {
            // we have a valid voxel, check intersections
            const auto grid_ind = m_gridIndices[getIndicesFlat(idx)];
            for (auto flat_idx = grid_ind.begin; flat_idx < grid_ind.end; ++flat_idx) {
                const auto& tet = m_tets[m_tetIndices[flat_idx]];
                if constexpr (COLLECTION == 65535) {
                    const auto res_cand = tet.intersect(p);
                    if (res_cand.valid() && res_cand.intersection <= tmax[dimension] && res_cand.intersection < res.intersection) {
                        res.intersection = res_cand.intersection;
                        res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                        res.item = &tet;
                        if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                            return res;
                    }
                } else {
                    if (tet.collection() == COLLECTION) {
                        const auto res_cand = tet.intersect(p);
                        if (res_cand.valid() && res_cand.intersection <= tmax[dimension] && res_cand.intersection < res.intersection) {
                            res.intersection = res_cand.intersection;
                            res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                            res.item = &tet;
                            if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                                return res;
                        }
                    }
                }
            }
            if (res.valid()) {
                cont = false;
            } else {
                idx[dimension] += step[dimension];
                if (0 <= idx[dimension] && idx[dimension] < m_N[dimension]) {
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
        for (const auto& tet : m_tets) {
            for (const auto& v : tet.vertices()) {
                for (std::uint_fast32_t i = 0; i < 3; ++i) {
                    m_aabb[i] = std::min(m_aabb[i], v[i]);
                    m_aabb[i + 3] = std::max(m_aabb[i + 3], v[i]);
                }
            }
        }
        for (std::uint_fast32_t i = 0; i < 3; ++i)
            m_spacing[i] = (m_aabb[i + 3] - m_aabb[i]) / m_N[i];
    }

    void assignGrid()
    {
        const auto [start, stop] = vectormath::splice(m_aabb);
        const auto size = std::reduce(m_N.cbegin(), m_N.cend(), 1, std::multiplies {});
        std::vector<std::vector<std::uint32_t>> grid(size);

        const std::array<int, 3> N = { m_N[0] - 1, m_N[1] - 1, m_N[2] - 1 };
        auto caster = [N](const std::array<double, 3>& v) -> std::array<int, 3> {
            std::array<int, 3> vi = {
                std::clamp(static_cast<int>(v[0]), 0, N[0]),
                std::clamp(static_cast<int>(v[1]), 0, N[1]),
                std::clamp(static_cast<int>(v[2]), 0, N[2])
            };
            return vi;
        };
        const std::array<double, 3> inv_spacing = { 1 / m_spacing[0], 1 / m_spacing[1], 1 / m_spacing[2] };
        for (std::uint32_t i = 0; i < m_tets.size(); ++i) {
            const auto tet_aabb = m_tets[i].AABB();
            const auto [tet_start, tet_stop] = vectormath::splice(tet_aabb);
            const auto start_ind = caster(vectormath::scale(vectormath::subtract(tet_start, start), inv_spacing));
            const auto stop_ind = caster(vectormath::scale(vectormath::subtract(tet_stop, start), inv_spacing));
            for (int z = start_ind[2]; z <= stop_ind[2]; ++z)
                for (int y = start_ind[1]; y <= stop_ind[1]; ++y)
                    for (int x = start_ind[0]; x <= stop_ind[0]; ++x) {
                        const int idx = getIndicesFlat(x, y, z);
                        grid[idx].push_back(i);
                    }
        }
        std::for_each(std::execution::par_unseq, grid.begin(), grid.end(), [](auto& v) { v.shrink_to_fit(); });

        // assigning indices
        std::uint32_t begin = 0;
        std::uint32_t end = 0;
        m_gridIndices.resize(grid.size());
        for (std::size_t gIdx = 0; gIdx < grid.size(); ++gIdx) {
            const auto& gVec = grid[gIdx];
            end = begin + static_cast<std::uint32_t>(gVec.size());
            m_gridIndices[gIdx] = { .begin = begin, .end = end };
            for (auto i : gVec)
                m_tetIndices.push_back(i);
            begin = end;
            grid[gIdx].clear();
            grid[gIdx].shrink_to_fit();
        }

        std::for_each(std::execution::par_unseq, m_gridIndices.cbegin(), m_gridIndices.cend(), [this](const auto& gIdx) {
            if (gIdx.begin < gIdx.end) {
                auto begin = m_tetIndices.begin() + gIdx.begin;
                auto end = m_tetIndices.begin() + gIdx.end;
                std::sort(begin, end);
            }
        });

        m_tetIndices.shrink_to_fit();
    }

private:
    struct GridIdx {
        std::uint32_t begin, end;
    };
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<double, 3> m_spacing = { 1, 1, 1 };
    std::array<int, 3> m_N = { 8, 8, 8 };
    std::vector<GridIdx> m_gridIndices; // same size as grid
    std::vector<std::uint32_t> m_tetIndices;
    std::vector<Tetrahedron> m_tets;
};
}
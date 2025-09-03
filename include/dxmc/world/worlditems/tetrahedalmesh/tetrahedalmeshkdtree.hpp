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

Copyright 2024 Erlend Andersen
*/

#pragma once
#include "dxmc/particle.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <optional>
#include <vector>

namespace dxmc {

// This is far worse in terms of speed than a structured grid structure

class TetrahedralMeshKDTree {
public:
    TetrahedralMeshKDTree() {};
    TetrahedralMeshKDTree(const std::vector<Tetrahedron>& tets, std::uint32_t max_depth = 8)
        : m_items(tets)
    {
        build(max_depth);
    }

    void setData(const std::vector<Tetrahedron>& tets, std::uint32_t max_depth = 8)
    {
        m_items = tets;
        build(max_depth);
    }

    std::vector<Tetrahedron>& tetrahedrons()
    {
        return m_items;
    }

    const std::vector<Tetrahedron>& tetrahedrons() const
    {
        return m_items;
    }

    void clearDoseScored()
    {
        std::for_each(std::execution::par_unseq, m_items.begin(), m_items.end(), [&](auto& tri) {
            tri.clearDoseScored();
        });
    }

    void clearEnergyScored()
    {
        std::for_each(std::execution::par_unseq, m_items.begin(), m_items.end(), [&](auto& tri) {
            tri.clearEnergyScored();
        });
    }

    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    std::size_t maxThetrahedronsVoxelCount() const
    {
        std::size_t max = 0;
        auto it = std::max_element(m_nodes.cbegin(), m_nodes.cend(), [](const auto& lh, const auto& rh) {
            if (lh.isLeaf() && rh.isLeaf())
                return lh.split_nelements.nelements < lh.split_nelements.nelements;
            else if (!lh.isLeaf() && !rh.isLeaf())
                return false;
            else
                return rh.isLeaf() ? true : false;
        });
        return static_cast<std::size_t>(it->split_nelements.nelements);
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void translate(const std::array<double, 3>& dist)
    {
        std::for_each(std::execution::par_unseq, m_items.begin(), m_items.end(), [&](auto& tri) {
            tri.translate(dist);
        });
        for (auto& node : m_nodes) {
            if (node.isLeaf()) {
                const auto dim = node.dim();
                node.split_nelements.split += dist[dim];
            }
        }
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<const Tetrahedron> intersect(const ParticleType auto& particle) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<true>(particle, m_aabb);
        return inter ? intersect<COLLECTION>(particle, *inter) : KDTreeIntersectionResult<const Tetrahedron> {};
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<const Tetrahedron> intersect(const ParticleType auto& particle, const std::array<double, 2>& tboxAABB) const
    {
        struct Stack {
            struct Element {
                Node node;
                std::array<double, 2> tbox;
            };
            std::array<Element, 32> items;
            std::uint32_t n_items = 0;

            Stack(const Node& node, const std::array<double, 2>& tbox)
            {
                addItem(node, tbox[0], tbox[1]);
            }
            void takeItem(Node& node, std::array<double, 2>& tbox)
            {
                node = items[n_items - 1].node;
                tbox = items[n_items - 1].tbox;
                --n_items;
            }
            void addItem(const Node& node, double tmin, double tmax)
            {
                items[n_items].node = node;
                items[n_items].tbox = { tmin, tmax };
                ++n_items;
            }
            bool isEmpty() const
            {
                return n_items == 0;
            }
            void clear()
            {
                n_items = 0;
            }
        };

        Stack stack(m_nodes[0], tboxAABB);
        KDTreeIntersectionResult<const Tetrahedron> res;
        res.intersection = std::numeric_limits<double>::max();
        Node node;
        std::array<double, 2> tbox;
        while (!stack.isEmpty()) {
            stack.takeItem(node, tbox);
            while (!node.isLeaf()) {
                const auto dim = node.dim();
                const auto split = node.split_nelements.split;

                const auto left_offset = node.offset();
                const auto right_offset = node.offset() + 1;

                // test for parallell beam
                if (std::abs(particle.dir[dim]) <= std::numeric_limits<double>::epsilon()) {
                    node = m_nodes[left_offset];
                    stack.addItem(m_nodes[right_offset], tbox[0], tbox[1]);
                } else {
                    const auto d = (split - particle.pos[dim]) / particle.dir[dim];
                    auto frontchild = particle.dir[dim] > 0 ? m_nodes[left_offset] : m_nodes[right_offset];
                    auto backchild = particle.dir[dim] > 0 ? m_nodes[right_offset] : m_nodes[left_offset];

                    if (d <= tbox[0]) {
                        // find back node
                        node = backchild;
                    } else if (d >= tbox[1]) { // find front node
                        node = frontchild;
                    } else {
                        stack.addItem(backchild, d, tbox[1]);
                        tbox[1] = d;
                        node = frontchild;
                    }
                }
            }
            // we have a leaf
            const auto startIdx = node.offset();
            const auto stopIdx = startIdx + node.split_nelements.nelements;
            for (auto idx = startIdx; idx < stopIdx; ++idx) {
                const auto& item = m_items[m_indices[idx]];
                if constexpr (COLLECTION == 65535) {
                    auto t_cand = item.intersect(particle);
                    if (t_cand.valid()) {
                        if (t_cand.intersection < res.intersection) {
                            if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                                res.intersection = t_cand.intersection;
                                res.item = &item;
                                res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                                if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                                    return res;
                            }
                        }
                    }
                } else {
                    if (item.collection() == COLLECTION) {
                        auto t_cand = item.intersect(particle);
                        if (t_cand.valid()) {
                            if (t_cand.intersection < res.intersection) {
                                if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                                    res.intersection = t_cand.intersection;
                                    res.item = &item;
                                    res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                                    if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                                        return res;
                                }
                            }
                        }
                    }
                }
            }
            if (res.valid())
                stack.clear();
        }
        return res;
    }

    Tetrahedron* pointInside(const std::array<double, 3>& pos)
    {
        Node node = m_nodes[0];
        while (!node.isLeaf()) {
            const auto dim = node.dim();
            const auto split = node.split_nelements.split;
            node = pos[dim] <= split ? m_nodes[node.offset()] : m_nodes[node.offset() + 1];
        }
        // we have a leaf
        // test all leaf tets
        const auto startIdx = node.offset();
        const auto stopIdx = startIdx + node.split_nelements.nelements;
        for (auto idx = startIdx; idx < stopIdx; ++idx) {
            auto& item = m_items[m_indices[idx]];
            if (item.pointInside(pos)) {
                return &item;
            }
        }
        return nullptr;
    }

    std::uint32_t pointInsideIndex(const std::array<double, 3>& pos) const
    {
        Node node = m_nodes[0];
        while (!node.isLeaf()) {
            const auto dim = node.dim();
            const auto split = node.split_nelements.split;
            node = pos[dim] <= split ? m_nodes[node.offset()] : m_nodes[node.offset() + 1];
        }
        // we have a leaf
        // test all leaf tets
        const auto startIdx = node.offset();
        const auto stopIdx = startIdx + node.split_nelements.nelements;
        for (auto idx = startIdx; idx < stopIdx; ++idx) {
            const auto& item = m_items[m_indices[idx]];
            if (item.pointInside(pos)) {
                return idx;
            }
        }
        return std::numeric_limits<std::uint32_t>::max();
    }

protected:
    void sortTetrahedrons()
    {
        std::sort(std::execution::par_unseq, m_items.begin(), m_items.end(), [](const auto& lh, const auto& rh) {
            constexpr std::array<double, 3> v = { 1, 1, 1 };
            const auto l = vectormath::dot(lh.center(), v);
            const auto r = vectormath::dot(rh.center(), v);
            return l < r;
        });
    }

    std::array<double, 6> calculateAABB(const std::vector<std::uint32_t>& indices) const
    {

        std::array<double, 6> aabb = { 0, 0, 0, 0, 0, 0 };

        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] = std::transform_reduce(
                std::execution::par_unseq, indices.cbegin(), indices.cend(), std::numeric_limits<double>::max(), [](const auto& lh, const auto& rh) { return std::min(lh, rh); }, [i, this](const auto& t) {                
                const auto& v = this->m_items[t].vertices();
                return std::min({ v[0][i], v[1][i], v[2][i], v[3][i] }); });
            aabb[i + 3] = std::transform_reduce(
                std::execution::par_unseq, indices.cbegin(), indices.cend(), std::numeric_limits<double>::lowest(), [](const auto& lh, const auto& rh) { return std::max(lh, rh); }, [i, this](const auto& t) {                
                const auto& v = this->m_items[t].vertices();
                return std::max({ v[0][i], v[1][i], v[2][i], v[3][i] }); });
        }
        return aabb;
    }

    void build(std::uint32_t max_depth = 8)
    {
        sortTetrahedrons();
        m_indices.clear();
        m_max_depth = max_depth;

        struct NodeTemplate {
            Node node;
            std::vector<std::uint32_t> indices;
        };

        // init nodes
        std::vector<NodeTemplate> nodes;
        nodes.reserve(m_items.size());
        nodes.push_back({});
        nodes[0].indices.resize(m_items.size());
        std::iota(nodes[0].indices.begin(), nodes[0].indices.end(), 0);

        // calculate aabb
        m_aabb = calculateAABB(nodes[0].indices);

        std::size_t currentNodeIdx = 0;

        const auto max_leafs_number = ((1 << (max_depth + 1)) + 1) / 2; // 2^(max_depth + 1)
        auto number_of_leafs = 0;

        while (currentNodeIdx < nodes.size()) {
            auto& cnode = nodes[currentNodeIdx].node;
            auto& cind = nodes[currentNodeIdx].indices;
            auto fom = splitNode(cnode, cind);

            if (fom == cind.size() || cind.size() < 2 || number_of_leafs >= max_leafs_number) {
                // leaf
                cnode.setLeaf();
                cnode.split_nelements.nelements = static_cast<std::uint32_t>(cind.size());
                cnode.setOffset(static_cast<std::uint32_t>(m_indices.size()));
                for (auto i : cind)
                    m_indices.push_back(i);
                cind.clear();
                cind.shrink_to_fit();
                ++number_of_leafs;
            } else {
                cnode.setOffset(nodes.size());
                NodeTemplate left, right;
                for (auto idx : cind) {
                    const auto& item = m_items[idx];
                    auto side = planeSide(item, cnode.dim(), cnode.split_nelements.split);
                    if (side <= 0)
                        left.indices.push_back(idx);
                    if (side >= 0)
                        right.indices.push_back(idx);
                }
                // purge memory
                cind.clear();
                cind.shrink_to_fit();
                nodes.push_back(left);
                nodes.push_back(right);
            }
            ++currentNodeIdx;
        }
        m_indices.shrink_to_fit();
        m_nodes.clear();
        m_nodes.resize(nodes.size());
        std::transform(std::execution::par_unseq, nodes.cbegin(), nodes.cend(), m_nodes.begin(), [](const auto& nr) {
            return nr.node;
        });
    }

    std::uint32_t splitAxis(const std::vector<Tetrahedron>& items, const std::vector<std::uint32_t>& indices)
    {
        // finding aabb
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };

        for (auto idx : indices) {
            const auto& item = items[idx];
            const auto aabb_tri = item.AABB();
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };
        return vectormath::argmax3<std::uint32_t>(extent);
    }

    float splitPlane(const std::vector<Tetrahedron>& items, const std::vector<std::uint32_t>& indices, const std::uint32_t dim)
    {
        const auto N = indices.size();
        std::vector<float> vals;
        vals.reserve(N);

        for (auto idx : indices) {
            const auto& item = items[idx];
            const auto v = item.center();
            vals.push_back(static_cast<float>(v[dim]));
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) / 2;
        }
    }

    static int planeSide(const Tetrahedron& item, const std::uint32_t D, const float plane)
    {
        const auto aabb = item.AABB();
        const auto min = aabb[D];
        const auto max = aabb[D + 3];

        if (max - plane < -epsilon())
            return -1;
        if (min - plane > epsilon())
            return 1;
        return 0;
    }

    int figureOfMerit(const std::vector<Tetrahedron>& items, const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx : indices) {
            const auto& item = items[idx];
            const auto side = planeSide(item, dim, planesep);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

    constexpr static double epsilon()
    {
        // Huristic epsilon
        return 11 * std::numeric_limits<double>::epsilon();
    }

private:
    struct Node {
        struct {
            std::uint32_t dim : 2; // dimensjon for branch
            std::uint32_t offset : 29; // offset to first child (branch) or to first element (leaf)
            std::uint32_t flag : 1; // Node is leaf (1) or branch (0)
        } dim_offset_flag;

        union {
            float split = 0; // Union split is for branches
            std::uint32_t nelements; // Union nelements is for number of items
        } split_nelements;

        Node()
        {
            dim_offset_flag.dim = 0;
            dim_offset_flag.offset = 0;
            dim_offset_flag.flag = 0;
        }

        std::uint32_t dim() const
        {
            return dim_offset_flag.dim;
        }
        void setDim(std::uint32_t dim)
        {
            dim_offset_flag.dim = dim;
        }

        void setLeaf()
        {
            dim_offset_flag.flag = std::uint32_t { 1 };
        }
        std::uint32_t isLeaf() const
        {
            return dim_offset_flag.flag;
        }

        void setOffset(std::uint32_t offset)
        {
            dim_offset_flag.offset = offset;
        }
        std::uint32_t offset() const
        {
            return dim_offset_flag.offset;
        }
    };

    int splitNode(Node& node, const std::vector<std::uint32_t>& indices) const
    {

        const auto aabb = calculateAABB(indices);
        const std::array<double, 3> extent = {
            aabb[3] - aabb[0],
            aabb[4] - aabb[1],
            aabb[5] - aabb[2],
        };
        const auto split_axis = vectormath::argmax3<std::uint32_t>(extent);
        node.setDim(split_axis);

        // split value
        std::vector<float> vals(indices.size());
        std::transform(std::execution::par_unseq, indices.cbegin(), indices.cend(), vals.begin(), [this, split_axis](auto i) {
            const auto& item = this->m_items[i];
            return static_cast<float>(item.vertices()[0][split_axis]);
        });
        std::sort(std::execution::par_unseq, vals.begin(), vals.end());
        const auto N = vals.size();
        float split_val = 0;
        if (N % 2 == 1)
            split_val = vals[N / 2];
        else
            split_val = (vals[N / 2] + vals[N / 2 - 1]) * 0.5;

        node.split_nelements.split = split_val;

        std::atomic<int> shared = 0;
        std::atomic<int> side = 0;
        std::for_each(std::execution::par_unseq, indices.cbegin(), indices.cend(), [this, &shared, &side, split_axis, split_val](const auto idx) {
            const auto& item = m_items[idx];
            const auto s = planeSide(item, split_axis, split_val);
            side += s;
            if (s == 0)
                shared++;
        });
        int fom = shared + std::abs(side);
        return fom;
    }

    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::uint32_t m_max_depth = 8;
    std::vector<std::uint32_t> m_indices;
    std::vector<Tetrahedron> m_items;
    std::vector<Node> m_nodes;
};
}
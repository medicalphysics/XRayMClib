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

#pragma once
#include "xraymc/particle.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/kdtreeintersectionresult.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreetype.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <memory>
#include <optional>
#include <vector>

namespace xraymc {

template <MeshKDTreeType U>
class MeshKDTreeFlat {
public:
    MeshKDTreeFlat() { };
    MeshKDTreeFlat(const std::vector<U>& triangles, std::uint32_t max_depth = 8)
    {
        setData(triangles, max_depth);
    }
    void setData(const std::vector<U>& triangles, std::uint32_t max_depth = 8)
    {
        m_indices.clear();
        m_nodes.clear();
        build(triangles, max_depth);
    }

    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (auto& node : m_nodes) {
            if (!node.isLeaf()) {
                const auto dim = node.dim();
                node.split_nelements.split += dist[dim];
            }
        }
    }
    void scale(double s)
    {
        for (auto& node : m_nodes) {
            if (!node.isLeaf()) {
                node.split_nelements.split *= s;
            }
        }
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::vector<U>& triangles, const std::array<double, 6>& aabb) const
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, triangles, *inter) : KDTreeIntersectionResult<const U> {};
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::vector<U>& items, const std::array<double, 2>& tboxAABB) const
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
        KDTreeIntersectionResult<const U> res;
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
                const auto& item = items[m_indices[idx]];
                auto t_cand = item.intersect(particle);
                if (t_cand) {
                    if (0 <= *t_cand && *t_cand < res.intersection) {
                        if (tbox[0] <= *t_cand && *t_cand <= tbox[1]) {
                            res.intersection = *t_cand;
                            res.item = &item;
                            res.rayOriginIsInsideItem = vectormath::dot(particle.dir, item.planeVector()) > 0;
                        }
                    }
                }
            }
            if (res.valid())
                stack.clear();
        }
        return res;
    }

protected:
    void build(const std::vector<U>& items, std::uint32_t max_depth = 8)
    {
        m_max_depth = max_depth;
        std::vector<std::uint32_t> indices(items.size());
        std::iota(indices.begin(), indices.end(), 0);

        struct NodeTemplate {
            Node node;
            std::vector<std::uint32_t> indices;
        };

        std::vector<NodeTemplate> nodes(1);
        nodes.reserve(indices.size());

        nodes[0].indices = indices;
        std::size_t currentNodeIdx = 0;

        const auto max_leafs_number = ((1 << (max_depth + 1)) + 1) / 2; // 2^(max_depth + 1)
        auto number_of_leafs = 0;

        while (currentNodeIdx < nodes.size()) {
            auto& cnode = nodes[currentNodeIdx].node;
            auto& cind = nodes[currentNodeIdx].indices;
            auto split_dim = splitAxis(items, cind);
            auto split_val = splitPlane(items, cind, split_dim);
            auto fom = figureOfMerit(items, cind, split_dim, split_val);
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
                // branch
                cnode.setDim(split_dim);
                cnode.split_nelements.split = split_val;
                cnode.setOffset(nodes.size());
                NodeTemplate left, right;
                for (auto idx : cind) {
                    const auto& item = items[idx];
                    auto side = planeSide(item, split_dim, split_val);
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
        m_nodes.clear();
        m_nodes.reserve(nodes.size());
        for (const auto& n : nodes)
            m_nodes.push_back(n.node);
    }

    std::uint32_t splitAxis(const std::vector<U>& items, const std::vector<std::uint32_t>& indices)
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

    float splitPlane(const std::vector<U>& items, const std::vector<std::uint32_t>& indices, const std::uint32_t dim)
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

    int planeSide(const U& item, const std::uint32_t D, const float plane)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        const auto aabb = item.AABB();

        min = aabb[D];
        max = aabb[D + 3];

        if (max - plane < -epsilon())
            return -1;
        if (min - plane > epsilon())
            return 1;
        return 0;
    }

    int figureOfMerit(const std::vector<U>& items, const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
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

        std::uint32_t dim()
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
        std::uint32_t isLeaf()
        {
            return dim_offset_flag.flag;
        }

        void setOffset(std::uint32_t offset)
        {
            dim_offset_flag.offset = offset;
        }
        std::uint32_t offset()
        {
            return dim_offset_flag.offset;
        }
    };
    std::uint32_t m_max_depth = 8;
    std::vector<std::uint32_t> m_indices;
    std::vector<U> m_items;
    std::vector<Node> m_nodes;
};
}
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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/kdtreeintersectionresult.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worlditems/worlditemtype.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <variant>
#include <vector>

namespace xraymc {

template <WorldItemType... Us>
class KDTreeFlat {
public:
    KDTreeFlat() { };
    KDTreeFlat(std::vector<std::variant<Us...>*>& items, std::uint32_t max_depth = 8)
    {
        setData(items, max_depth);
    }

    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    void setData(std::vector<std::variant<Us...>*>& items, std::uint32_t max_depth = 8)
    {
        std::vector<std::uint32_t> indices(items.size());
        std::iota(indices.begin(), indices.end(), 0);
        m_items = items;
        m_indices.clear();
        m_nodes.clear();
        build(indices, max_depth);
    };

    std::array<double, 6> AABB() const
    {
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };

        for (const auto* item : m_items) {
            const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
                aabb[i + 3] = std::max(aabb[i + 3], aabb_tri[i + 3]);
            }
        }
        return aabb;
    }

    VisualizationIntersectionResult<std::variant<Us...>> intersectVisualization(const ParticleType auto& particle, const std::array<double, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersectVisualization(particle, *inter) : VisualizationIntersectionResult<std::variant<Us...>> {};
    }

    VisualizationIntersectionResult<std::variant<Us...>> intersectVisualization(const ParticleType auto& particle, const std::array<double, 2>& tboxAABB) const
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
        VisualizationIntersectionResult<std::variant<Us...>> res;
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
                const auto* item = m_items[m_indices[idx]];
                auto t_cand = std::visit([&particle](const auto& it) { return it.template intersectVisualization<std::variant<Us...>>(particle); }, *item);
                if (t_cand.valid()) {
                    if (t_cand.intersection <= res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res = t_cand;
                            res.item = item;
                        }
                    }
                }
            }
            if (res.valid())
                stack.clear();
        }
        return res;
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb)
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<std::variant<Us...>> {};
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 2>& tboxAABB)
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
        KDTreeIntersectionResult<std::variant<Us...>> res = { .item = nullptr, .intersection = std::numeric_limits<double>::max() };
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
                auto* item = m_items[m_indices[idx]];
                auto t_cand = std::visit([&particle](const auto& it) { return it.intersect(particle); }, *item);
                if (t_cand.valid()) {
                    if (t_cand.intersection <= res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res.intersection = t_cand.intersection;
                            res.item = item;
                            res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
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
    void build(std::vector<std::uint32_t>& indices, std::uint32_t max_depth = 8)
    {
        m_max_depth = max_depth;
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
            const auto [split_dim, split_val] = splitAxisPlane(cind);
            auto fom = figureOfMerit(cind, split_dim, split_val);
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
                    auto* item = m_items[idx];
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

    std::pair<std::uint32_t, float> splitAxisPlane(const std::vector<std::uint32_t>& indices)
    {
        if (indices.size() == 0)
            return std::make_pair(std::uint32_t { 0 }, float { 0 });

        // split where we find best separation between objects
        // finding AA segments
        std::array<std::vector<std::pair<float, float>>, 3> segs;
        for (auto idx : indices) {
            const auto& u = m_items[idx];
            const auto aabb = std::visit([](const auto& item) { return item.AABB(); }, *u);
            for (std::uint32_t i = 0; i < 3; ++i) {
                const auto seg = std::make_pair(static_cast<float>(aabb[i]), static_cast<float>(aabb[i + 3]));
                segs[i].push_back(seg);
            }
        }

        // sorting segments
        std::uint32_t best_dim = 0;
        float best_plane = 0;
        int best_fom = static_cast<int>(indices.size());

        for (std::uint32_t i = 0; i < 3; ++i) {
            std::sort(segs[i].begin(), segs[i].end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });
            float max = std::numeric_limits<float>::lowest();
            float min = std::numeric_limits<float>::max();
            if (segs[i].size() > 1) {
                for (std::size_t idx = 0; idx < segs[i].size() - 1; ++idx) {
                    auto left = segs[i][idx].second;
                    auto right = segs[i][idx + 1].first;
                    if (left > right) {
                        left = segs[i][idx].first;
                        right = segs[i][idx + 1].second;
                    }
                    auto plane = (left + right) / 2;
                    auto cfom = figureOfMerit(indices, i, plane);
                    if (cfom < best_fom) {
                        best_fom = cfom;
                        best_plane = plane;
                        best_dim = i;
                    }
                    max = std::max(max, right);
                    min = std::min(min, left);
                }
            } else {
                max = std::max(segs[i][0].first, segs[i][0].second);
                min = std::min(segs[i][0].first, segs[i][0].second);
                float plane = (segs[i][0].first + segs[i][0].second) / 2;
                auto cfom = figureOfMerit(indices, i, plane);
                if (cfom < best_fom) {
                    best_fom = cfom;
                    best_plane = plane;
                    best_dim = i;
                }
            }
            if (best_fom > 0) {
                // lets also test the naive middle point
                float plane = (max + min) / 2;
                auto cfom = figureOfMerit(indices, i, plane);
                if (cfom < best_fom) {
                    best_fom = cfom;
                    best_plane = plane;
                    best_dim = i;
                }
            }
        }

        return std::make_pair(best_dim, best_plane);
    }

    int planeSide(std::variant<Us...>* item, const std::uint32_t D, const float plane)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        const auto aabb = std::visit([](const auto& it) { return it.AABB(); }, *item);

        min = aabb[D];
        max = aabb[D + 3];

        if (max - plane < -epsilon())
            return -1;
        if (min - plane > epsilon())
            return 1;
        return 0;
    }

    int figureOfMerit(const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx : indices) {
            auto* item = m_items[idx];
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
    std::vector<std::variant<Us...>*> m_items;
    std::vector<Node> m_nodes;
};
}
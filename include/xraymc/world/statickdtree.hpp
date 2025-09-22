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

#include <array>
#include <concepts>
#include <execution>
#include <numeric>
#include <optional>
#include <vector>

namespace xraymc {

template <typename U>
concept StaticKDTreeType = requires(U u, std::array<double, 3> vec, Particle p) {
    u.translate(vec);
    {
        u.intersect(p)
    } -> std::same_as<WorldIntersectionResult>;
    {
        u.center()
    } -> std::convertible_to<std::array<double, 3>>;
    {
        u.AABB()
    } -> std::convertible_to<std::array<double, 6>>;
};

template <StaticKDTreeType U>
class StaticKDTree {
public:
    StaticKDTree() { }
    StaticKDTree(const std::vector<U>& items, std::uint32_t max_depth = 8)
    {
        setData(items, max_depth);
    }

    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    void setData(const std::vector<U>& items, std::uint32_t max_depth = 8)
    {
        m_items = items;
        m_indices.clear();
        m_nodes.clear();

        build(max_depth);
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (auto& item : m_items)
            item.translate(dist);
        for (auto& node : m_nodes) {
            if (node.isLeaf()) {
                const auto dim = node.dim();
                node.split_nelements.split += dist[dim];
            }
        }
    }

    std::array<double, 6> AABB() const
    {
        std::array<double, 6> aabb = { 0, 0, 0, 0, 0, 0 };
        if (m_items.size() > 0) {
            aabb = m_items[0].AABB();
            for (auto& item : m_items) {
                const auto aabb_cand = item.AABB();
                for (int i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], aabb_cand[i]);
                    aabb[i + 3] = std::min(aabb[i + 3], aabb_cand[i + 3]);
                }
            }
        }
        return aabb;
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb) const
    {
        const auto tbox = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return tbox ? intersect(particle, *tbox) : KDTreeIntersectionResult<const U> {};
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::array<double, 2>& tboxAABB) const
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
                const auto& item = m_items[m_indices[idx]];
                auto t_cand = item.intersect(particle);
                if (t_cand.valid()) {
                    if (t_cand.intersection <= res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res.intersection = t_cand.intersection;
                            res.item = &item;
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
    void build(std::uint32_t max_depth = 8)
    {
        m_max_depth = max_depth;
        std::vector<std::uint32_t> indices(m_items.size());
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
                    const auto& item = m_items[idx];
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
            const auto aabb = u.AABB();
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

    int figureOfMerit(const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx : indices) {
            const auto& item = m_items[idx];
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

// Depricated

/*
template <int Depth, StaticKDTreeType U>
class StaticKDTree {
public:
    StaticKDTree() { }
    StaticKDTree(std::vector<U>& data)
    {
        if (data.size() == 0)
            return;

        const auto [D, plane] = planeSplit(data);

        m_D = D;
        m_plane = plane;

        // moving object depending on plane splits
        std::vector<U> left, right;
        for (std::size_t i = 0; i < data.size(); ++i) {
            const auto side = planeSide(data[i], m_D, m_plane);
            if (side == -1) {
                left.push_back(data[i]);
            } else if (side == 1) {
                right.push_back(data[i]);
            } else {
                left.push_back(data[i]);
                right.push_back(data[i]);
            }
        }

        m_left = StaticKDTree<Depth - 1, U>(left);
        m_right = StaticKDTree<Depth - 1, U>(right);
    }

    void translate(const std::array<double, 3>& dir)
    {
        m_plane += dir[m_D];
        m_left.translate(dir);
        m_right.translate(dir);
    }

    std::array<double, 6> AABB() const
    {
        auto left = m_left.AABB();
        const auto right = m_right.AABB();
        for (std::size_t i = 0; i < 3; ++i) {
            left[i] = std::min(left[i], right[i]);
        }
        for (std::size_t i = 3; i < 6; ++i) {
            left[i] = std::max(left[i], right[i]);
        }
        return left;
    }
    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb) const
    {
        const auto tbox = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return tbox ? intersect(particle, *tbox) : KDTreeIntersectionResult<const U> {};
    }
    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::array<double, 2>& tbox) const
    {
        // test for parallell beam, if parallell we must test both sides.
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<double>::epsilon()) {
            const auto hit_left = m_left.intersect(particle, tbox);
            const auto hit_right = m_right.intersect(particle, tbox);
            if (hit_left.valid() && hit_right.valid())
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            return hit_left.valid() ? hit_left : hit_right;
        }

        const auto [front, back] = particle.dir[m_D] > 0 ? std::make_pair(&m_left, &m_right) : std::make_pair(&m_right, &m_left);

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect(particle, tbox);
        }

        // both directions (start with front)
        const std::array<double, 2> t_front { tbox[0], t };
        const auto hit = front->intersect(particle, t_front);
        if (hit.valid()) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<double, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

protected:
    static std::pair<std::uint_fast32_t, double> planeSplit(const std::vector<U>& data)
    {
        // split where we find best separation between objects

        std::array<std::vector<std::pair<double, double>>, 3> segs;
        for (const auto& u : data) {
            const auto aabb = u.AABB();
            for (std::size_t i = 0; i < 3; ++i) {
                const auto seg = std::make_pair(aabb[i], aabb[i + 3]);
                bool added = false;
                std::size_t teller = 0;
                while (!added && teller != segs[i].size()) {
                    auto& r_seg = segs[i][teller];
                    if (r_seg.second > seg.first || r_seg.first > seg.second) {
                        r_seg.first = std::min(r_seg.first, seg.first);
                        r_seg.second = std::min(r_seg.second, seg.second);
                        added = true;
                    }
                    ++teller;
                }
                if (!added) {
                    segs[i].push_back(seg);
                }
            }
        }

        const std::array<std::size_t, 3> n_splits = { segs[0].size(), segs[1].size(), segs[2].size() };
        const auto n_splits_max = std::max(n_splits[0], std::max(n_splits[1], n_splits[2]));

        if (n_splits_max == 1) {
            const std::uint_fast32_t D = 0;
            const auto plane = std::nextafter(segs[D][0].second, std::numeric_limits<double>::max());
            return std::make_pair(D, plane);
        } else {
            // finding dim with min 2 splits and max extent
            std::uint_fast32_t D = 0;
            double extent = 0;
            for (std::uint_fast32_t i = 0; i < 3; ++i) {
                if (n_splits[i] == n_splits_max) {
                    auto min = std::numeric_limits<double>::max();
                    auto max = std::numeric_limits<double>::lowest();
                    for (const auto& [fi, se] : segs[i]) {
                        min = std::min(min, fi);
                        max = std::max(max, se);
                    }
                    const auto ex_cand = max - min;
                    if (ex_cand > extent) {
                        D = i;
                        extent = ex_cand;
                    }
                }
            }

            const auto ind_split = n_splits_max / 2;
            const auto plane = (segs[D][ind_split - 1].second + segs[D][ind_split].first) * 0.5;
            return std::make_pair(D, plane);
        }
    }

    static int planeSide(const U& obj, std::uint_fast32_t D, const double planesep)
    {
        const auto& aabb = obj.AABB();
        const auto min = aabb[D];
        const auto max = aabb[D + 3];

        if (max < planesep)
            return -1;
        if (min > planesep)
            return 1;
        return 0;
    }

private:
    std::uint_fast32_t m_D = 0;
    double m_plane = 0;
    StaticKDTree<Depth - 1, U> m_left;
    StaticKDTree<Depth - 1, U> m_right;
};

template <StaticKDTreeType U>
class StaticKDTree<0, U> {
public:
    StaticKDTree() { }
    StaticKDTree(std::vector<U>& data)
        : m_data(data)
    {
    }
    void translate(const std::array<double, 3>& dir)
    {
        for (auto& u : m_data)
            u.translate(dir);
    }
    std::array<double, 6> AABB() const noexcept
    {

        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (const auto& u : m_data) {
            const auto& aabb_obj = u.AABB();
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_obj[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_obj[i]);
            }
        }

        return aabb;
    }

    static constexpr std::size_t depth()
    {
        return 0;
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::array<double, 2>& tbox) const noexcept
    {
        KDTreeIntersectionResult<const U> res { .item = nullptr, .intersection = std::numeric_limits<double>::max() };

        for (const auto& u : m_data) {
            const auto t_cand = u.intersect(particle);
            if (t_cand.valid())
                if (t_cand.intersection < res.intersection) {
                    res.intersection = t_cand.intersection;
                    res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                    res.item = &u;
                }
        }
        return res;
    }

private:
    std::vector<U> m_data;
};*/
}
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

/**
 * @brief Flat-array KD-tree accelerator over heterogeneous world items for transport
 *        and visualization intersection queries.
 *
 * Organises a collection of world item variant pointers into a binary spatial tree
 * stored as a contiguous node array. Each internal node packs a split axis (2 bits),
 * a child offset (29 bits), and a leaf/branch flag (1 bit) into 32 bits, with the
 * split coordinate or leaf element count in a parallel 32-bit union. Both transport
 * (`intersect`) and visualization (`intersectVisualization`) queries use the same
 * iterative front-to-back traversal with a fixed-size stack (depth ≤ 32).
 *
 * The split plane is chosen by evaluating candidate planes derived from AABB
 * boundaries between adjacent items along all three axes, selecting the plane that
 * minimises the figure of merit (imbalance + straddle count). This differs from the
 * pure median-cut used by MeshKDTreeFlat and can produce better partitions for
 * spatially irregular world items.
 *
 * @tparam Us... World item types, each satisfying WorldItemType. Items are held as
 *              pointers to `std::variant<Us...>` and must remain valid for the
 *              lifetime of the tree.
 */
template <WorldItemType... Us>
class KDTreeFlat {
public:
    /// @brief Constructs an empty tree with no items.
    KDTreeFlat() { };

    /**
     * @brief Constructs a tree from a vector of world item pointers.
     * @param items     Pointers to world item variants to index.
     * @param max_depth Controls the maximum number of leaves (up to 2^(max_depth+1)).
     */
    KDTreeFlat(std::vector<std::variant<Us...>*>& items, std::uint32_t max_depth = 8)
    {
        setData(items, max_depth);
    }

    /// @brief Returns the maximum depth used when the tree was last built.
    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    /**
     * @brief Rebuilds the tree from @p items.
     *
     * Copies the item pointers into m_items, resets the index and node arrays, and
     * calls build() to construct the flat node array.
     * @param items     Pointers to world item variants to index.
     * @param max_depth Controls the maximum number of leaves (up to 2^(max_depth+1)).
     */
    void setData(std::vector<std::variant<Us...>*>& items, std::uint32_t max_depth = 8)
    {
        std::vector<std::uint32_t> indices(items.size());
        std::iota(indices.begin(), indices.end(), 0);
        m_items = items;
        m_indices.clear();
        m_nodes.clear();
        build(indices, max_depth);
    };

    /// @brief Returns the axis-aligned bounding box of all items as {xmin,ymin,zmin,xmax,ymax,zmax} in cm.
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

    /**
     * @brief Tests a particle ray for visualization intersection, guarded by an AABB pre-test.
     *
     * Returns an invalid result immediately if the ray misses @p aabb; otherwise
     * the AABB hit interval is forwarded to the iterative visualization traversal.
     * @param particle Particle whose position and direction define the ray.
     * @param aabb     Bounding box of all items as {xmin,ymin,zmin,xmax,ymax,zmax} in cm.
     * @return Closest visualization intersection result, or an invalid result on miss.
     */
    VisualizationIntersectionResult<std::variant<Us...>> intersectVisualization(const ParticleType auto& particle, const std::array<double, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersectVisualization(particle, *inter) : VisualizationIntersectionResult<std::variant<Us...>> {};
    }

    /**
     * @brief Iterative ray–tree traversal for visualization intersection within @p tboxAABB.
     *
     * Uses a fixed-size stack (depth ≤ 32). At each internal node the front child
     * (relative to the ray direction) is traversed first; the back child is pushed
     * onto the stack. Rays parallel to the split axis visit both children. At leaf
     * nodes each item's `intersectVisualization()` is called via `std::visit` and the
     * closest hit within the current interval is recorded. Once a hit is found the
     * stack is cleared.
     * @param particle  Particle whose position and direction define the ray.
     * @param tboxAABB  Initial ray interval {t_near, t_far} from the AABB intersection.
     * @return Closest visualization intersection result, or an invalid result.
     */
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

    /**
     * @brief Tests a particle ray for transport intersection, guarded by an AABB pre-test.
     *
     * Returns an invalid result immediately if the ray misses @p aabb; otherwise
     * the AABB hit interval is forwarded to the iterative transport traversal.
     * @param particle Particle whose position and direction define the ray.
     * @param aabb     Bounding box of all items as {xmin,ymin,zmin,xmax,ymax,zmax} in cm.
     * @return Closest transport intersection result, or an invalid result on miss.
     */
    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb)
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<std::variant<Us...>> {};
    }

    /**
     * @brief Iterative ray–tree traversal for transport intersection within @p tboxAABB.
     *
     * Identical in structure to the visualization overload but calls each item's
     * `intersect()` via `std::visit`, which returns a transport intersection result
     * including the `rayOriginIsInsideItem` flag.
     * @param particle  Particle whose position and direction define the ray.
     * @param tboxAABB  Initial ray interval {t_near, t_far} from the AABB intersection.
     * @return Closest transport intersection result, or an invalid result.
     */
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
    /**
     * @brief Builds the flat node array from an index list using a breadth-first strategy.
     *
     * Iterates over a working list of NodeTemplate entries. For each entry
     * `splitAxisPlane()` selects the best axis and plane; if the figure of merit
     * equals the subset size, the subset has fewer than two items, or the leaf budget
     * is exhausted, a leaf is created and its indices are appended to m_indices.
     * Otherwise a branch is created and two child NodeTemplates are appended.
     * @param indices   Indices into m_items for the full item set.
     * @param max_depth Controls the maximum number of leaves (up to 2^(max_depth+1)).
     */
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
                cnode.setOffset(static_cast<std::uint32_t>(nodes.size()));
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

    /**
     * @brief Selects the best split axis and plane by evaluating AABB-boundary candidates.
     *
     * For each of the three axes, item AABBs are sorted by their minimum coordinate.
     * Candidate planes are placed midway between adjacent AABB boundaries (left-max /
     * right-min pairs), and the midpoint of the overall extent is also tested. The
     * candidate with the lowest figure-of-merit across all axes is returned.
     * @param indices Indices into m_items for the current partition.
     * @return Pair of {best_axis (0–2), best_split_plane} in cm.
     */
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

    /**
     * @brief Determines which side of a split plane an item lies on, using its AABB.
     * @param item  Pointer to the world item variant to test.
     * @param D     Split axis index: 0 = x, 1 = y, 2 = z.
     * @param plane Split plane coordinate along @p D.
     * @return -1 if the item is entirely left of the plane, +1 if entirely right,
     *         or 0 if it straddles the plane (within epsilon).
     */
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

    /**
     * @brief Evaluates the quality of a candidate split plane.
     *
     * Sums signed side values (+1 right, −1 left) and adds the count of straddling
     * items. A lower value indicates a more balanced partition. If the result equals
     * indices.size() the split offers no benefit.
     * @param indices   Indices into m_items for the current partition.
     * @param dim       Split axis index.
     * @param planesep  Candidate split plane coordinate.
     * @return Figure-of-merit score (lower is better).
     */
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

    /// @brief Returns the heuristic epsilon used for plane-side comparisons (11 × machine epsilon).
    constexpr static double epsilon()
    {
        // Huristic epsilon
        return 11 * std::numeric_limits<double>::epsilon();
    }

private:
    /**
     * @brief Packed KD-tree node storing branch or leaf data in 64 bits.
     *
     * The first 32-bit word is a bitfield: dim (2 bits) holds the split axis for
     * branch nodes, offset (29 bits) is the index of the first child node (branch)
     * or the start of the element index range in m_indices (leaf), and flag (1 bit)
     * distinguishes leaf (1) from branch (0).
     *
     * The second 32-bit word is a union: `split` (float) holds the split plane
     * coordinate for branch nodes; `nelements` (uint32) holds the element count
     * for leaf nodes.
     */
    struct Node {
        struct {
            std::uint32_t dim : 2;     ///< Split axis index (0–2); valid for branch nodes only.
            std::uint32_t offset : 29; ///< Index of first child (branch) or first index entry (leaf).
            std::uint32_t flag : 1;    ///< 1 = leaf node, 0 = branch node.
        } dim_offset_flag;

        union {
            float split = 0;           ///< Split plane coordinate along dim; used by branch nodes.
            std::uint32_t nelements;   ///< Number of item indices in this leaf.
        } split_nelements;

        /// @brief Constructs a default branch node (dim=0, offset=0, flag=0).
        Node()
        {
            dim_offset_flag.dim = 0;
            dim_offset_flag.offset = 0;
            dim_offset_flag.flag = 0;
        }

        /// @brief Returns the split axis index stored in this branch node.
        std::uint32_t dim()
        {
            return dim_offset_flag.dim;
        }

        /// @brief Sets the split axis index for this branch node.
        void setDim(std::uint32_t dim)
        {
            dim_offset_flag.dim = dim;
        }

        /// @brief Marks this node as a leaf.
        void setLeaf()
        {
            dim_offset_flag.flag = std::uint32_t { 1 };
        }

        /// @brief Returns non-zero if this node is a leaf.
        std::uint32_t isLeaf()
        {
            return dim_offset_flag.flag;
        }

        /// @brief Sets the child/element offset for this node.
        void setOffset(std::uint32_t offset)
        {
            dim_offset_flag.offset = offset;
        }

        /// @brief Returns the child/element offset stored in this node.
        std::uint32_t offset()
        {
            return dim_offset_flag.offset;
        }
    };

    std::uint32_t m_max_depth = 8;                      ///< Maximum depth used when the tree was built.
    std::vector<std::uint32_t> m_indices;               ///< Flat list of item indices referenced by leaf nodes.
    std::vector<std::variant<Us...>*> m_items;          ///< Pointers to all world item variants; owned by the caller.
    std::vector<Node> m_nodes;                          ///< Flat node array; node 0 is the root.
};
}
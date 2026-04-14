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

/**
 * @brief Flat-array KD-tree accelerator for ray–mesh intersection queries.
 *
 * Organises a set of triangles (or any type satisfying MeshKDTreeType) into a binary
 * spatial tree stored as a contiguous node array rather than a pointer-linked structure.
 * Each internal node records a split axis (2 bits), a child offset (29 bits), and a
 * leaf/branch flag (1 bit) packed into a single 32-bit word, with the split coordinate
 * or leaf element count stored in a parallel 32-bit union. Traversal is iterative and
 * uses a fixed-size stack, avoiding recursion overhead.
 *
 * The tree is built with a median-cut strategy. A leaf is created when further splitting
 * would not reduce the figure of merit, the subset has fewer than two items, or the
 * maximum leaf budget (2^(max_depth+1) leaves) is exhausted.
 *
 * @tparam U Triangle-like type satisfying the MeshKDTreeType concept
 *           (must expose AABB(), center(), intersect(), and planeVector()).
 */
template <MeshKDTreeType U>
class MeshKDTreeFlat {
public:
    /// @brief Constructs an empty tree with no triangles.
    MeshKDTreeFlat() { };

    /**
     * @brief Constructs a tree over all triangles in @p triangles.
     * @param triangles Triangle collection; must remain valid for the lifetime of queries.
     * @param max_depth Controls the maximum number of leaves (up to 2^(max_depth+1)).
     */
    MeshKDTreeFlat(const std::vector<U>& triangles, std::uint32_t max_depth = 8)
    {
        setData(triangles, max_depth);
    }

    /**
     * @brief Rebuilds the tree over all triangles in @p triangles.
     * @param triangles Triangle collection.
     * @param max_depth Controls the maximum number of leaves (up to 2^(max_depth+1)).
     */
    void setData(const std::vector<U>& triangles, std::uint32_t max_depth = 8)
    {
        m_indices.clear();
        m_nodes.clear();
        build(triangles, max_depth);
    }

    /// @brief Returns the maximum depth used when the tree was last built.
    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    /**
     * @brief Translates the split-plane coordinates of all internal nodes by @p dist.
     *
     * Only the component along each node's stored split axis is applied; leaf nodes
     * are not modified.
     * @param dist Displacement vector in cm.
     */
    void translate(const std::array<double, 3>& dist)
    {
        for (auto& node : m_nodes) {
            if (!node.isLeaf()) {
                const auto dim = node.dim();
                node.split_nelements.split += dist[dim];
            }
        }
    }

    /**
     * @brief Uniformly scales the split-plane coordinates of all internal nodes by @p s.
     * @param s Scale factor applied to every split-plane value.
     */
    void scale(double s)
    {
        for (auto& node : m_nodes) {
            if (!node.isLeaf()) {
                node.split_nelements.split *= s;
            }
        }
    }

    /**
     * @brief Tests a particle ray against the mesh, guarded by an AABB pre-test.
     *
     * Returns an invalid result immediately if the ray misses @p aabb; otherwise
     * forwards the AABB hit interval to the iterative tree traversal.
     * @param particle  Particle whose position and direction define the ray.
     * @param triangles Triangle collection indexed by the tree's leaf nodes.
     * @param aabb      Axis-aligned bounding box of the full mesh as
     *                  {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     * @return Closest intersection within the AABB, or an invalid result.
     */
    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::vector<U>& triangles, const std::array<double, 6>& aabb) const
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, triangles, *inter) : KDTreeIntersectionResult<const U> {};
    }

    /**
     * @brief Iterative ray–tree traversal within the ray interval @p tboxAABB.
     *
     * Uses a fixed-size stack (depth ≤ 32) to avoid recursion. Internal nodes are
     * traversed front-to-back relative to the ray direction; the back child is pushed
     * onto the stack. Rays parallel to the split axis visit both children. At each
     * leaf the stored triangles are tested and the closest hit within the current
     * interval is recorded. Once a hit is found the stack is cleared to skip
     * remaining nodes.
     * @param particle  Particle whose position and direction define the ray.
     * @param items     Triangle collection indexed by the tree's leaf nodes.
     * @param tboxAABB  Initial ray interval {t_near, t_far} from the AABB intersection.
     * @return Closest intersection within the interval, or an invalid result.
     */
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
    /**
     * @brief Builds the flat node array from @p items using a breadth-first median-cut strategy.
     *
     * Iterates over a working list of `NodeTemplate` entries. For each entry the split
     * axis and plane are computed; if the figure of merit equals the subset size, the
     * subset has fewer than two items, or the leaf budget (2^(max_depth+1) / 2 leaves)
     * is exhausted, a leaf node is created and its indices are appended to m_indices.
     * Otherwise a branch node is created and two child NodeTemplates are appended.
     * After all nodes are processed the final `Node` values are copied into m_nodes.
     * @param items     Triangle collection to index.
     * @param max_depth Controls the maximum number of leaves.
     */
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

    /**
     * @brief Selects the split axis as the dimension with the largest AABB extent.
     * @param items   Triangle collection.
     * @param indices Indices of the triangles in the current partition.
     * @return Axis index (0 = x, 1 = y, 2 = z) with the widest span.
     */
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

    /**
     * @brief Computes the median-cut split plane along @p dim.
     *
     * Collects each triangle's centroid coordinate along @p dim, sorts them, and
     * returns the median (average of the two middle values for even N). The result
     * is stored as a `float` to match the packed Node representation.
     * @param items   Triangle collection.
     * @param indices Indices of the triangles in the current partition.
     * @param dim     Split axis index.
     * @return Median centroid coordinate along @p dim.
     */
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

    /**
     * @brief Determines which side of a split plane a triangle lies on, using its AABB.
     * @param item  Triangle to test.
     * @param D     Split axis index: 0 = x, 1 = y, 2 = z.
     * @param plane Split plane coordinate along @p D.
     * @return -1 if the triangle is entirely left of the plane, +1 if entirely right,
     *         or 0 if it straddles the plane (within epsilon).
     */
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

    /**
     * @brief Evaluates the quality of a candidate split plane.
     *
     * Sums the signed side values (+1 right, −1 left) and adds the count of
     * straddling triangles. A lower value indicates a more balanced partition.
     * If the result equals indices.size() the split is not useful.
     * @param items    Triangle collection.
     * @param indices  Indices of the triangles in the current partition.
     * @param dim      Split axis index.
     * @param planesep Candidate split plane coordinate.
     * @return Figure-of-merit score (lower is better).
     */
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
            std::uint32_t dim : 2;    ///< Split axis index (0–2); valid for branch nodes only.
            std::uint32_t offset : 29; ///< Index of first child (branch) or first index entry (leaf).
            std::uint32_t flag : 1;   ///< 1 = leaf node, 0 = branch node.
        } dim_offset_flag;

        union {
            float split = 0;           ///< Split plane coordinate along dim; used by branch nodes.
            std::uint32_t nelements;   ///< Number of triangle indices in this leaf.
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
    std::uint32_t m_max_depth = 8;      ///< Maximum depth used when the tree was built.
    std::vector<std::uint32_t> m_indices; ///< Flat list of triangle indices referenced by leaf nodes.
    std::vector<U> m_items;              ///< (Unused storage; triangles are owned by the caller.)
    std::vector<Node> m_nodes;           ///< Flat node array; node 0 is the root.
};
}
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
#include "xraymc/world/basicshapes/tetrahedron.hpp"
#include "xraymc/world/kdtreeintersectionresult.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <memory>
#include <optional>
#include <vector>

namespace xraymc {

/**
 * @brief A KD-tree over the outer-contour triangles of a tetrahedral mesh for fast
 *        ray–surface intersection.
 *
 * This tree accelerates the ray-entry test for a TetrahedalMesh: only the triangular
 * faces that form the outer surface (contour) of the mesh are stored. When a particle
 * ray enters the mesh AABB, the tree finds the first contour triangle hit and determines
 * whether the ray origin is inside or outside the mesh (via the face normal direction).
 * The tree is built by median-cut spatial partitioning along the longest axis of the
 * current node's bounding box. Each node stores either a split plane (branch) or a list
 * of triangle indices (leaf). Indices and nodes are stored in flat arrays for cache
 * efficiency.
 */
class TetrahedalMeshContourKDTree {
public:
    /// @brief Default constructor. Creates an empty tree.
    TetrahedalMeshContourKDTree() { };

    /**
     * @brief Constructs the KD-tree immediately from the provided surface triangles.
     * @param vertices  Vertex position array (indexed by element vertex indices).
     * @param elements  Triangle index triples defining the contour faces.
     * @param max_depth Maximum depth of the KD-tree; deeper trees are faster to query
     *                  but use more memory.
     */
    TetrahedalMeshContourKDTree(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, std::uint32_t max_depth = 8)
    {
        setData(vertices, elements, max_depth);
    }

    /**
     * @brief (Re-)builds the KD-tree from the provided surface triangles.
     *        Any previously stored tree data is discarded.
     * @param vertices  Vertex position array.
     * @param elements  Triangle index triples defining the contour faces.
     * @param max_depth Maximum depth of the KD-tree.
     */
    void setData(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, std::uint32_t max_depth = 8)
    {
        m_indices.clear();
        m_nodes.clear();
        build(vertices, elements, max_depth);
    }

    /// @brief Returns the maximum depth the tree was built with.
    std::uint32_t maxDepth() const
    {
        return m_max_depth;
    }

    /**
     * @brief Translates all branch split planes by @p dist.
     *
     * Only the split value of internal (branch) nodes is updated; leaf node data
     * (triangle indices) stays in place because it refers to the external vertex array.
     * @param dist Displacement vector in cm along {x, y, z}.
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
     * @brief AABB-filtered ray intersection: rejects the query immediately if the ray
     *        misses the overall bounding box, then delegates to the interval overload.
     * @param particle  Particle whose position and direction define the ray.
     * @param vertices  Vertex position array shared with the owning mesh.
     * @param elements  Triangle index triples for the contour faces.
     * @param aabb      Axis-aligned bounding box of the full mesh as {xmin,ymin,zmin,xmax,ymax,zmax}.
     * @return The nearest hit triangle and distance, or an invalid result if the ray misses.
     */
    KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> intersect(const ParticleType auto& particle, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::array<double, 6>& aabb) const
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, vertices, elements, *inter) : KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> {};
    }

    /**
     * @brief KD-tree traversal ray intersection within a pre-computed ray interval.
     *
     * Traverses the tree using an explicit fixed-size stack, testing each leaf's
     * triangles with a Möller–Trumbore ray–triangle test. The first valid hit within
     * the interval @p tboxAABB is returned. The rayOriginIsInsideItem flag is set by
     * comparing the hit triangle's normal against the ray direction.
     * @param particle   Particle whose position and direction define the ray.
     * @param vertices   Vertex position array.
     * @param elements   Triangle index triples for the contour faces.
     * @param tboxAABB   {tmin, tmax} ray interval from the AABB pre-filter.
     * @return The nearest hit triangle and distance within the interval, or invalid.
     */
    KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> intersect(const ParticleType auto& particle, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::array<double, 2>& tboxAABB) const
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
        KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> res;
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
                const auto& el = elements[m_indices[idx]];
                auto t_cand = basicshape::tetrahedron::intersectTriangle(vertices[el[0]], vertices[el[1]], vertices[el[2]], particle);
                if (t_cand) {
                    if (0 <= *t_cand && *t_cand < res.intersection) {
                        constexpr double epsilon = 1E-8;
                        if (tbox[0] - epsilon <= *t_cand && *t_cand <= tbox[1] + epsilon) {
                            res.intersection = *t_cand;
                            res.item = &el;
                            const auto normal = basicshape::tetrahedron::normalVector<false>(vertices[el[0]], vertices[el[1]], vertices[el[2]]);
                            res.rayOriginIsInsideItem = vectormath::dot(normal, particle.dir) > 0.0;
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
     * @brief Recursively builds the KD-tree using median-cut spatial partitioning.
     *
     * Each node is split along the longest axis of its bounding box at the median
     * triangle centroid. Triangles that straddle the split plane are duplicated into
     * both children. A node becomes a leaf when splitting would not reduce element
     * count, the node holds fewer than two triangles, or the maximum depth is reached.
     * @param vertices  Vertex position array.
     * @param elements  Triangle index triples for the contour faces.
     * @param max_depth Maximum recursion depth.
     */
    void build(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, std::uint32_t max_depth = 8)
    {
        m_max_depth = max_depth;
        std::vector<std::uint32_t> indices(elements.size());
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
            auto split_dim = splitAxis(vertices, elements, cind);
            auto split_val = splitPlane(vertices, elements, cind, split_dim);
            auto fom = figureOfMerit(vertices, elements, cind, split_dim, split_val);
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
                    const auto& item = elements[idx];
                    auto side = planeSide(vertices, item, split_dim, split_val);
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
     * @brief Selects the split axis as the longest spatial extent of the node's AABB.
     * @param vertices Vertex position array.
     * @param elements Triangle index triples.
     * @param indices  Indices of the triangles in the current node.
     * @return Axis index (0=x, 1=y, 2=z) with the greatest extent.
     */
    std::uint32_t splitAxis(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::vector<std::uint32_t>& indices)
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
            const auto& item = elements[idx];
            for (const auto vIdx : item) {
                for (std::size_t i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], vertices[vIdx][i]);
                    aabb[i + 3] = std::max(aabb[i + 3], vertices[vIdx][i]);
                }
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };
        return vectormath::argmax3<std::uint32_t>(extent);
    }

    /**
     * @brief Computes the median split position along @p dim from triangle centroid coordinates.
     * @param vertices Vertex position array.
     * @param elements Triangle index triples.
     * @param indices  Indices of the triangles in the current node.
     * @param dim      Axis index along which to split.
     * @return Median centroid coordinate along @p dim as a split plane position.
     */
    float splitPlane(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::vector<std::uint32_t>& indices, const std::uint32_t dim)
    {
        const auto N = indices.size();
        std::vector<float> vals;
        vals.reserve(N);

        for (auto idx : indices) {
            double pos = 0;
            for (auto vIdx : elements[idx])
                pos += vertices[vIdx][dim];
            vals.push_back(static_cast<float>(pos / 3.0));
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) / 2;
        }
    }

    /**
     * @brief Classifies a triangle relative to a split plane along axis @p D.
     * @param vertices Vertex position array.
     * @param element  Triangle vertex indices.
     * @param D        Axis index (0=x, 1=y, 2=z).
     * @param plane    Split plane coordinate along @p D.
     * @return -1 if the triangle is fully on the negative side, +1 if fully on the
     *         positive side, or 0 if it straddles the plane (within epsilon()).
     */
    int planeSide(const std::vector<std::array<double, 3>>& vertices, const std::array<std::uint32_t, 3>& element, const std::uint32_t D, const float plane)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        for (auto vIdx : element) {
            max = std::max(max, vertices[vIdx][D]);
            min = std::min(min, vertices[vIdx][D]);
        }

        if (max - plane < -epsilon())
            return -1;
        if (min - plane > epsilon())
            return 1;
        return 0;
    }

    /**
     * @brief Evaluates the quality of a proposed split as a figure of merit.
     *
     * Lower values indicate a more balanced split. The metric is the absolute imbalance
     * (|left_count - right_count|) plus the number of triangles that straddle the plane,
     * since straddling triangles must be duplicated into both children. A split that
     * produces no improvement returns the original triangle count.
     * @param vertices  Vertex position array.
     * @param elements  Triangle index triples.
     * @param indices   Indices of the triangles in the current node.
     * @param dim       Axis index of the candidate split plane.
     * @param planesep  Position of the candidate split plane along @p dim.
     * @return Figure-of-merit value; equal to indices.size() if the split is useless.
     */
    int figureOfMerit(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx : indices) {
            const auto side = planeSide(vertices, elements[idx], dim, planesep);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

    /**
     * @brief Returns the heuristic epsilon used for plane-side classification.
     *
     * Set to 11× machine epsilon to account for floating-point rounding when
     * comparing vertex coordinates against a split plane.
     * @return Tolerance value for plane-side tests.
     */
    constexpr static double epsilon()
    {
        // Huristic epsilon
        return 11 * std::numeric_limits<double>::epsilon();
    }

private:
    /**
     * @brief Compact KD-tree node packed into two 32-bit words.
     *
     * The first word encodes three fields in a bitfield:
     * - `dim`    (2 bits): split axis for branch nodes (0=x, 1=y, 2=z).
     * - `offset` (29 bits): index of the first child node (branch) or first
     *   element index in m_indices (leaf).
     * - `flag`   (1 bit): 1 if this node is a leaf, 0 if it is a branch.
     *
     * The second word is a union:
     * - `split` (float): split plane coordinate along `dim` for branch nodes.
     * - `nelements` (uint32): number of triangle indices stored in the leaf.
     */
    struct Node {
        struct {
            std::uint32_t dim : 2;    ///< Split axis index (branch only): 0=x, 1=y, 2=z.
            std::uint32_t offset : 29; ///< Child offset (branch) or index-array offset (leaf).
            std::uint32_t flag : 1;   ///< 1 = leaf node, 0 = branch node.
        } dim_offset_flag;

        union {
            float split = 0;           ///< Split plane coordinate along dim (branch nodes).
            std::uint32_t nelements;   ///< Number of triangle indices in this leaf.
        } split_nelements;

        /// @brief Default constructor. Initializes as a branch node with dim=0 and offset=0.
        Node()
        {
            dim_offset_flag.dim = 0;
            dim_offset_flag.offset = 0;
            dim_offset_flag.flag = 0;
        }

        /// @brief Returns the split axis index for this branch node.
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

        /// @brief Returns non-zero if this node is a leaf, zero if it is a branch.
        std::uint32_t isLeaf()
        {
            return dim_offset_flag.flag;
        }

        /// @brief Sets the child or element-array offset for this node.
        void setOffset(std::uint32_t offset)
        {
            dim_offset_flag.offset = offset;
        }

        /// @brief Returns the child or element-array offset for this node.
        std::uint32_t offset()
        {
            return dim_offset_flag.offset;
        }
    };
    std::uint32_t m_max_depth = 8;
    std::vector<std::uint32_t> m_indices;
    std::vector<Node> m_nodes;
};
}
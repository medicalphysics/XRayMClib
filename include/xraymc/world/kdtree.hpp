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

// DEPRICATED

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
 * @brief Recursive pointer-based KD-tree over heterogeneous world items.
 *
 * @deprecated Prefer the flat-array accelerator MeshKDTreeFlat for new code.
 *
 * Organises a collection of world item variant pointers into a binary spatial tree
 * using a median-cut strategy along the widest AABB axis. Each internal node stores
 * a splitting plane; leaf nodes hold the item pointers directly. The tree supports
 * both transport intersection queries (intersect()) and visualization queries
 * (intersectVisualization()). Both use the same front-to-back recursive traversal.
 *
 * @tparam Us... World item types, each satisfying WorldItemType. Items are held as
 *              pointers to `std::variant<Us...>` and must remain valid for the
 *              lifetime of the tree.
 */
template <WorldItemType... Us>
class KDTree {
public:
    /// @brief Constructs an empty tree with no items.
    KDTree() { }

    /**
     * @brief Constructs a tree from a vector of world item pointers.
     * @param items     Pointers to world item variants to index.
     * @param max_depth Maximum recursion depth (default 8).
     */
    KDTree(std::vector<std::variant<Us...>*>& items, const std::size_t max_depth = 8)
    {
        setData(items, max_depth);
    }

    /**
     * @brief Rebuilds the tree from @p items using a median-cut strategy.
     *
     * Computes the AABB of all items, selects the widest axis, finds the median-cut
     * plane, and partitions items into left and right children. Recursion stops when
     * the figure of merit equals the item count, @p max_depth reaches 1, or there is
     * only one item; those items become a leaf.
     * @param items     Pointers to world item variants to index.
     * @param max_depth Maximum remaining recursion depth.
     */
    void setData(std::vector<std::variant<Us...>*>& items, const std::size_t max_depth = 8)
    {
        m_items.clear();
        m_left = nullptr;
        m_right = nullptr;
        if (items.size() < 2) {
            for (const auto& item : items)
                m_items.push_back(item);
            return;
        }
        // finding aabb
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (auto& item : items) {
            const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<std::uint_fast32_t>(extent);

        const auto split = planeSplit(items);

        const auto fom = figureOfMerit(items, split);

        if (fom == items.size() || max_depth <= 1 || items.size() <= 1) {
            m_items = items;
        } else {
            m_plane = split;
            std::vector<std::variant<Us...>*> left;
            std::vector<std::variant<Us...>*> right;
            for (const auto& item : items) {
                const auto side = planeSide(item, m_plane, m_D);
                if (side <= 0)
                    left.push_back(item);
                if (side >= 0)
                    right.push_back(item);
            }
            m_left = std::make_unique<KDTree>(left, max_depth - 1);
            m_right = std::make_unique<KDTree>(right, max_depth - 1);
        }
    }
    /// @brief Returns the axis-aligned bounding box of all items in the tree as {xmin,ymin,zmin,xmax,ymax,zmax} in cm.
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
        AABB_iterator(aabb);
        return aabb;
    }

    /// @brief Returns the current depth of the tree (levels from root to deepest leaf).
    std::size_t depth() const
    {
        std::size_t teller = 0;
        depth_iterator(teller);
        return teller;
    }

    /**
     * @brief Translates the splitting planes of all nodes by @p dist.
     *
     * Only the component along each node's split axis is applied.
     * @param dist Displacement vector in cm.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_plane += dist[m_D];
        if (m_left) {
            m_left->translate(dist);
            m_right->translate(dist);
        }
    }

    /**
     * @brief Tests a particle ray against the tree, guarded by an AABB pre-test.
     *
     * Returns an invalid result immediately if the ray misses @p aabb; otherwise
     * the AABB hit interval is forwarded to the recursive traversal.
     * @param particle Particle whose position and direction define the ray.
     * @param aabb     Bounding box of the full set of items.
     * @return Closest intersection result, or an invalid result on miss.
     */
    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb)
    {
        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<std::variant<Us...>> { };
    }

    /**
     * @brief Recursive ray–tree traversal for transport intersection within @p tbox.
     *
     * At leaf nodes each stored item is tested via `std::visit` and the closest hit
     * within [tbox[0], tbox[1]] is returned. At internal nodes the split plane
     * determines front/back child order; the back child is only visited when the
     * front child yields no hit closer than the split plane.
     * @param particle Particle whose position and direction define the ray.
     * @param tbox     Active ray interval {t_near, t_far}.
     * @return Closest intersection result within the interval, or an invalid result.
     */
    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 2>& tbox)
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            KDTreeIntersectionResult<std::variant<Us...>> res = { .item = nullptr, .intersection = std::numeric_limits<double>::max() };
            for (auto& item : m_items) {
                // const auto t_cand = item->intersect(particle);
                const auto t_cand = std::visit([&particle](const auto& it) { return it.intersect(particle); }, *item);
                if (t_cand.valid()) {
                    const auto border_delta = t_cand.rayOriginIsInsideItem ? -GEOMETRIC_ERROR() : GEOMETRIC_ERROR();
                    if (t_cand.intersection + border_delta < res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res.intersection = t_cand.intersection;
                            res.item = item;
                            res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<double>::epsilon()) {
            const auto hit_left = m_left->intersect(particle, tbox);
            const auto hit_right = m_right->intersect(particle, tbox);
            if (hit_left.item && hit_right.item)
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            if (!hit_left.item)
                return hit_right;
            return hit_left;
        }

        auto front = particle.dir[m_D] > 0 ? m_left.get() : m_right.get();
        auto back = particle.dir[m_D] > 0 ? m_right.get() : m_left.get();

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
        if (hit.item) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<double, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

    /**
     * @brief Tests a particle ray for visualization intersection, guarded by an AABB pre-test.
     *
     * Returns an invalid result immediately if the ray misses @p aabb; otherwise
     * the AABB hit interval is forwarded to the recursive visualization traversal.
     * @param particle Particle whose position and direction define the ray.
     * @param aabb     Bounding box of the full set of items.
     * @return Closest visualization intersection result, or an invalid result on miss.
     */
    VisualizationIntersectionResult<std::variant<Us...>> intersectVisualization(const ParticleType auto& particle, const std::array<double, 6>& aabb) const
    {
        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersectVisualization(particle, *inter) : VisualizationIntersectionResult<std::variant<Us...>> { };
    }

    /**
     * @brief Recursive ray–tree traversal for visualization intersection within @p tbox.
     *
     * Identical in structure to the transport intersect() overload but calls each item's
     * `intersectVisualization()` via `std::visit`, which additionally returns a surface
     * normal and a scalar value (e.g. dose) for shading.
     * @param particle Particle whose position and direction define the ray.
     * @param tbox     Active ray interval {t_near, t_far}.
     * @return Closest visualization intersection result, or an invalid result.
     */
    VisualizationIntersectionResult<std::variant<Us...>> intersectVisualization(const ParticleType auto& particle, const std::array<double, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            VisualizationIntersectionResult<std::variant<Us...>> res;
            res.intersection = std::numeric_limits<double>::max();
            for (auto& item : m_items) {
                auto t_cand = std::visit([&particle](const auto& it) { return it.template intersectVisualization<std::variant<Us...>>(particle); }, *item);
                if (t_cand.valid()) {
                    if (t_cand.intersection < res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res = t_cand;
                            res.item = item;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<double>::epsilon()) {
            const auto hit_left = m_left->intersectVisualization(particle, tbox);
            const auto hit_right = m_right->intersectVisualization(particle, tbox);
            if (hit_left.item && hit_right.item)
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            if (!hit_left.item)
                return hit_right;
            return hit_left;
        }

        auto front = particle.dir[m_D] > 0 ? m_left.get() : m_right.get();
        auto back = particle.dir[m_D] > 0 ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersectVisualization(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersectVisualization(particle, tbox);
        }

        // both directions (start with front)
        const std::array<double, 2> t_front { tbox[0], t };
        const auto hit = front->intersectVisualization(particle, t_front);
        if (hit.item) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<double, 2> t_back { t, tbox[1] };
        return back->intersectVisualization(particle, t_back);
    }

protected:
    /**
     * @brief Computes the median-cut split plane along the current split axis.
     *
     * Collects the centroid coordinate along m_D for every item (via std::visit),
     * sorts them, and returns the median value (average of the two middle values
     * for even N).
     * @param items Items in the current partition.
     * @return Split plane coordinate along m_D in cm.
     */
    double planeSplit(const std::vector<std::variant<Us...>*>& items) const
    {
        const auto N = items.size();
        std::vector<double> vals;
        vals.reserve(N);

        for (const auto& item : items) {
            const auto v = std::visit([](const auto& it) { return it.center(); }, *item);
            vals.push_back(v[m_D]);
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) * 0.5;
        }
    }

    /**
     * @brief Evaluates the quality of a candidate split plane.
     *
     * Sums signed side values (+1 right, −1 left) over all items and adds the count
     * of straddling items. A lower value indicates a more balanced partition.
     * If the result equals items.size() the split offers no benefit.
     * @param items    Items in the current partition.
     * @param planesep Candidate split plane coordinate along m_D.
     * @return Figure-of-merit score (lower is better).
     */
    int figureOfMerit(const std::vector<std::variant<Us...>*>& items, const double planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto& item : items) {
            const auto side = planeSide(item, planesep, m_D);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

    /**
     * @brief Determines which side of a split plane an item lies on, using its AABB.
     * @param item  Pointer to the world item variant to test.
     * @param plane Split plane coordinate along axis @p D.
     * @param D     Split axis index: 0 = x, 1 = y, 2 = z.
     * @return -1 if the item is entirely left of the plane, +1 if entirely right,
     *         or 0 if it straddles the plane.
     */
    static int planeSide(std::variant<Us...>* item, const double plane, const unsigned int D)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        const auto aabb = std::visit([](const auto& it) { return it.AABB(); }, *item);

        min = aabb[D];
        max = aabb[D + 3];

        if (lessOrEqual(max, plane))
            // if (max - plane <= epsilon())
            return -1;
        if (greaterOrEqual(min, plane))
            // if (plane - min <= epsilon())
            return 1;
        return 0;
    }

    /// @brief Increments @p teller once per level by recursing into the left child.
    void depth_iterator(std::size_t& teller) const
    {
        teller++;
        if (m_left)
            m_left->depth_iterator(teller);
    }

    /**
     * @brief Collects all item pointers from leaf nodes into @p all (const overload).
     * @param all Output vector; leaf item pointers are appended via back_inserter.
     */
    void item_iterator(std::vector<const std::variant<Us...>*>& all) const
    {
        if (m_left) {
            m_left->item_iterator(all);
            m_right->item_iterator(all);
        } else {
            std::copy(m_items.cbegin(), m_items.cend(), std::back_inserter(all));
        }
    }

    /**
     * @brief Collects all item pointers from leaf nodes into @p all (mutable overload).
     * @param all Output vector; leaf item pointers are appended via back_inserter.
     */
    void item_iterator(std::vector<std::variant<Us...>*>& all)
    {
        if (m_left) {
            m_left->item_iterator(all);
            m_right->item_iterator(all);
        } else {
            std::copy(m_items.cbegin(), m_items.cend(), std::back_inserter(all));
        }
    }

    /**
     * @brief Recursively expands @p aabb to encompass all leaf items in the subtree.
     * @param aabb Running bounding box updated in place; must be pre-initialized to
     *             {max, max, max, lowest, lowest, lowest}.
     */
    void AABB_iterator(std::array<double, 6>& aabb) const
    {
        if (!m_left) {
            for (const auto& item : m_items) {
                const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
                for (std::size_t i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], aabb_tri[i]);
                }
                for (std::size_t i = 3; i < 6; ++i) {
                    aabb[i] = std::max(aabb[i], aabb_tri[i]);
                }
            }
        } else {
            m_left->AABB_iterator(aabb);
            m_right->AABB_iterator(aabb);
        }
    }
    /// @brief Returns the heuristic epsilon used for plane-side comparisons (11 × machine epsilon).
    constexpr static double epsilon()
    {
        // Huristic epsilon for triangle intersections
        return 11 * std::numeric_limits<double>::epsilon();
    }

    /// @brief Returns true if @p a ≤ @p b within the heuristic epsilon tolerance.
    constexpr static bool lessOrEqual(double a, double b)
    {
        return a - b <= epsilon() * a;
    }

    /// @brief Returns true if @p a ≥ @p b within the heuristic epsilon tolerance.
    constexpr static bool greaterOrEqual(double a, double b)
    {
        return b - a <= epsilon() * a;
    }

private:
    std::uint_fast32_t m_D = 0;
    double m_plane = 0;
    std::vector<std::variant<Us...>*> m_items;
    std::unique_ptr<KDTree> m_left = nullptr;
    std::unique_ptr<KDTree> m_right = nullptr;
};
}
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

// DEPRECATED

#pragma once

#include "dxmc/particle.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreetype.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <memory>
#include <optional>
#include <vector>

namespace dxmc {

template <MeshKDTreeType U>
class MeshKDTree {
public:
    MeshKDTree() {};
    MeshKDTree(MeshKDTree<U>& other) = delete;
    MeshKDTree(const MeshKDTree<U>& other) = delete;

    MeshKDTree(MeshKDTree<U>&& other) noexcept
    {
        m_D = other.m_D;
        m_plane = other.m_plane;
        m_trianglesIdx = std::move(other.m_trianglesIdx);
        if (other.m_left) {
            m_left = std::unique_ptr(std::move(other.m_left));
        }
        if (other.m_right) {
            m_right = std::unique_ptr(std::move(other.m_right));
        }
    }

    MeshKDTree(const std::vector<U>& triangles, const std::size_t max_depth = 8)
    {
        std::vector<std::uint32_t> idx(triangles.size());
        std::iota(idx.begin(), idx.end(), 0);
        setData(triangles, idx, max_depth);
    }

    MeshKDTree(const std::vector<U>& triangles, const std::vector<std::uint32_t>& idx, const std::size_t max_depth = 8)
    {
        setData(triangles, idx, max_depth);
    }

    void setData(const std::vector<U>& triangles, const std::size_t max_depth = 8)
    {
        std::vector<std::uint32_t> idx(triangles.size());
        std::iota(idx.begin(), idx.end(), 0);
        setData(triangles, idx, max_depth);
    }

    void setData(const std::vector<U>& triangles, const std::vector<std::uint32_t>& idx, const std::size_t max_depth = 8)
    {
        if (idx.size() == 0)
            return;
        // finding aabb
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };

        for (const auto tidx : idx) {
            // for (const auto& tri : triangles) {
            const auto& tri = triangles[tidx];
            const auto aabb_tri = tri.AABB();

            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<std::uint_fast32_t, double>(extent);

        const auto split = planeSplit(triangles, idx);

        const auto fom = figureOfMerit(triangles, idx, split);

        if (fom == idx.size() || max_depth <= 1 || idx.size() <= 1) {
            m_trianglesIdx = idx;
            m_left = nullptr;
            m_right = nullptr;
        } else {
            m_plane = split;
            std::vector<std::uint32_t> left;
            std::vector<std::uint32_t> right;
            for (const auto tidx : idx) {
                const auto& triangle = triangles[tidx];
                // for (const auto& triangle : triangles) {
                const auto side = planeSide(triangle, m_plane, m_D);
                if (side <= 0)
                    left.push_back(tidx);
                if (side >= 0)
                    right.push_back(tidx);
            }
            m_left = std::make_unique<MeshKDTree<U>>(triangles, left, max_depth - 1);
            m_right = std::make_unique<MeshKDTree<U>>(triangles, right, max_depth - 1);
        }
    }

    std::size_t depth() const
    {
        std::size_t teller = 0;
        depth_iterator(teller);
        return teller;
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_plane += dist[m_D];
        if (m_left) {
            m_left->translate(dist);
            m_right->translate(dist);
        }
    }

    void scale(double s)
    {
        m_plane *= s;
        if (m_left) {
            m_left->scale(s);
            m_right->scale(s);
        }
    }

    void mirror(const std::array<double, 3>& point)
    {
        m_plane += -2 * (m_plane - point[m_D]);
        if (m_left) {
            m_left->mirror(point);
            m_right->mirror(point);
        }
    }

    void mirror(const double value, const std::uint_fast32_t dim)
    {
        if (m_D == dim)
            m_plane += -2 * (m_plane - value);
        if (m_left) {
            m_left->mirror(value, dim);
            m_right->mirror(value, dim);
        }
    }

    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::vector<U>& triangles, const std::array<double, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, triangles, *inter) : KDTreeIntersectionResult<const U> {};
    }

protected:
    KDTreeIntersectionResult<const U> intersect(const ParticleType auto& particle, const std::vector<U>& triangles, const std::array<double, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;
            KDTreeIntersectionResult<const U> res;
            res.intersection = std::numeric_limits<double>::max();

            for (const auto triIdx : m_trianglesIdx) {
                // for (const U& triangle : m_triangles) {
                const auto& triangle = triangles[triIdx];
                const auto t_cand = triangle.intersect(particle);
                if (t_cand) {
                    if (0 <= *t_cand && *t_cand < res.intersection) {
                        if (tbox[0] <= *t_cand && *t_cand <= tbox[1]) {
                            res.intersection = *t_cand;
                            res.item = &triangle;
                            res.rayOriginIsInsideItem = vectormath::dot(particle.dir, triangle.planeVector()) > 0;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<double>::epsilon()) {
            auto hit_left = m_left->intersect(particle, triangles, tbox);
            auto hit_right = m_right->intersect(particle, triangles, tbox);
            if (hit_left.valid() && hit_right.valid())
                return hit_left.intersection > hit_right.intersection ? hit_right : hit_left;
            if (hit_right.valid())
                return hit_right;
            return hit_left;
        }

        const MeshKDTree<U>* const front = particle.dir[m_D] > 0 ? m_left.get() : m_right.get();
        const MeshKDTree<U>* const back = particle.dir[m_D] > 0 ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect(particle, triangles, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect(particle, triangles, tbox);
        }

        // both directions (start with front)
        const std::array<double, 2> t_front { tbox[0], t };
        auto hit = front->intersect(particle, triangles, t_front);
        if (hit.valid()) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<double, 2> t_back { t, tbox[1] };
        return back->intersect(particle, triangles, t_back);
    }

    double planeSplit(const std::vector<U>& triangles, const std::vector<std::uint32_t>& idx) const
    {
        const auto N = idx.size();
        std::vector<double> vals;
        vals.reserve(N);

        for (const auto triIdx : idx) {
            const auto& triangle = triangles[triIdx];
            const auto v = triangle.center()[m_D];
            vals.push_back(v);
        }
        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) * 0.5;
        }
    }

    int figureOfMerit(const std::vector<U>& triangles, const std::vector<std::uint32_t>& idx, const double planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto triIdx : idx) {
            const auto& triangle = triangles[triIdx];
            const auto side = planeSide(triangle, planesep, m_D);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

    static std::int_fast8_t planeSide(const U& triangle, const double plane, const std::uint_fast32_t D)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        const auto& aabb = triangle.AABB();
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

    void depth_iterator(std::size_t& teller) const
    {
        teller++;
        if (m_left)
            m_left->depth_iterator(teller);
    }

    constexpr static double epsilon()
    {
        // Huristic epsilon for triangle intersections
        return 11 * std::numeric_limits<double>::epsilon();
    }

    constexpr static bool lessOrEqual(double a, double b)
    {
        return a - b <= epsilon() * a;
    }
    constexpr static bool greaterOrEqual(double a, double b)
    {
        return b - a <= epsilon() * a;
    }

private:
    double m_plane = 0;
    std::unique_ptr<MeshKDTree<U>> m_left = nullptr;
    std::unique_ptr<MeshKDTree<U>> m_right = nullptr;
    std::vector<std::uint32_t> m_trianglesIdx;
    std::uint_fast32_t m_D = 0;
};

}
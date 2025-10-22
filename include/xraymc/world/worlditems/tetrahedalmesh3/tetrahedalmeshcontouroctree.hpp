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

Copyright 2025 Erlend Andersen
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

class TetrahedalMeshContourOcTree {
public:
    TetrahedalMeshContourOcTree() { };
    TetrahedalMeshContourOcTree(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, std::uint32_t max_depth = 8)
    {
        setData(vertices, elements, max_depth);
    }
    void setData(const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, std::uint32_t max_depth = 8)
    {
    }

    std::uint32_t maxDepth() const
    {
        return 0;
    }

    void translate(const std::array<double, 3>& dist)
    {
    }

    KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> intersect(const ParticleType auto& particle, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::array<double, 6>& aabb) const
    {
        auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, vertices, elements, *inter) : KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> {};
    }

    KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> intersect(const ParticleType auto& particle, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<std::uint32_t, 3>>& elements, const std::array<double, 2>& tboxAABB) const
    {
        return KDTreeIntersectionResult<const std::array<std::uint32_t, 3>> {};
    }

protected:
    inline auto childIndex(const std::array<double, 3>& pos, const std::array<double, 3>& nodeCenter)
    {
        return ((pos[0] > nodeCenter[0]) * 4) + ((pos[1] > nodeCenter[1]) * 2) + (pos[2] > nodeCenter[2]);
    }

private:
};
}
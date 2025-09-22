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

#include "xraymc/interactions.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/material/nistmaterials.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangle.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreeflat.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"

#include <array>
#include <execution>
#include <numeric>
#include <vector>

namespace xraymc {

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TriangulatedMesh {
public:
    TriangulatedMesh()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    TriangulatedMesh(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 4)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_triangles(triangles)
    {
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 4)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        m_triangles = reader(path);
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    TriangulatedMesh(const std::string& path, double scale, const std::size_t max_tree_dept = 4)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        m_triangles = reader(path);
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [=](auto& tri) {
            tri.scale(scale);
        });
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    void setMaterial(const Material<NMaterialShells>& mat)
    {
        m_material = mat;
    }

    void setMaterialDensity(double dens)
    {
        m_materialDensity = std::max(std::abs(dens), 0.00000001);
    }

    void setMaterial(const Material<NMaterialShells>& mat, double dens)
    {
        m_material = mat;
        setMaterialDensity(dens);
    }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials::density(nist_name);
            return true;
        }
        return false;
    }

    const MeshKDTreeFlat<Triangle>& kdtree() const
    {
        return m_kdtree;
    }

    const std::vector<Triangle>& triangles() const
    {
        return m_triangles;
    }

    void translate(const std::array<double, 3>& dist)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.translate(dist);
        });

        m_kdtree.translate(dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    void scale(double s)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.scale(s);
        });

        m_kdtree.scale(s);
        for (auto& e : m_aabb) {
            e *= s;
        }
    }

    void mirror(const std::array<double, 3>& point)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.mirror(point);
        });
        const auto depth = m_kdtree.maxDepth();
        m_kdtree.setData(m_triangles, depth);
        calculateAABB();
    }

    void mirror(const double value, const std::uint_fast32_t dim)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.mirror(value, dim);
        });
        const auto depth = m_kdtree.maxDepth();
        m_kdtree.setData(m_triangles, depth);
        calculateAABB();
    }

    void mirror(const std::uint_fast32_t dim)
    {
        const auto value = (m_aabb[dim] + m_aabb[dim + 3]) * 0.5;
        mirror(value, dim);
    }

    void rotate(const double angle, const std::array<double, 3>& axis)
    {
        constexpr std::array<double, 3> offset = { 0, 0, 0 };
        rotate(angle, axis, offset);
    }

    void rotate(const double angle, const std::array<double, 3>& axis, const std::array<double, 3>& point)
    {
        const auto depth = m_kdtree.maxDepth();

        const std::array point_neg = { -point[0], -point[1], -point[2] };

        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.translate(point_neg);
            tri.rotate(angle, axis);
            tri.translate(point);
        });
        m_kdtree.setData(m_triangles, depth);
        calculateAABB();
    }

    void setData(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 8)
    {
        m_triangles = triangles;
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    std::array<double, 3> center() const
    {
        std::array<double, 3> center { 0, 0, 0 };
        std::for_each(std::execution::unseq, m_triangles.cbegin(), m_triangles.cend(), [&](const auto& tri) {
            const auto c = tri.center();
            for (std::size_t i = 0; i < 3; ++i) {
                center[i] += c[i];
            }
        });
        for (std::size_t i = 0; i < 3; ++i) {
            center[i] /= m_triangles.size();
        }
        return center;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        return WorldIntersectionResult { .intersection = res.intersection, .rayOriginIsInsideItem = res.rayOriginIsInsideItem, .intersectionValid = res.item != nullptr };
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        VisualizationIntersectionResult<U> res_int;
        if (res.valid()) {
            res_int.normal = res.item->planeVector();
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            if (res_int.rayOriginIsInsideItem)
                res_int.normal = vectormath::scale(res_int.normal, -1.0);
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        bool cont = true;
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            const auto intM = intersect(p);
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
            if (intM.intersection < stepLen) {
                // transport to edge
                p.border_translate(intM.intersection);
                cont = false;
            } else {
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
                if (intRes.particleEnergyChanged) {
                    m_energyScored.scoreEnergy(intRes.energyImparted);
                }
            }
        }
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        // not defined for fleunce counter
        const auto vol = volume();
        m_dose.addScoredEnergy(m_energyScored, vol, m_materialDensity, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    void clearDoseScored()
    {
        m_dose.clear();
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    [[nodiscard]] double volume() const
    {
        // using signed volume of a thetrahedron by a point and the three vertices (Gauss theorem of divergence)
        // sign_vol = v1 * (v2 x v3) / 6 for three vectors defining a triangle, i.e triple product
        // volume is then the sum of signed volumes over all triangles

        // finding a neat reference point instead of {0,0,0}
        const std::array center = {
            (m_aabb[0] + m_aabb[3]) / 2,
            (m_aabb[1] + m_aabb[4]) / 2,
            (m_aabb[2] + m_aabb[5]) / 2
        };

        const auto vol = std::transform_reduce(std::execution::par_unseq, m_triangles.cbegin(), m_triangles.cend(), 0.0, std::plus<>(), [&center](const auto& tri) {
            const auto& v = tri.vertices();
            const auto v1 = vectormath::subtract(v[0], center);
            const auto v2 = vectormath::subtract(v[1], center);
            const auto v3 = vectormath::subtract(v[2], center);
            return vectormath::tripleProduct(v1, v2, v3);
        });
        return std::abs(vol) / 6;
    }

protected:
    void calculateAABB()
    {
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };

        for (const auto& tri : m_triangles) {
            const auto aabb_tri = tri.AABB();

            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }

        // Extending AABB
        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] -= GEOMETRIC_ERROR();
            aabb[i + 3] += GEOMETRIC_ERROR();
        }
        m_aabb = aabb;
    }

private:
    double m_materialDensity = 1;
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore m_energyScored;
    DoseScore m_dose;
    std::vector<Triangle> m_triangles;
    MeshKDTreeFlat<Triangle> m_kdtree;
    Material<NMaterialShells> m_material;
};
}
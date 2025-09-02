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

Copyright 2023 Erlend Andersen
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
class TriangulatedOpenSurface {
public:
    TriangulatedOpenSurface()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    TriangulatedOpenSurface(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    TriangulatedOpenSurface(const std::string& path, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        auto triangles = reader(path);
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    double surfaceThickness() const
    {
        return m_thickness;
    }
    void setSurfaceThickness(double cm)
    {
        m_thickness = std::abs(cm);
    }

    void setMaterial(const Material<NMaterialShells>& mat)
    {
        m_material = mat;
    }

    void setMaterialDensity(double dens)
    {
        m_materialDensity = std::abs(dens);
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

    const std::vector<Triangle>& getTriangles() const
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
        m_kdtree.mirror(point);
        calculateAABB();
    }

    void mirror(const double value, const std::uint_fast32_t dim)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.mirror(value, dim);
        });
        m_kdtree.mirror(value, dim);
        calculateAABB();
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

    void setData(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
    {

        m_thickness = std::abs(surfaceThickness);
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
        // Not to be used for internal intersections
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        WorldIntersectionResult wres;
        if (res.valid()) {
            const auto normal = res.item->planeVector();
            const auto length_scale = std::abs(vectormath::dot(p.dir, normal));
            const auto length = m_thickness / (2 * length_scale);

            wres.intersection = res.intersection - length;
            wres.rayOriginIsInsideItem = false;
            wres.intersectionValid = true;
        }
        return wres;
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        VisualizationIntersectionResult<U> res_int;
        if (res.valid()) {
            res_int.normal = res.item->planeVector();
            const auto length_scale = vectormath::dot(p.dir, res_int.normal);
            if (length_scale > 0.0) // fix normal vector
                res_int.normal = vectormath::scale(res_int.normal, -1.0);
            const auto length = m_thickness / (2 * std::abs(length_scale));
            res_int.intersection = res.intersection - length;
            res_int.rayOriginIsInsideItem = false;
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        const auto intersection = m_kdtree.intersect(p, m_triangles, m_aabb);

        const auto normal = intersection.item->planeVector();

        const auto plane_point_plus = vectormath::add(intersection.item->vertices()[0], vectormath::scale(normal, m_thickness * 0.5));
        const auto plane_point_minus = vectormath::add(intersection.item->vertices()[0], vectormath::scale(normal, -m_thickness * 0.5));

        bool inside;
        do {
            const auto nd = vectormath::dot(p.dir, normal);
            if (nd > GEOMETRIC_ERROR<double>()) {
                const auto pv = vectormath::subtract(p.pos, plane_point_plus);
                const auto lenght = -vectormath::dot(pv, normal) / nd;

                const auto att = m_material.attenuationValues(p.energy);
                const auto attSumInv = 1.0 / (att.sum() * m_materialDensity);
                const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
                if (stepLen > lenght) {
                    p.border_translate(lenght);
                    inside = false;
                } else {
                    // interact
                    p.translate(stepLen);
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    if (intRes.particleEnergyChanged) {
                        m_energyScored.scoreEnergy(intRes.energyImparted);
                    }
                    inside = intRes.particleAlive;
                }
            } else if (nd < GEOMETRIC_ERROR<double>()) {
                const auto pv = vectormath::subtract(p.pos, plane_point_minus);
                const auto lenght = -vectormath::dot(pv, normal) / nd;

                const auto att = m_material.attenuationValues(p.energy);
                const auto attSumInv = 1.0 / (att.sum() * m_materialDensity);
                const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
                if (stepLen > lenght) {
                    p.border_translate(lenght);
                    inside = false;
                } else {
                    // interact
                    p.translate(stepLen);
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    if (intRes.particleEnergyChanged) {
                        m_energyScored.scoreEnergy(intRes.energyImparted);
                    }
                    inside = intRes.particleAlive;
                }
            } else {
                // we do an interaction if in plane
                const auto att = m_material.attenuationValues(p.energy);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                if (intRes.particleEnergyChanged) {
                    m_energyScored.scoreEnergy(intRes.energyImparted);
                }
                inside = intRes.particleAlive;
            }
        } while (inside);
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
        const auto volume = calculateVolume(m_triangles, m_thickness);
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
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

    /*std::optional<double> intersectOffsetTriangle(const ParticleType auto& p const Triangle* tri)
    {
        Triangle off = *tri;
        const auto dist = vectormath::scale(off.planeVector(), m_thickness);
        off.translate(dist);
        return off.intersect(p)
    }

    std::optional<double> intersectPrismPlanes(const ParticleType auto& p, const Triangle* tri)
    {
        const auto& verts = tri->vertices();
        const auto v_norm = tri->planeVector();

        double min_inter = -1;
        // first plane
        {
            const auto n = vectormath::cross(v_norm, vectormath::subtract(verts[0], verts[1]));
            const auto dn = vectormath::dot(p.dir, n);
            if (std::abs(dn) > GEOMETRIC_ERROR<double>()) {
                const auto p = vectormath::subtract(p.pos, verts[0]);
                min_inter = vectormath::dot(p, b) / dn;
            }
        }
        // second plane
        {
            const auto n = vectormath::cross(v_norm, vectormath::subtract(verts[1], verts[2]));
            const auto dn = vectormath::dot(p.dir, n);
            if (std::abs(dn) > GEOMETRIC_ERROR<double>()) {
                const auto p = vectormath::subtract(p.pos, verts[1]);
                const auto cand = vectormath::dot(p, b) / dn;
                if (cand > 0) {
                    if (min_inter > 0)
                        min_inter = std::min(cand, min_inter);
                    else
                        min_inter = cand;
                }
            }
        }
        // third plane
        {
            const auto n = vectormath::cross(v_norm, vectormath::subtract(verts[2], verts[0]));
            const auto dn = vectormath::dot(p.dir, n);
            if (std::abs(dn) > GEOMETRIC_ERROR<double>()) {
                const auto p = vectormath::subtract(p.pos, verts[2]);
                const auto cand = vectormath::dot(p, b) / dn;
                if (cand > 0) {
                    if (min_inter > 0)
                        min_inter = std::min(cand, min_inter);
                    else
                        min_inter = cand;
                }
            }
        }
        return min_inter < 0 ? std::nullopt : min_inter;
    }*/

    static double calculateVolume(const std::vector<Triangle>& triangles, double thickness)
    {
        const auto area = std::transform_reduce(std::execution::par_unseq, triangles.cbegin(), triangles.cend(), 0.0, std::plus<>(), [](const auto& tri) {
            return tri.area();
        });
        return area * thickness;
    }

private:
    double m_materialDensity = 1;
    double m_thickness = 0.035; // cm
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore m_energyScored;
    DoseScore m_dose;
    MeshKDTreeFlat<Triangle> m_kdtree;
    std::vector<Triangle> m_triangles;
    Material<NMaterialShells> m_material;
};
}
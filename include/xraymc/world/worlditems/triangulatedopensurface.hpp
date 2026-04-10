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
#include "xraymc/serializer.hpp"
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

/**
 * @brief A triangulated open-surface geometry that models a thin shell of material.
 *
 * Unlike a closed mesh, this class does not enclose a volume. Instead, each triangle
 * represents a patch of a surface with a configurable finite thickness. Particles
 * entering a triangle face travel through a slab of the given thickness perpendicular
 * to the triangle plane before exiting the other side. Ray–triangle intersections are
 * accelerated by a flat KD-tree (MeshKDTreeFlat). Absorbed energy is accumulated in
 * an EnergyScore and can be converted to dose via addEnergyScoredToDoseScore(), using
 * the total surface area times thickness as the volume estimate.
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 */
template <int NMaterialShells = 16, int LOWENERGYCORRECTION = 2>
class TriangulatedOpenSurface {
public:
    /**
     * @brief Default constructor. Initializes the surface with no triangles and
     *        dry air material at standard sea-level density.
     */
    TriangulatedOpenSurface()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    /**
     * @brief Constructs an open surface from an existing triangle list.
     * @param triangles        Collection of triangles defining the surface patches.
     * @param surfaceThickness Thickness of the material slab at each triangle in cm.
     * @param max_tree_dept    Maximum depth of the KD-tree used for intersection acceleration.
     */
    TriangulatedOpenSurface(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    /**
     * @brief Constructs an open surface by reading triangles from an STL file.
     * @param path             Path to the STL file.
     * @param surfaceThickness Thickness of the material slab at each triangle in cm.
     * @param max_tree_dept    Maximum depth of the KD-tree used for intersection acceleration.
     */
    TriangulatedOpenSurface(const std::string& path, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        auto triangles = reader(path);
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    /// @brief Returns the surface slab thickness in cm.
    double surfaceThickness() const
    {
        return m_thickness;
    }

    /**
     * @brief Sets the surface slab thickness.
     * @param cm New thickness in cm; absolute value is used.
     */
    void setSurfaceThickness(double cm)
    {
        m_thickness = std::abs(cm);
    }

    /**
     * @brief Sets the material used for attenuation and interaction sampling.
     * @param mat Material cross-section data.
     */
    void setMaterial(const Material<NMaterialShells>& mat)
    {
        m_material = mat;
    }

    /**
     * @brief Sets the material mass density.
     * @param dens Density in g/cm³; absolute value is used.
     */
    void setMaterialDensity(double dens)
    {
        m_materialDensity = std::abs(dens);
    }

    /**
     * @brief Sets both the material and its mass density in one call.
     * @param mat  Material cross-section data.
     * @param dens Density in g/cm³; absolute value is used.
     */
    void setMaterial(const Material<NMaterialShells>& mat, double dens)
    {
        m_material = mat;
        setMaterialDensity(dens);
    }

    /**
     * @brief Sets the material and density from the NIST material database by name.
     * @param nist_name Name of the NIST material (e.g., "Water, Liquid").
     * @return true if the name was found and the material was set; false otherwise.
     */
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

    /// @brief Returns the flat KD-tree used to accelerate ray–triangle intersection tests.
    const MeshKDTreeFlat<Triangle>& kdtree() const
    {
        return m_kdtree;
    }

    /// @brief Returns the list of triangles defining the open surface.
    const std::vector<Triangle>& getTriangles() const
    {
        return m_triangles;
    }

    /**
     * @brief Translates all triangles, the KD-tree, and the AABB by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
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

    /**
     * @brief Uniformly scales all triangle vertices, the KD-tree, and the AABB by factor @p s.
     * @param s Scale factor applied to all vertex coordinates.
     */
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

    /**
     * @brief Rotates all triangles about @p point by @p angle radians around @p axis,
     *        then rebuilds the KD-tree and AABB.
     * @param angle Rotation angle in radians.
     * @param axis  Rotation axis (need not be normalized).
     * @param point World-space pivot point for the rotation.
     */
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

    /**
     * @brief Replaces the current triangle data, sets the surface thickness, rebuilds the
     *        KD-tree, and recalculates the AABB.
     * @param triangles        New collection of triangles defining the surface.
     * @param surfaceThickness Thickness of the material slab at each triangle in cm.
     * @param max_tree_dept    Maximum depth of the KD-tree.
     */
    void setData(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
    {

        m_thickness = std::abs(surfaceThickness);
        m_triangles = triangles;
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    /**
     * @brief Returns the centroid of the surface, computed as the mean of all triangle centers.
     * @return World-space centroid in cm.
     */
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

    /**
     * @brief Tests a particle ray against the surface using the KD-tree accelerator.
     *
     * The returned intersection distance is adjusted inward by half the slab thickness
     * (projected along the ray) so the particle enters the material before reaching the
     * triangle plane. Not intended for particles originating inside the surface.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the adjusted entry distance, always with
     *         rayOriginIsInsideItem = false.
     */
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

    /**
     * @brief Like intersect(), but also returns the triangle surface normal (oriented
     *        toward the incoming ray) and the accumulated dose value for visualization.
     * @tparam U Scalar type used for the value field (set to the current mean dose).
     * @param p  Particle whose position and direction define the ray.
     */
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

    /**
     * @brief Transports a particle through the surface slab using analog Monte Carlo sampling.
     *
     * The particle travels through a slab of thickness m_thickness centered on the
     * intersected triangle plane. Free-path sampling and photon interactions continue
     * until the particle exits the slab or is absorbed. The exit boundary is determined
     * by projecting the slab half-thickness along the particle direction. Imparted energy
     * is recorded in the EnergyScore accumulator.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Returns the energy-score accumulator for the surface.
     * @param index Unused; present for interface consistency.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    /// @brief Resets the energy-score accumulator to zero.
    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    /**
     * @brief Converts the accumulated energy score to a dose score using the surface
     *        area times thickness as the effective volume.
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = calculateVolume(m_triangles, m_thickness);
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the surface.
     * @param index Unused; present for interface consistency.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    /// @brief Resets the dose-score accumulator to zero.
    void clearDoseScored()
    {
        m_dose.clear();
    }

    /// @brief Returns the axis-aligned bounding box of the surface as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "TriOpenSurface1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells);
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether a raw data buffer begins with the expected magic identifier.
     * @param data Buffer to inspect; must be at least 32 bytes for a positive result.
     * @return true if the first 32 bytes match magicID().
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the surface (triangle vertex data, material density, surface thickness,
     *        material composition, and dose score) to a byte vector that can be restored via
     *        deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        std::vector<double> tris_points;
        tris_points.reserve(m_triangles.size() * 3 * 3); // three vertices and three coords per vertice
        for (const auto& tri : m_triangles) {
            for (const auto& v : tri.vertices()) {
                for (const auto& p : v) {
                    tris_points.push_back(p);
                }
            }
        }
        Serializer::serialize(tris_points, buffer);

        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serialize(m_thickness, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a surface from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed surface, or std::nullopt if the data is malformed (triangle count
     *         not a multiple of 9, or material composition cannot be resolved).
     */
    static std::optional<TriangulatedOpenSurface<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
    {
        std::vector<double> points;
        buffer = Serializer::deserialize(points, buffer);

        // Test if number of points is making triangles
        if (points.size() % 9 != 0)
            return std::nullopt;

        std::vector<Triangle> triangles;
        triangles.reserve(points.size() / (3 * 3));
        for (std::size_t i = 0; i < points.size(); i = i + 9) {
            std::array<std::array<double, 3>, 3> data;
            std::size_t k = i;
            for (auto& v : data) {
                for (auto& p : v) {
                    p = points[k++];
                }
            }
            triangles.emplace_back(Triangle { data });
        }

        TriangulatedOpenSurface<NMaterialShells, LOWENERGYCORRECTION> item(triangles);

        buffer = Serializer::deserialize(item.m_materialDensity, buffer);
        buffer = Serializer::deserialize(item.m_thickness, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        auto material_opt = Material<NMaterialShells>::byWeight(mat_weights);
        if (material_opt) {
            item.m_material = material_opt.value();
        } else {
            return std::nullopt;
        }

        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        return item;
    }

protected:
    /// @brief Recomputes the AABB by iterating over all triangle AABBs and extending by GEOMETRIC_ERROR().
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

    /**
     * @brief Estimates the effective volume as total surface area times @p thickness.
     * @param triangles Collection of triangles whose areas are summed.
     * @param thickness Surface slab thickness in cm.
     * @return Effective volume in cm³ used for dose normalization.
     */
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
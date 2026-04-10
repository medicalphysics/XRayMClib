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
 * @brief A triangulated surface mesh geometry for Monte Carlo particle transport.
 *
 * The mesh is defined by a set of triangles (typically loaded from an STL file) that
 * enclose a homogeneous volume with a single material and density. Ray–triangle
 * intersections are accelerated by a flat KD-tree (MeshKDTreeFlat). Particles are
 * tracked inside the enclosed volume until they exit through a triangle face.
 * Absorbed energy is accumulated in an EnergyScore and can be converted to dose
 * via addEnergyScoredToDoseScore().
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 */
template <int NMaterialShells = 16, int LOWENERGYCORRECTION = 2>
class TriangulatedMesh {
public:
    /**
     * @brief Default constructor. Initializes the mesh with no triangles and
     *        dry air material at standard sea-level density.
     */
    TriangulatedMesh()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    /**
     * @brief Constructs a mesh from an existing triangle list.
     * @param triangles     Collection of triangles defining the surface.
     * @param max_tree_dept Maximum depth of the KD-tree used for intersection acceleration.
     */
    TriangulatedMesh(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 4)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_triangles(triangles)
    {
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    /**
     * @brief Constructs a mesh by reading triangles from an STL file.
     * @param path          Path to the STL file.
     * @param max_tree_dept Maximum depth of the KD-tree used for intersection acceleration.
     */
    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 4)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        m_triangles = reader(path);
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    /**
     * @brief Constructs a mesh by reading triangles from an STL file and uniformly scaling them.
     * @param path          Path to the STL file.
     * @param scale         Uniform scale factor applied to all vertex coordinates.
     * @param max_tree_dept Maximum depth of the KD-tree used for intersection acceleration.
     */
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
     * @param dens Density in g/cm³; absolute value is used, clamped to at least 1e-8 g/cm³.
     */
    void setMaterialDensity(double dens)
    {
        m_materialDensity = std::max(std::abs(dens), 0.00000001);
    }

    /**
     * @brief Sets both the material and its mass density in one call.
     * @param mat  Material cross-section data.
     * @param dens Density in g/cm³; clamped to at least 1e-8 g/cm³.
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

    /// @brief Returns the list of triangles defining the mesh surface.
    const std::vector<Triangle>& triangles() const
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
     * @brief Mirrors the mesh through a given world-space point and rebuilds the KD-tree.
     * @param point The point through which each triangle vertex is reflected.
     */
    void mirror(const std::array<double, 3>& point)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.mirror(point);
        });
        const auto depth = m_kdtree.maxDepth();
        m_kdtree.setData(m_triangles, depth);
        calculateAABB();
    }

    /**
     * @brief Mirrors the mesh about a plane perpendicular to axis @p dim at position @p value,
     *        then rebuilds the KD-tree and AABB.
     * @param value Coordinate of the mirror plane along @p dim.
     * @param dim   Axis index: 0 = x, 1 = y, 2 = z.
     */
    void mirror(const double value, const std::uint_fast32_t dim)
    {
        std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
            tri.mirror(value, dim);
        });
        const auto depth = m_kdtree.maxDepth();
        m_kdtree.setData(m_triangles, depth);
        calculateAABB();
    }

    /**
     * @brief Mirrors the mesh about the mid-plane of the AABB along axis @p dim.
     * @param dim Axis index: 0 = x, 1 = y, 2 = z.
     */
    void mirror(const std::uint_fast32_t dim)
    {
        const auto value = (m_aabb[dim] + m_aabb[dim + 3]) * 0.5;
        mirror(value, dim);
    }

    /**
     * @brief Rotates all triangles about the world origin by @p angle radians around @p axis,
     *        then rebuilds the KD-tree and AABB.
     * @param angle Rotation angle in radians.
     * @param axis  Rotation axis (need not be normalized).
     */
    void rotate(const double angle, const std::array<double, 3>& axis)
    {
        constexpr std::array<double, 3> offset = { 0, 0, 0 };
        rotate(angle, axis, offset);
    }

    /**
     * @brief Rotates all triangles about an arbitrary @p point by @p angle radians around @p axis,
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
     * @brief Replaces the current triangle data, rebuilds the KD-tree, and recalculates the AABB.
     * @param triangles     New collection of triangles defining the surface.
     * @param max_tree_dept Maximum depth of the KD-tree.
     */
    void setData(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 8)
    {
        m_triangles = triangles;
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    /**
     * @brief Returns the centroid of the mesh, computed as the mean of all triangle centers.
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
     * @brief Tests a particle ray against the mesh surface using the KD-tree accelerator.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with distance to the nearest triangle face, including
     *         whether the ray origin is inside the mesh.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        return WorldIntersectionResult { .intersection = res.intersection, .rayOriginIsInsideItem = res.rayOriginIsInsideItem, .intersectionValid = res.item != nullptr };
    }

    /**
     * @brief Like intersect(), but also returns the triangle surface normal and the
     *        accumulated dose value for visualization.
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
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            if (res_int.rayOriginIsInsideItem)
                res_int.normal = vectormath::scale(res_int.normal, -1.0);
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    /**
     * @brief Transports a particle through the mesh using analog Monte Carlo sampling.
     *
     * The particle undergoes repeated free-path sampling and photon interactions until
     * it exits the mesh surface or is absorbed. Imparted energy is recorded in the
     * EnergyScore accumulator.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Returns the energy-score accumulator for the mesh.
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
     * @brief Converts the accumulated energy score to a dose score using the mesh volume and density.
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        // not defined for fleunce counter
        const auto vol = volume();
        m_dose.addScoredEnergy(m_energyScored, vol, m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the mesh.
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

    /// @brief Returns the axis-aligned bounding box of the mesh as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Computes the enclosed volume of the mesh using the divergence theorem
     *        (signed tetrahedron volumes summed over all triangles).
     * @return Enclosed volume in cm³.
     */
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

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "TriMesh1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells);
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
     * @brief Serializes the mesh (triangle vertex data, material density, material composition,
     *        and dose score) to a byte vector that can be restored via deserialize().
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
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a mesh from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed mesh, or std::nullopt if the data is malformed (triangle count
     *         not a multiple of 9, or material composition cannot be resolved).
     */
    static std::optional<TriangulatedMesh<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
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

        TriangulatedMesh<NMaterialShells, LOWENERGYCORRECTION> item(triangles);

        buffer = Serializer::deserialize(item.m_materialDensity, buffer);

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
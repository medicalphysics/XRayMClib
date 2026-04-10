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
#include "xraymc/particletracker.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcollection.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcontourkdtree.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <optional>

namespace xraymc {

/**
 * @brief An unstructured tetrahedral mesh geometry for Monte Carlo particle transport.
 *
 * The mesh is partitioned into named collections, each with its own material and
 * density.  Particles are tracked cell-by-cell using a Siddon-style cumulative step
 * algorithm for numerical stability.  The outer surface is accelerated by a KD-tree
 * built from boundary triangles.  Per-tetrahedron energy and dose accumulators are
 * provided along with volume-weighted collection summaries.
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 * @tparam FORCEDINTERACTION   When true, forces photoelectric interactions at each
 *                             tetrahedron boundary for variance-reduction.
 */
template <int NMaterialShells = 16, int LOWENERGYCORRECTION = 2, bool FORCEDINTERACTION = false>
class TetrahedalMesh {
    // static constexpr auto EPSILON = std::numeric_limits<double>::epsilon() * 1000;
    static constexpr double EPSILON = 1E-10;
    // static constexpr double EPSILON = GEOMETRIC_ERROR();

public:
    /// @brief Default constructor. Creates an empty, invalid mesh.
    TetrahedalMesh() { };

    /**
     * @brief Constructs a mesh and immediately loads geometry and material data.
     * @param data Validated mesh data; see setData() for details.
     */
    TetrahedalMesh(const TetrahedalMeshData& data)
    {
        setData(data);
    }

    /**
     * @brief Loads mesh geometry and collection data, building all internal structures.
     * @param data Mesh data containing nodes, elements, collection indices, densities,
     *             material compositions, and collection names.
     * @return true on success; false if @p data is invalid or any material cannot be
     *         constructed from its elemental weight map.
     */
    bool setData(const TetrahedalMeshData& data)
    {
        if (!data.valid())
            return false;
        makeStructure(data);
        calculateAABB();
        buildContourKDTree();

        m_doseScore.resize(m_tetrahedrons.size());
        m_energyScore.resize(m_tetrahedrons.size());

        m_collectionIdx = data.collectionIndices;
        m_collectionDensities = data.collectionDensities;
        m_collectionNames = data.collectionNames;

        std::vector<std::optional<Material<NMaterialShells>>> mats(data.collectionMaterialComposition.size());
        std::transform(std::execution::par_unseq, data.collectionMaterialComposition.cbegin(), data.collectionMaterialComposition.cend(), mats.begin(), [](const auto& w) {
            return Material<NMaterialShells>::byWeight(w);
        });
        m_collectionMaterials.clear();
        m_collectionMaterials.reserve(data.collectionMaterialComposition.size());
        for (auto& m : mats) {
            if (m)
                m_collectionMaterials.push_back(m.value());
            else
                return false;
        }

        return true;
    }

    /**
     * @brief Translates all vertices, the AABB, and the contour KD-tree by @p vec.
     * @param vec Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& vec)
    {
        std::transform(std::execution::par_unseq, m_vertices.begin(), m_vertices.end(), m_vertices.begin(), [=](const auto& v) {
            return vectormath::add(v, vec);
        });
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = m_aabb[i] + vec[i];
            m_aabb[i + 3] = m_aabb[i + 3] + vec[i];
        }
        m_kdtree.translate(vec);
    }

    /**
     * @brief Rotates all vertices about @p axis (through the mesh center) by @p angle,
     *        then rebuilds the AABB and contour KD-tree.
     * @param axis  Rotation axis (need not be normalized).
     * @param angle Rotation angle in radians.
     */
    void rotate(const std::array<double, 3>& axis, double angle)
    {
        const std::array center = { (m_aabb[3] + m_aabb[0]) / 2, (m_aabb[4] + m_aabb[1]) / 2, (m_aabb[5] + m_aabb[2]) / 2 };
        std::transform(std::execution::par_unseq, m_vertices.begin(), m_vertices.end(), m_vertices.begin(), [&](const auto& v) {
            return vectormath::add(vectormath::rotate(vectormath::subtract(v, center), axis, angle), center);
        });
        calculateAABB();
        m_kdtree.setData(m_vertices, m_outer_triangles, 8);
    }

    /**
     * @brief Restricts visualization to tetrahedrons belonging to the named collection.
     * @param collectionName Name to match against the collection list.
     * @return true if the name was found and the filter was set; false otherwise.
     */
    bool setDisplayCollectionIndexFilter(const std::string& collectionName)
    {
        std::uint32_t teller = 0;
        bool found = false;
        while (teller < m_collectionNames.size() && !found) {
            found = collectionName.compare(m_collectionNames[teller]) == 0;
        }
        return setDisplayCollectionIndexFilter(teller);
    }

    /**
     * @brief Restricts visualization to tetrahedrons belonging to collection @p idx.
     *        Pass max uint32 (or an out-of-range index) to show all collections.
     * @param idx Zero-based collection index.
     * @return true if @p idx is valid; false otherwise (filter is cleared).
     */
    bool setDisplayCollectionIndexFilter(std::uint32_t idx)
    {
        if (idx < m_collectionMaterials.size()) {
            m_collectionIndexForDisplay = idx;
            return true;
        }
        m_collectionIndexForDisplay = std::numeric_limits<std::uint32_t>::max();
        return false;
    }

    /// @brief Returns the axis-aligned bounding box of the mesh as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /// @brief Returns the geometric center of the AABB in world coordinates.
    std::array<double, 3> center() const
    {
        std::array<double, 3> res = {
            (m_aabb[0] + m_aabb[3]) * 0.5,
            (m_aabb[1] + m_aabb[4]) * 0.5,
            (m_aabb[2] + m_aabb[5]) * 0.5
        };
        return res;
    }

    /**
     * @brief Returns the volume of a single tetrahedron in cm³.
     * @param index Zero-based tetrahedron index.
     */
    double tetrahedalVolume(std::uint32_t index) const
    {
        const auto& tet = m_tetrahedrons.at(index).verticeIdx;
        const auto& v0 = m_vertices[tet[0]];
        const auto& v1 = m_vertices[tet[1]];
        const auto& v2 = m_vertices[tet[2]];
        const auto& v3 = m_vertices[tet[3]];
        return basicshape::tetrahedron::volume(v0, v1, v2, v3);
    }

    /// @brief Returns the total number of tetrahedrons in the mesh.
    std::uint32_t numberOfThetrahedrons() const
    {
        return static_cast<std::uint32_t>(m_tetrahedrons.size());
    }

    /// @brief Returns the names of all material collections in the mesh.
    const std::vector<std::string> collectionNames() const
    {
        return m_collectionNames;
    }

    /**
     * @brief Looks up the collection index for a given collection name.
     * @param name Collection name to search for.
     * @return The zero-based index on success, or std::nullopt if not found.
     */
    std::optional<std::uint32_t> tetrahedronCollectionIndex(const std::string& name) const
    {
        std::uint32_t teller = 0;
        while (teller < m_collectionNames.size()) {
            if (name.compare(m_collectionNames[teller]) == 0)
                return teller;
            ++teller;
        }
        return std::nullopt;
    }

    /**
     * @brief Returns the collection index for a given tetrahedron.
     * @param idx Zero-based tetrahedron index.
     */
    std::uint32_t tetrahedronCollectionIndex(std::uint32_t idx) const
    {
        return m_collectionIdx.at(idx);
    }

    /**
     * @brief Returns the four vertex positions of a tetrahedron.
     * @param ind Zero-based tetrahedron index.
     * @return Array of four {x, y, z} vertex positions in cm.
     */
    std::array<std::array<double, 3>, 4> tetrahedron(std::uint32_t ind) const
    {
        const auto& t = m_tetrahedrons.at(ind);
        std::array<std::array<double, 3>, 4> tet;
        for (std::uint32_t i = 0; i < 4; ++i) {
            tet[i] = m_vertices[t.verticeIdx[i]];
        }
        return tet;
    }

    /// @brief Returns the tetrahedron index for each outer-contour triangle (parallel to the triangle list).
    const std::vector<std::uint32_t>& outerContourTetrahedronIndices() const
    {
        return m_outerTriangleTetMembership;
    }

    /// @brief Returns a read-only reference to the particle tracker.
    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    /// @brief Returns a mutable reference to the particle tracker.
    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

    /**
     * @brief Aggregates per-tetrahedron dose scores into volume-weighted collection summaries.
     * @return One TetrahedalMeshCollection entry per collection, containing total volume,
     *         volume-weighted mean dose, dose variance, event count, density, and name.
     */
    std::vector<TetrahedalMeshCollection> collectionData() const
    {
        std::vector<TetrahedalMeshCollection> data(m_collectionDensities.size());
        for (std::uint32_t i = 0; i < m_tetrahedrons.size(); ++i) {
            const auto& tet = m_tetrahedrons[i];
            const auto collIdx = m_collectionIdx[i];
            const auto vol = tetrahedalVolume(i);
            auto& dat = data[collIdx];
            dat.volume += vol;
            dat.dose += m_doseScore[i].dose() * vol;
            dat.doseVariance += m_doseScore[i].variance() * vol * vol;
            dat.numberOfEvents += m_doseScore[i].numberOfEvents();
        }

        for (auto& d : data) {
            d.dose /= d.volume;
            d.doseVariance /= d.volume * d.volume;
        }
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i].density = m_collectionDensities[i];
            data[i].name = m_collectionNames[i];
        }
        return data;
    }

    /**
     * @brief Tests a particle ray against the mesh outer surface (AABB pre-filter + KD-tree).
     * @param particle Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the outer surface, or invalid if missed.
     */
    WorldIntersectionResult intersect(const ParticleType auto& particle) const
    {
        WorldIntersectionResult res { };
        if (auto inter = basicshape::AABB::intersectForwardInterval(particle, m_aabb); inter) {
            auto kdres = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, *inter);
            if (kdres.valid()) {
                res.intersection = kdres.intersection;
                res.rayOriginIsInsideItem = kdres.rayOriginIsInsideItem;
                res.intersectionValid = true;
            }
        };
        return res;
    }

    /**
     * @brief Like intersect(), but also returns the surface normal and per-tetrahedron dose
     *        for rendering.  When a display collection filter is active, the ray walks
     *        forward through the mesh until a tetrahedron of the selected collection is hit.
     * @tparam U Scalar type used for the dose value in the visualization result.
     * @param particle Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& particle) const
    {
        VisualizationIntersectionResult<U> res { };
        if (auto inter = basicshape::AABB::intersectForwardInterval(particle, m_aabb); inter) {
            auto kdres = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, *inter);
            if (kdres.valid()) {
                if (m_collectionIndexForDisplay == std::numeric_limits<std::uint32_t>::max()) {
                    res.intersection = kdres.intersection;
                    res.rayOriginIsInsideItem = kdres.rayOriginIsInsideItem;
                    res.intersectionValid = true;

                    // normal
                    const auto& vIdx = *(kdres.item);
                    res.normal = basicshape::tetrahedron::normalVector(m_vertices[vIdx[0]], m_vertices[vIdx[1]], m_vertices[vIdx[2]]);
                    // value
                    const auto faceIdx = std::distance(&(m_outer_triangles[0]), kdres.item);
                    const auto tetIdx = m_outerTriangleTetMembership[faceIdx];
                    res.value = m_doseScore[tetIdx].dose();
                } else {
                    const auto faceIdx = std::distance(&(m_outer_triangles[0]), kdres.item);
                    auto currentTetIdx = m_outerTriangleTetMembership[faceIdx];

                    double intersection = kdres.intersection;
                    bool cont = true;
                    do {
                        if (m_collectionIdx[currentTetIdx] == m_collectionIndexForDisplay) {
                            cont = false;
                        } else {
                            const auto [lenght, nextTetIdx] = nextTetrahedron(currentTetIdx, particle);
                            if (nextTetIdx == currentTetIdx) {
                                // we hit the end, are there others behind the bondary?
                                std::array<double, 2> interseg = { lenght + GEOMETRIC_ERROR<>() * 2, (*inter)[1] };
                                auto kdseg = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, interseg);
                                if (kdseg.valid()) {
                                    const auto segfaceIdx = std::distance(&(m_outer_triangles[0]), kdseg.item);
                                    currentTetIdx = m_outerTriangleTetMembership[segfaceIdx];
                                    intersection = kdseg.intersection;
                                } else {
                                    cont = false;
                                }
                            } else {
                                intersection = lenght;
                                currentTetIdx = nextTetIdx;
                            }
                        }

                    } while (cont);

                    if (m_collectionIdx[currentTetIdx] == m_collectionIndexForDisplay) {
                        res.intersection = intersection;
                        res.rayOriginIsInsideItem = kdres.rayOriginIsInsideItem;
                        res.intersectionValid = true;
                        const auto point = vectormath::add(particle.pos, vectormath::scale(particle.dir, res.intersection));
                        const auto& v = m_tetrahedrons[currentTetIdx].verticeIdx;
                        res.normal = basicshape::tetrahedron::closestNormalToPoint(m_vertices[v[0]], m_vertices[v[1]], m_vertices[v[2]], m_vertices[v[3]], point);
                        res.value = m_doseScore[currentTetIdx].dose();
                    }
                }
            }
        }
        return res;
    }

    /**
     * @brief Returns a copy of the energy-score accumulator for a single tetrahedron.
     * @param index Zero-based tetrahedron index; defaults to the first tetrahedron.
     */
    EnergyScore energyScored(std::size_t index = 0) const
    {
        return m_energyScore.at(index);
    }

    /**
     * @brief Returns a copy of the dose-score accumulator for a single tetrahedron.
     * @param index Zero-based tetrahedron index; defaults to the first tetrahedron.
     */
    DoseScore doseScored(std::size_t index = 0) const
    {
        return m_doseScore.at(index);
    }

    /// @brief Resets all per-tetrahedron dose-score accumulators to zero.
    void clearDoseScored()
    {
        std::for_each(std::execution::par_unseq, m_doseScore.begin(), m_doseScore.end(), [](auto& d) {
            d.clear();
        });
    }

    /// @brief Resets all per-tetrahedron energy-score accumulators to zero.
    void clearEnergyScored()
    {
        std::for_each(std::execution::par_unseq, m_energyScore.begin(), m_energyScore.end(), [](auto& d) {
            d.clear();
        });
    }

    /**
     * @brief Accumulates per-tetrahedron energy scores into dose scores using each
     *        tetrahedron's volume, its collection's density, and a calibration factor.
     * @param factor Optional scaling factor applied to each scored energy value.
     */
    void addEnergyScoredToDoseScore(double factor = 1)
    {
        std::vector<double> volumes(m_tetrahedrons.size());
        std::transform(std::execution::par_unseq, m_tetrahedrons.cbegin(), m_tetrahedrons.cend(), volumes.begin(), [&](const auto& tet) {
            const auto& vIdx = tet.verticeIdx;
            const auto& v0 = m_vertices[vIdx[0]];
            const auto& v1 = m_vertices[vIdx[1]];
            const auto& v2 = m_vertices[vIdx[2]];
            const auto& v3 = m_vertices[vIdx[3]];
            return basicshape::tetrahedron::volume(v0, v1, v2, v3);
        });

        std::vector<double> density(m_tetrahedrons.size());
        std::transform(std::execution::par_unseq, m_collectionIdx.cbegin(), m_collectionIdx.cend(), density.begin(), [&](auto idx) {
            return m_collectionDensities[idx];
        });
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            m_doseScore[i].addScoredEnergy(m_energyScore[i], volumes[i], density[i], factor);
        }
    }

    /**
     * @brief Transports a particle through the mesh until it exits or is absorbed,
     *        scoring deposited energy per tetrahedron.  Registers the track when
     *        @p P is ParticleTrack.
     * @param particle Particle to transport; modified in place.
     * @param state    Random number generator state.
     */
    template <ParticleType P>
    void transport(P& particle, RandomState& state)
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            m_tracker.registerParticle(particle);
        }

        siddonTransportStep(particle, state);
    }

    /**
     * @brief Reconstructs a TetrahedalMeshData snapshot from the current internal state
     *        (vertices, elements, collections, materials, and names).
     * @return A fully populated TetrahedalMeshData that can be passed to setData().
     */
    TetrahedalMeshData copyData() const
    {
        TetrahedalMeshData data;
        data.nodes = m_vertices;

        data.elements.resize(m_tetrahedrons.size());
        std::transform(std::execution::unseq, m_tetrahedrons.cbegin(), m_tetrahedrons.cend(), data.elements.begin(), [](const auto& t) { return t.verticeIdx; });

        data.collectionIndices = m_collectionIdx;
        data.collectionDensities = m_collectionDensities;
        data.collectionMaterialComposition.resize(m_collectionMaterials.size());
        std::transform(std::execution::unseq, m_collectionMaterials.cbegin(), m_collectionMaterials.cend(), data.collectionMaterialComposition.begin(), [](const auto& m) { return m.composition(); });
        data.collectionNames = m_collectionNames;

        return data;
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "TetMesh1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells) + std::to_string(FORCEDINTERACTION);
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
     * @brief Serializes the mesh (geometry, collections, materials, and dose scores)
     *        to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {

        auto buffer = Serializer::getEmptyBuffer();
        {
            std::vector<double> nodes_flat;
            nodes_flat.reserve(m_vertices.size() * 3);
            for (const auto& n : m_vertices)
                for (const auto& v : n)
                    nodes_flat.push_back(v);
            Serializer::serialize(nodes_flat, buffer);
        }
        {
            std::vector<std::uint32_t> elements_flat;
            elements_flat.reserve(m_tetrahedrons.size() * 4);
            for (const auto& n : m_tetrahedrons)
                for (const auto& v : n.verticeIdx)
                    elements_flat.push_back(v);
            Serializer::serialize(elements_flat, buffer);
        }
        Serializer::serialize(m_collectionIdx, buffer);
        Serializer::serialize(m_collectionDensities, buffer);
        {
            std::vector<std::map<std::uint8_t, double>> material_weights(m_collectionMaterials.size());
            std::transform(std::execution::unseq, m_collectionMaterials.cbegin(), m_collectionMaterials.cend(),
                material_weights.begin(), [](const auto& m) {
                    return m.composition();
                });
            Serializer::serializeMaterialWeights(material_weights, buffer);
        }
        Serializer::serialize(m_collectionNames, buffer);

        Serializer::serializeDoseScore(m_doseScore, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a mesh from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed mesh on success, or std::nullopt if setData() fails or
     *         the dose score vector size does not match the tetrahedron count.
     */
    static std::optional<TetrahedalMesh<NMaterialShells, LOWENERGYCORRECTION, FORCEDINTERACTION>> deserialize(std::span<const char> buffer)
    {

        TetrahedalMeshData data;
        { // nodes
            std::vector<double> nodes_flat;
            buffer = Serializer::deserialize(nodes_flat, buffer);
            data.nodes.resize(nodes_flat.size() / 3);
            for (std::size_t i = 0; i < data.nodes.size(); ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    data.nodes[i][j] = nodes_flat[i * 3 + j];
                }
            }
        }
        { // elements
            std::vector<std::uint32_t> elements_flat;
            buffer = Serializer::deserialize(elements_flat, buffer);
            data.elements.resize(elements_flat.size() / 4);
            for (std::size_t i = 0; i < data.elements.size(); ++i) {
                for (std::size_t j = 0; j < 4; ++j) {
                    data.elements[i][j] = elements_flat[i * 4 + j];
                }
            }
        }
        buffer = Serializer::deserialize(data.collectionIndices, buffer);
        buffer = Serializer::deserialize(data.collectionDensities, buffer);
        buffer = Serializer::deserializeMaterialWeights(data.collectionMaterialComposition, buffer);
        buffer = Serializer::deserialize(data.collectionNames, buffer);

        TetrahedalMesh<NMaterialShells, LOWENERGYCORRECTION, FORCEDINTERACTION> item;

        if (item.setData(data)) {
            buffer = Serializer::deserializeDoseScore(item.m_doseScore, buffer);
            if (item.m_doseScore.size() == item.m_tetrahedrons.size())
                return item;
        }
        return std::nullopt;
    }

protected:
    /// @brief Returns the index (0–2) of the component with the largest absolute value.
    static auto absargmax3(const std::array<double, 3>& arr)
    {
        const auto a = std::abs(arr[0]);
        const auto b = std::abs(arr[1]);
        const auto c = std::abs(arr[2]);
        return (b > a && b > c) ? 1 : ((c > a) ? 2 : 0);
    }

    /**
     * @brief Nudges the particle position by one ULP along its dominant direction to
     *        escape a face boundary.  Currently unused but kept for reference.
     * @tparam FORWARD When true, nudges forward along the direction; false nudges backward.
     */
    template <bool FORWARD = true>
    static void nudgeParticle(ParticleType auto& particle)
    {
        // Currently not used
        constexpr auto min = std::numeric_limits<double>::lowest();
        constexpr auto max = std::numeric_limits<double>::max();
        const auto i = absargmax3(particle.dir);
        if constexpr (FORWARD) {
            const auto dir = particle.dir[i] < 0.0 ? min : max;
            particle.pos[i] = std::nextafter(particle.pos[i], dir);
        } else {
            const auto dir = particle.dir[i] < 0.0 ? max : min;
            particle.pos[i] = std::nextafter(particle.pos[i], dir);
        }
    }

    /**
     * @brief Finds the tetrahedron that contains the particle's current position by
     *        shooting a reverse ray to the outer surface, then walking forward.
     *        May advance the particle slightly to escape a face if needed.
     * @param particle Particle whose position is to be located; may be translated.
     * @return Index of the tetrahedron containing the particle.
     */
    std::uint32_t intersectedTetrahedron(ParticleType auto& particle) const
    {
        // copy particle
        Particle p;
        p.pos = particle.pos;
        p.dir = vectormath::scale(particle.dir, -1.0);
        auto kdres = m_kdtree.intersect(p, m_vertices, m_outer_triangles, m_aabb);
        if (!kdres.valid()) {
            kdres = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, m_aabb);
            particle.translate(kdres.intersection + EPSILON);
        }
        const auto faceIdx = std::distance(&(m_outer_triangles[0]), kdres.item);
        auto currentTetIdx = m_outerTriangleTetMembership[faceIdx];
        // test if we are inside our tet;
        auto tet = walkTetrahedronLine(currentTetIdx, particle);
        return tet;
    }

    /**
     * @brief Walks the tetrahedral connectivity along the particle direction until the
     *        tetrahedron that straddles the current position is found.
     *        The particle may be translated by a small epsilon to resolve boundary cases.
     * @param currentIdx Starting tetrahedron index.
     * @param particle   Particle whose position defines the search point; may be translated.
     * @return Index of the tetrahedron that contains the particle.
     */
    std::uint32_t walkTetrahedronLine(std::uint32_t currentIdx, ParticleType auto& particle) const
    {
        std::array<double, 2> t = { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() };
        std::array<std::uint32_t, 2> faces;
        bool found = false;
        std::uint32_t oldIdx = std::numeric_limits<std::uint32_t>::max();
        while (!found) {
            const auto& tet = m_tetrahedrons[currentIdx];
            const auto& vIdx = tet.verticeIdx;
            const auto& v0 = m_vertices[vIdx[0]];
            const auto& v1 = m_vertices[vIdx[1]];
            const auto& v2 = m_vertices[vIdx[2]];
            const auto& v3 = m_vertices[vIdx[3]];
            std::uint32_t hit_counter = 0;
            const std::array<std::optional<double>, 4> hits = {
                basicshape::tetrahedron::intersectTriangle(v0, v1, v2, particle),
                basicshape::tetrahedron::intersectTriangle(v1, v0, v3, particle),
                basicshape::tetrahedron::intersectTriangle(v2, v3, v0, particle),
                basicshape::tetrahedron::intersectTriangle(v3, v2, v1, particle)
            };

            for (std::uint32_t i = 0; i < 4; ++i) {
                const auto& h = hits[i];
                if (h) {
                    if (*h < t[0]) {
                        t[0] = *h;
                        faces[0] = i;
                    }
                    if (*h > t[1]) {
                        t[1] = *h;
                        faces[1] = i;
                    }
                }
            }

            found = t[0] <= 0 && t[1] >= 0.0;
            if (!found) {
                // directions
                if (t[0] <= 0.0) {
                    if (currentIdx == tet.neighborIdx[faces[1]]) {
                        const auto dist = t[0] - EPSILON;
                        particle.translate(dist);
                    } else {
                        found = tet.neighborIdx[faces[1]] == oldIdx;
                        oldIdx = currentIdx;
                        currentIdx = tet.neighborIdx[faces[1]];
                    }
                } else { // (t[0] >= 0.0)
                    if (currentIdx == tet.neighborIdx[faces[0]]) {
                        const auto dist = t[0] + EPSILON;
                        particle.translate(dist);
                    } else {
                        found = tet.neighborIdx[faces[0]] == oldIdx;
                        oldIdx = currentIdx;
                        currentIdx = tet.neighborIdx[faces[0]];
                    }
                }
            }
        }
        return currentIdx;
    }

    /**
     * @brief Finds the exit face of the current tetrahedron and returns the path length
     *        to it and the index of the neighboring tetrahedron on the other side.
     *        Returns (length, currentIdx) when no neighbor exists (outer boundary).
     * @param currentIdx Index of the tetrahedron the particle is currently inside.
     * @param particle   Particle defining the ray direction.
     * @return {distance to exit face, index of next tetrahedron}.
     */
    std::tuple<double, std::uint32_t> nextTetrahedron(std::uint32_t currentIdx, const ParticleType auto& particle) const
    {
        const auto& tet = m_tetrahedrons[currentIdx];
        const auto& vIdx = tet.verticeIdx;
        const auto& v0 = m_vertices[vIdx[0]];
        const auto& v1 = m_vertices[vIdx[1]];
        const auto& v2 = m_vertices[vIdx[2]];
        const auto& v3 = m_vertices[vIdx[3]];

        // tet[0], tet[1], tet[2]
        // tet[1], tet[0], tet[3]
        // tet[2], tet[3], tet[0]
        // tet[3], tet[2], tet[1]

        double lenght = std::numeric_limits<double>::lowest();
        auto nextIdx = currentIdx;

        // optimize: we really only can hit one, early exit?
        const std::array<std::optional<double>, 4> hits = {
            basicshape::tetrahedron::intersectTriangle(v0, v1, v2, particle),
            basicshape::tetrahedron::intersectTriangle(v1, v0, v3, particle),
            basicshape::tetrahedron::intersectTriangle(v2, v3, v0, particle),
            basicshape::tetrahedron::intersectTriangle(v3, v2, v1, particle)
        };

        for (std::uint32_t i = 0; i < 4; ++i) {
            const auto& h = hits[i];
            if (h) {
                if (*h > lenght) {
                    lenght = *h;
                    nextIdx = tet.neighborIdx[i];
                }
            }
        }
        return std::make_pair(lenght, nextIdx);
    }

    /**
     * @brief Per-tetrahedron analog transport: samples a free path in each tetrahedron
     *        and either interacts or crosses to the next cell.
     *        Uses cumulative probability steps for numerical stability on small tetrahedrons.
     *        Unused but kept for reference
     * @param particle Particle to transport; modified in place.
     * @param state    Random number generator state.
     */
    void siddonTransport(ParticleType auto& particle, RandomState& state)
    {
        // we do cummulative steps instead of stepping on each tet due to numerical stability of tet intersections
        // when tets become very small

        auto currentTetIdx = intersectedTetrahedron(particle);
        auto collIdx = m_collectionIdx[currentTetIdx];
        bool still_inside = true;
        bool updateAtt = true;
        AttenuationValues attenuation;

        do {
            if (updateAtt) {
                attenuation = m_collectionMaterials[collIdx].attenuationValues(particle.energy);
                updateAtt = false;
            }
            const auto [borderlen, nextTetIdx] = nextTetrahedron(currentTetIdx, particle);
            const auto attSum = attenuation.sum();
            const auto steplen = -std::log(state.randomUniform()) / (attSum * m_collectionDensities[collIdx]);

            if constexpr (FORCEDINTERACTION) {
                // Forced photoelectric effect
                const auto relativePeProbability = attenuation.photoelectric / attSum;
                const auto probNotInteraction = std::exp(-attSum * m_collectionDensities[collIdx] * borderlen);
                m_energyScore[currentTetIdx].scoreEnergy(particle.energy * particle.weight * (1 - probNotInteraction) * relativePeProbability);
            }

            if (steplen < borderlen) {
                // interaction
                particle.translate(steplen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(attenuation, particle, m_collectionMaterials[collIdx], state);
                if constexpr (FORCEDINTERACTION) {
                    if (intRes.interactionWasIncoherent) // We have already scored PE events
                        m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                } else {
                    m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                }
                still_inside = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            } else {
                // translate
                still_inside = nextTetIdx != currentTetIdx;
                if (still_inside && borderlen > 0.0) {
                    particle.translate(borderlen);
                    if (collIdx != m_collectionIdx[nextTetIdx]) {
                        updateAtt = true;
                        collIdx = m_collectionIdx[nextTetIdx];
                    }
                    currentTetIdx = nextTetIdx;
                } else {
                    particle.border_translate(borderlen);
                }
            }
        } while (still_inside);
    }

    /**
     * @brief Cumulative-probability Siddon transport: accumulates survival probability
     *        across tetrahedron boundaries and resolves an interaction only when the
     *        threshold is crossed, improving stability for very small cells.
     *        This is the method called by transport().
     * @param particle Particle to transport; modified in place.
     * @param state    Random number generator state.
     */
    void siddonTransportStep(ParticleType auto& particle, RandomState& state)
    {
        // we do cummulative steps instead of stepping on each tet due to numerical stability of tet intersections
        // when tets become very small

        auto currentTetIdx = intersectedTetrahedron(particle);
        auto collIdx = m_collectionIdx[currentTetIdx];
        bool still_inside = true;
        bool updateAtt = true;
        AttenuationValues attenuation;
        double attSum;

        double probability_thres = state.randomUniform();
        double interaction_probability = 1.0;
        double steplenght = 0.0;

        do {
            if (updateAtt) {
                attenuation = m_collectionMaterials[collIdx].attenuationValues(particle.energy);
                attSum = attenuation.sum() * m_collectionDensities[collIdx];
                updateAtt = false;
            }

            const auto [borderlen, nextTetIdx] = nextTetrahedron(currentTetIdx, particle);
            if (borderlen > steplenght) {
                // valid border test

                // update next prob
                const auto delta_prob = std::exp(-attSum * (borderlen - steplenght));
                interaction_probability *= delta_prob;

                if constexpr (FORCEDINTERACTION) {
                    // Forced photoelectric effect
                    const auto relativePeProbability = attenuation.photoelectric * m_collectionDensities[collIdx] / attSum;
                    m_energyScore[currentTetIdx].scoreEnergy(particle.energy * particle.weight * (1 - delta_prob) * relativePeProbability);
                }

                if (interaction_probability < probability_thres) { // interaction
                    const auto delta_dist = std::log(probability_thres / interaction_probability) / attSum; // correcting for overspilling
                    const auto dist = borderlen - delta_dist;
                    // moving particle
                    particle.translate(dist);
                    // resetting
                    probability_thres = state.randomUniform();
                    interaction_probability = 1.0;
                    steplenght = 0.0;
                    // Do interaction
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(attenuation, particle, m_collectionMaterials[collIdx], state);
                    if constexpr (FORCEDINTERACTION) {
                        if (intRes.interactionWasIncoherent) // We have already scored PE events
                            m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                    } else {
                        m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                    }
                    still_inside = intRes.particleAlive;
                    updateAtt = intRes.particleEnergyChanged;
                } else {
                    steplenght = borderlen;
                    still_inside = currentTetIdx != nextTetIdx;
                    if (still_inside) {
                        currentTetIdx = nextTetIdx;
                        updateAtt = collIdx != m_collectionIdx[nextTetIdx];
                    } else {
                        particle.border_translate(steplenght);
                    }
                }
            } else {
                // something is wrong, advance
                // still_inside = currentTetIdx != nextTetIdx;
                // if (still_inside) {
                //    currentTetIdx = nextTetIdx;
                //    updateAtt = collIdx != m_collectionIdx[nextTetIdx];
                //} else {
                //    particle.border_translate(steplenght);
                //}
                // per now simply kill particle, happens once per 1E7 photons
                particle.energy = 0;
                still_inside = false;
            }
        } while (still_inside);
    }

    /// @brief Recomputes the AABB by scanning all vertices.
    void calculateAABB()
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = std::numeric_limits<double>::max();
            m_aabb[i + 3] = std::numeric_limits<double>::lowest();
        }

        for (const auto& v : m_vertices) {
            for (std::size_t i = 0; i < 3; ++i) {
                m_aabb[i] = std::min(m_aabb[i], v[i]);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], v[i]);
            }
        }
    }

    /**
     * @brief Builds the internal tetrahedron list and computes neighbor connectivity
     *        by matching shared faces between adjacent tetrahedrons.
     * @param data Source mesh data; nodes and elements are consumed directly.
     */
    void makeStructure(const TetrahedalMeshData& data)
    {
        m_vertices = data.nodes;

        m_tetrahedrons.resize(data.elements.size());
        m_tetrahedrons.shrink_to_fit();

        std::vector<std::uint32_t> indices(m_tetrahedrons.size());
        std::iota(indices.begin(), indices.end(), 0);

        std::transform(std::execution::par_unseq, data.elements.cbegin(), data.elements.cend(), indices.cbegin(), m_tetrahedrons.begin(), [](const auto& tInds, auto i) {
            return Tetrahedron { .verticeIdx = tInds, .neighborIdx = { i, i, i, i } };
        });

        if (m_tetrahedrons.size() < 2) {
            return;
        }

        // finding neighbors
        struct Face {
            Face() { }
            Face(std::uint32_t x, std::uint32_t y, std::uint32_t z, std::uint32_t i, std::uint32_t num)
            {
                verticeIdx = { x, y, z };
                tetIdx = i;
                faceNumber = num;
            }
            void sort() noexcept
            {
                if (verticeIdx[0] > verticeIdx[1])
                    std::swap(verticeIdx[0], verticeIdx[1]);
                if (verticeIdx[1] > verticeIdx[2])
                    std::swap(verticeIdx[1], verticeIdx[2]);
                if (verticeIdx[0] > verticeIdx[1])
                    std::swap(verticeIdx[0], verticeIdx[1]);
            }
            auto operator<=>(const Face& other) const
            {
                return verticeIdx <=> other.verticeIdx;
            }
            bool operator==(const Face& other) const
            {
                return verticeIdx == other.verticeIdx;
            }
            std::array<std::uint32_t, 3> verticeIdx;
            std::uint32_t tetIdx = 0;
            std::uint32_t faceNumber = 0;
        };

        // vector of faces
        std::vector<Face> faces;
        faces.reserve(m_tetrahedrons.size() * 4);
        for (std::uint32_t i = 0; i < m_tetrahedrons.size(); ++i) {
            const auto& tet = m_tetrahedrons[i].verticeIdx;
            faces.push_back({ tet[0], tet[1], tet[2], i, std::uint32_t { 0 } });
            faces.push_back({ tet[1], tet[0], tet[3], i, std::uint32_t { 1 } });
            faces.push_back({ tet[2], tet[3], tet[0], i, std::uint32_t { 2 } });
            faces.push_back({ tet[3], tet[2], tet[1], i, std::uint32_t { 3 } });
        }
        // sorting faces idx for comparison
        std::for_each(std::execution::par_unseq, faces.begin(), faces.end(), [](auto& f) { f.sort(); });
        std::sort(std::execution::par_unseq, faces.begin(), faces.end());

        // most faces will have one and only one neighbor, others zero
        auto first = faces.begin();
        auto second = first + 1;
        do {
            if (*first == *second) {
                m_tetrahedrons[first->tetIdx].neighborIdx[first->faceNumber] = second->tetIdx;
                m_tetrahedrons[second->tetIdx].neighborIdx[second->faceNumber] = first->tetIdx;
            }
            first = second++;
        } while (second != faces.end());
    }

    /**
     * @brief Identifies all boundary (outer) triangular faces — those whose neighbor
     *        index equals their own tetrahedron — and builds the contour KD-tree from them.
     */
    void buildContourKDTree()
    {
        m_outer_triangles.clear();
        m_outerTriangleTetMembership.clear();
        m_outer_triangles.reserve(m_tetrahedrons.size());
        m_outerTriangleTetMembership.reserve(m_tetrahedrons.size());
        for (std::uint32_t i = 0; i < m_tetrahedrons.size(); ++i) {
            const auto& tet = m_tetrahedrons[i].verticeIdx;
            if (m_tetrahedrons[i].neighborIdx[0] == i) {
                m_outer_triangles.push_back({ tet[0], tet[1], tet[2] });
                m_outerTriangleTetMembership.push_back(i);
            }
            if (m_tetrahedrons[i].neighborIdx[1] == i) {
                m_outer_triangles.push_back({ tet[1], tet[0], tet[3] });
                m_outerTriangleTetMembership.push_back(i);
            }
            if (m_tetrahedrons[i].neighborIdx[2] == i) {
                m_outer_triangles.push_back({ tet[2], tet[3], tet[0] });
                m_outerTriangleTetMembership.push_back(i);
            }
            if (m_tetrahedrons[i].neighborIdx[3] == i) {
                m_outer_triangles.push_back({ tet[3], tet[2], tet[1] });
                m_outerTriangleTetMembership.push_back(i);
            }
        }
        m_outer_triangles.shrink_to_fit();
        m_outerTriangleTetMembership.shrink_to_fit();
        m_kdtree.setData(m_vertices, m_outer_triangles, 8);
    }

private:
    struct Tetrahedron {
        std::array<std::uint32_t, 4> verticeIdx = { 0, 0, 0, 0 };
        std::array<std::uint32_t, 4> neighborIdx = { 0, 0, 0, 0 };
    };

    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };

    // tetrahedrons
    std::vector<std::array<double, 3>> m_vertices;
    std::vector<Tetrahedron> m_tetrahedrons;

    // structure for outer contour
    std::vector<std::array<std::uint32_t, 3>> m_outer_triangles; // outer triangles
    std::vector<std::uint32_t> m_outerTriangleTetMembership; // same size as outer triangles
    TetrahedalMeshContourKDTree m_kdtree;

    std::vector<std::uint32_t> m_collectionIdx; // same size as tetrahedrons
    std::vector<EnergyScore> m_energyScore; // same size as tetrahedrons
    std::vector<DoseScore> m_doseScore; // same size as tetrahedrons

    std::vector<double> m_collectionDensities; // same size as n collections
    std::vector<Material<NMaterialShells>> m_collectionMaterials; // same size as n collections
    std::vector<std::string> m_collectionNames; // same size as n collections
    ParticleTracker m_tracker;
    std::uint32_t m_collectionIndexForDisplay = std::numeric_limits<std::uint32_t>::max();
};
}
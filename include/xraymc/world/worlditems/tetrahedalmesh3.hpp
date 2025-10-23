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
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcollection.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh3/tetrahedalmeshcontourkdtree.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <optional>

namespace xraymc {

template <int NMaterialShells = 36, int LOWENERGYCORRECTION = 2, bool FORCEDINTERACTION = false>
class TetrahedalMesh3 {
    // static constexpr auto EPSILON = std::numeric_limits<double>::epsilon() * 1000;
    static constexpr double EPSILON = 1E-10;
    // static constexpr double EPSILON = GEOMETRIC_ERROR();

public:
    TetrahedalMesh3() { };
    TetrahedalMesh3(const TetrahedalMeshData& data)
    {
        setData(data);
    }

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

    void rotate(const std::array<double, 3>& axis, double angle)
    {
        std::transform(std::execution::par_unseq, m_vertices.begin(), m_vertices.end(), m_vertices.begin(), [=](const auto& v) {
            return vectormath::rotate(v, axis, angle);
        });
        calculateAABB();
        m_kdtree.setData(m_vertices, m_outer_triangles, 8);
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    std::array<double, 3> center() const
    {
        std::array<double, 3> res = {
            (m_aabb[0] + m_aabb[3]) * 0.5,
            (m_aabb[1] + m_aabb[4]) * 0.5,
            (m_aabb[2] + m_aabb[5]) * 0.5
        };
        return res;
    }

    double tetrahedalVolume(std::uint32_t index) const
    {
        const auto& tet = m_tetrahedrons.at(index).verticeIdx;
        const auto& v0 = m_vertices[tet[0]];
        const auto& v1 = m_vertices[tet[1]];
        const auto& v2 = m_vertices[tet[2]];
        const auto& v3 = m_vertices[tet[3]];
        return basicshape::tetrahedron::volume(v0, v1, v2, v3);
    }

    const std::vector<std::uint32_t>& outerContourTetrahedronIndices() const
    {
        return m_outerTriangleTetMembership;
    }

    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

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

    WorldIntersectionResult intersect(const ParticleType auto& particle) const
    {
        WorldIntersectionResult res {};
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

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& particle) const
    {
        VisualizationIntersectionResult<U> res {};
        if (auto inter = basicshape::AABB::intersectForwardInterval(particle, m_aabb); inter) {
            auto kdres = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, *inter);
            if (kdres.valid()) {
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
            }
        };
        return res;
    }

    std::uint32_t numberOfThetrahedrons() const
    {
        return static_cast<std::uint32_t>(m_tetrahedrons.size());
    }

    EnergyScore energyScored(std::size_t index) const
    {
        return m_energyScore.at(index);
    }

    DoseScore doseScored(std::size_t index) const
    {
        return m_doseScore.at(index);
    }

    void clearDoseScored()
    {
        std::for_each(std::execution::par_unseq, m_doseScore.begin(), m_doseScore.end(), [](auto& d) {
            d.clear();
        });
    }

    void clearEnergyScored()
    {
        std::for_each(std::execution::par_unseq, m_energyScore.begin(), m_energyScore.end(), [](auto& d) {
            d.clear();
        });
    }

    void addEnergyScoredToDoseScore(double factor)
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

    template <ParticleType P>
    void transport(P& particle, RandomState& state)
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            m_tracker.registerParticle(particle);
        }
        siddonTransport(particle, state);
    }

protected:
    static auto absargmax3(const std::array<double, 3>& arr)
    {
        const auto a = std::abs(arr[0]);
        const auto b = std::abs(arr[1]);
        const auto c = std::abs(arr[2]);
        return (b > a && b > c) ? 1 : ((c > a) ? 2 : 0);
    }

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

    std::uint32_t intersectedTetrahedron(ParticleType auto& particle) const
    {
        // copy particle
        Particle p;
        p.pos = particle.pos;
        p.dir = vectormath::scale(particle.dir, -1.0);
        auto kdres = m_kdtree.intersect(p, m_vertices, m_outer_triangles, m_aabb);
        if (!kdres.valid()) {
            kdres = m_kdtree.intersect(particle, m_vertices, m_outer_triangles, m_aabb);
            particle.border_translate(kdres.intersection);
        }
        const auto faceIdx = std::distance(&(m_outer_triangles[0]), kdres.item);
        auto currentTetIdx = m_outerTriangleTetMembership[faceIdx];
        // test if we are inside our tet;
        auto tet = walkTetrahedronLine(currentTetIdx, particle);
        return tet;
    }

    std::uint32_t walkTetrahedronLine(std::uint32_t currentIdx, ParticleType auto& particle) const
    {
        std::array<double, 2> t;
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
                    faces[hit_counter] = i;
                    t[hit_counter++] = *h;
                }
            }
            if (t[0] > t[1]) {
                t = { t[1], t[0] };
                faces = { faces[1], faces[0] };
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

        auto lenght = std::numeric_limits<double>::lowest();
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

    void siddonTransport(ParticleType auto& particle, RandomState& state)
    {
        // we do cummulative steps instead of stepping on each tet due to numerical stability of tet intersections
        // when tets become very small

        std::uint32_t currentTetIdx = intersectedTetrahedron(particle);
        bool still_inside = true;
        double attSum;
        double relativePEprobability;
        AttenuationValues attenuation;
        std::uint32_t collIdx;
        bool updateAtt = true;
        double prob_thres = state.randomUniform();
        double prob = 1;
        double travel_distance = 0;

        while (still_inside) {
            if (updateAtt) {
                if constexpr (FORCEDINTERACTION) {
                    collIdx = m_collectionIdx[currentTetIdx];
                    attenuation = m_collectionMaterials[collIdx].attenuationValues(particle.energy);
                    const auto attTot = attenuation.sum();
                    relativePEprobability = attenuation.photoelectric / attTot;
                    attSum = attTot * m_collectionDensities[collIdx];
                    updateAtt = false;
                } else {
                    collIdx = m_collectionIdx[currentTetIdx];
                    attenuation = m_collectionMaterials[collIdx].attenuationValues(particle.energy);
                    attSum = attenuation.sum() * m_collectionDensities[collIdx];
                    updateAtt = false;
                }
            }

            const auto [borderlen, nextTetIdx] = nextTetrahedron(currentTetIdx, particle);
            const auto delta_prob = std::exp(-attSum * (borderlen - travel_distance));
            prob *= delta_prob; // accumulating probability like the Heart of Gold

            if constexpr (FORCEDINTERACTION) {
                // Forced photoel
                const auto pe_prob = (1.0 - delta_prob) * relativePEprobability;
                m_energyScore[currentTetIdx].scoreEnergy(particle.energy * particle.weight * pe_prob);
            }

            if (prob <= prob_thres) { // interaction
                const auto delta_dist = -std::log(prob / prob_thres) / attSum; // correcting for overspilling
                const auto dist = borderlen - delta_dist;
                // moving particle
                particle.translate(dist);
                travel_distance = 0; // resetting travelled distance
                // Do interaction
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(attenuation, particle, m_collectionMaterials[collIdx], state);
                if constexpr (FORCEDINTERACTION) {
                    if (!intRes.interactionWasPhotoelectric && intRes.energyImparted > 0.0)
                        // We have already scored photoelectric energy, skipping scoring photoelectric to prevent bias
                        m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                } else {
                    if (intRes.energyImparted > 0.0)
                        m_energyScore[currentTetIdx].scoreEnergy(intRes.energyImparted);
                }
                if (intRes.particleAlive) {
                    updateAtt = intRes.particleEnergyChanged;
                    // updating prob threshold and prob  for next round
                    prob_thres = state.randomUniform();
                    prob = 1;
                } else {
                    still_inside = false;
                }
            } else {
                updateAtt = collIdx != m_collectionIdx[nextTetIdx];
                still_inside = currentTetIdx != nextTetIdx;
                currentTetIdx = nextTetIdx;
                travel_distance = borderlen;
            }
        }
        // Particle  either dead or we exit mesh
        if (particle.energy != 0.0)
            particle.border_translate(travel_distance);
    }

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

    void makeStructure(const TetrahedalMeshData& data)
    {
        m_vertices = data.nodes;

        m_tetrahedrons.resize(data.elements.size());
        m_tetrahedrons.shrink_to_fit();
        for (std::uint32_t i = 0; i < data.elements.size(); ++i) {
            m_tetrahedrons[i].verticeIdx = data.elements[i];
            // set all neighbors to self
            m_tetrahedrons[i].neighborIdx = { i, i, i, i };
        }

        if (m_tetrahedrons.size() < 2) {
            return;
        }

        // finding neighbors
        struct Face {
            Face(std::uint32_t x, std::uint32_t y, std::uint32_t z, std::uint32_t i, std::uint32_t num)
            {
                verticeIdx = { x, y, z };
                tetIdx = i;
                faceNumber = num;
            }
            void sort()
            {
                std::sort(verticeIdx.begin(), verticeIdx.end());
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

    void buildContourKDTree()
    {
        m_outer_triangles.clear();
        m_outerTriangleTetMembership.clear();
        for (std::uint32_t i = 0; i < m_tetrahedrons.size(); ++i) {
            const auto& tet = m_tetrahedrons[i].verticeIdx;
            if (m_tetrahedrons[i].neighborIdx[0] == i) {
                m_outer_triangles.push_back({ tet[0], tet[1], tet[2] });
                m_outerTriangleTetMembership.push_back(i);
            } else if (m_tetrahedrons[i].neighborIdx[1] == i) {
                m_outer_triangles.push_back({ tet[1], tet[0], tet[3] });
                m_outerTriangleTetMembership.push_back(i);
            } else if (m_tetrahedrons[i].neighborIdx[2] == i) {
                m_outer_triangles.push_back({ tet[2], tet[3], tet[0] });
                m_outerTriangleTetMembership.push_back(i);
            } else if (m_tetrahedrons[i].neighborIdx[3] == i) {
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
};
}
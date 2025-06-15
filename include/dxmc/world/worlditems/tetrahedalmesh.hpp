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

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/particletracker.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshgrid.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

struct TetrahedalMeshCollection {
    double density = 0;
    double volume = 0;
    double dose = 0;
    std::string name;
};

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FLUENCESCORING = true>
class TetrahedalMesh {
public:
    TetrahedalMesh()
    {
    }

    TetrahedalMesh(const TetrahedalMeshData& data, std::array<int, 3> dimensions = { 8, 8, 8 })
    {
        setData(data, dimensions);
    }

    bool setData(const TetrahedalMeshData& data, std::array<int, 3> dimensions)
    {
        if (!data.valid())
            return false;

        std::vector<Tetrahedron> tets(data.elements.size());
        std::transform(std::execution::par_unseq, data.elements.cbegin(), data.elements.cend(), data.collectionIndices.cbegin(), tets.begin(), [&data](const auto& e, const auto idx) {
            const auto& v0 = data.nodes[e[0]];
            const auto& v1 = data.nodes[e[1]];
            const auto& v2 = data.nodes[e[2]];
            const auto& v3 = data.nodes[e[3]];
            return Tetrahedron { v0, v1, v2, v3, static_cast<std::uint16_t>(idx) };
        });

        m_collectionDensity = data.collectionDensities;

        m_collectionMaterial.clear();
        m_collectionMaterial.reserve(data.collectionMaterialComposition.size());
        for (const auto& weights : data.collectionMaterialComposition) {
            auto mat_opt = dxmc::Material<NMaterialShells>::byWeight(weights);
            if (mat_opt)
                m_collectionMaterial.push_back(mat_opt.value());
            else
                return false;
        }

        m_grid.setData(tets, dimensions);

        if constexpr (!FLUENCESCORING) {
            generateWoodcockStepTable();
        }
        return true;
    }

    bool setData(const std::vector<Tetrahedron>& tets, const std::vector<double>& collectionDensities, const std::vector<Material<NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {}, std::array<int, 3> depth = { 8, 8, 8 })
    {
        // finding max collectionIdx and Material index
        const auto maxCollectionIdx = std::transform_reduce(
            std::execution::par_unseq, tets.cbegin(), tets.cend(), std::uint16_t { 0 },
            [](const auto lh, const auto rh) { return std::max(lh, rh); },
            [](const auto& t) { return t.collection(); });
        const auto maxMaterialIdx = std::transform_reduce(
            std::execution::par_unseq, tets.cbegin(), tets.cend(), std::uint16_t { 0 },
            [](const auto lh, const auto rh) { return std::max(lh, rh); },
            [](const auto& t) { return t.materialIndex(); });
        if (collectionDensities.size() <= maxCollectionIdx || materials.size() <= maxMaterialIdx)
            return false;
        m_collectionMaterial = materials;
        m_collectionDensity.reserve(collectionDensities.size());

        std::vector<std::atomic<double>> volumes(maxCollectionIdx + 1);
        std::for_each(std::execution::par_unseq, volumes.begin(), volumes.end(), [](auto& v) { v.store(double { 0 }); });
        std::for_each(std::execution::par_unseq, tets.cbegin(), tets.cend(), [&volumes](const auto& tet) {
            const auto idx = tet.collection();
            volumes[idx].fetch_add(tet.volume());
        });

        for (std::size_t i = 0; i <= maxCollectionIdx; ++i) {
            const auto d = collectionDensities[i];
            const auto v = volumes[i].load();
            m_collectionDensity.emplace_back(d, v);
        }

        if (collectionNames.size() == m_collectionDensity.size())
            m_collectionNames = collectionNames;
        else
            m_collectionNames.resize(m_collectionDensity.size());

        m_grid.setData(tets, depth);

        if constexpr (!FLUENCESCORING) {
            generateWoodcockStepTable();
        }

        return true;
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_grid.translate(dist);
    }

    std::array<double, 3> center() const
    {
        const auto& aabb = m_grid.AABB();
        const auto [low, high] = vectormath::splice(aabb);
        auto c = vectormath::add(low, high);
        return vectormath::scale(c, 0.5);
    }

    const std::array<double, 6>& AABB() const
    {
        return m_grid.AABB();
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        if (const auto kres = m_grid.intersect(p); kres.valid()) {
            res.intersection = kres.intersection;
            res.rayOriginIsInsideItem = kres.rayOriginIsInsideItem;
            res.intersectionValid = kres.valid();
        }
        return res;
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        VisualizationIntersectionResult<U> res;
        if (const auto kres = m_grid.intersect(p); kres.valid()) {
            res.intersection = kres.intersection;
            res.rayOriginIsInsideItem = kres.rayOriginIsInsideItem;
            res.intersectionValid = kres.valid();
            res.value = kres.item->doseScored().dose();
            const auto hit_pos = vectormath::add(p.pos, vectormath::scale(p.dir, kres.intersection));
            res.normal = kres.item->normal(hit_pos);
        }
        return res;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        const auto tets = m_grid.tetrahedrons();
        return tets.at(index).energyScored();
    }

    void clearEnergyScored()
    {
        m_grid.clearEnergyScored();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        auto& tets = m_grid.tetrahedrons();
        std::for_each(std::execution::par_unseq, tets.begin(), tets.end(), [calibration_factor, this](auto& tet) {
            const auto cidx = tet.collection();
            const auto dens = this->m_collectionDensity[cidx];
            tet.addEnergyScoredToDoseScore(dens, calibration_factor);
        });
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        const auto& tets = m_grid.tetrahedrons();
        return tets.at(index).doseScored();
    }

    void clearDoseScored()
    {
        m_grid.clearDoseScored();
    }

    template <ParticleType P>
    void transport(P& p, RandomState& state)
    {
        if constexpr (std::is_same_v<P, ParticleTrack>) {
            m_tracker.registerParticle(p);
        }
        if constexpr (FLUENCESCORING)
            transportSiddonForced(p, state);
        else
            transportWoodcock(p, state);
    }

    std::size_t numberOfCollections() const { return m_collectionDensity.size(); }
    std::size_t numberOfTetrahedra() const { return m_grid.tetrahedrons().size(); }
    const std::vector<Tetrahedron>& tetrahedrons() const
    {
        return m_grid.tetrahedrons();
    }
    std::size_t maxThetrahedronsVoxelCount() const
    {
        return m_grid.maxThetrahedronsVoxelCount();
    }

    std::vector<TetrahedalMeshCollection> collectionData() const
    {
        std::vector<TetrahedalMeshCollection> data(m_collectionDensity.size());
        const auto& tets = m_grid.tetrahedrons();
        for (std::size_t i = 0; i < m_collectionDensity.size(); ++i) {
            data[i].density = m_collectionDensity[i].density;
            data[i].name = m_collectionNames[i];
            const auto collectionIdx = static_cast<std::uint16_t>(i);
            data[i].volume = std::transform_reduce(std::execution::par_unseq, tets.cbegin(), tets.cend(), 0.0, std::plus {}, [collectionIdx](const auto& t) -> double {
                return t.collection() == collectionIdx ? t.volume() : 0.0;
            });
            const auto cdens = data[i].density;
            if (data[i].density * data[i].volume > 0.0) {
                data[i].dose = std::transform_reduce(std::execution::par_unseq, tets.cbegin(), tets.cend(), 0.0, std::plus {}, [collectionIdx, cdens](const auto& t) -> double {
                    if (t.collection() == collectionIdx) {
                        const auto tetmass = t.volume() * cdens;
                        const auto energyImparted = t.doseScored().dose() * tetmass;
                        return energyImparted;
                    } else
                        return 0.0;
                });
                data[i].dose /= data[i].density * data[i].volume;
            } else {
                data[i].dose = 0;
            }
        }
        return data;
    }

    void setMaterial(const Material<NMaterialShells>& material, std::size_t index)
    {
        if (index < m_collectionMaterial.size())
            m_collectionMaterial[index] = material;
    }

    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

protected:
    void generateWoodcockStepTable()
    {
        std::vector<double> energy;
        {
            auto e = std::log(MIN_ENERGY());
            const auto emax = std::log(MAX_ENERGY());
            const auto estep = (emax - e) / 10;
            while (e <= emax) {
                energy.push_back(e);
                e += estep;
            }
        }
        // adding edges;
        for (const auto& mat : m_collectionMaterial) {
            for (std::size_t i = 0; i < mat.numberOfShells(); ++i) {
                const auto& shell = mat.shell(i);
                const auto e = shell.bindingEnergy + 0.01;
                if (e > MIN_ENERGY()) {
                    energy.push_back(std::log(e));
                }
            }
        }
        std::sort(energy.begin(), energy.end());
        auto remove = std::unique(energy.begin(), energy.end());
        energy.erase(remove, energy.end());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), energy.begin(), [](const auto e) { return std::exp(e); });

        // finding max attenuation for each energy
        std::vector<double> att(energy.size(), 0.0);
        for (std::size_t mIdx = 0; mIdx < m_collectionMaterial.size(); ++mIdx) {
            const auto& mat = m_collectionMaterial[mIdx];
            const auto d = m_collectionDensity[mIdx];
            for (std::size_t i = 0; i < energy.size(); ++i) {
                const auto aval = mat.attenuationValues(energy[i]);
                att[i] = std::max(aval.sum() * d, att[i]);
            }
        }
        std::vector<std::pair<double, double>> data(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.cbegin(), data.begin(), [](const auto e, const auto a) {
            return std::make_pair(e, a);
        });

        m_woodcockStepTableLin = data;
    }

    void transportWoodcock(ParticleType auto& p, RandomState& state)
    {
        bool still_inside = true;
        double attMaxInv;
        bool updateAtt = true;
        while (still_inside) {
            if (updateAtt) {
                attMaxInv = 1 / interpolate(m_woodcockStepTableLin, p.energy);
                updateAtt = false;
            }

            // making interaction step
            const auto steplen = -std::log(state.randomUniform()) * attMaxInv;
            p.translate(steplen);

            // finding current tet
            auto* currentTet = m_grid.pointInside(p.pos);

            if (currentTet) { // is interaction virtual?
                const auto collectionIdx = currentTet->collection();
                const auto attenuation = m_collectionMaterial[collectionIdx].attenuationValues(p.energy);
                const auto attSum = attenuation.sum() * m_collectionDensity[collectionIdx];
                if (state.randomUniform() < attSum * attMaxInv) {
                    // we have a real interaction
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(attenuation, p, m_collectionMaterial[collectionIdx], state);
                    currentTet->scoreEnergy(intRes.energyImparted);
                    still_inside = intRes.particleAlive;
                    updateAtt = intRes.particleEnergyChanged;
                }
            } else {
                // we are outside item, backtrack and return
                const Particle pback = { .pos = p.pos, .dir = vectormath::scale(p.dir, -1.0) };
                const auto inter_back = intersect(pback);
                if (inter_back.valid()) {
                    p.border_translate(-inter_back.intersection);
                }
                still_inside = false;
            }
        }
    }

    void transportSiddonForced(ParticleType auto& p, RandomState& state)
    {
        Tetrahedron* tet = m_grid.pointInside(p.pos);

        while (tet) {
            const auto& material = m_collectionMaterial[tet->collection()];
            const auto& density = m_collectionDensity[tet->collection()];
            const auto intLen = tet->intersect(p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, density, p, material, state);
            tet->scoreEnergy(intRes.energyImparted);
            tet = intRes.particleAlive ? m_grid.pointInside(p.pos) : nullptr;
        }
    }

    void transportSiddon(ParticleType auto& p, RandomState& state)
    {
        Tetrahedron* tet = m_grid.pointInside(p.pos);
        std::uint16_t currentCollection;
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (tet) {
            if (updateAtt) {
                currentCollection = tet->collection();
                att = m_collectionMaterial[currentCollection].attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_collectionDensity[currentCollection]);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = tet->intersect(p).intersection;
            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto& material = m_collectionMaterial[currentCollection];
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, material, state);
                tet->scoreEnergy(intRes.energyImparted);
                if (intRes.particleAlive)
                    updateAtt = intRes.particleEnergyChanged;
                else
                    tet = nullptr;
            } else {
                // transport to border
                p.border_translate(intLen);
                tet = m_grid.pointInside(p.pos);
                if (tet)
                    updateAtt = currentCollection != tet->collection();
            }
        }
    }

private:
    TetrahedalMeshGrid m_grid;
    std::vector<std::pair<double, double>> m_woodcockStepTableLin;
    std::vector<double> m_collectionDensity;
    std::vector<Material<NMaterialShells>> m_collectionMaterial;
    std::vector<std::string> m_collectionNames;
    ParticleTracker m_tracker;
};
}

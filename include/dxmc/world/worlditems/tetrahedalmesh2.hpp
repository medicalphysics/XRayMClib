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
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh2/tetrahedalmeshgrid2.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FLUENCESCORING = true>
class TetrahedalMesh2 {
public:
    TetrahedalMesh2()
    {
    }
    TetrahedalMesh2(const TetrahedalMeshData& data, std::array<uint32_t, 3> dimensions = { 128, 128, 128 })
    {
        setData(data, dimensions);
    }

    bool setData(const TetrahedalMeshData& data, std::array<std::uint32_t, 3> dimensions)
    {
        if (!data.valid())
            return false;

        m_collectionMaterials.clear();
        m_collectionMaterials.reserve(data.collectionMaterialComposition.size());
        for (const auto& weights : data.collectionMaterialComposition) {
            auto mat_opt = dxmc::Material<NMaterialShells>::byWeight(weights);
            if (mat_opt)
                m_collectionMaterials.push_back(mat_opt.value());
            else
                return false;
        }

        m_grid.setData(data, dimensions);
        m_collectionIdx = data.collectionIndices;
        m_collectionDensities = data.collectionDensities;
        m_collectionNames = data.collectionNames;
        m_doseScore.resize(m_collectionIdx.size());
        m_energyScore.resize(m_collectionIdx.size());
        return true;
    }

    void translate(const std::array<double, 3>& vec)
    {
        m_grid.translate(vec);
    }

    std::array<double, 3> center() const
    {
        const auto& aabb = m_grid.AABB();
        std::array<double, 3> c = {
            (aabb[0] + aabb[3]) * 0.5,
            (aabb[1] + aabb[4]) * 0.5,
            (aabb[2] + aabb[5]) * 0.5
        };
        return c;
    }
    const std::array<double, 6>& AABB() const
    {
        return m_grid.AABB();
    }

    WorldIntersectionResult intersect(ParticleType auto& p) const
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
            const auto index = *kres.item;
            const auto& dose = m_doseScore[index];
            res.value = dose.dose();
            const auto hit_pos = vectormath::add(p.pos, vectormath::scale(p.dir, kres.intersection));

            const auto& elements = m_grid.elements();
            const auto& nodes = m_grid.nodes();
            const auto& v0 = nodes[elements[index][0]];
            const auto& v1 = nodes[elements[index][1]];
            const auto& v2 = nodes[elements[index][2]];
            const auto& v3 = nodes[elements[index][3]];
            res.normal = basicshape::tetrahedron::closestNormalToPoint(v0, v1, v2, v3, hit_pos);
        }
        return res;
    }

    const EnergyScore& energyScored(std::size_t index=0) const
    {
        return m_energyScore.at(index);
    }
    const DoseScore& doseScored(std::size_t index=0) const
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
        auto elementVolumes = m_grid.volumes();
        std::vector<double> elementDensity(elementVolumes.size());
        std::transform(std::execution::par_unseq, m_collectionIdx.cbegin(), m_collectionIdx.cend(), elementDensity.begin(), [&](const auto collIdx) {
            return m_collectionDensities[collIdx];
        });

        for (std::size_t i = 0; i < elementVolumes.size(); ++i) {
            m_doseScore[i].addScoredEnergy(m_energyScore[i], elementVolumes[i], elementDensity[i], factor);
        }
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

    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

    void rotate(const std::array<double, 3>& axis, double angle)
    {
        m_grid.rotate(axis, angle);
        std::transform(std::execution::par_unseq, m_data.nodes.cbegin(), m_data.nodes.cend(), m_data.nodes.begin(), [angle, &axis](const auto& v) {
            return vectormath::rotate(v, axis, angle);
        });
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
        for (const auto& mat : m_collectionMaterials) {
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
        for (std::size_t mIdx = 0; mIdx < m_collectionMaterials.size(); ++mIdx) {
            const auto& mat = m_collectionMaterials[mIdx];
            const auto d = m_collectionDensities[mIdx];
            for (std::size_t i = 0; i < energy.size(); ++i) {
                const auto aval = mat.attenuationValues(energy[i]);
                att[i] = std::max(aval.sum() * d, att[i]);
            }
        }

        m_woodcockStepTableLin.resize(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.cbegin(), m_woodcockStepTableLin.begin(), [](const auto e, const auto a) {
            return std::make_pair(e, a);
        });
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
            const auto tetIdx = m_grid.pointInside(p.pos);

            if (tetIdx) { // is interaction virtual?
                const auto collIdx = m_collectionIdx[tetIdx.value()];
                const auto attenuation = m_collectionMaterials[collIdx].attenuationValues(p.energy);
                const auto attSum = attenuation.sum() * m_collectionDensities[collIdx];
                if (state.randomUniform() < attSum * attMaxInv) {
                    // we have a real interaction
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(attenuation, p, m_collectionMaterials[collIdx], state);
                    m_energyScore[tetIdx].scoreEnergy(intRes.energyImparted);
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
        auto tetIdx_cand = m_grid.pointInside(p.pos);

        const auto& elements = m_grid.elements();
        const auto& nodes = m_grid.nodes();

        while (tetIdx_cand) {
            const auto tetIdx = tetIdx_cand.value();
            const auto collIdx = m_collectionIdx[tetIdx];
            const auto& material = m_collectionMaterials[collIdx];
            const auto& density = m_collectionDensities[collIdx];

            const auto& v0 = nodes[elements[tetIdx][0]];
            const auto& v1 = nodes[elements[tetIdx][1]];
            const auto& v2 = nodes[elements[tetIdx][2]];
            const auto& v3 = nodes[elements[tetIdx][3]];

            const auto intLen = basicshape::tetrahedron::intersect(v0, v1, v2, v3, p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, density, p, material, state);
            m_energyScore[tetIdx].scoreEnergy(intRes.energyImparted);
            tetIdx_cand = intRes.particleAlive ? m_grid.pointInside(p.pos) : std::nullopt;
        }
    }

private:
    TetrahedalMeshGrid2 m_grid;
    std::vector<std::uint32_t> m_collectionIdx;
    std::vector<EnergyScore> m_energyScore;
    std::vector<DoseScore> m_doseScore;

    std::vector<std::pair<double, double>> m_woodcockStepTableLin;

    std::vector<double> m_collectionDensities;
    std::vector<Material<NMaterialShells>> m_collectionMaterials;
    std::vector<std::string> m_collectionNames;

    ParticleTracker m_tracker;
};
}
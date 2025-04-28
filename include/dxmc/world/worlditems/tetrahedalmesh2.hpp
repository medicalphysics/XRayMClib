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

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

struct TetrahedalMeshData {
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::array<std::size_t, 4>> elements;
    std::vector<std::size_t> collectionIndices; // same size as elements

    // specify collection/organ properties
    std::vector<double> collectionDensities;
    std::vector<std::map<std::size_t, double>> collectionMaterialComposition;
    std::vector<std::string> collectionNames;
};

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FLUENCESCORING = true>
class TetrahedalMesh2 {
public:
    TetrahedalMesh2()
    {
    }
    TetrahedalMesh2(const TetrahedalMeshData& data, std::array<std::size_t, 3> spacing = { 128, 128, 128 })
    {
        setData(data, spacing);
    }

    bool setData(const TetrahedalMeshData& data, std::array<std::size_t, 3> spacing)
    {
    }

    void translate(const std::array<double, 3>& vec)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += vec[i];
            m_aabb[i + 3] += vec[i];
        }
        std::for_each(std::execution::par_unseq(m_nodes.begin(), m_nodes.end(), [&vec](auto& n){
            n = vectormath::add(n, vec);
        });
    }

    std::array<double, 3> center() const
    {
        std::array<double, 3> c = {
            (m_aabb[0] + m_aabb[3]) * 0.5,
            (m_aabb[1] + m_aabb[4]) * 0.5,
            (m_aabb[2] + m_aabb[5]) * 0.5
        };
        return c;
    }
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    WorldIntersectionResult intersect(auto Particle& p) const
    { // not implemented
        return WorldIntersectionResult {};
    }

    VisualizationIntersectionResult intersectVisualization(auto Particle& p) const
    { // not implemented
        return VisualizationIntersectionResult {};
    }
    const EnergyScore& energyScored(std::size_t index) const
    {
        return m_energyScore.at(index);
    }
    const DoseScore& doseScored(std::size_t index) const
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
    }

    void transport(auto Particle& p, RandomState& state)
    { // not implemented
    }

protected:
    void computeAABB()
    {
        m_aabb = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (const auto& n : m_nodes) {
            for (std::size_t i = 0; i < 3; ++i) {
                m_aabb[i] = std::min(m_aabb[i], n[i]);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], n[i]);
            }
        }
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] -= GEOMETRIC_ERROR();
            m_aabb[i + 3] += GEOMETRIC_ERROR();
        }
    }

private:
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };

    std::vector<std::array<double, 3>> m_nodes;
    std::vector<std::array<std::size_t, 4>> m_elements;
    std::vector<std::uint16_t> m_collectionIndices;
    std::vector<EnergyScore> m_energyScore;
    std::vector<DoseScore> m_doseScore;

    std::vector<double> m_collectionDensities;
    std::vector<Material<NMaterialShells>> m_collectionMaterials;
    std::vector<std::string> m_collectionNames;
};
}
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
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <cmath>
#include <limits>
#include <optional>

namespace xraymc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class WorldBoxGrid {
public:
    WorldBoxGrid(const std::array<double, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : m_aabb(aabb)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
        setVoxelDimensions({ 1, 1, 1 });
    }

    WorldBoxGrid(double aabb_size, std::array<double, 3> pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
        setVoxelDimensions({ 1, 1, 1 });
    }

    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }
    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        m_material = material;
        setMaterialDensity(density);
    }

    void setMaterialDensity(double density) { m_materialDensity = std::abs(density); }

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

    void setVoxelDimensions(const std::array<std::uint_fast32_t, 3>& dim)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_voxelDim[i] = std::min(std::max(std::uint_fast32_t { 1 }, dim[i]), std::uint_fast32_t { 1000 });

        for (std::size_t i = 0; i < 3; ++i) {
            m_voxelSize[i] = (m_aabb[i + 3] - m_aabb[i]) / m_voxelDim[i];
            m_voxelSizeInv[i] = 1.0 / m_voxelSize[i];
        }
        const auto ndim = m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
        m_energyScored.resize(ndim);
        m_dose.resize(ndim);
    }

    std::uint_fast32_t totalNumberOfVoxels() const
    {
        return m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
    }

    const std::array<std::uint_fast32_t, 3>& voxelDimensions() const
    {
        return m_voxelDim;
    }

    const std::array<double, 3>& voxelSpacing() const
    {
        return m_voxelSize;
    }

    template <bool BOUNDSCHECK = true>
    std::uint_fast32_t gridIndex(const std::array<double, 3>& pos) const noexcept
    {
        if constexpr (BOUNDSCHECK) {
            const auto x = static_cast<std::uint_fast32_t>(std::clamp((pos[0] - m_aabb[0]) * m_voxelSizeInv[0], double { 0 }, static_cast<double>(m_voxelDim[0] - 1)));
            const auto y = static_cast<std::uint_fast32_t>(std::clamp((pos[1] - m_aabb[1]) * m_voxelSizeInv[1], double { 0 }, static_cast<double>(m_voxelDim[1] - 1)));
            const auto z = static_cast<std::uint_fast32_t>(std::clamp((pos[2] - m_aabb[2]) * m_voxelSizeInv[2], double { 0 }, static_cast<double>(m_voxelDim[2] - 1)));
            return x + (y + z * m_voxelDim[1]) * m_voxelDim[0];
        } else {
            const auto x = static_cast<std::uint_fast32_t>((pos[0] - m_aabb[0]) * m_voxelSizeInv[0]);
            const auto y = static_cast<std::uint_fast32_t>((pos[1] - m_aabb[1]) * m_voxelSizeInv[1]);
            const auto z = static_cast<std::uint_fast32_t>((pos[2] - m_aabb[2]) * m_voxelSizeInv[2]);
            return x + (y + z * m_voxelDim[1]) * m_voxelDim[0];
        }
    }

    void translate(const std::array<double, 3>& dist) noexcept
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<double, 3> center() const noexcept
    {
        std::array<double, 3> c {
            (m_aabb[0] + m_aabb[3]) * 0.5,
            (m_aabb[1] + m_aabb[4]) * 0.5,
            (m_aabb[2] + m_aabb[5]) * 0.5,
        };
        return c;
    }

    const std::array<double, 6>& AABB() const noexcept
    {
        return m_aabb;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::AABB::intersect(p, m_aabb);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::AABB::template intersectVisualization<U>(p, m_aabb);
        if (inter.valid()) {
            auto p_int = p;
            p_int.translate(inter.intersection);
            const auto ind = gridIndex(p_int.pos);
            inter.value = m_dose[ind].dose();
        }
        return inter;
    }

    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::AABB::pointInside(p.pos, m_aabb);
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                const auto scoreIdx = gridIndex<false>(p.pos);
                m_energyScored[scoreIdx].scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored.at(index);
    }

    void clearEnergyScored()
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = m_voxelSize[0] * m_voxelSize[1] * m_voxelSize[2];
        for (std::size_t i = 0; i < m_energyScored.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScored[i], volume, m_materialDensity, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose.at(index);
    }

    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

protected:
    void correctAABB()
    {
        auto test = [](const auto& aabb) -> bool {
            bool ok = true;
            for (std::size_t i = 0; i < 3; ++i)
                ok = ok && aabb[i] < aabb[i + 3];
            return ok;
        };
        if (!test(m_aabb)) {
            for (std::size_t i = 0; i < 3; ++i) {
                if (std::abs(m_aabb[i] - m_aabb[i + 3]) < GEOMETRIC_ERROR()) {
                    if (m_aabb[i] > m_aabb[i + 3]) {
                        std::swap(m_aabb[i], m_aabb[i + 3]);
                    }
                    m_aabb[i] -= GEOMETRIC_ERROR();
                    m_aabb[i + 3] += GEOMETRIC_ERROR();
                } else if (m_aabb[i] > m_aabb[i + 3]) {
                    std::swap(m_aabb[i], m_aabb[i + 3]);
                }
            }
        }
    }

private:
    std::array<double, 6> m_aabb;
    double m_materialDensity = 1;
    std::array<double, 3> m_voxelSize = { 1, 1, 1 };
    std::array<double, 3> m_voxelSizeInv = { 1, 1, 1 };
    std::array<std::uint_fast32_t, 3> m_voxelDim = { 1, 1, 1 };
    Material<NMaterialShells> m_material;
    std::vector<EnergyScore> m_energyScored;
    std::vector<DoseScore> m_dose;
};
}

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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FORCE_INTERACTIONS = false>
class WorldBox {
public:
    WorldBox(const std::array<double, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : m_aabb(aabb)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
    }

    WorldBox(double aabb_size, std::array<double, 3> pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        correctAABB();
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
    }

    void setAABB(const std::array<double, 6>& aabb)
    {
        m_aabb = aabb;
        correctAABB();
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

    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
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
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        if constexpr (FORCE_INTERACTIONS) {
            transportForced(p, state);
        } else {
            transportRandom(p, state);
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
        const auto [l, h] = vectormath::splice(m_aabb);
        const auto sides = vectormath::subtract(h, l);
        const auto volume = std::reduce(sides.cbegin(), sides.cend(), double { 1 }, std::multiplies<>());

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

protected:
    void transportForced(ParticleType auto& p, RandomState& state) noexcept
    {

        bool cont = basicshape::AABB::pointInside(p.pos, m_aabb);
        while (cont) {
            const auto intLen = intersect(p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, m_materialDensity, p, m_material, state);
            m_energyScored.scoreEnergy(intRes.energyImparted);
            cont = intRes.particleAlive && basicshape::AABB::pointInside(p.pos, m_aabb);
        }
    }
    void transportRandom(ParticleType auto& p, RandomState& state) noexcept
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
                m_energyScored.scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

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
                if (m_aabb[i] == m_aabb[i + 3]) {
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
    Material<NMaterialShells> m_material;
    double m_materialDensity = 1;
    EnergyScore m_energyScored;
    DoseScore m_dose;
};

}

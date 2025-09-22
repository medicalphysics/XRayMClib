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

#include "xraymc/floating.hpp"
#include "xraymc/interactions.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TG195World42 {
protected:
    static bool insideChild(const Particle& p, const basicshape::cylinder::Cylinder& child, const std::array<double, 6>& aabb)
    {
        if (basicshape::AABB::pointInside(p.pos, aabb)) {
            return basicshape::cylinder::pointInside(p.pos, child);
        }
        return false;
    }

public:
    TG195World42(double radius = 16, double height = 10, const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_cylinder.center = pos;
        m_cylinder.direction = { 0, 0, 1 };
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");

        m_centerChild.center = m_cylinder.center;
        m_centerChild.direction = m_cylinder.direction;
        m_centerChild.radius = 0.5;
        m_centerChild.half_height = 5;
        m_periferyChild.center = { m_cylinder.center[0] - (m_cylinder.radius - 1), m_cylinder.center[1], m_cylinder.center[2] };
        m_periferyChild.direction = m_cylinder.direction;
        m_periferyChild.radius = 0.5;
        m_periferyChild.half_height = 5;
        for (std::size_t i = 0; i < 3; ++i) {
            m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.radius;
            m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.radius;
            m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.radius;
            m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.radius;
            if (i == 2) {
                m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.half_height;
                m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.half_height;
                m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.half_height;
                m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.half_height;
            }
        }
    }

    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(double density) { m_materialDensity = density; }

    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        m_centerChild.center = vectormath::add(m_centerChild.center, dist);
        m_periferyChild.center = vectormath::add(m_periferyChild.center, dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_periferyChild_aabb[i] += dist[i];
            m_periferyChild_aabb[i + 3] += dist[i];
            m_centerChild_aabb[i] += dist[i];
            m_centerChild_aabb[i + 3] += dist[i];
        }
    }

    std::array<double, 3> center() const
    {
        return m_cylinder.center;
    }

    std::array<double, 6> AABB() const
    {
        return xraymc::basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    WorldIntersectionResult intersect(const Particle& p) const noexcept
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const Particle& p) const noexcept
    {
        return basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
    }

    void transport(Particle& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1.0 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            const auto stepLen = -std::log(state.randomUniform()) * attSumInv;

            const auto intLen = intersect(p);
            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                if (insideChild(p, m_centerChild, m_centerChild_aabb)) {
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    m_centerChild_energyScored.scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else if (insideChild(p, m_periferyChild, m_periferyChild_aabb)) {
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    m_periferyChild_energyScored.scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else {
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                }

            } else {
                // transport to border
                p.border_translate(intLen.intersection);
                cont = false;
            }
        }
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_centerChild_energyScored.clear();
        m_periferyChild_energyScored.clear();
    }
    const EnergyScore& energyScoredCenterCylinder() const
    {
        return m_centerChild_energyScored;
    }
    const EnergyScore& energyScoredPeriferyCylinder() const
    {
        return m_periferyChild_energyScored;
    }
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_centerChild_energyScored;
        else if (index == 1)
            return m_periferyChild_energyScored;
        return m_energyScored;
    }
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto c_vol = m_centerChild.volume();
        const auto p_vol = m_periferyChild.volume();
        const auto vol = m_cylinder.volume() - c_vol - p_vol;
        m_doseScored.addScoredEnergy(m_energyScored, vol, m_materialDensity, calibration_factor);
        m_centerChild_doseScored.addScoredEnergy(m_centerChild_energyScored, c_vol, m_materialDensity, calibration_factor);
        m_periferyChild_doseScored.addScoredEnergy(m_periferyChild_energyScored, p_vol, m_materialDensity, calibration_factor);
    };
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_centerChild_doseScored;
        else if (index == 1)
            return m_periferyChild_doseScored;
        return m_doseScored;
    }
    void clearDoseScored()
    {
        m_doseScored.clear();
        m_centerChild_doseScored.clear();
        m_periferyChild_doseScored.clear();
    }

private:
    basicshape::cylinder::Cylinder m_cylinder;
    double m_materialDensity = 1;
    Material<NMaterialShells> m_material;
    EnergyScore m_energyScored;
    EnergyScore m_periferyChild_energyScored;
    EnergyScore m_centerChild_energyScored;
    DoseScore m_doseScored;
    DoseScore m_periferyChild_doseScored;
    DoseScore m_centerChild_doseScored;
    basicshape::cylinder::Cylinder m_centerChild;
    basicshape::cylinder::Cylinder m_periferyChild;
    std::array<double, 6> m_centerChild_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<double, 6> m_periferyChild_aabb = { 0, 0, 0, 0, 0, 0 };
};
}

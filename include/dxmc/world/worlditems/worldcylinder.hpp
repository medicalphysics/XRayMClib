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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class WorldCylinder {
public:
    WorldCylinder(double radius = 16, double height = 10, const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 })
        : m_material(Material<NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
    {
        m_cylinder.radius = std::abs(radius);
        m_cylinder.half_height = std::abs(height) / 2;
        m_cylinder.center = center;
        m_cylinder.direction = vectormath::normalized(dir);
        m_materialDensity = NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");
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

    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        updateAABB();
    }

    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    void setCenter(const std::array<double, 3>& c)
    {
        m_cylinder.center = c;
        updateAABB();
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void setRadius(double r)
    {
        m_cylinder.radius = std::abs(r);
        updateAABB();
    }
    double radius() const
    {
        return m_cylinder.radius;
    }

    void setHeight(double h)
    {
        m_cylinder.half_height = std::abs(h / 2);
        updateAABB();
    }
    double height() const
    {
        return m_cylinder.half_height * 2;
    }

    void setDirection(const std::array<double, 3>& dir)
    {
        m_cylinder.direction = vectormath::normalized(dir);
        updateAABB();
    }

    const std::array<double, 3>& direction() const
    {
        return m_cylinder.direction;
    }

    WorldIntersectionResult
    intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
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
            const auto intLen = intersect(p);

            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                m_energyScored.scoreEnergy(intRes.energyImparted);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
            } else {
                // transport to border
                p.border_translate(intLen.intersection);
                cont = false;
            }
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
        m_dose.addScoredEnergy(m_energyScored, m_cylinder.volume(), m_materialDensity, calibration_factor);
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
    void updateAABB()
    {
        m_aabb = basicshape::cylinder::cylinderAABB(m_cylinder);
    }

private:
    basicshape::cylinder::Cylinder m_cylinder;
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    double m_materialDensity = 1;
    Material<NMaterialShells> m_material;
    EnergyScore m_energyScored;
    DoseScore m_dose;
};

}

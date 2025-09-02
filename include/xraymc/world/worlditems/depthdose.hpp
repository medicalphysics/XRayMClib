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
#include "xraymc/particle.hpp"
#include "xraymc/particletracker.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class DepthDose {
public:
    DepthDose(double radius = 16, double height = 10, std::size_t resolution = 100, const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_cylinder.center = pos;
        m_cylinder.direction = dir;
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;
        updateAABB();

        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        m_energyScored.resize(resolution);
        m_dose.resize(resolution);
    }

    double length() const { return m_cylinder.half_height * 2; }
    double radius() const { return m_cylinder.radius; }

    std::size_t resolution() const
    {
        return m_energyScored.size();
    }

    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
    }

    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        setMaterial(material);
        setMaterialDensity(density);
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

    const Material<NMaterialShells>& material() const
    {
        return m_material;
    }

    double density() const { return m_materialDensity; }

    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        updateAABB();
    }

    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        auto inter = basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
        if (inter.valid()) {
            auto p_int = p;
            p_int.translate(inter.intersection);
            const auto dose_ind_f = (p_int.pos[2] - (m_cylinder.center[2] - m_cylinder.half_height)) * m_energyScored.size() / (m_cylinder.half_height * 2);
            const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_energyScored.size() - 1);
            inter.value = m_dose[ind].dose();
        }
        return inter;
    }

    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }
        if constexpr (FORCEINTERACTIONS)
            transportForced(p, state);
        else
            transportRandom(p, state);
    }

    const std::vector<std::pair<double, EnergyScore>> depthEnergyScored() const
    {
        std::vector<std::pair<double, EnergyScore>> depth;
        depth.reserve(m_energyScored.size());

        const auto step = (2 * m_cylinder.half_height) / m_energyScored.size();
        const auto start = m_cylinder.center[2] - m_cylinder.half_height + step / 2;
        for (std::size_t i = 0; i < m_energyScored.size(); ++i) {
            depth.push_back(std::make_pair(start + step * i, m_energyScored[i]));
        }
        return depth;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored[index];
    }

    void clearEnergyScored()
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto totalVolume = m_cylinder.volume();
        const auto partVolume = totalVolume / m_dose.size();
        for (std::size_t i = 0; i < m_dose.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScored[i], partVolume, m_materialDensity, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose[index];
    }

    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
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
    void updateAABB()
    {
        m_aabb = basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    template <bool BOUNDS_CHECK = true>
    std::size_t cylinderIndex(const std::array<double, 3>& pos) const
    {
        const auto cstart = vectormath::subtract(m_cylinder.center, vectormath::scale(m_cylinder.direction, m_cylinder.half_height));
        const auto cdelta = vectormath::subtract(pos, cstart);
        const auto dz = vectormath::dot(cdelta, m_cylinder.direction);
        if constexpr (BOUNDS_CHECK) {
            const auto ind_f = std::clamp(dz * m_energyScored.size() / (2 * m_cylinder.half_height), 0.0, static_cast<double>(m_energyScored.size() - 1));
            return static_cast<std::size_t>(ind_f);
        } else {
            const auto ind_f = dz * m_energyScored.size() / (2 * m_cylinder.half_height);
            return static_cast<std::size_t>(ind_f);
        }
    }

    void transportRandom(ParticleType auto& p, RandomState& state)
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
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
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;

                const auto ind = cylinderIndex<true>(p.pos);
                m_energyScored[ind].scoreEnergy(intRes.energyImparted);
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    void transportForced(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        while (cont) {
            const auto inter = intersect(p);
            auto intLen = inter.intersection; // this must be valid
            const auto ind = cylinderIndex<true>(p.pos);

            // do we hit a seperating plane?
            {
                const auto p_cyl_proj = vectormath::dot(m_cylinder.direction, p.dir);
                if (std::abs(p_cyl_proj) > GEOMETRIC_ERROR<>()) {
                    const auto plane_point = vectormath::add(m_cylinder.center, vectormath::scale(m_cylinder.direction, m_cylinder.half_height * (2.0 * (ind + (p_cyl_proj > 0.0 ? 1 : 0)) / resolution() - 1.0)));
                    const auto t_plane = vectormath::dot(vectormath::subtract(plane_point, p.pos), m_cylinder.direction) / p_cyl_proj;
                    intLen = std::min(intLen, t_plane);
                }
            }
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, m_materialDensity, p, m_material, state);

            m_energyScored[ind].scoreEnergy(intRes.energyImparted);
            cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
        }
    }

private:
    std::array<double, 6> m_aabb;
    xraymc::basicshape::cylinder::Cylinder m_cylinder;
    double m_materialDensity = 1;
    Material<NMaterialShells> m_material;
    ParticleTracker m_tracker;
    std::vector<EnergyScore> m_energyScored;
    std::vector<DoseScore> m_dose;
};

}

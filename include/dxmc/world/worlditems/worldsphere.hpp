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
#include "dxmc/particletracker.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/sphere.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class WorldSphere {
public:
    WorldSphere(double radius = 16, const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_radius(std::abs(radius))
        , m_center(pos)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
    }

    void setRadius(double r)
    {
        m_radius = std::abs(r);
    }

    auto radius() const { return m_radius; }

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

    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
    }

    std::array<double, 3> center() const
    {
        return m_center;
    }

    std::array<double, 6> AABB() const
    {
        std::array<double, 6> aabb {
            m_center[0] - m_radius,
            m_center[1] - m_radius,
            m_center[2] - m_radius,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_radius
        };
        return aabb;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::sphere::intersect(p, m_center, m_radius);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::sphere::template intersectVisualization<U>(p, m_center, m_radius);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        m_meanEnergyAccumulator.addEnergy(p.energy);

        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }
        if constexpr (FORCEINTERACTIONS)
            transportForced(p, state);
        else
            transportRandom(p, state);
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
        const auto volume = (4 * std::numbers::pi_v<double> * m_radius * m_radius * m_radius) / 3;
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

    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

protected:
    void transportRandom(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::sphere::pointInside(p.pos, m_center, m_radius);
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
            const auto intLen = intersect(p).intersection; // this must be valid

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

    void transportForced(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::sphere::pointInside(p.pos, m_center, m_radius);
        while (cont) {
            const auto intLen = intersect(p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, m_materialDensity, p, m_material, state);
            m_energyScored.scoreEnergy(intRes.energyImparted);
            cont = intRes.particleAlive && basicshape::sphere::pointInside(p.pos, m_center, m_radius);
        }
    }

private:
    struct MeanEnergyAccumulator {
        double energysum = 0;
        std::uint64_t n_sum = 0;
        void addEnergy(double e)
        {
            auto eref = std::atomic_ref(energysum);
            eref.fetch_add(e, std::memory_order_relaxed);
            auto nref = std::atomic_ref(n_sum);
            nref.fetch_add(n_sum, std::memory_order_relaxed);
        }
        double meanEnergy() const
        {
            return n_sum > 0 ? energysum / n_sum : 0.0;
        }
    };
    double m_radius = 0;
    std::array<double, 3> m_center;
    Material<NMaterialShells> m_material;
    double m_materialDensity = 1;
    EnergyScore m_energyScored;
    DoseScore m_dose;
    MeanEnergyAccumulator m_meanEnergyAccumulator;
    ParticleTracker m_tracker;
};
}

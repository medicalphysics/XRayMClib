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
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinderz.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TG195World3Breast {
public:
    TG195World3Breast()
        : m_skin_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_tissue_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        translateBox(m_dose_boxes[0], { double { 5 }, double { 5 }, double { 0 } });
        translateBox(m_dose_boxes[1], { double { 2 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[2], { double { 5 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[3], { double { 8 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[4], { double { 5 }, double { -5 }, double { 0 } });
        translateBox(m_dose_boxes[5], { double { 5 }, double { 0 }, double { -1.5f } });
        translateBox(m_dose_boxes[6], { double { 5 }, double { 0 }, double { 1.5f } });
    }

    void setTissueMaterial(const Material<NMaterialShells>& material, double dens)
    {
        m_tissue_material = material;
        m_tissue_density = std::abs(dens);
    }

    void setSkinMaterial(const Material<NMaterialShells>& material, double dens)
    {
        m_skin_material = material;
        m_skin_density = std::abs(dens);
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
        for (auto& b : m_dose_boxes) {
            for (std::size_t i = 0; i < 3; ++i) {
                b.aabb[i] += dist[i];
                b.aabb[i + 3] += dist[i];
            }
        }
    }

    std::array<double, 3> center() const
    {
        return m_center;
    }

    std::array<double, 6> AABB() const
    {
        std::array aabb {
            m_center[0],
            m_center[1] - m_radius,
            m_center[2] - m_halfHeight,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_halfHeight
        };
        return aabb;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_skin_energy.clear();
        for (auto& b : m_dose_boxes) {
            b.energyScored.clear();
        }
    }
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index < m_dose_boxes.size())
            return m_dose_boxes[index].energyScored;
        else if (index == m_dose_boxes.size())
            return m_skin_energy;
        return m_energyScored;
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        auto skin_volume = std::numbers::pi_v<double> * (m_halfHeight - m_skin_thick) * 2 * (m_radius * m_radius - (m_radius - m_skin_thick) * (m_radius - m_skin_thick));
        skin_volume += std::numbers::pi_v<double> * 2 * m_skin_thick * m_radius * m_radius;
        // only half cylinder
        skin_volume /= 2;

        std::vector<double> box_volumes(m_dose_boxes.size());
        std::transform(m_dose_boxes.cbegin(), m_dose_boxes.cend(), box_volumes.begin(), [](const auto& box) {
            double vol = 1;
            for (std::size_t i = 0; i < 3; ++i) {
                vol *= box.aabb[i + 3] - box.aabb[i];
            }
            return vol;
        });

        const auto box_volume_total = std::reduce(box_volumes.cbegin(), box_volumes.cend());
        const auto total_volume = std::numbers::pi_v<double> * 2 * m_halfHeight * m_radius * m_radius / 2 - skin_volume - box_volume_total;

        m_doseScored.addScoredEnergy(m_energyScored, total_volume, m_tissue_density, calibration_factor);
        m_skin_dose.addScoredEnergy(m_skin_energy, skin_volume, m_skin_density, calibration_factor);
        for (std::size_t i = 0; i < m_dose_boxes.size(); ++i) {
            auto& bd = m_dose_boxes[i].doseScored;
            auto& be = m_dose_boxes[i].energyScored;
            bd.addScoredEnergy(be, box_volumes[i], m_tissue_density, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index < m_dose_boxes.size())
            return m_dose_boxes[index].doseScored;
        else if (index == m_dose_boxes.size())
            return m_skin_dose;
        return m_doseScored;
    }

    void clearDoseScored()
    {
        m_doseScored.clear();
        m_skin_dose.clear();
        for (auto& b : m_dose_boxes) {
            b.doseScored.clear();
        }
    }

    WorldIntersectionResult intersect(const Particle& p) const noexcept
    {
        return intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const Particle& p) const noexcept
    {
        const auto aabb = AABB();
        auto cyl = basicshape::cylinderZ::template intersectVisualization<U>(p, m_center, m_radius, m_halfHeight);
        const auto box = basicshape::AABB::template intersectVisualization<U>(p, aabb);
        if (cyl.valid() && box.valid()) {
            if (basicshape::AABB::pointInside(p.pos, aabb)) {
                cyl.intersection = std::min(cyl.intersection, box.intersection);
                cyl.rayOriginIsInsideItem = true;
            } else {
                if (cyl.intersection < box.intersection) {
                    cyl.normal = { -1, 0, 0 };
                    cyl.intersection = box.intersection;
                }
                cyl.rayOriginIsInsideItem = false;
            }
        } else {
            cyl.intersectionValid = false;
        }

        return cyl;
    }

    void transport(Particle& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinderZ::pointInside(p.pos, m_center, m_radius, m_halfHeight) && basicshape::AABB::pointInside(p.pos, AABB());
        while (cont) {
            const auto intBreast = intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
            if (intBreast.valid()) {
                const auto intTissue = intersectHalfCylindar(p, m_center, m_radius - m_skin_thick, m_halfHeight - m_skin_thick);
                if (intTissue.valid()) {
                    if (intTissue.rayOriginIsInsideItem) {
                        // only inside tissue
                        const auto att = m_tissue_material.attenuationValues(p.energy);
                        const auto attSumInv = 1 / (att.sum() * m_tissue_density);
                        const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
                        if (stepLen < intTissue.intersection) {
                            p.translate(stepLen);
                            const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_tissue_material, state);
                            cont = intRes.particleAlive;
                            scoreEnergyImparted(p, intRes.energyImparted);
                        } else {
                            p.border_translate(intTissue.intersection);
                            cont = basicshape::cylinderZ::pointInside(p.pos, m_center, m_radius, m_halfHeight) && basicshape::AABB::pointInside(p.pos, AABB());
                        }
                    } else {
                        // starts in skin and goes to tissue
                        const auto att = m_skin_material.attenuationValues(p.energy);
                        const auto attSumInv = 1 / (att.sum() * m_skin_density);
                        const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
                        if (stepLen < intTissue.intersection) {
                            p.translate(stepLen);
                            const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_skin_material, state);
                            cont = intRes.particleAlive;
                            m_skin_energy.scoreEnergy(intRes.energyImparted);
                        } else {
                            p.border_translate(intTissue.intersection);
                            cont = basicshape::cylinderZ::pointInside(p.pos, m_center, m_radius, m_halfHeight) && basicshape::AABB::pointInside(p.pos, AABB());
                        }
                    }
                } else {
                    // only intersects skin
                    const auto att = m_skin_material.attenuationValues(p.energy);
                    const auto attSumInv = 1 / (att.sum() * m_skin_density);
                    const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
                    if (stepLen < intBreast.intersection) {
                        p.translate(stepLen);
                        const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_skin_material, state);
                        cont = intRes.particleAlive;
                        m_skin_energy.scoreEnergy(intRes.energyImparted);
                    } else {
                        p.border_translate(intBreast.intersection);
                        cont = false;
                    }
                }
            } else {
                cont = false;
            }
        }
    }

protected:
    struct ScoreBoxChild {
        std::array<double, 6> aabb = { -1, -1, -.5, 1, 1, .5 };
        EnergyScore energyScored;
        DoseScore doseScored;
    };
    static void translateBox(ScoreBoxChild& box, const std::array<double, 3>& vec)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            box.aabb[i] += vec[i];
            box.aabb[i + 3] += vec[i];
        }
    }

    void scoreEnergyImparted(const Particle& p, double energy)
    {
        bool hit = false;
        for (auto& b : m_dose_boxes) {
            if (!hit && basicshape::AABB::pointInside(p.pos, b.aabb)) {
                b.energyScored.scoreEnergy(energy);
                hit = true;
            }
        }
        m_energyScored.scoreEnergy(energy);
    }

    static WorldIntersectionResult intersectHalfCylindar(const Particle& p, const std::array<double, 3>& center, const double radius, const double halfHeight)
    {
        const std::array aabb = {
            center[0],
            center[1] - radius,
            center[2] - halfHeight,
            center[0] + radius,
            center[1] + radius,
            center[2] + halfHeight
        };
        auto box = basicshape::AABB::intersect(p, aabb);
        if (box.valid()) {
            const auto cyl = basicshape::cylinderZ::intersect(p, center, radius, halfHeight);
            if (cyl.valid()) {
                box.rayOriginIsInsideItem = box.rayOriginIsInsideItem && cyl.rayOriginIsInsideItem;
                if (box.rayOriginIsInsideItem) {
                    box.intersection = std::min(cyl.intersection, box.intersection);
                } else {
                    box.intersection = std::max(cyl.intersection, box.intersection);
                }
            } else {
                box.intersectionValid = false;
                box.intersection = 0;
                box.rayOriginIsInsideItem = false;
            }
        }
        return box;
    }

private:
    double m_radius = 10;
    double m_skin_thick = 0.2;
    double m_halfHeight = 2.5;
    std::array<double, 3> m_center = { 0, 0, 0 };
    double m_skin_density = 1;
    double m_tissue_density = 1;
    Material<NMaterialShells> m_skin_material;
    Material<NMaterialShells> m_tissue_material;
    EnergyScore m_energyScored;
    DoseScore m_doseScored;
    EnergyScore m_skin_energy;
    DoseScore m_skin_dose;
    std::array<ScoreBoxChild, 7> m_dose_boxes;
};
}

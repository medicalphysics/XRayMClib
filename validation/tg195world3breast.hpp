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

#include "xraymc/floating.hpp"
#include "xraymc/interactions.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinderz.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

/**
 * @brief TG195 Case 3 breast phantom: a half-cylinder with a skin layer and embedded dose scoring boxes.
 *
 * Models the compressed breast geometry from AAPM TG-195. The phantom is a half-cylinder
 * (flat face toward -X, curved face toward +X) with a uniform skin shell of configurable
 * thickness surrounding glandular/adipose tissue. Seven 2x2x1 cm dose-scoring boxes are
 * pre-positioned inside the phantom at locations defined by TG-195.
 *
 * Dose scoring index convention (used by energyScored() and doseScored()):
 *   - 0–6  : individual dose boxes
 *   - 7    : skin shell
 *   - 8    : remaining tissue volume
 *
 * @tparam NMaterialShells  Number of electron shells used in the material model (default 5).
 * @tparam LOWENERGYCORRECTION  Low-energy correction mode passed to the interaction kernel (default 2).
 */
template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TG195World3Breast {
public:
    /**
     * @brief Constructs the phantom with default air materials and positions the seven TG-195 dose boxes.
     *
     * Materials default to dry air so that setTissueMaterial() / setSkinMaterial() must be called
     * before a meaningful simulation is run.
     */
    TG195World3Breast()
        : m_skin_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_tissue_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        // TG-195 Table 3 dose-box centre positions (x, y, z) in cm relative to phantom centre.
        translateBox(m_dose_boxes[0], { double { 5 }, double { 5 }, double { 0 } });
        translateBox(m_dose_boxes[1], { double { 2 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[2], { double { 5 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[3], { double { 8 }, double { 0 }, double { 0 } });
        translateBox(m_dose_boxes[4], { double { 5 }, double { -5 }, double { 0 } });
        translateBox(m_dose_boxes[5], { double { 5 }, double { 0 }, double { -1.5f } });
        translateBox(m_dose_boxes[6], { double { 5 }, double { 0 }, double { 1.5f } });
    }

    /**
     * @brief Sets the material and mass density for the inner tissue volume.
     * @param material  Material definition (e.g. glandular tissue or adipose).
     * @param dens      Mass density in g/cm³; the absolute value is used.
     */
    void setTissueMaterial(const Material<NMaterialShells>& material, double dens)
    {
        m_tissue_material = material;
        m_tissue_density = std::abs(dens);
    }

    /**
     * @brief Sets the material and mass density for the skin shell.
     * @param material  Material definition for the skin layer.
     * @param dens      Mass density in g/cm³; the absolute value is used.
     */
    void setSkinMaterial(const Material<NMaterialShells>& material, double dens)
    {
        m_skin_material = material;
        m_skin_density = std::abs(dens);
    }

    /**
     * @brief Translates the entire phantom (centre + all dose boxes) by a displacement vector.
     * @param dist  Translation vector {dx, dy, dz} in cm.
     */
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

    /**
     * @brief Returns the phantom centre position.
     * @return Centre coordinates {x, y, z} in cm.
     */
    std::array<double, 3> center() const
    {
        return m_center;
    }

    /**
     * @brief Returns the axis-aligned bounding box of the full phantom.
     *
     * The flat face lies at x = center[0]; the curved face extends to x = center[0] + radius.
     *
     * @return AABB as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     */
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

    /**
     * @brief Resets all energy-scored accumulators (dose boxes, skin, and bulk tissue).
     */
    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_skin_energy.clear();
        for (auto& b : m_dose_boxes) {
            b.energyScored.clear();
        }
    }

    /**
     * @brief Returns the energy scorer for a given scoring region.
     * @param index  0–6: individual dose box; 7: skin shell; anything else: bulk tissue.
     * @return Const reference to the corresponding EnergyScore.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index < m_dose_boxes.size())
            return m_dose_boxes[index].energyScored;
        else if (index == m_dose_boxes.size())
            return m_skin_energy;
        return m_energyScored;
    }

    /**
     * @brief Converts accumulated energy scores to absorbed dose for all scoring regions.
     *
     * Volumes are computed analytically. The skin-shell volume accounts for the curved side
     * and the two flat end caps; the bulk-tissue volume subtracts the skin shell and all dose
     * boxes from the total half-cylinder volume. Each region's dose accumulator is updated via
     * DoseScore::addScoredEnergy().
     *
     * @param calibration_factor  Optional scale factor applied to all dose values (default 1).
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        // Curved side walls of the skin shell (half-cylinder annulus, both z-faces).
        auto skin_volume = std::numbers::pi_v<double> * (m_halfHeight - m_skin_thick) * 2 * (m_radius * m_radius - (m_radius - m_skin_thick) * (m_radius - m_skin_thick));
        // Add the two flat end-cap annuli of the skin shell.
        skin_volume += std::numbers::pi_v<double> * 2 * m_skin_thick * m_radius * m_radius;
        // The phantom is a half-cylinder, so halve the full-cylinder skin volume.
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
        // Bulk tissue = total half-cylinder - skin shell - all dose boxes.
        const auto total_volume = std::numbers::pi_v<double> * 2 * m_halfHeight * m_radius * m_radius / 2 - skin_volume - box_volume_total;

        m_doseScored.addScoredEnergy(m_energyScored, total_volume, m_tissue_density, calibration_factor);
        m_skin_dose.addScoredEnergy(m_skin_energy, skin_volume, m_skin_density, calibration_factor);
        for (std::size_t i = 0; i < m_dose_boxes.size(); ++i) {
            auto& bd = m_dose_boxes[i].doseScored;
            auto& be = m_dose_boxes[i].energyScored;
            bd.addScoredEnergy(be, box_volumes[i], m_tissue_density, calibration_factor);
        }
    }

    /**
     * @brief Returns the dose scorer for a given scoring region.
     * @param index  0–6: individual dose box; 7: skin shell; anything else: bulk tissue.
     * @return Const reference to the corresponding DoseScore.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index < m_dose_boxes.size())
            return m_dose_boxes[index].doseScored;
        else if (index == m_dose_boxes.size())
            return m_skin_dose;
        return m_doseScored;
    }

    /**
     * @brief Resets all dose accumulators (dose boxes, skin, and bulk tissue).
     */
    void clearDoseScored()
    {
        m_doseScored.clear();
        m_skin_dose.clear();
        for (auto& b : m_dose_boxes) {
            b.doseScored.clear();
        }
    }

    /**
     * @brief Computes the intersection of a particle ray with the half-cylinder phantom boundary.
     * @param p  Particle whose position and direction define the ray.
     * @return   WorldIntersectionResult describing the next boundary crossing.
     */
    WorldIntersectionResult intersect(const Particle& p) const noexcept
    {
        return intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
    }

    /**
     * @brief Computes the intersection for visualization purposes, returning surface normal information.
     *
     * Clips the full cylinder against the flat bounding box to recover the half-cylinder geometry.
     * The flat face normal is forced to {-1, 0, 0} when the ray enters through it.
     *
     * @tparam U  Scalar type used for intersection distances.
     * @param p   Particle whose position and direction define the ray.
     * @return    VisualizationIntersectionResult with distance, normal, and inside-flag.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const Particle& p) const noexcept
    {
        const auto aabb = AABB();
        auto cyl = basicshape::cylinderZ::template intersectVisualization<U>(p, m_center, m_radius, m_halfHeight);
        const auto box = basicshape::AABB::template intersectVisualization<U>(p, aabb);
        if (cyl.valid() && box.valid()) {
            if (basicshape::AABB::pointInside(p.pos, aabb)) {
                // Ray originates inside the phantom: take the nearer of the two exits.
                cyl.intersection = std::min(cyl.intersection, box.intersection);
                cyl.rayOriginIsInsideItem = true;
            } else {
                // Ray originates outside: the flat face (AABB) may be the entry surface.
                if (cyl.intersection < box.intersection) {
                    cyl.normal = { -1, 0, 0 }; // flat face normal
                    cyl.intersection = box.intersection;
                }
                cyl.rayOriginIsInsideItem = false;
            }
        } else {
            cyl.intersectionValid = false;
        }

        return cyl;
    }

    /**
     * @brief Transports a particle through the phantom, sampling interactions until it exits or is absorbed.
     *
     * Handles three spatial cases at each step:
     *  1. Particle inside the inner tissue volume – tissue attenuation applied.
     *  2. Particle inside the skin shell heading toward tissue – skin attenuation applied.
     *  3. Particle traversing only the skin shell – skin attenuation applied; exits phantom on miss.
     *
     * Energy imparted to tissue is routed to the appropriate dose box (if inside one) and to
     * the bulk-tissue accumulator. Energy imparted to skin goes to the skin accumulator only.
     *
     * @param p      Particle to transport; modified in place (position, direction, energy, alive flag).
     * @param state  Random number generator state.
     */
    void transport(Particle& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinderZ::pointInside(p.pos, m_center, m_radius, m_halfHeight) && basicshape::AABB::pointInside(p.pos, AABB());
        while (cont) {
            const auto intBreast = intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
            if (intBreast.valid()) {
                const auto intTissue = intersectHalfCylindar(p, m_center, m_radius - m_skin_thick, m_halfHeight - m_skin_thick);
                if (intTissue.valid()) {
                    if (intTissue.rayOriginIsInsideItem) {
                        // Particle is inside the tissue core.
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
                        }
                    } else {
                        // Particle is in the skin shell, ray crosses into tissue.
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
                        }
                    }
                } else {
                    // Ray passes through skin only (e.g. grazing trajectory).
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
    /**
     * @brief A small AABB dose-scoring voxel embedded inside the phantom.
     *
     * Default extents are a 2×2×1 cm box centred at the origin; translateBox() is used
     * in the constructor to move each box to its TG-195 position.
     */
    struct ScoreBoxChild {
        std::array<double, 6> aabb = { -1, -1, -.5, 1, 1, .5 };
        EnergyScore energyScored;
        DoseScore doseScored;
    };

    /**
     * @brief Translates a ScoreBoxChild by a displacement vector.
     * @param box  Box to translate (modified in place).
     * @param vec  Translation {dx, dy, dz} in cm.
     */
    static void translateBox(ScoreBoxChild& box, const std::array<double, 3>& vec)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            box.aabb[i] += vec[i];
            box.aabb[i + 3] += vec[i];
        }
    }

    /**
     * @brief Scores energy imparted to tissue in the appropriate dose box and the bulk accumulator.
     *
     * Only the first matching dose box receives the energy to avoid double-counting if boxes overlap.
     * The bulk accumulator always receives the energy regardless of box membership.
     *
     * @param p       Particle at the interaction site (position is used for box lookup).
     * @param energy  Energy imparted at the interaction site in keV.
     */
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

    /**
     * @brief Computes the intersection of a particle ray with a half-cylinder (flat face at x = center[0]).
     *
     * The half-cylinder is formed by the intersection of a full Z-aligned cylinder and an AABB
     * that clips everything at x < center[0]. Intersection logic covers three cases:
     *  - Inside both:  take the nearer exit surface.
     *  - Inside AABB only (between flat face and cylinder surface): take nearer exit, mark as outside.
     *  - Outside both: take the farther entry surface (both must be entered to be inside).
     *
     * @param p           Particle defining the ray.
     * @param center      Centre of the cylinder base circle.
     * @param radius      Cylinder radius in cm.
     * @param halfHeight  Half-height of the cylinder along Z in cm.
     * @return            WorldIntersectionResult for the next boundary crossing.
     */
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
                if (box.rayOriginIsInsideItem && cyl.rayOriginIsInsideItem) {
                    // Inside both shapes: exit at whichever boundary is closer.
                    box.intersection = std::min(cyl.intersection, box.intersection);
                } else if (box.rayOriginIsInsideItem && !cyl.rayOriginIsInsideItem) {
                    // Between the flat face and the cylinder surface: not truly inside the half-cylinder.
                    box.intersection = std::min(cyl.intersection, box.intersection);
                    box.rayOriginIsInsideItem = false;
                } else {
                    // Outside: entry is at whichever boundary the ray hits last.
                    box.intersection = std::max(cyl.intersection, box.intersection);
                }
            } else {
                // Ray misses the cylinder entirely — no valid intersection.
                box.intersectionValid = false;
                box.intersection = 0;
                box.rayOriginIsInsideItem = false;
            }
        }
        return box;
    }

private:
    double m_radius = 10;           ///< Outer radius of the breast phantom in cm.
    double m_skin_thick = 0.2;      ///< Skin shell thickness in cm.
    double m_halfHeight = 2.5;      ///< Half-height (half compression thickness) in cm.
    std::array<double, 3> m_center = { 0, 0, 0 }; ///< Centre of the flat face of the half-cylinder.
    double m_skin_density = 1;      ///< Skin mass density in g/cm³.
    double m_tissue_density = 1;    ///< Tissue mass density in g/cm³.
    Material<NMaterialShells> m_skin_material;   ///< Skin material (attenuation data).
    Material<NMaterialShells> m_tissue_material; ///< Tissue material (attenuation data).
    EnergyScore m_energyScored;  ///< Bulk tissue energy accumulator.
    DoseScore m_doseScored;      ///< Bulk tissue dose accumulator.
    EnergyScore m_skin_energy;   ///< Skin shell energy accumulator.
    DoseScore m_skin_dose;       ///< Skin shell dose accumulator.
    std::array<ScoreBoxChild, 7> m_dose_boxes; ///< Seven TG-195 dose-scoring voxels.
};
}

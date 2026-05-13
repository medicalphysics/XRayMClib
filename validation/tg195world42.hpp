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

/**
 * @brief TG195 Case 4.2 phantom: a homogeneous cylinder with two embedded dose-scoring cylinders.
 *
 * Models the geometry from AAPM TG-195 Case 4.2. The phantom is a single large cylinder
 * (default radius 16 cm, height 10 cm) filled with a uniform material. Two small cylindrical
 * scoring volumes (radius 0.5 cm, half-height 5 cm, axis along Z) are embedded inside:
 *  - **Centre child**: coaxial with the phantom cylinder.
 *  - **Periphery child**: offset 1 cm inside the outer wall along -X.
 *
 * Dose scoring index convention (used by energyScored() and doseScored()):
 *   - 0 : centre child cylinder
 *   - 1 : periphery child cylinder
 *   - anything else : remaining bulk volume
 *
 * @tparam NMaterialShells      Number of electron shells in the material model (default 5).
 * @tparam LOWENERGYCORRECTION  Low-energy correction mode for the interaction kernel (default 2).
 */
template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TG195World42 {
protected:
    /**
     * @brief Tests whether a particle is inside a child cylinder using a fast AABB pre-check.
     *
     * The AABB rejection avoids the more expensive cylinder test for most particles.
     *
     * @param p     Particle whose position is tested.
     * @param child Cylinder geometry of the child scoring volume.
     * @param aabb  Pre-computed AABB of @p child.
     * @return      True if the particle position is inside @p child.
     */
    static bool insideChild(const Particle& p, const basicshape::cylinder::Cylinder& child, const std::array<double, 6>& aabb)
    {
        if (basicshape::AABB::pointInside(p.pos, aabb)) {
            return basicshape::cylinder::pointInside(p.pos, child);
        }
        return false;
    }

public:
    /**
     * @brief Constructs the TG-195 Case 4.2 phantom and positions the two scoring cylinders.
     *
     * The phantom axis is aligned with Z. The periphery scoring cylinder is placed at
     * x = center[0] - (radius - 1), i.e. 1 cm inside the outer wall.
     * Both child cylinders are given tight AABBs for fast ray-box pre-rejection.
     * The bulk material defaults to dry air; call setMaterial() before simulating.
     *
     * @param radius  Outer radius of the phantom cylinder in cm (default 16).
     * @param height  Full height of the phantom cylinder in cm (default 10).
     * @param pos     Centre position of the phantom in world coordinates (default origin).
     */
    TG195World42(double radius = 16, double height = 10, const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_cylinder.center = pos;
        m_cylinder.direction = { 0, 0, 1 };
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");

        // Centre scoring cylinder: coaxial with the phantom, radius 0.5 cm, half-height 5 cm.
        m_centerChild.center = m_cylinder.center;
        m_centerChild.direction = m_cylinder.direction;
        m_centerChild.radius = 0.5;
        m_centerChild.half_height = 5;

        // Periphery scoring cylinder: offset to 1 cm inside the outer wall along -X.
        m_periferyChild.center = { m_cylinder.center[0] - (m_cylinder.radius - 1), m_cylinder.center[1], m_cylinder.center[2] };
        m_periferyChild.direction = m_cylinder.direction;
        m_periferyChild.radius = 0.5;
        m_periferyChild.half_height = 5;

        // Build tight AABBs for the two child cylinders (XY extents from radius, Z from half_height).
        for (std::size_t i = 0; i < 3; ++i) {
            m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.radius;
            m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.radius;
            m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.radius;
            m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.radius;
            if (i == 2) {
                // Override Z extents with the (larger) half_height instead of radius.
                m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.half_height;
                m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.half_height;
                m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.half_height;
                m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.half_height;
            }
        }
    }

    /**
     * @brief Sets the bulk material for the entire phantom volume (including child scoring cylinders).
     * @param material  Material definition to use for attenuation lookups.
     */
    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    /**
     * @brief Sets the mass density of the bulk material.
     * @param density  Mass density in g/cm³.
     */
    void setMaterialDensity(double density) { m_materialDensity = density; }

    /**
     * @brief Translates the entire phantom (outer cylinder and both child cylinders) by a vector.
     * @param dist  Translation vector {dx, dy, dz} in cm.
     */
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

    /**
     * @brief Returns the centre of the outer phantom cylinder.
     * @return Centre coordinates {x, y, z} in cm.
     */
    std::array<double, 3> center() const
    {
        return m_cylinder.center;
    }

    /**
     * @brief Returns the axis-aligned bounding box of the outer phantom cylinder.
     * @return AABB as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     */
    std::array<double, 6> AABB() const
    {
        return xraymc::basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    /**
     * @brief Computes the intersection of a particle ray with the outer cylinder boundary.
     * @param p  Particle whose position and direction define the ray.
     * @return   WorldIntersectionResult for the next boundary crossing.
     */
    WorldIntersectionResult intersect(const Particle& p) const noexcept
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    /**
     * @brief Computes the intersection for visualization, returning surface normal information.
     * @tparam U  Scalar type for intersection distances.
     * @param p   Particle defining the ray.
     * @return    VisualizationIntersectionResult with distance, normal, and inside-flag.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const Particle& p) const noexcept
    {
        return basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
    }

    /**
     * @brief Transports a particle through the phantom, sampling interactions until it exits or is absorbed.
     *
     * Attenuation values are cached between steps and only recomputed when the particle energy changes
     * (e.g. after a Compton scatter), avoiding redundant material lookups on photo-electric absorption
     * or coherent scatter steps where energy is unchanged.
     *
     * After each interaction the particle position is tested against the two child cylinders:
     * energy is routed to the matching child scorer, or discarded (bulk volume) if neither matches.
     *
     * @param p      Particle to transport; modified in place.
     * @param state  Random number generator state.
     */
    void transport(Particle& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            // Recompute attenuation only when particle energy changed.
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1.0 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            const auto stepLen = -std::log(state.randomUniform()) * attSumInv;

            const auto intLen = intersect(p);
            if (stepLen < intLen.intersection) {
                // Interaction occurs before the boundary — advance and interact.
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
                    // Interaction in bulk volume — energy not scored separately.
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                }
            } else {
                // Step reaches the boundary — transport to the surface and exit.
                p.border_translate(intLen.intersection);
                cont = false;
            }
        }
    }

    /**
     * @brief Resets all energy-scored accumulators (bulk, centre child, periphery child).
     */
    void clearEnergyScored()
    {
        m_energyScored.clear();
        m_centerChild_energyScored.clear();
        m_periferyChild_energyScored.clear();
    }

    /**
     * @brief Returns the energy scorer for the centre child cylinder.
     * @return Const reference to the centre child EnergyScore.
     */
    const EnergyScore& energyScoredCenterCylinder() const
    {
        return m_centerChild_energyScored;
    }

    /**
     * @brief Returns the energy scorer for the periphery child cylinder.
     * @return Const reference to the periphery child EnergyScore.
     */
    const EnergyScore& energyScoredPeriferyCylinder() const
    {
        return m_periferyChild_energyScored;
    }

    /**
     * @brief Returns the energy scorer for a given scoring region by index.
     * @param index  0: centre child; 1: periphery child; anything else: bulk volume.
     * @return Const reference to the corresponding EnergyScore.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_centerChild_energyScored;
        else if (index == 1)
            return m_periferyChild_energyScored;
        return m_energyScored;
    }

    /**
     * @brief Converts accumulated energy scores to absorbed dose for all scoring regions.
     *
     * The bulk volume is the total phantom volume minus the two child cylinder volumes.
     * All three regions share the same material density.
     *
     * @param calibration_factor  Optional scale factor applied to all dose values (default 1).
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto c_vol = m_centerChild.volume();
        const auto p_vol = m_periferyChild.volume();
        // Subtract child volumes so energy is not double-counted in the bulk dose.
        const auto vol = m_cylinder.volume() - c_vol - p_vol;
        m_doseScored.addScoredEnergy(m_energyScored, vol, m_materialDensity, calibration_factor);
        m_centerChild_doseScored.addScoredEnergy(m_centerChild_energyScored, c_vol, m_materialDensity, calibration_factor);
        m_periferyChild_doseScored.addScoredEnergy(m_periferyChild_energyScored, p_vol, m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose scorer for a given scoring region by index.
     * @param index  0: centre child; 1: periphery child; anything else: bulk volume.
     * @return Const reference to the corresponding DoseScore.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_centerChild_doseScored;
        else if (index == 1)
            return m_periferyChild_doseScored;
        return m_doseScored;
    }

    /**
     * @brief Resets all dose accumulators (bulk, centre child, periphery child).
     */
    void clearDoseScored()
    {
        m_doseScored.clear();
        m_centerChild_doseScored.clear();
        m_periferyChild_doseScored.clear();
    }

private:
    basicshape::cylinder::Cylinder m_cylinder;          ///< Outer phantom cylinder geometry.
    double m_materialDensity = 1;                       ///< Bulk material mass density in g/cm³.
    Material<NMaterialShells> m_material;               ///< Bulk material (attenuation data).
    EnergyScore m_energyScored;                         ///< Bulk volume energy accumulator.
    EnergyScore m_periferyChild_energyScored;           ///< Periphery child energy accumulator.
    EnergyScore m_centerChild_energyScored;             ///< Centre child energy accumulator.
    DoseScore m_doseScored;                             ///< Bulk volume dose accumulator.
    DoseScore m_periferyChild_doseScored;               ///< Periphery child dose accumulator.
    DoseScore m_centerChild_doseScored;                 ///< Centre child dose accumulator.
    basicshape::cylinder::Cylinder m_centerChild;       ///< Centre scoring cylinder geometry.
    basicshape::cylinder::Cylinder m_periferyChild;     ///< Periphery scoring cylinder geometry.
    std::array<double, 6> m_centerChild_aabb = { 0, 0, 0, 0, 0, 0 };   ///< Tight AABB for centre child (fast rejection).
    std::array<double, 6> m_periferyChild_aabb = { 0, 0, 0, 0, 0, 0 }; ///< Tight AABB for periphery child (fast rejection).
};
}

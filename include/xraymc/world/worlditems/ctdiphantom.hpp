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
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/statickdtree.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief A CTDI (Computed Tomography Dose Index) phantom for CT dosimetry simulations.
 *
 * The phantom models the standard PMMA cylindrical body used in CT dose measurements:
 * a solid PMMA cylinder with five cylindrical air holes — one at the center and four
 * at the periphery (0°, 90°, 180°, 270°). Each air hole acts as an independent dose
 * scorer. Particles entering the PMMA body undergo forced photoelectric interactions
 * by default (FORCEDINTERACTIONS = true) for variance reduction; particles entering
 * an air hole are transported through it without interaction. Dose to each hole
 * can be retrieved individually via doseScored().
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 * @tparam FORCEDINTERACTIONS  When true (default), forces photoelectric interactions in
 *                             the PMMA body for variance reduction.
 */
template <int NMaterialShells = 16, int LOWENERGYCORRECTION = 2, bool FORCEDINTERACTIONS = true>
class CTDIPhantom {
public:
    /**
     * @brief Constructs a CTDI phantom: a PMMA cylinder with five air holes (one
     *        central, four peripheral at 0°/90°/180°/270°) for dose measurement.
     * @param radius    Outer radius of the PMMA cylinder in cm; clamped to at least
     *                  six times the hole radius.
     * @param height    Full height of the cylinder in cm; clamped to at least holeHeight().
     * @param pos       World-space center of the cylinder.
     * @param direction Axis direction of the cylinder (unit vector).
     */
    CTDIPhantom(double radius = 16, double height = 15, const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& direction = { 0, 0, 1 })
        : m_pmma(Material<NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
        , m_air(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        radius = std::max(std::abs(radius), holeRadii() * 6);
        height = std::max(std::abs(height), holeHeight());
        m_cylinder.center = pos;
        m_cylinder.direction = direction;
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_pmma_density = NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");
        m_air_density = NISTMaterials::density("Air, Dry (near sea level)");

        std::vector<CTDIAirHole> holes;
        holes.reserve(5);
        constexpr std::array<double, 10> positions = { 0, 0, 0, 1, 1, 0, 0, -1, -1, 0 };
        std::uint8_t index = 0;
        for (std::size_t i = 0; i < 10; i = i + 2) {
            const auto x = positions[i];
            const auto y = positions[i + 1];
            CTDIAirHole hole;
            hole.cylinder.center = m_cylinder.center;
            hole.cylinder.half_height = holeHeight() / 2;
            hole.cylinder.direction = m_cylinder.direction;
            hole.cylinder.radius = holeRadii();
            hole.index = index++;
            hole.cylinder.center[0] += x * (radius - holeEdgeDistance());
            hole.cylinder.center[1] += y * (radius - holeEdgeDistance());
            holes.push_back(hole);
        }
        m_kdtree.setData(holes);
        return;
    }

    /**
     * @brief Translates the phantom (cylinder and all air holes) by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        m_kdtree.translate(dist);
    }

    /// @brief Returns the world-space center of the phantom cylinder.
    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    /**
     * @brief Overrides the material filling the measurement holes (default: dry air).
     * @param nistName NIST material name; silently ignored if not found in the database.
     * @param density  Mass density to use for the hole material in g/cm³.
     */
    void setHoleMaterial(const std::string& nistName, double density)
    {
        auto m = Material<NMaterialShells>::byNistName(nistName);
        if (m) {
            m_air = m.value();
            m_air_density = density;
        }
    }

    /**
     * @brief Returns the energy-score accumulator for one measurement hole.
     * @param index Hole index: 0 = center, 1–4 = peripheral holes; defaults to center.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScore.at(index);
    }

    /// @brief Resets all five hole energy-score accumulators to zero.
    void clearEnergyScored()
    {
        for (auto& d : m_energyScore) {
            d.clear();
        }
    }

    /**
     * @brief Accumulates the current per-hole energy scores into the dose scores,
     *        converting energy to dose using the hole volume, hole material density,
     *        and a calibration factor.
     * @param calibration_factor Optional scaling factor applied to each scored energy value.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto holeVolume = holeHeight() * std::numbers::pi_v<double> * holeRadii() * holeRadii();

        for (std::size_t i = 0; i < m_energyScore.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScore[i], holeVolume, m_air_density, calibration_factor);
        }
    }

    /**
     * @brief Returns the dose-score accumulator for one measurement hole.
     * @param index Hole index: 0 = center, 1–4 = peripheral holes; defaults to center.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose.at(index);
    }

    /// @brief Returns the absorbed dose (Gy) scored in the central hole.
    const double centerDoseScored() const
    {
        return m_dose[0].dose();
    }

    /// @brief Returns the mean absorbed dose (Gy) averaged over the four peripheral holes.
    const double pheriferyDoseScored() const
    {
        return (m_dose[1].dose() + m_dose[2].dose() + m_dose[3].dose() + m_dose[4].dose()) / 4;
    }

    /// @brief Resets all five hole dose-score accumulators to zero.
    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    /// @brief Returns the axis-aligned bounding box of the phantom as {xmin, ymin, zmin, xmax, ymax, zmax}.
    std::array<double, 6> AABB() const
    {
        return basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    /**
     * @brief Tests a particle ray against the outer PMMA cylinder.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result containing the distance to the cylinder surface.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    /**
     * @brief Like intersect(), but returns a visualization result that includes the
     *        dose at the nearest hole (peripheral or central) and the surface normal.
     * @tparam U Scalar type used for the dose value in the visualization result.
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto res = basicshape::cylinder::intersectVisualization<U>(p, m_cylinder);
        if (res.valid()) {
            const std::array tbox = { res.intersection, res.intersection + m_cylinder.radius };
            auto holes = m_kdtree.intersect(p, tbox);
            if (holes.valid()) {
                res.value = m_dose[holes.item->index].dose();
            } else {
                res.value = m_dose[0].dose();
            }
        }
        return res;
    }

    /**
     * @brief Transports a particle through the phantom until it exits or is absorbed.
     *        Uses forced interactions inside air holes when FORCEDINTERACTIONS is true,
     *        otherwise samples analog interactions.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transport(ParticleType auto& p, RandomState& state)
    {
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        while (cont) {
            if (updateAtt) {
                att = m_pmma.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_pmma_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intCTDI = intersect(p); // this can not be nullopt
            const auto intLen = intCTDI.intersection;
            const std::array<double, 2> tbox { 0.0, intLen };
            const auto intHoles = m_kdtree.intersect(p, tbox);

            if (intHoles.valid()) {
                // We intersect a hole
                if (intHoles.rayOriginIsInsideItem) {
                    // We are inside a hole
                    if constexpr (FORCEDINTERACTIONS) {
                        const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intHoles.intersection, m_air_density, p, m_air, state);
                        m_energyScore[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        cont = intRes.particleAlive;
                    } else {
                        const auto airAtt = m_air.attenuationValues(p.energy);
                        const auto airStepLen = -std::log(state.randomUniform()) / (airAtt.sum() * m_air_density);
                        if (airStepLen < intHoles.intersection) {
                            const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(airAtt, p, m_air, state);
                            m_energyScore[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                            updateAtt = intRes.particleEnergyChanged;
                            cont = intRes.particleAlive;
                        } else {
                            p.border_translate(intHoles.intersection);
                        }
                    }
                } else {
                    if (stepLen < intHoles.intersection) {
                        // interaction happends before particle hit hole
                        p.translate(stepLen);
                        const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_pmma, state);
                        updateAtt = intRes.particleEnergyChanged;
                        cont = intRes.particleAlive;
                    } else {
                        // transport particle to hole
                        p.border_translate(intHoles.intersection);
                        if constexpr (FORCEDINTERACTIONS) {
                            // find distance of hole crossing
                            const auto dist = intHoles.item->intersect(p).intersection;
                            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(dist, m_air_density, p, m_air, state);
                            m_energyScore[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                            updateAtt = intRes.particleEnergyChanged;
                            cont = intRes.particleAlive;
                        }
                    }
                }
            } else {
                if (stepLen < intLen) {
                    // interaction happends
                    p.translate(stepLen);
                    const auto intRes = interactions::interact(att, p, m_pmma, state);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else {
                    // transport to border
                    p.border_translate(intLen);
                    cont = false;
                }
            }
        }
    }

    /// @brief Returns the radius of each measurement hole in cm (0.5 cm).
    static constexpr double holeRadii() noexcept
    {
        return 0.5;
    }
    /// @brief Returns the full height of each measurement hole in cm (10 cm).
    static constexpr double holeHeight() noexcept
    {
        return 10.0;
    }

    /// @brief Returns the full height of the phantom cylinder in cm.
    double height() const
    {
        return m_cylinder.half_height * 2;
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "CTDIPhantom1" + std::to_string(NMaterialShells) + std::to_string(LOWENERGYCORRECTION) + std::to_string(FORCEDINTERACTIONS);
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether a raw data buffer begins with the expected magic identifier.
     * @param data Buffer to inspect; must be at least 32 bytes for a positive result.
     * @return true if the first 32 bytes match magicID().
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the phantom (geometry, hole material, and dose scores) to a
     *        byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_cylinder.center, buffer);
        Serializer::serialize(m_cylinder.radius, buffer);
        Serializer::serialize(m_cylinder.direction, buffer);
        Serializer::serialize(m_cylinder.half_height, buffer);

        Serializer::serialize(m_air_density, buffer);

        Serializer::serializeMaterialWeights(m_air.composition(), buffer);

        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a phantom from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed phantom on success, or std::nullopt if the hole material
     *         weight set cannot be parsed.
     */
    static std::optional<CTDIPhantom<NMaterialShells, LOWENERGYCORRECTION, FORCEDINTERACTIONS>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> center;
        buffer = Serializer::deserialize(center, buffer);
        double radius;
        buffer = Serializer::deserialize(radius, buffer);
        std::array<double, 3> direction;
        buffer = Serializer::deserialize(direction, buffer);
        double half_height;
        buffer = Serializer::deserialize(half_height, buffer);

        CTDIPhantom<NMaterialShells, LOWENERGYCORRECTION, FORCEDINTERACTIONS> item(radius, half_height * 2, center, direction);
        buffer = Serializer::deserialize(item.m_air_density, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        if (auto mat_opt = Material<NMaterialShells>::byWeight(mat_weights); mat_opt) {
            item.m_air = mat_opt.value();
        } else {
            return std::nullopt;
        }
        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        return item;
    }

protected:
    struct CTDIAirHole {
        basicshape::cylinder::Cylinder cylinder;
        std::uint8_t index = 0;

        /// @brief Translates the hole cylinder center by @p d.
        void translate(const std::array<double, 3>& d) noexcept
        {
            for (std::size_t i = 0; i < 3; ++i)
                cylinder.center[i] += d[i];
        }

        /// @brief Returns the world-space center of the hole cylinder.
        const std::array<double, 3>& center() const noexcept { return cylinder.center; }

        /// @brief Returns the axis-aligned bounding box of the hole cylinder.
        std::array<double, 6> AABB() const noexcept
        {
            return basicshape::cylinder::cylinderAABB(cylinder);
        }

        /// @brief Tests a particle ray against the hole cylinder and returns the intersection result.
        auto intersect(const ParticleType auto& p) const noexcept
        {
            return basicshape::cylinder::intersect(p, cylinder);
        }
    };

    /// @brief Returns the distance from a peripheral hole center to the outer cylinder edge in cm (1 cm).
    static constexpr double holeEdgeDistance() noexcept
    {
        return 1.0;
    }

private:
    basicshape::cylinder::Cylinder m_cylinder;
    double m_pmma_density = 0;
    double m_air_density = 0;
    std::array<EnergyScore, 5> m_energyScore;
    std::array<DoseScore, 5> m_dose;
    StaticKDTree<CTDIAirHole> m_kdtree;
    Material<NMaterialShells> m_pmma;
    Material<NMaterialShells> m_air;
};
}

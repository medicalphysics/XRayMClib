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
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

/**
 * @brief A hollow box geometry representing a shielded room: a uniform-material wall
 *        of configurable thickness surrounds an air-filled inner cavity.
 *        Energy is scored in the wall volume only; the inner cavity is transparent.
 *
 * @tparam NMaterialShells    Number of electron shells used in material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 */
template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2>
class EnclosedRoom {
public:
    /**
     * @brief Constructs a room from an explicit inner AABB and wall thickness.
     * @param wallthickness Wall thickness in cm; clamped to at least 0.001 cm.
     * @param inner_aabb    Inner cavity as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     */
    EnclosedRoom(double wallthickness, const std::array<double, 6>& inner_aabb)
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_wallThickness = std::max(std::abs(wallthickness), 0.001);
        setInnerRoomAABB(inner_aabb);
        m_density = NISTMaterials::density("Air, Dry (near sea level)");
    }

    /**
     * @brief Constructs a room centered at the origin from per-axis inner dimensions.
     * @param wallthickness   Wall thickness in cm; clamped to at least 0.001 cm.
     * @param inner_room_size Inner cavity extents {dx, dy, dz} in cm; centered at origin.
     */
    EnclosedRoom(double wallthickness, const std::array<double, 3>& inner_room_size)
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        std::array<double, 6> inner_aabb;
        for (std::size_t i = 0; i < 3; ++i) {
            inner_aabb[i] = -inner_room_size[i] * 0.5;
            inner_aabb[i + 3] = inner_room_size[i] * 0.5;
        }
        m_wallThickness = std::max(std::abs(wallthickness), 0.001);
        setInnerRoomAABB(inner_aabb);
        m_density = NISTMaterials::density("Air, Dry (near sea level)");
    }

    /**
     * @brief Constructs a cubic room centered at the origin with equal inner side lengths.
     * @param wallthickness  Wall thickness in cm; clamped to at least 0.001 cm.
     * @param inner_room_size Inner cavity side length in cm; same for all three axes.
     */
    EnclosedRoom(double wallthickness = 10, double inner_room_size = 2)
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        std::array<double, 6> inner_aabb;
        for (std::size_t i = 0; i < 3; ++i) {
            inner_aabb[i] = -inner_room_size * 0.5;
            inner_aabb[i + 3] = inner_room_size * 0.5;
        }
        m_wallThickness = std::max(std::abs(wallthickness), 0.001);
        setInnerRoomAABB(inner_aabb);
        m_density = NISTMaterials::density("Air, Dry (near sea level)");
    }

    /// @brief Returns the current wall thickness in cm.
    const double wallThickness() const
    {
        return m_wallThickness;
    }

    /**
     * @brief Sets the wall thickness and recomputes the outer AABB.
     * @param cm New thickness in cm; clamped to at least 0.001 cm.
     */
    void setWallThickness(double cm)
    {
        m_wallThickness = std::max(std::abs(cm), 0.001);
        m_outerAABB = outerAABfromInner(m_innerAABB, m_wallThickness);
    }

    /// @brief Sets the wall material without changing the density.
    void setMaterial(const Material<NMaterialShells>& material) { m_material = material; }

    /**
     * @brief Sets the wall material and its mass density together.
     * @param material Wall material definition.
     * @param density  Mass density in g/cm³.
     */
    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        m_material = material;
        setDensity(density);
    }

    /**
     * @brief Sets the wall material mass density in g/cm³.
     * @param dens Absolute value is used.
     */
    void setDensity(double dens)
    {
        m_density = std::abs(dens);
    }

    /**
     * @brief Sets the inner cavity AABB and recomputes the outer AABB.
     * @param aabb Inner cavity as {xmin, ymin, zmin, xmax, ymax, zmax}.
     */
    void setInnerRoomAABB(const std::array<double, 6>& aabb)
    {
        m_innerAABB = aabb;
        m_outerAABB = outerAABfromInner(m_innerAABB, m_wallThickness);
    }

    /**
     * @brief Translates both inner and outer AABBs by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_innerAABB[i] += dist[i];
            m_innerAABB[i + 3] += dist[i];
        }
    }

    /// @brief Returns the geometric center of the inner cavity.
    std::array<double, 3> center() const
    {
        std::array center = {
            (m_innerAABB[0] + m_innerAABB[3]) / 2,
            (m_innerAABB[1] + m_innerAABB[4]) / 2,
            (m_innerAABB[2] + m_innerAABB[5]) / 2
        };
        return center;
    }

    /// @brief Returns the outer AABB (wall exterior) as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_outerAABB;
    }

    /// @brief Returns the inner cavity AABB as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& innerRoomAABB() const
    {
        return m_innerAABB;
    }

    /**
     * @brief Returns the distance to the nearest wall surface along the particle ray.
     *        Rays originating inside the cavity are directed toward the inner wall;
     *        rays inside the wall material are directed toward whichever surface
     *        (inner or outer) is closer.
     * @param p Particle whose position and direction define the ray.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const bool is_inside_inner = basicshape::AABB::pointInside(p.pos, m_innerAABB);
        WorldIntersectionResult intersect;
        if (is_inside_inner) {
            intersect = basicshape::AABB::intersect(p, m_innerAABB);
            intersect.rayOriginIsInsideItem = false;
        } else {
            intersect = basicshape::AABB::intersect(p, m_outerAABB);
            if (intersect.rayOriginIsInsideItem) {
                // test if we intersect inner aabb
                const auto intersect_inner = basicshape::AABB::intersect(p, m_innerAABB);
                if (intersect_inner.valid() && intersect_inner.intersection < intersect.intersection) {
                    intersect = intersect_inner;
                }
            }
        }
        return intersect;
    }

    /**
     * @brief Like intersect(), but skips the near wall faces so that the interior is
     *        visible during rendering, and returns the surface normal and dose value.
     * @tparam U Scalar type used for the dose value in the visualization result.
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        VisualizationIntersectionResult<U> intersection = basicshape::AABB::template intersectVisualization<U>(p, m_innerAABB);
        if (intersection.valid()) {
            if (intersection.rayOriginIsInsideItem) {
                intersection.rayOriginIsInsideItem = false;
            } else {
                // we intersect inner box from outside and want to render closest walls invisible
                auto p_copy = p;
                p_copy.border_translate(intersection.intersection);
                const auto past_wall_intersection = basicshape::AABB::template intersectVisualization<U>(p_copy, m_innerAABB);
                intersection.intersection += past_wall_intersection.intersection;
                intersection.normal = past_wall_intersection.normal;
            }
        } else {
            intersection = basicshape::AABB::template intersectVisualization<U>(p, m_outerAABB);
        }
        return intersection;
    }

    /**
     * @brief Returns the energy-score accumulator for the wall volume.
     * @param index Unused; present for interface consistency.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScore;
    }

    /// @brief Resets the wall energy-score accumulator to zero.
    void clearEnergyScored()
    {
        m_energyScore.clear();
    }

    /**
     * @brief Accumulates the wall energy score into the dose score, using the wall
     *        volume (outer minus inner), wall density, and a calibration factor.
     * @param calibration_factor Optional scaling factor applied to the scored energy.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        std::array<double, 3> inner_sides;
        std::array<double, 3> outer_sides;
        for (std::size_t i = 0; i < 3; ++i) {
            inner_sides[i] = m_innerAABB[i + 3] - m_innerAABB[i];
            outer_sides[i] = m_outerAABB[i + 3] - m_outerAABB[i];
        }
        const auto inner_volume = std::reduce(inner_sides.cbegin(), inner_sides.cend(), 1.0, std::multiplies<>());
        const auto outer_volume = std::reduce(outer_sides.cbegin(), outer_sides.cend(), 1.0, std::multiplies<>());

        m_dose.addScoredEnergy(m_energyScore, outer_volume - inner_volume, m_density, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the wall volume.
     * @param index Unused; present for interface consistency.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    /// @brief Resets the wall dose-score accumulator to zero.
    void clearDoseScored()
    {
        m_dose.clear();
    }

    /**
     * @brief Transports a particle through the wall material until it exits or is
     *        absorbed, scoring deposited energy in the wall accumulator.
     *        Particles inside the inner cavity are not transported in this method.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = pointInside(p.pos);
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                m_energyScore.scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "EnclosedRoom1" + std::to_string(NMaterialShells) + std::to_string(LOWENERGYCORRECTION);
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
     * @brief Serializes the room (geometry, wall material, density, and dose score) to
     *        a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_wallThickness, buffer);
        Serializer::serialize(m_innerAABB, buffer);
        Serializer::serialize(m_density, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a room from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed room on success, or std::nullopt if the material
     *         weight set cannot be parsed.
     */
    static std::optional<EnclosedRoom<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
    {
        EnclosedRoom<NMaterialShells, LOWENERGYCORRECTION> item;

        buffer = Serializer::deserialize(item.m_wallThickness, buffer);
        buffer = Serializer::deserialize(item.m_innerAABB, buffer);
        item.m_outerAABB = item.outerAABfromInner(item.m_innerAABB, item.m_wallThickness);
        buffer = Serializer::deserialize(item.m_density, buffer);
        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        if (auto mat_opt = Material<NMaterialShells>::byWeight(mat_weights); mat_opt) {
            item.m_material = mat_opt.value();
        } else {
            return std::nullopt;
        }
        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        return item;
    }

protected:
    /**
     * @brief Computes the outer AABB by expanding @p inner by @p wall_thickness on all sides.
     * @param inner          Inner cavity AABB as {xmin, ymin, zmin, xmax, ymax, zmax}.
     * @param wall_thickness Thickness to add to each face in cm.
     * @return Outer AABB.
     */
    static std::array<double, 6> outerAABfromInner(const std::array<double, 6>& inner, double wall_thickness)
    {
        std::array<double, 6> aabb = {
            inner[0] - wall_thickness,
            inner[1] - wall_thickness,
            inner[2] - wall_thickness,
            inner[3] + wall_thickness,
            inner[4] + wall_thickness,
            inner[5] + wall_thickness
        };
        return aabb;
    }

    /**
     * @brief Returns true when @p p is inside the wall material (inside the outer AABB
     *        but outside the inner AABB).
     */
    bool pointInside(const std::array<double, 3>& p) const noexcept
    {
        return basicshape::AABB::pointInside(p, m_outerAABB) && !basicshape::AABB::pointInside(p, m_innerAABB);
    }

private:
    double m_wallThickness = 10;
    std::array<double, 6> m_innerAABB = { -1, -1, -1, 1, 1, 1 };
    std::array<double, 6> m_outerAABB = { -1, -1, -1, 1, 1, 1 };
    double m_density = 1;
    Material<NMaterialShells> m_material;
    EnergyScore m_energyScore;
    DoseScore m_dose;
};
}
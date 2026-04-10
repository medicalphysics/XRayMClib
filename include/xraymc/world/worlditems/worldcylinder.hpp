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
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

/**
 * @brief An arbitrarily oriented solid cylinder geometry for Monte Carlo particle transport.
 *
 * The cylinder is parameterised by a center, an axis direction, a radius, and a full height.
 * It is filled with a single homogeneous material, defaulting to PMMA (Lucite/Perspex).
 * Particles are tracked inside the volume using analog Monte Carlo free-path sampling until
 * they exit through the curved side or one of the two end caps. Absorbed energy is
 * accumulated in an EnergyScore and can be converted to dose via addEnergyScoredToDoseScore().
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 */
template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2>
class WorldCylinder {
public:
    /**
     * @brief Constructs a cylinder, defaulting to PMMA at its standard NIST density.
     * @param radius Cylinder radius in cm; absolute value is used.
     * @param height Full height of the cylinder in cm; absolute value is used.
     * @param center World-space center of the cylinder in cm.
     * @param dir    Axis direction of the cylinder; normalized internally.
     */
    WorldCylinder(double radius = 16, double height = 10, const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 })
        : m_material(Material<NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
    {
        m_cylinder.radius = std::abs(radius);
        m_cylinder.half_height = std::abs(height) / 2;
        m_cylinder.center = center;
        m_cylinder.direction = vectormath::normalized(dir);
        m_materialDensity = NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");
    }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const WorldCylinder<NMaterialShells, LOWENERGYCORRECTION>&) const = default;

    /**
     * @brief Sets the material used for attenuation and interaction sampling.
     * @param material Material cross-section data.
     */
    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    /**
     * @brief Sets both the material and its mass density in one call.
     * @param material Material cross-section data.
     * @param density  Density in g/cm³; absolute value is used.
     */
    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        m_material = material;
        setMaterialDensity(density);
    }

    /**
     * @brief Sets the material mass density.
     * @param density Density in g/cm³; absolute value is used.
     */
    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
    }

    /**
     * @brief Sets the material and density from the NIST material database by name.
     * @param nist_name Name of the NIST material (e.g., "Water, Liquid").
     * @return true if the name was found and the material was set; false otherwise.
     */
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

    /**
     * @brief Translates the cylinder center by @p dist and updates the AABB.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        updateAABB();
    }

    /// @brief Returns the world-space center of the cylinder in cm.
    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    /**
     * @brief Sets the world-space center of the cylinder and updates the AABB.
     * @param c New center position in cm.
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_cylinder.center = c;
        updateAABB();
    }

    /// @brief Returns the axis-aligned bounding box of the cylinder as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Sets the cylinder radius and updates the AABB.
     * @param r Radius in cm; absolute value is used.
     */
    void setRadius(double r)
    {
        m_cylinder.radius = std::abs(r);
        updateAABB();
    }

    /// @brief Returns the cylinder radius in cm.
    double radius() const
    {
        return m_cylinder.radius;
    }

    /**
     * @brief Sets the full cylinder height and updates the AABB.
     * @param h Full height in cm; absolute value is used.
     */
    void setHeight(double h)
    {
        m_cylinder.half_height = std::abs(h / 2);
        updateAABB();
    }

    /// @brief Returns the full height of the cylinder in cm.
    double height() const
    {
        return m_cylinder.half_height * 2;
    }

    /**
     * @brief Sets the cylinder axis direction and updates the AABB.
     * @param dir Desired axis direction; normalized internally.
     */
    void setDirection(const std::array<double, 3>& dir)
    {
        m_cylinder.direction = vectormath::normalized(dir);
        updateAABB();
    }

    /// @brief Returns the normalized cylinder axis direction.
    const std::array<double, 3>& direction() const
    {
        return m_cylinder.direction;
    }

    /**
     * @brief Tests a particle ray against the cylinder (curved side and end caps).
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the nearest surface, including
     *         whether the ray origin is inside the cylinder.
     */
    WorldIntersectionResult
    intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    /**
     * @brief Like intersect(), but also returns the surface normal and the accumulated
     *        dose value for visualization.
     * @tparam U Scalar type used for the value field (set to the current mean dose).
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    /**
     * @brief Transports a particle through the cylinder using analog Monte Carlo sampling.
     *
     * Free-path lengths are sampled exponentially; if the sampled path is shorter than
     * the distance to the exit surface an interaction occurs, otherwise the particle is
     * advanced to the border. Continues until the particle exits or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Returns the energy-score accumulator for the cylinder.
     * @param index Unused; present for interface consistency.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    /// @brief Resets the energy-score accumulator to zero.
    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    /**
     * @brief Converts the accumulated energy score to a dose score using the cylinder volume.
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        m_dose.addScoredEnergy(m_energyScored, m_cylinder.volume(), m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the cylinder.
     * @param index Unused; present for interface consistency.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    /// @brief Resets the dose-score accumulator to zero.
    void clearDoseScored()
    {
        m_dose.clear();
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "Cylinder1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells);
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
     * @brief Serializes the cylinder (geometry, material density, material composition,
     *        and dose score) to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        // Serialize cylinder
        Serializer::serialize(m_cylinder.center, buffer);
        Serializer::serialize(m_cylinder.direction, buffer);
        Serializer::serialize(m_cylinder.radius, buffer);
        Serializer::serialize(m_cylinder.half_height, buffer);
        // materials
        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a cylinder from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed cylinder (always valid for well-formed input).
     */
    static std::optional<WorldCylinder<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
    {
        std::array<double, 3> center, direction;
        buffer = Serializer::deserialize(center, buffer);
        buffer = Serializer::deserialize(direction, buffer);

        double radius, half_heigh;
        buffer = Serializer::deserialize(radius, buffer);
        buffer = Serializer::deserialize(half_heigh, buffer);

        WorldCylinder<NMaterialShells, LOWENERGYCORRECTION> item(radius, half_heigh * 2, center, direction);

        double density;
        buffer = Serializer::deserialize(density, buffer);
        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);

        auto material_opt = Material<NMaterialShells>::byWeight(mat_weights);
        if (material_opt) {
            item.setMaterial(material_opt.value(), density);
        }

        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        return std::make_optional(item);
    }

protected:
    /// @brief Recomputes the AABB from the current cylinder geometry.
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

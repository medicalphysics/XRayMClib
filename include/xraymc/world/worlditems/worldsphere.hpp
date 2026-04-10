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
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/sphere.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

/**
 * @brief A solid sphere geometry for Monte Carlo particle transport.
 *
 * The sphere is parameterised by a center and a radius and filled with a single
 * homogeneous material. Particles are tracked inside the volume using analog or forced
 * Monte Carlo sampling until they exit through the spherical surface. Absorbed energy
 * is accumulated in an EnergyScore and can be converted to dose via
 * addEnergyScoredToDoseScore(). When FORCEINTERACTIONS is true, deterministic forced
 * photoelectric interactions are applied at each step for variance reduction.
 * An optional ParticleTracker records particle histories when ParticleTrack particles
 * are transported.
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 * @tparam FORCEINTERACTIONS   When true, forces photoelectric interactions at each step
 *                             for variance reduction; when false, uses analog sampling.
 */
template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class WorldSphere {
public:
    /**
     * @brief Constructs a sphere, initialized with dry air at standard sea-level density.
     * @param radius Sphere radius in cm; absolute value is used.
     * @param pos    World-space center of the sphere in cm.
     */
    WorldSphere(double radius = 16, const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_radius(std::abs(radius))
        , m_center(pos)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
    }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const WorldSphere<NMaterialShells, LOWENERGYCORRECTION, FORCEINTERACTIONS>& other) const = default;

    /**
     * @brief Sets the sphere radius.
     * @param r Radius in cm; absolute value is used.
     */
    void setRadius(double r)
    {
        m_radius = std::abs(r);
    }

    /// @brief Returns the sphere radius in cm.
    auto radius() const { return m_radius; }

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
    void setMaterialDensity(double density) { m_materialDensity = std::abs(density); }

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
     * @brief Translates the sphere center by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
    }

    /**
     * @brief Sets the world-space center of the sphere.
     * @param c New center position in cm.
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
    }

    /// @brief Returns the world-space center of the sphere in cm.
    std::array<double, 3> center() const
    {
        return m_center;
    }

    /// @brief Returns the axis-aligned bounding box of the sphere as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
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

    /**
     * @brief Tests a particle ray against the sphere.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the sphere surface, including
     *         whether the ray origin is inside the sphere.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::sphere::intersect(p, m_center, m_radius);
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
        auto inter = basicshape::sphere::template intersectVisualization<U>(p, m_center, m_radius);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    /**
     * @brief Transports a particle through the sphere, dispatching to forced or analog
     *        sampling based on the FORCEINTERACTIONS template parameter.
     *        If @p p is a ParticleTrack, it is registered with the particle tracker before transport.
     * @tparam P Particle type; if ParticleTrack, tracking is enabled.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Returns the energy-score accumulator for the sphere.
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
     * @brief Converts the accumulated energy score to a dose score using the sphere volume (4πr³/3).
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = (4 * std::numbers::pi_v<double> * m_radius * m_radius * m_radius) / 3;
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the sphere.
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

    /// @brief Returns the particle tracker (const overload).
    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    /// @brief Returns the particle tracker (mutable overload).
    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "Sphere1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(FORCEINTERACTIONS) + std::to_string(NMaterialShells);
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
     * @brief Serializes the sphere (radius, center, material density, material composition,
     *        and dose score) to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_radius, buffer);
        Serializer::serialize(m_center, buffer);
        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);

        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a sphere from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed sphere, or std::nullopt if the material composition cannot be resolved.
     */
    static std::optional<WorldSphere<NMaterialShells, LOWENERGYCORRECTION, FORCEINTERACTIONS>> deserialize(std::span<const char> buffer)
    {
        WorldSphere<NMaterialShells, LOWENERGYCORRECTION, FORCEINTERACTIONS> item;

        buffer = Serializer::deserialize(item.m_radius, buffer);
        buffer = Serializer::deserialize(item.m_center, buffer);
        buffer = Serializer::deserialize(item.m_materialDensity, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        if (auto matCand = Material<NMaterialShells>::byWeight(mat_weights); matCand)
            item.m_material = matCand.value();
        else
            return std::nullopt;

        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        return std::make_optional(item);
    }

protected:
    /**
     * @brief Transports a particle through the sphere using analog Monte Carlo sampling.
     *
     * Free-path lengths are sampled exponentially; if the sampled path is shorter than
     * the distance to the exit surface an interaction occurs, otherwise the particle is
     * advanced to the border. Continues until the particle exits or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Transports a particle through the sphere using forced (deterministic) interactions.
     *
     * At each step the distance to the exit surface is computed and a forced photoelectric
     * interaction is applied over that path length with the appropriate survival weight.
     * Continues until the particle exits or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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
    double m_radius = 0;
    std::array<double, 3> m_center;
    Material<NMaterialShells> m_material;
    double m_materialDensity = 1;
    EnergyScore m_energyScored;
    DoseScore m_dose;
    ParticleTracker m_tracker;
};
}

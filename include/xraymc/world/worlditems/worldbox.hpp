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
 * @brief An axis-aligned box (AABB) geometry for Monte Carlo particle transport.
 *
 * The box is defined by its axis-aligned bounding box and filled with a single
 * homogeneous material. Particles are tracked inside the volume until they exit
 * through one of the six faces. Absorbed energy is accumulated in an EnergyScore
 * and can be converted to dose via addEnergyScoredToDoseScore(). When
 * FORCE_INTERACTIONS is true, forced (deterministic) photoelectric interactions are
 * used at every step instead of analog random sampling, providing variance reduction.
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 * @tparam FORCE_INTERACTIONS  When true, forces photoelectric interactions at each step
 *                             for variance reduction; when false, uses analog sampling.
 */
template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2, bool FORCE_INTERACTIONS = false>
class WorldBox {
public:
    /**
     * @brief Constructs a box from an explicit AABB, initialized with dry air at
     *        standard sea-level density.
     * @param aabb Box extents as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     *             Min/max pairs are swapped automatically if out of order.
     */
    WorldBox(const std::array<double, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : m_aabb(aabb)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
    }

    /**
     * @brief Constructs a symmetric box of half-size @p aabb_size centered at @p pos,
     *        initialized with dry air at standard sea-level density.
     * @param aabb_size Half-extent of the box in cm; absolute value is used.
     * @param pos       World-space center of the box in cm.
     */
    WorldBox(double aabb_size, std::array<double, 3> pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        correctAABB();
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
    }

    /**
     * @brief Constructs a box from an explicit AABB with a given material and density.
     * @param aabb     Box extents as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     * @param material Material cross-section data.
     * @param density  Material mass density in g/cm³; absolute value is used.
     */
    WorldBox(const std::array<double, 6>& aabb, const Material<NMaterialShells>& material, double density)
        : m_aabb(aabb)
        , m_material(material)
    {
        m_materialDensity = std::abs(density);
        correctAABB();
    }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const WorldBox<NMaterialShells, LOWENERGYCORRECTION, FORCE_INTERACTIONS>&) const = default;

    /**
     * @brief Sets the box extents and corrects the AABB if min/max are out of order.
     * @param aabb New extents as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     */
    void setAABB(const std::array<double, 6>& aabb)
    {
        m_aabb = aabb;
        correctAABB();
    }

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
     * @brief Translates the box by @p dist by shifting all six AABB planes.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist) noexcept
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    /// @brief Returns the world-space center of the box in cm.
    std::array<double, 3> center() const noexcept
    {
        std::array<double, 3> c {
            (m_aabb[0] + m_aabb[3]) * 0.5,
            (m_aabb[1] + m_aabb[4]) * 0.5,
            (m_aabb[2] + m_aabb[5]) * 0.5,
        };
        return c;
    }

    /// @brief Returns the axis-aligned bounding box as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
    const std::array<double, 6>& AABB() const noexcept
    {
        return m_aabb;
    }

    /**
     * @brief Tests a particle ray against the box faces.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the entry or exit face, including
     *         whether the ray origin is inside the box.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const noexcept
    {
        return basicshape::AABB::intersect(p, m_aabb);
    }

    /**
     * @brief Like intersect(), but also returns the face normal and the accumulated
     *        dose value for visualization.
     * @tparam U Scalar type used for the value field (set to the current mean dose).
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::AABB::template intersectVisualization<U>(p, m_aabb);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    /**
     * @brief Dispatches to transportForced() or transportRandom() depending on
     *        the FORCE_INTERACTIONS template parameter.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        if constexpr (FORCE_INTERACTIONS) {
            transportForced(p, state);
        } else {
            transportRandom(p, state);
        }
    }

    /**
     * @brief Returns the energy-score accumulator for the box.
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
     * @brief Converts the accumulated energy score to a dose score using the box volume.
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto [l, h] = vectormath::splice(m_aabb);
        const auto sides = vectormath::subtract(h, l);
        const auto volume = std::reduce(sides.cbegin(), sides.cend(), double { 1 }, std::multiplies<>());

        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    /**
     * @brief Returns the dose-score accumulator for the box.
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
        std::string name = "WorldBox1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells) + std::to_string(FORCE_INTERACTIONS);
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
     * @brief Serializes the box (AABB, material density, material composition, and dose score)
     *        to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_aabb, buffer);
        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a box from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed box, or std::nullopt if the material composition cannot be resolved.
     */
    static std::optional<WorldBox<NMaterialShells, LOWENERGYCORRECTION, FORCE_INTERACTIONS>> deserialize(std::span<const char> buffer)
    {

        std::array<double, 6> aabb;
        buffer = Serializer::deserialize(aabb, buffer);

        double dens;
        buffer = Serializer::deserialize(dens, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);

        auto mat_opt = Material<NMaterialShells>::byWeight(mat_weights);

        WorldBox<NMaterialShells, LOWENERGYCORRECTION, FORCE_INTERACTIONS> item(aabb, mat_opt.value(), dens);

        if (mat_opt) {
            item.setMaterial(mat_opt.value(), dens);
        } else {
            return std::nullopt;
        }
        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);
        return item;
    }

protected:
    /**
     * @brief Transports a particle through the box using forced (deterministic) interactions.
     *
     * At each step the distance to the exit face is computed and a forced photoelectric
     * interaction is applied over that path length, with the appropriate survival weight.
     * Continues until the particle exits or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transportForced(ParticleType auto& p, RandomState& state) noexcept
    {

        bool cont = basicshape::AABB::pointInside(p.pos, m_aabb);
        while (cont) {
            const auto intLen = intersect(p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, m_materialDensity, p, m_material, state);
            m_energyScored.scoreEnergy(intRes.energyImparted);
            cont = intRes.particleAlive && basicshape::AABB::pointInside(p.pos, m_aabb);
        }
    }

    /**
     * @brief Transports a particle through the box using analog Monte Carlo sampling.
     *
     * Free-path lengths are sampled exponentially; if the sampled path is shorter than
     * the distance to the exit face an interaction occurs, otherwise the particle is
     * advanced to the border. Continues until the particle exits or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transportRandom(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::AABB::pointInside(p.pos, m_aabb);
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
            const auto intLen = intersect(p).intersection; // this can not be nullopt

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
     * @brief Ensures the AABB is valid by swapping any inverted min/max pairs and
     *        expanding degenerate (zero-extent) axes by GEOMETRIC_ERROR().
     */
    void correctAABB()
    {
        auto test = [](const auto& aabb) -> bool {
            bool ok = true;
            for (std::size_t i = 0; i < 3; ++i)
                ok = ok && aabb[i] < aabb[i + 3];
            return ok;
        };
        if (!test(m_aabb)) {
            for (std::size_t i = 0; i < 3; ++i) {
                if (m_aabb[i] == m_aabb[i + 3]) {
                    m_aabb[i] -= GEOMETRIC_ERROR();
                    m_aabb[i + 3] += GEOMETRIC_ERROR();
                } else if (m_aabb[i] > m_aabb[i + 3]) {
                    std::swap(m_aabb[i], m_aabb[i + 3]);
                }
            }
        }
    }

private:
    std::array<double, 6> m_aabb;
    double m_materialDensity = 1;
    Material<NMaterialShells> m_material;
    EnergyScore m_energyScored;
    DoseScore m_dose;
};

}

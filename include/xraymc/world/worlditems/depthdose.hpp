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
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>

namespace xraymc {

template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class DepthDose {
public:
    /**
     * @brief Constructs a depth-dose scorer: a uniformly filled cylinder divided
     *        into @p resolution axial slabs, defaulting to dry air at sea level.
     * @param radius     Cylinder radius in cm.
     * @param height     Full cylinder height in cm.
     * @param resolution Number of axial scoring bins along the cylinder axis.
     * @param pos        World-space center of the cylinder.
     * @param dir        Cylinder axis direction (unit vector).
     */
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

    /**
     * @brief Sets the cylinder radius in cm and updates the bounding box.
     * @param radius Absolute value is used; negative values are treated as positive.
     */
    void setRadius(double radius)
    {
        m_cylinder.radius = std::abs(radius);
        updateAABB();
    }

    /**
     * @brief Sets the cylinder axis direction and updates the bounding box.
     * @param direction Desired axis vector; normalized internally.
     */
    void setDirection(const std::array<double, 3>& direction)
    {
        m_cylinder.direction = vectormath::normalized(direction);
        updateAABB();
    }

    /**
     * @brief Sets the full cylinder height in cm and updates the bounding box.
     * @param lenght Absolute value is used; negative values are treated as positive.
     */
    void setLenght(double lenght)
    {
        m_cylinder.half_height = std::abs(lenght * 0.5);
        updateAABB();
    }

    /**
     * @brief Sets the world-space center of the cylinder and updates the bounding box.
     * @param position New center position in cm.
     */
    void setCenter(const std::array<double, 3>& position)
    {
        m_cylinder.center = position;
        updateAABB();
    }

    /**
     * @brief Changes the number of axial scoring bins, clearing all accumulated scores.
     * @param resolution New number of depth bins.
     */
    void setResolution(std::size_t resolution)
    {
        clearEnergyScored();
        clearDoseScored();
        m_energyScored.resize(resolution);
        m_dose.resize(resolution);
    }

    /// @brief Returns the full cylinder height in cm.
    double length() const { return m_cylinder.half_height * 2; }

    /// @brief Returns the cylinder radius in cm.
    double radius() const { return m_cylinder.radius; }

    /// @brief Returns the cylinder axis direction (unit vector).
    const std::array<double, 3>& direction() const
    {
        return m_cylinder.direction;
    }
    /// @brief Returns the world-space center of the cylinder.
    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    /// @brief Returns the number of axial scoring bins.
    std::size_t resolution() const
    {
        return m_energyScored.size();
    }

    /**
     * @brief Sets the fill material without changing the density.
     * @param material Material definition to use for attenuation and interaction.
     */
    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    /**
     * @brief Sets the fill material mass density in g/cm³.
     * @param density Absolute value is used.
     */
    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
    }

    /**
     * @brief Sets the fill material and its mass density together.
     * @param material Material definition.
     * @param density  Mass density in g/cm³.
     */
    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        setMaterial(material);
        setMaterialDensity(density);
    }

    /**
     * @brief Sets the fill material and its NIST-tabulated density by name.
     * @param nist_name NIST material name (e.g. "Water, Liquid").
     * @return true on success; false if the name is not found in the database.
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

    /// @brief Returns the current fill material.
    const Material<NMaterialShells>& material() const
    {
        return m_material;
    }

    /// @brief Returns the fill material mass density in g/cm³.
    double density() const { return m_materialDensity; }

    /**
     * @brief Translates the cylinder center by @p dist and updates the bounding box.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        updateAABB();
    }

    /// @brief Returns the axis-aligned bounding box of the cylinder as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Tests a particle ray against the cylinder.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result containing the distance to the cylinder surface.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    /**
     * @brief Like intersect(), but returns a visualization result that includes the
     *        dose of the depth bin at the hit point and the cylinder surface normal.
     * @tparam U Scalar type used for the dose value in the visualization result.
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        auto inter = basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
        if (inter.valid()) {
            auto p_int = p;
            p_int.translate(inter.intersection);
            const auto ind = cylinderIndex(p_int.pos);
            inter.value = m_dose[ind].dose();
        }
        return inter;
    }

    /**
     * @brief Transports a particle through the cylinder until it exits or is absorbed,
     *        scoring deposited energy in the corresponding depth bin.
     *        Registers the particle track when @p P is ParticleTrack.
     *        Uses forced interactions when FORCEINTERACTIONS is true.
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
     * @brief Returns per-bin energy scores paired with the axial position of each bin center.
     * @return Vector of (depth in cm, EnergyScore) for each scoring bin, ordered from
     *         the bottom to the top of the cylinder.
     */
    const std::vector<std::pair<double, EnergyScore>> depthEnergyScored() const
    {
        std::vector<std::pair<double, EnergyScore>> depth;
        depth.reserve(m_dose.size());

        const auto step = (2 * m_cylinder.half_height) / m_dose.size();
        const auto start = m_cylinder.center[2] - m_cylinder.half_height + step / 2;
        for (std::size_t i = 0; i < m_dose.size(); ++i) {
            depth.push_back(std::make_pair(start + step * i, m_energyScored[i]));
        }
        return depth;
    }

    /**
     * @brief Returns per-bin dose scores paired with the axial position of each bin center.
     * @return Vector of (depth in cm, DoseScore) for each scoring bin, ordered from
     *         the bottom to the top of the cylinder.
     */
    const std::vector<std::pair<double, DoseScore>> depthDoseScored() const
    {
        std::vector<std::pair<double, DoseScore>> depth;
        depth.reserve(m_dose.size());

        const auto step = (2 * m_cylinder.half_height) / m_dose.size();
        const auto start = m_cylinder.center[2] - m_cylinder.half_height + step / 2;
        for (std::size_t i = 0; i < m_dose.size(); ++i) {
            depth.push_back(std::make_pair(start + step * i, m_dose[i]));
        }
        return depth;
    }

    /**
     * @brief Returns the energy-score accumulator for a single depth bin.
     * @param index Bin index (0 = bottom bin); defaults to the first bin.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored[index];
    }

    /// @brief Resets all depth-bin energy-score accumulators to zero.
    void clearEnergyScored()
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    /**
     * @brief Accumulates the current per-bin energy scores into the dose scores,
     *        converting energy to dose using the bin volume, material density, and a
     *        calibration factor.
     * @param calibration_factor Optional scaling factor applied to each scored energy value.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto totalVolume = m_cylinder.volume();
        const auto partVolume = totalVolume / m_dose.size();
        for (std::size_t i = 0; i < m_dose.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScored[i], partVolume, m_materialDensity, calibration_factor);
        }
    }

    /**
     * @brief Returns the dose-score accumulator for a single depth bin.
     * @param index Bin index (0 = bottom bin); defaults to the first bin.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose[index];
    }

    /// @brief Resets all depth-bin dose-score accumulators to zero.
    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    /// @brief Returns a read-only reference to the particle tracker (populated when
    ///        transporting ParticleTrack particles).
    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    /// @brief Returns a mutable reference to the particle tracker.
    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "DepthDose1" + std::to_string(NMaterialShells) + std::to_string(LOWENERGYCORRECTION) + std::to_string(FORCEINTERACTIONS);
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
     * @brief Serializes the scorer (geometry, material, density, and dose scores) to
     *        a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_cylinder.center, buffer);
        Serializer::serialize(m_cylinder.radius, buffer);
        Serializer::serialize(m_cylinder.direction, buffer);
        Serializer::serialize(m_cylinder.half_height, buffer);
        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a scorer from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed scorer on success, or std::nullopt if the material
     *         weight set cannot be parsed.
     */
    static std::optional<DepthDose<NMaterialShells, LOWENERGYCORRECTION, FORCEINTERACTIONS>> deserialize(std::span<const char> buffer)
    {
        DepthDose<NMaterialShells, LOWENERGYCORRECTION, FORCEINTERACTIONS> item;

        buffer = Serializer::deserialize(item.m_cylinder.center, buffer);
        buffer = Serializer::deserialize(item.m_cylinder.radius, buffer);
        buffer = Serializer::deserialize(item.m_cylinder.direction, buffer);
        buffer = Serializer::deserialize(item.m_cylinder.half_height, buffer);
        buffer = Serializer::deserialize(item.m_materialDensity, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        if (auto mat_opt = Material<NMaterialShells>::byWeight(mat_weights); mat_opt) {
            item.m_material = mat_opt.value();
        } else {
            return std::nullopt;
        }
        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

        item.updateAABB();

        return item;
    }

protected:
    /// @brief Recomputes the AABB from the current cylinder geometry.
    void updateAABB()
    {
        m_aabb = basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    /**
     * @brief Maps a world-space position to the corresponding depth-bin index.
     * @tparam BOUNDS_CHECK When true, clamps the result to [0, resolution-1];
     *                      when false, no clamping is performed.
     * @param pos World-space position; typically a point inside the cylinder.
     * @return Zero-based bin index along the cylinder axis.
     */
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

    /**
     * @brief Analog (random) transport: samples a free path and either interacts
     *        inside the current bin or exits the cylinder.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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

    /**
     * @brief Forced-interaction transport: guarantees an interaction within each depth
     *        bin (or at the cylinder boundary), improving scoring efficiency for thin
     *        or low-density geometries.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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
    std::array<double, 6> m_aabb;                  ///< Axis-aligned bounding box {xmin, ymin, zmin, xmax, ymax, zmax} [cm]; recomputed whenever geometry changes.
    xraymc::basicshape::cylinder::Cylinder m_cylinder; ///< Cylinder geometry (center, direction, radius, half-height) defining the scorer volume.
    double m_materialDensity = 1;                  ///< Mass density of the fill material [g/cm³].
    Material<NMaterialShells> m_material;          ///< Cross-section data for the fill material; defaults to dry air at sea level.
    ParticleTracker m_tracker;                     ///< Accumulates particle track segments when transporting ParticleTrack particles.
    std::vector<EnergyScore> m_energyScored;       ///< Per-bin energy imparted accumulators; size = resolution (one entry per axial slab).
    std::vector<DoseScore> m_dose;                 ///< Per-bin absorbed dose accumulators; parallel to m_energyScored.
};

}

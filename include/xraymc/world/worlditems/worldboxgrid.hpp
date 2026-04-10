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

#include <cmath>
#include <limits>
#include <optional>

namespace xraymc {

/**
 * @brief An axis-aligned box divided into a regular voxel grid for spatially resolved
 *        dose and energy scoring.
 *
 * The box is filled with a single homogeneous material, but its interior is subdivided
 * into a configurable NxMxK voxel grid. Each voxel has its own EnergyScore and
 * DoseScore accumulator, enabling 3-D dose maps. Particles are transported using
 * analog Monte Carlo free-path sampling, and energy deposited by each interaction is
 * credited to the voxel containing the interaction point.
 *
 * @tparam NMaterialShells     Number of electron shells for material cross-sections.
 * @tparam LOWENERGYCORRECTION Low-energy correction mode passed to interaction sampling.
 */
template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2>
class WorldBoxGrid {
public:
    /**
     * @brief Constructs a grid box from an explicit AABB with a single voxel,
     *        initialized with dry air at standard sea-level density.
     * @param aabb Box extents as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     *             Min/max pairs are corrected automatically if out of order.
     */
    WorldBoxGrid(const std::array<double, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : m_aabb(aabb)
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
        setVoxelDimensions({ 1, 1, 1 });
    }

    /**
     * @brief Constructs a grid box from an explicit AABB with a given material and density.
     * @param aabb     Box extents as {xmin, ymin, zmin, xmax, ymax, zmax} in cm.
     * @param material Material cross-section data.
     * @param density  Material mass density in g/cm³.
     */
    WorldBoxGrid(const std::array<double, 6>& aabb, const Material<NMaterialShells>& material, double density = 1)
        : m_aabb(aabb)
        , m_materialDensity(density)
        , m_material(material)
    {
        correctAABB();
        setVoxelDimensions({ 1, 1, 1 });
    }

    /**
     * @brief Constructs a symmetric grid box of half-size @p aabb_size centered at @p pos,
     *        initialized with dry air at standard sea-level density.
     * @param aabb_size Half-extent of the box in cm; absolute value is used.
     * @param pos       World-space center of the box in cm.
     */
    WorldBoxGrid(double aabb_size, std::array<double, 3> pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        correctAABB();
        setVoxelDimensions({ 1, 1, 1 });
    }

    /// @brief Defaulted equality comparison (compares all data members).
    bool operator==(const WorldBoxGrid<NMaterialShells, LOWENERGYCORRECTION>&) const = default;

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
     * @brief Sets the voxel grid dimensions and resizes the per-voxel score arrays.
     *
     * Each dimension is clamped to [1, 1000]. Voxel sizes and their reciprocals are
     * recomputed from the current AABB. All existing score data is discarded.
     * @param dim Number of voxels along {x, y, z}.
     */
    void setVoxelDimensions(const std::array<std::uint64_t, 3>& dim)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_voxelDim[i] = std::min(std::max(std::uint64_t { 1 }, dim[i]), std::uint64_t { 1000 });

        for (std::size_t i = 0; i < 3; ++i) {
            m_voxelSize[i] = (m_aabb[i + 3] - m_aabb[i]) / m_voxelDim[i];
            m_voxelSizeInv[i] = 1.0 / m_voxelSize[i];
        }
        const auto ndim = m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
        m_energyScored.resize(ndim);
        m_dose.resize(ndim);
    }

    /// @brief Returns the total number of voxels (product of all three dimensions).
    std::uint64_t totalNumberOfVoxels() const
    {
        return m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
    }

    /// @brief Returns the voxel grid dimensions as {nx, ny, nz}.
    const std::array<std::uint64_t, 3>& voxelDimensions() const
    {
        return m_voxelDim;
    }

    /// @brief Returns the voxel spacing in cm as {dx, dy, dz}.
    const std::array<double, 3>& voxelSpacing() const
    {
        return m_voxelSize;
    }

    /**
     * @brief Converts a world-space position to a flat voxel index using row-major (x-fastest) order.
     * @tparam BOUNDSCHECK When true (default), clamps each coordinate to the valid voxel range
     *                     before computing the index. When false, no bounds check is performed
     *                     (faster, but undefined behavior if @p pos is outside the box).
     * @param pos World-space position in cm.
     * @return Flat voxel index in [0, totalNumberOfVoxels()).
     */
    template <bool BOUNDSCHECK = true>
    std::uint64_t gridIndex(const std::array<double, 3>& pos) const noexcept
    {
        if constexpr (BOUNDSCHECK) {
            const auto x = static_cast<std::uint64_t>(std::clamp((pos[0] - m_aabb[0]) * m_voxelSizeInv[0], double { 0 }, static_cast<double>(m_voxelDim[0] - 1)));
            const auto y = static_cast<std::uint64_t>(std::clamp((pos[1] - m_aabb[1]) * m_voxelSizeInv[1], double { 0 }, static_cast<double>(m_voxelDim[1] - 1)));
            const auto z = static_cast<std::uint64_t>(std::clamp((pos[2] - m_aabb[2]) * m_voxelSizeInv[2], double { 0 }, static_cast<double>(m_voxelDim[2] - 1)));
            return x + (y + z * m_voxelDim[1]) * m_voxelDim[0];
        } else {
            const auto x = static_cast<std::uint64_t>((pos[0] - m_aabb[0]) * m_voxelSizeInv[0]);
            const auto y = static_cast<std::uint64_t>((pos[1] - m_aabb[1]) * m_voxelSizeInv[1]);
            const auto z = static_cast<std::uint64_t>((pos[2] - m_aabb[2]) * m_voxelSizeInv[2]);
            return x + (y + z * m_voxelDim[1]) * m_voxelDim[0];
        }
    }

    /**
     * @brief Translates the box by @p dist by shifting all six AABB planes.
     *        Note: voxel sizes are not recomputed; call setVoxelDimensions() afterwards
     *        if the AABB dimensions change.
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
     * @brief Like intersect(), but also returns the face normal and the dose value of the
     *        voxel at the intersection point for visualization.
     * @tparam U Scalar type used for the value field (set to the dose of the hit voxel).
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        auto inter = basicshape::AABB::template intersectVisualization<U>(p, m_aabb);
        if (inter.valid()) {
            auto p_int = p;
            p_int.translate(inter.intersection);
            const auto ind = gridIndex(p_int.pos);
            inter.value = m_dose[ind].dose();
        }
        return inter;
    }

    /**
     * @brief Transports a particle through the grid using analog Monte Carlo sampling.
     *
     * Free-path lengths are sampled exponentially; if the sampled path is shorter than
     * the distance to the exit face an interaction occurs and the deposited energy is
     * credited to the voxel containing the interaction point. Continues until the particle
     * exits the box or is absorbed.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
    void transport(ParticleType auto& p, RandomState& state) noexcept
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
                const auto scoreIdx = gridIndex<false>(p.pos);
                m_energyScored[scoreIdx].scoreEnergy(intRes.energyImparted);
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
     * @brief Returns the energy-score accumulator for the voxel at @p index.
     * @param index Flat voxel index; throws std::out_of_range if out of bounds.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored.at(index);
    }

    /// @brief Resets all per-voxel energy-score accumulators to zero.
    void clearEnergyScored()
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    /**
     * @brief Converts each per-voxel energy score to a dose score using the voxel volume.
     * @param calibration_factor Optional scaling factor applied during the conversion; defaults to 1.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = m_voxelSize[0] * m_voxelSize[1] * m_voxelSize[2];
        for (std::size_t i = 0; i < m_energyScored.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScored[i], volume, m_materialDensity, calibration_factor);
        }
    }

    /**
     * @brief Returns the dose-score accumulator for the voxel at @p index.
     * @param index Flat voxel index; throws std::out_of_range if out of bounds.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose.at(index);
    }

    /// @brief Resets all per-voxel dose-score accumulators to zero.
    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "WorldBoxGrid1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells);
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
     * @brief Serializes the grid (AABB, material density, voxel dimensions, material
     *        composition, and per-voxel dose scores) to a byte vector that can be
     *        restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_aabb, buffer);
        Serializer::serialize(m_materialDensity, buffer);
        Serializer::serialize(m_voxelDim, buffer);
        Serializer::serializeMaterialWeights(m_material.composition(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a grid box from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed grid, or std::nullopt if the material composition cannot be resolved.
     */
    static std::optional<WorldBoxGrid<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
    {

        std::array<double, 6> aabb;
        buffer = Serializer::deserialize(aabb, buffer);

        double dens;
        buffer = Serializer::deserialize(dens, buffer);
        std::array<std::uint64_t, 3> dim;
        buffer = Serializer::deserialize(dim, buffer);

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);

        auto mat_opt = Material<NMaterialShells>::byWeight(mat_weights);
        if (mat_opt) {
            WorldBoxGrid<NMaterialShells, LOWENERGYCORRECTION> item(aabb, mat_opt.value(), dens);
            item.setMaterialDensity(dens);
            item.setVoxelDimensions(dim);
            buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);

            return std::make_optional(item);
        }
        return std::nullopt;
    }

protected:
    /**
     * @brief Ensures the AABB is valid by swapping inverted min/max pairs and expanding
     *        near-degenerate axes (within GEOMETRIC_ERROR()) by GEOMETRIC_ERROR().
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
                if (std::abs(m_aabb[i] - m_aabb[i + 3]) < GEOMETRIC_ERROR()) {
                    if (m_aabb[i] > m_aabb[i + 3]) {
                        std::swap(m_aabb[i], m_aabb[i + 3]);
                    }
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
    std::array<double, 3> m_voxelSize = { 1, 1, 1 };
    std::array<double, 3> m_voxelSizeInv = { 1, 1, 1 };
    std::array<std::uint64_t, 3> m_voxelDim = { 1, 1, 1 };
    Material<NMaterialShells> m_material;
    std::vector<EnergyScore> m_energyScored;
    std::vector<DoseScore> m_dose;
};
}

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
#include <vector>

namespace xraymc {

/**
 * @brief A zero-thickness, infinitely thin ideal flat-panel detector.
 *
 * The detector is an oriented rectangle in 3-D space, parameterised by a center
 * position, two direction-cosine vectors (row and column), pixel spacing, and
 * pixel grid dimensions.  Arriving particles are absorbed immediately and their
 * energy (weighted by particle weight) is scored in the corresponding pixel.
 * Per-pixel energy and dose accumulators are provided, along with optional
 * particle-track recording for ParticleTrack particles.
 */
class FlatDetector {
public:
    /**
     * @brief Default constructor. Creates a 128×128 detector centered at the origin
     *        with 0.1 cm pixel spacing, lying in the XY plane.
     */
    FlatDetector()
    {
        const auto nPix = m_detector_dimensions[0] * m_detector_dimensions[1];
        m_energyScore.resize(nPix);
        m_doseScore.resize(nPix);
        m_aabb = calculateAABB();
    }

    /// @brief Returns the axis-aligned bounding box of the detector panel as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Translates the detector center by @p dist and updates the AABB.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
        m_aabb = calculateAABB();
    }

    /// @brief Returns the world-space center of the detector panel.
    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    /**
     * @brief Sets the world-space center of the detector panel and updates the AABB.
     * @param c New center position in cm.
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        m_aabb = calculateAABB();
    }

    /**
     * @brief Sets the row and column direction cosines that orient the detector panel.
     *        The vectors must be orthogonal; if not, the call is silently ignored.
     *        Both vectors are normalized internally.
     * @param dir_cosX Row direction (column-index increases along this vector).
     * @param dir_cosY Column direction (row-index increases along this vector).
     */
    void setDirectionCosines(const std::array<double, 3>& dir_cosX, const std::array<double, 3>& dir_cosY)
    {
        if (vectormath::dot(dir_cosX, dir_cosY) < 2 * GEOMETRIC_ERROR<>()) {
            m_direction_cosines[0] = vectormath::normalized(dir_cosX);
            m_direction_cosines[1] = vectormath::normalized(dir_cosY);
        }
        m_aabb = calculateAABB();
    }

    /// @brief Returns the two direction-cosine vectors {row direction, column direction}.
    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_direction_cosines;
    }

    /**
     * @brief Sets pixel spacing from an array and updates the AABB.
     * @param spacing {spacingX, spacingY} in cm; values below 0.00001 are clamped.
     */
    void setPixelSpacing(const std::array<double, 2>& spacing)
    {
        setPixelSpacing(spacing[0], spacing[1]);
    }

    /**
     * @brief Sets per-axis pixel spacing and updates the AABB.
     * @param spacingX Column pixel size in cm; clamped to at least 0.00001 cm.
     * @param spacingY Row pixel size in cm; clamped to at least 0.00001 cm.
     */
    void setPixelSpacing(double spacingX, double spacingY)
    {
        m_pixel_spacing[0] = std::max(spacingX, 0.00001);
        m_pixel_spacing[1] = std::max(spacingY, 0.00001);
        m_aabb = calculateAABB();
    }

    /// @brief Returns the pixel spacing as {spacingX, spacingY} in cm.
    std::array<double, 2> pixelSpacing() const
    {
        return m_pixel_spacing;
    }

    /**
     * @brief Sets the pixel grid dimensions from an array and updates the AABB.
     * @param dimensions {columns, rows}; each clamped to at least 1.
     */
    void setDetectorDimensions(const std::array<std::size_t, 2>& dimensions)
    {
        setDetectorDimensions(dimensions[0], dimensions[1]);
    }

    /**
     * @brief Sets the pixel grid dimensions and updates the AABB.
     *        Resizes the energy and dose score vectors accordingly.
     * @param dimX Number of columns; clamped to at least 1.
     * @param dimY Number of rows; clamped to at least 1.
     */
    void setDetectorDimensions(std::size_t dimX, std::size_t dimY)
    {
        constexpr std::size_t min = 1;
        m_detector_dimensions[0] = std::max(dimX, min);
        m_detector_dimensions[1] = std::max(dimY, min);
        m_energyScore.resize(dimX * dimY);
        m_doseScore.resize(dimX * dimY);
        m_aabb = calculateAABB();
    }

    /// @brief Returns the pixel grid dimensions as {columns, rows}.
    const std::array<std::size_t, 2>& detectorDimensions() const
    {
        return m_detector_dimensions;
    }

    /// @brief Returns the detector surface normal (cross product of the two direction cosines).
    std::array<double, 3> normal() const
    {
        return vectormath::cross(m_direction_cosines[0], m_direction_cosines[1]);
    }

    /**
     * @brief Tests a particle ray against the detector plane and pixel grid.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the detector plane, or invalid
     *         if the ray misses the pixel grid or travels away from it.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        const auto dir = vectormath::changeBasisInverse(m_direction_cosines[0], m_direction_cosines[1], normal(), p.dir);
        const auto pos = vectormath::changeBasisInverse(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto denom = dir[2];
        if (std::abs(denom) > GEOMETRIC_ERROR<>()) {
            const auto t = -pos[2] / denom;
            if (t > 0) {
                const auto x = pos[0] + dir[0] * t;
                const auto y = pos[1] + dir[1] * t;
                if (std::abs(x) <= m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5 && std::abs(y) <= m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) {
                    res.intersection = t;
                    res.rayOriginIsInsideItem = false;
                    res.intersectionValid = true;
                }
            }
        }
        return res;
    }

    /**
     * @brief Like intersect(), but also returns the surface normal and the dose value
     *        of the hit pixel for visualization purposes.
     * @tparam U Scalar type used for the dose value in the visualization result.
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        VisualizationIntersectionResult<U> res;
        const auto dir = vectormath::changeBasisInverse(m_direction_cosines[0], m_direction_cosines[1], normal(), p.dir);
        const auto pos = vectormath::changeBasisInverse(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto denom = dir[2];
        if (std::abs(denom) > GEOMETRIC_ERROR<>()) {
            const auto t = -pos[2] / denom;
            if (t > 0) {
                const auto x = pos[0] + dir[0] * t;
                const auto y = pos[1] + dir[1] * t;
                if (std::abs(x) <= m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5 && std::abs(y) <= m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) {
                    res.intersection = t;
                    res.rayOriginIsInsideItem = false;
                    res.intersectionValid = true;
                    res.normal = denom < 0 ? normal() : vectormath::scale(normal(), -1.0);
                    const auto indx = (x + m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5) / m_pixel_spacing[0];
                    const auto indy = (y + m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) / m_pixel_spacing[1];
                    const auto flatIdx = static_cast<std::size_t>(std::floor(indx)) * m_detector_dimensions[1] + static_cast<std::size_t>(std::floor(indy));
                    res.value = m_doseScore[flatIdx].dose();
                }
            }
        }
        return res;
    }

    /**
     * @brief Absorbs the particle into the detector pixel at its current position,
     *        scoring the weighted energy (energy × weight) and setting particle energy
     *        to zero.  Registers the track when @p P is ParticleTrack.
     * @param p     Particle to absorb; energy is zeroed after scoring.
     * @param state Random number generator state (unused, present for interface consistency).
     */
    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }

        const auto pos = vectormath::changeBasisInverse(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto x = pos[0];
        const auto y = pos[1];
        const auto indx = std::clamp((x + m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5) / m_pixel_spacing[0], 0.0, static_cast<double>(m_detector_dimensions[0] - 1));
        const auto indy = std::clamp((y + m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) / m_pixel_spacing[1], 0.0, static_cast<double>(m_detector_dimensions[1] - 1));
        const auto flatIdx = static_cast<std::size_t>(indx) * m_detector_dimensions[1] + static_cast<std::size_t>(indy);
        m_energyScore[flatIdx].scoreEnergy(p.energy * p.weight);
        p.energy = 0;
    }

    /**
     * @brief Returns the energy-score accumulator for a single pixel.
     * @param index Flat pixel index (row-major: row * columns + col); defaults to pixel 0.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScore[index];
    }

    /// @brief Resets all per-pixel energy-score accumulators to zero.
    void clearEnergyScored()
    {
        for (auto& d : m_energyScore) {
            d.clear();
        }
    }

    /**
     * @brief Accumulates per-pixel energy scores into dose scores using pixel area
     *        as the volume proxy (density = 1 g/cm²) and a calibration factor.
     * @param calibration_factor Optional scaling factor applied to each scored energy value.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = m_pixel_spacing[0] * m_pixel_spacing[1];
        constexpr double density = 1.0;
        for (std::size_t i = 0; i < m_doseScore.size(); ++i) {
            m_doseScore[i].addScoredEnergy(m_energyScore[i], volume, density, calibration_factor);
        }
    }

    /**
     * @brief Returns the dose-score accumulator for a single pixel.
     * @param index Flat pixel index (row-major); defaults to pixel 0.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_doseScore[index];
    }

    /// @brief Resets all per-pixel dose-score accumulators to zero.
    void clearDoseScored()
    {
        for (auto& d : m_doseScore) {
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
        std::string name = "FlatDetector1";
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
     * @brief Serializes the detector (geometry, orientation, and dose scores) to a byte
     *        vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_center, buffer);
        Serializer::serialize(m_pixel_spacing, buffer);
        Serializer::serialize(m_detector_dimensions, buffer);
        Serializer::serialize(m_direction_cosines[0], buffer);
        Serializer::serialize(m_direction_cosines[1], buffer);
        Serializer::serializeDoseScore(m_doseScore, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a detector from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed detector (always valid for well-formed input).
     */
    static std::optional<FlatDetector> deserialize(std::span<const char> buffer)
    {
        FlatDetector item;

        buffer = Serializer::deserialize(item.m_center, buffer);
        buffer = Serializer::deserialize(item.m_pixel_spacing, buffer);
        buffer = Serializer::deserialize(item.m_detector_dimensions, buffer);
        buffer = Serializer::deserialize(item.m_direction_cosines[0], buffer);
        buffer = Serializer::deserialize(item.m_direction_cosines[1], buffer);
        buffer = Serializer::deserializeDoseScore(item.m_doseScore, buffer);

        item.m_aabb = item.calculateAABB();

        return item;
    }

protected:
    /// @brief Returns {min, max} of two values as a pair.
    static std::pair<double, double> minmax(const double a, const double b)
    {
        return a < b ? std::pair { a, b } : std::pair { b, a };
    }

    /**
     * @brief Computes the AABB of the detector panel from the current center, direction
     *        cosines, pixel spacing, and grid dimensions. A thin normal-direction
     *        epsilon is added to give the panel non-zero volume.
     */
    std::array<double, 6> calculateAABB() const
    {
        const auto normal_vec = normal();
        const auto half_span_x = vectormath::scale(m_direction_cosines[0], m_pixel_spacing[0] * m_detector_dimensions[0] / 2);
        const auto half_span_y = vectormath::scale(m_direction_cosines[1], m_pixel_spacing[1] * m_detector_dimensions[1] / 2);
        const auto span_vecs = std::array { half_span_x, half_span_y, vectormath::scale(normal_vec, GEOMETRIC_ERROR<>()) };
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (const auto& s : span_vecs)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                aabb[i] = std::min(aabb[i], mi);
                aabb[i + 3] = std::max(aabb[i + 3], ma);
            }
        return aabb;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 };
    std::array<double, 2> m_pixel_spacing = { 0.1, 0.1 };
    std::array<std::size_t, 2> m_detector_dimensions = { 128, 128 };
    std::array<std::array<double, 3>, 2> m_direction_cosines = { 1, 0, 0, 0, 1, 0 };
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<EnergyScore> m_energyScore;
    std::vector<DoseScore> m_doseScore;
    ParticleTracker m_tracker;
};
};

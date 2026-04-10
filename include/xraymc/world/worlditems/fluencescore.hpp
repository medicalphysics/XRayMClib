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

#include "xraymc/constants.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief A circular disc scorer that tallies photon fluence and energy spectra.
 *
 * Particles crossing the disc plane within its radius are counted into an
 * energy-binned histogram (spectrum) and an EnergyScore accumulator.
 * Fluence (photons/cm²) and its variance can be derived from the raw counts
 * by normalising against the number of primary particles and the disc area.
 * Dose is not defined for this scorer; the corresponding methods are no-ops.
 */
class FluenceScore {
public:
    /**
     * @brief Constructs a fluence-scoring disc.
     * @param radius  Disc radius in cm.
     * @param center  World-space center of the disc.
     * @param normal  Surface normal of the disc plane (normalized internally).
     */
    FluenceScore(double radius = 16, const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 3>& normal = { 0, 0, 1 })
        : m_center(center)
        , m_radius(radius)
    {
        setPlaneNormal(normal);
        setEnergyStep(1);
        calculateAABB();
    }

    /**
     * @brief Sets the energy bin width for the spectrum histogram and resets all counts.
     * @param step Bin width in keV; clamped to at least 0.1 keV.
     */
    void setEnergyStep(double step)
    {
        m_energy_step = std::max(step, 0.1);
        const auto N = static_cast<std::uint64_t>(MAX_ENERGY() / m_energy_step) + 1;
        m_intensity.resize(N);
        std::fill(m_intensity.begin(), m_intensity.end(), 0);
    }

    /**
     * @brief Sets the world-space center of the disc and updates the AABB.
     * @param c New center position in cm.
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        calculateAABB();
    }

    /**
     * @brief Sets the disc radius and updates the AABB.
     * @param r Radius in cm; absolute value is used, clamped to at least 0.00001 cm.
     */
    void setRadius(double r)
    {
        m_radius = std::max(std::abs(r), 0.00001);
        calculateAABB();
    }

    /// @brief Returns the disc area in cm².
    double area() const
    {
        return std::numbers::pi_v<double> * m_radius * m_radius;
    }

    /**
     * @brief Translates the disc center and AABB by @p dist.
     * @param dist Displacement vector in cm along {x, y, z}.
     */
    void translate(const std::array<double, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_center[i] += dist[i];
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    /// @brief Returns the world-space center of the disc.
    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    /// @brief Returns the axis-aligned bounding box of the disc as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Sets the disc plane normal and updates the AABB.
     * @param normal Desired normal vector; normalized internally.
     */
    void setPlaneNormal(const std::array<double, 3>& normal)
    {
        m_normal = normal;
        vectormath::normalize(m_normal);
        calculateAABB();
    }

    /**
     * @brief Returns the raw photon count spectrum as (bin-center energy in keV, count) pairs.
     * @return Vector of (energy keV, photon count), one entry per energy bin.
     */
    std::vector<std::pair<double, std::uint64_t>> getSpecter() const
    {
        std::vector<std::pair<double, std::uint64_t>> spec(m_intensity.size());
        for (std::size_t i = 0; i < m_intensity.size(); ++i) {
            spec[i] = std::make_pair(m_energy_step * i + m_energy_step / 2, m_intensity[i]);
        }
        return spec;
    }

    /**
     * @brief Returns the differential fluence spectrum: photons per primary per cm² per bin.
     * @param N_primaries Number of primary particles used for normalization; defaults to 1.
     * @return Vector of (energy keV, fluence photons/cm²/primary), one entry per energy bin.
     */
    std::vector<std::pair<double, double>> getFluenceSpecter(std::uint64_t N_primaries = 1) const
    {
        const auto specter = getSpecter();
        std::vector<std::pair<double, double>> fluence(specter.size());
        const auto a = area();
        const auto norm = N_primaries > 0 ? static_cast<double>(N_primaries) : 1.0;
        const auto norm_inv = 1.0 / norm;
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), fluence.begin(), [a, norm_inv](const auto& s) {
            return std::make_pair(s.first, s.second * norm_inv / a);
        });
        return fluence;
    }

    /**
     * @brief Returns the per-bin variance of the fluence estimate using Bernoulli statistics.
     * @param N_primaries Number of primary particles; defaults to 1.
     * @return Vector of (energy keV, variance (photons/cm²/primary)²), one entry per energy bin.
     */
    std::vector<std::pair<double, double>> getFluenceVariance(std::uint64_t N_primaries = 1) const
    {
        const auto specter = getSpecter();
        std::vector<std::pair<double, double>> fluence(specter.size());
        const auto a = area();
        const auto norm = N_primaries > 0 ? static_cast<double>(N_primaries) : 1.0;
        const auto norm_inv = 1.0 / norm;
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), fluence.begin(), [a, norm_inv](const auto& s) {
            auto var = norm_inv * (s.second * norm_inv * (1 - s.second * norm_inv));
            return std::make_pair(s.first, var / (a * a));
        });
        return fluence;
    }

    /**
     * @brief Tests a particle ray against the disc (AABB pre-filter then exact disc test).
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the disc plane, or invalid if missed.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const auto aabb_inter = basicshape::AABB::intersectForwardInterval(p, m_aabb);
        return aabb_inter ? intersectDisc(p) : WorldIntersectionResult { };
    }

    /**
     * @brief Like intersect(), but also returns the correctly-oriented disc normal
     *        for shading during visualization.
     * @tparam U Scalar type used for the value field (unused; set to zero).
     * @param p  Particle whose position and direction define the ray.
     */
    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        const auto res = intersect(p);
        VisualizationIntersectionResult<U> w;
        if (res.valid()) {
            w.rayOriginIsInsideItem = false;
            w.intersection = res.intersection;
            w.intersectionValid = true;
            w.normal = vectormath::dot(p.dir, m_normal) <= 0 ? m_normal : vectormath::scale(m_normal, -1.0);
        }
        return w;
    }

    /**
     * @brief Returns the total energy-score accumulator for the disc.
     * @param index Unused; present for interface consistency.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    /// @brief No-op. Energy scoring is reset via clearDoseScored() for this scorer.
    void clearEnergyScored()
    {
        return;
    }

    /// @brief No-op. Dose is not defined for a fluence scorer.
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        // not defined for fluence counter
        return;
    }

    /**
     * @brief Returns an empty DoseScore. Dose is not defined for a fluence scorer.
     * @param index Unused; present for interface consistency.
     */
    const DoseScore doseScored(std::size_t index = 0) const
    {
        return DoseScore { };
    }

    /// @brief Resets the energy-score accumulator and all spectrum bin counts to zero.
    void clearDoseScored()
    {
        m_energyScored.clear();
        std::fill(m_intensity.begin(), m_intensity.end(), std::uint64_t { 0 });
        return;
    }

    /**
     * @brief Records a particle crossing: increments the energy-bin counter and the
     *        total energy-score accumulator, then advances the particle past the disc.
     *        Thread-safe for the bin counter via std::atomic_ref.
     * @param particle Particle arriving at the disc surface; border-translated in place.
     * @param state    Random number generator state (unused).
     */
    void transport(ParticleType auto& particle, RandomState& state)
    {
        // Assuming particle is on the disc
        const auto eIdx = static_cast<std::size_t>(particle.energy / m_energy_step);
        auto counter = std::atomic_ref(m_intensity[eIdx]);
        counter++;
        m_energyScored.scoreEnergy(particle.energy);
        particle.border_translate(0);
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "FluenceScore1";
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
     * @brief Serializes the scorer (geometry, spectrum histogram, and energy accumulator)
     *        to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_center, buffer);
        Serializer::serialize(m_normal, buffer);
        Serializer::serialize(m_radius, buffer);
        Serializer::serialize(m_energy_step, buffer);
        Serializer::serialize(m_intensity, buffer);
        Serializer::serialize(m_energyScored.energyImparted(), buffer);
        Serializer::serialize(m_energyScored.energyImpartedSquared(), buffer);
        Serializer::serialize(m_energyScored.numberOfEvents(), buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a scorer from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed scorer (always valid for well-formed input).
     */
    static std::optional<FluenceScore> deserialize(std::span<const char> buffer)
    {
        FluenceScore item;

        buffer = Serializer::deserialize(item.m_center, buffer);
        buffer = Serializer::deserialize(item.m_normal, buffer);
        buffer = Serializer::deserialize(item.m_radius, buffer);
        buffer = Serializer::deserialize(item.m_energy_step, buffer);
        buffer = Serializer::deserialize(item.m_intensity, buffer);

        double energy, energySq;
        buffer = Serializer::deserialize(energy, buffer);
        buffer = Serializer::deserialize(energySq, buffer);

        std::uint64_t nevents;
        buffer = Serializer::deserialize(nevents, buffer);

        item.m_energyScored.set(energy, energySq, nevents);

        item.calculateAABB();

        return std::make_optional(item);
    }

protected:
    /// @brief Returns {min, max} of two values as a pair.
    ///        Uses a hand-rolled comparison instead of std::minmax to avoid an MSVC /O2 bug or weird feature.
    static inline std::pair<double, double> minmax(double v1, double v2)
    {
        // use own minmax instead of std::minmax due to bug or weird feature of MSVC compiler with /O2
        return v1 <= v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
    }

    /// @brief Recomputes the AABB to tightly bound the disc, ensuring a minimum thickness
    ///        along any axis-aligned direction by adding a margin.
    void calculateAABB()
    {
        const std::array<std::array<double, 3>, 3> span = {
            vectormath::scale(vectormath::cross(m_normal, { 1, 0, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 1, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 0, 1 }), m_radius)
        };
        m_aabb = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };

        for (const auto& s : span)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                m_aabb[i] = std::min(m_aabb[i], mi);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], ma);
            }
        // ensure a min size;
        for (std::size_t i = 0; i < 3; ++i) {
            if (m_aabb[i + 3] - m_aabb[i] < 2 * GEOMETRIC_ERROR()) {
                m_aabb[i] -= GEOMETRIC_ERROR();
                m_aabb[i + 3] += GEOMETRIC_ERROR();
            }
        }
    }

    /**
     * @brief Exact ray–disc intersection test against the plane and radius.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the disc plane, or invalid if
     *         the ray is nearly parallel to the plane or the hit point is outside the radius.
     */
    WorldIntersectionResult intersectDisc(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        const auto D = vectormath::dot(p.dir, m_normal);
        constexpr double minOrt = 1E-6;

        if (D < -minOrt || D > minOrt) {
            res.intersection = vectormath::dot(vectormath::subtract(m_center, p.pos), m_normal) / D;
            if (res.intersection > 0) {
                // intersection point
                const auto p_int = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                // distance from center
                const auto c_dist = vectormath::subtract(m_center, p_int);
                // check if distance from center is less than radius
                res.intersectionValid = vectormath::length_sqr(c_dist) <= m_radius * m_radius;
            }
        }
        return res;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 };
    std::array<double, 3> m_normal = { 0, 0, 1 };
    double m_radius = 16;
    double m_energy_step = 1;
    std::array<double, 6> m_aabb;
    std::vector<std::uint64_t> m_intensity;
    EnergyScore m_energyScored;
};
}
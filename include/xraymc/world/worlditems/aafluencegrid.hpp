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
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <execution>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief An axis-aligned planar grid scorer that tallies a photon spectrum per pixel.
 *
 * The grid is a rectangle of pixels lying in a plane perpendicular to one of the
 * world axes. Each photon crossing the plane inside the grid is counted in the
 * energy bin of its pixel, giving a spectrum for every pixel. Only photons travelling
 * against the plane normal are scored; photons travelling along the normal pass
 * through unseen. Fluence (photons/cm²) can be derived from the raw counts by
 * normalising against the number of primary particles and pixelArea().
 *
 * After transport, addEnergyScoredToDoseScore() converts the per-pixel spectra to
 * air kerma (keV/g), which is available per pixel through doseScored() and is
 * shown as the pixel value in visualizations.
 *
 * The two in-plane axes are called first and second axis. They are the remaining
 * world axes in increasing order, e.g. {x, z} for a normal along y.
 *
 * @tparam NORMALAXIS Plane normal as a signed 1-based axis: ±1 for x, ±2 for y, ±3 for z.
 *                    The sign gives the normal direction, i.e -2 is a normal of {0, -1, 0}.
 */
template <int NORMALAXIS = 1>
    requires(NORMALAXIS > -4 && NORMALAXIS < 4 && NORMALAXIS != 0)
class AAFluenceGrid {
public:
    /**
     * @brief Constructs a fluence grid with 1 keV energy bins.
     * @param center     World-space center of the grid in cm.
     * @param spacing    Pixel size in cm along {first axis, second axis}.
     * @param resolution Number of pixels along {first axis, second axis}.
     */
    AAFluenceGrid(const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 2>& spacing = { 1, 1 }, const std::array<std::uint64_t, 2> resolution = { 100, 100 })
        : m_center(center)
    {
        setSpacing(spacing);
        setResolution(resolution);
        setEnergySpacing(1);
    }

    /// @brief Returns the 0-based world axis index of the plane normal.
    static constexpr int normalAxisIndex()
    {
        // std::abs is not constexpr before C++23
        return (NORMALAXIS > 0 ? NORMALAXIS : -NORMALAXIS) - 1;
    }

    /// @brief Returns the sign of the plane normal, +1 or -1.
    static constexpr int normalAxisDirection()
    {
        return NORMALAXIS > 0 ? 1 : -1;
    }

    /// @brief Returns the world axis index of the first in-plane axis.
    static constexpr int firstAxisIndex()
    {
        constexpr auto N = normalAxisIndex();
        if constexpr (N == 0)
            return 1;
        else
            return 0;
    }

    /// @brief Returns the world axis index of the second in-plane axis.
    static constexpr int secondAxisIndex()
    {
        constexpr auto N = normalAxisIndex();
        if constexpr (N != 2)
            return 2;
        else
            return 1;
    }

    /// @brief Returns the plane normal as a unit vector.
    static constexpr std::array<double, 3> planeNormal()
    {
        std::array<double, 3> n = { 0, 0, 0 };
        n[normalAxisIndex()] = normalAxisDirection();
        return n;
    }

    /**
     * @brief Sets the pixel size along {first axis, second axis} and updates the AABB.
     * @param spacing Pixel size in cm; each value is clamped to at least 0.001 cm.
     */
    void setSpacing(const std::array<double, 2>& spacing)
    {
        setSpacing(spacing[0], spacing[1]);
    }

    /**
     * @brief Sets a square pixel size and updates the AABB.
     * @param spacing Pixel size in cm; clamped to at least 0.001 cm.
     */
    void setSpacing(double spacing) { setSpacing(spacing, spacing); }

    /**
     * @brief Sets the pixel size and updates the AABB. Scored counts are kept.
     * @param first_axis  Pixel size in cm along the first axis; clamped to at least 0.001 cm.
     * @param second_axis Pixel size in cm along the second axis; clamped to at least 0.001 cm.
     */
    void setSpacing(double first_axis, double second_axis)
    {
        m_spacing[0] = std::max(first_axis, 0.001);
        m_spacing[1] = std::max(second_axis, 0.001);
        calculateAABB();
    }

    /// @brief Returns the pixel size in cm along {first axis, second axis}.
    const std::array<double, 2>& spacing() const
    {
        return m_spacing;
    }

    /**
     * @brief Sets the number of pixels along {first axis, second axis} and resets all counts.
     * @param res Number of pixels; zero is replaced by one.
     */
    void setResolution(const std::array<std::uint64_t, 2>& res)
    {
        setResolution(res[0], res[1]);
    }

    /**
     * @brief Sets the number of pixels, updates the AABB and resets all counts.
     * @param first_axis  Number of pixels along the first axis; zero is replaced by one.
     * @param second_axis Number of pixels along the second axis; zero is replaced by one.
     */
    void setResolution(std::uint64_t first_axis, std::uint64_t second_axis)
    {
        m_resolution[0] = first_axis > 0 ? first_axis : 1;
        m_resolution[1] = second_axis > 0 ? second_axis : 1;
        calculateAABB();
        gridChanged();
    }

    /// @brief Returns the number of pixels along {first axis, second axis}.
    const std::array<std::uint64_t, 2>& resolution() const
    {
        return m_resolution;
    }

    /**
     * @brief Sets the energy bin width of the per-pixel spectra and resets all counts.
     *        Bins cover [0, MAX_ENERGY()].
     * @param step Bin width in keV; clamped to at least 0.1 keV.
     */
    void setEnergySpacing(double step)
    {
        m_energy_step = std::max(step, 0.1);
        m_energy_resolution = static_cast<std::uint64_t>(MAX_ENERGY() / m_energy_step) + 1;
        gridChanged();
    }

    /// @brief Returns the energy bin width in keV.
    double energySpacing() const
    {
        return m_energy_step;
    }

    /// @brief Returns the bin-center energies in keV of the per-pixel spectra.
    std::vector<double> getSpectrumEnergy() const
    {
        std::vector<double> e(m_energy_resolution);
        const auto start = m_energy_step * 0.5;
        for (std::size_t i = 0; i < m_energy_resolution; ++i)
            e[i] = start + i * m_energy_step;
        return e;
    }

    /**
     * @brief Sets the world-space center of the grid and updates the AABB.
     * @param c New center position in cm.
     */
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        calculateAABB();
    }

    /// @brief Returns the area of a single pixel in cm².
    double pixelArea() const
    {
        return m_spacing[0] * m_spacing[1];
    }

    /**
     * @brief Translates the grid center and AABB by @p dist.
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

    /// @brief Returns the world-space center of the grid.
    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    /// @brief Returns the axis-aligned bounding box of the grid as {xmin, ymin, zmin, xmax, ymax, zmax}.
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /**
     * @brief Returns the raw photon count spectrum of a pixel.
     * @param first_index  Pixel index along the first axis, counted from the low end.
     * @param second_index Pixel index along the second axis, counted from the low end.
     * @return Photon count per energy bin; see getSpectrumEnergy() for bin energies.
     */
    std::vector<std::uint64_t> getSpecter(std::uint64_t first_index, std::uint64_t second_index) const
    {
        const auto index = first_index + m_resolution[0] * second_index;
        return getSpecter(index);
    }

    /**
     * @brief Returns the raw photon count spectrum of a pixel.
     * @param index Flat pixel index, first_index + resolution()[0] * second_index.
     * @return Photon count per energy bin; see getSpectrumEnergy() for bin energies.
     */
    std::vector<std::uint64_t> getSpecter(std::uint64_t index) const
    {
        std::vector<std::uint64_t> s(m_energy_resolution);
        const auto bIdx = index * m_energy_resolution;
        std::copy(m_intensity.cbegin() + bIdx, m_intensity.cbegin() + bIdx + m_energy_resolution, s.begin());
        return s;
    }

    /**
     * @brief Tests a particle ray against the grid plane and rectangle.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the grid plane, or invalid if the
     *         ray misses the grid, is parallel to it or travels along the plane normal.
     */
    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return intersectPlane(p);
    }

    /**
     * @brief Like intersect(), but also returns the plane normal for shading and the
     *        air kerma of the pixel hit as the value, for visualization.
     * @tparam U Scalar type used for the value field.
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
            w.normal = planeNormal();
            constexpr auto X = firstAxisIndex();
            constexpr auto Y = secondAxisIndex();
            const auto x = p.pos[X] + res.intersection * p.dir[X];
            const auto y = p.pos[Y] + res.intersection * p.dir[Y];
            const auto fIdx = binIndex((x - m_aabb[X]) / m_spacing[0], m_resolution[0]);
            const auto sIdx = binIndex((y - m_aabb[Y]) / m_spacing[1], m_resolution[1]);
            w.value = m_dose[fIdx + m_resolution[0] * sIdx].dose();
        }
        return w;
    }

    /**
     * @brief Returns the total energy-score accumulator for the grid.
     * @param index Unused; present for interface consistency.
     */
    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    /// @brief Resets the total energy-score accumulator. Spectrum bin counts are kept.
    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    /**
     * @brief Converts the per-pixel spectra to air kerma and stores it as the pixel dose.
     *
     * For each pixel the kerma is the sum over energy bins of
     * count * E * (mu_tr/rho)_air(E) / pixelArea(), evaluated at the bin-center energy E.
     * The stored dose is replaced, not accumulated, since it is computed from all counts
     * scored since the last clearDoseScored(). No variance is estimated.
     * @param calibration_factor Multiplicative factor applied to the kerma, e.g the beam calibration factor.
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto energy = getSpectrumEnergy();
        const auto air = Material<16>::byNistName("Air, Dry (near sea level)").value();

        std::vector<double> uen(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), uen.begin(), [&](const auto& e) { return air.massEnergyTransferAttenuation(e); });

        for (std::size_t sIdx = 0; sIdx < m_resolution[1]; ++sIdx) {
            for (std::size_t fIdx = 0; fIdx < m_resolution[0]; ++fIdx) {
                double dose = 0;
                const auto index = fIdx + sIdx * m_resolution[0];
                for (std::size_t eIdx = 0; eIdx < m_energy_resolution; ++eIdx) {
                    const auto dindex = index * m_energy_resolution + eIdx;
                    const auto& fluence = m_intensity[dindex];
                    // skip empty bins, uen may be NaN for bin centers below the attenuation tables
                    if (fluence > 0)
                        dose += energy[eIdx] * uen[eIdx] * fluence;
                }
                m_dose[index].set(dose * calibration_factor / pixelArea(), 0, 0);
            }
        }
        return;
    }

    /**
     * @brief Returns the air kerma of a pixel, as computed by addEnergyScoredToDoseScore().
     * @param index Flat pixel index, first_index + resolution()[0] * second_index.
     */
    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose[index];
    }

    /// @brief Resets the energy-score accumulator, all spectrum bin counts and the pixel doses to zero.
    void clearDoseScored()
    {
        m_energyScored.clear();
        std::fill(m_intensity.begin(), m_intensity.end(), std::uint64_t { 0 });
        for (auto& d : m_dose)
            d.clear();
        return;
    }

    /**
     * @brief Records a particle crossing: increments the energy bin counter of the pixel hit
     *        and the total energy-score accumulator, then advances the particle past the grid.
     *        Thread-safe for the bin counter via std::atomic_ref.
     * @param particle Particle arriving at the grid plane; border-translated in place.
     * @param state    Random number generator state (unused).
     */
    void transport(ParticleType auto& particle, RandomState& state)
    {
        // Assuming particle is on the plane
        constexpr auto fD = firstAxisIndex();
        constexpr auto sD = secondAxisIndex();

        const auto eIdx = binIndex(particle.energy / m_energy_step, m_energy_resolution);
        const auto fIdx = binIndex((particle.pos[fD] - m_aabb[fD]) / m_spacing[0], m_resolution[0]);
        const auto sIdx = binIndex((particle.pos[sD] - m_aabb[sD]) / m_spacing[1], m_resolution[1]);

        const auto index = sIdx * (m_energy_resolution * m_resolution[0]) + fIdx * m_energy_resolution + eIdx;
        auto counter = std::atomic_ref(m_intensity[index]);
        counter++;
        m_energyScored.scoreEnergy(particle.energy);
        particle.border_translate(0);
    }

    /// @brief Returns the 32-byte magic identifier used to tag serialized buffers.
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "AAFluenceGrid1";
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
     * @brief Serializes the scorer (geometry, spectrum histograms, energy accumulator and pixel doses)
     *        to a byte vector that can be restored via deserialize().
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_center, buffer);
        Serializer::serialize(m_spacing, buffer);
        Serializer::serialize(m_resolution, buffer);
        Serializer::serialize(m_energy_step, buffer);
        Serializer::serialize(m_energy_resolution, buffer);
        Serializer::serialize(m_intensity, buffer);        
        Serializer::serialize(m_energyScored.energyImparted(), buffer);
        Serializer::serialize(m_energyScored.energyImpartedSquared(), buffer);
        Serializer::serialize(m_energyScored.numberOfEvents(), buffer);
        Serializer::serializeDoseScore(m_dose, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a scorer from a byte buffer produced by serialize().
     * @param buffer Serialized data; the magic ID is expected to have been validated beforehand.
     * @return The reconstructed scorer (always valid for well-formed input).
     */
    static std::optional<AAFluenceGrid<NORMALAXIS>> deserialize(std::span<const char> buffer)
    {
        AAFluenceGrid<NORMALAXIS> item;

        buffer = Serializer::deserialize(item.m_center, buffer);
        buffer = Serializer::deserialize(item.m_spacing, buffer);
        buffer = Serializer::deserialize(item.m_resolution, buffer);
        buffer = Serializer::deserialize(item.m_energy_step, buffer);
        buffer = Serializer::deserialize(item.m_energy_resolution, buffer);
        buffer = Serializer::deserialize(item.m_intensity, buffer);

        double energy, energySq;
        buffer = Serializer::deserialize(energy, buffer);
        buffer = Serializer::deserialize(energySq, buffer);

        std::uint64_t nevents;
        buffer = Serializer::deserialize(nevents, buffer);

        item.m_energyScored.set(energy, energySq, nevents);

        buffer = Serializer::deserializeDoseScore(item.m_dose, buffer);
        
        item.calculateAABB();
        return std::make_optional(item);
    }

protected:
    /// @brief Resizes the count buffer to pixels × energy bins and the dose buffer to pixels, and zeroes both.
    void gridChanged()
    {
        const auto N = m_resolution[0] * m_resolution[1] * m_energy_resolution;
        m_intensity.resize(N);
        std::fill(m_intensity.begin(), m_intensity.end(), std::uint64_t { 0 });

        m_dose.resize(m_resolution[0] * m_resolution[1]);
        for (auto& d : m_dose)
            d.clear();
    }

    /// @brief Truncates a fractional bin coordinate to a bin index clamped to [0, size - 1].
    static inline std::uint64_t binIndex(double coordinate, std::uint64_t size)
    {
        if (coordinate <= 0)
            return 0;
        return std::min(static_cast<std::uint64_t>(coordinate), size - 1);
    }

    /// @brief Recomputes the AABB to bound the grid, padded by GEOMETRIC_ERROR() along the normal axis.
    void calculateAABB()
    {
        constexpr auto X = firstAxisIndex();
        constexpr auto Y = secondAxisIndex();
        constexpr auto Z = normalAxisIndex();
        m_aabb[X] = m_center[X] - m_spacing[0] * m_resolution[0] * 0.5;
        m_aabb[X + 3] = m_center[X] + m_spacing[0] * m_resolution[0] * 0.5;
        m_aabb[Y] = m_center[Y] - m_spacing[1] * m_resolution[1] * 0.5;
        m_aabb[Y + 3] = m_center[Y] + m_spacing[1] * m_resolution[1] * 0.5;
        m_aabb[Z] = m_center[Z] - GEOMETRIC_ERROR();
        m_aabb[Z + 3] = m_center[Z] + GEOMETRIC_ERROR();
    }

    /**
     * @brief Exact ray–plane intersection test against the grid rectangle.
     * @param p Particle whose position and direction define the ray.
     * @return Intersection result with the distance to the grid plane, or invalid if the ray
     *         is nearly parallel to the plane, travels along the normal, or hits outside the grid.
     */
    WorldIntersectionResult intersectPlane(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        constexpr auto N = normalAxisIndex();
        constexpr auto D = normalAxisDirection();
        if (std::abs(p.dir[N]) > GEOMETRIC_ERROR()) {
            if ((p.dir[N] < 0) == (D > 0)) { // only rays that are opposite direction of plane normal
                const auto t = (m_center[N] - p.pos[N]) / p.dir[N];
                if (t > 0.0) { // only forward intersect
                    constexpr int first = firstAxisIndex();
                    const auto pos_f = p.pos[first] + p.dir[first] * t;
                    if (m_aabb[first] < pos_f && pos_f < m_aabb[first + 3]) {
                        constexpr int second = secondAxisIndex();
                        const auto pos_s = p.pos[second] + p.dir[second] * t;
                        if (m_aabb[second] < pos_s && pos_s < m_aabb[second + 3]) {
                            res.intersectionValid = true;
                            res.intersection = t;
                        }
                    }
                }
            }
        }
        return res;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 };
    std::array<double, 2> m_spacing = { 1.0, 1.0 };
    std::array<std::uint64_t, 2> m_resolution = { 100, 100 };
    double m_energy_step = 1;
    std::uint64_t m_energy_resolution = 1;
    std::array<double, 6> m_aabb;
    std::vector<std::uint64_t> m_intensity; // counts, energy bin fastest, then first axis, then second axis
    EnergyScore m_energyScored;
    std::vector<DoseScore> m_dose; // air kerma per pixel
};
}

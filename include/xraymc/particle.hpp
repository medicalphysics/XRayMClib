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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "xraymc/constants.hpp"

#include <array>
#include <concepts>
#include <cstdint>

namespace xraymc {

/**
 * @brief Concept satisfied by any photon particle type usable in the transport loop.
 *
 * A conforming type must expose:
 * - `pos` — `std::array<double, 3>` position in cm.
 * - `dir` — `std::array<double, 3>` unit direction vector.
 * - `energy` — `double` photon energy in keV.
 * - `weight` — `double` statistical weight.
 * - `translate(double)` — advances the particle by the given distance along `dir`.
 * - `border_translate(double)` — same as `translate` but adds a small geometric
 *   offset to ensure the particle clears the boundary.
 *
 * Both `Particle` and `ParticleTrack` satisfy this concept.
 *
 * @tparam U The candidate particle type.
 */
template <typename U>
concept ParticleType = requires(U p, const U pc, std::array<double, 3> vec, double factor) {
    p.translate(factor);
    p.border_translate(factor);
    requires std::same_as<decltype(pc.pos), std::array<double, 3>>;
    requires std::same_as<decltype(pc.dir), std::array<double, 3>>;
    requires std::same_as<decltype(pc.energy), double>;
    requires std::same_as<decltype(pc.weight), double>;
};

/**
 * @brief Lightweight photon particle for standard Monte Carlo transport.
 *
 * Stores the minimal state required by the transport loop: position, direction,
 * energy, and statistical weight. No interaction history is recorded.
 * Satisfies the `ParticleType` concept.
 */
struct Particle {
    std::array<double, 3> pos; ///< Position in cm.
    std::array<double, 3> dir; ///< Unit direction vector.
    double energy; ///< Photon energy in keV.
    double weight; ///< Statistical weight (1.0 for unweighted histories).

    /**
     * @brief Advances the particle along its direction by @p dist cm.
     * @param dist Distance to translate in cm.
     */
    inline void translate(const double dist)
    {
        pos[0] += dir[0] * dist;
        pos[1] += dir[1] * dist;
        pos[2] += dir[2] * dist;
    }

    /**
     * @brief Returns the minimum extra offset added by `border_translate`.
     *
     * Equal to `GEOMETRIC_ERROR<double>()` (1×10⁻⁶ cm), ensuring the particle
     * clears any boundary surface after a step.
     *
     * @return Minimum translation offset in cm.
     */
    inline static constexpr double border_translate_minimum()
    {
        return GEOMETRIC_ERROR<double>();
    }

    /**
     * @brief Advances the particle to (and just beyond) a geometry boundary.
     *
     * Adds `border_translate_minimum()` to @p dist before translating, ensuring
     * the particle is placed on the far side of the boundary and will not
     * re-intersect the same surface on the next step.
     *
     * @param dist Nominal step length to the boundary in cm.
     */
    inline void border_translate(const double dist)
    {
        // Make sure we translate particle beyond any border we want to translate to
        // We simply add 100 nm to the distance, works for float and double
        translate(dist + border_translate_minimum());
    }
};

/**
 * @brief Photon particle that records the last N interaction positions.
 *
 * Extends the basic `Particle` state with a circular buffer of the N most recent
 * positions registered via `registerPosition()`. The interaction samplers call
 * `registerPosition()` before each scatter event when `P` is `ParticleTrack`,
 * allowing post-transport inspection of the particle trajectory.
 * Satisfies the `ParticleType` concept.
 */
struct ParticleTrack {

    std::array<double, 3> pos; ///< Position in cm.
    std::array<double, 3> dir; ///< Unit direction vector.
    double energy; ///< Photon energy in keV.
    double weight; ///< Statistical weight.

    /// @brief Maximum number of positions stored in the history buffer.
    static constexpr std::uint_fast32_t N = 5;

    std::array<std::array<double, 3>, N> m_history; ///< Circular position history buffer.
    std::uint_fast32_t m_index = 0; ///< Total number of positions registered.

    /**
     * @brief Advances the particle along its direction by @p dist cm.
     * @param dist Distance to translate in cm.
     */
    inline void translate(const double dist)
    {
        pos[0] += dir[0] * dist;
        pos[1] += dir[1] * dist;
        pos[2] += dir[2] * dist;
    }

    /**
     * @brief Returns the minimum extra offset added by `border_translate`.
     * @return Minimum translation offset in cm (`GEOMETRIC_ERROR<double>()`).
     */
    inline static constexpr double border_translate_minimum()
    {
        return GEOMETRIC_ERROR<double>();
    }

    /**
     * @brief Advances the particle to (and just beyond) a geometry boundary.
     * @param dist Nominal step length to the boundary in cm.
     */
    inline void border_translate(const double dist)
    {
        // Make sure we translate particle beyond any border we want to translate to
        // We simply add 100 nm to the distance, works for float and double
        translate(dist + border_translate_minimum());
    }

    /**
     * @brief Records the current position in the circular history buffer.
     *
     * Writes `pos` into `m_history[m_index % N]` and increments `m_index`.
     * Old entries are overwritten once more than N positions have been recorded.
     */
    void registerPosition()
    {
        // incrementing and circle around
        m_history[m_index % N] = pos;
        m_index++;
    }

    /**
     * @brief Returns the stored positions in chronological order.
     *
     * Reconstructs the correct temporal order from the circular buffer. If fewer
     * than N positions have been recorded, only those positions are meaningful;
     * call `getSize()` to determine the valid count.
     *
     * @return Array of N positions; entries beyond `getSize()` are uninitialised.
     */
    std::array<std::array<double, 3>, N> getHistory() const
    {
        std::array<std::array<double, 3>, N> r;
        if (m_index <= N) {
            for (std::uint_fast32_t i = 0; i < m_index; ++i)
                r[i] = m_history[i];
        } else {
            const std::uint_fast32_t first = m_index % N;
            for (std::uint_fast32_t i = 0; i < N; ++i)
                r[i] = m_history[(first + i) % N];
        }
        return r;
    }

    /**
     * @brief Returns the number of valid entries in the history buffer.
     * @return min(total registrations, N).
     */
    std::uint_fast32_t getSize() const
    {
        return std::min(N, m_index);
    }
};

}
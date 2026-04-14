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

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "xraymc/particle.hpp"

#include <array>
#include <atomic>
#include <vector>

/**
 * @brief Thread-safe collector of photon trajectory points across multiple transport histories.
 *
 * Maintains a pre-allocated flat buffer of `TrackPoint` entries. Each entry stores
 * a particle ID and a 3D position. When `registerParticle()` is called from
 * concurrent transport threads, an `std::atomic_ref` fetch-add on the write index
 * ensures each thread reserves a contiguous, non-overlapping slice of the buffer
 * without locking.
 *
 * Points for a specific particle can be retrieved in insertion order via `track()`.
 * If the buffer is full, further calls to `registerParticle()` are silently ignored.
 */
class ParticleTracker {
public:
    /**
     * @brief Constructs the tracker with a pre-allocated buffer capacity.
     * @param size Maximum number of `TrackPoint` entries to store (default 1024).
     */
    ParticleTracker(std::size_t size = 1024)
    {
        m_points.resize(size);
    }

    /// @brief Default equality comparison (compares buffer contents and indices).
    bool operator==(const ParticleTracker&) const = default;

    /**
     * @brief Resizes the internal point buffer.
     *
     * Replaces the current buffer with one of @p size default-constructed entries.
     * Should be called before transport begins; not safe to call concurrently with
     * `registerParticle()`.
     *
     * @param size New buffer capacity in number of `TrackPoint` entries.
     */
    void setNumberOfPoints(std::size_t size)
    {
        m_points.resize(size);
    }

    /**
     * @brief Records the interaction history of a `ParticleTrack` into the buffer.
     *
     * Atomically reserves `p.getSize() + 1` slots (history positions plus the
     * current position) and writes them with a unique particle ID. If the buffer
     * does not have enough remaining capacity, the particle is silently dropped.
     *
     * Safe to call from multiple threads simultaneously.
     *
     * @param p The particle whose history and current position are to be recorded.
     */
    void registerParticle(const ParticleTrack& p)
    {
        // Threadsafe particle register
        const auto tracksize = p.getSize() + 1;

        auto ain = std::atomic_ref(m_index);
        const auto currentindex = ain.fetch_add(tracksize);

        if ((currentindex + tracksize) < m_points.size()) {
            // we have space for particle
            auto aid = std::atomic_ref(m_currentId);
            const auto id = aid.fetch_add(std::uint64_t { 1 });
            const auto track = p.getHistory();
            for (std::size_t i = 0; i < tracksize - 1; ++i)
                m_points[currentindex + i] = { .particleID = id, .position = track[i] };
            // adding current position
            m_points[currentindex + tracksize - 1] = { .particleID = id, .position = p.pos };
        }
    }

    /**
     * @brief Returns the number of particles successfully registered so far.
     * @return Count of particles whose histories have been stored in the buffer.
     */
    std::uint64_t numberOfParticles() const
    {
        return m_currentId - 1;
    }

    /**
     * @brief Retrieves the stored positions for a specific particle in order.
     *
     * Scans the entire buffer and collects all entries whose particle ID matches
     * @p particleNumber (0-based). Returns them in the order they were written,
     * which corresponds to the sequence of interaction positions followed by the
     * particle's final position.
     *
     * @param particleNumber Zero-based index of the particle to retrieve.
     * @return Vector of 3D positions (in cm) for the requested particle.
     */
    std::vector<std::array<double, 3>> track(std::uint64_t particleNumber) const
    {
        const auto id = particleNumber + 1;
        std::vector<std::array<double, 3>> res;
        for (const auto& el : m_points)
            if (el.particleID == id)
                res.push_back(el.position);
        return res;
    }

private:
    /// @brief A single recorded point associated with a particle history.
    struct TrackPoint {
        std::uint64_t particleID = 0;      ///< 1-based ID of the owning particle.
        std::array<double, 3> position;    ///< Position in cm.
        bool operator==(const TrackPoint&) const = default;
    };
    std::vector<TrackPoint> m_points; ///< Pre-allocated flat point buffer.
    std::size_t m_index = 0;          ///< Next free write index (updated atomically).
    std::uint64_t m_currentId = 1;    ///< Next particle ID to assign (updated atomically).
};
}
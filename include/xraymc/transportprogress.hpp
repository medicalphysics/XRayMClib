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

#include <atomic>
#include <chrono>
#include <string>
#include <utility>

/**
 * @brief Thread-safe progress tracker for a Monte Carlo transport run.
 *
 * Tracks the number of simulated particles, elapsed wall-clock time, and an
 * early-termination flag. Worker threads call `addCompletedNumber()` after each
 * exposure; the UI or controlling thread polls `message()`, `progress()`, or
 * `continueSimulation()`.
 *
 * Counter updates use `std::atomic_ref` with relaxed ordering for low overhead.
 * The stop flag uses a sequentially-consistent store/load so that a `setStopSimulation()`
 * call from any thread is visible to all workers on their next `continueSimulation()` check.
 */
class TransportProgress {
public:
    /// @brief Default constructor — object must be initialised with `start()` before use.
    TransportProgress() { }

    /**
     * @brief Constructs and immediately starts the progress tracker.
     * @param n_particles Total number of particles to simulate.
     */
    TransportProgress(std::uint64_t n_particles)
    {
        start(n_particles);
    }

    /**
     * @brief (Re-)initialises the tracker for a new simulation run.
     *
     * Resets the particle counter and wall-clock start time, and re-enables
     * the continue flag. Safe to call between successive `Transport::run` calls.
     *
     * @param N Total number of particles to simulate (clamped to ≥ 1).
     */
    void start(std::uint64_t N)
    {
        m_continue_simulation_flag = true;
        m_nParticles = std::max(N, std::uint64_t { 1 });
        m_nParticleCount = 0;
        m_start = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Records that @p N additional particles have been transported.
     *
     * Atomically increments the completed-particle counter and updates the
     * elapsed-time snapshot used for ETA estimation. Records the end time when
     * the total reaches `m_nParticles`. Safe to call concurrently from multiple threads.
     *
     * @param N Number of particles completed in this batch.
     */
    void addCompletedNumber(std::uint64_t N)
    {
        auto aref = std::atomic_ref(m_nParticleCount);
        aref.fetch_add(N, std::memory_order_relaxed);

        auto dref = std::atomic_ref(m_elapsed);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_start);
        dref.exchange(duration);

        if (aref.load() >= m_nParticles) {
            m_end = std::chrono::high_resolution_clock::now();
        }
    }

    /**
     * @brief Returns whether transport should continue.
     *
     * Worker threads check this flag after each exposure. Returns false after
     * `setStopSimulation()` has been called.
     *
     * @return True if transport should continue, false if it has been cancelled.
     */
    bool continueSimulation() const
    {
        return m_continue_simulation_flag;
    }

    /**
     * @brief Returns the current and total particle counts.
     * @return A pair of {completed, total} particle counts.
     */
    std::pair<std::uint64_t, std::uint64_t> progress() const
    {
        return std::make_pair(m_nParticleCount, m_nParticles);
    }

    /**
     * @brief Requests early termination of the simulation.
     *
     * Sets the continue flag to false via `std::atomic_ref`. Worker threads
     * will observe this on their next `continueSimulation()` check and exit.
     */
    void setStopSimulation()
    {
        auto flag = std::atomic_ref(m_continue_simulation_flag);
        flag.store(false);
    }

    /**
     * @brief Returns the total simulation wall-clock time as a human-readable string.
     *
     * Returns the duration between `start()` and completion if all particles have
     * been transported; otherwise returns "Not done yet".
     *
     * @return Human-readable elapsed time string, or "Not done yet".
     */
    std::string humanTotalTime() const
    {
        if (m_nParticleCount == m_nParticles) {
            return human_time(std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start));
        }
        return "Not done yet";
    }

    /**
     * @brief Returns the total simulation wall-clock time in milliseconds.
     *
     * Only meaningful after all particles have been transported (i.e. after
     * `addCompletedNumber` has been called with the final batch).
     *
     * @return Duration from `start()` to completion.
     */
    std::chrono::milliseconds totalTime() const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start);
    }

    /**
     * @brief Generates a one-line progress message suitable for console display.
     *
     * Format: "<elapsed>, remaining <ETA> [<percent>%]"
     * ETA is estimated by linear extrapolation from the elapsed time and the
     * fraction of particles completed. Returns "NA [0%]" if no particles have
     * been completed yet.
     *
     * @return Progress string with elapsed time, estimated remaining time, and percentage.
     */
    std::string message() const
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_start);

        std::string message = human_time(elapsed) + ", remaining ";

        if (m_nParticleCount > 0) {
            const auto p = std::to_string((m_nParticleCount * 100) / m_nParticles);
            const auto remaining = (m_elapsed * (m_nParticles - m_nParticleCount)) / m_nParticleCount;
            message += human_time(remaining) + " [" + p + "%]";
        } else {
            message += "NA [0%]";
        }
        return message;
    }

    /**
     * @brief Formats a millisecond duration as a human-readable string.
     *
     * Selects the most appropriate unit:
     * - Hours if duration > 3 hours.
     * - Minutes if duration > 3 minutes.
     * - Seconds otherwise.
     *
     * @param time Duration to format.
     * @return String such as "42 sec", "7 min", or "2 hrs".
     */
    static std::string human_time(const std::chrono::milliseconds& time)
    {
        if (time > std::chrono::hours(3))
            return std::to_string(std::chrono::duration_cast<std::chrono::hours>(time).count()) + " hrs";
        else if (time > std::chrono::minutes(3))
            return std::to_string(std::chrono::duration_cast<std::chrono::minutes>(time).count()) + " min";
        else
            return std::to_string(std::chrono::duration_cast<std::chrono::seconds>(time).count()) + " sec";
    }

private:
    std::uint64_t m_nParticles = 1;          ///< Total number of particles to simulate.
    std::uint64_t m_nParticleCount = 0;       ///< Number of particles completed so far.
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start; ///< Wall-clock start time.
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;   ///< Wall-clock end time.
    std::chrono::milliseconds m_elapsed;      ///< Last-recorded elapsed time snapshot for ETA calculation.
    bool m_continue_simulation_flag = true;   ///< False after `setStopSimulation()` is called.
};
}

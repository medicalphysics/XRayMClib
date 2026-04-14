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

#include "xraymc/beams/beamtype.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <iostream>
#include <thread>

namespace xraymc {

/**
 * @brief Multi-threaded Monte Carlo photon transport driver.
 *
 * Distributes beam exposures across a configurable number of worker threads using
 * an atomic fetch-add counter. Each thread draws exposures, samples particles from
 * each exposure, and calls `world.transport()` per particle.
 *
 * After all threads complete, the beam's calibration factor is applied to convert
 * the accumulated energy scores to dose (unless `useBeamCalibration` is false).
 *
 * Typical usage:
 * @code
 * Transport transport;
 * transport(world, beam);           // operator()
 * Transport::run(world, beam, 4);   // static, 4 threads
 * Transport::runConsole(world, beam); // prints progress to stdout
 * @endcode
 */
class Transport {
public:
    /**
     * @brief Constructs a Transport object, defaulting to all available hardware threads.
     */
    Transport()
    {
        const std::uint64_t hvc = std::thread::hardware_concurrency();
        m_nThreads = std::max(hvc, std::uint64_t { 1 });
    }

    /// @brief Returns the number of worker threads that will be used by `operator()`.
    std::uint64_t numberOfThreads() const { return m_nThreads; }

    /**
     * @brief Sets the number of worker threads for subsequent `operator()` calls.
     * @param n Desired thread count; clamped to a minimum of 1.
     */
    void setNumberOfThreads(std::uint64_t n) { m_nThreads = std::max(n, std::uint64_t { 0 }); }

    /**
     * @brief Runs transport using the instance's thread count setting.
     *
     * Delegates to `Transport::run` with the stored `m_nThreads`.
     *
     * @tparam B   Beam type satisfying `BeamType`.
     * @tparam Ws  World item types.
     * @param world              The simulation world; energy scores are accumulated in place.
     * @param beam               Source of exposures and particles.
     * @param progress           Optional progress tracker (may be nullptr).
     * @param useBeamCalibration If true, applies `beam.calibrationFactor()` to convert
     *                           energy scores to dose after transport.
     */
    template <BeamType B, WorldItemType... Ws>
    auto operator()(World<Ws...>& world, const B& beam, TransportProgress* progress = nullptr, bool useBeamCalibration = true) const
    {
        return run(world, beam, m_nThreads, progress, useBeamCalibration);
    }

    /**
     * @brief Runs transport and prints periodic progress messages to stdout.
     *
     * Launches transport in a background thread and polls `TransportProgress::message()`
     * every @p update_ms milliseconds, overwriting the previous line with `\r`.
     *
     * @tparam B   Beam type satisfying `BeamType`.
     * @tparam Ws  World item types.
     * @param world              The simulation world.
     * @param beam               Source of exposures and particles.
     * @param nThreads           Number of worker threads (0 = hardware concurrency).
     * @param useBeamCalibration If true, applies beam calibration after transport.
     * @param update_ms          Console refresh interval in milliseconds (default 2000).
     * @return Total simulation wall-clock time as reported by `TransportProgress`.
     */
    template <BeamType B, WorldItemType... Ws>
    static auto runConsole(World<Ws...>& world, const B& beam, std::uint64_t nThreads = 0, bool useBeamCalibration = true, std::uint32_t update_ms = 2000)
    {
        xraymc::TransportProgress progress;

        bool running = true;
        if (nThreads == 0)
            nThreads = std::thread::hardware_concurrency();

        std::thread job([&]() {
            Transport::run<B, Ws...>(world, beam, nThreads, &progress, useBeamCalibration);
            running = false;
        });
        std::string message;
        while (running) {
            std::this_thread::sleep_for(std::chrono::milliseconds(update_ms));
            std::cout << std::string(message.length(), ' ') << "\r";
            message = progress.message();
            std::cout << message << std::flush << "\r";
        }
        job.join();
        std::cout << std::string(message.length(), ' ') << "\r";
        return progress.totalTime();
    }

    /**
     * @brief Runs multi-threaded Monte Carlo transport.
     *
     * Clears all energy scores, spawns `nThreads − 1` `std::jthread` workers plus
     * one call on the calling thread, and waits for all to finish. Exposures are
     * distributed with an `std::atomic<uint64_t>` fetch-add counter so no two
     * threads process the same exposure.
     *
     * If `useBeamCalibration` is true and the simulation was not cancelled via
     * @p progress, `beam.calibrationFactor()` is applied and energy scores are
     * converted to dose via `world.addEnergyScoredToDoseScore()`.
     *
     * @tparam B   Beam type satisfying `BeamType`.
     * @tparam Ws  World item types.
     * @param world              The simulation world; modified in place.
     * @param beam               Source of exposures and particles.
     * @param nThreads           Number of worker threads (clamped to ≥ 1).
     * @param progress           Optional progress tracker (may be nullptr).
     * @param useBeamCalibration If true, applies beam calibration after transport.
     */
    template <BeamType B, WorldItemType... Ws>
    static void run(World<Ws...>& world, const B& beam, std::uint64_t nThreads = 1, TransportProgress* progress = nullptr, bool useBeamCalibration = true)
    {
        // clearing scored energy before run
        world.clearEnergyScored();

        nThreads = std::max(nThreads, std::uint64_t { 1 });

        const auto nExposures = beam.numberOfExposures();
        std::vector<std::jthread> threads;
        threads.reserve(nThreads - 1);
        std::vector<RandomState> states(nThreads - 1);
        std::atomic<std::uint64_t> start(0);

        if (progress)
            progress->start(beam.numberOfParticles());

        for (std::size_t i = 0; i < nThreads - 1; ++i) {
            threads.emplace_back(Transport::template runWorker<B, Ws...>, std::ref(world), std::cref(beam), std::ref(states[i]), std::ref(start), nExposures, progress);
        }
        RandomState state;
        runWorker(world, beam, state, start, nExposures, progress);

        for (auto& thread : threads) {
            thread.join();
        }
        bool cont = progress ? progress->continueSimulation() : true;
        if (useBeamCalibration && cont) {
            const auto beamCalibrationFactor = beam.calibrationFactor(progress);
            world.addEnergyScoredToDoseScore(beamCalibrationFactor);
            world.clearEnergyScored();
        }
    }

protected:
    /**
     * @brief Per-thread worker loop that processes beam exposures until all are done.
     *
     * Atomically claims the next exposure index via `exposureStart.fetch_add(1)` and
     * transports all particles in that exposure. Continues until `exposureEnd` is
     * reached or `progress->continueSimulation()` returns false.
     *
     * @tparam B   Beam type satisfying `BeamType`.
     * @tparam Ws  World item types.
     * @param world         The simulation world; `transport()` is called per particle.
     * @param beam          Source of exposures and particles.
     * @param state         Per-thread PRNG state.
     * @param exposureStart Shared atomic counter of the next unstarted exposure index.
     * @param exposureEnd   Total number of exposures (exclusive upper bound).
     * @param progress      Optional progress tracker; signals early termination if set.
     */
    template <BeamType B, WorldItemType... Ws>
    static void runWorker(World<Ws...>& world, const B& beam, RandomState& state, std::atomic<std::uint64_t>& exposureStart, std::uint64_t exposureEnd, TransportProgress* progress = nullptr)
    {
        auto n = exposureStart.fetch_add(1);
        while (n < exposureEnd) {
            const auto exposure = beam.exposure(n);
            const auto nParticles = exposure.numberOfParticles();
            for (std::uint64_t i = 0; i != nParticles; ++i) {
                auto particle = exposure.sampleParticle(state);
                world.transport(particle, state);
            }
            n = exposureStart.fetch_add(1);
            if (progress) {
                progress->addCompletedNumber(nParticles);
                if (!progress->continueSimulation()) {
                    n = exposureEnd; // we stop simulation
                }
            }
        }
    }

private:
    std::uint64_t m_nThreads = 1;
};

}

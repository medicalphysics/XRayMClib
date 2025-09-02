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

class Transport {
public:
    Transport()
    {
        const std::uint64_t hvc = std::thread::hardware_concurrency();
        m_nThreads = std::max(hvc, std::uint64_t { 1 });
    }

    std::uint64_t numberOfThreads() const { return m_nThreads; }
    void setNumberOfThreads(std::uint64_t n) { m_nThreads = std::max(n, std::uint64_t { 1 }); }

    template <BeamType B, WorldItemType... Ws>
    auto operator()(World<Ws...>& world, const B& beam, TransportProgress* progress = nullptr, bool useBeamCalibration = true) const
    {
        return run(world, beam, m_nThreads, progress, useBeamCalibration);
    }

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
        if (useBeamCalibration) {
            const auto beamCalibrationFactor = beam.calibrationFactor(progress);
            world.addEnergyScoredToDoseScore(beamCalibrationFactor);
            world.clearEnergyScored();
        }
    }

protected:
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

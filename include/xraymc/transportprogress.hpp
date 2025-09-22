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

namespace xraymc {
class TransportProgress {
public:
    TransportProgress() { }
    TransportProgress(std::uint64_t n_particles)
    {
        start(n_particles);
    }

    void start(std::uint64_t N)
    {
        m_continue_simulation_flag = true;
        m_nParticles = std::max(N, std::uint64_t { 1 });
        m_nParticleCount = 0;
        m_start = std::chrono::high_resolution_clock::now();
    }

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

    bool continueSimulation() const
    {
        return m_continue_simulation_flag;
    }

    std::pair<std::uint64_t, std::uint64_t> progress() const
    {
        return std::make_pair(m_nParticleCount, m_nParticles);
    }

    void setStopSimulation()
    {
        auto flag = std::atomic_ref(m_continue_simulation_flag);
        flag.store(false);
    }

    std::string humanTotalTime() const
    {
        if (m_nParticleCount == m_nParticles) {
            return human_time(std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start));
        }
        return "Not done yet";
    }

    std::chrono::milliseconds totalTime() const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start);
    }

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
    std::uint64_t m_nParticles = 1;
    std::uint64_t m_nParticleCount = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
    std::chrono::milliseconds m_elapsed;
    bool m_continue_simulation_flag = true;
};
}

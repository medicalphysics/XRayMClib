/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "dxmc/particle.hpp"

#include <array>
#include <atomic>
#include <vector>

namespace dxmc {
class ParticleTracker {
public:
    ParticleTracker(std::size_t size = 1024)
    {
        m_points.resize(size);
    }

    void setNumberOfPoints(std::size_t size)
    {
        m_points.resize(size);
    }

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

    std::uint64_t numberOfParticles() const
    {
        return m_currentId - 1;
    }

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
    struct TrackPoint {
        std::uint64_t particleID = 0;
        std::array<double, 3> position;
    };
    std::vector<TrackPoint> m_points;
    std::size_t m_index = 0;
    std::uint64_t m_currentId = 1;
};
}
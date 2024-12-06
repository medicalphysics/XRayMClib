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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {
template <bool ENABLETRACKING = false>
class IsotropicMonoEnergyBeamCircleExposure {
public:
    IsotropicMonoEnergyBeamCircleExposure(const std::array<double, 3>& pos, const std::array<double, 3>& direction, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_dir(direction)
        , m_NParticles(N)
    {

        m_dirCosines = calculateDirectionCosines(m_dir);
    }

    const std::array<double, 3>& position() const { return m_pos; }

    void setCollimationHalfAngle(double angle)
    {
        m_r = std::tan(std::abs(angle));
    }

    void setEnergy(double e)
    {
        m_energy = e;
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {

        auto k_x = state.randomUniform() * 2 - 1;
        auto k_y = state.randomUniform() * 2 - 1;
        while (k_x * k_x + k_y * k_y > 1) {
            k_x = state.randomUniform() * 2 - 1;
            k_y = state.randomUniform() * 2 - 1;
        }

        const auto x = m_r * k_x;
        const auto y = m_r * k_y;

        const std::array<double, 3> dir = {
            m_dir[0] + x * m_dirCosines[0][0] + y * m_dirCosines[1][0],
            m_dir[1] + x * m_dirCosines[0][1] + y * m_dirCosines[1][1],
            m_dir[2] + x * m_dirCosines[0][2] + y * m_dirCosines[1][2]
        };

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = vectormath::normalized(dir),
                .energy = m_energy,
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = vectormath::normalized(dir),
                .energy = m_energy,
                .weight = 1
            };
            return p;
        }
    }

protected:
    static std::array<std::array<double, 3>, 2> calculateDirectionCosines(const std::array<double, 3>& dir)
    {
        std::array<std::array<double, 3>, 2> cos = { { { 0, 0, 0 }, { 0, 0, 0 } } };
        const auto minInd = vectormath::argmin3<std::uint_fast32_t, double>(dir);
        std::array<double, 3> k { 0, 0, 0 };
        k[minInd] = 1;

        const auto vec_xy_raw = vectormath::cross(dir, k);
        cos[0] = vectormath::normalized(vec_xy_raw);
        cos[1] = vectormath::cross(cos[0], dir);
        return cos;
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    double m_r = 0;
    std::uint64_t m_NParticles = 100;
    double m_energy = 60;
};

template <bool ENABLETRACKING = false>
class IsotropicMonoEnergyBeamCircle {
public:
    IsotropicMonoEnergyBeamCircle(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 1, 0, 0 })
        : m_pos(pos)
        , m_dir(dir)
    {
        vectormath::normalize(m_dir);
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }
    const std::array<double, 3>& position() const { return m_pos; }

    const std::array<double, 3>& direction() const
    {
        return m_dir;
    }

    void setDirection(const std::array<double, 3>& dir)
    {
        m_dir = vectormath::normalized(dir);
    }

    void setEnergy(double energy)
    {
        m_energy = std::clamp(std::abs(energy), MIN_ENERGY(), MAX_ENERGY());
    }

    void setCollimationHalfAngle(double ang)
    {
        m_collimationHalfAngle = std::clamp(std::abs(ang), 0.0, PI_VAL());
    }

    double collimationHalfAngle() const { return m_collimationHalfAngle; }

    IsotropicMonoEnergyBeamCircleExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicMonoEnergyBeamCircleExposure<ENABLETRACKING> exp(m_pos, m_dir, m_particlesPerExposure);
        exp.setCollimationHalfAngle(m_collimationHalfAngle);
        exp.setEnergy(m_energy);
        return exp;
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const noexcept
    {
        return 1;
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 1, 0, 0 };
    double m_collimationHalfAngle = 0;
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
    double m_energy = 60;
};
}
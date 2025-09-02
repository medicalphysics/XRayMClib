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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "xraymc/floating.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {
template <bool ENABLETRACKING = false>
class PencilBeamExposure {
public:
    PencilBeamExposure(const std::array<double, 3>& pos, const std::array<double, 3>& dir, double energy, double weight, std::uint64_t N)
        : m_energy(energy)
        , m_weight(weight)
        , m_pos(pos)
        , m_dir(dir)
        , m_NParticles(N)
    {
    }

    const std::array<double, 3>& position() const { return m_pos; }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {
        if constexpr (ENABLETRACKING) {
            ParticleTrack p = { .pos = m_pos,
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            p.registerPosition();
            return p;
        } else {
            Particle p = { .pos = m_pos,
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            return p;
        }
    }

private:
    double m_energy = 60;
    double m_weight = 1;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::uint64_t m_NParticles = 100;
};

template <bool ENABLETRACKING = false>
class PencilBeam {
public:
    PencilBeam(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 }, double energy = 60)
        : m_energy(energy)
        , m_pos(pos)
        , m_dir(dir)
    {
        vectormath::normalize(m_dir);
    }

    void setEnergy(double energy)
    {
        m_energy = std::abs(energy);
    }

    double energy() const { return m_energy; }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }

    void setPosition(const std::array<double, 3>& pos)
    {
        m_pos = pos;
    }

    const std::array<double, 3>& position() const
    {
        return m_pos;
    }

    void setDirection(const std::array<double, 3>& dir)
    {
        m_dir = dir;
        vectormath::normalize(m_dir);
    }

    const std::array<double, 3>& direction() const
    {
        return m_dir;
    }

    std::array<std::array<double, 3>, 2> directionCosines() const
    {
        std::array<double, 3> cand = { 0, 0, 0 };
        const auto minIdx = vectormath::argmin3(m_dir);
        cand[minIdx] = 1;

        std::array<std::array<double, 3>, 2> cos;
        cos[0] = vectormath::cross(m_dir, cand);
        vectormath::normalize(cos[0]);
        cos[1] = vectormath::cross(m_dir, cos[0]);
        return cos;
    }

    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dir = vectormath::cross(dir[0], dir[1]);
        vectormath::normalize(m_dir);
    }
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dir = vectormath::cross(xdir, ydir);
        vectormath::normalize(m_dir);
    }

    const std::array<double, 4> collimationAngles() const
    {
        return std::array { 0.0, 0.0, 0.0, 0.0 };
    }

    void setParticleWeight(double weight = 1)
    {
        m_weight = weight;
    }

    PencilBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        PencilBeamExposure<ENABLETRACKING> exp(m_pos, m_dir, m_energy, m_weight, m_particlesPerExposure);
        return exp;
    }

    void setAirKerma(double k)
    {
        m_airKerma = std::abs(k);
    }

    double airKerma() const { return m_airKerma; }

    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        auto air_cand = Material<5>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();
        const double kerma = numberOfParticles() * air.massEnergyTransferAttenuation(m_energy);
        return m_airKerma / kerma;
    }

private:
    double m_energy = 60;
    double m_weight = 1;
    double m_airKerma = 1;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
};
}
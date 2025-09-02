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

#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {
template <bool ENABLETRACKING = false>
class IsotropicCircularMonoEnergyBeamExposure {
public:
    IsotropicCircularMonoEnergyBeamExposure(const std::array<double, 3>& pos, double energy, double radius = 0, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_energy(energy)
        , m_radius(radius)
        , m_NParticles(N)
    {
    }

    const std::array<double, 3>& position() const { return m_pos; }

    void setCollimationHalfAngles(const std::array<double, 4>& angles)
    {
        m_directionSampler.setData(angles);
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {
        const auto angle = state.randomUniform(2 * PI_VAL());
        const auto sang = std::sin(angle);
        const auto cang = std::cos(angle);
        const auto X = vectormath::rotate({ 0, 1, 0 }, { 0, 0, 1 }, sang, cang);
        constexpr std::array<double, 3> Y = { 0, 0, 1 };
        const auto Z = vectormath::cross(X, Y);
        const auto pos = vectormath::add(vectormath::scale(Z, -m_radius), m_pos);

        const auto dirz = m_directionSampler(state);
        const auto dir = vectormath::changeBasis(X, Y, Z, dirz);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = pos,
                .dir = dir,
                .energy = m_energy,
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = pos,
                .dir = dir,
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
    SphereSamplingRectangularField m_directionSampler;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::uint64_t m_NParticles = 100;
    double m_radius = 0;
    double m_energy = 60;
};

template <bool ENABLETRACKING = false>
class IsotropicCircularMonoEnergyBeam {
public:
    IsotropicCircularMonoEnergyBeam(const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_pos(pos)
    {
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }
    const std::array<double, 3>& position() const { return m_pos; }

    void setRadius(double r)
    {
        m_radius = std::abs(r);
    }
    double radius() const { return m_radius; }

    void setEnergy(double energy) { m_energy = std::min(std::max(MIN_ENERGY(), energy), MAX_ENERGY()); }
    double energy() const { return m_energy; }

    const std::array<double, 4>& collimationHalfAngles() const { return m_collimationHalfAngles; }

    void setCollimationHalfAngles(const std::array<double, 4>& angles) { m_collimationHalfAngles = angles; }
    void setCollimationHalfAngles(double minX, double minY, double maxX, double maxY)
    {
        m_collimationHalfAngles[0] = minX;
        m_collimationHalfAngles[1] = minY;
        m_collimationHalfAngles[2] = maxX;
        m_collimationHalfAngles[3] = maxY;
    }

    IsotropicCircularMonoEnergyBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicCircularMonoEnergyBeamExposure<ENABLETRACKING> exp(m_pos, m_energy, m_radius, m_particlesPerExposure);
        exp.setCollimationHalfAngles(m_collimationHalfAngles);
        return exp;
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const noexcept
    {
        return 1;
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 4> m_collimationHalfAngles = { 0, 0, 0, 0 };
    double m_radius = 0;
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
    double m_energy = 60;
};
}
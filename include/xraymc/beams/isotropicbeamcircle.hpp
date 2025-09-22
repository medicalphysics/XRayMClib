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

#include "xraymc/beams/utilities/spheresamplingcircularfield.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {
template <bool ENABLETRACKING = false>
class IsotropicBeamCircleExposure {
public:
    IsotropicBeamCircleExposure(const std::array<double, 3>& pos, const std::array<double, 3>& direction, std::uint64_t N = 1E6)
        : m_pos(pos)
        , m_dir(direction)
        , m_NParticles(N)
    {

        m_dirCosines = calculateDirectionCosines(m_dir);
    }

    const std::array<double, 3>& position() const { return m_pos; }

    void setCollimationHalfAngle(double angle)
    {
        m_directionSampler.setData(angle);
    }

    void setSpecterDistribution(const SpecterDistribution<double>& s)
    {
        m_specterDist = s;
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {

        const auto dirz = m_directionSampler(state);
        const auto dir = vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dirz);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = dir,
                .energy = m_specterDist.sampleValue(state),
                .weight = 1
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = dir,
                .energy = m_specterDist.sampleValue(state),
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
    SphereSamplingCircularField m_directionSampler;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::uint64_t m_NParticles = 100;
    SpecterDistribution<double> m_specterDist;
};

template <bool ENABLETRACKING = false>
class IsotropicBeamCircle {
public:
    IsotropicBeamCircle(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& dir = { 1, 0, 0 })
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

    void setEnergySpecter(const std::vector<std::pair<double, double>>& specter)
    {
        m_specter = SpecterDistribution(specter);
    }

    void setCollimationHalfAngle(double ang)
    {
        m_collimationHalfAngle = std::clamp(std::abs(ang), 0.0, PI_VAL());
    }

    double collimationHalfAngle() const { return m_collimationHalfAngle; }

    IsotropicBeamCircleExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        IsotropicBeamCircleExposure<ENABLETRACKING> exp(m_pos, m_dir, m_particlesPerExposure);
        exp.setCollimationHalfAngle(m_collimationHalfAngle);
        exp.setSpecterDistribution(m_specter);
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
    SpecterDistribution<double> m_specter;
};
}
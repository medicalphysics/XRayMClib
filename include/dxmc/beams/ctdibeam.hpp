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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/beams/filters/bowtiefilter.hpp"
#include "dxmc/beams/filters/ctorganaecfilter.hpp"
#include "dxmc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {
template <bool ENABLETRACKING = false>
class CTDIBeamExposure {
public:
    CTDIBeamExposure(double angle, double SDD, std::uint64_t historiesPerExposure,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double>* specter, const BowtieFilter* bowtie, double weight = 1)
        : m_directionSampler(collimationHalfAngles)
        , m_Nparticles(historiesPerExposure)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_weight(weight)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
    {
        m_dirCosineX = vectormath::rotate(m_dirCosineX, m_dirCosineY, angle);
        m_dir = vectormath::cross(m_dirCosineX, m_dirCosineY);
        m_pos = vectormath::scale(m_dir, -SDD / 2);
    }

    CTDIBeamExposure() = delete;

    const std::array<double, 3>& position() const { return m_pos; }

    std::uint64_t numberOfParticles() const
    {
        return m_Nparticles;
    }

    auto sampleParticle(RandomState& state) const noexcept
    {
        const auto dir = m_directionSampler(state);
        const auto angx = std::atan(dir[0] / dir[2]);
        const auto bowtie_weight = m_bowtieFilter->operator()(angx);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = vectormath::changeBasis(m_dirCosineX, m_dirCosineY, m_dir, dir),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight * bowtie_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = vectormath::changeBasis(m_dirCosineX, m_dirCosineY, m_dir, dir),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight * bowtie_weight
            };
            return p;
        }
    }

protected:
private:
    SphereSamplingRectangularField m_directionSampler;
    std::uint64_t m_Nparticles = 1;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 1, 0 };
    std::array<double, 3> m_dirCosineX = { 1, 0, 0 };
    std::array<double, 3> m_dirCosineY = { 0, 0, 1 };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    double m_weight = 1;
    const SpecterDistribution<double>* m_specter = nullptr;
    const BowtieFilter* m_bowtieFilter = nullptr;
};

template <bool ENABLETRACKING = false>
class CTDIBeam {
public:
    CTDIBeam(double angleStep, double SDD, const std::array<double, 2>& collimationHalfAngles, std::uint64_t particlesPerExposure, const SpecterDistribution<double>& specter, const BowtieFilter& bowtie, const CTOrganAECFilter& organFilter, double weight = 1)
        : m_angleStep(angleStep)
        , m_sdd(SDD)
        , m_weight(weight)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_particlesPerExposure(particlesPerExposure)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
        , m_organFilter(organFilter)
    {
    }

    CTDIBeam(double angleStep, double SDD, const std::array<double, 2>& collimationHalfAngles, std::uint64_t particlesPerExposure, const SpecterDistribution<double>& specter, const BowtieFilter& bowtie, double weight = 1)
        : m_angleStep(angleStep)
        , m_sdd(SDD)
        , m_weight(weight)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_particlesPerExposure(particlesPerExposure)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
    {
    }

    std::uint64_t numberOfExposures() const
    {
        constexpr auto pi2 = PI_VAL() * 2;
        return static_cast<std::uint64_t>(pi2 / m_angleStep);
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }

    CTDIBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        const auto angle = i * m_angleStep;
        const auto organWeight = m_organFilter.useFilter() ? m_organFilter(angle) : 1.0;
        return CTDIBeamExposure<ENABLETRACKING>(angle, m_sdd, m_particlesPerExposure, m_collimationHalfAngles, &m_specter, &m_bowtieFilter, m_weight * organWeight);
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        return 1;
    }

private:
    double m_angleStep = 0;
    double m_sdd = 1;
    double m_weight = 1;
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    std::uint64_t m_particlesPerExposure = 1;
    SpecterDistribution<double> m_specter;
    BowtieFilter m_bowtieFilter;
    CTOrganAECFilter m_organFilter;
};
}
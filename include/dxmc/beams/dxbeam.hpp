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

#include "dxmc/beams/tube/tube.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {
template <bool ENABLETRACKING = false>
class DXBeamExposure {
public:
    DXBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double> specter)
        : m_pos(pos)
        , m_dirCosines(dircosines)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_NParticles(N)
        , m_weight(weight)
        , m_specter(specter)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    const std::array<double, 3>& position() const { return m_pos; }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {
        const auto angx = state.randomUniform(-m_collimationHalfAngles[0], m_collimationHalfAngles[0]);
        const auto angy = state.randomUniform(-m_collimationHalfAngles[1], m_collimationHalfAngles[1]);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = particleDirection(angx, angy),
                .energy = m_specter.sampleValue(state),
                .weight = m_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(angx, angy),
                .energy = m_specter.sampleValue(state),
                .weight = m_weight
            };
            return p;
        }
    }

protected:
    std::array<double, 3> particleDirection(double anglex, double angley) const
    {
        const auto dx = std::tan(anglex);
        const auto dy = std::tan(angley);
        const std::array pdir = {
            m_dirCosines[0][0] * dx + m_dirCosines[1][0] * dy + m_dir[0],
            m_dirCosines[0][1] * dx + m_dirCosines[1][1] * dy + m_dir[1],
            m_dirCosines[0][2] * dx + m_dirCosines[1][2] * dy + m_dir[2]
        };
        return vectormath::normalized(pdir);
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    std::uint64_t m_NParticles = 100;
    double m_weight = 1;
    SpecterDistribution<double> m_specter;
};

template <bool ENABLETRACKING = false>
class DXBeam {
public:
    DXBeam(
        const std::array<double, 3>& pos = { 0, 0, 0 },
        const std::array<std::array<double, 3>, 2>& dircosines = { { { 1, 0, 0 }, { 0, 1, 0 } } },
        const std::map<std::size_t, double>& filtrationMaterials = {})
        : m_pos(pos)
    {
        setDirectionCosines(dircosines);

        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    const std::array<double, 3>& position() const { return m_pos; }
    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }

    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dirCosines[0] = xdir;
        m_dirCosines[1] = ydir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    double DAPvalue() const { return m_measuredDAP; }
    void setDAPvalue(double dap) { m_measuredDAP = std::abs(dap); }

    const Tube& tube() const { return m_tube; }
    void setTube(const Tube& tube)
    {
        m_tube = tube;
        tubeChanged();
    }
    void setTubeVoltage(double voltage)
    {
        m_tube.setVoltage(voltage);
        tubeChanged();
    }
    void setTubeAnodeAngle(double ang)
    {
        m_tube.setAnodeAngle(ang);
        tubeChanged();
    }
    void setTubeAnodeAngleDeg(double ang)
    {
        m_tube.setAnodeAngleDeg(ang);
        tubeChanged();
    }
    void addTubeFiltrationMaterial(std::size_t Z, double mm)
    {
        m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }
    double tubeFiltration(std::size_t Z) const
    {
        return m_tube.filtration(Z);
    }
    void clearTubeFiltrationMaterials()
    {
        m_tube.clearFiltrationMaterials();
        tubeChanged();
    }
    void setTubeEnergyResolution(double energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    double tubeAlHalfValueLayer()
    {
        return m_tube.mmAlHalfValueLayer();
    }
    double tubeMeanSpecterEnergy()
    {
        return m_tube.meanSpecterEnergy();
    }
    const std::array<double, 2>& collimationAngles() const
    {
        return m_collimationHalfAngles;
    }
    void setCollimationHalfAngles(const std::array<double, 2>& angles)
    {
        setCollimationHalfAngles(angles[0], angles[1]);
    }
    void setCollimationHalfAngles(double X, double Y)
    {
        m_collimationHalfAngles[0] = std::clamp(X, 0.0, PI_VAL());
        m_collimationHalfAngles[1] = std::clamp(Y, 0.0, PI_VAL());
    }
    std::array<double, 2> collimationHalfAnglesDeg() const
    {
        auto d = m_collimationHalfAngles;
        d[0] *= RAD_TO_DEG();
        d[1] *= RAD_TO_DEG();
        return d;
    }
    void setCollimationHalfAnglesDeg(const std::array<double, 2>& angles)
    {
        setCollimationHalfAnglesDeg(angles[0], angles[1]);
    }
    void setCollimationHalfAnglesDeg(double X, double Y)
    {
        setCollimationHalfAngles(DEG_TO_RAD() * X, DEG_TO_RAD() * Y);
    }

    void setBeamSize(double beamSizeX, double beamSizeY, double sourceDetectorDistance)
    {
        if (sourceDetectorDistance > 0) {
            setCollimationHalfAngles(
                m_collimationHalfAngles[0] = std::atan(std::abs(beamSizeX) / (2 * sourceDetectorDistance)),
                m_collimationHalfAngles[1] = std::atan(std::abs(beamSizeY) / (2 * sourceDetectorDistance)));
        }
    }

    DXBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        DXBeamExposure<ENABLETRACKING> exp(m_pos, m_dirCosines, m_particlesPerExposure, m_weight, m_collimationHalfAngles, m_specter);
        return exp;
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, true);
        auto air_cand = Material<5>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();

        const auto kerma_per_history = std::transform_reduce(std::execution::par_unseq, energies.cbegin(), energies.cend(), weights.cbegin(), 0.0, std::plus<>(), [&](const auto e, const auto w) -> double {
            const auto uen = air.massEnergyTransferAttenuation(e);
            return w * e * uen;
        });

        const auto kerma_total = kerma_per_history * numberOfParticles();
        return m_measuredDAP / kerma_total;
    }

protected:
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
    double m_weight = 1;
    double m_measuredDAP = 1;
    Tube m_tube;
    SpecterDistribution<double> m_specter;
};
}
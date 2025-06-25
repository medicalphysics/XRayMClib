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
#include "dxmc/beams/utilities/spheresamplingrectangularfield.hpp"
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
class CBCTBeamExposure {
public:
    CBCTBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double>* specter)
        : m_directionSampler(collimationHalfAngles)
        , m_pos(pos)
        , m_dirCosines(dircosines)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_NParticles(N)
        , m_weight(weight)
        , m_specter(specter)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    CBCTBeamExposure() = delete;

    const std::array<double, 3>& position() const { return m_pos; }

    const std::array<std::array<double, 3>, 2>& directionCosines() const { return m_dirCosines; }

    const std::array<double, 2>& collimationHalfAngles() const
    {
        return m_collimationHalfAngles;
    }

    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

    auto sampleParticle(RandomState& state) const noexcept
    {
        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = particleDirection(state),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight
            };
            return p;
        }
    }

protected:
    std::array<double, 3> particleDirection(RandomState& state) const
    {
        const auto dir = m_directionSampler(state);
        return vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dir);
    }

private:
    SphereSamplingRectangularField m_directionSampler;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    std::uint64_t m_NParticles = 100;
    double m_weight = 1;
    const SpecterDistribution<double>* m_specter = nullptr;
};

template <bool ENABLETRACKING = false>
class CBCTBeam {
public:
    CBCTBeam(
        const std::array<double, 3>& isocenter = { 0, 0, 0 },
        const std::array<double, 3>& direction = { 0, 0, 1 },
        const std::map<std::size_t, double>& filtrationMaterials = {})
        : m_isocenter(isocenter)
    {
        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }

    std::uint64_t numberOfExposures() const
    {
        const auto dAngle = std::abs(m_angleStop - m_angleStart);
        const auto steps = std::max(1.0, dAngle / m_angleStep);
        return static_cast<std::uint64_t>(steps);
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = std::max(n, std::uint64_t { 1 }); }

    const std::array<double, 3>& isocenter() const { return m_isocenter; }
    void setIsocenter(const std::array<double, 3>& pos) { m_isocenter = pos; }

    const std::array<double, 3>& rotationAxis() const
    {
        return m_direction;
    }

    void setRotationAxis(const std::array<double, 3>& dir)
    {
        m_direction = vectormath::normalized(dir);
    }

    double startAngle() const { return m_angleStart; }
    void setStartAngle(double angle) { m_angleStart = angle; }
    double startAngleDeg() const { return m_angleStart * RAD_TO_DEG(); }
    void setStartAngleDeg(double angle) { m_angleStart = angle * DEG_TO_RAD(); }

    double stepAngle() const { return m_angleStep; }
    void setStepAngle(double angle)
    {
        m_angleStep = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }
    double stepAngleDeg() const { return m_angleStep * RAD_TO_DEG(); }
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    double stopAngle() const { return m_angleStop; }
    void setStopAngle(double angle) { m_angleStop = angle; }
    double stopAngleDeg() const { return m_angleStop * RAD_TO_DEG(); }
    void setStopAngleDeg(double angle) { m_angleStop = angle * DEG_TO_RAD(); }

    double DAPvalue() const { return m_measuredDAP; }
    void setDAPvalue(double dap) { m_measuredDAP = std::abs(dap); }

    const Tube& tube() const { return m_tube; }
    void setTube(const Tube&& tube)
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

    const std::array<double, 2>& collimationHalfAngles() const
    {
        return m_collimationHalfAngles;
    }
    void setCollimationHalfAngles(const std::array<double, 2>& angles)
    {
        setCollimationHalfAngles(angles[0], angles[1]);
    }
    void setCollimationHalfAngles(double X, double Y)
    {
        m_collimationHalfAngles[0] = std::clamp(X, 0.0, PI_VAL() * 0.5);
        m_collimationHalfAngles[1] = std::clamp(Y, 0.0, PI_VAL() * 0.5);
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

    double sourceDetectorDistance() const
    {
        return m_SDD;
    }
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    CBCTBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        auto angle = i * m_angleStep;
        if (m_angleStart > m_angleStop)
            angle = -angle;

        // finding normal vector to direction
        const auto normal_ind = vectormath::argmin3(m_direction);
        std::array<double, 3> normal = { 0, 0, 0 };
        normal[normal_ind] = 1;
        normal = vectormath::normalized(vectormath::cross(normal, m_direction));
        normal = vectormath::rotate(normal, m_direction, m_angleStart + angle);

        const auto beamdir = vectormath::cross(normal, m_direction);

        const std::array<std::array<double, 3>, 2> cosines = { normal, m_direction };

        auto pos = vectormath::add(m_isocenter, vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis

        CBCTBeamExposure<ENABLETRACKING> exp(pos, cosines, m_particlesPerExposure, m_weight, m_collimationHalfAngles, &m_specter);
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
    std::array<double, 3> m_isocenter = { 0, 0, 0 };
    std::array<double, 3> m_direction = { 0, 0, 1 };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    double m_angleStart = 0;
    double m_angleStop = PI_VAL();
    double m_angleStep = PI_VAL() / 180;
    double m_SDD = 50;
    std::uint64_t m_particlesPerExposure = 100000;
    double m_weight = 1;
    double m_measuredDAP = 1;
    Tube m_tube;
    SpecterDistribution<double> m_specter;
};
}
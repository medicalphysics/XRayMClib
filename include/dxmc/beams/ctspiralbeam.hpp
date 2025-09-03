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

#include "dxmc/beams/ctdibeam.hpp"
#include "dxmc/beams/filters/bowtiefilter.hpp"
#include "dxmc/beams/filters/ctaecfilter.hpp"
#include "dxmc/beams/filters/ctorganaecfilter.hpp"
#include "dxmc/beams/tube/tube.hpp"
#include "dxmc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"

#include <array>
#include <cmath>
#include <map>

namespace dxmc {
template <bool ENABLETRACKING = false>
class CTSpiralBeamExposure {
public:
    CTSpiralBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
        const std::array<double, 2>& collimationHalfAngles, const SpecterDistribution<double>* specter, const BowtieFilter* bowtie)
        : m_directionSampler(collimationHalfAngles)
        , m_pos(pos)
        , m_dirCosines(dircosines)
        , m_collimationHalfAngles(collimationHalfAngles)
        , m_NParticles(N)
        , m_weight(weight)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    CTSpiralBeamExposure() = delete;

    const std::array<double, 3>& position() const { return m_pos; }

    const std::array<std::array<double, 3>, 2>& directionCosines() const { return m_dirCosines; }

    const std::array<double, 2>& collimationHalfAngles() const { return m_collimationHalfAngles; }

    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

    auto sampleParticle(RandomState& state) const noexcept
    {
        const auto dir = m_directionSampler(state);
        const auto angx = std::atan(dir[0] / dir[2]);
        const auto bowtie_weight = m_bowtieFilter->operator()(angx);

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = {
                .pos = m_pos,
                .dir = vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dir),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight * bowtie_weight
            };
            p.registerPosition();
            return p;
        } else {
            Particle p = {
                .pos = m_pos,
                .dir = vectormath::changeBasis(m_dirCosines[0], m_dirCosines[1], m_dir, dir),
                .energy = m_specter->sampleValue(state),
                .weight = m_weight * bowtie_weight
            };
            return p;
        }
    }

    double weight() const
    {
        return m_weight;
    }

protected:
private:
    SphereSamplingRectangularField m_directionSampler;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };
    std::uint64_t m_NParticles = 100;
    double m_weight = 1;
    const SpecterDistribution<double>* m_specter = nullptr;
    const BowtieFilter* m_bowtieFilter = nullptr;
};

template <bool ENABLETRACKING = false>
class CTSpiralBeam {
public:
    CTSpiralBeam(
        const std::array<double, 3>& start_pos = { 0, 0, 0 },
        const std::array<double, 3>& stop_pos = { 0, 0, 0 },
        const std::map<std::size_t, double>& filtrationMaterials = {})
        : m_start(start_pos)
        , m_stop(stop_pos)
    {
        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    std::uint64_t numberOfExposures() const
    {
        const auto direction = vectormath::subtract(m_stop, m_start);
        const auto dz = m_pitch * m_collimation;
        const auto total_rot_angle = vectormath::length(direction) * (PI_VAL() * 2) / dz;
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles;
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    const std::array<double, 3>& startPosition() const { return m_start; }
    const std::array<double, 3>& stopPosition() const { return m_stop; }
    void setStartPosition(const std::array<double, 3>& start)
    {
        m_start = start;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }
    void setStopPosition(const std::array<double, 3>& stop)
    {
        m_stop = stop;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }
    void setStartStopPosition(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        m_start = start;
        m_stop = stop;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    double collimation() const { return m_collimation; }
    void setCollimation(double coll_cm)
    {
        // collimation must be larger than 1 mm (0.1 cm)
        m_collimation = std::max(std::abs(coll_cm), 0.1);
    }

    std::array<double, 2> collimationHalfAngles() const
    {
        std::array<double, 2> r = {
            std::atan(m_FOV / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    double sourceDetectorDistance() const
    {
        return m_SDD;
    }
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    double scanFieldOfView() const { return m_FOV; }
    void setScanFieldOfView(double fov_cm)
    {
        m_FOV = std::max(std::abs(fov_cm), 1.0);
    }

    const BowtieFilter& bowtieFilter() const
    {
        return m_bowtieFilter;
    }
    void setBowtieFilter(const BowtieFilter& filter)
    {
        m_bowtieFilter = filter;
    }

    double pitch() const { return m_pitch; }
    void setPitch(double p)
    {
        m_pitch = std::max(std::abs(p), 0.1);
    }

    void setCTDIvol(double ctdi) { m_CTDIvol = ctdi; }
    double CTDIvol() const { return m_CTDIvol; }
    void setCTDIdiameter(double d) { m_CTDIdiameter = d; }
    double CTDIdiameter() const { return m_CTDIdiameter; }

    double startAngle() const { return m_startAngle; }
    void setStartAngle(double angle) { m_startAngle = angle; }
    double startAngleDeg() const { return m_startAngle * RAD_TO_DEG(); }
    void setStartAngleDeg(double angle) { m_startAngle = angle * DEG_TO_RAD(); }

    double stepAngle() const { return m_stepAngle; }
    void setStepAngle(double angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }
    double stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG(); }
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

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
        auto success = m_tube.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    void clearTubeFiltrationMaterials()
    {
        m_tube.clearFiltrationMaterials();
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
    void setTubeEnergyResolution(double energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    void setAECFilter(const CTAECFilter& filter)
    {
        m_aecFilter = filter;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }
    void setAECFilterData(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& data)
    {
        m_aecFilter.setData(start, stop, data);
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }
    const CTAECFilter& AECFilter() const
    {
        return m_aecFilter;
    }

    CTOrganAECFilter& organAECFilter()
    {
        return m_organFilter;
    }
    const CTOrganAECFilter& organAECFilter() const
    {
        return m_organFilter;
    }

    CTSpiralBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        constexpr auto pi2 = PI_VAL() * 2;
        const auto angle = i * m_stepAngle;
        const auto dz = m_pitch * m_collimation * angle / pi2;
        const auto directionZ = vectormath::subtract(m_stop, m_start);
        const auto direction = vectormath::normalized(directionZ);

        // finding normal vector to direction
        const auto normal_ind = vectormath::argmin3(direction);
        std::array<double, 3> normal = { 0, 0, 0 };
        normal[normal_ind] = 1;
        normal = vectormath::normalized(vectormath::cross(normal, direction));
        normal = vectormath::rotate(normal, direction, m_startAngle + angle);

        const auto beamdir = vectormath::cross(normal, direction);

        const std::array<std::array<double, 3>, 2> cosines = { normal, direction };

        auto pos = vectormath::add(vectormath::add(m_start, vectormath::scale(direction, dz)), vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis
        const auto angx = std::atan(m_FOV / m_SDD);
        const auto angy = std::atan(0.5 * m_collimation / m_SDD);

        std::array<double, 2> angles = { angx, angy };

        const auto organWeight = m_organFilter.useFilter() ? m_organFilter(angle) : 1.0;

        const auto weight = m_weight * m_aecFilter(pos) * organWeight;

        CTSpiralBeamExposure<ENABLETRACKING> exp(pos, cosines, m_particlesPerExposure, weight, angles, &m_specter, &m_bowtieFilter);
        return exp;
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        // generating scoring world
        using Phantom = CTDIPhantom<5, 1>;
        World<Phantom> world;
        world.reserveNumberOfItems(1);
        const auto& ctdi = world.template addItem<Phantom>({ m_CTDIdiameter });
        world.build();

        // generating CTDIbeam
        const auto angx = std::atan(m_FOV / m_SDD);
        const auto angy = std::atan(0.5 * m_collimation / m_SDD);
        const std::array<double, 2> collimationHalfAngles = { angx, angy };
        CTDIBeam beam(m_stepAngle, m_SDD, collimationHalfAngles, m_particlesPerExposure, m_specter, m_bowtieFilter, m_organFilter);

        Transport transport;

        transport(world, beam, progress);

        const auto ctdiw_calc = (ctdi.centerDoseScored() + 2 * ctdi.pheriferyDoseScored()) * 10.0 / (3 * m_collimation);

        const auto ctdiw_beam = m_CTDIvol * m_pitch;
        return ctdiw_beam / ctdiw_calc;
    }

protected:
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 };
    std::array<double, 3> m_stop = { 0, 0, 0 };
    double m_FOV = 50;
    double m_SDD = 100;
    double m_collimation = 1; // cm
    double m_pitch = 1;
    double m_startAngle = 0;
    double m_stepAngle = 0.018; // about a degree;
    double m_weight = 1;
    double m_CTDIvol = 1;
    double m_CTDIdiameter = 32;
    std::uint64_t m_particlesPerExposure = 100;
    Tube m_tube;
    SpecterDistribution<double> m_specter;
    CTAECFilter m_aecFilter;
    BowtieFilter m_bowtieFilter;
    CTOrganAECFilter m_organFilter;
};
}
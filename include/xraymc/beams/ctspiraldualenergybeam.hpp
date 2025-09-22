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

#include "xraymc/beams/ctdibeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/filters/ctaecfilter.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <cmath>
#include <map>

namespace xraymc {
template <bool ENABLETRACKING = false>
class CTSpiralDualEnergyBeam {
public:
    CTSpiralDualEnergyBeam(
        const std::array<double, 3>& start_pos = { 0, 0, 0 },
        const std::array<double, 3>& stop_pos = { 0, 0, 0 },
        const std::map<std::size_t, double>& filtrationMaterials = {})
        : m_start(start_pos)
        , m_stop(stop_pos)
    {
        for (const auto [Z, mm] : filtrationMaterials) {
            m_tubeA.addFiltrationMaterial(Z, mm);
            m_tubeB.addFiltrationMaterial(Z, mm);
        }
        tubeChanged();
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    std::uint64_t numberOfExposures() const
    {
        const auto direction = vectormath::subtract(m_stop, m_start);
        const auto dz = m_pitch * m_collimation;
        const auto total_rot_angle = vectormath::length(direction) * (PI_VAL() * 2) / dz;
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles * 2;
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

    double sourceDetectorDistance() const { return m_SDD; }
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    double scanFieldOfViewB() const { return m_FOVB; }
    void setScanFieldOfViewB(double fov_cm)
    {
        m_FOVB = std::max(std::abs(fov_cm), 1.0);
    }
    double scanFieldOfViewA() const { return m_FOVA; }
    void setScanFieldOfViewA(double fov_cm)
    {
        m_FOVA = std::max(std::abs(fov_cm), 1.0);
    }

    std::array<double, 2> collimationHalfAnglesA() const
    {
        std::array r = {
            std::atan(m_FOVA / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    std::array<double, 2> collimationHalfAnglesB() const
    {
        std::array r = {
            std::atan(m_FOVB / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    const BowtieFilter& bowtieFilterA() const
    {
        return m_bowtieFilterA;
    }
    void setBowtieFilterA(const BowtieFilter& filter)
    {
        m_bowtieFilterA = filter;
    }
    const BowtieFilter& bowtieFilterB() const
    {
        return m_bowtieFilterB;
    }
    void setBowtieFilterB(const BowtieFilter& filter)
    {
        m_bowtieFilterB = filter;
    }

    double pitch() const { return m_pitch; }
    void setPitch(double p)
    {
        m_pitch = std::max(std::abs(p), 0.1);
    }

    void setCTDIvol(double ctdi) { m_CTDIvol = ctdi; }
    double CTDIvol() const { return m_CTDIvol; }
    void setCTDIdiameter(double d) { m_CTDIdiameter = std::max(d, 3.0); }
    double CTDIdiameter() const { return m_CTDIdiameter; }

    double startAngle() const { return m_startAngle; }
    void setStartAngle(double angle) { m_startAngle = angle; }
    double startAngleDeg() const { return m_startAngle * RAD_TO_DEG(); }
    void setStartAngleDeg(double angle) { m_startAngle = angle * DEG_TO_RAD(); }

    double tubeBoffsetAngle() const { return m_tubeBoffsetAngle; }
    void setTubeBoffsetAngle(double angle) { m_tubeBoffsetAngle = angle; }
    double tubeBoffsetAngleDeg() const { return m_tubeBoffsetAngle * RAD_TO_DEG(); }
    void setTubeBoffsetAngleDeg(double angle) { m_tubeBoffsetAngle = angle * DEG_TO_RAD(); }

    double stepAngle() const { return m_stepAngle; }
    void setStepAngle(double angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }
    double stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG(); }
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    const Tube& tubeA() const { return m_tubeA; }
    void setTubeA(const Tube&& tube)
    {
        m_tubeA = tube;
        tubeChanged();
    }
    const Tube& tubeB() const { return m_tubeB; }
    void setTubeB(const Tube&& tube)
    {
        m_tubeB = tube;
        tubeChanged();
    }
    void setTubeAVoltage(double voltage)
    {
        m_tubeA.setVoltage(voltage);
        tubeChanged();
    }
    void setTubeBVoltage(double voltage)
    {
        m_tubeB.setVoltage(voltage);
        tubeChanged();
    }
    void setTubesAnodeAngle(double ang)
    {
        m_tubeA.setAnodeAngle(ang);
        m_tubeB.setAnodeAngle(ang);
        tubeChanged();
    }
    void setTubesAnodeAngleDeg(double ang)
    {
        m_tubeA.setAnodeAngleDeg(ang);
        m_tubeB.setAnodeAngleDeg(ang);
        tubeChanged();
    }
    void addTubeAFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tubeA.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    void addTubeBFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tubeB.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    void clearTubesFiltrationMaterials()
    {
        m_tubeA.clearFiltrationMaterials();
        m_tubeB.clearFiltrationMaterials();
        tubeChanged();
    }
    double tubeAAlHalfValueLayer()
    {
        return m_tubeA.mmAlHalfValueLayer();
    }
    double tubeBAlHalfValueLayer()
    {
        return m_tubeB.mmAlHalfValueLayer();
    }
    double tubeAMeanSpecterEnergy()
    {
        return m_tubeA.meanSpecterEnergy();
    }
    double tubeBMeanSpecterEnergy()
    {
        return m_tubeB.meanSpecterEnergy();
    }
    void setTubesEnergyResolution(double energyResolution)
    {
        m_tubeA.setEnergyResolution(energyResolution);
        m_tubeB.setEnergyResolution(energyResolution);
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

    // helpers to set tube weights
    void setRelativeMasTubeA(double mas)
    {
        m_relativeMasA = mas;
        tubeChanged();
    }
    void setRelativeMasTubeB(double mas)
    {
        m_relativeMasB = mas;
        tubeChanged();
    }
    double tubeRelativeWeightA() const
    {
        return m_weightA;
    }
    double tubeRelativeWeightB() const
    {
        return m_weightB;
    }
    double relativeMasTubeA() const
    {
        return m_relativeMasA;
    }
    double relativeMasTubeB() const
    {
        return m_relativeMasB;
    }

    CTOrganAECFilter& organAECFilter()
    {
        return m_organFilter;
    }
    const CTOrganAECFilter& organAECFilter() const
    {
        return m_organFilter;
    }

    CTSpiralBeamExposure<ENABLETRACKING> exposure(std::size_t dualExposureIndex) const noexcept
    {
        const std::size_t i = dualExposureIndex / 2;
        const bool isTubeB = dualExposureIndex % 2 == 1;

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
        const auto rot_angle = isTubeB ? m_startAngle + angle + m_tubeBoffsetAngle : m_startAngle + angle;
        normal = vectormath::rotate(normal, direction, rot_angle);

        const auto beamdir = vectormath::cross(normal, direction);

        const std::array<std::array<double, 3>, 2> cosines = { normal, direction };

        auto pos = vectormath::add(vectormath::add(m_start, vectormath::scale(direction, dz)), vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis
        const auto angx = isTubeB ? std::atan(m_FOVB / m_SDD) : std::atan(m_FOVA / m_SDD);
        const auto angy = std::atan(0.5 * m_collimation / m_SDD);

        std::array<double, 2> angles = { angx, angy };

        auto weight = isTubeB ? m_weightB : m_weightA;
        weight *= m_aecFilter(pos);
        if (m_organFilter.useFilter())
            weight *= m_organFilter(angle);
        const SpecterDistribution<double>* const specter = isTubeB ? &m_specterB : &m_specterA;
        const BowtieFilter* const bowtie = isTubeB ? &m_bowtieFilterB : &m_bowtieFilterA;
        CTSpiralBeamExposure<ENABLETRACKING> exp(pos, cosines, m_particlesPerExposure, weight, angles, specter, bowtie);
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
        const auto angxA = std::atan(m_FOVA / m_SDD);
        const auto angy = std::atan(0.5 * m_collimation / m_SDD);
        const std::array<double, 2> collimationAnglesA = { angxA, angy };
        CTDIBeam beamA(m_stepAngle, m_SDD, collimationAnglesA, m_particlesPerExposure, m_specterA, m_bowtieFilterA, m_organFilter, m_weightA);
        Transport transport;
        transport(world, beamA, progress, false);

        const auto angxB = std::atan(m_FOVB / m_SDD);
        const std::array collimationAnglesB = { angxB, angy };
        CTDIBeam beamB(m_stepAngle, m_SDD, collimationAnglesB, m_particlesPerExposure, m_specterB, m_bowtieFilterB, m_organFilter, m_weightB);
        transport(world, beamB, progress, false);

        world.addEnergyScoredToDoseScore();

        const auto ctdiw_calc = (ctdi.centerDoseScored() + 2 * ctdi.pheriferyDoseScored()) * 10.0 / (3 * m_collimation);

        const auto ctdiw_beam = m_CTDIvol * m_pitch;
        return ctdiw_beam / ctdiw_calc;
    }

protected:
    void tubeChanged()
    {
        auto energiesA = m_tubeA.getEnergy();
        auto weightsA = m_tubeA.getSpecter(energiesA, false);
        m_specterA = SpecterDistribution(energiesA, weightsA);
        const auto weightA = m_relativeMasA * std::reduce(std::execution::par_unseq, weightsA.cbegin(), weightsA.cend(), 0.0);

        auto energiesB = m_tubeB.getEnergy();
        auto weightsB = m_tubeB.getSpecter(energiesB, false);
        m_specterB = SpecterDistribution(energiesB, weightsB);
        const auto weightB = m_relativeMasB * std::reduce(std::execution::par_unseq, weightsB.cbegin(), weightsB.cend(), 0.0);

        m_weightA = 2 * weightA / (weightA + weightB);
        m_weightB = 2 * weightB / (weightA + weightB);
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 };
    std::array<double, 3> m_stop = { 0, 0, 0 };
    double m_FOVA = 50;
    double m_FOVB = 50;
    double m_SDD = 100;
    double m_collimation = 1; // cm
    double m_pitch = 1;
    double m_startAngle = 0;
    double m_tubeBoffsetAngle = DEG_TO_RAD() * 90;
    double m_stepAngle = DEG_TO_RAD();
    double m_weightA = 1;
    double m_weightB = 1;
    double m_relativeMasA = 1;
    double m_relativeMasB = 1;
    double m_CTDIvol = 1;
    double m_CTDIdiameter = 32;
    std::uint64_t m_particlesPerExposure = 100;
    Tube m_tubeA;
    Tube m_tubeB;
    SpecterDistribution<double> m_specterA;
    SpecterDistribution<double> m_specterB;
    CTAECFilter m_aecFilter;
    BowtieFilter m_bowtieFilterA;
    BowtieFilter m_bowtieFilterB;
    CTOrganAECFilter m_organFilter;
};
}

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
#include "xraymc/beams/filters/bowtiefilter.hpp"
#include "xraymc/beams/filters/ctaecfilter.hpp"
#include "xraymc/beams/filters/ctorganaecfilter.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <cmath>
#include <map>
#include <span>

namespace xraymc {

/**
 * @brief A single gantry exposure used internally by `CTSpiralBeam`.
 *
 * Represents the X-ray tube at one position along the helical trajectory. Photon
 * directions are sampled from a `SphereSamplingRectangularField` over the fan/slice
 * collimation angles, then modulated by a `BowtieFilter` weight computed from the
 * in-plane fan angle. Energies are drawn from a `SpecterDistribution`. The combined
 * photon weight is `m_weight × bowtie_weight`, where `m_weight` already incorporates
 * the AEC and organ-AEC modulation applied by `CTSpiralBeam::exposure()`.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class CTSpiralBeamExposure {
public:
    /**
     * @brief Constructs a single CT spiral exposure.
     *
     * Derives the central beam direction as `dircosines[0] × dircosines[1]` and builds a
     * `SphereSamplingRectangularField` for the supplied collimation half-angles.
     *
     * @param pos                   Source position in world space [cm].
     * @param dircosines            Two orthonormal vectors `{cos_x, cos_y}` — the in-plane
     *                              normal and the helix axis direction.
     * @param N                     Number of photon histories for this exposure.
     * @param weight                Combined statistical weight (base × AEC × organ-AEC).
     * @param collimationHalfAngles `{half_fan, half_slice}` collimation half-angles [rad],
     *                              where half_fan = atan(FOV/SDD) and half_slice = atan(collimation/(2·SDD)).
     * @param specter               Non-owning pointer to the energy spectrum sampler.
     * @param bowtie                Non-owning pointer to the bowtie filter.
     */
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

    /// @brief Default constructor is deleted — an exposure must always have geometry and a spectrum.
    CTSpiralBeamExposure() = delete;

    /// @brief Returns the source position for this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the two direction cosines `{cos_x, cos_y}` for this exposure.
    const std::array<std::array<double, 3>, 2>& directionCosines() const { return m_dirCosines; }

    /// @brief Returns the collimation half-angles `{half_fan, half_slice}` [rad].
    const std::array<double, 2>& collimationHalfAngles() const { return m_collimationHalfAngles; }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * 1. Draws a fan-angle direction from `SphereSamplingRectangularField`.
     * 2. Derives the in-plane fan angle and looks up the bowtie weight.
     * 3. Transforms the direction to world space via `changeBasis`.
     * 4. Draws an energy from the spectrum.
     * 5. Returns a photon with weight = `m_weight × bowtie_weight`.
     *
     * @param state  Per-thread PRNG state.
     * @return `Particle` or `ParticleTrack` depending on `ENABLETRACKING`.
     */
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

    /// @brief Returns the combined photon weight (base × AEC × organ-AEC) for this exposure.
    double weight() const
    {
        return m_weight;
    }

private:
    SphereSamplingRectangularField m_directionSampler; ///< Samples directions within the fan/slice collimation field.
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Central beam direction (derived from direction cosines).
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< Orthonormal basis {in-plane normal, helix axis} perpendicular to m_dir.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 }; ///< Collimation half-angles {half_fan, half_slice} [rad].
    std::uint64_t m_NParticles = 100; ///< Number of photon histories.
    double m_weight = 1; ///< Combined photon weight (base × AEC × organ-AEC).
    const SpecterDistribution<double>* m_specter = nullptr; ///< Non-owning pointer to the energy spectrum.
    const BowtieFilter* m_bowtieFilter = nullptr; ///< Non-owning pointer to the bowtie filter.
};

/**
 * @brief A helical CT spiral beam, satisfying `BeamType`.
 *
 * Simulates a complete helical (spiral) CT scan between two positions. The X-ray tube
 * rotates around the patient axis while the table advances, tracing a helix from
 * `startPosition()` to `stopPosition()`. The source-to-isocentre distance is SDD/2.
 *
 * The number of exposures is derived from the helix geometry:
 *   N = floor(|stop − start| × 2π / (pitch × collimation) / stepAngle)
 *
 * Per-exposure photon weights incorporate:
 * - A base weight `m_weight`.
 * - An axial AEC weight from `CTAECFilter` (tube current modulation along the patient axis).
 * - An angular organ-AEC weight from `CTOrganAECFilter` (angular tube current modulation).
 * - An in-plane bowtie filter weight from `BowtieFilter`.
 *
 * `calibrationFactor()` performs an internal full-rotation CTDI simulation using a
 * `CTDIPhantom` of the configured diameter and computes the factor needed to scale
 * scored doses to match the prescribed CTDIvol.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CTSpiralBeam {
public:
    /**
     * @brief Constructs a CT spiral beam scanning from @p start_pos to @p stop_pos.
     *
     * Applies any supplied filtration materials to the internal tube, rebuilds the
     * spectrum cache, and normalises the AEC filter between the two endpoints.
     *
     * @param start_pos           Scan start position [cm]. Default: origin.
     * @param stop_pos            Scan stop position [cm]. Default: origin.
     * @param filtrationMaterials Map of `{atomic number Z → thickness [mm]}` added as
     *                            filtration to the internal `Tube`. Default: empty.
     */
    CTSpiralBeam(
        const std::array<double, 3>& start_pos = { 0, 0, 0 },
        const std::array<double, 3>& stop_pos = { 0, 0, 0 },
        const std::map<std::size_t, double>& filtrationMaterials = { })
        : m_start(start_pos)
        , m_stop(stop_pos)
    {
        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /**
     * @brief Returns the number of gantry exposures along the full helix.
     *
     * Computed as floor(|stop − start| × 2π / (pitch × collimation) / stepAngle).
     */
    std::uint64_t numberOfExposures() const
    {
        const auto direction = vectormath::subtract(m_stop, m_start);
        const auto dz = m_pitch * m_collimation;
        const auto total_rot_angle = vectormath::length(direction) * (PI_VAL() * 2) / dz;
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles;
    }

    /// @brief Returns the total number of photon histories across all exposures.
    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }

    /// @brief Returns the number of photon histories per exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    /**
     * @brief Sets the number of photon histories per exposure.
     * @param n  Histories per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /// @brief Returns the scan start position [cm].
    const std::array<double, 3>& startPosition() const { return m_start; }

    /// @brief Returns the scan stop position [cm].
    const std::array<double, 3>& stopPosition() const { return m_stop; }

    /**
     * @brief Sets the scan start position [cm]; renormalises the AEC filter.
     * @param start  New start position [cm].
     */
    void setStartPosition(const std::array<double, 3>& start)
    {
        m_start = start;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /**
     * @brief Sets the scan stop position [cm]; renormalises the AEC filter.
     * @param stop  New stop position [cm].
     */
    void setStopPosition(const std::array<double, 3>& stop)
    {
        m_stop = stop;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /**
     * @brief Sets both scan endpoints simultaneously [cm]; renormalises the AEC filter once.
     * @param start  New start position [cm].
     * @param stop   New stop position [cm].
     */
    void setStartStopPosition(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        m_start = start;
        m_stop = stop;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /// @brief Returns the total beam collimation width [cm].
    double collimation() const { return m_collimation; }

    /**
     * @brief Sets the total beam collimation width; clamped to ≥ 1 mm (0.1 cm).
     * @param coll_cm  Collimation width [cm].
     */
    void setCollimation(double coll_cm)
    {
        // collimation must be larger than 1 mm (0.1 cm)
        m_collimation = std::max(std::abs(coll_cm), 0.1);
    }

    /**
     * @brief Returns the current collimation half-angles as `{half_fan, half_slice}` [rad].
     *
     * Computed from the scan FOV and collimation:
     *   half_fan   = atan(FOV / SDD)
     *   half_slice = atan(collimation / (2 × SDD))
     */
    std::array<double, 2> collimationHalfAngles() const
    {
        std::array<double, 2> r = {
            std::atan(m_FOV / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    /// @brief Returns the source-to-detector distance [cm].
    double sourceDetectorDistance() const
    {
        return m_SDD;
    }

    /**
     * @brief Sets the source-to-detector distance [cm]; clamped to ≥ 1 cm.
     * @param SDD_cm  Source-to-detector distance [cm].
     */
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    /// @brief Returns the scan field of view radius [cm].
    double scanFieldOfView() const { return m_FOV; }

    /**
     * @brief Sets the scan field of view radius [cm]; clamped to ≥ 1 cm.
     * @param fov_cm  Scan FOV radius [cm].
     */
    void setScanFieldOfView(double fov_cm)
    {
        m_FOV = std::max(std::abs(fov_cm), 1.0);
    }

    /// @brief Returns a const reference to the bowtie filter.
    const BowtieFilter& bowtieFilter() const
    {
        return m_bowtieFilter;
    }

    /**
     * @brief Replaces the bowtie filter (copied).
     * @param filter  New `BowtieFilter` instance.
     */
    void setBowtieFilter(const BowtieFilter& filter)
    {
        m_bowtieFilter = filter;
    }

    /// @brief Returns the helical pitch (table advance per rotation / collimation width).
    double pitch() const { return m_pitch; }

    /**
     * @brief Sets the helical pitch; clamped to ≥ 0.1.
     * @param p  Pitch value (dimensionless).
     */
    void setPitch(double p)
    {
        m_pitch = std::max(std::abs(p), 0.1);
    }

    /**
     * @brief Sets the prescribed CTDIvol used for dose calibration [mGy].
     * @param ctdi  CTDIvol [mGy].
     */
    void setCTDIvol(double ctdi) { m_CTDIvol = ctdi; }

    /// @brief Returns the prescribed CTDIvol used for dose calibration [mGy].
    double CTDIvol() const { return m_CTDIvol; }

    /**
     * @brief Sets the CTDI phantom diameter for the internal calibration simulation [cm].
     * @param d  Phantom diameter [cm]. Typical values: 16 cm (head) or 32 cm (body).
     */
    void setCTDIdiameter(double d) { m_CTDIdiameter = d; }

    /// @brief Returns the CTDI phantom diameter used for calibration [cm].
    double CTDIdiameter() const { return m_CTDIdiameter; }

    /// @brief Returns the helix start angle [rad].
    double startAngle() const { return m_startAngle; }
    /// @brief Sets the helix start angle [rad].
    void setStartAngle(double angle) { m_startAngle = angle; }
    /// @brief Returns the helix start angle [deg].
    double startAngleDeg() const { return m_startAngle * RAD_TO_DEG(); }
    /// @brief Sets the helix start angle [deg].
    void setStartAngleDeg(double angle) { m_startAngle = angle * DEG_TO_RAD(); }

    /// @brief Returns the angular step between consecutive exposures [rad].
    double stepAngle() const { return m_stepAngle; }

    /**
     * @brief Sets the angular step between consecutive exposures [rad].
     *
     * Clamped to a minimum of 0.1° to avoid degenerate helices.
     *
     * @param angle  Desired step size [rad]; the absolute value is used.
     */
    void setStepAngle(double angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }

    /// @brief Returns the angular step between consecutive exposures [deg].
    double stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG(); }
    /// @brief Sets the angular step between consecutive exposures [deg].
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    /// @brief Returns a const reference to the internal X-ray tube model.
    const Tube& tube() const { return m_tube; }

    /**
     * @brief Replaces the internal tube model; rebuilds the spectrum cache.
     * @param tube  New `Tube` instance (moved in).
     */
    void setTube(const Tube&& tube)
    {
        m_tube = tube;
        tubeChanged();
    }

    /**
     * @brief Sets the tube voltage [kV]; rebuilds the spectrum cache.
     * @param voltage  Tube voltage [kV].
     */
    void setTubeVoltage(double voltage)
    {
        m_tube.setVoltage(voltage);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle [rad]; rebuilds the spectrum cache.
     * @param ang  Anode angle [rad].
     */
    void setTubeAnodeAngle(double ang)
    {
        m_tube.setAnodeAngle(ang);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle [deg]; rebuilds the spectrum cache.
     * @param ang  Anode angle [deg].
     */
    void setTubeAnodeAngleDeg(double ang)
    {
        m_tube.setAnodeAngleDeg(ang);
        tubeChanged();
    }

    /**
     * @brief Adds a filtration material; rebuilds the spectrum cache on success.
     * @param Z   Atomic number of the filter material.
     * @param mm  Filter thickness [mm].
     */
    void addTubeFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tube.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }

    /// @brief Removes all filtration materials; rebuilds the spectrum cache.
    void clearTubeFiltrationMaterials()
    {
        m_tube.clearFiltrationMaterials();
        tubeChanged();
    }

    /// @brief Returns the aluminium half-value layer of the tube spectrum [mm Al].
    double tubeAlHalfValueLayer()
    {
        return m_tube.mmAlHalfValueLayer();
    }

    /// @brief Returns the mean photon energy of the tube spectrum [keV].
    double tubeMeanSpecterEnergy()
    {
        return m_tube.meanSpecterEnergy();
    }

    /**
     * @brief Sets the energy resolution of the tube spectrum; rebuilds the spectrum cache.
     * @param energyResolution  Energy bin width [keV].
     */
    void setTubeEnergyResolution(double energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    /**
     * @brief Replaces the axial AEC filter; renormalises it between the current endpoints.
     * @param filter  New `CTAECFilter` instance.
     */
    void setAECFilter(const CTAECFilter& filter)
    {
        m_aecFilter = filter;
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /**
     * @brief Sets the axial AEC filter data from a raw weight profile; renormalises it.
     *
     * @param start  Start position of the weight profile [cm].
     * @param stop   Stop position of the weight profile [cm].
     * @param data   Per-position tube-current weights (normalised internally).
     */
    void setAECFilterData(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& data)
    {
        m_aecFilter.setData(start, stop, data);
        m_aecFilter.normalizeBetween(m_start, m_stop);
    }

    /// @brief Returns a const reference to the axial AEC filter.
    const CTAECFilter& AECFilter() const
    {
        return m_aecFilter;
    }

    /// @brief Returns a mutable reference to the organ AEC filter.
    CTOrganAECFilter& organAECFilter()
    {
        return m_organFilter;
    }

    /// @brief Returns a const reference to the organ AEC filter.
    const CTOrganAECFilter& organAECFilter() const
    {
        return m_organFilter;
    }

    /**
     * @brief Returns the `CTSpiralBeamExposure` for helix position @p i.
     *
     * Computes the gantry angle as `i × stepAngle` and the table advance as
     * `pitch × collimation × angle / (2π)`. The in-plane normal is derived by
     * rotating a seed vector around the helix axis by `startAngle + angle`. The
     * source is placed at `start + direction × dz − beamdir × SDD/2`. Collimation
     * half-angles are computed from the scan FOV and collimation width. The total
     * photon weight is `m_weight × aecFilter(pos) × organWeight`.
     *
     * @param i  Zero-based exposure index; should be in [0, numberOfExposures()).
     * @return A fully configured `CTSpiralBeamExposure` for the requested position.
     */
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

    /**
     * @brief Returns the dose calibration factor for this spiral beam.
     *
     * Runs an internal CTDI simulation: constructs a `CTDIPhantom` of diameter
     * `CTDIdiameter()`, performs a full 360° `CTDIBeam` transport, then computes:
     *
     *   factor = (CTDIvol × pitch) / CTDIw_calc
     *
     * where CTDIw_calc = (D_centre + 2 × D_periphery) × 10 / (3 × collimation).
     * The organ-AEC filter is forwarded to the `CTDIBeam` so angular modulation
     * is correctly accounted for.
     *
     * @param progress  Optional progress reporter forwarded to the internal transport.
     * @return Calibration factor to multiply accumulated energy scores by.
     */
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

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMCTSpiralBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMCTSpiralBeam";
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether @p data begins with the expected magic identifier.
     * @param data  Byte span to inspect; must be at least 32 bytes.
     * @return True if the first 32 bytes match `magicID()`, false otherwise.
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the beam configuration to a byte buffer.
     *
     * Writes start/stop positions, FOV, SDD, collimation, start angle, step angle,
     * weight, CTDIvol, CTDI phantom diameter, pitch, particles per exposure, and the
     * serialized `Tube`, `BowtieFilter`, `CTOrganAECFilter`, and `CTAECFilter`
     * sub-blocks using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        Serializer::serialize(m_start, buffer);
        Serializer::serialize(m_stop, buffer);
        Serializer::serialize(m_FOV, buffer);
        Serializer::serialize(m_SDD, buffer);
        Serializer::serialize(m_collimation, buffer);
        Serializer::serialize(m_startAngle, buffer);
        Serializer::serialize(m_stepAngle, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_CTDIvol, buffer);
        Serializer::serialize(m_CTDIdiameter, buffer);
        Serializer::serialize(m_pitch, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serializeItem(m_tube, buffer);
        Serializer::serializeItem(m_bowtieFilter, buffer);
        Serializer::serializeItem(m_organFilter, buffer);
        Serializer::serializeItem(m_aecFilter, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `CTSpiralBeam` from a serialized byte buffer.
     *
     * Reads all scalar fields then deserializes the `Tube`, `BowtieFilter`,
     * `CTOrganAECFilter`, and `CTAECFilter` sub-blocks in order. Returns
     * `std::nullopt` if any sub-block fails to deserialize.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CTSpiralBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<CTSpiralBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        CTSpiralBeam<ENABLETRACKING> item;

        buffer = Serializer::deserialize(item.m_start, buffer);
        buffer = Serializer::deserialize(item.m_stop, buffer);
        buffer = Serializer::deserialize(item.m_FOV, buffer);
        buffer = Serializer::deserialize(item.m_SDD, buffer);
        buffer = Serializer::deserialize(item.m_collimation, buffer);
        buffer = Serializer::deserialize(item.m_startAngle, buffer);
        buffer = Serializer::deserialize(item.m_stepAngle, buffer);
        buffer = Serializer::deserialize(item.m_weight, buffer);
        buffer = Serializer::deserialize(item.m_CTDIvol, buffer);
        buffer = Serializer::deserialize(item.m_CTDIdiameter, buffer);
        buffer = Serializer::deserialize(item.m_pitch, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        auto name = Serializer::getNameIDTemplate();
        auto item_buffer = Serializer::getEmptyBuffer();

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto tube_opt = Tube::deserialize(item_buffer);
        if (!tube_opt) {
            return std::nullopt;
        } else {
            item.m_tube = tube_opt.value();
            item.tubeChanged();
        }

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto bowtie_opt = BowtieFilter::deserialize(item_buffer);
        if (!bowtie_opt) {
            return std::nullopt;
        } else {
            item.m_bowtieFilter = bowtie_opt.value();
        }

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto aec_organ_opt = CTOrganAECFilter::deserialize(item_buffer);
        if (!aec_organ_opt) {
            return std::nullopt;
        } else {
            item.m_organFilter = aec_organ_opt.value();
        }

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto aec_opt = CTAECFilter::deserialize(item_buffer);
        if (!aec_opt) {
            return std::nullopt;
        } else {
            item.m_aecFilter = aec_opt.value();
        }

        return std::make_optional(item);
    }

protected:
    /**
     * @brief Rebuilds the `SpecterDistribution` cache from the current `Tube` state.
     *
     * Called automatically after any tube parameter change. Queries the tube for its
     * energy bins and unnormalised weights, then constructs a new `SpecterDistribution`
     * for use in `exposure()` and `calibrationFactor()`.
     */
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 }; ///< Scan start position [cm].
    std::array<double, 3> m_stop = { 0, 0, 0 }; ///< Scan stop position [cm].
    double m_FOV = 50; ///< Scan field-of-view radius [cm].
    double m_SDD = 100; ///< Source-to-detector distance [cm].
    double m_collimation = 1; ///< Total beam collimation width [cm].
    double m_pitch = 1; ///< Helical pitch (table advance / collimation per rotation).
    double m_startAngle = 0; ///< Helix start angle [rad].
    double m_stepAngle = 0.018; ///< Angular step between exposures [rad]. Default: ~1°.
    double m_weight = 1; ///< Base photon weight (before AEC and organ-AEC modulation).
    double m_CTDIvol = 1; ///< Prescribed CTDIvol for dose calibration [mGy].
    double m_CTDIdiameter = 32; ///< CTDI phantom diameter for calibration simulation [cm].
    std::uint64_t m_particlesPerExposure = 100; ///< Photon histories per gantry position.
    Tube m_tube; ///< X-ray tube model (owned copy).
    SpecterDistribution<double> m_specter; ///< Spectrum cache rebuilt on tube changes.
    CTAECFilter m_aecFilter; ///< Axial AEC filter (tube current modulation along patient axis).
    BowtieFilter m_bowtieFilter; ///< Bowtie filter (angle-dependent intensity modulation).
    CTOrganAECFilter m_organFilter; ///< Organ AEC filter (angular tube current modulation).
};
}
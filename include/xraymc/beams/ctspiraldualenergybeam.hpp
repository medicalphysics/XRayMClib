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
 * @brief A dual-energy helical CT spiral beam, satisfying `BeamType`.
 *
 * Simulates a dual-energy CT scan using two X-ray tubes (A and B) interleaved along a
 * single helical trajectory. Each physical gantry angle produces two exposures — one from
 * tube A and one from tube B — so `numberOfExposures()` is twice the number of helix steps.
 * Even-indexed exposures (0, 2, 4, …) belong to tube A; odd-indexed exposures (1, 3, 5, …)
 * belong to tube B, which is rotated by `m_tubeBoffsetAngle` (default: 90°) relative to A.
 *
 * Relative tube weights are computed from user-supplied mAs values and scan FOV radii so
 * that the photon density per unit area is equalised between tubes:
 *   weight_A = 2 × (mAs_A × FOV_A × spectrumSum_A) / (weight_A + weight_B)
 *   weight_B = 2 × (mAs_B × FOV_B × spectrumSum_B) / (weight_A + weight_B)
 *
 * `calibrationFactor()` performs an internal dual-rotation CTDI simulation using both
 * tubes and computes the factor needed to scale scored doses to the prescribed CTDIvol.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CTSpiralDualEnergyBeam {
public:
    /**
     * @brief Constructs a dual-energy CT spiral beam scanning from @p start_pos to @p stop_pos.
     *
     * Applies @p filtrationMaterials to both tubes, rebuilds the spectrum caches, and
     * normalises the AEC filter between the two endpoints.
     *
     * @param start_pos           Scan start position [cm]. Default: origin.
     * @param stop_pos            Scan stop position [cm]. Default: origin.
     * @param filtrationMaterials Map of `{atomic number Z → thickness [mm]}` applied to
     *                            both tubes. Default: empty.
     */
    CTSpiralDualEnergyBeam(
        const std::array<double, 3>& start_pos = { 0, 0, 0 },
        const std::array<double, 3>& stop_pos = { 0, 0, 0 },
        const std::map<std::size_t, double>& filtrationMaterials = { })
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

    /**
     * @brief Returns the total number of exposures across both tubes.
     *
     * Computed as 2 × floor(|stop − start| × 2π / (pitch × collimation) / stepAngle).
     * Each helix step yields one tube-A exposure followed by one tube-B exposure.
     */
    std::uint64_t numberOfExposures() const
    {
        const auto direction = vectormath::subtract(m_stop, m_start);
        const auto dz = m_pitch * m_collimation;
        const auto total_rot_angle = vectormath::length(direction) * (PI_VAL() * 2) / dz;
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles * 2;
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

    /// @brief Returns the source-to-detector distance [cm].
    double sourceDetectorDistance() const { return m_SDD; }

    /**
     * @brief Sets the source-to-detector distance [cm]; clamped to ≥ 1 cm.
     * @param SDD_cm  Source-to-detector distance [cm].
     */
    void setSourceDetectorDistance(double SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), 1.0);
    }

    /// @brief Returns the scan field-of-view radius for tube B [cm].
    double scanFieldOfViewB() const { return m_FOVB; }

    /**
     * @brief Sets the scan FOV radius for tube B [cm]; clamped to ≥ 1 cm. Rebuilds spectrum weights.
     * @param fov_cm  Scan FOV radius [cm].
     */
    void setScanFieldOfViewB(double fov_cm)
    {
        m_FOVB = std::max(std::abs(fov_cm), 1.0);
        tubeChanged();
    }

    /// @brief Returns the scan field-of-view radius for tube A [cm].
    double scanFieldOfViewA() const { return m_FOVA; }

    /**
     * @brief Sets the scan FOV radius for tube A [cm]; clamped to ≥ 1 cm. Rebuilds spectrum weights.
     * @param fov_cm  Scan FOV radius [cm].
     */
    void setScanFieldOfViewA(double fov_cm)
    {
        m_FOVA = std::max(std::abs(fov_cm), 1.0);
        tubeChanged();
    }

    /**
     * @brief Returns the collimation half-angles for tube A as `{half_fan_A, half_slice}` [rad].
     *
     * half_fan_A = atan(FOV_A / SDD), half_slice = atan(collimation / (2 × SDD)).
     */
    std::array<double, 2> collimationHalfAnglesA() const
    {
        std::array r = {
            std::atan(m_FOVA / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    /**
     * @brief Returns the collimation half-angles for tube B as `{half_fan_B, half_slice}` [rad].
     *
     * half_fan_B = atan(FOV_B / SDD), half_slice = atan(collimation / (2 × SDD)).
     */
    std::array<double, 2> collimationHalfAnglesB() const
    {
        std::array r = {
            std::atan(m_FOVB / m_SDD),
            std::atan(0.5 * m_collimation / m_SDD)
        };
        return r;
    }

    /// @brief Returns a const reference to the bowtie filter for tube A.
    const BowtieFilter& bowtieFilterA() const
    {
        return m_bowtieFilterA;
    }

    /**
     * @brief Replaces the bowtie filter for tube A (copied).
     * @param filter  New `BowtieFilter` instance for tube A.
     */
    void setBowtieFilterA(const BowtieFilter& filter)
    {
        m_bowtieFilterA = filter;
    }

    /// @brief Returns a const reference to the bowtie filter for tube B.
    const BowtieFilter& bowtieFilterB() const
    {
        return m_bowtieFilterB;
    }

    /**
     * @brief Replaces the bowtie filter for tube B (copied).
     * @param filter  New `BowtieFilter` instance for tube B.
     */
    void setBowtieFilterB(const BowtieFilter& filter)
    {
        m_bowtieFilterB = filter;
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
     * @brief Sets the CTDI phantom diameter for the internal calibration simulation [cm]; clamped to ≥ 3 cm.
     * @param d  Phantom diameter [cm]. Typical values: 16 cm (head) or 32 cm (body).
     */
    void setCTDIdiameter(double d) { m_CTDIdiameter = std::max(d, 3.0); }

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

    /// @brief Returns the angular offset of tube B relative to tube A [rad]. Default: 90°.
    double tubeBoffsetAngle() const { return m_tubeBoffsetAngle; }
    /// @brief Sets the angular offset of tube B relative to tube A [rad].
    void setTubeBoffsetAngle(double angle) { m_tubeBoffsetAngle = angle; }
    /// @brief Returns the angular offset of tube B relative to tube A [deg].
    double tubeBoffsetAngleDeg() const { return m_tubeBoffsetAngle * RAD_TO_DEG(); }
    /// @brief Sets the angular offset of tube B relative to tube A [deg].
    void setTubeBoffsetAngleDeg(double angle) { m_tubeBoffsetAngle = angle * DEG_TO_RAD(); }

    /// @brief Returns the angular step between consecutive helix positions [rad].
    double stepAngle() const { return m_stepAngle; }

    /**
     * @brief Sets the angular step between consecutive helix positions [rad].
     *
     * Clamped to a minimum of 0.1° to avoid degenerate helices.
     *
     * @param angle  Desired step size [rad]; the absolute value is used.
     */
    void setStepAngle(double angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }

    /// @brief Returns the angular step between consecutive helix positions [deg].
    double stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG(); }
    /// @brief Sets the angular step between consecutive helix positions [deg].
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    /// @brief Returns a const reference to tube A.
    const Tube& tubeA() const { return m_tubeA; }

    /**
     * @brief Replaces tube A; rebuilds the spectrum caches and relative weights.
     * @param tube  New `Tube` instance (moved in).
     */
    void setTubeA(const Tube&& tube)
    {
        m_tubeA = tube;
        tubeChanged();
    }

    /// @brief Returns a const reference to tube B.
    const Tube& tubeB() const { return m_tubeB; }

    /**
     * @brief Replaces tube B; rebuilds the spectrum caches and relative weights.
     * @param tube  New `Tube` instance (moved in).
     */
    void setTubeB(const Tube&& tube)
    {
        m_tubeB = tube;
        tubeChanged();
    }

    /**
     * @brief Sets the voltage of tube A [kV]; rebuilds the spectrum caches.
     * @param voltage  Tube A voltage [kV].
     */
    void setTubeAVoltage(double voltage)
    {
        m_tubeA.setVoltage(voltage);
        tubeChanged();
    }

    /**
     * @brief Sets the voltage of tube B [kV]; rebuilds the spectrum caches.
     * @param voltage  Tube B voltage [kV].
     */
    void setTubeBVoltage(double voltage)
    {
        m_tubeB.setVoltage(voltage);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle for both tubes simultaneously [rad]; rebuilds the spectrum caches.
     * @param ang  Anode angle [rad].
     */
    void setTubesAnodeAngle(double ang)
    {
        m_tubeA.setAnodeAngle(ang);
        m_tubeB.setAnodeAngle(ang);
        tubeChanged();
    }

    /**
     * @brief Sets the anode angle for both tubes simultaneously [deg]; rebuilds the spectrum caches.
     * @param ang  Anode angle [deg].
     */
    void setTubesAnodeAngleDeg(double ang)
    {
        m_tubeA.setAnodeAngleDeg(ang);
        m_tubeB.setAnodeAngleDeg(ang);
        tubeChanged();
    }

    /**
     * @brief Adds a filtration material to tube A; rebuilds the spectrum caches on success.
     * @param Z   Atomic number of the filter material.
     * @param mm  Filter thickness [mm].
     */
    void addTubeAFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tubeA.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }

    /**
     * @brief Adds a filtration material to tube B; rebuilds the spectrum caches on success.
     * @param Z   Atomic number of the filter material.
     * @param mm  Filter thickness [mm].
     */
    void addTubeBFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tubeB.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }

    /// @brief Removes all filtration materials from both tubes; rebuilds the spectrum caches.
    void clearTubesFiltrationMaterials()
    {
        m_tubeA.clearFiltrationMaterials();
        m_tubeB.clearFiltrationMaterials();
        tubeChanged();
    }

    /// @brief Returns the aluminium half-value layer of tube A's spectrum [mm Al].
    double tubeAAlHalfValueLayer()
    {
        return m_tubeA.mmAlHalfValueLayer();
    }

    /// @brief Returns the aluminium half-value layer of tube B's spectrum [mm Al].
    double tubeBAlHalfValueLayer()
    {
        return m_tubeB.mmAlHalfValueLayer();
    }

    /// @brief Returns the mean photon energy of tube A's spectrum [keV].
    double tubeAMeanSpecterEnergy()
    {
        return m_tubeA.meanSpecterEnergy();
    }

    /// @brief Returns the mean photon energy of tube B's spectrum [keV].
    double tubeBMeanSpecterEnergy()
    {
        return m_tubeB.meanSpecterEnergy();
    }

    /**
     * @brief Sets the energy resolution for both tube spectra; rebuilds the spectrum caches.
     * @param energyResolution  Energy bin width [keV].
     */
    void setTubesEnergyResolution(double energyResolution)
    {
        m_tubeA.setEnergyResolution(energyResolution);
        m_tubeB.setEnergyResolution(energyResolution);
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

    /**
     * @brief Sets the relative mAs for tube A; rebuilds the relative weight normalization.
     *
     * The relative mAs scales the tube's spectral output before the photon-density
     * equalisation between tubes A and B.
     *
     * @param mas  Relative mAs for tube A (dimensionless scale factor).
     */
    void setRelativeMasTubeA(double mas)
    {
        m_relativeMasA = mas;
        tubeChanged();
    }

    /**
     * @brief Sets the relative mAs for tube B; rebuilds the relative weight normalization.
     * @param mas  Relative mAs for tube B (dimensionless scale factor).
     */
    void setRelativeMasTubeB(double mas)
    {
        m_relativeMasB = mas;
        tubeChanged();
    }

    /**
     * @brief Returns the computed relative photon weight for tube A.
     *
     * Normalised so that `weightA + weightB = 2`, preserving the total number of
     * simulated photons while reflecting the mAs and field-size ratio.
     */
    double tubeRelativeWeightA() const
    {
        return m_weightA;
    }

    /**
     * @brief Returns the computed relative photon weight for tube B.
     *
     * Normalised so that `weightA + weightB = 2`.
     */
    double tubeRelativeWeightB() const
    {
        return m_weightB;
    }

    /// @brief Returns the user-supplied relative mAs for tube A.
    double relativeMasTubeA() const
    {
        return m_relativeMasA;
    }

    /// @brief Returns the user-supplied relative mAs for tube B.
    double relativeMasTubeB() const
    {
        return m_relativeMasB;
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
     * @brief Returns a `CTSpiralBeamExposure` for the given dual-tube exposure index.
     *
     * The mapping from @p dualExposureIndex to helix position and tube is:
     * - `i = dualExposureIndex / 2` — helix step index.
     * - `isTubeB = dualExposureIndex % 2 == 1` — even indices → tube A, odd → tube B.
     *
     * Tube B's gantry angle is offset by `m_tubeBoffsetAngle` relative to tube A.
     * The source is placed at `start + direction × dz − beamdir × SDD/2`. The fan
     * half-angle uses the respective tube's FOV; the slice half-angle is shared.
     * Photon weight = `tube_weight × aecFilter(pos) × organWeight`.
     *
     * @param dualExposureIndex  Zero-based index into the interleaved A/B exposure sequence.
     * @return A fully configured `CTSpiralBeamExposure` for the requested position and tube.
     */
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

    /**
     * @brief Returns the dose calibration factor for this dual-energy beam.
     *
     * Runs two sequential internal CTDI simulations — one for tube A and one for tube B —
     * on a shared `CTDIPhantom` of diameter `CTDIdiameter()`, accumulating energy scores
     * without intermediate dose conversion. Then computes:
     *
     *   factor = (CTDIvol × pitch) / CTDIw_calc
     *
     * where CTDIw_calc = (D_centre + 2 × D_periphery) × 10 / (3 × collimation).
     * Each tube uses its own spectrum, bowtie filter, per-tube weight, and the shared
     * organ-AEC filter.
     *
     * @param progress  Optional progress reporter forwarded to the internal transports.
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

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMCTSpiralDualEnergyBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMCTSpiralDualEnergyBeam";
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
     * Writes all scalar fields (start/stop, FOVs, SDD, collimation, angles, weights, mAs,
     * CTDIvol, phantom diameter, pitch, particles per exposure) followed by serialized
     * sub-blocks for tubeA, tubeB, bowtieFilterA, bowtieFilterB, organFilter, and aecFilter.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        Serializer::serialize(m_start, buffer);
        Serializer::serialize(m_stop, buffer);
        Serializer::serialize(m_FOVA, buffer);
        Serializer::serialize(m_FOVB, buffer);
        Serializer::serialize(m_SDD, buffer);
        Serializer::serialize(m_collimation, buffer);
        Serializer::serialize(m_startAngle, buffer);
        Serializer::serialize(m_tubeBoffsetAngle, buffer);
        Serializer::serialize(m_stepAngle, buffer);
        Serializer::serialize(m_weightA, buffer);
        Serializer::serialize(m_weightB, buffer);
        Serializer::serialize(m_relativeMasA, buffer);
        Serializer::serialize(m_relativeMasB, buffer);
        Serializer::serialize(m_CTDIvol, buffer);
        Serializer::serialize(m_CTDIdiameter, buffer);
        Serializer::serialize(m_pitch, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serializeItem(m_tubeA, buffer);
        Serializer::serializeItem(m_tubeB, buffer);
        Serializer::serializeItem(m_bowtieFilterA, buffer);
        Serializer::serializeItem(m_bowtieFilterB, buffer);
        Serializer::serializeItem(m_organFilter, buffer);
        Serializer::serializeItem(m_aecFilter, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `CTSpiralDualEnergyBeam` from a serialized byte buffer.
     *
     * Reads all scalar fields then deserializes the six sub-blocks (tubeA, tubeB,
     * bowtieFilterA, bowtieFilterB, organFilter, aecFilter) in order. Calls
     * `tubeChanged()` after all sub-blocks are restored to rebuild the spectrum caches
     * and relative weights. Returns `std::nullopt` if any sub-block fails.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CTSpiralDualEnergyBeam>` on success, or `std::nullopt` on failure.
     */
    static std::optional<CTSpiralDualEnergyBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        CTSpiralDualEnergyBeam<ENABLETRACKING> item;

        buffer = Serializer::deserialize(item.m_start, buffer);
        buffer = Serializer::deserialize(item.m_stop, buffer);
        buffer = Serializer::deserialize(item.m_FOVA, buffer);
        buffer = Serializer::deserialize(item.m_FOVB, buffer);
        buffer = Serializer::deserialize(item.m_SDD, buffer);
        buffer = Serializer::deserialize(item.m_collimation, buffer);
        buffer = Serializer::deserialize(item.m_startAngle, buffer);
        buffer = Serializer::deserialize(item.m_tubeBoffsetAngle, buffer);
        buffer = Serializer::deserialize(item.m_stepAngle, buffer);
        buffer = Serializer::deserialize(item.m_weightA, buffer);
        buffer = Serializer::deserialize(item.m_weightB, buffer);
        buffer = Serializer::deserialize(item.m_relativeMasA, buffer);
        buffer = Serializer::deserialize(item.m_relativeMasB, buffer);
        buffer = Serializer::deserialize(item.m_CTDIvol, buffer);
        buffer = Serializer::deserialize(item.m_CTDIdiameter, buffer);
        buffer = Serializer::deserialize(item.m_pitch, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        auto name = Serializer::getNameIDTemplate();
        auto item_buffer = Serializer::getEmptyBuffer();

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto tube_optA = Tube::deserialize(item_buffer);
        if (!tube_optA) {
            return std::nullopt;
        } else {
            item.m_tubeA = tube_optA.value();
        }

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto tube_optB = Tube::deserialize(item_buffer);
        if (!tube_optB) {
            return std::nullopt;
        } else {
            item.m_tubeB = tube_optB.value();
        }

        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto bowtie_optA = BowtieFilter::deserialize(item_buffer);
        if (!bowtie_optA) {
            return std::nullopt;
        } else {
            item.m_bowtieFilterA = bowtie_optA.value();
        }
        buffer = Serializer::deserializeItem(name, item_buffer, buffer);
        auto bowtie_optB = BowtieFilter::deserialize(item_buffer);
        if (!bowtie_optB) {
            return std::nullopt;
        } else {
            item.m_bowtieFilterB = bowtie_optB.value();
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
        item.tubeChanged();
        return std::make_optional(item);
    }

protected:
    /**
     * @brief Rebuilds the spectrum caches and normalised relative weights for both tubes.
     *
     * For each tube, queries the current energy bins and spectral weights, constructs a
     * new `SpecterDistribution`, and computes a raw weight:
     *   raw_weight = mAs × spectralSum × FOV
     *
     * The factor FOV corrects for different field sizes (photon density per unit area is
     * equalised). The final normalised weights are:
     *   m_weightA = 2 × raw_A / (raw_A + raw_B)
     *   m_weightB = 2 × raw_B / (raw_A + raw_B)
     *
     * Called automatically after any tube, mAs, or FOV parameter change.
     */
    void tubeChanged()
    {
        // Correcting for different field sizes since photon density per area should be equal
        auto energiesA = m_tubeA.getEnergy();
        auto weightsA = m_tubeA.getSpecter(energiesA, false);
        m_specterA = SpecterDistribution(energiesA, weightsA);
        const double relativeFieldSizeA = m_FOVA; // Collimation is the same for A and B tubes
        const auto weightA = m_relativeMasA * std::reduce(std::execution::par_unseq, weightsA.cbegin(), weightsA.cend(), 0.0) * relativeFieldSizeA;

        auto energiesB = m_tubeB.getEnergy();
        auto weightsB = m_tubeB.getSpecter(energiesB, false);
        m_specterB = SpecterDistribution(energiesB, weightsB);
        const double relativeFieldSizeB = m_FOVB; // Collimation is the same for A and B tubes
        const auto weightB = m_relativeMasB * std::reduce(std::execution::par_unseq, weightsB.cbegin(), weightsB.cend(), 0.0) * relativeFieldSizeB;

        m_weightA = 2 * weightA / (weightA + weightB);
        m_weightB = 2 * weightB / (weightA + weightB);
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 }; ///< Scan start position [cm].
    std::array<double, 3> m_stop = { 0, 0, 0 }; ///< Scan stop position [cm].
    double m_FOVA = 50; ///< Scan field-of-view radius for tube A [cm].
    double m_FOVB = 50; ///< Scan field-of-view radius for tube B [cm].
    double m_SDD = 100; ///< Source-to-detector distance [cm].
    double m_collimation = 1; ///< Total beam collimation width [cm].
    double m_pitch = 1; ///< Helical pitch (table advance / collimation per rotation).
    double m_startAngle = 0; ///< Helix start angle [rad].
    double m_tubeBoffsetAngle = DEG_TO_RAD() * 90; ///< Angular offset of tube B relative to tube A [rad]. Default: 90°.
    double m_stepAngle = DEG_TO_RAD(); ///< Angular step between consecutive helix positions [rad]. Default: 1°.
    double m_weightA = 1; ///< Normalised relative weight for tube A (recomputed by tubeChanged()).
    double m_weightB = 1; ///< Normalised relative weight for tube B (recomputed by tubeChanged()).
    double m_relativeMasA = 1; ///< User-supplied relative mAs for tube A.
    double m_relativeMasB = 1; ///< User-supplied relative mAs for tube B.
    double m_CTDIvol = 1; ///< Prescribed CTDIvol for dose calibration [mGy].
    double m_CTDIdiameter = 32; ///< CTDI phantom diameter for calibration simulation [cm].
    std::uint64_t m_particlesPerExposure = 100; ///< Photon histories per exposure.
    Tube m_tubeA; ///< X-ray tube A model (owned copy).
    Tube m_tubeB; ///< X-ray tube B model (owned copy).
    SpecterDistribution<double> m_specterA; ///< Spectrum cache for tube A, rebuilt by tubeChanged().
    SpecterDistribution<double> m_specterB; ///< Spectrum cache for tube B, rebuilt by tubeChanged().
    CTAECFilter m_aecFilter; ///< Axial AEC filter shared by both tubes.
    BowtieFilter m_bowtieFilterA; ///< Bowtie filter for tube A.
    BowtieFilter m_bowtieFilterB; ///< Bowtie filter for tube B.
    CTOrganAECFilter m_organFilter; ///< Organ AEC filter shared by both tubes.
};
}

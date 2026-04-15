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
 * @brief A single gantry-angle exposure in a sequential (step-and-shoot) CT scan.
 *
 * Returned by `CTSequentialBeam::exposure()`. Each exposure represents one gantry
 * position within a discrete axial rotation. Photons are sampled via
 * `SphereSamplingRectangularField` over the fan-beam solid angle, transformed to
 * world space, then weighted by the combined base weight and the `BowtieFilter`
 * attenuation. Photon energy is drawn from the supplied `SpecterDistribution`.
 *
 * The direction cosines `{normal, scanNormal}` define the in-plane and along-axis
 * directions; `m_dir = cross(normal, scanNormal)` is the central beam direction
 * pointing from source to isocentre.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with
 *                         the start position registered; otherwise returns `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CTSequentialBeamExposure {
public:
    /**
     * @brief Constructs an exposure at the given source position with all sampling parameters.
     *
     * Computes the central beam direction as `cross(dircosines[0], dircosines[1])` and
     * initialises the rectangular-field direction sampler with @p collimationHalfAngles.
     *
     * @param pos                   Source position in world space [cm].
     * @param dircosines            Orthonormal pair `{normal, scanNormal}` defining the
     *                              in-plane and along-axis directions.
     * @param N                     Number of photon histories for this exposure.
     * @param weight                Combined base weight (bowtie weight is applied on top
     *                              during `sampleParticle`).
     * @param collimationHalfAngles `{half_fan_angle, half_cone_angle}` [rad]; fan angle
     *                              from `atan(FOV/SDD)`, cone angle from
     *                              `atan(collimation/(2×SDD))`.
     * @param specter               Non-owning pointer to the photon energy spectrum.
     * @param bowtie                Non-owning pointer to the bowtie filter.
     */
    CTSequentialBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, std::uint64_t N, double weight,
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

    /// @brief Default construction is disabled; all parameters must be supplied.
    CTSequentialBeamExposure() = delete;

    /// @brief Returns the source position of this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the orthonormal direction-cosine pair `{normal, scanNormal}`.
    const std::array<std::array<double, 3>, 2>& directionCosines() const { return m_dirCosines; }

    /// @brief Returns the collimation half-angles `{half_fan, half_cone}` [rad].
    const std::array<double, 2> collimationHalfAngles() const { return m_collimationHalfAngles; }

    /// @brief Returns the number of photon histories for this exposure.
    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * Draws a direction from `SphereSamplingRectangularField`, computes the fan
     * angle `angx = atan(dir_x / dir_z)` in local space, evaluates the bowtie
     * weight at that angle, then transforms the direction to world space via
     * `changeBasis`. Photon energy is drawn from the spectrum. The final weight
     * is `m_weight × bowtie_weight`.
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

    /// @brief Returns the base photon weight for this exposure (before bowtie modulation).
    double weight() const
    {
        return m_weight;
    }

protected:
private:
    SphereSamplingRectangularField m_directionSampler;                                          ///< Samples directions within the fan-beam solid angle.
    std::array<double, 3> m_pos = { 0, 0, 0 };                                                 ///< Source position [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 };                                                 ///< Central beam direction (cross product of direction cosines).
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };     ///< Orthonormal basis: {in-plane normal, scan axis}.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };                                  ///< Half-angles {fan, cone} [rad].
    std::uint64_t m_NParticles = 100;                                                            ///< Number of photon histories.
    double m_weight = 1;                                                                          ///< Base photon weight.
    const SpecterDistribution<double>* m_specter = nullptr;                                      ///< Non-owning pointer to the photon energy spectrum.
    const BowtieFilter* m_bowtieFilter = nullptr;                                                ///< Non-owning pointer to the bowtie filter.
};

/**
 * @brief A sequential (step-and-shoot) CT beam satisfying `BeamType` and `SerializeItemType`.
 *
 * Models a CT scanner that acquires one full 360° axial rotation at each table position
 * before stepping the couch to the next slice. There is no helical pitch; the table
 * advances by `m_sliceSpacing` between rotations.
 *
 * **Geometry**: the source orbits an isocentre defined by `m_position`. At gantry angle
 * `θ = i × stepAngle` the in-plane normal is `rotate(n₀, scanNormal, startAngle + θ)`
 * and the source position is `isocentre + scanNormal × (sliceSpacing × sliceNumber)
 * − beamdir × SDD/2`, where `sliceNumber = floor(θ / 2π)`.
 *
 * **Exposure count**: `numberOfExposures = floor(numberOfSlices × 2π / stepAngle)`.
 * Each exposure rotates the gantry by one `stepAngle`; a new slice begins every
 * `floor(2π / stepAngle)` exposures.
 *
 * **Weighting**: each exposure weight = `m_weight × organAECFilter(θ)`. Bowtie
 * attenuation is applied inside `CTSequentialBeamExposure::sampleParticle`.
 *
 * **Calibration** (`calibrationFactor()`): builds a `CTDIPhantom`, transports a
 * `CTDIBeam`, and returns `CTDIw_target / CTDIw_simulated` where
 * `CTDIw = (D_centre + 2 × D_periphery) × 10 / (3 × collimation)`.
 *
 * The `SpecterDistribution` cache is rebuilt by `tubeChanged()` whenever any tube
 * parameter is modified.
 *
 * Satisfies the `BeamType` concept and the `SerializeItemType` concept.
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CTSequentialBeam {
public:
    /**
     * @brief Constructs a sequential CT beam at the given start position and scan axis.
     *
     * @p scan_normal is normalised internally. Each entry in @p filtrationMaterials
     * is added to the X-ray tube via `addFiltrationMaterial`; `tubeChanged()` is
     * called once after construction to build the initial `SpecterDistribution`.
     *
     * @param start_pos           Isocentre position for the first slice [cm]. Default: origin.
     * @param scan_normal         Unit vector along the scan (couch-travel) axis. Default: +z.
     * @param filtrationMaterials Map of `{atomic number Z → thickness [mm]}` added to the tube.
     */
    CTSequentialBeam(
        const std::array<double, 3>& start_pos = { 0, 0, 0 },
        const std::array<double, 3>& scan_normal = { 0, 0, 1 },
        const std::map<std::size_t, double>& filtrationMaterials = { })
        : m_position(start_pos)
        , m_scanNormal(scan_normal)
    {
        vectormath::normalize(m_scanNormal);
        for (const auto [Z, mm] : filtrationMaterials)
            m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }

    /**
     * @brief Returns the total number of gantry-angle exposures across all slices.
     *
     * Computed as `floor(numberOfSlices × 2π / stepAngle)`. Each full rotation
     * spans `floor(2π / stepAngle)` exposures; the table advances by `sliceSpacing`
     * every rotation.
     */
    std::uint64_t numberOfExposures() const
    {
        const auto total_rot_angle = m_numberOfSlices * PI_VAL() * 2;
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles;
    }

    /// @brief Returns the total number of photon histories across all exposures.
    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }
    /// @brief Returns the number of photon histories per gantry-angle exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }
    /**
     * @brief Sets the number of photon histories per gantry-angle exposure.
     * @param n  Histories per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /// @brief Returns the isocentre position of the first slice [cm].
    const std::array<double, 3>& position() const { return m_position; }
    /// @brief Returns the scan-axis unit vector.
    const std::array<double, 3>& scanNormal() const { return m_scanNormal; }
    /**
     * @brief Sets the isocentre position of the first slice [cm].
     * @param position  3-D position in world space [cm].
     */
    void setPosition(const std::array<double, 3>& position)
    {
        m_position = position;
    }
    /**
     * @brief Sets the scan (couch-travel) axis; normalised internally.
     * @param normal  Direction vector (need not be a unit vector).
     */
    void setScanNormal(const std::array<double, 3>& normal)
    {
        m_scanNormal = vectormath::normalized(normal);
    }
    /// @brief Returns the couch step between consecutive slices [cm].
    double sliceSpacing() const
    {
        return m_sliceSpacing;
    }
    /**
     * @brief Sets the couch step between consecutive slices [cm]; clamped to ≥ 0.
     * @param spacing  Couch step [cm].
     */
    void setSliceSpacing(double spacing)
    {
        m_sliceSpacing = std::max(0.0, spacing);
    }
    /// @brief Returns the number of discrete axial slices to acquire.
    std::uint64_t numberOfSlices() const { return m_numberOfSlices; }
    /**
     * @brief Sets the number of discrete axial slices; clamped to ≥ 1.
     * @param N  Number of slices.
     */
    void setNumberOfSlices(std::uint64_t N)
    {
        m_numberOfSlices = std::max(N, std::uint64_t { 1 });
    }

    /// @brief Returns the beam collimation (total detector row width) [cm].
    double collimation() const
    {
        return m_collimation;
    }
    /**
     * @brief Sets the beam collimation [cm]; minimum 1 mm (0.1 cm).
     * @param coll_cm  Collimation [cm]; absolute value is clamped to ≥ 0.1.
     */
    void setCollimation(double coll_cm)
    {
        // collimation must be larger than 1 mm (0.1 cm)
        m_collimation = std::max(std::abs(coll_cm), 0.1);
    }

    /**
     * @brief Returns the collimation half-angles `{half_fan, half_cone}` [rad].
     *
     * Computed as `{atan(FOV/SDD), atan(collimation/(2×SDD))}`.
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

    /// @brief Returns the scan field-of-view radius [cm].
    double scanFieldOfView() const { return m_FOV; }
    /**
     * @brief Sets the scan field-of-view radius [cm]; clamped to ≥ 1 cm.
     * @param fov_cm  FOV radius [cm].
     */
    void setScanFieldOfView(double fov_cm)
    {
        m_FOV = std::max(std::abs(fov_cm), 1.0);
    }

    /// @brief Returns the bowtie filter applied to all exposures.
    const BowtieFilter& bowtieFilter() const
    {
        return m_bowtieFilter;
    }
    /**
     * @brief Replaces the bowtie filter.
     * @param filter  New `BowtieFilter` (copied by value).
     */
    void setBowtieFilter(const BowtieFilter& filter)
    {
        m_bowtieFilter = filter;
    }

    /**
     * @brief Sets the target CTDIw used for dose calibration [mGy].
     * @param ctdi  Target CTDIw value [mGy].
     */
    void setCTDIw(double ctdi) { m_CTDIw = ctdi; }
    /// @brief Returns the target CTDIw [mGy].
    double CTDIw() const { return m_CTDIw; }
    /**
     * @brief Sets the CTDI phantom diameter used in calibration [cm].
     * @param d  Phantom diameter [cm] (typically 16 cm head or 32 cm body).
     */
    void setCTDIdiameter(double d) { m_CTDIdiameter = d; }
    /// @brief Returns the CTDI phantom diameter [cm].
    double CTDIdiameter() const { return m_CTDIdiameter; }

    /// @brief Returns the gantry start angle [rad].
    double startAngle() const { return m_startAngle; }
    /**
     * @brief Sets the gantry start angle [rad].
     * @param angle  Start angle [rad].
     */
    void setStartAngle(double angle) { m_startAngle = angle; }
    /// @brief Returns the gantry start angle [deg].
    double startAngleDeg() const { return m_startAngle * RAD_TO_DEG(); }
    /**
     * @brief Sets the gantry start angle [deg].
     * @param angle  Start angle [deg].
     */
    void setStartAngleDeg(double angle) { m_startAngle = angle * DEG_TO_RAD(); }

    /// @brief Returns the angular step between consecutive exposures [rad].
    double stepAngle() const { return m_stepAngle; }
    /**
     * @brief Sets the angular step between consecutive exposures [rad];
     *        clamped to ≥ 0.1°.
     * @param angle  Step angle [rad].
     */
    void setStepAngle(double angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD() / 10);
    }
    /// @brief Returns the angular step between consecutive exposures [deg].
    double stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG(); }
    /**
     * @brief Sets the angular step between consecutive exposures [deg].
     * @param angle  Step angle [deg].
     */
    void setStepAngleDeg(double angle) { setStepAngle(angle * DEG_TO_RAD()); }

    /// @brief Returns the X-ray tube configuration.
    const Tube& tube() const { return m_tube; }
    /**
     * @brief Replaces the X-ray tube and rebuilds the spectrum cache.
     * @param tube  New `Tube` (moved); `tubeChanged()` is called afterwards.
     */
    void setTube(const Tube&& tube)
    {
        m_tube = tube;
        tubeChanged();
    }
    /**
     * @brief Sets the tube peak voltage [kV] and rebuilds the spectrum cache.
     * @param voltage  Peak voltage [kV].
     */
    void setTubeVoltage(double voltage)
    {
        m_tube.setVoltage(voltage);
        tubeChanged();
    }
    /**
     * @brief Sets the anode angle [rad] and rebuilds the spectrum cache.
     * @param ang  Anode angle [rad].
     */
    void setTubeAnodeAngle(double ang)
    {
        m_tube.setAnodeAngle(ang);
        tubeChanged();
    }
    /**
     * @brief Sets the anode angle [deg] and rebuilds the spectrum cache.
     * @param ang  Anode angle [deg].
     */
    void setTubeAnodeAngleDeg(double ang)
    {
        m_tube.setAnodeAngleDeg(ang);
        tubeChanged();
    }
    /**
     * @brief Adds a filtration material and rebuilds the spectrum cache if accepted.
     * @param Z   Atomic number of the filter material.
     * @param mm  Filter thickness [mm].
     */
    void addTubeFiltrationMaterial(std::size_t Z, double mm)
    {
        auto success = m_tube.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    /**
     * @brief Removes all filtration materials and rebuilds the spectrum cache.
     */
    void clearTubeFiltrationMaterials()
    {
        m_tube.clearFiltrationMaterials();
        tubeChanged();
    }
    /**
     * @brief Returns the aluminium half-value layer of the current tube spectrum [mm Al].
     * @return Half-value layer [mm Al].
     */
    double tubeAlHalfValueLayer()
    {
        return m_tube.mmAlHalfValueLayer();
    }
    /**
     * @brief Returns the mean photon energy of the current tube spectrum [keV].
     * @return Mean spectrum energy [keV].
     */
    double tubeMeanSpecterEnergy()
    {
        return m_tube.meanSpecterEnergy();
    }
    /**
     * @brief Sets the tube energy resolution and rebuilds the spectrum cache.
     * @param energyResolution  Energy bin resolution [keV].
     */
    void setTubeEnergyResolution(double energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    /**
     * @brief Returns a mutable reference to the organ angular AEC filter.
     *
     * The organ AEC filter applies three-zone angular modulation (anterior,
     * lateral, posterior) to the photon weight within each rotation.
     */
    CTOrganAECFilter& organAECFilter()
    {
        return m_organFilter;
    }
    /// @brief Returns a const reference to the organ angular AEC filter.
    const CTOrganAECFilter& organAECFilter() const
    {
        return m_organFilter;
    }

    /**
     * @brief Returns a `CTSequentialBeamExposure` for the given gantry-angle index.
     *
     * The slice number is derived from `floor(i × stepAngle / 2π)`; the source is
     * offset along the scan axis by `sliceSpacing × sliceNumber`. The in-plane
     * normal is rotated by `startAngle + i × stepAngle` around `scanNormal`.
     * The organ AEC weight is evaluated at the cumulative rotation angle.
     *
     * @param i  Exposure index in [0, numberOfExposures()).
     * @return A fully configured `CTSequentialBeamExposure`.
     */
    CTSequentialBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        constexpr auto pi2 = PI_VAL() * 2;
        const auto angle = i * m_stepAngle;
        const auto sliceNumber = static_cast<int>(angle / pi2);

        // finding normal vector to direction
        const auto normal_ind = vectormath::argmin3(m_scanNormal);
        std::array<double, 3> normal = { 0, 0, 0 };
        normal[normal_ind] = 1;
        normal = vectormath::normalized(vectormath::cross(normal, m_scanNormal));
        normal = vectormath::rotate(normal, m_scanNormal, m_startAngle + angle);

        const auto beamdir = vectormath::cross(normal, m_scanNormal);

        const std::array<std::array<double, 3>, 2> cosines = { normal, m_scanNormal };

        auto pos = vectormath::add(vectormath::add(m_position, vectormath::scale(m_scanNormal, m_sliceSpacing * sliceNumber)), vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis
        const auto angx = std::atan(m_FOV / m_SDD);
        const auto angy = std::atan(0.5 * m_collimation / m_SDD);

        std::array<double, 2> angles = { angx, angy };

        const auto organWeight = m_organFilter.useFilter() ? m_organFilter(angle) : 1.0;

        CTSequentialBeamExposure<ENABLETRACKING> exp(pos, cosines, m_particlesPerExposure, m_weight * organWeight, angles, &m_specter, &m_bowtieFilter);
        return exp;
    }

    /**
     * @brief Computes the dose calibration factor via CTDI phantom simulation.
     *
     * Builds a `CTDIPhantom<5,1>` of diameter `m_CTDIdiameter`, transports a
     * `CTDIBeam` (one full 360° rotation) through it, then computes:
     *
     *   CTDIw_calc = (D_centre + 2 × D_periphery) × 10 / (3 × collimation)
     *
     * Returns `CTDIw_target / CTDIw_calc` so that multiplying all photon weights
     * by this factor reproduces the user-specified `m_CTDIw`.
     *
     * @param progress  Optional progress reporter; passed through to `Transport`.
     * @return Calibration factor = `CTDIw_target / CTDIw_simulated`.
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

        return m_CTDIw / ctdiw_calc;
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMCTSeqBeam" padded with spaces.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMCTSeqBeam";
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
     * Writes scalar parameters (position, scanNormal, FOV, SDD, collimation,
     * startAngle, stepAngle, weight, CTDIw, CTDIdiameter, sliceSpacing,
     * numberOfSlices, particlesPerExposure) followed by three sub-item blocks:
     * `Tube`, `BowtieFilter`, and `CTOrganAECFilter`.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        Serializer::serialize(m_position, buffer);
        Serializer::serialize(m_scanNormal, buffer);
        Serializer::serialize(m_FOV, buffer);
        Serializer::serialize(m_SDD, buffer);
        Serializer::serialize(m_collimation, buffer);
        Serializer::serialize(m_startAngle, buffer);
        Serializer::serialize(m_stepAngle, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_CTDIw, buffer);
        Serializer::serialize(m_CTDIdiameter, buffer);
        Serializer::serialize(m_sliceSpacing, buffer);
        Serializer::serialize(m_numberOfSlices, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);
        Serializer::serializeItem(m_tube, buffer);
        Serializer::serializeItem(m_bowtieFilter, buffer);
        Serializer::serializeItem(m_organFilter, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `CTSequentialBeam` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`, calls `tubeChanged()` after
     * restoring the tube to rebuild the spectrum cache. Returns `nullopt` if
     * any sub-item (Tube, BowtieFilter, or CTOrganAECFilter) fails to deserialize.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CTSequentialBeam>` on success, `nullopt` on failure.
     */
    static std::optional<CTSequentialBeam<ENABLETRACKING>> deserialize(std::span<const char> buffer)
    {
        CTSequentialBeam<ENABLETRACKING> item;

        buffer = Serializer::deserialize(item.m_position, buffer);
        buffer = Serializer::deserialize(item.m_scanNormal, buffer);
        buffer = Serializer::deserialize(item.m_FOV, buffer);
        buffer = Serializer::deserialize(item.m_SDD, buffer);
        buffer = Serializer::deserialize(item.m_collimation, buffer);
        buffer = Serializer::deserialize(item.m_startAngle, buffer);
        buffer = Serializer::deserialize(item.m_stepAngle, buffer);
        buffer = Serializer::deserialize(item.m_weight, buffer);
        buffer = Serializer::deserialize(item.m_CTDIw, buffer);
        buffer = Serializer::deserialize(item.m_CTDIdiameter, buffer);
        buffer = Serializer::deserialize(item.m_sliceSpacing, buffer);
        buffer = Serializer::deserialize(item.m_numberOfSlices, buffer);
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
        auto aec_opt = CTOrganAECFilter::deserialize(item_buffer);
        if (!aec_opt) {
            return std::nullopt;
        } else {
            item.m_organFilter = aec_opt.value();
        }

        return std::make_optional(item);
    }

protected:
    /**
     * @brief Rebuilds the `SpecterDistribution` cache after any tube parameter change.
     *
     * Queries the tube for its energy bins and unfiltered spectral weights, then
     * constructs a new `SpecterDistribution` for sampling in `exposure()`.
     * Called by all `setTube*` methods and the constructor.
     */
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<double, 3> m_position = { 0, 0, 0 };   ///< Isocentre position of the first slice [cm].
    std::array<double, 3> m_scanNormal = { 0, 0, 1 };  ///< Unit vector along the scan (couch-travel) axis.
    double m_FOV = 50;                                  ///< Scan field-of-view radius [cm].
    double m_SDD = 100;                                 ///< Source-to-detector distance [cm].
    double m_collimation = 1;                           ///< Total beam collimation (detector row width) [cm].
    double m_startAngle = 0;                            ///< Gantry start angle [rad].
    double m_stepAngle = 0.018;                         ///< Angular step between exposures [rad] (~1°).
    double m_weight = 1;                                ///< Base photon weight applied to all exposures.
    double m_CTDIw = 1;                                 ///< Target CTDIw for dose calibration [mGy].
    double m_CTDIdiameter = 32;                         ///< CTDI phantom diameter [cm].
    double m_sliceSpacing = 0;                          ///< Couch step between consecutive slices [cm].
    std::uint64_t m_numberOfSlices = 1;                 ///< Number of discrete axial slices to acquire.
    std::uint64_t m_particlesPerExposure = 100;         ///< Photon histories per gantry-angle exposure.
    Tube m_tube;                                        ///< X-ray tube configuration.
    SpecterDistribution<double> m_specter;              ///< Cached photon energy spectrum (rebuilt by tubeChanged()).
    BowtieFilter m_bowtieFilter;                        ///< Bowtie filter for fan-angle intensity modulation.
    CTOrganAECFilter m_organFilter;                     ///< Organ angular AEC filter for angular weight modulation.
};
}
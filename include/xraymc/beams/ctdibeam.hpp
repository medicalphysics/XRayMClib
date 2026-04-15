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

#include "xraymc/beams/filters/bowtiefilter.hpp"
#include "xraymc/beams/filters/ctorganaecfilter.hpp"
#include "xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {

/**
 * @brief A single CT gantry exposure used internally by `CTDIBeam`.
 *
 * Represents the X-ray tube at one gantry angle of a full-rotation CTDI scan. The
 * source position is placed at distance SDD/2 from the isocentre along the beam
 * direction. Photon energies are drawn from a `SpecterDistribution`, fan-angle
 * directions from a `SphereSamplingRectangularField`, and per-photon weights are
 * modulated by a `BowtieFilter`.
 *
 * @tparam ENABLETRACKING  If true, `sampleParticle()` returns `ParticleTrack` with the
 *                         start position registered; otherwise returns `Particle`. Default: false.
 */
template <bool ENABLETRACKING = false>
class CTDIBeamExposure {
public:
    /**
     * @brief Constructs a single CT exposure at the given gantry angle.
     *
     * Rotates the default in-plane direction cosine by @p angle around the z-axis to
     * orient the source, then places it at SDD/2 from the isocentre.
     *
     * @param angle                  Gantry angle for this exposure [rad].
     * @param SDD                    Source-to-detector distance [cm]; source is at SDD/2 from isocentre.
     * @param historiesPerExposure   Number of photon histories to simulate in this exposure.
     * @param collimationHalfAngles  `{half_angle_x, half_angle_y}` — rectangular field
     *                               half-angles [rad] passed to `SphereSamplingRectangularField`.
     * @param specter                Non-owning pointer to the energy spectrum sampler.
     * @param bowtie                 Non-owning pointer to the bowtie filter.
     * @param weight                 Base statistical weight for each sampled photon. Default: 1.
     */
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

    /// @brief Default constructor is deleted — an exposure must always have geometry.
    CTDIBeamExposure() = delete;

    /// @brief Returns the source position for this exposure [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the number of photon histories in this exposure.
    std::uint64_t numberOfParticles() const
    {
        return m_Nparticles;
    }

    /**
     * @brief Samples a single photon from this exposure.
     *
     * 1. Draws a fan-angle direction from `SphereSamplingRectangularField`.
     * 2. Derives the in-plane fan angle and looks up the bowtie weight.
     * 3. Transforms the local direction into world space via `changeBasis`.
     * 4. Draws an energy from the spectrum sampler.
     * 5. Returns a `Particle` (or `ParticleTrack` when `ENABLETRACKING` is true)
     *    with weight = `m_weight × bowtie_weight`.
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

private:
    SphereSamplingRectangularField m_directionSampler;          ///< Samples fan-angle directions within the collimated field.
    std::uint64_t m_Nparticles = 1;                             ///< Number of photon histories for this exposure.
    std::array<double, 3> m_pos = { 0, 0, 0 };                 ///< Source position in world space [cm].
    std::array<double, 3> m_dir = { 0, 1, 0 };                 ///< Central beam direction (towards isocentre).
    std::array<double, 3> m_dirCosineX = { 1, 0, 0 };          ///< In-plane direction cosine (rotated by gantry angle).
    std::array<double, 3> m_dirCosineY = { 0, 0, 1 };          ///< Slice-direction (z-axis) direction cosine.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 };  ///< Rectangular collimation half-angles {x, y} [rad].
    double m_weight = 1;                                         ///< Base photon weight (before bowtie modulation).
    const SpecterDistribution<double>* m_specter = nullptr;     ///< Non-owning pointer to the energy spectrum.
    const BowtieFilter* m_bowtieFilter = nullptr;               ///< Non-owning pointer to the bowtie filter.
};

/**
 * @brief A full-rotation CT beam for CTDI dosimetry, satisfying `BeamType`.
 *
 * Simulates a complete 360° axial CT rotation by distributing exposures evenly around
 * the gantry. Each exposure is a `CTDIBeamExposure` at angle `i × angleStep`, with
 * photon energies drawn from a `SpecterDistribution`, fan-angle directions sampled by
 * `SphereSamplingRectangularField`, and per-photon weights modulated by a `BowtieFilter`.
 *
 * An optional `CTOrganAECFilter` can apply angle-dependent tube-current modulation:
 * when `organFilter.useFilter()` is true, the exposure weight is multiplied by the
 * filter weight at that gantry angle.
 *
 * `calibrationFactor()` returns 1 — dose conversion is assumed to be handled
 * externally (e.g. by scaling to a measured CTDI value).
 *
 * @tparam ENABLETRACKING  If true, exposures return `ParticleTrack`; otherwise `Particle`.
 *                         Default: false.
 */
template <bool ENABLETRACKING = false>
class CTDIBeam {
public:
    /**
     * @brief Constructs a full-rotation CT beam with organ AEC modulation.
     *
     * @param angleStep              Angular step between consecutive exposures [rad].
     *                               `numberOfExposures()` = floor(2π / angleStep).
     * @param SDD                    Source-to-detector distance [cm].
     * @param collimationHalfAngles  `{half_x, half_y}` collimation half-angles [rad].
     * @param particlesPerExposure   Number of photon histories per gantry position.
     * @param specter                Energy spectrum sampler (copied into the beam).
     * @param bowtie                 Bowtie filter (copied into the beam).
     * @param organFilter            Organ AEC filter for angle-dependent weight modulation (copied).
     * @param weight                 Base statistical weight. Default: 1.
     */
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

    /**
     * @brief Constructs a full-rotation CT beam without organ AEC modulation.
     *
     * Uses a default-constructed (disabled) `CTOrganAECFilter`.
     *
     * @param angleStep              Angular step between consecutive exposures [rad].
     * @param SDD                    Source-to-detector distance [cm].
     * @param collimationHalfAngles  `{half_x, half_y}` collimation half-angles [rad].
     * @param particlesPerExposure   Number of photon histories per gantry position.
     * @param specter                Energy spectrum sampler (copied into the beam).
     * @param bowtie                 Bowtie filter (copied into the beam).
     * @param weight                 Base statistical weight. Default: 1.
     */
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

    /**
     * @brief Returns the number of gantry exposures in one full rotation.
     *
     * Computed as floor(2π / angleStep).
     */
    std::uint64_t numberOfExposures() const
    {
        constexpr auto pi2 = PI_VAL() * 2;
        return static_cast<std::uint64_t>(pi2 / m_angleStep);
    }

    /// @brief Returns the total number of photon histories across all exposures.
    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }

    /**
     * @brief Returns the `CTDIBeamExposure` for gantry position @p i.
     *
     * Computes the gantry angle as `i × m_angleStep`. If the organ AEC filter is
     * enabled, the exposure weight is scaled by the filter weight at that angle.
     *
     * @param i  Zero-based exposure index; should be in [0, numberOfExposures()).
     * @return A fully configured `CTDIBeamExposure` for the requested gantry position.
     */
    CTDIBeamExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        const auto angle = i * m_angleStep;
        const auto organWeight = m_organFilter.useFilter() ? m_organFilter(angle) : 1.0;
        return CTDIBeamExposure<ENABLETRACKING>(angle, m_sdd, m_particlesPerExposure, m_collimationHalfAngles, &m_specter, &m_bowtieFilter, m_weight * organWeight);
    }

    /**
     * @brief Returns the dose calibration factor (always 1 for CTDI beams).
     *
     * CTDI beams are intended for relative dosimetry; the calibration to absolute
     * dose units is applied externally. The `progress` parameter is accepted for
     * `BeamType` interface compatibility but is not used.
     *
     * @param progress  Unused.
     * @return 1.0.
     */
    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        return 1;
    }

private:
    double m_angleStep = 0;                                    ///< Angular step between exposures [rad].
    double m_sdd = 1;                                          ///< Source-to-detector distance [cm].
    double m_weight = 1;                                       ///< Base photon weight.
    std::array<double, 2> m_collimationHalfAngles = { 0, 0 }; ///< Rectangular collimation half-angles {x, y} [rad].
    std::uint64_t m_particlesPerExposure = 1;                  ///< Photon histories per gantry position.
    SpecterDistribution<double> m_specter;                     ///< Energy spectrum sampler (owned copy).
    BowtieFilter m_bowtieFilter;                               ///< Bowtie filter (owned copy).
    CTOrganAECFilter m_organFilter;                            ///< Organ AEC filter (owned copy; disabled by default).
};
}
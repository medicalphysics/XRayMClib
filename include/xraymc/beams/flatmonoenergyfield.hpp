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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "xraymc/constants.hpp"
#include "xraymc/floating.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/transportprogress.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <span>

namespace xraymc {

/**
 * @brief A single exposure event for a flat, monoenergetic radiation field.
 *
 * Represents one snapshot of the beam geometry and energy used to sample
 * photon particles uniformly across a rectangular field. Created by
 * `FlatMonoEnergyField::exposure()`.
 *
 * @tparam ENABLETRACKING  When true, sampled particles carry full track
 *                         history (`ParticleTrack`); otherwise plain
 *                         `Particle` objects are returned.
 */
template <bool ENABLETRACKING = false>
class FlatMonoEnergyFieldExposure {
public:
    /**
     * @brief Constructs an exposure from explicit beam parameters.
     *
     * @param pos      Centre of the field in world coordinates [cm].
     * @param cosines  Two orthonormal in-plane direction cosines {x-axis, y-axis}.
     * @param lenghtX  Half-width of the field along the x-axis [cm].
     * @param lenghtY  Half-width of the field along the y-axis [cm].
     * @param energy   Photon energy [keV].
     * @param weight   Statistical weight assigned to every sampled particle.
     * @param N        Number of particles to emit in this exposure.
     */
    FlatMonoEnergyFieldExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& cosines, double lenghtX, double lenghtY, double energy, double weight, std::uint64_t N)
        : m_energy(energy)
        , m_weight(weight)
        , m_pos(pos)
        , m_NParticles(N)
    {
        m_dir = vectormath::cross(cosines);
        m_dirXY[0] = vectormath::scale(cosines[0], lenghtX);
        m_dirXY[1] = vectormath::scale(cosines[1], lenghtY);
    }

    /// @brief Returns the centre of the field in world coordinates [cm].
    const std::array<double, 3>& position() const { return m_pos; }

    /// @brief Returns the number of particles to emit in this exposure.
    std::uint64_t numberOfParticles() const { return m_NParticles; }

    /**
     * @brief Samples a single photon particle from a uniform distribution over the field.
     *
     * Draws two uniform random numbers in [-1, 1] to select a position within the
     * rectangular field, then constructs a particle (or `ParticleTrack` when tracking
     * is enabled) with the fixed beam direction and energy.
     *
     * @param state  Random-number generator state (modified in place).
     * @return A `Particle` or `ParticleTrack` with position, direction, energy, and weight set.
     */
    auto sampleParticle(RandomState& state) const noexcept
    {

        const auto r1 = state.randomUniform(-1.0, 1.0);
        const auto r2 = state.randomUniform(-1.0, 1.0);

        const auto pos = vectormath::add(vectormath::scale(m_dirXY[0], r1), vectormath::scale(m_dirXY[1], r2));

        if constexpr (ENABLETRACKING) {
            ParticleTrack p = { .pos = vectormath::add(pos, m_pos),
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            p.registerPosition();
            return p;
        } else {
            Particle p = { .pos = vectormath::add(pos, m_pos),
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            return p;
        }
    }

private:
    double m_energy = 60; ///< Photon energy [keV].
    double m_weight = 1; ///< Statistical weight of each particle.
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirXY = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction vectors scaled by half-widths [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Beam direction unit vector.
    std::uint64_t m_NParticles = 100; ///< Number of particles in this exposure.
};

/**
 * @brief A flat, monoenergetic photon beam source with a rectangular field.
 *
 * Models a parallel beam of photons with a single energy that are emitted
 * uniformly from a rectangular aperture. The beam is defined by a centre
 * position, two in-plane direction cosines, field half-widths in X and Y,
 * and a photon energy. Multiple exposures can be generated; each exposure
 * independently samples particles from the same geometry.
 *
 * Dose calibration is supported via a prescribed air-KERMA value.
 *
 * @tparam ENABLETRACKING  When true, all sampled particles carry full track
 *                         history; otherwise plain `Particle` objects are used.
 */
template <bool ENABLETRACKING = false>
class FlatMonoEnergyField {
public:
    /**
     * @brief Constructs a beam with default settings.
     *
     * @param pos     Centre of the field in world coordinates [cm]. Defaults to origin.
     * @param energy  Photon energy [keV]. Defaults to 60 keV.
     */
    FlatMonoEnergyField(const std::array<double, 3>& pos = { 0, 0, 0 }, double energy = 60)
        : m_energy(energy)
        , m_pos(pos)
    {
    }

    /**
     * @brief Sets the photon energy, clamped to the valid simulation range.
     * @param energy  Desired energy [keV]; the absolute value is used.
     */
    void setEnergy(double energy)
    {
        m_energy = std::clamp(std::abs(energy), MIN_ENERGY(), MAX_ENERGY());
    }

    /// @brief Returns the photon energy [keV].
    double energy() const { return m_energy; }

    /// @brief Returns the number of exposures.
    std::uint64_t numberOfExposures() const { return m_Nexposures; }

    /**
     * @brief Sets the number of exposures; enforces a minimum of one.
     * @param n  Desired exposure count.
     */
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    /**
     * @brief Sets the number of particles emitted per exposure.
     * @param n  Desired particle count per exposure.
     */
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    /// @brief Returns the number of particles emitted per exposure.
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    /// @brief Returns the total number of particles across all exposures.
    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }

    /**
     * @brief Sets the centre position of the field.
     * @param pos  New source position in world coordinates [cm].
     */
    void setPosition(const std::array<double, 3>& pos)
    {
        m_pos = pos;
    }

    /// @brief Returns the centre position of the field in world coordinates [cm].
    const std::array<double, 3>& position() const
    {
        return m_pos;
    }

    /// @brief Returns the beam propagation direction as the cross product of the two in-plane cosines.
    const std::array<double, 3> direction() const
    {
        return vectormath::cross(m_dirCosines[0], m_dirCosines[1]);
    }

    /// @brief Returns the two normalised in-plane direction cosines {x-axis, y-axis}.
    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

    /**
     * @brief Sets both in-plane direction cosines and normalises them.
     * @param dir  Two-element array containing {x-cosine, y-cosine}.
     */
    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    /**
     * @brief Sets both in-plane direction cosines from separate vectors.
     * @param xdir  Direction along the field x-axis (will be normalised).
     * @param ydir  Direction along the field y-axis (will be normalised).
     */
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dirCosines[0] = vectormath::normalized(xdir);
        m_dirCosines[1] = vectormath::normalized(ydir);
    }

    /**
     * @brief Sets the field half-width along the x-axis.
     * @param x  Half-width [cm]; the absolute value is stored.
     */
    void setLenghtX(double x)
    {
        m_lenghtX = std::abs(x);
    }

    /**
     * @brief Sets the field half-width along the y-axis.
     * @param y  Half-width [cm]; the absolute value is stored.
     */
    void setLenghtY(double y)
    {
        m_lenghtY = std::abs(y);
    }

    /**
     * @brief Sets equal half-widths in both x and y directions.
     * @param l  Half-width [cm] applied to both axes.
     */
    void setLenght(double l)
    {
        setLenghtX(l);
        setLenghtY(l);
    }

    /**
     * @brief Sets independent half-widths in x and y directions.
     * @param x  Half-width along the x-axis [cm].
     * @param y  Half-width along the y-axis [cm].
     */
    void setLenght(double x, double y)
    {
        setLenghtX(x);
        setLenghtY(y);
    }

    /// @brief Returns the field half-width along the x-axis [cm].
    double lenghtX() const
    {
        return m_lenghtX;
    }

    /// @brief Returns the field half-width along the y-axis [cm].
    double lenghtY() const
    {
        return m_lenghtY;
    }

    /**
     * @brief Sets the statistical weight assigned to every sampled particle.
     * @param weight  Particle weight. Defaults to 1.
     */
    void setParticleWeight(double weight = 1)
    {
        m_weight = weight;
    }

    /**
     * @brief Creates an exposure object for the given exposure index.
     *
     * All exposures share the same geometry; the index parameter is accepted
     * for interface compatibility but is not used.
     *
     * @param i  Exposure index (unused).
     * @return A `FlatMonoEnergyFieldExposure` configured with the current beam parameters.
     */
    FlatMonoEnergyFieldExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        FlatMonoEnergyFieldExposure<ENABLETRACKING> exp(m_pos, m_dirCosines, m_lenghtX, m_lenghtY, m_energy, m_weight, m_particlesPerExposure);
        return exp;
    }

    /**
     * @brief Sets the prescribed air-KERMA for dose calibration.
     * @param k  Air-KERMA value [mGy]; the absolute value is stored.
     */
    void setAirKerma(double k)
    {
        m_airKerma = std::abs(k);
    }

    /// @brief Returns the prescribed air-KERMA value used for dose calibration.
    double airKerma() const { return m_airKerma; }

    /**
     * @brief Returns the dose calibration factor for this beam.
     *
     * Computes the ratio of the prescribed air KERMA to the expected air KERMA from
     * the simulated fluence:
     *
     *   factor = airKerma / (numberOfParticles · μ_tr/ρ(E) · E)
     *
     * where μ_tr/ρ is the NIST mass energy-transfer coefficient of dry air at the
     * beam energy. Returns 0 if the NIST air material cannot be loaded.
     *
     * @param progress  Unused; accepted for `BeamType` interface compatibility.
     * @return Calibration factor to multiply accumulated energy scores by.
     */
    double calibrationFactor(TransportProgress* progress = nullptr) const
    {
        auto air_cand = Material<5>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();
        const auto area = m_lenghtX * m_lenghtY;
        const double kerma = m_energy * numberOfParticles() * air.massEnergyTransferAttenuation(m_energy) / area; // keV/g/cm2
        return m_airKerma / kerma; //(mGy/kg)(g/keV/cm2)
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMFlatMono" padded with spaces, used by the
     *         serializer to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BEAMFlatMono";
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
     * Writes energy, weight, air KERMA, position, direction, exposure count, and
     * particles-per-exposure using the `Serializer` format.
     *
     * @return Byte buffer containing the complete beam state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_energy, buffer);
        Serializer::serialize(m_weight, buffer);
        Serializer::serialize(m_airKerma, buffer);
        Serializer::serialize(m_pos, buffer);
        Serializer::serialize(m_dirCosines[0], buffer);
        Serializer::serialize(m_dirCosines[1], buffer);
        Serializer::serialize(m_lenghtX, buffer);
        Serializer::serialize(m_lenghtY, buffer);
        Serializer::serialize(m_Nexposures, buffer);
        Serializer::serialize(m_particlesPerExposure, buffer);

        return buffer;
    }

    /**
     * @brief Reconstructs a `FlatMonoEnergyField` from a serialized byte buffer.
     *
     * Reads the fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<FlatMonoEnergyField>` containing the restored beam.
     */
    static std::optional<FlatMonoEnergyField<ENABLETRACKING>> deserialize1(std::span<const char> buffer)
    {
        FlatMonoEnergyField<ENABLETRACKING> item;
        buffer = Serializer::deserialize(item.m_energy, buffer);
        buffer = Serializer::deserialize(item.m_weight, buffer);
        buffer = Serializer::deserialize(item.m_airKerma, buffer);
        buffer = Serializer::deserialize(item.m_pos, buffer);
        buffer = Serializer::deserialize(item.m_dirCosines[0], buffer);
        buffer = Serializer::deserialize(item.m_dirCosines[1], buffer);
        buffer = Serializer::deserialize(item.m_lenghtX, buffer);
        buffer = Serializer::deserialize(item.m_lenghtY, buffer);
        buffer = Serializer::deserialize(item.m_Nexposures, buffer);
        buffer = Serializer::deserialize(item.m_particlesPerExposure, buffer);

        return std::make_optional(item);
    }

private:
    double m_energy = 60; ///< Photon energy [keV].
    double m_weight = 1; ///< Statistical weight of each particle.
    double m_airKerma = 1; ///< Prescribed air KERMA for dose calibration [mGy].
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction cosines {cos_x, cos_y}, normalised.
    double m_lenghtX = 1; ///< Field half-width along the x-axis [cm].
    double m_lenghtY = 1; ///< Field half-width along the y-axis [cm].
    std::uint64_t m_Nexposures = 100; ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100; ///< Particles emitted per exposure.
};
}

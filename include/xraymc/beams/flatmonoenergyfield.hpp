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

template <bool ENABLETRACKING = false>
class FlatMonoEnergyFieldExposure {
public:
    FlatMonoEnergyFieldExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& cosines, double lenghtX, double lenghtY, double energy, double weight, std::uint64_t N)
        : m_energy(energy)
        , m_weight(weight)
        , m_pos(pos)
        , m_NParticles(N)
    {
        m_dir = vectormath::cross(m_dirCosines);
        m_dirXY[0] = vectormath::scale(cosines[0], lenghtX);
        m_dirXY[1] = vectormath::scale(cosines[1], lenghtY);
    }

    const std::array<double, 3>& position() const { return m_pos; }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    auto sampleParticle(RandomState& state) const noexcept
    {

        const auto r1 = 2 * state.randomUniform() - 1;
        const auto r2 = 2 * state.randomUniform() - 1;

        const auto pos = vectormath::add(vectormath::scale(m_dirXY[0], r1), vectormath::scale(m_dirXY[1], r2));
        if constexpr (ENABLETRACKING) {
            ParticleTrack p = { .pos = pos,
                .dir = m_dir,
                .energy = m_energy,
                .weight = m_weight };
            p.registerPosition();
            return p;
        } else {
            Particle p = { .pos = pos,
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
    std::array<std::array<double, 3>, 2> m_dirXY = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction unormalised.
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Beam direction unit vector.
    std::uint64_t m_NParticles = 100; ///< Number of particles in this exposure.
};

template <bool ENABLETRACKING = false>
class FlatMonoEnergyField {
public:
    FlatMonoEnergyField(const std::array<double, 3>& pos = { 0, 0, 0 }, double energy = 60)
        : m_pos(pos)
        , m_energy(energy)
    {
    }

    void setEnergy(double energy)
    {
        m_energy = std::clamp(std::abs(energy), MIN_ENERGY(), MAX_ENERGY());
    }

    double energy() const { return m_energy; }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }

    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }

    void setPosition(const std::array<double, 3>& pos)
    {
        m_pos = pos;
    }

    const std::array<double, 3>& position() const
    {
        return m_pos;
    }

    const std::array<double, 3>& direction() const
    {
        return m_dir;
    }

    std::array<std::array<double, 3>, 2>& directionCosines() const
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
        m_dir = vectormath::cross(xdir, ydir);
        vectormath::normalize(m_dir);
    }

    void setLenghtX(double x)
    {
        m_lenghtX = std::abs(x);
    }

    void setLenghtY(double y)
    {
        m_lenghtY = std::abs(y);
    }

    void setLenght(double l)
    {
        setLenghtX(l);
        setLenghtY(l);
    }
    void setLenght(double x.double y)
    {
        setLenghtX(x);
        setLenghtY(y);
    }

    double lenghtX() const
    {
        return m_lenghtX;
    }
    double lenghtY() const
    {
        return m_lenghtY;
    }

    void setParticleWeight(double weight = 1)
    {
        m_weight = weight;
    }

    FlatMonoEnergyFieldExposure<ENABLETRACKING> exposure(std::size_t i) const noexcept
    {
        FlatMonoEnergyFieldExposure<ENABLETRACKING> exp(m_pos, m_dir, m_energy, m_weight, m_particlesPerExposure);
        return exp;
    }

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
        const double kerma = m_energy * numberOfParticles() * air.massEnergyTransferAttenuation(m_energy); // keV/g
        return m_airKerma / kerma; //(mGy/kg)(g/keV)
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BEAMPencilBeam" padded with spaces, used by the
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
     * @brief Reconstructs a `PencilBeam` from a serialized byte buffer.
     *
     * Reads the seven fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<PencilBeam>` containing the restored beam.
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
    double m_airKerma = 1; ///< Prescribed air KERMA for dose calibration.
    std::array<double, 3> m_pos = { 0, 0, 0 }; ///< Source position [cm].
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }; ///< In-plane direction cosines {cos_x, cos_y}, normalised.
    double lenghtX = 1;
    double lenghtY = 1;
    std::uint64_t m_Nexposures = 100; ///< Number of exposures.
    std::uint64_t m_particlesPerExposure = 100; ///< Particles emitted per exposure.
};
}
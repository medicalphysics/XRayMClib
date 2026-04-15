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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "xraymc/beams/tube/betheHeitlerCrossSection.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/material/atomhandler.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/serializer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <utility>
#include <vector>

namespace xraymc {

/**
 * @brief X-ray tube model for computing bremsstrahlung and characteristic spectra.
 *
 * Computes the photon energy spectrum emitted by a tungsten-anode X-ray tube using
 * the semi-relativistic Bethe-Heitler model (see `BetheHeitlerCrossSection`). The
 * spectrum is filtered by user-configured added filtration materials (Al, Cu, Sn, Ag
 * or any element by Z) and includes the tungsten K-edge characteristic lines.
 *
 * Derived quantities — half-value layer (HVL) in mm Al and mean spectrum energy — are
 * computed on demand and cached until any parameter changes.
 *
 * Satisfies the `SerializeItemType` concept via `magicID()`, `validMagicID()`,
 * `serialize()`, and `deserialize()`.
 *
 * Valid tube voltage range: 50–150 kV (see `minVoltage()` / `maxVoltage()`).
 */
class Tube {
public:
    /**
     * @brief Constructs a Tube with the given operating parameters.
     *
     * @param tubeVoltage      Peak tube voltage (kVp) [kV]. Clamped to [50, 150].
     * @param anodeAngleDeg    Anode take-off angle [degrees]. Controls self-filtration
     *                         and characteristic-line intensities. Default: 12°.
     * @param energyResolution Energy bin width used when sampling the spectrum [keV].
     *                         Clamped to [0.1, 10]. Default: 1.0 keV.
     */
    Tube(double tubeVoltage = 120.0, double anodeAngleDeg = 12.0, double energyResolution = 1.0)
        : m_voltage(tubeVoltage)
        , m_energyResolution(energyResolution)
    {
        setAnodeAngleDeg(anodeAngleDeg);
        m_hasCachedSpecterCalculations = false;
    }

    /// @brief Returns the maximum supported tube voltage [kV] (150 kV).
    static constexpr double maxVoltage() { return 150; }
    /// @brief Returns the minimum supported tube voltage [kV] (50 kV).
    static constexpr double minVoltage() { return 50; }

    /// @brief Returns the current tube voltage [kV].
    double voltage() const { return m_voltage; }

    /**
     * @brief Sets the tube voltage [kV].
     *
     * Invalidates any cached HVL and mean-energy values.
     *
     * @param voltage  Desired tube voltage [kV]; clamped to [minVoltage(), maxVoltage()].
     */
    void setVoltage(double voltage)
    {
        m_voltage = std::clamp(voltage, minVoltage(), maxVoltage());
        m_hasCachedSpecterCalculations = false;
    }

    /// @brief Returns the anode take-off angle [rad].
    double anodeAngle() const { return m_anodeAngle; }
    /// @brief Returns the anode take-off angle [degrees].
    double anodeAngleDeg() const { return m_anodeAngle * RAD_TO_DEG(); }

    /**
     * @brief Sets the anode take-off angle [rad].
     *
     * The angle is clamped to [0, π/2]. Invalidates cached spectrum values.
     *
     * @param angle  Anode angle in radians; the absolute value is used.
     */
    void setAnodeAngle(double angle)
    {
        auto a = std::abs(angle);
        if (a > PI_VAL() * 0.5)
            a = PI_VAL() * 0.5;
        m_anodeAngle = a;
        m_hasCachedSpecterCalculations = false;
    }

    /**
     * @brief Sets the anode take-off angle [degrees].
     * @param angle  Anode angle in degrees; converted to radians internally.
     */
    void setAnodeAngleDeg(double angle)
    {
        setAnodeAngle(angle * DEG_TO_RAD());
    }

    /**
     * @brief Adds or replaces an added-filtration material by atomic number.
     *
     * The material is identified by its atomic number Z and stored as a thickness
     * in mm. If Z is already present its thickness is overwritten. Invalidates
     * cached spectrum values.
     *
     * @param Z   Atomic number of the filtration element (must exist in `AtomHandler`).
     * @param mm  Filtration thickness [mm]; the absolute value is used.
     * @return True if Z is a known element and the filtration was set; false otherwise.
     */
    bool addFiltrationMaterial(const std::integral auto Z, const double mm)
    {
        if (AtomHandler::atomExists(Z)) {
            m_filtrationMaterials[static_cast<std::uint8_t>(Z)] = std::abs(mm);
            m_hasCachedSpecterCalculations = false;
            return true;
        }
        return false;
    }

    /// @brief Returns a read-only view of all added-filtration materials (Z → mm).
    const std::map<std::uint8_t, double>& filtrationMaterials() const { return m_filtrationMaterials; }

    /// @brief Convenience setter — adds aluminium (Z = 13) filtration of @p mm mm.
    void setAlFiltration(double mm)
    {
        addFiltrationMaterial(13, mm);
    }
    /// @brief Convenience setter — adds copper (Z = 29) filtration of @p mm mm.
    void setCuFiltration(double mm)
    {
        addFiltrationMaterial(29, mm);
    }
    /// @brief Convenience setter — adds tin (Z = 50) filtration of @p mm mm.
    void setSnFiltration(double mm)
    {
        addFiltrationMaterial(50, mm);
    }
    /// @brief Convenience setter — adds silver (Z = 47) filtration of @p mm mm.
    void setAgFiltration(double mm)
    {
        addFiltrationMaterial(47, mm);
    }

    /**
     * @brief Returns the added-filtration thickness for element Z [mm].
     *
     * @param Z  Atomic number of the filtration element.
     * @return Filtration thickness [mm], or 0 if Z has not been set.
     */
    double filtration(std::integral auto Z) const
    {
        const auto z8 = static_cast<std::uint8_t>(Z);
        if (m_filtrationMaterials.contains(z8))
            return m_filtrationMaterials.at(z8);
        return 0;
    }

    /// @brief Returns the aluminium (Z = 13) added filtration [mm].
    double AlFiltration() const
    {
        return filtration(13);
    }
    /// @brief Returns the copper (Z = 29) added filtration [mm].
    double CuFiltration() const
    {
        return filtration(29);
    }
    /// @brief Returns the tin (Z = 50) added filtration [mm].
    double SnFiltration() const
    {
        return filtration(50);
    }

    /**
     * @brief Removes all added-filtration materials.
     *
     * Invalidates cached HVL and mean-energy values.
     */
    void clearFiltrationMaterials()
    {
        m_filtrationMaterials.clear();
        m_hasCachedSpecterCalculations = false;
    }

    /**
     * @brief Sets the energy bin width used when computing the spectrum [keV].
     *
     * Smaller values give finer spectral resolution at the cost of computation time.
     * Invalidates cached HVL and mean-energy values.
     *
     * @param energyResolution  Bin width [keV]; clamped to [0.1, 10.0].
     */
    void setEnergyResolution(double energyResolution)
    {
        m_energyResolution = std::clamp(energyResolution, 0.1, 10.0);
        m_hasCachedSpecterCalculations = false;
    }

    /// @brief Returns the current energy bin width [keV].
    double energyResolution() const { return m_energyResolution; }

    /**
     * @brief Returns the energy bin centres used for spectrum sampling [keV].
     *
     * Generates energies from `energyResolution()` up to (and including) `voltage()`
     * in steps of `energyResolution()`.
     *
     * @return Vector of energy values [keV].
     */
    std::vector<double> getEnergy() const
    {
        std::vector<double> energies;
        auto hv = m_energyResolution;
        const auto n_elem = static_cast<std::size_t>(std::ceil(m_voltage / m_energyResolution));
        energies.reserve(n_elem);
        while (hv <= m_voltage) {
            energies.push_back(hv);
            hv = hv + m_energyResolution;
        }
        return energies;
    }

    /**
     * @brief Returns the full spectrum as a vector of (energy [keV], weight) pairs.
     *
     * Uses the internal energy grid from `getEnergy()` and delegates to the
     * vector overload of `getSpecter`.
     *
     * @param normalize  If true (default), weights sum to 1.0.
     * @return Vector of `{energy [keV], weight}` pairs.
     */
    std::vector<std::pair<double, double>> getSpecter(bool normalize = true) const
    {
        auto energies = getEnergy();

        auto specter = this->getSpecter(energies, normalize);

        std::vector<std::pair<double, double>> map;
        map.reserve(specter.size());
        for (std::size_t i = 0; i < specter.size(); ++i)
            map.push_back(std::make_pair(energies[i], specter[i]));
        return map;
    }

    /**
     * @brief Returns the spectrum weights for the given energy bins at an explicit anode angle.
     *
     * Computes bremsstrahlung intensities via the Bethe-Heitler model, adds K-edge
     * characteristic radiation, applies added-filtration attenuation, and optionally
     * normalises the result.
     *
     * @param energies    Energy bin centres [keV].
     * @param anodeAngle  Take-off angle from the anode surface [rad]; overrides `m_anodeAngle`.
     * @param normalize   If true (default), the returned weights sum to 1.0.
     * @return Vector of spectral weights, one per element of @p energies.
     */
    std::vector<double> getSpecter(const std::vector<double>& energies, const double anodeAngle, bool normalize = true) const
    {
        std::vector<double> specter(energies.size());
        const auto kVp = voltage();
        std::transform(std::execution::par_unseq, energies.begin(), energies.end(), specter.begin(), [kVp, anodeAngle](auto hv) -> double {
            return BetheHeitlerCrossSection::betheHeitlerSpectra(kVp, hv, anodeAngle);
        });

        // adding characteristic radiation
        addCharacteristicEnergy(energies, specter);
        filterSpecter(energies, specter);
        if (normalize) {
            normalizeSpecter(specter);
        }
        return specter;
    }
    /**
     * @brief Returns the spectrum weights for the given energy bins using the configured anode angle.
     *
     * Convenience overload that delegates to `getSpecter(energies, anodeAngle, normalize)`
     * with the stored `m_anodeAngle`.
     *
     * @param energies   Energy bin centres [keV].
     * @param normalize  If true (default), the returned weights sum to 1.0.
     * @return Vector of spectral weights, one per element of @p energies.
     */
    std::vector<double> getSpecter(const std::vector<double>& energies, bool normalize = true) const
    {
        return getSpecter(energies, m_anodeAngle, normalize);
    }

    /**
     * @brief Returns the first half-value layer in aluminium [mm].
     *
     * Computed from the air-KERMA-weighted spectrum using a boosted gradient-descent
     * search on the Beer-Lambert transmission through aluminium. Result is cached;
     * the non-const overload also caches the mean energy in the same pass.
     *
     * @return HVL in mm Al.
     */
    double mmAlHalfValueLayer()
    {
        if (!m_hasCachedSpecterCalculations) {
            const auto res = calculateSpecterValues();
            m_cachedHVL = res.HVL;
            m_cachedMeanEnergy = res.meanEnergy;
            m_hasCachedSpecterCalculations = true;
        }
        return m_cachedHVL;
    }

    /**
     * @brief Returns the first half-value layer in aluminium [mm] (const overload).
     *
     * Uses the cache when available, otherwise computes on-the-fly without caching.
     *
     * @return HVL in mm Al.
     */
    [[nodiscard]] double mmAlHalfValueLayer() const
    {
        if (!m_hasCachedSpecterCalculations) {
            return calculateSpecterValues().HVL;
        }
        return m_cachedHVL;
    }

    /**
     * @brief Returns the fluence-weighted mean photon energy of the spectrum [keV].
     *
     * Computed as Σ(E_i · w_i) over the normalised spectrum. Result is cached;
     * the non-const overload also caches the HVL in the same pass.
     *
     * @return Mean photon energy [keV].
     */
    double meanSpecterEnergy()
    {
        if (!m_hasCachedSpecterCalculations) {
            const auto res = calculateSpecterValues();
            m_cachedHVL = res.HVL;
            m_cachedMeanEnergy = res.meanEnergy;
            m_hasCachedSpecterCalculations = true;
        }
        return m_cachedMeanEnergy;
    }

    /**
     * @brief Returns the fluence-weighted mean photon energy [keV] (const overload).
     *
     * Uses the cache when available, otherwise computes on-the-fly without caching.
     *
     * @return Mean photon energy [keV].
     */
    [[nodiscard]] double meanSpecterEnergy() const
    {
        if (!m_hasCachedSpecterCalculations) {
            return calculateSpecterValues().meanEnergy;
        }
        return m_cachedMeanEnergy;
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "XrayTube" padded with spaces, used by the serializer
     *         to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "XrayTube";
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
     * @brief Serializes the tube configuration to a byte buffer.
     *
     * Writes voltage, energy resolution, anode angle, and filtration materials
     * using the `Serializer` format.
     *
     * @return Byte buffer containing the complete tube state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_voltage, buffer);
        Serializer::serialize(m_energyResolution, buffer);
        Serializer::serialize(m_anodeAngle, buffer);
        Serializer::serializeMaterialWeights(m_filtrationMaterials, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `Tube` from a serialized byte buffer.
     *
     * Reads the four fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<Tube>` containing the restored tube configuration.
     */
    static std::optional<Tube> deserialize(std::span<const char> buffer)
    {
        Tube item;

        buffer = Serializer::deserialize(item.m_voltage, buffer);
        buffer = Serializer::deserialize(item.m_energyResolution, buffer);
        buffer = Serializer::deserialize(item.m_anodeAngle, buffer);
        buffer = Serializer::deserializeMaterialWeights(item.m_filtrationMaterials, buffer);
        return item;
    }

protected:
    /**
     * @brief POD holding the two lazily computed spectrum summary statistics.
     */
    struct SpecterCalculatedValues {
        double HVL;        ///< First half-value layer in aluminium [mm].
        double meanEnergy; ///< Fluence-weighted mean photon energy [keV].
    };

    /**
     * @brief Adds tungsten K-edge characteristic line intensities to a spectrum.
     *
     * For each of the four dominant K-edge lines, finds the nearest energy bin
     * (within 2 keV) in @p energy and increments the corresponding bin in @p specter
     * by the line intensity returned by `BetheHeitlerCrossSection::characteristicTungstenKedge`.
     *
     * @param energy   Energy bin centres [keV].
     * @param specter  Spectral weights to modify in place.
     */
    void addCharacteristicEnergy(const std::vector<double>& energy, std::vector<double>& specter) const
    {
        auto energyBegin = energy.begin();
        auto energyEnd = energy.end();

        const auto voltage_d = voltage();
        const auto anode_d = anodeAngle();

        const auto kEdge = BetheHeitlerCrossSection::characteristicTungstenKedge(voltage_d, anode_d);
        for (const auto& [e, n] : kEdge) {
            // find closest energy
            auto eIdx = std::lower_bound(energyBegin, energyEnd, e);
            if (eIdx != energyEnd) {

                if (std::abs(e - *eIdx) <= 2.0) {
                    // we only add characteristic radiation if specter energy is closer than 2 keV from the K edge
                    auto nIdx = specter.begin();
                    std::advance(nIdx, std::distance(energyBegin, eIdx));
                    *nIdx = *nIdx + n; // adding characteristic k edge intensity to specter
                }
            }
        }
    }

    /**
     * @brief Applies Beer-Lambert attenuation from all added-filtration materials.
     *
     * Accumulates the total linear attenuation (photoelectric + Compton + Rayleigh)
     * for each filtration element at each energy, then multiplies the spectral weights
     * by exp(−∑ μ_i · ρ_i · t_i).
     *
     * @param energies  Energy bin centres [keV].
     * @param specter   Spectral weights to attenuate in place.
     */
    void filterSpecter(const std::vector<double>& energies, std::vector<double>& specter) const
    {
        std::vector<double> totAtt(energies.size(), 0.0);
        for (auto const& [Z, mm] : m_filtrationMaterials) {
            const auto& atom = AtomHandler::Atom(Z);
            const auto cm = mm * 0.1; // for mm -> cm
            const auto dist = atom.standardDensity * cm;
            auto p = interpolate(atom.photoel, energies);
            auto in = interpolate(atom.incoherent, energies);
            auto co = interpolate(atom.coherent, energies);
            for (std::size_t i = 0; i < energies.size(); ++i) {
                totAtt[i] += (p[i] + in[i] + co[i]) * dist;
            }
        }
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), totAtt.cbegin(), specter.begin(),
            [&](const auto s, const auto att) { return s * std::exp(-att); });
    }

    /**
     * @brief Normalises spectral weights so they sum to 1.
     * @param specter  Spectral weights to normalise in place.
     */
    void normalizeSpecter(std::vector<double>& specter) const
    {
        const auto sum = std::reduce(std::execution::par_unseq, specter.cbegin(), specter.cend());
        std::for_each(std::execution::par_unseq, specter.begin(), specter.end(), [sum](auto& n) { n = n / sum; });
    }

    /**
     * @brief Computes and returns both HVL and mean energy from the current spectrum.
     *
     * Generates the normalised spectrum, computes the fluence-weighted mean energy,
     * then finds the HVL in mm Al using a boosted gradient-descent search on the
     * air-KERMA transmission through aluminium slabs.
     *
     * @return `SpecterCalculatedValues` containing HVL [mm Al] and mean energy [keV].
     */
    [[nodiscard]] SpecterCalculatedValues calculateSpecterValues() const
    {
        const auto energy = getEnergy();
        const auto specter = getSpecter(energy, true);

        SpecterCalculatedValues res;

        // Mean energy
        res.meanEnergy = std::transform_reduce(std::execution::par_unseq, energy.cbegin(), energy.cend(), specter.cbegin(), double { 0 }, std::plus { }, std::multiplies { });

        // HVL
        // Al att vector
        const auto& Al = AtomHandler::Atom(13);
        const auto photo = interpolate(Al.photoel, energy);
        const auto incoherent = interpolate(Al.incoherent, energy);
        const auto coherent = interpolate(Al.coherent, energy);
        auto att = addVectors(photo, incoherent, coherent);
        const auto dens = Al.standardDensity;
        std::for_each(std::execution::par_unseq, att.begin(), att.end(), [dens](auto& a) { a *= dens; });

        // Calculating air KERMA
        auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
        std::vector<double> kerma(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), specter.cbegin(), kerma.begin(), [&air](auto e, auto s) {
            return e * s * air.massEnergyTransferAttenuation(e);
        });
        const auto norm = std::reduce(std::execution::par_unseq, kerma.cbegin(), kerma.cend(), 0.0);
        std::for_each(std::execution::par_unseq, kerma.begin(), kerma.end(), [norm](auto& k) { k /= norm; });

        // Boosted gradient decent for finding half value layer
        auto x = 0.5;
        double step = 0;
        int iter = 0;
        double g;
        do {
            x = std::max(x + step, 0.0001);
            g = std::transform_reduce(std::execution::par_unseq, kerma.cbegin(), kerma.cend(), att.cbegin(), 0.0, std::plus<>(), [x](auto s, auto a) -> double { return s * std::exp(-a * x); });
            step = (g - 0.5) * std::max(5 - iter, 1);
        } while ((std::abs(g - 0.5) > 0.005 && iter++ < 10));

        res.HVL = x * 10; // cm -> mm
        return res;
    }

private:
    double m_voltage = 120;                          ///< Tube voltage [kV].
    double m_energyResolution = 1;                   ///< Energy bin width [keV].
    double m_anodeAngle = 0.21;                      ///< Anode take-off angle [rad] (~12°).
    double m_cachedHVL = 0;                          ///< Cached first HVL in aluminium [mm].
    double m_cachedMeanEnergy = 0;                   ///< Cached fluence-weighted mean energy [keV].
    std::map<std::uint8_t, double> m_filtrationMaterials; ///< Added filtration: Z → thickness [mm].
    bool m_hasCachedSpecterCalculations = false;     ///< True when cached HVL and mean energy are valid.
};
}
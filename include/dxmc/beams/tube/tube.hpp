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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/beams/tube/betheHeitlerCrossSection.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/material.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <utility>
#include <vector>

namespace dxmc {

class Tube {
public:
    Tube(double tubeVoltage = 120.0, double anodeAngleDeg = 12.0, double energyResolution = 1.0)
        : m_voltage(tubeVoltage)
        , m_energyResolution(energyResolution)
    {
        setAnodeAngleDeg(anodeAngleDeg);
        m_hasCachedSpecterCalculations = false;
    }

    static constexpr double maxVoltage() { return 150; }
    static constexpr double minVoltage() { return 50; }

    double voltage() const { return m_voltage; }
    void setVoltage(double voltage)
    {
        m_voltage = std::clamp(voltage, minVoltage(), maxVoltage());
        m_hasCachedSpecterCalculations = false;
    }

    double anodeAngle() const { return m_anodeAngle; }
    double anodeAngleDeg() const { return m_anodeAngle * RAD_TO_DEG(); }
    void setAnodeAngle(double angle)
    {
        auto a = std::abs(angle);
        if (a > PI_VAL() * 0.5)
            a = PI_VAL() * 0.5;
        m_anodeAngle = a;
        m_hasCachedSpecterCalculations = false;
    }
    void setAnodeAngleDeg(double angle)
    {
        setAnodeAngle(angle * DEG_TO_RAD());
    }

    bool addFiltrationMaterial(const std::size_t Z, const double mm)
    {
        if (AtomHandler::atomExists(Z)) {
            m_filtrationMaterials[Z] = std::abs(mm);
            m_hasCachedSpecterCalculations = false;
            return true;
        }
        return false;
    }

    const std::map<std::size_t, double>& filtrationMaterials() const { return m_filtrationMaterials; }

    void setAlFiltration(double mm)
    {
        addFiltrationMaterial(13, mm);
    }
    void setCuFiltration(double mm)
    {
        addFiltrationMaterial(29, mm);
    }
    void setSnFiltration(double mm)
    {
        addFiltrationMaterial(50, mm);
    }
    void setAgFiltrationMaterial(double mm)
    {
        addFiltrationMaterial(47, mm);
    }

    double filtration(std::size_t Z) const
    {
        if (m_filtrationMaterials.contains(Z))
            return m_filtrationMaterials.at(Z);
        return 0;
    }

    double AlFiltration() const
    {
        return filtration(13);
    }
    double CuFiltration() const
    {
        return filtration(29);
    }
    double SnFiltration() const
    {
        return filtration(50);
    }

    void clearFiltrationMaterials()
    {
        m_filtrationMaterials.clear();
        m_hasCachedSpecterCalculations = false;
    }

    void setEnergyResolution(double energyResolution)
    {
        m_energyResolution = std::clamp(energyResolution, 0.1, 10.0);
        m_hasCachedSpecterCalculations = false;
    }

    double energyResolution() const { return m_energyResolution; }

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
    std::vector<double> getSpecter(const std::vector<double>& energies, bool normalize = true) const
    {
        return getSpecter(energies, m_anodeAngle, normalize);
    }

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

    [[nodiscard]] double mmAlHalfValueLayer() const
    {
        if (!m_hasCachedSpecterCalculations) {
            return calculateSpecterValues().HVL;
        }
        return m_cachedHVL;
    }

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

    [[nodiscard]] double meanSpecterEnergy() const
    {
        if (!m_hasCachedSpecterCalculations) {
            return calculateSpecterValues().meanEnergy;
        }
        return m_cachedMeanEnergy;
    }

protected:
    struct SpecterCalculatedValues {
        double HVL, meanEnergy;
    };
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

    void normalizeSpecter(std::vector<double>& specter) const
    {
        const auto sum = std::reduce(std::execution::par_unseq, specter.cbegin(), specter.cend());
        std::for_each(std::execution::par_unseq, specter.begin(), specter.end(), [sum](auto& n) { n = n / sum; });
    }

    [[nodiscard]] SpecterCalculatedValues calculateSpecterValues() const
    {
        const auto energy = getEnergy();
        const auto specter = getSpecter(energy, true);

        SpecterCalculatedValues res;

        // Mean energy
        res.meanEnergy = std::transform_reduce(std::execution::par_unseq, energy.cbegin(), energy.cend(), specter.cbegin(), double { 0 }, std::plus {}, std::multiplies {});

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
    double m_voltage = 120;
    double m_energyResolution = 1;
    double m_anodeAngle = 0.21; // about 12 degrees
    double m_cachedHVL = 0;
    double m_cachedMeanEnergy = 0;
    std::map<std::size_t, double> m_filtrationMaterials;
    bool m_hasCachedSpecterCalculations = false;
};
}
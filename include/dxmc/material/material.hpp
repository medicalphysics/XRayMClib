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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/nistmaterials.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <execution>

namespace dxmc {

struct MaterialShell {
    double numberOfElectronsFraction = 0;
    double bindingEnergy = 0;
    double HartreeFockOrbital_0 = 0;
    double numberOfPhotonsPerInitVacancy = 0;
    double energyOfPhotonsPerInitVacancy = 0;
};

struct AttenuationValues {
    double photoelectric;
    double incoherent;
    double coherent;
    double sum() const noexcept { return photoelectric + incoherent + coherent; }
};

template <std::size_t N = 5>
class Material {
public:
    static std::optional<Material<N>> byZ(std::size_t Z)
    {
        auto a = AtomHandler::Atom(Z);
        if (a.Z == Z) {
            std::map<std::size_t, double> w;
            w[Z] = 1;
            return byWeight(w);
        }
        return std::nullopt;
    }
    static std::optional<Material<N>> byWeight(const std::map<std::size_t, double>& weights)
    {
        auto m = constructMaterial(weights);
        return m;
    }

    static std::optional<Material<N>> byChemicalFormula(const std::string& str)
    {
        auto numberDensCompound = parseCompoundStr(str);
        if (numberDensCompound.size() == 0)
            return std::nullopt;
        std::map<std::size_t, double> weight;
        for (const auto [Z, numDens] : numberDensCompound) {
            if (Z == 0 || Z > 99)
                return std::nullopt;
            const auto& atom = AtomHandler::Atom(Z);
            weight[Z] = numDens * atom.atomicWeight;
        }
        return Material<N>::byWeight(weight);
    }

    static std::optional<Material<N>> byNistName(const std::string& name)
    {
        const auto& w = NISTMaterials::Composition(name);
        if (!w.empty())
            return Material<N>::byWeight(w);

        return std::nullopt;
    }

    static std::vector<std::string> listNistCompoundNames()
    {
        return NISTMaterials::listNames();
    }

    inline double effectiveZ() const { return m_effectiveZ; }

    inline double scatterFactor(double momentumTransfer) const
    {

        return 0.0;
    }

    inline double formFactor(const double momentumTransfer) const
    {
        return 0.0;
    }

    inline double sampleSquaredMomentumTransferFromFormFactorSquared(double qsquared_max, RandomState& state) const
    {
        return m_formFactorInvSamp(qsquared_max, state);
    }

    inline double meanIncoherentScatterEnergy(double energy) const
    {

        return 0.0;
    }

    inline double massEnergyTransferAttenuation(double energy) const
    {
        const auto att = attenuationValues(energy);
        return massEnergyTransferAttenuation(att, energy);
    }

    inline double massEnergyTransferAttenuation(const AttenuationValues& att, double energy) const
    {
        double avg_fluro_energy = 0;
        for (std::uint_fast8_t i = 0; i < m_numberOfShells; ++i) {
            if (m_shells[i].bindingEnergy > MIN_ENERGY())
                avg_fluro_energy += m_shells[i].energyOfPhotonsPerInitVacancy * attenuationPhotoelectricShell(i, energy);
        }
        avg_fluro_energy /= att.photoelectric;

        const auto avg_fluro_frac = avg_fluro_energy / energy;
        const auto incoherent_scatter_energy = meanIncoherentScatterEnergy(energy);
        const auto f_pe = 1 - avg_fluro_frac;
        const auto f_inc = 1 - (incoherent_scatter_energy + avg_fluro_energy) / energy;
        return att.photoelectric * f_pe + att.incoherent * f_inc;
    }

    // photo, coherent, incoherent
    inline AttenuationValues attenuationValues(double energy) const
    {

        AttenuationValues att {
            .photoelectric = 0,
            .incoherent = 0,
            .coherent = 0
        };
        return att;
    }
    // photo, coherent, incoherent
    std::vector<AttenuationValues> attenuationValues(const std::vector<double>& energy) const
    {
        std::vector<AttenuationValues> att(energy.size());
        return att;
    }
    inline std::vector<double> totalAttenuationValue(const std::vector<double>& energy) const
    {
        const auto att = attenuationValues(energy);
        std::vector<double> sum(att.size());
        std::transform(std::execution::par_unseq, att.cbegin(), att.cend(), sum.begin(), [](const auto& a) { return a.sum(); });
        return sum;
    }

    inline double attenuationPhotoeletric(double energy) const
    {
        return 0.0;
    }

    inline double attenuationPhotoelectricShell(std::uint8_t shell, double energy) const
    {
        return 0.0;
    }

    inline std::uint8_t numberOfShells() const
    {
        return m_numberOfShells;
    }

    inline const MaterialShell& shell(std::size_t shell) const
    {
        return m_shells[shell];
    }

    inline const auto& shells() const { return m_shells; }

    static double momentumTransfer(double energy, double angle)
    {
        return momentumTransferMax(energy) * std::sin(angle * 0.5); // per Å
    }

    static double momentumTransferCosAngle(double energy, double cosAngle)
    {
        return momentumTransferMax(energy) * std::sqrt((1 - cosAngle) * 0.5); // per Å
    }

    static constexpr double momentumTransferMax(double energy)
    {
        static constexpr double hc_si = 1.239841193E-6; // ev*m
        static constexpr double m2A = 1E10; // meters to ångstrøm
        static constexpr double eV2keV = 1E-3; // eV to keV
        static constexpr double hc = hc_si * m2A * eV2keV; // kev*Å
        static constexpr double hc_inv = 1 / hc;
        return energy * hc_inv; // per Å
    }

    static bool validCompoundString(const std::string& cmp)
    {
        const auto m = parseCompoundStr(cmp);
        bool valid = true;
        for (const auto& [Z, w] : m) {
            valid = valid && Z > 0 && Z <= 99 && w > 0;
        }
        return valid;
    }

    static bool validNistName(const std::string& name)
    {
        const auto& w = NISTMaterials::Composition(name);
        return !w.empty();
    }

    /**
     * This function parses a chemical formula string and returns a map of elements
     * Z and number density (not normalized). It's kinda messy but supports parenthesis
     * in the expression.
     */
    static std::map<std::uint64_t, double> parseCompoundStr(const std::string& str)
    {
        static std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };

        std::map<std::uint64_t, double> parsed;

        auto begin = str.begin();
        std::vector<std::pair<std::string, std::string>> vals;
        while (begin != str.end()) {
            const char l = *begin;
            if (std::isupper(l)) {
                auto v = std::pair<std::string, std::string>();
                v.first = std::string(&l, 1);
                vals.push_back(v);
            } else if (std::islower(l)) {
                if (vals.size() > 0)
                    vals.back().first.append(&l, 1);
            } else if ((l == '.') || std::isdigit(l)) {
                if (vals.size() > 0)
                    vals.back().second.append(&l, 1);
            } else if (l == '(') {
                auto startP = begin + 1;
                auto endP = startP;
                while (endP != str.end()) {
                    const char l1 = *endP;
                    if (l1 == ')') {
                        std::string strP(startP, endP);
                        auto parseP = Material<N>::parseCompoundStr(strP);
                        auto endPdig = endP + 1;
                        double w = 1;
                        std::string str_w;
                        while (endPdig != str.end()) {
                            const char l2 = *endPdig;
                            if (std::isdigit(l2) || l2 == '.') {
                                str_w.append(&l2, 1);
                                begin = endPdig;
                                endPdig++;
                            } else {
                                endPdig = str.end();
                            }
                        }
                        try {
                            w = std::stod(str_w);
                        } catch (const std::invalid_argument& e) {
                        }
                        for (const auto& [el, frac] : parseP) {
                            if (!parsed.contains(el)) {
                                parsed[el] = 0;
                            }
                            parsed[el] += (frac * w);
                        }
                        endP = str.end();
                    } else {
                        endP++;
                    }
                }
            }
            begin++;
        }
        for (const auto& [el, num] : vals) {
            auto pos = std::find(S.begin(), S.end(), el);
            if (pos != S.end()) {
                std::uint64_t Z = std::distance(S.begin(), pos) + 1;
                double w = 1.0;
                try {
                    w = std::stod(num);
                } catch (const std::invalid_argument& e) {
                }
                if (!parsed.contains(Z)) {
                    parsed[Z] = 0;
                }
                parsed[Z] += w;
            }
        }
        return parsed;
    }

protected:
    Material()
    {
    }

    enum class LUTType {
        photoelectric,
        coherent,
        incoherent,
        formfactor,
        scatterfactor,
        incoherentenergy
    };

    static std::optional<Material<N>> constructMaterial(const std::map<std::size_t, double>& compositionByWeight)
    {
        for (const auto& [Z, w] : compositionByWeight) {
            const auto& a = AtomHandler::Atom(Z);
            if (a.Z != Z)
                return std::nullopt;
        }

        auto weight = compositionByWeight;
        const auto totalWeight = std::transform_reduce(weight.cbegin(), weight.cend(), 0.0, std::plus<>(), [](const auto& right) -> double { return right.second; });
        std::for_each(weight.begin(), weight.end(), [=](auto& w) { w.second /= totalWeight; });

        Material<N> m;
        m.m_effectiveZ = std::transform_reduce(weight.cbegin(), weight.cend(), 0.0, std::plus<>(), [](const auto& pair) -> double { return pair.first * pair.second; });

        std::vector<double> photoel, coherent, incoherent;

        for (const auto& [Z, w] : weight) {
            const auto& a = AtomHandler::Atom(Z);
            for (const auto& p : a.photoel)
                photoel.push_back(p.first);
            for (const auto& p : a.coherent)
                coherent.push_back(p.first);
            for (const auto& p : a.incoherent)
                incoherent.push_back(p.first);
        }
        std::sort(photoel.begin(), photoel.end());
        std::sort(coherent.begin(), coherent.end());
        std::sort(incoherent.begin(), incoherent.end());

        auto AddScale = [](std::vector<double>& res, const std::vector<double>& adder, double scale) {
            std::transform(std::execution::par_unseq, adder.cbegin() adder.cend(), res.cbegin(), res.begin(), [scale](double a, double res) { return res + a * scale; });
        };
        std::vector<double> photoel_val(photoel.size(), 0.0);
        std::vector<double> coherent_val(coherent.size(), 0.0);
        std::vector<double> incoherent_val(incoherent.size(), 0.0);
        for (const auto& [Z, w] : weight) {
            const auto& a = AtomHandler::Atom(Z);
            const auto p = interpolate(a.photoel, photoel);
            AddScale(photoel_val, p, w);
            const auto in = interpolate(a.incoherent, incoherent);
            AddScale(incoherent_val,in, w);
            const auto co = interpolate(a.coherent, coherent);
            AddScale(coherent_val, co, w);
        }
        


        // createMaterialAtomicShells(m, weight, offset);

        generateFormFactorInverseSampling(m);

        return m;
    }

    static void generateFormFactorInverseSampling(Material<N>& material)
    {
        const auto qmax = Material<N>::momentumTransferMax(MAX_ENERGY());
        constexpr double qmin = 0.001;
        auto func = [&material](double qsquared) -> double {
            const auto q = std::sqrt(qsquared);
            const auto f = material.formFactor(q);
            return f * f;
        };

        material.m_formFactorInvSamp = CPDFSampling<double, 20>(qmin * qmin, qmax * qmax, func);
    }

    static void createMaterialAtomicShells(Material<N>& material, const std::map<std::size_t, double>& normalizedWeight)
    {
    }

private:
    double m_effectiveZ = 0;
    std::vector<std::pair<double, double>> m_photoel;
    std::vector<std::pair<double, double>> m_coherent;
    std::vector<std::pair<double, double>> m_incoherent;
    CPDFSampling<double, 20> m_formFactorInvSamp;
    std::array<MaterialShell, N + 1> m_shells;
    std::uint8_t m_numberOfShells = 0;
};
}
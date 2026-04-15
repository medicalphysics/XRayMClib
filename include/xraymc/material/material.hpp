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
#include "xraymc/interpolation.hpp"
#include "xraymc/material/atomhandler.hpp"
#include "xraymc/material/atomicelement.hpp"
#include "xraymc/material/nistmaterials.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <concepts>
#include <execution>

namespace xraymc {

/**
 * @brief Per-subshell data mixed down from the constituent elements of a `Material`.
 *
 * Stores the weighted-average orbital parameters needed to sample photoelectric
 * interactions and fluorescence emission within a compound material. At most `N`
 * explicit shells are tracked; any remaining shells are collapsed into a single
 * remainder shell by `Material::createMaterialAtomicShells`.
 */
struct MaterialShell {
    double numberOfElectronsFraction = 0;       ///< Fraction of total material electrons in this shell (normalised to 1 across all shells).
    double bindingEnergy = 0;                    ///< Shell binding energy [keV].
    double HartreeFockOrbital_0 = 0;             ///< Hartree–Fock orbital momentum parameter p₀ [a.u.] for Compton profile sampling.
    double numberOfPhotonsPerInitVacancy = 0;    ///< Average fluorescence photons emitted per initial vacancy in this shell.
    double energyOfPhotonsPerInitVacancy = 0;    ///< Average total fluorescence energy per initial vacancy [keV].

    /// @brief Equality comparison (all fields).
    bool operator==(const MaterialShell& other) const = default;
};

/**
 * @brief Mass attenuation coefficients [cm²/g] for the three photon interaction types at a given energy.
 */
struct AttenuationValues {
    double photoelectric;   ///< Photoelectric absorption mass attenuation coefficient [cm²/g].
    double incoherent;      ///< Incoherent (Compton) scattering mass attenuation coefficient [cm²/g].
    double coherent;        ///< Coherent (Rayleigh) scattering mass attenuation coefficient [cm²/g].
    /// @brief Returns the total mass attenuation coefficient (sum of all three interactions) [cm²/g].
    double sum() const noexcept { return photoelectric + incoherent + coherent; }
};

/**
 * @brief A photon-transport material built from weighted elemental cross-section data.
 *
 * Stores cubic least-squares spline interpolation tables (in log-log space) for
 * photoelectric absorption, incoherent (Compton) scattering, coherent (Rayleigh)
 * scattering, form factors, incoherent scattering functions, and mean Compton scatter
 * energies. Per-subshell photoelectric cross-sections are kept separately for up to `N`
 * shells; excess shells are collapsed into a single remainder shell.
 *
 * A `Material` cannot be default-constructed publicly; use one of the static factory
 * methods:
 * - `byZ(Z)`                    — pure element by atomic number.
 * - `byWeight(weights)`         — mixture by elemental mass fractions.
 * - `byChemicalFormula(str)`    — compound from a chemical formula string (e.g. "H2O").
 * - `byNistName(name)`          — compound from the NIST material registry.
 *
 * @tparam N  Maximum number of explicit atomic subshells tracked per material.
 *            Shells beyond `N` are averaged into a remainder shell. Default: 16.
 */
template <std::size_t N = 16>
class Material {
public:
    /// @brief Equality comparison across all internal tables and composition data.
    bool operator==(const Material<N>& other) const = default;

    /**
     * @brief Constructs a pure-element material for atomic number @p Z.
     *
     * @param Z  Atomic number (integral type).
     * @return An engaged `optional<Material>` if @p Z is in the loaded physics data,
     *         `nullopt` otherwise.
     */
    static std::optional<Material<N>> byZ(std::integral auto Z)
    {
        auto a = AtomHandler::Atom(Z);
        if (a.Z == static_cast<std::uint64_t>(Z)) {
            std::map<std::uint8_t, double> w;
            w[static_cast<std::uint8_t>(Z)] = 1;
            return byWeight(w);
        }
        return std::nullopt;
    }

    /**
     * @brief Constructs a compound or mixture material from elemental mass fractions.
     *
     * Weights are normalised internally so they do not need to sum to 1. Entries
     * whose Z is not found in the loaded physics data cause `nullopt` to be returned.
     *
     * @tparam I  Integral key type for the weight map.
     * @param weights  Map of `{Z → mass fraction}` (unnormalised).
     * @return An engaged `optional<Material>` on success, `nullopt` if any Z is invalid.
     */
    template <std::integral I>
    static std::optional<Material<N>> byWeight(const std::map<I, double>& weights)
    {
        if constexpr (std::is_same<I, std::uint8_t>::value) {
            auto m = constructMaterial(weights);
            return m;
        } else {
            std::map<std::uint8_t, double> wt;
            for (const auto& [Z, w] : weights) {
                wt[static_cast<std::uint8_t>(Z)] = w;
            }
            auto m = constructMaterial(wt);
            return m;
        }
    }

    /**
     * @brief Constructs a material from a chemical formula string (e.g. "H2O", "Ca3(PO4)2").
     *
     * Parses element symbols and stoichiometric counts (including parenthesised groups)
     * via `parseCompoundStr`, converts atom counts to mass fractions using tabulated
     * atomic weights, then delegates to `byWeight`. Returns `nullopt` if the string
     * cannot be parsed or contains unsupported elements (Z outside [1, 99]).
     *
     * @param str  Chemical formula string.
     * @return An engaged `optional<Material>` on success, `nullopt` on parse failure.
     */
    static std::optional<Material<N>> byChemicalFormula(const std::string& str)
    {
        auto numberDensCompound = parseCompoundStr(str);
        if (numberDensCompound.size() == 0)
            return std::nullopt;
        std::map<std::uint8_t, double> weight;
        for (const auto [Z, numDens] : numberDensCompound) {
            if (Z == 0 || Z > 99)
                return std::nullopt;
            const auto& atom = AtomHandler::Atom(Z);
            weight[Z] = numDens * atom.atomicWeight;
        }
        return Material<N>::byWeight(weight);
    }

    /**
     * @brief Constructs a material from a NIST compound name.
     *
     * Looks up the mass fractions in `NISTMaterials` and delegates to `byWeight`.
     * Returns `nullopt` if @p name is not found in the NIST registry.
     *
     * @param name  NIST compound name as returned by `listNistCompoundNames()`.
     * @return An engaged `optional<Material>` on success, `nullopt` if not found.
     */
    static std::optional<Material<N>> byNistName(const std::string& name)
    {
        const auto& w = NISTMaterials::Composition(name);
        if (!w.empty())
            return Material<N>::byWeight(w);

        return std::nullopt;
    }

    /**
     * @brief Returns the names of all materials available in the NIST registry.
     * @return Vector of NIST material name strings.
     */
    static std::vector<std::string> listNistCompoundNames()
    {
        return NISTMaterials::listNames();
    }

    /**
     * @brief Returns the effective atomic number Z_eff of the material.
     *
     * Computed as the mass-fraction-weighted mean of the constituent atomic numbers:
     * Z_eff = Σ Z_i × w_i.
     *
     * @return Effective atomic number (dimensionless).
     */
    double effectiveZ() const
    {
        double Zeff = std::transform_reduce(m_composition.cbegin(), m_composition.cend(), 0.0, std::plus<>(), [](const auto& pair) -> double { return pair.first * pair.second; });
        return Zeff;
    }

    /**
     * @brief Returns the normalised elemental composition as mass fractions.
     * @return Const reference to the `{Z → mass fraction}` map (fractions sum to 1).
     */
    const std::map<std::uint8_t, double>& composition() const
    {
        return m_composition;
    }

    /**
     * @brief Returns the incoherent scattering function S(x) at the given momentum transfer.
     *
     * Evaluates the cubic spline fit to the mass-fraction-weighted incoherent
     * scattering function in log-log space.
     *
     * @param momentumTransfer  Momentum transfer x [Å⁻¹].
     * @return S(x) — incoherent scattering function (dimensionless, 0–Z).
     */
    inline double scatterFactor(double momentumTransfer) const
    {
        const auto logmomt = std::log(momentumTransfer);
        const auto begin = m_attenuationTable.cbegin() + m_attenuationTableOffset[4];
        const auto end = m_attenuationTable.cbegin() + m_attenuationTableOffset[5];
        return std::exp(CubicLSInterpolator<double>::evaluateSpline(logmomt, begin, end));
    }

    /**
     * @brief Returns the coherent scattering form factor F(x) at the given momentum transfer.
     *
     * Evaluates the cubic spline fit to the mass-fraction-weighted form factor in
     * log-log space.
     *
     * @param momentumTransfer  Momentum transfer x [Å⁻¹].
     * @return F(x) — form factor (dimensionless, 0–Z).
     */
    inline double formFactor(const double momentumTransfer) const
    {
        const auto logmomt = std::log(momentumTransfer);
        const auto begin = m_attenuationTable.cbegin() + m_attenuationTableOffset[3];
        const auto end = m_attenuationTable.cbegin() + m_attenuationTableOffset[4];
        return std::exp(CubicLSInterpolator<double>::evaluateSpline(logmomt, begin, end));
    }

    /**
     * @brief Samples a squared momentum transfer q² from the form-factor-squared distribution.
     *
     * Uses the precomputed CPDF inverse-sampling table to draw q² ∈ [q²_min, @p qsquared_max]
     * proportional to F²(q) — the squared form factor — for coherent (Rayleigh) scatter angle
     * sampling.
     *
     * @param qsquared_max  Upper bound for q² [Å⁻²], typically `momentumTransferMax(energy)²`.
     * @param state         Per-thread PRNG state.
     * @return Sampled q² [Å⁻²].
     */
    inline double sampleSquaredMomentumTransferFromFormFactorSquared(double qsquared_max, RandomState& state) const
    {
        return m_formFactorInvSamp(qsquared_max, state);
    }

    /**
     * @brief Returns the mean scattered photon energy for Compton interactions [keV].
     *
     * @param energy  Incident photon energy [keV].
     * @return Mean energy of the scattered photon [keV].
     */
    inline double meanIncoherentScatterEnergy(double energy) const
    {
        const auto begin = m_attenuationTable.cbegin() + m_attenuationTableOffset[5];
        const auto end = m_attenuationTable.cbegin() + m_attenuationTableOffset[6];
        return std::exp(CubicLSInterpolator<double>::evaluateSpline(std::log(energy), begin, end));
    }

    /**
     * @brief Returns the mean scattered photon energy for Compton interactions [keV],
     *        with the incident energy supplied as its natural logarithm.
     *
     * Avoids a redundant `log()` call when the log-energy is already available.
     *
     * @param logenergy  Natural logarithm of the incident photon energy [ln(keV)].
     * @return Mean energy of the scattered photon [keV].
     */
    inline double meanIncoherentScatterLogEnergy(double logenergy) const
    {
        const auto begin = m_attenuationTable.cbegin() + m_attenuationTableOffset[5];
        const auto end = m_attenuationTable.cbegin() + m_attenuationTableOffset[6];
        return std::exp(CubicLSInterpolator<double>::evaluateSpline(logenergy, begin, end));
    }

    /**
     * @brief Returns the mass energy-transfer attenuation coefficient μ_tr/ρ [cm²/g].
     *
     * Computes the attenuation values at @p energy then delegates to the overload
     * that accepts pre-computed `AttenuationValues`.
     *
     * @param energy  Photon energy [keV].
     * @return μ_tr/ρ [cm²/g].
     */
    inline double massEnergyTransferAttenuation(double energy) const
    {
        const auto att = attenuationValues(energy);
        return massEnergyTransferAttenuation(att, energy);
    }

    /**
     * @brief Returns the mass energy-transfer attenuation coefficient μ_tr/ρ [cm²/g].
     *
     * Accounts for the fluorescence escape fraction in photoelectric interactions and
     * the fraction of Compton energy retained by the recoil electron:
     *
     *   μ_tr/ρ = μ_pe/ρ × (1 − f_fluorescence) + μ_inc/ρ × (1 − E_scatter/E)
     *
     * @param att     Pre-computed `AttenuationValues` at @p energy [cm²/g each].
     * @param energy  Incident photon energy [keV].
     * @return μ_tr/ρ [cm²/g].
     */
    inline double massEnergyTransferAttenuation(const AttenuationValues& att, double energy) const
    {
        const auto logEnergy = std::log(energy);
        double avg_fluro_energy = 0;
        for (std::uint_fast8_t i = 0; i < m_numberOfShells; ++i) {
            if (m_shells[i].bindingEnergy > MIN_ENERGY())
                avg_fluro_energy += m_shells[i].energyOfPhotonsPerInitVacancy * attenuationPhotoelectricShell_logEnergy(i, logEnergy);
        }
        avg_fluro_energy /= att.photoelectric;

        const auto avg_fluro_frac = avg_fluro_energy / energy;
        const auto incoherent_scatter_energy = meanIncoherentScatterLogEnergy(logEnergy);
        const auto f_pe = 1 - avg_fluro_frac;
        const auto f_inc = 1 - (incoherent_scatter_energy + avg_fluro_energy) / energy;
        return att.photoelectric * f_pe + att.incoherent * f_inc;
    }

    /**
     * @brief Returns the photoelectric, incoherent, and coherent mass attenuation
     *        coefficients at a single energy [cm²/g].
     *
     * @param energy  Photon energy [keV].
     * @return `AttenuationValues` with photoelectric, incoherent, and coherent coefficients.
     */
    inline AttenuationValues attenuationValues(double energy) const
    {
        const auto logEnergy = std::log(energy);

        const auto begin_p = m_attenuationTable.cbegin();
        const auto begin_i = begin_p + m_attenuationTableOffset[1];
        const auto begin_c = begin_p + m_attenuationTableOffset[2];
        const auto end_c = begin_p + m_attenuationTableOffset[3];
        AttenuationValues att {
            .photoelectric = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_p, begin_i)),
            .incoherent = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_i, begin_c)),
            .coherent = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_c, end_c))
        };
        return att;
    }
    /**
     * @brief Returns the photoelectric, incoherent, and coherent mass attenuation
     *        coefficients for a vector of energies [cm²/g], evaluated in parallel.
     *
     * @param energy  Vector of photon energies [keV].
     * @return Vector of `AttenuationValues`, one entry per input energy.
     */
    std::vector<AttenuationValues> attenuationValues(const std::vector<double>& energy) const
    {
        std::vector<AttenuationValues> att(energy.size());

        const auto begin_p = m_attenuationTable.cbegin();
        const auto begin_i = begin_p + m_attenuationTableOffset[1];
        const auto begin_c = begin_p + m_attenuationTableOffset[2];
        const auto end_c = begin_p + m_attenuationTableOffset[3];

        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.begin(), [&](const auto e) {
            const auto logEnergy = std::log(e);
            AttenuationValues a {
                .photoelectric = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_p, begin_i)),
                .incoherent = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_i, begin_c)),
                .coherent = std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin_c, end_c))
            };
            return a; });
        return att;
    }
    /**
     * @brief Returns the total mass attenuation coefficient μ/ρ [cm²/g] for a vector of energies.
     *
     * Equivalent to summing the three components of `attenuationValues(energy)`, evaluated
     * in parallel.
     *
     * @param energy  Vector of photon energies [keV].
     * @return Vector of total mass attenuation coefficients [cm²/g].
     */
    inline std::vector<double> totalAttenuationValue(const std::vector<double>& energy) const
    {
        const auto att = attenuationValues(energy);
        std::vector<double> sum(att.size());
        std::transform(std::execution::par_unseq, att.cbegin(), att.cend(), sum.begin(), [](const auto& a) { return a.sum(); });
        return sum;
    }

    /**
     * @brief Returns the total photoelectric mass attenuation coefficient μ_pe/ρ [cm²/g].
     *
     * @param energy  Photon energy [keV].
     * @return Photoelectric mass attenuation coefficient [cm²/g].
     */
    inline double attenuationPhotoeletric(double energy) const
    {
        const auto logEnergy = std::log(energy);
        const auto begin = m_attenuationTable.cbegin();
        const auto end = m_attenuationTable.cbegin() + m_attenuationTableOffset[1];
        return std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin, end));
    }

    /**
     * @brief Returns the per-subshell photoelectric mass attenuation coefficient [cm²/g].
     *
     * @param shell   Subshell index in [0, numberOfShells()).
     * @param energy  Photon energy [keV].
     * @return Subshell photoelectric mass attenuation coefficient [cm²/g], or 0 if
     *         @p energy is below the subshell binding energy or @p shell is out of range.
     */
    inline double attenuationPhotoelectricShell(std::uint8_t shell, double energy) const
    {
        const auto logEnergy = std::log(energy);
        return attenuationPhotoelectricShell_logEnergy(shell, logEnergy);
    }

    /// @brief Returns the number of tracked subshells (at most `N + 1` including the remainder shell).
    inline std::uint8_t numberOfShells() const
    {
        return m_numberOfShells;
    }

    /**
     * @brief Returns the `MaterialShell` data for the given subshell index.
     * @param shell  Index in [0, numberOfShells()).
     * @return Const reference to the `MaterialShell`.
     */
    inline const MaterialShell& shell(std::size_t shell) const
    {
        return m_shells[shell];
    }

    /// @brief Returns the full array of all tracked subshells.
    inline const auto& shells() const { return m_shells; }

    /**
     * @brief Returns the momentum transfer x for a given photon energy and scattering angle.
     *
     * x = (E / hc) × sin(θ/2) [Å⁻¹]
     *
     * @param energy  Incident photon energy [keV].
     * @param angle   Scattering angle θ [rad].
     * @return Momentum transfer x [Å⁻¹].
     */
    static double momentumTransfer(double energy, double angle)
    {
        return momentumTransferMax(energy) * std::sin(angle * 0.5); // per Å
    }

    /**
     * @brief Returns the momentum transfer x for a given photon energy and cosine of scattering angle.
     *
     * x = (E / hc) × sqrt((1 − cosθ) / 2) [Å⁻¹]
     *
     * @param energy    Incident photon energy [keV].
     * @param cosAngle  Cosine of the scattering angle.
     * @return Momentum transfer x [Å⁻¹].
     */
    static double momentumTransferCosAngle(double energy, double cosAngle)
    {
        return momentumTransferMax(energy) * std::sqrt((1 - cosAngle) * 0.5); // per Å
    }

    /**
     * @brief Returns the maximum possible momentum transfer for back-scattering (θ = π) [Å⁻¹].
     *
     * x_max = E / hc [Å⁻¹]
     *
     * @param energy  Incident photon energy [keV].
     * @return Maximum momentum transfer [Å⁻¹].
     */
    static constexpr double momentumTransferMax(double energy)
    {
        static constexpr double hc_si = 1.239841193E-6; // ev*m
        static constexpr double m2A = 1E10; // meters to ångstrøm
        static constexpr double eV2keV = 1E-3; // eV to keV
        static constexpr double hc = hc_si * m2A * eV2keV; // kev*Å
        static constexpr double hc_inv = 1 / hc;
        return energy * hc_inv; // per Å
    }

    /**
     * @brief Returns true if @p cmp is a parseable chemical formula with supported elements.
     *
     * A string is valid if `parseCompoundStr` produces at least one entry and every
     * entry has Z in (0, 99] with a positive stoichiometric count.
     *
     * @param cmp  Chemical formula string to validate.
     * @return True if valid; false otherwise.
     */
    static bool validCompoundString(const std::string& cmp)
    {
        const auto m = parseCompoundStr(cmp);
        bool valid = true;
        for (const auto& [Z, w] : m) {
            valid = valid && Z > 0 && Z <= 99 && w > 0;
        }
        return valid;
    }

    /**
     * @brief Returns true if @p name is a known NIST material name.
     * @param name  Material name to look up.
     * @return True if the name exists in the NIST registry; false otherwise.
     */
    static bool validNistName(const std::string& name)
    {
        const auto& w = NISTMaterials::Composition(name);
        return !w.empty();
    }

    /**
     * @brief Parses a chemical formula string into a map of atomic number to stoichiometric count.
     *
     * Supports element symbols, integer and decimal stoichiometric coefficients, and
     * nested parentheses (e.g. "Ca3(PO4)2"). Counts are **not** normalised to mass
     * fractions — conversion to mass fractions is done by `byChemicalFormula` using
     * tabulated atomic weights.
     *
     * @param str  Chemical formula string (e.g. "H2O", "Ca3(PO4)2").
     * @return Map of `{Z → stoichiometric count}`. Empty if the string cannot be parsed.
     */
    static std::map<std::uint8_t, double> parseCompoundStr(const std::string& str)
    {
        static std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };

        std::map<std::uint8_t, double> parsed;

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
                std::uint8_t Z = std::distance(S.begin(), pos) + 1;
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
    /// @brief Default constructor; accessible only to factory methods.
    Material()
    {
    }

    /**
     * @brief Returns the per-subshell photoelectric mass attenuation coefficient [cm²/g],
     *        with the energy supplied as its natural logarithm.
     *
     * Avoids a redundant `log()` call in hot inner loops. Returns 0 if @p shell is
     * out of range or if @p logEnergy is below the subshell binding-energy threshold.
     *
     * @param shell      Subshell index in [0, numberOfShells()).
     * @param logEnergy  Natural logarithm of the photon energy [ln(keV)].
     * @return Subshell photoelectric mass attenuation coefficient [cm²/g].
     */
    inline double attenuationPhotoelectricShell_logEnergy(std::uint8_t shell, double logEnergy) const
    {
        if (shell < m_numberOfShells) {
            const auto start = m_attenuationTableOffset[6 + shell];
            const auto stop = m_attenuationTableOffset[6 + shell + 1];
            const auto begin = m_attenuationTable.cbegin() + start;
            const auto end = m_attenuationTable.cbegin() + stop;
            if (logEnergy >= (*begin)[0]) {
                return std::exp(CubicLSInterpolator<double>::evaluateSpline(logEnergy, begin, end));
            }
        }
        return 0;
    }

    /**
     * @brief Identifies which cross-section or scattering-function table to build.
     *
     * Used internally by `constructSplineInterpolator` to select the correct
     * data array from each `AtomicElement` and tune the spline knot count and
     * discontinuity handling.
     */
    enum class LUTType {
        photoelectric,    ///< Total photoelectric cross-section table.
        coherent,         ///< Coherent (Rayleigh) scattering cross-section table.
        incoherent,       ///< Incoherent (Compton) scattering cross-section table.
        formfactor,       ///< Coherent form factor F(x) table.
        scatterfactor,    ///< Incoherent scattering function S(x) table.
        incoherentenergy  ///< Mean Compton scatter energy table.
    };

    /**
     * @brief Builds a weighted least-squares cubic spline from multiple element data arrays.
     *
     * Merges all element data into a single array, applies mass-fraction weights, sorts
     * and de-duplicates the combined points, optionally converts to log-log space, then
     * fits a `CubicLSInterpolator` with @p nknots knots.
     *
     * @param data           One data array per element: `{(x, y)}` pairs (linear space).
     * @param weights        Mass fractions corresponding to each data array.
     * @param loglog         If true, both axes are log-transformed before fitting.
     * @param nknots         Number of spline knots.
     * @param maybe_discont  If true, the fitter is allowed to place knots at discontinuities
     *                       (used for photoelectric edges).
     * @return Fitted `CubicLSInterpolator`.
     */
    static CubicLSInterpolator<double> constructSplineInterpolator(const std::vector<std::vector<std::pair<double, double>>>& data, const std::vector<double>& weights, bool loglog = false, std::size_t nknots = 15, bool maybe_discont = true)
    {
        const auto Nvec = std::transform_reduce(data.cbegin(), data.cend(), std::size_t { 0 }, std::plus<>(), [](const auto& rh) -> std::size_t { return rh.size(); });

        std::vector<double> w(data.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data.size() != weights.size())
                w[i] = 1;
            else
                w[i] = weights[i];
        }
        const auto sum_weights = std::reduce(w.cbegin(), w.cend(), 0.0);
        std::transform(w.cbegin(), w.cend(), w.begin(), [sum_weights](const auto& ww) { return ww / sum_weights; });

        std::vector<std::pair<double, double>> arr(Nvec);
        std::vector<std::size_t> idx(Nvec);
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            std::copy(std::execution::par_unseq, data[i].cbegin(), data[i].cend(), arr.begin() + start);
            std::fill(idx.begin() + start, idx.begin() + start + data[i].size(), i);
            start += data[i].size();
        }
        // applying weights for first array in data
        std::transform(std::execution::par_unseq, arr.cbegin(), arr.cend(), idx.cbegin(), arr.begin(), [&](const auto& pair, const auto i) {
            return std::make_pair(pair.first, pair.second * w[i]);
        });

        for (std::size_t i = 0; i < data.size(); ++i) {
            std::transform(std::execution::par, arr.begin(), arr.end(), idx.cbegin(), arr.begin(), [&](const auto& pair, const auto index) {
                if (index == i)
                    return pair;
                else
                    return std::make_pair(pair.first, pair.second + w[i] * interpolate(data[i], pair.first));
            });
        }

        // sorting for generating interpolation lut
        std::sort(arr.begin(), arr.end(), [](const auto& lh, const auto& rh) {
            if (lh.first < rh.first)
                return true;
            else if (lh.first == rh.first)
                return lh.second < rh.second;
            return false;
        });

        // erasing duplicate items
        auto erase_from = std::unique(arr.begin(), arr.end(), [](const auto& lh, const auto& rh) {
            constexpr auto e = std::numeric_limits<double>::epsilon();
            if (lh.first == rh.first)
                return std::abs(lh.second - rh.second) <= e;
            else
                return std::abs(lh.first - rh.first) <= e;
        });
        if (std::distance(erase_from, arr.end()) != 0)
            arr.erase(erase_from, arr.end());

        if (loglog) {
            // removing zero and negative items
            auto last = std::remove_if(arr.begin(), arr.end(), [](const auto& pair) -> bool {
                return pair.first <= 0.0 || pair.second <= 0.0;
            });
            arr.erase(last, arr.end());

            std::for_each(std::execution::par_unseq, arr.begin(), arr.end(), [](auto& v) {
                v.first = std::log(v.first);
                v.second = std::log(v.second);
            });
        }

        auto interpolator = CubicLSInterpolator(arr, nknots, maybe_discont);
        return interpolator;
    }

    /**
     * @brief Builds a mass-fraction-weighted least-squares spline for the given `LUTType`.
     *
     * Gathers the appropriate raw data array from each constituent element, applies
     * type-specific preprocessing (coherent edge insertion, incoherent point densification),
     * selects the knot count and discontinuity flag, then delegates to the data-array
     * overload of `constructSplineInterpolator`.
     *
     * @param normalizedWeight  Normalised mass fractions `{Z → fraction}` (sum to 1).
     * @param type              Which cross-section or function table to build.
     * @return Fitted `CubicLSInterpolator` in log-log space.
     */
    static CubicLSInterpolator<double> constructSplineInterpolator(const std::map<std::uint8_t, double>& normalizedWeight, LUTType type)
    {
        auto getAtomArr = [&](const AtomicElement& atom, LUTType type = LUTType::photoelectric) -> const std::vector<std::pair<double, double>>& {
            if (type == LUTType::photoelectric)
                return atom.photoel;
            else if (type == LUTType::incoherent)
                return atom.incoherent;
            else if (type == LUTType::coherent)
                return atom.coherent;
            else if (type == LUTType::formfactor)
                return atom.formFactor;
            else if (type == LUTType::scatterfactor)
                return atom.incoherentSF;
            else if (type == LUTType::incoherentenergy)
                return atom.incoherentMeanScatterEnergy;
            return atom.photoel;
        };

        std::vector<std::vector<std::pair<double, double>>> data;
        std::vector<double> weights;
        for (const auto& [Z, w] : normalizedWeight) {
            const auto& a = AtomHandler::Atom(Z);
            if (type == LUTType::coherent) {
                // here we are finding the highest dip, energy wise, in coherent data and introduces an discontinuity
                auto arr = getAtomArr(a, type);
                auto max_binding_energy = MAX_ENERGY();
                for (const auto& [sidx, shell] : a.shells) {
                    const auto binding_energy = shell.bindingEnergy;
                    if (max_binding_energy > MIN_ENERGY() && binding_energy < max_binding_energy) {
                        auto pos = std::upper_bound(arr.begin(), arr.end(), shell.bindingEnergy, [](const auto v, const auto& el) { return v < el.first; });
                        if (pos != arr.end() && pos != arr.cbegin()) {
                            auto prev = pos - 1;
                            arr.insert(pos, std::make_pair(prev->first, pos->second));
                        }
                        max_binding_energy = binding_energy - 5;
                    }
                }
                data.push_back(arr);

            } else if (type == LUTType::incoherent) {
                // We increse number of datapoints for incoherent to construct a cubic spline
                auto arr = getAtomArr(a, type);
                std::vector<std::pair<double, double>> addon;
                addon.reserve(arr.size());
                for (std::size_t i = 1; i < arr.size(); ++i) {
                    const auto& first = arr[i - 1];
                    const auto& last = arr[i];
                    const auto diff = last.first - first.first;
                    constexpr std::array<double, 3> fac = {
                        0.15,
                        0.5,
                        0.85,
                    };
                    for (const auto f : fac) {
                        const auto x = first.first + diff * f;
                        const auto y = interp(first, last, x);
                        addon.push_back(std::make_pair(x, y));
                    }
                }
                for (auto& a : addon)
                    arr.push_back(a);
                data.push_back(arr);
            } else {
                data.push_back(getAtomArr(a, type));
            }

            weights.push_back(w);
        };

        std::size_t nknots = 10;
        if (type == LUTType::photoelectric) {
            nknots = 16;
        } else if (type == LUTType::coherent) {
            nknots = 32;
        } else if (type == LUTType::incoherent) {
            nknots = 32;
        } else if (type == LUTType::scatterfactor) {
            nknots = 20;
        } else if (type == LUTType::formfactor) {
            nknots = 20;
        }

        bool discont = false;
        if (type == LUTType::photoelectric || type == LUTType::coherent)
            discont = true;

        constexpr auto loglog = true;
        auto lut = constructSplineInterpolator(data, weights, loglog, nknots, discont);
        return lut;
    }

    /**
     * @brief Core factory: builds all interpolation tables and shell data for the given composition.
     *
     * Normalises the weight fractions, fits six spline LUTs (photoelectric, incoherent,
     * coherent, form factor, scatter factor, mean Compton energy), packs them into the
     * flat `m_attenuationTable` array, builds per-subshell photoelectric tables via
     * `createMaterialAtomicShells`, and precomputes the form-factor inverse-sampling
     * table via `generateFormFactorInverseSampling`.
     *
     * @param compositionByWeight  Map of `{Z → mass fraction}` (need not be normalised).
     * @return An engaged `optional<Material>` on success, `nullopt` if any Z is unknown.
     */
    static std::optional<Material<N>> constructMaterial(const std::map<std::uint8_t, double>& compositionByWeight)
    {
        for (const auto& [Z, w] : compositionByWeight) {
            const auto& a = AtomHandler::Atom(Z);
            if (a.Z != Z)
                return std::nullopt;
        }

        auto weight = compositionByWeight;
        const auto totalWeight = std::transform_reduce(weight.cbegin(), weight.cend(), 0.0, std::plus<>(), [](const auto& right) -> double { return right.second; });
        std::for_each(weight.begin(), weight.end(), [=](auto& w) { w.second /= totalWeight; });

        std::array<CubicLSInterpolator<double>, 6> attenuation = {
            constructSplineInterpolator(weight, LUTType::photoelectric),
            constructSplineInterpolator(weight, LUTType::incoherent),
            constructSplineInterpolator(weight, LUTType::coherent),
            constructSplineInterpolator(weight, LUTType::formfactor),
            constructSplineInterpolator(weight, LUTType::scatterfactor),
            constructSplineInterpolator(weight, LUTType::incoherentenergy),
        };
        Material<N> m;
        m.m_composition = weight; // setting composition

        m.m_attenuationTable.clear();
        std::array<std::size_t, attenuation.size() + N> offset;
        for (std::size_t i = 0; i < attenuation.size(); ++i) {
            const auto& table = attenuation[i].getDataTable();
            auto begin = m.m_attenuationTable.insert(m.m_attenuationTable.end(), table.cbegin(), table.cend());
            offset[i] = std::distance(m.m_attenuationTable.begin(), begin);
        }
        createMaterialAtomicShells(m, weight, offset);

        for (std::size_t i = 0; i < offset.size(); ++i) {
            if (i < attenuation.size() + std::min(std::uint8_t { N }, m.numberOfShells()))
                m.m_attenuationTableOffset[i] = static_cast<std::uint_fast32_t>(std::distance(m.m_attenuationTable.begin(), m.m_attenuationTable.begin() + offset[i]));
            else
                m.m_attenuationTableOffset[i] = static_cast<std::uint_fast32_t>(std::distance(m.m_attenuationTable.begin(), m.m_attenuationTable.end()));
        }
        m.m_attenuationTableOffset[offset.size()] = static_cast<std::uint_fast32_t>(std::distance(m.m_attenuationTable.begin(), m.m_attenuationTable.end()));

        // creating lookuptable for inverse sampling of formfactor
        generateFormFactorInverseSampling(m);

        return m;
    }

    /**
     * @brief Precomputes the CPDF inverse-sampling table for the coherent form factor squared.
     *
     * Builds a `CPDFSampling<double, 20>` over q² ∈ [q²_min, q²_max] proportional to
     * F²(q), enabling rejection-free sampling of the coherent scattering angle.
     *
     * @param material  Material to update in-place (sets `m_formFactorInvSamp`).
     */
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

    /**
     * @brief Builds the per-subshell photoelectric tables and `MaterialShell` records.
     *
     * Collects all subshells from all constituent elements, sorts them by binding energy
     * (highest first), keeps the top `N` explicitly (each with its own spline), and
     * collapses the remainder into a single averaged shell. Also normalises the
     * `numberOfElectronsFraction` across all shells.
     *
     * @param material          Material to update in-place.
     * @param normalizedWeight  Normalised mass fractions `{Z → fraction}`.
     * @param offset            Output array of table offsets updated by this function
     *                          for indices 6 through 6+N.
     */
    static void createMaterialAtomicShells(Material<N>& material, const std::map<std::uint8_t, double>& normalizedWeight, std::array<std::size_t, 6 + N>& offset)
    {
        struct Shell {
            std::uint64_t Z = 0;
            std::uint64_t S = 0;
            double weight = 0;
            double bindingEnergy = 0;
        };

        std::vector<Shell> shells;
        for (const auto& [Z, w] : normalizedWeight) {
            const auto& atom = AtomHandler::Atom(Z);
            for (const auto& [S, shell] : atom.shells) {
                shells.push_back({ .Z = Z, .S = S, .weight = w, .bindingEnergy = shell.bindingEnergy });
            }
        }
        std::sort(shells.begin(), shells.end(), [](const auto& lh, const auto& rh) {
            return lh.bindingEnergy > rh.bindingEnergy;
        });

        const auto sum_weight = std::transform_reduce(shells.cbegin(), shells.cend(), 0.0, std::plus<>(), [](const auto& s) { return s.weight; });

        const auto Nshells = std::min(shells.size(), N);
        material.m_numberOfShells = static_cast<std::uint8_t>(Nshells);
        for (std::size_t i = 0; i < Nshells; ++i) {
            const auto& shell = AtomHandler::Atom(shells[i].Z).shells.at(shells[i].S);
            std::vector<std::pair<double, double>> photolog(shell.photoel.size());
            std::transform(std::execution::par_unseq, shell.photoel.cbegin(), shell.photoel.cend(), photolog.begin(),
                [=](const auto& p) {
                    return std::make_pair(std::log(p.first), std::log(p.second * shells[i].weight));
                });
            CubicLSInterpolator<double> inter(photolog, 5, false);

            auto begin = inter.getDataTable().begin();
            auto end = inter.getDataTable().end();
            auto table_beg = material.m_attenuationTable.insert(material.m_attenuationTable.end(), begin, end);
            offset[i + 6] = std::distance(material.m_attenuationTable.begin(), table_beg);

            auto& materialshell = material.m_shells[i];
            materialshell.numberOfElectronsFraction = shells[i].weight * shell.numberOfElectrons / sum_weight;
            materialshell.bindingEnergy = shell.bindingEnergy;
            materialshell.HartreeFockOrbital_0 = shell.HartreeFockOrbital_0;
            materialshell.numberOfPhotonsPerInitVacancy = shell.numberOfPhotonsPerInitVacancy;
            materialshell.energyOfPhotonsPerInitVacancy = shell.energyOfPhotonsPerInitVacancy;
        }
        // Filling remainder shell
        if (shells.size() > Nshells) {
            material.m_numberOfShells++;
            auto& materialshell = material.m_shells[Nshells];
            const auto mean_fac = 1.0 / (shells.size() - Nshells);
            for (std::size_t i = Nshells; i < shells.size(); ++i) {
                const auto& shell = AtomHandler::Atom(shells[i].Z).shells.at(shells[i].S);
                const auto w = shells[i].weight;
                materialshell.bindingEnergy += shell.bindingEnergy * mean_fac;
                materialshell.numberOfElectronsFraction += shell.numberOfElectrons * w / sum_weight;
                materialshell.HartreeFockOrbital_0 += shell.HartreeFockOrbital_0 * mean_fac;
                materialshell.numberOfPhotonsPerInitVacancy += shell.numberOfPhotonsPerInitVacancy * mean_fac;
                materialshell.energyOfPhotonsPerInitVacancy += shell.energyOfPhotonsPerInitVacancy * mean_fac;
            }
        }

        // normalize number of electrons fraction
        const auto sumElFraction = std::transform_reduce(material.m_shells.cbegin(), material.m_shells.cend(), 0.0, std::plus<>(), [](const auto& s) { return s.numberOfElectronsFraction; });
        std::for_each(material.m_shells.begin(), material.m_shells.end(), [sumElFraction](auto& s) { s.numberOfElectronsFraction /= sumElFraction; });
    }

private:
    std::array<std::uint_fast32_t, 6 + N + 1> m_attenuationTableOffset;  ///< Start indices into `m_attenuationTable` for each of the 6 LUTs + N shell tables + sentinel.
    CPDFSampling<double, 20> m_formFactorInvSamp;                          ///< Inverse-CDF sampler for coherent scattering angle via F²(q).
    std::array<MaterialShell, N + 1> m_shells;                             ///< Per-subshell data (N explicit shells + 1 remainder shell).
    std::uint8_t m_numberOfShells = 0;                                     ///< Actual number of shells populated (≤ N + 1).
    std::vector<std::array<double, 3>> m_attenuationTable;                 ///< Flat packed spline knot table: {log_x, log_y, dy/dx} for all LUTs.
    std::map<std::uint8_t, double> m_composition;                          ///< Normalised elemental mass fractions {Z → fraction}.
};

}
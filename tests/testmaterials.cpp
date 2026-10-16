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

#include "testutils.hpp"

#include "xraymc/constants.hpp"
#include "xraymc/material/atomhandler.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/material/nistmaterials.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

// Reference data from NIST XCOM / Hubbell & Seltzer tables:
// energy [keV], mass attenuation coefficient mu/rho [cm2/g], mass energy absorption coefficient mu_en/rho [cm2/g]
struct NISTData {
    double energy, att, atten;
};

std::vector<NISTData> nistWater()
{
    return {
        { 10, 5.329, 4.944 },
        { 15, 1.673, 1.374 },
        { 20, 0.8096, 0.5503 },
        { 30, 0.3756, 0.1557 },
        { 40, 0.2683, 0.06947 },
        { 50, 0.2269, 0.04223 },
        { 60, 0.2059, 0.03190 },
        { 80, 0.1837, 0.02597 },
        { 100, 0.1707, 0.02546 },
        { 150, 0.1505, 0.02764 }
    };
}

std::vector<NISTData> nistAir()
{
    // Argon K-edge at 3.2029 keV is excluded
    return {
        { 1.5, 1.191E+03, 1.188E+03 },
        { 2.0, 5.279E+02, 5.262E+02 },
        { 3.0, 1.625E+02, 1.614E+02 },
        { 4.0, 7.788E+01, 7.636E+01 },
        { 5.0, 4.027E+01, 3.931E+01 },
        { 6.0, 2.341E+01, 2.270E+01 },
        { 8.0, 9.921E+00, 9.446E+00 },
        { 10, 5.120E+00, 4.742E+00 },
        { 15, 1.614E+00, 1.334E+00 },
        { 20, 7.779E-01, 5.389E-01 },
        { 30, 3.538E-01, 1.537E-01 },
        { 40, 2.485E-01, 6.833E-02 },
        { 50, 2.080E-01, 4.098E-02 },
        { 60, 1.875E-01, 3.041E-02 },
        { 80, 1.662E-01, 2.407E-02 },
        { 100, 1.541E-01, 2.325E-02 },
        { 150, 1.356E-01, 2.496E-02 }
    };
}

std::vector<NISTData> nistLead()
{
    // energies close to absorption edges are excluded, mu_en not used
    return {
        { 10, 130.6, 0 },
        { 20, 86.36, 0 },
        { 30, 30.32, 0 },
        { 40, 14.36, 0 },
        { 50, 8.041, 0 },
        { 60, 5.021, 0 },
        { 80, 2.419, 0 },
        { 100, 5.549, 0 },
        { 150, 2.014, 0 }
    };
}

void testAtomHandler()
{
    const auto& atoms = xraymc::AtomHandler::allAtoms();
    check(atoms.size() == 100, "AtomHandler: physics data contains elements 1 to 100 (" + std::to_string(atoms.size()) + ")");
    check(xraymc::AtomHandler::atomExists(1) && xraymc::AtomHandler::atomExists(100), "AtomHandler: first and last element exist");
    check(!xraymc::AtomHandler::atomExists(0) && !xraymc::AtomHandler::atomExists(101), "AtomHandler: Z = 0 and Z = 101 do not exist");
    check(xraymc::AtomHandler::toSymbol(82) == "Pb" && xraymc::AtomHandler::toSymbol(1) == "H", "AtomHandler: element symbols");
    check(xraymc::AtomHandler::toSymbol(0).empty(), "AtomHandler: no symbol for Z = 0");

    const auto& pb = xraymc::AtomHandler::Atom(82);
    check(pb.Z == 82, "AtomHandler: lead has Z = 82");
    checkRelative(pb.atomicWeight, 207.2, 1E-3, "AtomHandler: lead atomic weight");
    checkRelative(pb.standardDensity, 11.35, 1E-2, "AtomHandler: lead standard density");
    check(xraymc::AtomHandler::Atom(0).Z == 0, "AtomHandler: unknown element returns dummy element");

    bool tablesOk = true;
    for (const auto& [Z, atom] : atoms) {
        tablesOk = tablesOk && atom.Z == Z && !atom.photoel.empty() && !atom.incoherent.empty() && !atom.coherent.empty() && !atom.shells.empty();
        tablesOk = tablesOk && atom.atomicWeight > 0;
    }
    check(tablesOk, "AtomHandler: all elements have cross section tables and shells");

    // K-shell binding energy of lead is 88.0 keV
    double maxBinding = 0;
    for (const auto& [id, shell] : pb.shells)
        maxBinding = std::max(maxBinding, shell.bindingEnergy);
    checkNear(maxBinding, 88.0, 0.1, "AtomHandler: lead K-shell binding energy");
}

void testNISTMaterials()
{
    const auto names = xraymc::NISTMaterials::listNames();
    check(names.size() > 100, "NISTMaterials: registry contains materials (" + std::to_string(names.size()) + ")");
    bool compOk = true;
    bool densOk = true;
    std::string bad;
    for (const auto& name : names) {
        const auto comp = xraymc::NISTMaterials::Composition(name);
        const auto sum = std::transform_reduce(comp.cbegin(), comp.cend(), 0.0, std::plus<>(), [](const auto& p) { return p.second; });
        if (comp.empty() || std::abs(sum - 1) > 2E-3) {
            compOk = false;
            bad = name;
        }
        densOk = densOk && xraymc::NISTMaterials::density(name) > 0;
    }
    check(compOk, "NISTMaterials: all mass fractions sum to one" + (bad.empty() ? "" : " (" + bad + ")"));
    check(densOk, "NISTMaterials: all densities positive");
    checkNear(xraymc::NISTMaterials::density("Water, Liquid"), 1.0, 1E-6, "NISTMaterials: water density");
    check(xraymc::NISTMaterials::density("No such material") == 0, "NISTMaterials: unknown material density is 0");
    check(xraymc::NISTMaterials::Composition("No such material").empty(), "NISTMaterials: unknown material has empty composition");
}

void testConstruction()
{
    using Material = xraymc::Material<5>;

    bool allZ = true;
    for (std::size_t Z = 1; Z <= 100; ++Z)
        allZ = allZ && Material::byZ(Z).has_value();
    check(allZ, "Material: byZ constructs all elements 1 to 100");
    check(!Material::byZ(0).has_value(), "Material: byZ(0) returns nullopt");
    check(!Material::byZ(101).has_value(), "Material: byZ(101) returns nullopt");
    check(!Material::byZ(300).has_value(), "Material: byZ(300) returns nullopt (no narrowing to valid Z)");
    check(!Material::byZ(-1).has_value(), "Material: byZ(-1) returns nullopt");

    bool allNist = true;
    std::string bad;
    for (const auto& name : Material::listNistCompoundNames()) {
        if (!Material::byNistName(name)) {
            allNist = false;
            bad = name;
        }
    }
    check(allNist, "Material: byNistName constructs all NIST materials" + (bad.empty() ? "" : " (" + bad + ")"));
    check(!Material::byNistName("No such material").has_value(), "Material: byNistName with unknown name returns nullopt");
    check(Material::validNistName("Water, Liquid") && !Material::validNistName("Water"), "Material: validNistName");

    const auto water = Material::byNistName("Water, Liquid").value();
    const auto waterFormula = Material::byChemicalFormula("H2O").value();
    const auto& wc = waterFormula.composition();
    check(wc.size() == 2, "Material: H2O has two elements");
    checkNear(wc.at(1), 0.111894, 1E-4, "Material: H2O hydrogen mass fraction");
    checkNear(wc.at(8), 0.888106, 1E-4, "Material: H2O oxygen mass fraction");
    checkRelative(waterFormula.attenuationValues(60).sum(), water.attenuationValues(60).sum(), 1E-3, "Material: H2O by formula equals NIST water");
    checkNear(water.effectiveZ(), 0.111894 * 1 + 0.888106 * 8, 1E-4, "Material: effectiveZ is mass weighted Z");

    // byWeight normalizes weights
    const auto m1 = Material::byWeight(std::map<int, double> { { 1, 2.0 }, { 8, 6.0 } }).value();
    checkNear(m1.composition().at(1), 0.25, 1E-12, "Material: byWeight normalizes weights");
    checkNear(m1.composition().at(8), 0.75, 1E-12, "Material: byWeight normalizes weights (second element)");
    check(!Material::byWeight(std::map<int, double> {}).has_value(), "Material: byWeight with empty composition returns nullopt");
    check(!Material::byWeight(std::map<int, double> { { 0, 1.0 } }).has_value(), "Material: byWeight with Z = 0 returns nullopt");
    check(!Material::byWeight(std::map<int, double> { { 1, 0.0 } }).has_value(), "Material: byWeight with zero total weight returns nullopt");
    check(!Material::byChemicalFormula("Xx2").has_value(), "Material: byChemicalFormula with unknown element returns nullopt");
    check(!Material::byChemicalFormula("").has_value(), "Material: byChemicalFormula with empty string returns nullopt");

    check(Material::byZ(13).value() == Material::byZ(13).value(), "Material: equal materials compare equal");
    check(!(Material::byZ(13).value() == Material::byZ(14).value()), "Material: different materials compare unequal");
}

void testParseCompound()
{
    using Material = xraymc::Material<5>;
    auto eq = [](const std::map<std::uint8_t, double>& a, const std::map<std::uint8_t, double>& b) {
        if (a.size() != b.size())
            return false;
        for (const auto& [Z, n] : a)
            if (!b.contains(Z) || std::abs(b.at(Z) - n) > 1E-12)
                return false;
        return true;
    };
    check(eq(Material::parseCompoundStr("H2O"), { { 1, 2 }, { 8, 1 } }), "parseCompoundStr: H2O");
    check(eq(Material::parseCompoundStr("C6H12O6"), { { 6, 6 }, { 1, 12 }, { 8, 6 } }), "parseCompoundStr: C6H12O6");
    check(eq(Material::parseCompoundStr("NaCl"), { { 11, 1 }, { 17, 1 } }), "parseCompoundStr: two letter symbols");
    check(eq(Material::parseCompoundStr("H0.5O"), { { 1, 0.5 }, { 8, 1 } }), "parseCompoundStr: fractional number");
    check(eq(Material::parseCompoundStr("Ca10(PO4)6(OH)2"), { { 20, 10 }, { 15, 6 }, { 8, 26 }, { 1, 2 } }), "parseCompoundStr: hydroxyapatite with parenthesis");
    check(eq(Material::parseCompoundStr("K(OH)"), { { 19, 1 }, { 8, 1 }, { 1, 1 } }), "parseCompoundStr: parenthesis without multiplier");
    check(eq(Material::parseCompoundStr("HH"), { { 1, 2 } }), "parseCompoundStr: repeated element is summed");
    check(Material::parseCompoundStr("Xx2").empty(), "parseCompoundStr: unknown element gives empty result");
    check(Material::parseCompoundStr("XxO").empty(), "parseCompoundStr: unknown element among known gives empty result");
    check(Material::validCompoundString("H2O"), "validCompoundString: H2O is valid");
    check(!Material::validCompoundString("Xx2"), "validCompoundString: unknown element is invalid");
    check(!Material::validCompoundString(""), "validCompoundString: empty string is invalid");
}

void testReferenceAttenuation()
{
    using Material = xraymc::Material<16>;
    const auto water = Material::byNistName("Water, Liquid").value();
    for (const auto& d : nistWater()) {
        checkRelative(water.attenuationValues(d.energy).sum(), d.att, 0.01, "Water: total attenuation vs NIST at " + std::to_string(int(d.energy)) + " keV");
        checkRelative(water.massEnergyTransferAttenuation(d.energy), d.atten, 0.01, "Water: energy transfer vs NIST mu_en at " + std::to_string(int(d.energy)) + " keV");
    }
    const auto air = Material::byNistName("Air, Dry (near sea level)").value();
    for (const auto& d : nistAir()) {
        checkRelative(air.attenuationValues(d.energy).sum(), d.att, 0.02, "Air: total attenuation vs NIST at " + std::to_string(d.energy) + " keV");
        checkRelative(air.massEnergyTransferAttenuation(d.energy), d.atten, 0.02, "Air: energy transfer vs NIST mu_en at " + std::to_string(d.energy) + " keV");
    }
    const auto lead = Material::byZ(82).value();
    for (const auto& d : nistLead()) {
        checkRelative(lead.attenuationValues(d.energy).sum(), d.att, 0.03, "Lead: total attenuation vs NIST at " + std::to_string(int(d.energy)) + " keV");
    }
    const auto below = lead.attenuationValues(87.9).sum();
    const auto above = lead.attenuationValues(88.1).sum();
    check(above / below > 3.5, "Lead: K-edge jump at 88 keV (ratio " + std::to_string(above / below) + ")");
}

void testAttenuationProperties()
{
    using Material = xraymc::Material<16>;
    const auto water = Material::byNistName("Water, Liquid").value();

    bool positive = true;
    bool peDecreasing = true;
    bool totalDecreasing = true;
    bool mutrLess = true;
    double lastPe = std::numeric_limits<double>::max();
    double lastTot = std::numeric_limits<double>::max();
    std::vector<double> energies;
    for (double e = 1.0; e <= 150.0; e += 0.5) {
        energies.push_back(e);
        const auto a = water.attenuationValues(e);
        positive = positive && a.photoelectric > 0 && a.incoherent > 0 && a.coherent > 0;
        peDecreasing = peDecreasing && a.photoelectric < lastPe;
        totalDecreasing = totalDecreasing && a.sum() < lastTot;
        mutrLess = mutrLess && water.massEnergyTransferAttenuation(e) < a.sum() && water.massEnergyTransferAttenuation(e) > 0;
        lastPe = a.photoelectric;
        lastTot = a.sum();
    }
    check(positive, "Water: all attenuation components positive 1-150 keV");
    check(peDecreasing, "Water: photoelectric attenuation decreases with energy (no edges)");
    check(totalDecreasing, "Water: total attenuation decreases with energy");
    check(mutrLess, "Water: 0 < mu_tr < mu");

    const auto vec = water.attenuationValues(energies);
    const auto tot = water.totalAttenuationValue(energies);
    bool vecOk = vec.size() == energies.size() && tot.size() == energies.size();
    for (std::size_t i = 0; i < energies.size() && vecOk; ++i) {
        const auto a = water.attenuationValues(energies[i]);
        vecOk = std::abs(vec[i].photoelectric - a.photoelectric) <= 1E-12 * a.photoelectric && std::abs(tot[i] - a.sum()) <= 1E-12 * a.sum();
    }
    check(vecOk, "Water: vector attenuation equals scalar attenuation");
    checkRelative(water.attenuationPhotoeletric(30), water.attenuationValues(30).photoelectric, 1E-12, "Water: attenuationPhotoeletric equals photoelectric component");

    // shell data
    const auto lead = Material::byZ(82).value();
    double frac = 0;
    bool sorted = true;
    for (std::size_t i = 0; i < lead.numberOfShells(); ++i) {
        frac += lead.shell(i).numberOfElectronsFraction;
        if (i > 0)
            sorted = sorted && lead.shell(i).bindingEnergy <= lead.shell(i - 1).bindingEnergy;
    }
    checkNear(frac, 1.0, 1E-6, "Lead: shell electron fractions sum to one");
    check(sorted, "Lead: shells sorted by decreasing binding energy");
    checkNear(lead.shell(0).bindingEnergy, 88.0, 0.1, "Lead: first shell is K-shell");

    // 16 explicit shells covers all but the outermost shells of lead, the remainder shell has no table
    check(lead.numberOfShells() == 17, "Lead: 16 explicit shells and a remainder shell");
    checkNear(lead.attenuationPhotoelectricShell(16, 60.0), 0.0, 0.0, "Lead: remainder shell has no photoelectric table");
    for (double e : { 20.0, 60.0, 100.0, 140.0 }) {
        double shellSum = 0;
        for (std::size_t i = 0; i < lead.numberOfShells(); ++i)
            if (lead.shell(i).bindingEnergy < e)
                shellSum += lead.attenuationPhotoelectricShell(static_cast<std::uint8_t>(i), e);
        checkRelative(shellSum, lead.attenuationValues(e).photoelectric, 0.02, "Lead: sum of shell photoelectric cross sections equals total at " + std::to_string(int(e)) + " keV");
    }
    const auto lead5 = xraymc::Material<5>::byZ(82).value();
    bool finite = true;
    for (std::size_t i = 0; i <= 5; ++i)
        finite = finite && std::isfinite(lead5.attenuationPhotoelectricShell(static_cast<std::uint8_t>(i), 60.0));
    check(finite, "Lead with 5 shells: all shell cross sections finite including remainder");
    checkRelative(lead5.massEnergyTransferAttenuation(60), lead.massEnergyTransferAttenuation(60), 1E-3, "Lead: energy transfer independent of number of explicit shells");
}

void testScatterFunctions()
{
    using Material = xraymc::Material<16>;
    const auto water = Material::byNistName("Water, Liquid").value();

    // momentum transfer x = sin(theta/2) / lambda in inverse Angstrom
    const double E = 60;
    checkRelative(Material::momentumTransferMax(E), E / 12.398520, 1E-5, "momentumTransferMax is E / hc");
    checkRelative(Material::momentumTransfer(E, std::numbers::pi), Material::momentumTransferMax(E), 1E-12, "momentumTransfer at 180 deg is maximum");
    checkNear(Material::momentumTransfer(E, 0), 0, 1E-12, "momentumTransfer at 0 deg is zero");
    checkRelative(Material::momentumTransferCosAngle(E, std::cos(0.7)), Material::momentumTransfer(E, 0.7), 1E-10, "momentumTransferCosAngle agrees with momentumTransfer");

    bool ffDecreasing = true;
    bool sfIncreasing = true;
    double lastFF = std::numeric_limits<double>::max();
    double lastSF = 0;
    for (double x = 0.01; x < 3; x *= 1.2) {
        const auto ff = water.formFactor(x);
        const auto sf = water.scatterFactor(x);
        ffDecreasing = ffDecreasing && ff <= lastFF * (1 + 1E-6);
        sfIncreasing = sfIncreasing && sf >= lastSF * (1 - 1E-6);
        lastFF = ff;
        lastSF = sf;
    }
    check(ffDecreasing, "Water: form factor decreases with momentum transfer");
    check(sfIncreasing, "Water: incoherent scatter function increases with momentum transfer");
    // Both functions are mass fraction weighted sums over elements and tend to sum(w_i Z_i)
    checkRelative(water.formFactor(1E-3), water.effectiveZ(), 0.01, "Water: form factor at zero momentum transfer equals effective Z");
    checkRelative(water.scatterFactor(9.0), water.effectiveZ(), 0.01, "Water: scatter function at large momentum transfer equals effective Z");
    checkNear(water.scatterFactor(1E-3), 0, 0.01, "Water: scatter function vanishes at zero momentum transfer");

    bool meanOk = true;
    for (double e = 5; e <= 150; e += 5) {
        const auto k = e / xraymc::ELECTRON_REST_MASS();
        const auto emin = e / (1 + 2 * k);
        const auto mean = water.meanIncoherentScatterEnergy(e);
        meanOk = meanOk && mean < e && mean > emin;
    }
    check(meanOk, "Water: mean Compton scattered energy between kinematic limits");

    auto state = xraymctest::seededState();
    const auto qmax = Material::momentumTransferMax(60);
    bool qOk = true;
    for (int i = 0; i < 100000; ++i) {
        const auto q2 = water.sampleSquaredMomentumTransferFromFormFactorSquared(qmax * qmax, state);
        qOk = qOk && q2 >= 0 && q2 <= qmax * qmax;
    }
    check(qOk, "Water: sampled squared momentum transfer within [0, qmax^2]");
}

int main()
{
    std::cout << "Testing materials and attenuation data\n";
    testAtomHandler();
    testNISTMaterials();
    testConstruction();
    testParseCompound();
    testReferenceAttenuation();
    testAttenuationProperties();
    testScatterFunctions();
    return xraymctest::summary("testmaterials");
}

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

#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/material/atomhandler.hpp"
#include "xraymc/material/material.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;
using xraymctest::checkRelative;

// Air KERMA transmission through mm of aluminium, computed independently of Tube
double kermaTransmission(const Tube& t, double mmAl)
{
    const auto e = t.getEnergy();
    const auto s = t.getSpecter(e);
    const auto& al = AtomHandler::Atom(13);
    const auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
    double I0 = 0, I1 = 0;
    for (std::size_t i = 0; i < e.size(); ++i) {
        if (e[i] < MIN_ENERGY())
            continue;
        const auto att = (interpolate(al.photoel, e[i]) + interpolate(al.incoherent, e[i]) + interpolate(al.coherent, e[i])) * al.standardDensity;
        const auto k = e[i] * s[i] * air.massEnergyTransferAttenuation(e[i]);
        I0 += k;
        I1 += k * std::exp(-att * mmAl * 0.1);
    }
    return I1 / I0;
}

void testSpecter()
{
    Tube t(120, 12, 1.0);
    t.setAlFiltration(2.5);
    const auto energy = t.getEnergy();
    check(!energy.empty() && std::abs(energy.front() - 1.0) < 1E-12 && std::abs(energy.back() - 120.0) < 1E-9, "Tube: energy grid from resolution to tube voltage");
    bool spaced = true;
    for (std::size_t i = 1; i < energy.size(); ++i)
        spaced = spaced && std::abs(energy[i] - energy[i - 1] - 1.0) < 1E-9;
    check(spaced, "Tube: energy grid spacing equals resolution");

    const auto s = t.getSpecter();
    const auto sum = std::transform_reduce(s.cbegin(), s.cend(), 0.0, std::plus<>(), [](const auto& p) { return p.second; });
    checkNear(sum, 1.0, 1E-9, "Tube: normalized spectrum sums to one");
    bool nonNegative = true;
    for (const auto& [e, n] : s)
        nonNegative = nonNegative && n >= 0 && e <= 120;
    check(nonNegative, "Tube: spectrum is non-negative and below tube voltage");
    check(s.front().second < 1E-6, "Tube: low energies are removed by filtration");

    const auto sUnnorm = t.getSpecter(false);
    const auto sumU = std::transform_reduce(sUnnorm.cbegin(), sUnnorm.cend(), 0.0, std::plus<>(), [](const auto& p) { return p.second; });
    bool proportional = sUnnorm.size() == s.size() && sumU > 0;
    for (std::size_t i = 0; i < s.size() && proportional; ++i)
        proportional = std::abs(sUnnorm[i].second / sumU - s[i].second) < 1E-9;
    check(proportional, "Tube: unnormalized spectrum is proportional to normalized");

    // Tungsten characteristic radiation: K-alpha at 58-59 keV and K-beta at 67 keV, K-edge at 69.5 keV
    auto value = [&](const auto& sp, double e) {
        for (const auto& [ei, n] : sp)
            if (std::abs(ei - e) < 1E-6)
                return n;
        return 0.0;
    };
    check(value(s, 58) > 1.5 * value(s, 57) && value(s, 60) > 1.5 * value(s, 61), "Tube: tungsten K-alpha lines present at 120 kV");
    check(value(s, 68) > 1.5 * value(s, 69), "Tube: tungsten K-beta line present at 120 kV");

    Tube t55(55, 12, 1.0);
    t55.setAlFiltration(2.5);
    const auto s55 = t55.getSpecter();
    check(value(s55, 54) < value(s55, 53), "Tube: no tungsten K lines below the K-edge (55 kV)");
}

void testVoltageAndFiltration()
{
    double lastMean = 0, lastHVL = 0;
    bool increasing = true;
    for (double kv : { 50.0, 70.0, 90.0, 110.0, 130.0, 150.0 }) {
        Tube t(kv);
        t.setAlFiltration(2.5);
        increasing = increasing && t.meanSpecterEnergy() > lastMean && t.mmAlHalfValueLayer() > lastHVL;
        lastMean = t.meanSpecterEnergy();
        lastHVL = t.mmAlHalfValueLayer();
    }
    check(increasing, "Tube: mean energy and HVL increase with tube voltage");

    lastMean = 0;
    lastHVL = 0;
    increasing = true;
    for (double mm : { 0.5, 1.0, 2.0, 4.0, 8.0 }) {
        Tube t(100);
        t.setAlFiltration(mm);
        increasing = increasing && t.meanSpecterEnergy() > lastMean && t.mmAlHalfValueLayer() > lastHVL;
        lastMean = t.meanSpecterEnergy();
        lastHVL = t.mmAlHalfValueLayer();
    }
    check(increasing, "Tube: mean energy and HVL increase with aluminium filtration");

    Tube tc(100);
    tc.setAlFiltration(2.5);
    const auto hvlAl = tc.mmAlHalfValueLayer();
    tc.setCuFiltration(0.2);
    check(tc.mmAlHalfValueLayer() > hvlAl, "Tube: added copper filtration hardens the beam");
    checkNear(tc.CuFiltration(), 0.2, 1E-12, "Tube: copper filtration is stored");
    checkNear(tc.AlFiltration(), 2.5, 1E-12, "Tube: aluminium filtration is kept when adding copper");
    tc.clearFiltrationMaterials();
    check(tc.AlFiltration() == 0 && tc.CuFiltration() == 0, "Tube: clearFiltrationMaterials removes all filtration");
    check(!tc.addFiltrationMaterial(0, 1.0), "Tube: invalid filtration material is rejected");

    Tube tv;
    tv.setVoltage(500);
    check(tv.voltage() == Tube::maxVoltage(), "Tube: voltage is clamped to maximum");
    tv.setVoltage(10);
    check(tv.voltage() == Tube::minVoltage(), "Tube: voltage is clamped to minimum");
    tv.setAnodeAngleDeg(15);
    checkNear(tv.anodeAngleDeg(), 15, 1E-12, "Tube: anode angle in degrees");
}

void testHalfValueLayer()
{
    for (double kv : { 60.0, 80.0, 120.0 }) {
        for (double res : { 1.0, 0.5 }) {
            Tube t(kv, 12, res);
            t.setAlFiltration(3.0);
            const auto hvl = t.mmAlHalfValueLayer();
            const auto label = "Tube " + std::to_string(int(kv)) + " kV (resolution " + std::to_string(res) + " keV)";
            checkNear(kermaTransmission(t, hvl), 0.5, 0.01, label + ": air KERMA is halved by one HVL");
            check(hvl > 1.0 && hvl < 10.0, label + ": HVL in plausible range (" + std::to_string(hvl) + " mm Al)");
        }
    }
    // HVL must not depend on the energy resolution of the spectrum
    Tube a(80, 12, 1.0), b(80, 12, 0.5);
    a.setAlFiltration(3.0);
    b.setAlFiltration(3.0);
    checkRelative(b.mmAlHalfValueLayer(), a.mmAlHalfValueLayer(), 0.03, "Tube: HVL independent of energy resolution");
    checkRelative(b.meanSpecterEnergy(), a.meanSpecterEnergy(), 0.02, "Tube: mean energy independent of energy resolution");

    // const and non-const overloads agree
    const Tube c(90);
    Tube d(90);
    checkNear(c.mmAlHalfValueLayer(), d.mmAlHalfValueLayer(), 1E-12, "Tube: const and cached HVL agree");
    checkNear(c.meanSpecterEnergy(), d.meanSpecterEnergy(), 1E-12, "Tube: const and cached mean energy agree");
}

void testSerialize()
{
    Tube t(85, 10, 0.5);
    t.setAlFiltration(2.0);
    t.setCuFiltration(0.1);
    const auto buffer = t.serialize();
    const auto r = Tube::deserialize(buffer);
    check(r.has_value(), "Tube: deserialize succeeds");
    if (r) {
        const auto& u = r.value();
        check(u.voltage() == t.voltage() && u.anodeAngle() == t.anodeAngle() && u.AlFiltration() == t.AlFiltration() && u.CuFiltration() == t.CuFiltration(), "Tube: deserialized parameters equal original");
        // spectra are normalized with parallel reductions, allow for rounding differences
        const auto su = u.getSpecter();
        const auto st = t.getSpecter();
        bool same = su.size() == st.size();
        for (std::size_t i = 0; i < su.size() && same; ++i)
            same = su[i].first == st[i].first && std::abs(su[i].second - st[i].second) <= 1E-12;
        check(same, "Tube: deserialized spectrum equals original");
    }
}

int main()
{
    std::cout << "Testing x-ray tube\n";
    testSpecter();
    testVoltageAndFiltration();
    testHalfValueLayer();
    testSerialize();
    return xraymctest::summary("testtube");
}

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

#include "atomicelementhandler.hpp"

#include <charconv>
#include <cmath>
#include <optional>
#include <string_view>

#ifdef USE_NIST
#include "xraylib++.h"
#endif

std::size_t lowerIdx(const std::vector<double>& data, double limit)
{
    for (std::size_t i = 0; i < data.size() - 2; i = i + 2) {
        if (data[i + 2] > limit)
            return i;
    }
    return 0;
}
std::size_t upperIdx(const std::vector<double>& data, double limit)
{
    for (int i = data.size() - 2; i > 0; i = i - 2) {
        if (data[i - 2] < limit)
            return i;
    }
    return data.size() - 2;
}

AtomicElementHandler::AtomicElementHandler(std::uint64_t Z)
{
    m_atom.Z = Z;
}

void AtomicElementHandler::setCoherentData(const std::vector<double>& data)
{
    auto start = lowerIdx(data, minPhotonEnergy() * keVToMeV());
    auto end = upperIdx(data, maxPhotonEnergy() * keVToMeV());
    m_atom.coherent.clear();
    m_atom.coherent.reserve((end - start) / 2);
    for (std::size_t i = start; i <= end; i = i + 2) {
        const double e = data[i] * MeVTokeV();
#ifdef USE_NIST
        double a;
        try {
            a = xrlpp::CS_Rayl(static_cast<int>(m_atom.Z), e);
        } catch (std::invalid_argument& e) {
            a = data[i + 1] * barnToAtt();
        }
#else
        const double a = data[i + 1] * barnToAtt();
#endif
        if (e > 0 && a > 0)
            m_atom.coherent.push_back(std::make_pair(e, a));
    }
    m_atom.coherent.shrink_to_fit();
}
void AtomicElementHandler::setIncoherentData(const std::vector<double>& data)
{
    auto start = lowerIdx(data, minPhotonEnergy() * keVToMeV());
    auto end = upperIdx(data, maxPhotonEnergy() * keVToMeV());
    m_atom.incoherent.clear();
    m_atom.incoherent.reserve((end - start) / 2);
    for (std::size_t i = start; i <= end; i = i + 2) {
        const double e = data[i] * MeVTokeV();
#ifdef USE_NIST
        double a;
        try {
            a = xrlpp::CS_Compt(static_cast<int>(m_atom.Z), e);
        } catch (std::invalid_argument& e) {
            a = data[i + 1] * barnToAtt();
        }
#else
        const double a = data[i + 1] * barnToAtt();
#endif
        if (e > 0 && a > 0)
            m_atom.incoherent.push_back(std::make_pair(e, a));
    }
    m_atom.incoherent.shrink_to_fit();
}

void AtomicElementHandler::setIncoherentAvgEnergyScatteredPhoton(const std::vector<double>& data)
{
    auto start = lowerIdx(data, minPhotonEnergy() * keVToMeV());
    auto end = upperIdx(data, maxPhotonEnergy() * keVToMeV());
    m_atom.incoherentMeanScatterEnergy.clear();
    m_atom.incoherentMeanScatterEnergy.reserve((end - start) / 2);
    for (std::size_t i = start; i <= end; i = i + 2) {
        const double e = data[i] * MeVTokeV();
        const double a = data[i + 1] * MeVTokeV();
        if (e > 0 && a > 0)
            m_atom.incoherentMeanScatterEnergy.push_back(std::make_pair(e, a));
    }
    m_atom.incoherentMeanScatterEnergy.shrink_to_fit();
}

void AtomicElementHandler::setFormFactor(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.formFactor.clear();
    m_atom.formFactor.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6; // for units of per �
        const auto f = data[i + 1];
        if (x <= momentumTransferMax() && x > 0.0) {
            m_atom.formFactor.push_back(std::make_pair(x, f));
        }
    }
}
/*
void AtomicElementHandler::setImaginaryAnomalousSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_imagAnomSF.clear();
    m_imagAnomSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E3;
        const auto f = data[i + 1];
        m_imagAnomSF.push_back(std::make_pair(x, f));
    }
}

void AtomicElementHandler::setRealAnomalousSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_realAnomSF.clear();
    m_realAnomSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E3;
        const auto f = data[i + 1];
        m_realAnomSF.push_back(std::make_pair(x, f));
    }
}
*/
void AtomicElementHandler::setIncoherentSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.incoherentSF.clear();
    m_atom.incoherentSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6; // for units of per �
        const auto f = data[i + 1];
        if (x <= momentumTransferMax() && x > 0.0) {
            m_atom.incoherentSF.push_back(std::make_pair(x, f));
        }
    }
}

void AtomicElementHandler::setPhotoelectricData(const std::vector<double>& data)
{
    auto start = lowerIdx(data, minPhotonEnergy() * keVToMeV());
    auto end = upperIdx(data, maxPhotonEnergy() * keVToMeV());
    m_atom.photoel.clear();
    m_atom.photoel.reserve((end - start) / 2);
    for (std::size_t i = start; i <= end; i = i + 2) {
        const double e = data[i] * MeVTokeV();
#ifdef USE_NIST
        double a;
        try {
            a = xrlpp::CS_Photo(static_cast<int>(m_atom.Z), e);
        } catch (std::invalid_argument& e) {
            a = data[i + 1] * barnToAtt();
        }
#else
        const double a = data[i + 1] * barnToAtt();
#endif
        if (e > 0 && a > 0)
            m_atom.photoel.push_back(std::make_pair(e, a));
    }
    m_atom.photoel.shrink_to_fit();
}

void AtomicElementHandler::setShellPhotoelectricData(std::uint64_t shell, const std::vector<double>& data)
{
    if (!m_atom.shells.contains(shell)) {
        m_atom.shells[shell] = xraymc::AtomicShell(shell);
    }
    auto start = lowerIdx(data, minPhotonEnergy() * keVToMeV());
    auto end = upperIdx(data, maxPhotonEnergy() * keVToMeV());
    std::vector<std::pair<double, double>> photoel;
    photoel.reserve((end - start) / 2);

    for (std::size_t i = start; i <= end; i = i + 2) {
        const double e = data[i] * MeVTokeV();
        const double a = data[i + 1] * barnToAtt();
        if (e > 0 && a > 0)
            photoel.push_back(std::make_pair(e, a));
    }
    photoel.shrink_to_fit();
    m_atom.shells[shell].photoel = photoel;
}

void AtomicElementHandler::setRadiativeTransitionProbabilities(std::uint64_t shell, const std::vector<double>& data)
{
    if (!m_atom.shells.contains(shell)) {
        m_atom.shells[shell] = xraymc::AtomicShell(shell);
    }
    const auto N = data.size() / 3;
    std::vector<xraymc::AtomicShellRadiativeEmission> trans(N);
    for (std::size_t i = 0; i < N; ++i) {
        trans[i].vacancy = static_cast<std::uint64_t>(data[i * 3 + 0]);
        trans[i].probability = (data[i * 3 + 1]);
        trans[i].energy = (data[i * 3 + 2]) * MeVTokeV();
    }
    m_atom.shells.at(shell).radiativeEmissions = trans;
    return;
}

void AtomicElementHandler::setShellBindingEnergy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = xraymc::AtomicShell(shell);
        }
        m_atom.shells[shell].bindingEnergy = E;
    }
}

void AtomicElementHandler::setShellKineticEnergy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = xraymc::AtomicShell(shell);
        }
        m_atom.shells[shell].kineticEnergy = E;
    }
}

void AtomicElementHandler::setShellNumberOfElectrons(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = xraymc::AtomicShell(shell);
        }
        m_atom.shells[shell].numberOfElectrons = N;
    }
}

void AtomicElementHandler::setShellNumberOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = xraymc::AtomicShell(shell);
        }
        m_atom.shells[shell].numberOfPhotonsPerInitVacancy = N;
    }
}

void AtomicElementHandler::setShellEnergyOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = xraymc::AtomicShell(shell);
        }
        m_atom.shells[shell].energyOfPhotonsPerInitVacancy = E;
    }
}

bool AtomicElementHandler::operator==(const AtomicElementHandler& other) const
{
    bool eq = m_atom.Z == other.m_atom.Z;
    eq = eq && m_atom.atomicWeight == other.m_atom.atomicWeight;
    eq = eq && m_atom.photoel == other.m_atom.photoel;
    eq = eq && m_atom.coherent == other.m_atom.coherent;
    eq = eq && m_atom.incoherent == other.m_atom.incoherent;
    eq = eq && m_atom.incoherentSF == other.m_atom.incoherentSF;
    eq = eq && m_atom.formFactor == other.m_atom.formFactor;
    eq = eq && m_atom.standardDensity == other.m_atom.standardDensity;

    for (const auto& [key, shell] : m_atom.shells) {
        if (!other.m_atom.shells.contains(key))
            return false;
        const auto& lshell = m_atom.shells.at(key);
        const auto& rshell = other.m_atom.shells.at(key);
        eq = eq && lshell.shell == rshell.shell;
        eq = eq && lshell.bindingEnergy == rshell.bindingEnergy;
        eq = eq && lshell.kineticEnergy == rshell.kineticEnergy;
        eq = eq && lshell.energyOfPhotonsPerInitVacancy == rshell.energyOfPhotonsPerInitVacancy;
        eq = eq && lshell.HartreeFockOrbital_0 == rshell.HartreeFockOrbital_0;
        eq = eq && lshell.numberOfElectrons == rshell.numberOfElectrons;
        eq = eq && lshell.numberOfPhotonsPerInitVacancy == rshell.numberOfPhotonsPerInitVacancy;
        eq = eq && lshell.photoel == rshell.photoel;
    }
    return eq;
}
constexpr double AtomicElementHandler::maxPhotonEnergy()
{
    auto to_double = [](const std::string_view s) -> std::optional<double> {
        if (double value; std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
            return value;
        else
            return std::nullopt;
    };
    return std::min(to_double(XRAYMCLIB_MAXENERGY).value_or(500.0), 500.0);
}
constexpr double AtomicElementHandler::minPhotonEnergy()
{
    auto to_double = [](const std::string_view s) -> std::optional<double> {
        if (double value; std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
            return value;
        else
            return std::nullopt;
    };
    return std::min(to_double(XRAYMCLIB_MINTABLEENERGY).value_or(1.0), 1.0);
}

void AtomicElementHandler::setShellHartreeFockProfile_0(std::uint64_t shell, double J)
{
    if (!m_atom.shells.contains(shell)) {
        m_atom.shells[shell] = xraymc::AtomicShell(shell);
    }
    m_atom.shells[shell].HartreeFockOrbital_0 = J;
}
void AtomicElementHandler::setStandardDensity(double dens)
{
    m_atom.standardDensity = dens;
}
constexpr double AtomicElementHandler::momentumTransferMax()
{
    constexpr double hc_si = 1.239841193E-6; // ev*m
    constexpr double m2A = 1E10; // meters to �ngstr�m
    constexpr double eV2keV = 1E-3; // eV to keV
    constexpr double hc = hc_si * m2A * eV2keV; // kev*�
    constexpr double hc_inv = 1.0 / hc;
    return maxPhotonEnergy() * hc_inv; // per �
}
double AtomicElementHandler::momentumTransfer(double energy, double angle)
{
    constexpr double hc_si = 1.239841193E-6; // ev*m
    constexpr double m2A = 1E10; // meters to �ngstr�m
    constexpr double eV2keV = 1E-3; // eV to keV
    constexpr double hc = hc_si * m2A * eV2keV; // kev*�
    constexpr double hc_inv = 1.0 / hc;
    return energy * std::sin(angle * 0.5) * hc_inv; // per �
}
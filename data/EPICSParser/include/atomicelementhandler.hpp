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
#include "dxmc/material/atomicelement.hpp"

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

class AtomicElementHandler {
public:
    AtomicElementHandler() {};
    AtomicElementHandler(std::uint64_t Z);
    AtomicElementHandler(const dxmc::AtomicElement& atom)
        : m_atom(atom)
    {
    }

    bool operator==(const AtomicElementHandler& other) const;

    dxmc::AtomicElement atom() const { return m_atom; }
    void setAtom(const dxmc::AtomicElement& atom) { m_atom = atom; }

    void setZ(std::uint64_t Z) { m_atom.Z = Z; }
    void setAtomicWeight(double AW) { m_atom.atomicWeight = AW; }

    void setPhotoelectricData(const std::vector<double>& data);
    void setCoherentData(const std::vector<double>& data);
    void setIncoherentData(const std::vector<double>& data);
    void setIncoherentAvgEnergyScatteredPhoton(const std::vector<double>& data);

    void setFormFactor(const std::vector<double>& data);
    // void setImaginaryAnomalousSF(const std::vector<double>& data);
    // void setRealAnomalousSF(const std::vector<double>& data);
    void setIncoherentSF(const std::vector<double>& data);

    void setShellBindingEnergy(const std::vector<double>& data);
    void setShellKineticEnergy(const std::vector<double>& data);
    void setShellPhotoelectricData(const std::uint64_t shell, const std::vector<double>& data);
    void setShellNumberOfElectrons(const std::vector<double>& data);
    void setShellNumberOfPhotonsPerInitVacancy(const std::vector<double>& data);
    void setShellEnergyOfPhotonsPerInitVacancy(const std::vector<double>& data);
    void setRadiativeTransitionProbabilities(std::uint64_t shell, const std::vector<double>& data);

    void setShellHartreeFockProfile_0(std::uint64_t shell, double J);
    void setStandardDensity(double dens);

    auto& getShells() { return m_atom.shells; }

    static constexpr double momentumTransferMax();
    static double momentumTransfer(double energy, double angle);

    static constexpr double maxPhotonEnergy();
    static constexpr double minPhotonEnergy();
    static constexpr double MeVTokeV() { return 1000; }
    static constexpr double keVToMeV() { return 1.0 / MeVTokeV(); }
    double barnToAtt()
    {
        constexpr double u = 1.6605402;
        return 1.0 / (u * m_atom.atomicWeight);
    }

protected:
private:
    dxmc::AtomicElement m_atom;
};

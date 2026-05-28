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

Copyright 2026 Erlend Andersen
*/

#include "xraymc/beams/pencilbeam.hpp"

#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <std::size_t N = 15, int L = 2>
bool testfluencescore()
{
    using Sphere = xraymc::WorldSphere<N, L, false>;
    using FluenceScore = xraymc::FluenceScore;
    using World = xraymc::World<Sphere, FluenceScore>;
    using Material = xraymc::Material<N>;

    World world;
    world.reserveNumberOfItems(2);
    auto& sphere = world.template addItem<Sphere>();
    sphere.setRadius(1);

    auto& fscore = world.template addItem<FluenceScore>();
    fscore.setRadius(1);
    fscore.setCenter({ 0, 2, 0 });
    fscore.setPlaneNormal({ 0, -1, 0 });
    fscore.setEnergyStep(1.0);

    auto material = Material::byChemicalFormula("PbBiSi").value();

    const double density = 11.0;
    sphere.setMaterial(material, density);

    world.build();

    xraymc::PencilBeam<false> beam({ -2, 0, 0 }, { 1, 0, 0 }, 135);
    beam.setNumberOfExposures(50);
    beam.setNumberOfParticlesPerExposure(1E6);

    xraymc::Transport transport;
    transport.runConsole(world, beam);

    xraymc::VisualizeWorld viz(world);

    viz.addLineProp({ -2, 0, 0 }, { 1, 0, 0 }, 1, 0.1);
    int height = 528 * 2;
    int width = 528 * 2;
    std::vector<double> buffer(height * width * 4, 1);
    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(50);
        viz.setPolarAngleDeg(30);
        viz.setAzimuthalAngleDeg(i * 360.0 / 12.0);
        viz.suggestFOV(3);
        viz.generate(world, buffer, width, height);
        std::string name = "fluence_score" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer, width, height);
    }

    const auto specter = fscore.getFluenceSpecter();
    std::cout << "Energy keV, I (" << N << ")" << std::endl;
    for (const auto [energy, I] : specter) {
        if (energy < beam.energy() + 1)
            std::cout << energy << "," << I << std::endl;
    }
    return true;
}

struct Transmission {
    std::uint64_t from = 0;
    std::uint64_t to = 0;
    double energy = 0;
    double totalProbability = 0;
};

void radTransWorker(std::vector<Transmission>& trans, const xraymc::AtomicElement& atom, double probability, std::uint64_t shell)
{
    for (const auto& [shIdx, sh] : atom.shells) {
        if (sh.shell == shell) {
            for (const auto& rad_em : sh.radiativeEmissions) {
                Transmission t;
                t.from = shell;
                t.to = rad_em.vacancy;
                t.energy = rad_em.energy;
                t.totalProbability = probability * rad_em.probability;
                trans.push_back(t);

                // adding next cascade
                radTransWorker(trans, atom, rad_em.probability, rad_em.vacancy);
            }
        }
    }
}

template <int N = 32>
bool testRadiativeTransitionsApproximation(std::size_t Z = 53, double photonEnergy = 60, bool add = false)
{
    // iodine: 53
    const auto atom = xraymc::AtomHandler::Atom(Z);

    std::vector<Transmission> transmissions;

    std::map<std::uint64_t, double> shellpeatt;
    double cum_shellpeatt = 0;
    for (const auto& [shIdx, sh] : atom.shells) {
        if (sh.photoel.front().first < photonEnergy) {
            std::pair<double, double> pe_int = std::make_pair(photonEnergy, 0.0);
            auto f1 = std::upper_bound(sh.photoel.cbegin(), sh.photoel.cend(), pe_int, [](const auto& lh, const auto& rh) {
                return lh.first < rh.first;
            });
            if (f1 != sh.photoel.cend()) {
                // linear interpolation
                auto f0 = f1 - 1;
                double att = xraymc::interp(*f0, *f1, photonEnergy);
                shellpeatt[sh.shell] = att;
                cum_shellpeatt += att;
            }
        }
    }

    for (const auto& [shIdx, sh] : atom.shells) {
        if (shellpeatt.contains(sh.shell)) {
            radTransWorker(transmissions, atom, shellpeatt.at(sh.shell) / cum_shellpeatt, sh.shell);
        }
    }

    std::ofstream file2;

    if (add) {
        file2.open("shell_transmission.txt", std::ios::app);
        if (!file2.is_open())
            file2.open("shell_transmission.txt");
    } else {
        file2.open("shell_transmission.txt");
        file2 << "From,To,Energy,Probability,Z,IncomingEnergy\n";
    }

    for (const auto& t : transmissions) {
        file2 << t.from << ',' << t.to << ',' << t.energy << ',' << t.totalProbability;
        file2 << ',' << Z << ',' << photonEnergy << std::endl;
    }
    file2.close();

    std::map<double, double> enp;
    for (const auto& t : transmissions) {
        if (enp.contains(t.energy)) {
            enp[t.energy] += t.totalProbability;
        } else {
            enp[t.energy] = t.totalProbability;
        }
    }

    std::ofstream file;
    if (add) {
        file.open("transmissions.txt", std::ios::app);
        if (!file.is_open())
            file.open("transmissions.txt");
    } else {
        file.open("transmissions.txt");
    }

    if (!add) {
        file << "IncomingEnergy [keV],Energy [keV],Probability,Type,Z\n";
    }
    for (const auto& [key, value] : enp) {
        file << photonEnergy << ',';
        file << key << "," << value;
        file << "," << "Cascade," << Z << std::endl;
    }

    // xraymc probs
    constexpr std::size_t samples = 1E7;
    const auto mat = xraymc::Material<N>::byZ(Z).value();
    const double xraymcTotPEcross = mat.attenuationValues(photonEnergy).photoelectric;
    xraymc::RandomState state;
    std::map<double, double> xraymcprobs;

    for (std::size_t i = 0; i < samples; ++i) {
        xraymc::Particle p;
        p.energy = photonEnergy;
        p.weight = 1;
        double absorbed_energy = xraymc::interactions::photoelectricEffectIA(xraymcTotPEcross, p, mat, state);
        if (p.energy > 0) {
            if (xraymcprobs.contains(p.energy)) {
                xraymcprobs[p.energy] += p.weight;
            } else {
                xraymcprobs[p.energy] = p.weight;
            }
        }
    }
    for (const auto& [key, value] : xraymcprobs) {
        file << photonEnergy << ',';
        file << key << "," << value / samples;
        file << "," << "XrayMCsampled," << Z << std::endl;
    }

    for (const auto& [shIdx, sh] : atom.shells) {
        if (sh.energyOfPhotonsPerInitVacancy > 0) {
            file << photonEnergy << ',';
            file << sh.energyOfPhotonsPerInitVacancy << ',' << sh.numberOfPhotonsPerInitVacancy * shellpeatt[sh.shell] / cum_shellpeatt;
            file << "," << "XrayMCmodel," << Z << std::endl;
        }
    }
    file.close();
    return true;
}

int main()
{
    std::cout << "Testing high Z material scatter\n";
    bool success = true;

    testRadiativeTransitionsApproximation<31>(82, 100);
    testRadiativeTransitionsApproximation<31>(53, 100, true);

    return 0;
    success = success && testfluencescore<63, 2>();
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
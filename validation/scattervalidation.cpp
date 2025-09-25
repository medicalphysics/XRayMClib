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

Copyright 2024 Erlend Andersen
*/

#include "xraymc/interactions.hpp"

#include <algorithm>
#include <atomic>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>

constexpr int NSHELLS = 36;
constexpr std::size_t N = 5E7;

struct Histogram {
    double start = 0;
    double stop = 1;
    int N = 100;
    std::vector<std::uint64_t> intensity;
    Histogram(int n = 100, double start_val = 0, double stop_val = 1)
    {
        N = std::clamp(n, 2, 1000);
        intensity.resize(N);
        start = std::min(start_val, stop_val);
        stop = std::max(start_val, stop_val);
    }
    void operator()(double v)
    {
        double fidx = ((v - start) / (stop - start)) * (N - 1);
        if (fidx > 0 && fidx < N - 1) {
            auto idx = static_cast<int>(fidx);
            intensity[idx]++;
        }
    }

    std::vector<double> getBins(bool midpoint = true) const
    {
        if (midpoint) {
            std::vector<double> bins(N);
            const auto step = (stop - start) / N;
            for (std::size_t i = 0; i < N; ++i) {
                bins[i] = start + (i + 0.5) * step;
            }
            return bins;
        } else {
            std::vector<double> bins(N + 1);
            const auto step = (stop - start) / N;
            for (std::size_t i = 0; i < N + 1; ++i) {
                bins[i] = start + i * step;
            }
            return bins;
        }
    }
    std::vector<double> getValues() const
    {
        std::vector<double> vals(N, 0.0);
        const auto sum = static_cast<double>(std::accumulate(intensity.cbegin(), intensity.cend(), std::size_t { 0 }));
        if (sum > 0)
            std::transform(intensity.cbegin(), intensity.cend(), vals.begin(), [=](auto i) { return i / sum; });
        return vals;
    }
    void save(const std::string& fname) const
    {
        std::ofstream myfile;
        myfile.open(fname, std::ios::out);
        auto bins = getBins();
        auto vals = getValues();
        for (std::size_t i = 0; i < vals.size(); ++i) {
            myfile << bins[i] << "," << vals[i] << '\n';
        }
        myfile.close();
    }
};

class ResultPrint {
public:
    ResultPrint()
    {
        m_myfile.open("validationScatterTable.txt", std::ios::trunc);
    }
    ResultPrint(const std::string& fname)
    {
        m_myfile.open(fname, std::ios::out);
    }

    ~ResultPrint()
    {
        m_myfile.close();
    }
    void header()
    {
        m_myfile << "Model;InteractionType;x;y;Energy;Material\n";
    }
    void operator()(const std::string& model, const std::string& type, double x, double y, double energy, const std::string& matname = "")
    {
        const std::lock_guard<std::mutex> lock(m_mutex);
        m_myfile << model << ";";
        m_myfile << type << ";";
        m_myfile << x << ";";
        m_myfile << y << ";";
        m_myfile << energy << ";";
        m_myfile << matname << '\n';
    }
    void operator()(const Histogram& hist, const std::string& model, const std::string& type, double energy, const std::string& matname = "")
    {
        const auto b = hist.getBins();
        const auto v = hist.getValues();
        const std::lock_guard<std::mutex> lock(m_mutex);
        for (std::size_t i = 0; i < v.size(); ++i) {
            m_myfile << model << ';';
            m_myfile << type << ';';
            m_myfile << b[i] << ';';
            m_myfile << v[i] << ';';
            m_myfile << energy << ';';
            m_myfile << matname << '\n';
        }
    }

private:
    std::ofstream m_myfile;
    std::mutex m_mutex;
};

template <int Model = 1>
Histogram comptonScatterAngle(const xraymc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    xraymc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        xraymc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        xraymc::interactions::comptonScatter<NSHELLS, Model>(p, material, state);
        const auto angle = xraymc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * xraymc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int Model = 1>
Histogram comptonScatterEnergy(const xraymc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    xraymc::RandomState state;
    Histogram h(1000, .8 / (1 + energy / xraymc::ELECTRON_REST_MASS() * 2), 1);
    for (std::size_t i = 0; i < N; ++i) {
        xraymc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        xraymc::interactions::comptonScatter<NSHELLS, Model>(p, material, state);
        h(p.energy / energy);
    }
    return h;
}

template <int Model = 1>
Histogram photoElectricEnergyIA(const xraymc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    xraymc::RandomState state;
    double max_en = 0;
    for (const auto& shell : material.shells()) {
        max_en = std::max(max_en, shell.energyOfPhotonsPerInitVacancy);
    }
    Histogram h(1000, 0, max_en + 0.1);

    const auto att = material.attenuationPhotoeletric(energy);

    for (std::size_t i = 0; i < N; ++i) {
        xraymc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        const auto E = xraymc::interactions::photoelectricEffect<NSHELLS, Model>(att, p, material, state);
        if (p.energy > 0)
            h(p.energy);
    }
    return h;
}

template <int Model = 1>
Histogram rayleightScatterAngle(const xraymc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    xraymc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        xraymc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        xraymc::interactions::rayleightScatter<NSHELLS, Model>(p, material, state);
        const auto angle = xraymc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * xraymc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int I = 0>
void saveHist(ResultPrint& p, const xraymc::Material<NSHELLS>& material, double energy, const std::string& matname, std::size_t N = 1E6)
{
    auto h_en = comptonScatterEnergy<I>(material, energy, N);
    auto h_ang = comptonScatterAngle<I>(material, energy, N);
    auto hr_ang = rayleightScatterAngle<I>(material, energy, N);
    std::string model = "NoneLC";
    if (I == 1)
        model = "Livermore";
    if (I == 2)
        model = "IA";

    p(h_ang, model, "ComptonAngle", energy, matname);
    p(h_en, model, "ComptonEnergy", energy, matname);
    p(hr_ang, model, "RayleighAngle", energy, matname);
    if constexpr (I == 2) {
        auto h_photo = photoElectricEnergyIA<I>(material, energy, N);
        p(h_photo, model, "PhotoElectricIA", energy, matname);
    }
}

auto TG195_breast_tissue()
{
    std::map<std::size_t, double> adipose_w;
    adipose_w[1] = 11.2;
    adipose_w[6] = 61.9;
    adipose_w[7] = 1.7;
    adipose_w[8] = 25.1;
    adipose_w[15] = 0.025;
    adipose_w[16] = 0.025;
    adipose_w[19] = 0.025;
    adipose_w[20] = 0.025;

    const double adipose_d = 0.93;

    std::map<std::size_t, double> gland_w;
    gland_w[1] = 10.2;
    gland_w[6] = 18.4;
    gland_w[7] = 3.2;
    gland_w[8] = 67.7;
    gland_w[15] = 0.125;
    gland_w[16] = 0.125;
    gland_w[19] = 0.125;
    gland_w[20] = 0.125;

    const double gland_d = 1.04;

    // weighted 20% gland 80% adipose
    std::map<std::size_t, double> w;
    for (const auto [Z, n] : adipose_w) {
        if (!w.contains(Z))
            w[Z] = 0;
        w[Z] += n * 0.8;
    }
    for (const auto [Z, n] : gland_w) {
        if (!w.contains(Z))
            w[Z] = 0;
        w[Z] += n * 0.2;
    }
    const auto d = adipose_d * 0.8 + gland_d * 0.2;

    return xraymc::Material<NSHELLS>::byWeight(w).value();
}

int main()
{

    ResultPrint p;
    p.header();

    std::vector<double> energies;
    // energies.push_back(15);
    energies.push_back(16.8);
    energies.push_back(30);
    //energies.push_back(50);
    energies.push_back(120);

    std::vector<std::string> material_names;
    //material_names.push_back("Water, Liquid");
    //material_names.push_back("Polymethyl Methacralate (Lucite, Perspex)");
    material_names.push_back("Gold");
    material_names.push_back("TG195Breast");

    std::vector<xraymc::Material<NSHELLS>> materials;
    //materials.push_back(xraymc::Material<NSHELLS>::byNistName("Water, Liquid").value());
    //materials.push_back(xraymc::Material<NSHELLS>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value());
    materials.push_back(xraymc::Material<NSHELLS>::byZ(79).value());
    materials.push_back(TG195_breast_tissue());

    std::vector<std::jthread> threads;
    threads.reserve(materials.size() * energies.size());

    for (std::size_t i = 0; i < materials.size(); ++i) {
        const auto& material_name = material_names[i];
        const auto& material = materials[i];
        for (auto energy : energies) {
            threads.emplace_back(saveHist<0>, std::ref(p), std::cref(material), energy, std::cref(material_name), N);
            threads.emplace_back(saveHist<1>, std::ref(p), std::cref(material), energy, std::cref(material_name), N);
            threads.emplace_back(saveHist<2>, std::ref(p), std::cref(material), energy, std::cref(material_name), N);
        }
    }

    return 0;
}
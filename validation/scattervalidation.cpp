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

Copyright 2024 Erlend Andersen
*/

#include "dxmc/interactions.hpp"

#include <algorithm>
#include <atomic>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>

constexpr int NSHELLS = 15;

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
        std::vector<double> vals(N);
        const auto sum = static_cast<double>(std::accumulate(intensity.cbegin(), intensity.cend(), std::size_t { 0 }));
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
Histogram comptonScatterAngle(const dxmc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::comptonScatter<NSHELLS, Model>(p, material, state);
        const auto angle = dxmc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * dxmc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int Model = 1>
Histogram comptonScatterEnergy(const dxmc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(1000, .9 / (1 + energy / dxmc::ELECTRON_REST_MASS() * 2), 1);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::comptonScatter<NSHELLS, Model>(p, material, state);
        h(p.energy / energy);
    }
    return h;
}

template <int Model = 1>
Histogram rayleightScatterAngle(const dxmc::Material<NSHELLS>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::rayleightScatter<NSHELLS, Model>(p, material, state);
        const auto angle = dxmc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * dxmc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int I = 0>
void saveHist(ResultPrint& p, const dxmc::Material<NSHELLS>& material, double energy, const std::string& matname, std::size_t N = 1E6)
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
}

int main()
{

    constexpr std::size_t N = 5E6;

    ResultPrint p;
    p.header();

    std::vector<double> energies;
    // energies.push_back(15);
    // energies.push_back(30);
    energies.push_back(50);

    std::vector<std::string> material_names;
    // material_names.push_back("Water, Liquid");
    // material_names.push_back("Polymethyl Methacralate (Lucite, Perspex)");
    material_names.push_back("Gold");

    std::vector<dxmc::Material<NSHELLS>> materials;
    // materials.push_back(dxmc::Material<NSHELLS>::byNistName("Water, Liquid").value());
    // materials.push_back(dxmc::Material<NSHELLS>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value());
    materials.push_back(dxmc::Material<NSHELLS>::byZ(79).value());

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
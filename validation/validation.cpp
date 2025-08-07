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

Copyright 2023 Erlend Andersen
*/

#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/isotropicbeamcircle.hpp"
#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/beams/isotropicmonoenergybeamcircle.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/fluencescore.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldboxgrid.hpp"
#include "dxmc/world/worlditems/worldcylinder.hpp"

#include "tg195world3breast.hpp"
#include "tg195world42.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

using namespace dxmc;

// Set this to true for a reduced number of photons (for testing)
constexpr bool SAMPLE_RUN = false;
constexpr std::size_t NShells = 12;
constexpr double Sigma = 1.96;

struct ResultKeys {
    std::string rCase = "unknown";
    std::string volume = "unknown";
    std::string specter = "unknown";
    std::string modus = "unknown";
    std::string model = "unknown";
    double TG195Result = 0;
    double result = 0;
    double result_std = 0;
    std::uint64_t nEvents = 0;
    std::uint64_t nMilliseconds = 0;
};

class ResultPrint {
private:
    std::ofstream m_myfile;
    static std::string m_filename;

public:
    ResultPrint()
    {
        if (m_filename.empty())
            m_filename = "validationTable.txt";
        m_myfile.open(m_filename, std::ios::out | std::ios::app);
    }
    ResultPrint(std::string_view fname)
    {
        if (fname.empty())
            m_filename = "validationTable.txt";
        else
            m_filename = std::string(fname);
        m_myfile.open(m_filename, std::ios::out | std::ios::app);
    }
    ~ResultPrint()
    {
        m_myfile.close();
    }
    void header(bool print_only_if_empty = true)
    {
        if (print_only_if_empty) {
            m_myfile.close();
            std::ifstream file(m_filename, std::ios::binary | std::ios::ate);
            bool empty = file.tellg() == 0;
            file.close();
            m_myfile.open(m_filename, std::ios::out | std::ios::app);
            if (!empty)
                return;
        }
        m_myfile << "Case,Volume,Specter,Model,Mode,Result,Uncertainty,nEvents,SimulationTime\n";
    }

    void operator()(const ResultKeys& r, bool terminal = true, const std::string& units = "ev/hist")
    {
        print(r, terminal, units);
    }

    void print(const ResultKeys& r, bool terminal = true, const std::string& units = "ev/hist")
    {
        m_myfile << r.rCase << ",";
        m_myfile << r.volume << ",";
        m_myfile << r.specter << ",";
        m_myfile << r.model << ",";
        m_myfile << r.modus << ",";
        m_myfile << r.result << ",";
        m_myfile << r.result_std << ",";
        m_myfile << r.nEvents << ",";
        m_myfile << r.nMilliseconds;
        m_myfile << std::endl;

        if (terminal) {
            std::cout << "VOI: " << r.volume;
            std::cout << " " << units << ": " << r.result;
            std::cout << " TG195: " << r.TG195Result;
            std::cout << " difference[%]: " << (r.result / r.TG195Result - 1) * 100 << std::endl;
        }
    }
};
// static initialization
std::string ResultPrint::m_filename;

template <typename W, typename B>
void saveImageOfWorld(const std::string& name, W& world, B& beam, double polarAngle = 90, double azimuthAngle = 90, double dist = 100, double zoom = 1, double linelenght = 250, double linethick = 0.2)
{
    dxmc::VisualizeWorld viz(world);
    viz.setAzimuthalAngleDeg(azimuthAngle);
    viz.setPolarAngleDeg(polarAngle);
    viz.setDistance(dist);
    viz.suggestFOV(zoom);
    auto buffer = viz.template createBuffer<double>(2048 * 2, 2048 * 2);
    if constexpr (std::is_same_v<B, dxmc::IsotropicMonoEnergyBeamCircle<false>> || std::is_same_v<B, dxmc::IsotropicBeamCircle<false>>) {
        const auto& p = beam.position();
        const auto& d = beam.direction();
        const auto a = beam.collimationHalfAngle();
        const auto x = [&d, a]() {
            const auto ind = dxmc::vectormath::argmin3<std::size_t>(d);
            std::array<double, 3> nr = { 0, 0, 0 };
            nr[ind] = 1;
            return dxmc::vectormath::scale(dxmc::vectormath::normalized(dxmc::vectormath::cross(d, nr)), std::atan(a));
        }();
        for (int ang = 0; ang < 360; ang = ang + 30) {
            const auto xa = dxmc::vectormath::rotate(x, d, static_cast<double>(ang));
            const auto dir = dxmc::vectormath::normalized(dxmc::vectormath::add(d, xa));
            viz.addLineProp(p, dir, linelenght, linethick);
        }
    } else {
        viz.addLineProp(beam, linelenght, linethick);
    }

    std::array<std::uint8_t, 3> color;
    color = { 255, 195, 170 }; // skin
    viz.setColorOfItem(world.getItemPointerFromName("Tissue"), color);

    viz.generate(world, buffer);
    viz.savePNG(name, buffer);
    return;
}

// energy weighs pair for spectre
/*RQR-8
W/Al
100 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.708 mm
Mean Energy: 50.6 keV
HVL: 3.950 mm Al
QVL: 9.840 mm Al
*/
const std::vector<double> TG195_100KV_raw({ 16.25, 1.423E-04, 16.75, 2.157E-04, 17.25, 3.102E-04, 17.75, 4.324E-04, 18.25, 5.840E-04, 18.75, 7.644E-04, 19.25, 9.784E-04, 19.75, 1.222E-03, 20.25, 1.491E-03, 20.75, 1.803E-03, 21.25, 2.129E-03, 21.75, 2.490E-03, 22.25, 2.863E-03, 22.75, 3.263E-03, 23.25, 3.658E-03, 23.75, 4.093E-03, 24.25, 4.504E-03, 24.75, 4.912E-03, 25.25, 5.347E-03, 25.75, 5.769E-03, 26.25, 6.168E-03, 26.75, 6.582E-03, 27.25, 6.965E-03, 27.75, 7.360E-03, 28.25, 7.710E-03, 28.75, 8.067E-03, 29.25, 8.368E-03, 29.75, 8.671E-03, 30.25, 8.975E-03, 30.75, 9.213E-03, 31.25, 9.476E-03, 31.75, 9.694E-03, 32.25, 9.903E-03, 32.75, 1.009E-02, 33.25, 1.025E-02, 33.75, 1.040E-02, 34.25, 1.053E-02, 34.75, 1.063E-02, 35.25, 1.073E-02, 35.75, 1.081E-02, 36.25, 1.087E-02, 36.75, 1.092E-02, 37.25, 1.096E-02, 37.75, 1.099E-02, 38.25, 1.100E-02, 38.75, 1.100E-02, 39.25, 1.099E-02, 39.75, 1.098E-02, 40.25, 1.095E-02, 40.75, 1.091E-02, 41.25, 1.086E-02, 41.75, 1.081E-02, 42.25, 1.076E-02, 42.75, 1.069E-02, 43.25, 1.063E-02, 43.75, 1.055E-02, 44.25, 1.048E-02, 44.75, 1.039E-02, 45.25, 1.031E-02, 45.75, 1.022E-02, 46.25, 1.012E-02, 46.75, 1.003E-02, 47.25, 9.933E-03, 47.75, 9.828E-03, 48.25, 9.732E-03, 48.75, 9.628E-03, 49.25, 9.516E-03, 49.75, 9.412E-03, 50.25, 9.302E-03, 50.75, 9.193E-03, 51.25, 9.084E-03, 51.75, 8.970E-03, 52.25, 8.862E-03, 52.75, 8.749E-03, 53.25, 8.637E-03, 53.75, 8.526E-03, 54.25, 8.409E-03, 54.75, 8.300E-03, 55.25, 8.185E-03, 55.75, 8.072E-03, 56.25, 7.959E-03, 56.75, 7.847E-03, 57.25, 7.737E-03, 57.75, 2.568E-02, 58.25, 7.513E-03, 58.75, 7.405E-03, 59.25, 3.920E-02, 59.75, 7.181E-03, 60.25, 7.071E-03, 60.75, 6.962E-03, 61.25, 6.854E-03, 61.75, 6.746E-03, 62.25, 6.640E-03, 62.75, 6.530E-03, 63.25, 6.425E-03, 63.75, 6.321E-03, 64.25, 6.214E-03, 64.75, 6.107E-03, 65.25, 6.006E-03, 65.75, 5.901E-03, 66.25, 5.797E-03, 66.75, 1.673E-02, 67.25, 5.592E-03, 67.75, 5.491E-03, 68.25, 5.390E-03, 68.75, 8.223E-03, 69.25, 5.055E-03, 69.75, 4.296E-03, 70.25, 4.236E-03, 70.75, 4.171E-03, 71.25, 4.110E-03, 71.75, 4.048E-03, 72.25, 3.982E-03, 72.75, 3.919E-03, 73.25, 3.852E-03, 73.75, 3.787E-03, 74.25, 3.719E-03, 74.75, 3.654E-03, 75.25, 3.585E-03, 75.75, 3.516E-03, 76.25, 3.449E-03, 76.75, 3.379E-03, 77.25, 3.308E-03, 77.75, 3.240E-03, 78.25, 3.169E-03, 78.75, 3.098E-03, 79.25, 3.026E-03, 79.75, 2.954E-03, 80.25, 2.882E-03, 80.75, 2.809E-03, 81.25, 2.736E-03, 81.75, 2.665E-03, 82.25, 2.592E-03, 82.75, 2.519E-03, 83.25, 2.445E-03, 83.75, 2.370E-03, 84.25, 2.296E-03, 84.75, 2.222E-03, 85.25, 2.148E-03, 85.75, 2.073E-03, 86.25, 1.999E-03, 86.75, 1.925E-03, 87.25, 1.850E-03, 87.75, 1.776E-03, 88.25, 1.700E-03, 88.75, 1.625E-03, 89.25, 1.550E-03, 89.75, 1.476E-03, 90.25, 1.400E-03, 90.75, 1.326E-03, 91.25, 1.251E-03, 91.75, 1.177E-03, 92.25, 1.101E-03, 92.75, 1.027E-03, 93.25, 9.529E-04, 93.75, 8.781E-04, 94.25, 8.041E-04, 94.75, 7.302E-04, 95.25, 6.559E-04, 95.75, 5.823E-04, 96.25, 5.089E-04, 96.75, 4.353E-04, 97.25, 3.623E-04, 97.75, 2.892E-04, 98.25, 2.166E-04, 98.75, 1.441E-04, 99.25, 7.193E-05, 99.75, 5.990E-06 });

/*RQR-9
W/Al
120 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.861 mm
Mean Energy: 56.4 keV
HVL: 5.010 mm Al
*/
const std::vector<double> TG195_120KV_raw({ 16.75, 1.107E-04, 17.25, 1.625E-04, 17.75, 2.308E-04, 18.25, 3.172E-04, 18.75, 4.220E-04, 19.25, 5.486E-04, 19.75, 6.956E-04, 20.25, 8.610E-04, 20.75, 1.056E-03, 21.25, 1.264E-03, 21.75, 1.499E-03, 22.25, 1.748E-03, 22.75, 2.019E-03, 23.25, 2.293E-03, 23.75, 2.601E-03, 24.25, 2.900E-03, 24.75, 3.203E-03, 25.25, 3.531E-03, 25.75, 3.858E-03, 26.25, 4.176E-03, 26.75, 4.511E-03, 27.25, 4.830E-03, 27.75, 5.163E-03, 28.25, 5.469E-03, 28.75, 5.786E-03, 29.25, 6.065E-03, 29.75, 6.349E-03, 30.25, 6.638E-03, 30.75, 6.879E-03, 31.25, 7.143E-03, 31.75, 7.372E-03, 32.25, 7.597E-03, 32.75, 7.804E-03, 33.25, 7.994E-03, 33.75, 8.171E-03, 34.25, 8.339E-03, 34.75, 8.483E-03, 35.25, 8.622E-03, 35.75, 8.745E-03, 36.25, 8.849E-03, 36.75, 8.949E-03, 37.25, 9.031E-03, 37.75, 9.109E-03, 38.25, 9.170E-03, 38.75, 9.219E-03, 39.25, 9.264E-03, 39.75, 9.297E-03, 40.25, 9.319E-03, 40.75, 9.332E-03, 41.25, 9.333E-03, 41.75, 9.332E-03, 42.25, 9.327E-03, 42.75, 9.307E-03, 43.25, 9.292E-03, 43.75, 9.259E-03, 44.25, 9.229E-03, 44.75, 9.187E-03, 45.25, 9.149E-03, 45.75, 9.101E-03, 46.25, 9.044E-03, 46.75, 8.996E-03, 47.25, 8.937E-03, 47.75, 8.871E-03, 48.25, 8.813E-03, 48.75, 8.747E-03, 49.25, 8.672E-03, 49.75, 8.605E-03, 50.25, 8.530E-03, 50.75, 8.456E-03, 51.25, 8.381E-03, 51.75, 8.300E-03, 52.25, 8.226E-03, 52.75, 8.145E-03, 53.25, 8.065E-03, 53.75, 7.985E-03, 54.25, 7.899E-03, 54.75, 7.820E-03, 55.25, 7.736E-03, 55.75, 7.652E-03, 56.25, 7.568E-03, 56.75, 7.486E-03, 57.25, 7.403E-03, 57.75, 3.335E-02, 58.25, 7.236E-03, 58.75, 7.155E-03, 59.25, 5.339E-02, 59.75, 6.986E-03, 60.25, 6.903E-03, 60.75, 6.821E-03, 61.25, 6.739E-03, 61.75, 6.658E-03, 62.25, 6.578E-03, 62.75, 6.494E-03, 63.25, 6.415E-03, 63.75, 6.338E-03, 64.25, 6.256E-03, 64.75, 6.175E-03, 65.25, 6.100E-03, 65.75, 6.021E-03, 66.25, 5.942E-03, 66.75, 2.242E-02, 67.25, 5.788E-03, 67.75, 5.712E-03, 68.25, 5.637E-03, 68.75, 9.988E-03, 69.25, 5.257E-03, 69.75, 4.045E-03, 70.25, 4.019E-03, 70.75, 3.988E-03, 71.25, 3.960E-03, 71.75, 3.932E-03, 72.25, 3.900E-03, 72.75, 3.871E-03, 73.25, 3.838E-03, 73.75, 3.808E-03, 74.25, 3.774E-03, 74.75, 3.743E-03, 75.25, 3.709E-03, 75.75, 3.674E-03, 76.25, 3.641E-03, 76.75, 3.606E-03, 77.25, 3.570E-03, 77.75, 3.537E-03, 78.25, 3.500E-03, 78.75, 3.463E-03, 79.25, 3.426E-03, 79.75, 3.389E-03, 80.25, 3.351E-03, 80.75, 3.313E-03, 81.25, 3.274E-03, 81.75, 3.238E-03, 82.25, 3.200E-03, 82.75, 3.160E-03, 83.25, 3.121E-03, 83.75, 3.079E-03, 84.25, 3.039E-03, 84.75, 3.000E-03, 85.25, 2.959E-03, 85.75, 2.919E-03, 86.25, 2.878E-03, 86.75, 2.838E-03, 87.25, 2.797E-03, 87.75, 2.756E-03, 88.25, 2.712E-03, 88.75, 2.671E-03, 89.25, 2.629E-03, 89.75, 2.588E-03, 90.25, 2.544E-03, 90.75, 2.502E-03, 91.25, 2.460E-03, 91.75, 2.418E-03, 92.25, 2.374E-03, 92.75, 2.331E-03, 93.25, 2.289E-03, 93.75, 2.244E-03, 94.25, 2.202E-03, 94.75, 2.159E-03, 95.25, 2.115E-03, 95.75, 2.072E-03, 96.25, 2.029E-03, 96.75, 1.984E-03, 97.25, 1.941E-03, 97.75, 1.896E-03, 98.25, 1.853E-03, 98.75, 1.809E-03, 99.25, 1.765E-03, 99.75, 1.722E-03, 100.25, 1.677E-03, 100.75, 1.634E-03, 101.25, 1.589E-03, 101.75, 1.546E-03, 102.25, 1.501E-03, 102.75, 1.458E-03, 103.25, 1.414E-03, 103.75, 1.370E-03, 104.25, 1.326E-03, 104.75, 1.282E-03, 105.25, 1.238E-03, 105.75, 1.195E-03, 106.25, 1.151E-03, 106.75, 1.107E-03, 107.25, 1.063E-03, 107.75, 1.019E-03, 108.25, 9.761E-04, 108.75, 9.323E-04, 109.25, 8.893E-04, 109.75, 8.456E-04, 110.25, 8.027E-04, 110.75, 7.592E-04, 111.25, 7.158E-04, 111.75, 6.731E-04, 112.25, 6.300E-04, 112.75, 5.874E-04, 113.25, 5.445E-04, 113.75, 5.017E-04, 114.25, 4.594E-04, 114.75, 4.168E-04, 115.25, 3.747E-04, 115.75, 3.324E-04, 116.25, 2.903E-04, 116.75, 2.485E-04, 117.25, 2.067E-04, 117.75, 1.650E-04, 118.25, 1.236E-04, 118.75, 8.222E-05, 119.25, 4.102E-05, 119.75, 3.417E-06 });

/*
RQR-M3
Mo/Mo
30 kVp
15 deg anode angle
0% Ripple
Mo filter thickness: 0.0386 mm
Mean Energy: 16.8 keV
HVL: 0.3431 mm Al
QVL: 0.7663 mm Al
*/
const std::vector<double> TG195_30KV_raw({ 7.25, 1.551E-04, 7.75, 4.691E-04, 8.25, 1.199E-03, 8.75, 2.405E-03, 9.25, 4.263E-03, 9.75, 6.797E-03, 10.25, 9.761E-03, 10.75, 1.314E-02, 11.25, 1.666E-02, 11.75, 2.013E-02, 12.25, 2.349E-02, 12.75, 2.666E-02, 13.25, 2.933E-02, 13.75, 3.167E-02, 14.25, 3.365E-02, 14.75, 3.534E-02, 15.25, 3.644E-02, 15.75, 3.741E-02, 16.25, 3.796E-02, 16.75, 3.823E-02, 17.25, 3.445E-01, 17.75, 3.770E-02, 18.25, 3.704E-02, 18.75, 3.639E-02, 19.25, 9.200E-02, 19.75, 2.178E-03, 20.25, 2.048E-03, 20.75, 2.043E-03, 21.25, 2.098E-03, 21.75, 2.193E-03, 22.25, 2.327E-03, 22.75, 2.471E-03, 23.25, 2.625E-03, 23.75, 2.770E-03, 24.25, 2.907E-03, 24.75, 3.000E-03, 25.25, 3.062E-03, 25.75, 3.058E-03, 26.25, 2.988E-03, 26.75, 2.823E-03, 27.25, 2.575E-03, 27.75, 2.233E-03, 28.25, 1.815E-03, 28.75, 1.290E-03, 29.25, 6.696E-04, 29.75, 4.086E-05 });

std::vector<std::pair<double, double>> TG195_specter(const std::vector<double>& raw)
{
    std::vector<std::pair<double, double>> s;
    s.reserve(raw.size() / 2);
    for (std::size_t i = 0; i < raw.size(); i = i + 2) {
        s.push_back(std::make_pair(static_cast<double>(raw[i] - 0.25), static_cast<double>(raw[i + 1])));
    }
    return s;
}

auto TG195_120KV()
{
    return TG195_specter(TG195_120KV_raw);
}

auto TG195_100KV()
{
    return TG195_specter(TG195_100KV_raw);
}

auto TG195_30KV()
{
    return TG195_specter(TG195_30KV_raw);
}

// mass energy abs coeff for air
std::vector<std::pair<double, double>> TG195_mass_en_abs_air()
{
    static std::vector<double> raw({ 0.25, 0.000E+00, 0.75, 0.000E+00, 1.25, 0.000E+00, 1.75, 0.000E+00, 2.25, 0.000E+00, 2.75, 0.000E+00, 3.25, 0.000E+00, 3.75, 3.815E+01, 4.25, 6.506E+01, 4.75, 4.657E+01, 5.25, 3.441E+01, 5.75, 2.609E+01, 6.25, 2.027E+01, 6.75, 1.606E+01, 7.25, 1.292E+01, 7.75, 1.050E+01, 8.25, 8.649E+00, 8.75, 7.219E+00, 9.25, 6.073E+00, 9.75, 5.154E+00, 10.25, 4.411E+00, 10.75, 3.804E+00, 11.25, 3.301E+00, 11.75, 2.880E+00, 12.25, 2.528E+00, 12.75, 2.230E+00, 13.25, 1.976E+00, 13.75, 1.759E+00, 14.25, 1.572E+00, 14.75, 1.410E+00, 15.25, 1.270E+00, 15.75, 1.146E+00, 16.25, 1.038E+00, 16.75, 9.434E-01, 17.25, 8.599E-01, 17.75, 7.864E-01, 18.25, 7.204E-01, 18.75, 6.611E-01, 19.25, 6.083E-01, 19.75, 5.612E-01, 20.25, 5.211E-01, 20.75, 4.845E-01, 21.25, 4.496E-01, 21.75, 4.177E-01, 22.25, 3.901E-01, 22.75, 3.650E-01, 23.25, 3.411E-01, 23.75, 3.181E-01, 24.25, 2.974E-01, 24.75, 2.788E-01, 25.25, 2.621E-01, 25.75, 2.469E-01, 26.25, 2.325E-01, 26.75, 2.193E-01, 27.25, 2.070E-01, 27.75, 1.956E-01, 28.25, 1.851E-01, 28.75, 1.753E-01, 29.25, 1.663E-01, 29.75, 1.579E-01, 30.25, 1.501E-01, 30.75, 1.428E-01, 31.25, 1.360E-01, 31.75, 1.298E-01, 32.25, 1.240E-01, 32.75, 1.185E-01, 33.25, 1.133E-01, 33.75, 1.086E-01, 34.25, 1.041E-01, 34.75, 9.979E-02, 35.25, 9.581E-02, 35.75, 9.211E-02, 36.25, 8.868E-02, 36.75, 8.544E-02, 37.25, 8.236E-02, 37.75, 7.945E-02, 38.25, 7.670E-02, 38.75, 7.410E-02, 39.25, 7.165E-02, 39.75, 6.933E-02, 40.25, 6.715E-02, 40.75, 6.510E-02, 41.25, 6.310E-02, 41.75, 6.122E-02, 42.25, 5.947E-02, 42.75, 5.782E-02, 43.25, 5.629E-02, 43.75, 5.479E-02, 44.25, 5.332E-02, 44.75, 5.196E-02, 45.25, 5.069E-02, 45.75, 4.946E-02, 46.25, 4.832E-02, 46.75, 4.720E-02, 47.25, 4.610E-02, 47.75, 4.510E-02, 48.25, 4.411E-02, 48.75, 4.321E-02, 49.25, 4.233E-02, 49.75, 4.147E-02, 50.25, 4.068E-02, 50.75, 3.991E-02, 51.25, 3.916E-02, 51.75, 3.848E-02, 52.25, 3.781E-02, 52.75, 3.715E-02, 53.25, 3.657E-02, 53.75, 3.594E-02, 54.25, 3.532E-02, 54.75, 3.483E-02, 55.25, 3.434E-02, 55.75, 3.381E-02, 56.25, 3.334E-02, 56.75, 3.289E-02, 57.25, 3.243E-02, 57.75, 3.205E-02, 58.25, 3.166E-02, 58.75, 3.124E-02, 59.25, 3.081E-02, 59.75, 3.050E-02, 60.25, 3.020E-02, 60.75, 2.985E-02, 61.25, 2.950E-02, 61.75, 2.921E-02, 62.25, 2.893E-02, 62.75, 2.865E-02, 63.25, 2.842E-02, 63.75, 2.815E-02, 64.25, 2.788E-02, 64.75, 2.767E-02, 65.25, 2.745E-02, 65.75, 2.724E-02, 66.25, 2.704E-02, 66.75, 2.684E-02, 67.25, 2.664E-02, 67.75, 2.644E-02, 68.25, 2.629E-02, 68.75, 2.615E-02, 69.25, 2.596E-02, 69.75, 2.582E-02, 70.25, 2.568E-02, 70.75, 2.554E-02, 71.25, 2.545E-02, 71.75, 2.531E-02, 72.25, 2.518E-02, 72.75, 2.509E-02, 73.25, 2.501E-02, 73.75, 2.488E-02, 74.25, 2.476E-02, 74.75, 2.467E-02, 75.25, 2.459E-02, 75.75, 2.455E-02, 76.25, 2.452E-02, 76.75, 2.444E-02, 77.25, 2.436E-02, 77.75, 2.428E-02, 78.25, 2.425E-02, 78.75, 2.421E-02, 79.25, 2.418E-02, 79.75, 2.414E-02, 80.25, 2.407E-02, 80.75, 2.400E-02, 81.25, 2.397E-02, 81.75, 2.394E-02, 82.25, 2.390E-02, 82.75, 2.387E-02, 83.25, 2.384E-02, 83.75, 2.381E-02, 84.25, 2.378E-02, 84.75, 2.375E-02, 85.25, 2.372E-02, 85.75, 2.369E-02, 86.25, 2.366E-02, 86.75, 2.364E-02, 87.25, 2.361E-02, 87.75, 2.358E-02, 88.25, 2.355E-02, 88.75, 2.352E-02, 89.25, 2.350E-02, 89.75, 2.351E-02, 90.25, 2.348E-02, 90.75, 2.345E-02, 91.25, 2.343E-02, 91.75, 2.340E-02, 92.25, 2.341E-02, 92.75, 2.338E-02, 93.25, 2.336E-02, 93.75, 2.337E-02, 94.25, 2.334E-02, 94.75, 2.332E-02, 95.25, 2.333E-02, 95.75, 2.334E-02, 96.25, 2.331E-02, 96.75, 2.329E-02, 97.25, 2.330E-02, 97.75, 2.327E-02, 98.25, 2.325E-02, 98.75, 2.326E-02, 99.25, 2.327E-02, 99.75, 2.328E-02 });
    const auto N = raw.size() / 2;
    std::vector<std::pair<double, double>> r(N);
    for (std::size_t i = 0; i < N; ++i) {
        r[i].first = raw[2 * i];
        r[i].second = raw[2 * i + 1];
    }
    return r;
}

std::pair<double, std::map<std::size_t, double>> TG195_soft_tissue()
{
    std::map<std::size_t, double> Zs;
    Zs[1] = 10.5;
    Zs[6] = 25.6;
    Zs[7] = 2.7;
    Zs[8] = 60.2;
    Zs[11] = 0.1;
    Zs[15] = 0.2;
    Zs[16] = 0.3;
    Zs[17] = 0.2;
    Zs[19] = 0.2;

    return std::make_pair(1.03, Zs);
}

std::pair<double, std::map<std::size_t, double>> TG195_pmma()
{
    std::map<std::size_t, double> pmma_w;
    pmma_w[1] = 8.0541;
    pmma_w[6] = 59.9846;
    pmma_w[8] = 31.9613;

    return std::make_pair(1.19, pmma_w);
}

std::pair<double, std::map<std::size_t, double>> TG195_skin()
{
    std::map<std::size_t, double> w;
    w[1] = 9.8;
    w[6] = 17.8;
    w[7] = 5;
    w[8] = 66.7;
    w[15] = .175;
    w[16] = .175;
    w[19] = .175;
    w[20] = .175;

    return std::make_pair(1.09, w);
}

std::pair<double, std::map<std::size_t, double>> TG195_air()
{
    std::map<std::size_t, double> w;
    w[6] = 0.0124;
    w[7] = 75.5268;
    w[8] = 23.1781;
    w[18] = 1.2827;

    return std::make_pair(0.001205, w);
}

std::pair<double, std::map<std::size_t, double>> TG195_water()
{
    std::map<std::size_t, double> w;
    w[1] = 11.1898;
    w[8] = 88.8102;
    return std::make_pair(1.0, w);
}

std::pair<double, std::map<std::size_t, double>> TG195_breast_tissue()
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

    return std::make_pair(d, w);
}

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    dxmc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress, false);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << std::flush << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <BeamType Beam, int LOWENERGYCORRECTION = 2>
    requires(std::same_as<Beam, IsotropicBeamCircle<>> || std::same_as<Beam, IsotropicMonoEnergyBeamCircle<>>)
bool TG195Case1Fluence(std::uint32_t N_threads, bool mammo = false)
{
    constexpr std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 32 : 256;
    constexpr std::uint64_t N_HISTORIES = SAMPLE_RUN ? 1000000 : 1000000;
    constexpr std::uint64_t TOTAL_HIST = N_EXPOSURES * N_HISTORIES;

    using Filter = WorldCylinder<NShells, LOWENERGYCORRECTION>;
    using Scoring = FluenceScore;
    using Material = dxmc::Material<NShells>;

    ResultKeys res;
    res.rCase = "Case 1";
    if (LOWENERGYCORRECTION == 0)
        res.model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        res.model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        res.model = "IA";

    World<Filter, Scoring> world;
    auto [air_density, air_composition] = TG195_air();
    world.setMaterialByWeight(air_composition, air_density);

    world.reserveNumberOfItems(2);
    auto& scoring = world.template addItem<Scoring>({ 0.5, { 0, 0, -100 }, { 0, 0, 1 } });
    scoring.setEnergyStep(0.5);

    // Without filter
    world.build(200);

    Beam beam;
    beam.setNumberOfExposures(N_EXPOSURES);
    beam.setNumberOfParticlesPerExposure(N_HISTORIES);
    beam.setPosition({ 0, 0, 0 });
    beam.setDirection({ 0, 0, -1 });
    beam.setCollimationHalfAngle(std::atan(0.5 / 100.0));

    std::array<double, 2> filter_thickness = { 0, 0 };
    std::array<double, 3> TG195result = { 0, 0, 0 };

    if constexpr (std::same_as<Beam, IsotropicBeamCircle<>>) { // we have a specter
        res.modus = "polyenergetic";
        if (mammo) {
            const auto specter = TG195_30KV();
            beam.setEnergySpecter(specter);
            res.specter = "30 kVp";
            filter_thickness = { 0.03431, 0.07663 };
            TG195result = { 1.857530E-01, 9.281901E-02, 4.642190E-02 };
        } else {
            const auto specter = TG195_100KV();
            beam.setEnergySpecter(specter);
            res.specter = "100 kVp";
            filter_thickness = { 0.3950, 0.9840 };
            TG195result = { 3.452395E-02, 1.725231E-02, 8.621281E-03 };
        }
    } else {
        res.modus = "monoenergetic";
        if (mammo) {
            beam.setEnergy(29.9999);
            res.specter = "30 keV";
            filter_thickness = { 0.2273, 0.4546 };
            TG195result = { 5.732642E-02, 2.868487E-02, 1.435717E-02 };

        } else {
            beam.setEnergy(99.9999);
            res.specter = "100 keV";
            filter_thickness = { 1.511, 3.022 };
            TG195result = { 2.902460E-02, 1.448180E-02, 7.226802E-03 };
        }
    }

    auto AirKerma = [](const auto& spec) {
        const auto u = TG195_mass_en_abs_air();
        double kerma = 0;
        for (std::size_t i = 0; i < u.size(); ++i) {
            kerma += u[i].first * u[i].second * spec[i].second / 100;
        }
        return kerma;
    };

    std::cout << "TG195 Case 1 for " << res.modus << " " << res.specter << " photons with low en model: " << res.model << std::endl;

    ResultPrint print;
    Transport transport;
    transport.setNumberOfThreads(N_threads);

    // None Filter
    auto time_elapsed1 = runDispatcher(transport, world, beam);
    auto fluenceNone = scoring.getFluenceSpecter(TOTAL_HIST);

    res.volume = "NoneFilter";
    res.result = AirKerma(fluenceNone);
    res.result_std = 0;
    res.TG195Result = TG195result[0];
    res.nEvents = scoring.energyScored().numberOfEvents();
    res.nMilliseconds = time_elapsed1.count();
    world.clearDoseScored();
    print(res, true, "Kerma/mm");

    // adding filter
    auto aluminum = Material::byZ(13).value();
    const auto aluminum_dens = AtomHandler::Atom(13).standardDensity;
    auto& filter = world.template addItem<Filter>({ 2, 1, { 0, 0, 0 }, { 0, 0, 1 } });
    filter.setMaterial(aluminum, aluminum_dens);

    // HVL filter
    filter.setHeight(filter_thickness[0]);
    filter.setCenter({ 0, 0, -10 - filter_thickness[0] / 2 });
    world.build(200);
    auto time_elapsed2 = runDispatcher(transport, world, beam);
    auto fluenceHVL = scoring.getFluenceSpecter(TOTAL_HIST);
    res.volume = "HVLFilter";
    res.result = AirKerma(fluenceHVL);
    res.result_std = 0;
    res.TG195Result = TG195result[1];
    res.nEvents = scoring.energyScored().numberOfEvents();
    res.nMilliseconds = time_elapsed2.count();
    world.clearDoseScored();
    print(res, true, "Kerma/mm");

    // QVL filter
    filter.setHeight(filter_thickness[1]);
    filter.setCenter({ 0, 0, -10 - filter_thickness[1] / 2 });
    world.build(200);
    auto time_elapsed3 = runDispatcher(transport, world, beam);
    auto fluenceQVL = scoring.getFluenceSpecter(TOTAL_HIST);
    res.volume = "QVLFilter";
    res.result = AirKerma(fluenceQVL);
    res.result_std = 0;
    res.TG195Result = TG195result[2];
    res.nEvents = scoring.energyScored().numberOfEvents();
    res.nMilliseconds = time_elapsed3.count();
    world.clearDoseScored();
    print(res, true, "Kerma/mm");

    if (LOWENERGYCORRECTION == 1 && std::same_as<Beam, IsotropicMonoEnergyBeamCircle<>> && !mammo) {
        saveImageOfWorld("Case1world.png", world, beam, 60, 120, 300, 5, 100, 0.05);
    }

    return true;
}

template <BeamType Beam, int LOWENERGYCORRECTION = 2>
    requires(std::same_as<Beam, IsotropicBeam<>> || std::same_as<Beam, IsotropicMonoEnergyBeam<>>)
bool TG195Case2AbsorbedEnergy(std::uint32_t N_threads, bool tomo = false)
{
    constexpr std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 32 : 1024;
    constexpr std::uint64_t N_HISTORIES = SAMPLE_RUN ? 1000000 : 1000000;

    using SimpleBox = dxmc::WorldBox<NShells, LOWENERGYCORRECTION>;
    using Box = WorldBoxGrid<NShells, LOWENERGYCORRECTION>;
    using Material = Material<NShells>;

    auto [mat_dens, mat_weights] = TG195_soft_tissue();
    auto mat = Material::byWeight(mat_weights).value();

    World<Box, SimpleBox> world;
    auto [air_density, air_composition] = TG195_air();
    world.setMaterialByWeight(air_composition, air_density);

    const auto box_halfside = 39.0 / 2;
    const double box_height = 20;
    const double box_zbegin = 155;
    const std::array box_aabb = { -box_halfside, -box_halfside, box_zbegin, box_halfside, box_halfside, box_zbegin + box_height };
    world.reserveNumberOfItems(2);
    auto& box = world.template addItem<Box>({ box_aabb });
    box.setVoxelDimensions({ 78, 78, 40 });
    box.setMaterial(mat);
    box.setMaterialDensity(mat_dens);

    auto& scoring_plane = world.template addItem<SimpleBox>({ { -box_halfside, -box_halfside, 180, box_halfside, box_halfside, 180.1 } });
    scoring_plane.setMaterial(dxmc::Material<NShells>::byWeight(air_composition).value(), world.fillMaterialDensity());

    world.build(180);

    Beam beam;

    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        const auto specter = TG195_120KV();
        beam.setEnergySpecter(specter);
    } else {
        beam.setEnergy(56.4);
    }
    beam.setNumberOfExposures(N_EXPOSURES);
    beam.setNumberOfParticlesPerExposure(N_HISTORIES);

    if (tomo) {
        constexpr auto alpha = 15 * DEG_TO_RAD();
        const auto h = 180 * std::tan(alpha);
        beam.setPosition({ 0, -h, 0 });
        const auto y_ang_min = std::atan((h - 39.0 / 2) / 180);
        const auto y_ang_max = std::atan((h + 39.0 / 2) / 180);
        const auto x_ang = std::atan(39.0 / (2 * 180));
        beam.setCollimationHalfAngles(-x_ang, y_ang_min, x_ang, y_ang_max);
    } else {
        const auto collangle = std::atan(39.0 / (2 * 180));
        beam.setCollimationHalfAngles(-collangle, -collangle, collangle, collangle);
    }

    ResultKeys res;
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "120 kVp" : "56.4 keV";
    res.rCase = "Case 2";
    res.modus = tomo ? "tomosynthesis" : "radiography";
    if (LOWENERGYCORRECTION == 0)
        res.model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        res.model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        res.model = "IA";

    std::cout << "TG195 Case 2 for " << res.modus << " orientation and " << res.specter << " photons with low en model: " << res.model << std::endl;

    double TG195_value, TG195_value_min, TG195_value_max;
    std::array<double, 9> TG195_voi_values;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        if (tomo) {
            TG195_value = 30923.125;
            TG195_value_min = 30896.40;
            TG195_value_max = 30967.00;
            TG195_voi_values = { 30.3499448, 23.5176818, 31.6384508, 23.5156038, 8.8986803, 70.5268025, 47.7388610, 20.3142663, 12.509331 };

        } else {
            TG195_value = 33125.975;
            TG195_value_min = 33078.00;
            TG195_value_max = 33189.00;
            TG195_voi_values = { 24.9681848, 24.9503038, 33.5210200, 24.9559023, 24.9661163, 72.7011165, 49.9880020, 21.7295070, 13.4818183 };
        }
    } else {
        if (tomo) {
            TG195_value = 30883.825;
            TG195_value_min = 30870.00;
            TG195_value_max = 30912.80;
            TG195_voi_values = { 33.0807985, 25.4752720, 34.6257073, 25.5054213, 9.7906903, 70.8049988, 51.0616275, 22.2764985, 13.5443103 };

        } else {
            TG195_value = 33171.400;
            TG195_value_min = 33134.90;
            TG195_value_max = 33205.40;
            TG195_voi_values = { 27.010506, 27.003098, 36.667360, 27.008295, 27.007973, 72.858791, 53.345097, 23.832066, 14.595186 };
        }
    }

    if constexpr (LOWENERGYCORRECTION == 1 && std::same_as<Beam, IsotropicMonoEnergyBeam<>>) {
        if (tomo) {
            saveImageOfWorld("Case2worldTomo.png", world, beam, 110, 120, 400, 2);
        } else {
            saveImageOfWorld("Case2world.png", world, beam, 110, 120, 400, 2);
        }
    }

    Transport transport;
    transport.setNumberOfThreads(N_threads);
    auto time_elapsed = runDispatcher(transport, world, beam);

    const auto total_hist = static_cast<double>(N_EXPOSURES * N_HISTORIES);

    double total_ev = 0;
    double total_ev_var = 0;
    std::vector<double> ev_vector(box.totalNumberOfVoxels());
    std::vector<double> ev_var_vector(box.totalNumberOfVoxels());
    std::vector<std::uint64_t> ev_events_vector(box.totalNumberOfVoxels());
    std::uint64_t total_number_events = 0;
    for (std::size_t i = 0; i < box.totalNumberOfVoxels(); ++i) {
        total_ev += box.energyScored(i).energyImparted();
        ev_vector[i] = box.energyScored(i).energyImparted();
        ev_var_vector[i] = box.energyScored(i).variance();
        total_ev_var += box.energyScored(i).variance();
        total_number_events += box.energyScored(i).numberOfEvents();
        ev_events_vector[i] = box.energyScored(i).numberOfEvents();
    }

    res.volume = "Total body";
    res.TG195Result = TG195_value;
    res.result = total_ev / (total_hist / 1000);
    res.result_std = Sigma * std::sqrt(total_ev_var) / (total_hist / 1000) / res.result;
    res.nEvents = total_number_events;
    res.nMilliseconds = time_elapsed.count();
    ResultPrint print;
    print(res, true);

    std::map<std::size_t, std::array<double, 3>> vois;
    vois[3] = { 0, 0, 0 };
    vois[8] = { 0, 0, 3 };
    vois[9] = { 0, 0, 6 };
    vois[7] = { 0, 0, -3 };
    vois[6] = { 0, 0, -6 };
    vois[1] = { 0, -15, 0 };
    vois[5] = { 0, 15, 0 };
    vois[2] = { -15, 0, 0 };
    vois[4] = { 15, 0, 0 };
    for (const auto& [ind, pos] : vois) {
        res.volume = "VOI " + std::to_string(ind);
        res.result = 0;
        res.nEvents = 0;
        const auto ds = 1.5;
        const auto& spacing = box.voxelSpacing();
        const auto& dim = box.voxelDimensions();
        for (std::size_t z = 0; z < dim[2]; ++z) {
            const auto posz = -(dim[2] * spacing[2]) / 2 + z * spacing[2] + spacing[2] / 2;
            if (pos[2] - ds <= posz && posz <= pos[2] + ds) {
                for (std::size_t y = 0; y < dim[1]; ++y) {
                    const auto posy = -(dim[1] * spacing[1]) / 2 + y * spacing[1] + spacing[1] / 2;
                    if (pos[1] - ds <= posy && posy <= pos[1] + ds) {
                        for (std::size_t x = 0; x < dim[0]; ++x) {
                            const auto posx = -(dim[0] * spacing[0]) / 2 + x * spacing[0] + spacing[0] / 2;
                            if (pos[0] - ds <= posx && posx <= pos[0] + ds) {
                                std::array<double, 3> vpos = { posx, posy, posz };
                                const auto boxpos = vectormath::add(vpos, box.center());
                                const auto dIdx = box.gridIndex(boxpos);
                                res.result += ev_vector[dIdx];
                                res.result_std += ev_var_vector[dIdx];
                                res.nEvents += ev_events_vector[dIdx];
                            }
                        }
                    }
                }
            }
        }
        res.TG195Result = TG195_voi_values[ind - 1];
        res.result /= (total_hist / 1000);
        res.result_std = Sigma * std::sqrt(res.result_std) / (total_hist / 1000) / res.result;
        print(res, true);
    }

    return true;
}

template <BeamType Beam, int LOWENERGYCORRECTION = 2>
    requires(std::same_as<Beam, IsotropicBeam<>> || std::same_as<Beam, IsotropicMonoEnergyBeam<>>)
bool TG195Case3AbsorbedEnergy(std::uint32_t N_threads, bool tomo = false)
{
    ResultKeys res;
    res.rCase = "Case 3";
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "30 kVp" : "16.8 keV";
    res.modus = tomo ? "tomosynthesis" : "radiography";
    std::string model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        model = "IA";

    std::cout << "TG195 Case 3 for " << res.modus << " orientation and " << res.specter << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 1024;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 1000000 : 1000000;

    using Box = WorldBox<NShells, LOWENERGYCORRECTION>;
    using Breast = TG195World3Breast<NShells, LOWENERGYCORRECTION>;
    using World = World<Box, Breast>;
    using Material = Material<NShells>;

    const auto [water_d, water_w] = TG195_water();
    const auto [pmma_d, pmma_w] = TG195_pmma();
    const auto [air_d, air_w] = TG195_air();
    const auto [breast_d, breast_w] = TG195_breast_tissue();
    const auto [skin_d, skin_w] = TG195_skin();

    const auto water = Material::byWeight(water_w).value();
    const auto pmma = Material::byWeight(pmma_w).value();
    const auto air = Material::byWeight(air_w).value();
    const auto breasttissue = Material::byWeight(breast_w).value();
    const auto skin = Material::byWeight(skin_w).value();

    World world;
    world.setMaterialByWeight(air_w, air_d);
    world.reserveNumberOfItems(5);

    auto& body = world.template addItem<Box>({ { -17 - dxmc::GEOMETRIC_ERROR(), -15, -15, 0 - dxmc::GEOMETRIC_ERROR(), 15, 15 } });
    body.setMaterial(water, water_d);
    auto& uplate = world.template addItem<Box>({ { 0, -13, 2.5 + dxmc::GEOMETRIC_ERROR(), 14, 13, 2.7 + dxmc::GEOMETRIC_ERROR() } });
    uplate.setMaterial(pmma, pmma_d);
    auto& lplate = world.template addItem<Box>({ { 0, -13, -2.7 - dxmc::GEOMETRIC_ERROR(), 14, 13, -2.5 - dxmc::GEOMETRIC_ERROR() } });
    lplate.setMaterial(pmma, pmma_d);
    auto& breast = world.template addItem<Breast>();
    breast.setSkinMaterial(skin, skin_d);
    breast.setTissueMaterial(breasttissue, breast_d);

    auto& scoringplane = world.template addItem<Box>({ { 0, -13, -2.5 - 1.5 - .1, 14, 13, -2.5 - 1.5 } });
    scoringplane.setMaterial(air, air_d);

    world.build(70);

    Beam beam;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        const auto specter = TG195_30KV();
        beam.setEnergySpecter(specter);
    } else {
        beam.setEnergy(16.8);
    }

    constexpr double source_plane = 66;
    constexpr double source_height = source_plane - 2.5 - 1.5;

    beam.setNumberOfExposures(N_EXPOSURES);
    beam.setNumberOfParticlesPerExposure(N_HISTORIES);
    if (tomo) {
        constexpr auto alpha = 15 * DEG_TO_RAD();
        auto beampos = vectormath::rotate<double>({ 0, 0, source_height }, { 1, 0, 0 }, alpha);
        constexpr double plane_sizey = 13.0;
        const auto height = beampos[2] + 2.5 + 1.5;
        const auto angy_max = std::atan((plane_sizey - beampos[1]) / height);
        const auto angy_min = std::atan((-plane_sizey - beampos[1]) / height);
        constexpr double plane_sizex = 14.0;
        const auto angx_max = std::atan(plane_sizex / height);
        constexpr double angx_min = 0;
        beam.setCollimationHalfAngles(angx_min, -angy_max, angx_max, -angy_min); // since we have ycosine = {0,-1,0}

        beam.setPosition(beampos);
        const std::array<double, 3> cosinex = { 1, 0, 0 };
        const std::array<double, 3> cosiney = { 0, -1, 0 };
        beam.setDirectionCosines(cosinex, cosiney);
    } else {
        beam.setPosition({ 0, 0, source_height });
        beam.setDirectionCosines({ 1, 0, 0 }, { 0, -1, 0 });
        const auto collanglex = std::atan(14.0 / source_plane);
        const auto collangley = std::atan(13.0 / source_plane);
        beam.setCollimationHalfAngles(0, -collangley, collanglex, collangley);
    }

    if constexpr (LOWENERGYCORRECTION == 1 && std::same_as<Beam, IsotropicMonoEnergyBeam<>>) {
        std::string tp;
        if (tomo)
            tp = "Tomo";
        else
            tp = "";
        std::string name = "Case3world" + tp + ".png";
        saveImageOfWorld(name, world, beam, 240, 60, 500, 6, 100, 0.1);
        name = "Case3world2" + tp + ".png";
        saveImageOfWorld(name, world, beam, 180, 60, 500, 6, 100, 0.1);
    }

    Transport transport;
    transport.setNumberOfThreads(N_threads);
    auto time_elapsed = runDispatcher(transport, world, beam);

    ResultPrint print;

    double sim_ev, sim_ev_min, sim_ev_max;
    std::array<double, 7> sim_subvol;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        // specter
        if (tomo) {
            sim_ev = 4188.833;
            sim_ev_min = 4176.5;
            sim_ev_max = 4198.2;
            sim_subvol = { 14.390, 15.825, 15.972, 15.445, 17.171, 5.619, 49.022 };
        } else {
            sim_ev = 4293.433;
            sim_ev_min = 4287.4;
            sim_ev_max = 4296.8;
            sim_subvol = { 16.502, 16.658, 16.814, 16.249, 16.521, 6.041, 50.041 };
        }
    } else {
        if (tomo) {
            sim_ev = 4577.743;
            sim_ev_min = 4565.2;
            sim_ev_max = 4588.2;
            sim_subvol = { 15.217, 16.836, 16.943, 16.431, 18.370, 5.043, 54.974 };
        } else {
            sim_ev = 4697.333;
            sim_ev_min = 4691.0;
            sim_ev_max = 4700.5;
            sim_subvol = { 17.692, 18.070, 17.865, 17.262, 17.768, 5.417, 56.017 };
        }
    }

    res.model = model;

    constexpr auto evNormal = 1000.0 / (N_HISTORIES * N_EXPOSURES);

    res.volume = "Total body";
    res.TG195Result = sim_ev;
    res.result = breast.energyScored(8).energyImparted() * evNormal;
    res.result_std = Sigma * breast.energyScored(8).standardDeviation() * evNormal / res.result;
    res.nEvents = breast.energyScored(8).numberOfEvents();
    res.nMilliseconds = time_elapsed.count();
    print(res, true);

    for (int i = 0; i < 7; ++i) {
        res.volume = "VOI " + std::to_string(i + 1);
        res.TG195Result = sim_subvol[i];
        res.result = breast.energyScored(i).energyImparted() * evNormal;
        res.result_std = Sigma * breast.energyScored(i).standardDeviation() * evNormal / res.result;
        res.nEvents = breast.energyScored(i).numberOfEvents();
        print(res, true);
    }

    return true;
}

template <int LOWENERGYCORRECTION = 2>
bool TG195Case41AbsorbedEnergy(std::uint32_t N_threads, bool specter = false, bool large_collimation = false)
{
    std::string model;
    if (LOWENERGYCORRECTION == 0)
        model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        model = "IA";

    std::cout << "TG195 Case 4.1 for ";
    if (large_collimation)
        std::cout << "80mm collimation and ";
    else
        std::cout << "10mm collimation and ";
    if (specter)
        std::cout << "120 kVp";
    else
        std::cout << "56.4 keV";
    std::cout << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 1024;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 100000 : 1000000;

    using Cylindar = DepthDose<NShells, LOWENERGYCORRECTION>;
    using World = World<Cylindar>;
    using Material = Material<NShells>;

    auto [mat_dens, mat_weights] = TG195_pmma();

    auto mat = Material::byWeight(mat_weights).value();

    World world;
    auto [air_density, air_composition] = TG195_air();
    world.setMaterialByWeight(air_composition, air_density);

    auto& cylinder = world.template addItem<Cylindar>({ 16, 300, 600 });
    world.build(60);
    cylinder.setMaterial(mat);
    cylinder.setMaterialDensity(mat_dens);
    std::chrono::milliseconds time_elapsed;
    if (specter) {
        using Beam = IsotropicBeam<>;
        Beam beam({ -60, 0, 0 }, { { { 0, 1, 0 }, { 0, 0, 1 } } });
        const auto specter = TG195_120KV();
        beam.setEnergySpecter(specter);
        const auto collangle_y = std::atan(16.0 / 60);
        const auto collangle_z = large_collimation ? std::atan(4.0 / 60) : std::atan(0.5 / 60);
        beam.setCollimationHalfAngles({ -collangle_y, -collangle_z, collangle_y, collangle_z });
        beam.setNumberOfExposures(N_EXPOSURES);
        beam.setNumberOfParticlesPerExposure(N_HISTORIES);

        Transport transport;
        transport.setNumberOfThreads(N_threads);
        time_elapsed = runDispatcher(transport, world, beam);
    } else {
        using Beam = IsotropicMonoEnergyBeam<>;
        Beam beam({ -60, 0, 0 }, { { { 0, 1, 0 }, { 0, 0, 1 } } }, 56.4);
        const auto collangle_y = std::atan(16.0 / 60);
        const auto collangle_z = large_collimation ? std::atan(4.0 / 60) : std::atan(0.5 / 60);
        beam.setCollimationHalfAngles({ -collangle_y, -collangle_z, collangle_y, collangle_z });
        beam.setNumberOfExposures(N_EXPOSURES);
        beam.setNumberOfParticlesPerExposure(N_HISTORIES);

        if constexpr (LOWENERGYCORRECTION == 1) {
            if (!specter) {
                std::string tp;
                if (large_collimation)
                    tp = "80";
                else
                    tp = "10";
                const auto lenght = std::sqrt(120.0 * 120.0 + 120 * std::tan(collangle_y) * 120 * std::tan(collangle_y));
                std::string name = "Case41world_" + tp + ".png";
                saveImageOfWorld(name, world, beam, 270, 90, 500, 1.6, lenght);
                name = "Case41world2_" + tp + ".png";
                saveImageOfWorld(name, world, beam, 270, 30, 500, 1.6, lenght);
            }
        }

        Transport transport;
        transport.setNumberOfThreads(N_threads);
        time_elapsed = runDispatcher(transport, world, beam);
    }
    const std::array<double, 4> voi_locations = { 0, 1, 2, 3 };
    std::array<double, 4> ev_history, ev_history_var;
    ev_history.fill(0);
    ev_history_var.fill(0);
    std::array<std::uint64_t, 4> ev_events = { 0, 0, 0, 0 };

    for (const auto& [z, energyScored] : cylinder.depthEnergyScored()) {
        for (std::size_t i = 0; i < voi_locations.size(); ++i) {
            const auto zmin = voi_locations[i] - 0.5;
            const auto zmax = voi_locations[i] + 0.5;
            if (zmin < z && z < zmax) {
                ev_history[i] += energyScored.energyImparted();
                ev_history_var[i] += energyScored.variance();
                ev_events[i] += energyScored.numberOfEvents();
            }
        }
    }

    for (std::size_t i = 0; i < voi_locations.size(); ++i) {
        ev_history[i] /= ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        ev_history_var[i] = Sigma * std::sqrt(ev_history_var[i]) / ((N_HISTORIES * N_EXPOSURES) / 1000.0) / ev_history[i];
    }

    ResultPrint print;
    ResultKeys res;
    res.specter = specter ? "120 kVp" : "56.4 keV";
    res.rCase = "Case 4.1";
    res.model = model;
    res.modus = large_collimation ? "80mm collimation" : "10mm collimation";
    res.nMilliseconds = time_elapsed.count();

    std::array<double, 4> tg195;
    if (specter && large_collimation) {
        tg195 = { 3586.59, 3537.84, 3378.99, 2672.21 };
    }
    if (specter && !large_collimation) {
        tg195 = { 13137.02, 2585.47, 1706.86, 1250.61 };
    }
    if (!specter && large_collimation) {
        tg195 = { 3380.39, 3332.64, 3176.44, 2559.58 };
    }
    if (!specter && !large_collimation) {
        tg195 = { 11592.27, 2576.72, 1766.85, 1330.53 };
    }

    for (std::size_t i = 0; i < voi_locations.size(); ++i) {
        res.volume = "VOI " + std::to_string(i + 1);
        res.TG195Result = tg195[i];
        res.result = ev_history[i];
        res.result_std = ev_history_var[i];
        res.nEvents = ev_events[i];
        print(res, true);
    }
    return true;
}

template <BeamType Beam, int LOWENERGYCORRECTION = 2>
    requires(std::same_as<Beam, IsotropicBeam<>> || std::same_as<Beam, IsotropicMonoEnergyBeam<>>)
bool TG195Case42AbsorbedEnergy(std::uint32_t N_threads, bool large_collimation = false)
{
    std::string model;
    if (LOWENERGYCORRECTION == 0)
        model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        model = "IA";

    std::cout << "TG195 Case 4.2 for ";
    if (large_collimation)
        std::cout << "80mm collimation and ";
    else
        std::cout << "10mm collimation and ";
    if constexpr (std::same_as<Beam, IsotropicBeam<>>)
        std::cout << "120 kVp";
    else
        std::cout << "56.4 keV";
    std::cout << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 1024;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 10000 : 2000000;

    using Cylindar = TG195World42<NShells, LOWENERGYCORRECTION>;
    using World = World<Cylindar>;
    using Material = Material<NShells>;
    auto [mat_dens, mat_weights] = TG195_pmma();
    auto mat = Material::byWeight(mat_weights).value();

    World world;
    auto [air_density, air_composition] = TG195_air();
    world.setMaterialByWeight(air_composition, air_density);

    auto& cylinder = world.template addItem<Cylindar>({ 16, 600 });
    world.build(90);
    cylinder.setMaterial(mat);
    cylinder.setMaterialDensity(mat_dens);

    ResultPrint print;
    ResultKeys res;
    if (LOWENERGYCORRECTION == 0)
        res.model = "NoneLC";
    else if (LOWENERGYCORRECTION == 1)
        res.model = "Livermore";
    else
        res.model = "IA";
    res.modus = large_collimation ? "80mm collimation" : "10mm collimation";
    res.rCase = "Case 4.2";
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "120 kVp" : "56.4 keV";

    Beam beam({ -60, 0, 0 }, { { { 0, 1, 0 }, { 0, 0, 1 } } });
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        const auto specter = TG195_120KV();
        beam.setEnergySpecter(specter);
    } else {
        beam.setEnergy(56.4);
    }
    const auto collangle_y = std::atan(16.0 / 60);
    const auto collangle_z = large_collimation ? std::atan(4.0 / 60) : std::atan(0.5 / 60);
    beam.setCollimationHalfAngles({ -collangle_y, -collangle_z, collangle_y, collangle_z });
    beam.setNumberOfExposures(N_EXPOSURES);
    beam.setNumberOfParticlesPerExposure(N_HISTORIES);

    // setting up benchmarking values
    std::array<double, 36> sim_ev_center, sim_ev_pher;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        if (large_collimation) {
            sim_ev_center = { 10.878025, 10.9243, 10.884625, 10.89795, 10.87265, 10.902675, 10.8994, 10.880875, 10.875475, 10.8862, 10.895975, 10.88105, 10.899600, 10.886225, 10.893400, 10.894200, 10.879025, 10.885500, 10.894125, 10.8898, 10.8916, 10.895875, 10.889525, 10.889775, 10.89365, 10.901875, 10.894475, 10.906975, 10.888025, 10.877475, 10.883325, 10.875925, 10.8881, 10.886775, 10.88975, 10.900075 };
            sim_ev_pher = { 115.34325, 113.76275, 109.16925, 101.706, 91.562975, 78.39105, 61.388325, 40.08625, 22.471075, 11.781725, 6.14551, 3.42218, 2.05605, 1.35319, 0.96088275, 0.743808, 0.61922025, 0.55457575, 0.5309405, 0.55428325, 0.6219885, 0.74445025, 0.96480125, 1.3481875, 2.0611025, 3.4154, 6.15532, 11.7854, 22.461525, 40.13715, 61.42595, 78.328975, 91.481375, 101.61325, 109.10425, 113.8365 };
        } else {
            sim_ev_center = { 11.35475, 11.399475, 11.38755, 11.402175, 11.400225, 11.370625, 11.402625, 11.37715, 11.385375, 11.4096, 11.399825, 11.376675, 11.369800, 11.361300, 11.388650, 11.37793, 11.3692, 11.371525, 11.3807, 11.36, 11.37645, 11.379075, 11.379975, 11.368975, 11.377675, 11.38375, 11.3871, 11.3844, 11.37395, 11.3821, 11.371825, 11.395575, 11.379075, 11.372125, 11.4001, 11.39895 };
            sim_ev_pher = { 116.7915, 115.31075, 110.532, 103.1535, 92.9892, 79.7158, 62.57135, 40.92065, 22.966025, 12.169, 6.4216075, 3.5657225, 2.1609475, 1.418235, 1.014725, 0.780394, 0.64618, 0.5763575, 0.559597, 0.5776625, 0.648646, 0.7762825, 1.008508, 1.412125, 2.1584, 3.577945, 6.4296425, 12.1714, 23.039025, 40.926825, 62.451325, 79.759575, 93.074725, 103.269, 110.595, 115.27 };
        }
    } else {
        if (large_collimation) {
            sim_ev_center = { 11.630625, 11.632925, 11.617175, 11.624825, 11.633075, 11.600225, 11.61655, 11.6235, 11.592875, 11.6258, 11.612, 11.612775, 11.608675, 11.623975, 11.611325, 11.6174, 11.6234, 11.627975, 11.60745, 11.632875, 11.628275, 11.6239, 11.61645, 11.617375, 11.621775, 11.6178, 11.6444, 11.61515, 11.626375, 11.64605, 11.63335, 11.628425, 11.622, 11.6198, 11.59835, 11.609925 };
            sim_ev_pher = { 99.665175, 98.300175, 94.509325, 88.48595, 80.1113, 69.261125, 55.124425, 37.2351, 21.754, 11.767, 6.269845, 3.460205, 2.073845, 1.3435025, 0.94542875, 0.71714775, 0.58643225, 0.52293525, 0.49996925, 0.5225535, 0.5875545, 0.719903, 0.94029425, 1.3461175, 2.07283, 3.4740625, 6.2371725, 11.79715, 21.73405, 37.28175, 55.1853, 69.2588, 80.036275, 88.397, 94.640025, 98.332825 };
        } else {
            sim_ev_center = { 12.168675, 12.112550, 12.163000, 12.135800, 12.088100, 12.102750, 12.135925, 12.124250, 12.149150, 12.149575, 12.148200, 12.159550, 12.130775, 12.148475, 12.153925, 12.138300, 12.134750, 12.148675, 12.157600, 12.145200, 12.156650, 12.162250, 12.156625, 12.156875, 12.139050, 12.156100, 12.159325, 12.137025, 12.149725, 12.095750, 12.147100, 12.123100, 12.121100, 12.136800, 12.134825, 12.139550 };
            sim_ev_pher = { 101.293750, 99.802475, 96.108725, 89.966275, 81.369475, 70.542400, 56.283500, 38.128300, 22.348775, 12.166625, 6.515755, 3.643485, 2.180365, 1.413790, 0.984371, 0.751024, 0.616579, 0.545221, 0.521630, 0.545158, 0.616436, 0.750186, 0.991487, 1.413838, 2.174915, 3.629998, 6.514243, 12.200725, 22.286500, 38.146275, 56.287325, 70.641650, 81.424625, 89.820150, 96.169950, 99.924350 };
        }
    }

    Transport transport;
    transport.setNumberOfThreads(N_threads);

    const std::array<double, 3> co_x = { 0, 1, 0 };
    const std::array<double, 3> co_y = { 0, 0, 1 };
    const std::array<double, 3> pos = { -60, 0, 0 };
    for (std::size_t i = 0; i < 36; ++i) {
        const auto angInt = i * 10;
        const auto angle = angInt * DEG_TO_RAD();
        auto x = vectormath::rotate(co_x, { 0, 0, 1 }, angle);
        auto p_ang = vectormath::rotate(pos, { 0, 0, 1 }, angle);
        beam.setPosition(p_ang);
        beam.setDirectionCosines(x, co_y);

        auto time_elepased = runDispatcher(transport, world, beam);
        res.nMilliseconds = time_elepased.count();

        std::cout << "Angle " << angInt;
        res.modus = large_collimation ? "Pherifery 80mm collimation" : "Pherifery 10mm collimation";
        res.volume = std::to_string(angInt);
        res.TG195Result = sim_ev_pher[i];
        res.result = cylinder.energyScoredPeriferyCylinder().energyImparted() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        res.result_std = Sigma * cylinder.energyScoredPeriferyCylinder().standardDeviation() / ((N_HISTORIES * N_EXPOSURES) / 1000.0) / res.result;
        res.nEvents = cylinder.energyScoredPeriferyCylinder().numberOfEvents();
        std::cout << " Pherifery: " << res.result << " sim/TG195: [" << (res.result / sim_ev_pher[i] - 1) * 100 << "%]";
        print(res, false);
        res.modus = large_collimation ? "Center 80mm collimation" : "Center 10mm collimation";
        res.volume = std::to_string(angInt);
        res.TG195Result = sim_ev_center[i];
        res.result = cylinder.energyScoredCenterCylinder().energyImparted() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        res.result_std = Sigma * cylinder.energyScoredCenterCylinder().standardDeviation() / ((N_HISTORIES * N_EXPOSURES) / 1000.0) / res.result;
        res.nEvents = cylinder.energyScoredCenterCylinder().numberOfEvents();
        std::cout << " Center: " << res.result << " sim/TG195: [" << (res.result / sim_ev_center[i] - 1) * 100 << "%]" << std::endl;
        print(res, false);
    }

    return true;
}

template <typename T>
std::vector<T> readBinaryArray(const std::string& path, std::size_t array_size)
{
    std::vector<T> buffer;

    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if (!ifs) {
        return buffer;
    }

    auto end = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    auto buffer_size = std::size_t(end - ifs.tellg());
    auto dim_size = array_size * sizeof(T);
    if (dim_size != buffer_size) {
        return buffer;
    }

    if (buffer_size == 0) { // avoid undefined behavior
        return buffer;
    }

    buffer.resize(array_size, 0);

    if (!ifs.read(reinterpret_cast<char*>(buffer.data()), buffer_size)) {
        return buffer;
    }
    return buffer;
}

template <std::size_t NMATSHELLS = 5, int LOWENERGYCORRECTION = 2, int TRANSPARENTVOXEL = 255>
std::pair<AAVoxelGrid<NMATSHELLS, LOWENERGYCORRECTION, TRANSPARENTVOXEL>, std::vector<std::pair<double, std::string>>> generateTG195World5()
{

    const std::vector<std::map<std::size_t, double>> matWeights = {
        {
            // Air
            { 6, 0.0124 },
            { 7, 75.5268 },
            { 8, 23.1781 },
            { 18, 1.2827 },
        },
        {
            // Cushion
            { 1, 7.8 },
            { 6, 64.7 },
            { 7, 8.4 },
            { 8, 19.1 },
        },
        {
            // Carbon
            { 6, 100 },
        },
        {
            // Soft tissue
            { 1, 10.5 },
            { 6, 25.6 },
            { 7, 2.7 },
            { 8, 60.2 },
            { 11, 0.1 },
            { 15, 0.2 },
            { 16, 0.3 },
            { 17, 0.2 },
            { 19, 0.2 },
        },
        {
            // Heart
            { 1, 10.4 },
            { 6, 13.9 },
            { 7, 2.9 },
            { 8, 71.8 },
            { 11, 0.1 },
            { 15, 0.2 },
            { 16, 0.2 },
            { 17, 0.2 },
            { 19, 0.3 },
        },
        {
            // Lung
            { 1, 10.3 },
            { 6, 10.5 },
            { 7, 3.1 },
            { 8, 74.9 },
            { 11, 0.2 },
            { 15, 0.2 },
            { 16, 0.3 },
            { 17, 0.3 },
            { 19, 0.2 },
        },
        {
            // Liver
            { 1, 10.2 },
            { 6, 13.9 },
            { 7, 3 },
            { 8, 71.6 },
            { 11, 0.2 },
            { 15, 0.3 },
            { 16, 0.3 },
            { 17, 0.2 },
            { 19, 0.3 },
        },
        {
            // Gallbladder
            { 1, 10.5 },
            { 6, 25.6 },
            { 7, 2.7 },
            { 8, 60.2 },
            { 11, 0.1 },
            { 15, 0.2 },
            { 16, 0.3 },
            { 17, 0.2 },
            { 19, 0.2 },
        },
        {
            // Spleen
            { 1, 10.3 },
            { 6, 11.3 },
            { 7, 3.2 },
            { 8, 74.1 },
            { 11, 0.1 },
            { 15, 0.3 },
            { 16, 0.2 },
            { 17, 0.2 },
            { 19, 0.3 },
        },
        {
            // Stomach
            { 1, 10.6 },
            { 6, 11.5 },
            { 7, 2.2 },
            { 8, 75.1 },
            { 11, 0.1 },
            { 15, 0.1 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.1 },
        },
        {
            // Large intestine
            { 1, 10.6 },
            { 6, 11.5 },
            { 7, 2.2 },
            { 8, 75.1 },
            { 11, 0.1 },
            { 15, 0.1 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.1 },
        },
        {
            // Pancreas
            { 1, 10.6 },
            { 6, 16.9 },
            { 7, 2.2 },
            { 8, 69.4 },
            { 11, 0.2 },
            { 15, 0.2 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.2 },
        },
        {
            // Adrenal
            { 1, 10.5 },
            { 6, 25.6 },
            { 7, 2.7 },
            { 8, 60.2 },
            { 11, 0.1 },
            { 15, 0.2 },
            { 16, 0.3 },
            { 17, 0.2 },
            { 19, 0.2 },
        },
        {
            // Thyroid
            { 1, 10.4 },
            { 6, 11.9 },
            { 7, 2.4 },
            { 8, 74.5 },
            { 11, 0.2 },
            { 15, 0.1 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.1 },
            { 53, 0.1 },
        },
        {
            // Thymus
            { 1, 10.5 },
            { 6, 25.6 },
            { 7, 2.7 },
            { 8, 60.2 },
            { 11, 0.1 },
            { 15, 0.2 },
            { 16, 0.3 },
            { 17, 0.2 },
            { 19, 0.2 },
        },
        {
            // Small Intestine
            { 1, 10.6 },
            { 6, 11.5 },
            { 7, 2.2 },
            { 8, 75.1 },
            { 11, 0.1 },
            { 15, 0.1 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.1 },
        },
        {
            // Esophagus
            { 1, 10.6 },
            { 6, 11.5 },
            { 7, 2.2 },
            { 8, 75.1 },
            { 11, 0.1 },
            { 15, 0.1 },
            { 16, 0.1 },
            { 17, 0.2 },
            { 19, 0.1 },
        },
        {
            // Skin
            { 1, 10 },
            { 6, 20.4 },
            { 7, 4.2 },
            { 8, 64.5 },
            { 11, 0.2 },
            { 15, 0.1 },
            { 16, 0.2 },
            { 17, 0.3 },
            { 19, 0.1 },
        },
        {
            // Breast
            { 1, 11.2 },
            { 6, 61.9 },
            { 7, 1.7 },
            { 8, 25.1 },
            { 15, 0.025 },
            { 16, 0.025 },
            { 19, 0.025 },
            { 20, 0.025 },
        },
        {
            // Cortical bone
            { 1, 3.4 },
            { 6, 15.5 },
            { 7, 4.2 },
            { 8, 43.5 },
            { 11, 0.1 },
            { 12, 0.2 },
            { 15, 10.3 },
            { 16, 0.3 },
            { 20, 22.5 },
        },
    };

    std::vector<Material<NMATSHELLS>> materials;

    std::transform(matWeights.cbegin(), matWeights.cend(), std::back_inserter(materials), [=](const auto& f) {
        auto mat_cand = Material<NMATSHELLS>::byWeight(f);
        return mat_cand.value();
    });

    constexpr std::array<std::size_t, 3> dim = { 500, 320, 260 };
    constexpr auto size = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    constexpr std::array<double, 3> spacing = { 0.1f, 0.1f, 0.1f };

    AAVoxelGrid<NMATSHELLS, LOWENERGYCORRECTION, TRANSPARENTVOXEL> grid;

    auto matArray = readBinaryArray<std::uint8_t>("case5world.bin", size);

    std::vector<std::pair<double, std::string>> matInfo;
    matInfo.push_back(std::make_pair(.001205, "Air"));
    matInfo.push_back(std::make_pair(.075, "Cushion Foam"));
    matInfo.push_back(std::make_pair(1.20, "Carbon fiber"));
    matInfo.push_back(std::make_pair(1.03, "Soft tissue"));
    matInfo.push_back(std::make_pair(1.05, "Heart"));
    matInfo.push_back(std::make_pair(0.26, "Lung"));
    matInfo.push_back(std::make_pair(1.06, "Liver"));
    matInfo.push_back(std::make_pair(1.03, "Gallbladder"));
    matInfo.push_back(std::make_pair(1.06, "Spleen"));
    matInfo.push_back(std::make_pair(1.03, "Stomach"));
    matInfo.push_back(std::make_pair(1.03, "Large Intestine"));
    matInfo.push_back(std::make_pair(1.04, "Pancreas"));
    matInfo.push_back(std::make_pair(1.03, "Adrenal"));
    matInfo.push_back(std::make_pair(1.05, "Thyroid"));
    matInfo.push_back(std::make_pair(1.03, "Thymus"));
    matInfo.push_back(std::make_pair(1.03, "Small Intestine"));
    matInfo.push_back(std::make_pair(1.03, "Esophagus"));
    matInfo.push_back(std::make_pair(1.09, "Skin"));
    matInfo.push_back(std::make_pair(0.93, "Breast"));
    matInfo.push_back(std::make_pair(1.92, "Cortical Bone"));

    auto res = std::make_pair(grid, matInfo);

    std::vector<double> density(size, 0.0);
    std::transform(std::execution::par_unseq, matArray.cbegin(), matArray.cend(), density.begin(), [&](const auto ind) {
        return matInfo[ind].first;
    });

    res.first.setData(dim, density, matArray, materials);
    res.first.setSpacing(spacing);
    return res;
}

template <BeamType B, int LOWENERGYCORRECTION = 2>
    requires(std::same_as<B, IsotropicBeam<>> || std::same_as<B, IsotropicMonoEnergyBeam<>>)
bool TG195Case5AbsorbedEnergy(std::uint32_t N_threads)
{
    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 1024;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 1000000 : 1000000;

    constexpr int TRANSPARENTVOXELS = 255;
    using World = World<AAVoxelGrid<NShells, LOWENERGYCORRECTION, TRANSPARENTVOXELS>>;
    auto [grid_object, matInf] = generateTG195World5<NShells, LOWENERGYCORRECTION, TRANSPARENTVOXELS>();

    World world;
    auto& grid = world.addItem(grid_object, "Tissue");
    world.build(60);

    ResultPrint print;
    ResultKeys res;
    if (LOWENERGYCORRECTION == 0)
        res.model = "NoneLC";
    else if (LOWENERGYCORRECTION == 1)
        res.model = "Livermore";
    else
        res.model = "IA";
    res.rCase = "Case 5";

    B beam;
    if constexpr (std::same_as<B, IsotropicBeam<>>) {
        const auto specter = TG195_120KV();
        beam.setEnergySpecter(specter);
        res.specter = "120 kVp";
    } else {
        beam.setEnergy(56.4);
        res.specter = "56.4 keV";
    }

    std::array<std::array<double, 17>, 8> tg195_doses;
    if constexpr (std::same_as<B, IsotropicBeam<>>) {
        tg195_doses[0] = { 12374.975000, 2917.750000, 1275.860000, 612.314000, 5.779858, 16.682150, 121.044500, 15.161925, 8.171573, 0.150686, 1.653300, 40.660375, 9.781930, 33.369275, 559.765000, 21.488925, 7727.772500 };
        tg195_doses[1] = { 12594.500000, 1801.817500, 1007.282500, 612.421000, 5.742833, 9.973140, 76.976000, 8.730155, 5.864173, 0.132321, 1.388425, 33.380250, 7.117933, 30.519400, 538.165250, 17.763675, 6631.552500 };
        tg195_doses[2] = { 10648.025000, 737.539250, 640.750500, 447.051750, 3.697750, 6.575955, 36.069700, 3.290558, 3.191215, 0.104060, 1.033353, 15.512825, 3.375828, 21.384350, 348.679250, 6.815653, 4814.567500 };
        tg195_doses[3] = { 10137.700000, 730.785250, 716.213750, 389.105000, 2.617823, 17.390300, 56.726925, 5.195065, 4.383270, 0.144614, 1.239305, 11.312000, 3.718875, 24.866200, 223.370250, 2.190193, 7437.655000 };
        tg195_doses[4] = { 10250.375000, 1211.352500, 1043.237500, 385.775000, 2.573733, 31.087675, 102.259750, 9.739023, 6.914105, 0.178788, 1.508200, 10.662950, 5.453338, 30.245775, 143.224000, 2.957988, 9718.545000 };
        tg195_doses[5] = { 10069.325000, 1121.222500, 687.389500, 243.211000, 1.624610, 33.913700, 107.841250, 11.186650, 6.924788, 0.153553, 1.333120, 8.570210, 5.341703, 29.445475, 232.172500, 2.136605, 7265.825000 };
        tg195_doses[6] = { 10666.100000, 1503.582500, 601.060250, 164.330750, 1.323968, 30.147475, 122.857750, 16.103400, 7.157005, 0.120301, 1.202235, 10.514325, 6.824010, 22.474550, 365.628500, 7.175110, 5072.787500 };
        tg195_doses[7] = { 12488.675000, 2558.895000, 876.999750, 375.597250, 3.763275, 24.311350, 138.639750, 18.520350, 8.427903, 0.134507, 1.485195, 26.485200, 9.596425, 27.676950, 562.674250, 18.287350, 6436.260000 };
    } else {
        tg195_doses[0] = { 11574.275000, 3086.415000, 1301.172500, 679.467000, 6.366870, 17.572300, 134.399750, 16.729425, 8.792378, 0.152669, 1.742475, 44.221675, 10.718650, 36.896025, 456.355000, 21.679725, 8761.232500 };
        tg195_doses[1] = { 11761.250000, 1932.042500, 1045.652500, 683.628250, 6.328008, 9.911785, 82.980825, 9.224485, 6.182205, 0.130073, 1.472253, 36.387000, 7.644288, 33.457325, 437.274750, 17.972150, 7669.282500 };
        tg195_doses[2] = { 9975.200000, 786.733750, 679.416750, 495.341250, 3.943625, 6.274665, 36.212625, 3.145198, 3.146200, 0.101546, 1.056510, 16.899150, 3.379468, 23.030800, 285.866500, 6.795723, 5611.520000 };
        tg195_doses[3] = { 9581.672500, 765.345000, 755.185250, 421.935500, 2.624548, 18.052325, 58.050950, 5.069595, 4.381183, 0.146320, 1.267433, 12.025825, 3.612820, 26.356000, 189.367000, 2.067310, 8510.072500 };
        tg195_doses[4] = { 9704.905000, 1265.142500, 1085.652500, 411.344500, 2.523603, 33.788550, 109.307750, 10.178175, 7.151828, 0.183346, 1.565265, 11.223200, 5.506005, 32.106700, 129.805500, 2.818055, 11003.750000 };
        tg195_doses[5] = { 9487.977500, 1202.942500, 721.263500, 251.049250, 1.516695, 37.435375, 117.358250, 12.021875, 7.309943, 0.158866, 1.379200, 8.953935, 5.484705, 31.534025, 195.788500, 1.996303, 8367.435000 };
        tg195_doses[6] = { 9970.395000, 1625.792500, 636.469500, 168.172000, 1.245140, 33.399550, 136.692250, 17.883625, 7.677150, 0.123194, 1.262395, 11.258675, 7.289940, 24.393100, 298.963500, 7.180418, 5876.965000 };
        tg195_doses[7] = { 11649.900000, 2725.147500, 916.402500, 409.546000, 3.993910, 26.541075, 155.303000, 20.681400, 9.132600, 0.134053, 1.550688, 28.927825, 10.487675, 30.430800, 455.474000, 18.393875, 7391.310000 };
    }

    const auto collangle_y = std::atan(25.0 / 60);
    const auto collangle_z = std::atan(0.5 / 60);
    beam.setCollimationHalfAngles({ -collangle_y, -collangle_z, collangle_y, collangle_z });
    beam.setNumberOfExposures(N_EXPOSURES);
    beam.setNumberOfParticlesPerExposure(N_HISTORIES);

    Transport transport;
    transport.setNumberOfThreads(N_threads);
    const std::array<double, 3> co_x = { -1, 0, 0 };
    const std::array<double, 3> co_y = { 0, 0, 1 };
    const std::array<double, 3> pos = { 0, -60, 0 };
    for (std::size_t angInt = 0; angInt < 8; angInt = ++angInt) {
        const auto angle = -(angInt * 45 * DEG_TO_RAD());
        auto x = vectormath::rotate(co_x, { 0, 0, 1 }, angle);
        auto p_ang = vectormath::rotate(pos, { 0, 0, 1 }, angle);
        beam.setPosition(p_ang);
        beam.setDirectionCosines(x, co_y);
        res.modus = std::to_string(angInt * 45);

        if constexpr (LOWENERGYCORRECTION == 1 && std::same_as<B, IsotropicMonoEnergyBeam<>>) {
            std::string name = "Case5world_" + std::to_string(angInt) + ".png";
            saveImageOfWorld(name, world, beam, 60, 60, 500, 2, 100, 0.1);
        }

        std::cout << "Case 5 specter: " << res.specter << " angle: " << res.modus << " model: " << res.model << std::endl;
        auto time_elapsed = runDispatcher(transport, world, beam);

        res.nMilliseconds = time_elapsed.count();

        const auto doseScore = grid.getEnergyScores();
        const auto materialIndex = grid.getMaterialIndex();

        std::uint8_t matIdx = 0;
        for (const auto& [density, material_name] : matInf) {
            if (matIdx > 2) {
                res.volume = material_name;
                const auto tg195_idx = matIdx - 3;
                res.TG195Result = tg195_doses[angInt][tg195_idx];
                const auto ei_tot = std::transform_reduce(std::execution::par_unseq, doseScore.cbegin(), doseScore.cend(), materialIndex.cbegin(), 0.0, std::plus<>(), [matIdx](const auto& energyScored, auto ind) -> double {
                    if (ind == matIdx)
                        return energyScored.energyImparted();
                    else
                        return 0;
                });
                res.result = ei_tot / ((N_HISTORIES * N_EXPOSURES) / 1000);

                const auto ei_var = std::transform_reduce(std::execution::par_unseq, doseScore.cbegin(), doseScore.cend(), materialIndex.cbegin(), 0.0, std::plus<>(), [matIdx](const auto& energyScored, auto ind) -> double {
                    if (ind == matIdx)
                        return energyScored.variance();
                    else
                        return 0;
                });
                res.result_std = Sigma * std::sqrt(ei_var) / ((N_HISTORIES * N_EXPOSURES) / 1000) / res.result;
                const auto events = std::transform_reduce(std::execution::par_unseq, doseScore.cbegin(), doseScore.cend(), materialIndex.cbegin(), 0.0, std::plus<>(), [matIdx](const auto& energyScored, auto ind) -> double {
                    if (ind == matIdx)
                        return energyScored.numberOfEvents();
                    else
                        return 0;
                });
                res.nEvents = events;
                print(res, true);
            }
            matIdx++;
        }
        world.clearEnergyScored();
        world.clearDoseScored();
    }

    return true;
}

template <int LOWENERGYCORRECTION>
bool runAll(std::uint32_t N_threads)
{
    auto success = true;

    success = success && TG195Case1Fluence<IsotropicMonoEnergyBeamCircle<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case1Fluence<IsotropicMonoEnergyBeamCircle<>, LOWENERGYCORRECTION>(N_threads, true);
    success = success && TG195Case1Fluence<IsotropicBeamCircle<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case1Fluence<IsotropicBeamCircle<>, LOWENERGYCORRECTION>(N_threads, true);

    success = success && TG195Case2AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case2AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, true);
    success = success && TG195Case2AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case2AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, true);

    success = success && TG195Case3AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case3AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, true);
    success = success && TG195Case3AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case3AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, true);

    success = success && TG195Case41AbsorbedEnergy<LOWENERGYCORRECTION>(N_threads, false, false);
    success = success && TG195Case41AbsorbedEnergy<LOWENERGYCORRECTION>(N_threads, false, true);
    success = success && TG195Case41AbsorbedEnergy<LOWENERGYCORRECTION>(N_threads, true, false);
    success = success && TG195Case41AbsorbedEnergy<LOWENERGYCORRECTION>(N_threads, true, true);

    success = success && TG195Case42AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case42AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads, true);
    success = success && TG195Case42AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, false);
    success = success && TG195Case42AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads, true);

    success = success && TG195Case5AbsorbedEnergy<IsotropicMonoEnergyBeam<>, LOWENERGYCORRECTION>(N_threads);
    success = success && TG195Case5AbsorbedEnergy<IsotropicBeam<>, LOWENERGYCORRECTION>(N_threads);

    return success;
}

struct Arguments {
    std::uint32_t N_threads = 0;
    enum class LECorrection {
        All,
        None,
        Livermore,
        IA
    };
    LECorrection correction = LECorrection::All;
    bool exit = false;
    std::string filename = "validationTable.txt";

    std::string correctionString() const
    {
        switch (correction) {
        case LECorrection::None:
            return "None";
        case LECorrection::Livermore:
            return "Livermore";
        case LECorrection::IA:
            return "IA";
        }
        return "All";
    }
};

Arguments argparse(int argc, char* argv[])
{
    Arguments args;
    const std::vector<std::string_view> args_str(argv + 1, argv + argc);

    std::size_t i = 0;
    while (i < args_str.size()) {
        auto arg = args_str[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << "Help for validation run of dxmclib\n";
            std::cout << "-j or --jobs [N] where N is number of threads, default all available threads\n";
            std::cout << "-b or --binding [All | None | Livermore | IA]  where options is number of electron binding energy correction type\n";
            std::cout << "-f or --filename [filename]  specify output filename, default validationTable.txt\n";
            args.exit = true;
            return args;
        } else if (arg == "-j" || arg == "--jobs") {
            i++;
            arg = args_str[i];
            if (std::from_chars(arg.data(), arg.data() + arg.size(), args.N_threads).ec != std::errc {}) {
                std::cout << "Number of jobs must be specified, see --help\n";
                args.exit = true;
                return args;
            }
        } else if (arg == "-b" || arg == "--binding") {
            i++;
            arg = args_str[i];
            std::string arg_str(arg);
            std::transform(arg.begin(), arg.end(), arg_str.begin(), [](auto& c) { return std::tolower(c); });
            if (arg_str == "all")
                args.correction = Arguments::LECorrection::All;
            else if (arg_str == "ia")
                args.correction = Arguments::LECorrection::IA;
            else if (arg_str == "livermore")
                args.correction = Arguments::LECorrection::Livermore;
            else if (arg_str == "none")
                args.correction = Arguments::LECorrection::None;
            else {
                args.exit = true;
                std::cout << "Binding correction type must be specified, see --help\n";
                return args;
            }
        } else if (arg == "-f" || arg == "--filename") {
            i++;
            arg = args_str[i];
            if (!arg.empty()) {
                args.filename = arg;
            } else {
                args.exit = true;
                std::cout << "Filename must be specified, see --help\n";
                return args;
            }
        }

        i++;
    }
    return args;
}

void printStart(const Arguments& args)
{
    std::cout << "Validation run of dxmclib\n";
    std::cout << "Number of threads: " << args.N_threads << std::endl;
    std::cout << "Binding correction: " << args.correctionString() << std::endl;
    std::cout << "Output file: " << args.filename << std::endl;

    ResultPrint resPrint(args.filename);
    resPrint.header();
}

int main(int argc, char* argv[])
{
    auto args = argparse(argc, argv);

    if (args.exit) {
        return EXIT_SUCCESS;
    }

    if (args.N_threads == 0)
        args.N_threads = std::thread::hardware_concurrency();

    auto success = true;
    printStart(args);

    if (args.correction == Arguments::LECorrection::None) {
        success = runAll<0>(args.N_threads);
    } else if (args.correction == Arguments::LECorrection::Livermore) {
        success = runAll<1>(args.N_threads);
    } else if (args.correction == Arguments::LECorrection::IA) {
        success = runAll<2>(args.N_threads);
    } else {
        success = runAll<0>(args.N_threads);
        success = runAll<1>(args.N_threads);
        success = runAll<2>(args.N_threads);
    }

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
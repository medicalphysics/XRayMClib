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

struct ResultKeys {
    std::string rCase = "unknown";
    std::string volume = "unknown";
    std::string specter = "unknown";
    std::string modus = "unknown";
    std::string model = "unknown";
    double TG195Result = 0;
    double TG195Result_min = 0;
    double TG195Result_max = 0;
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
        m_myfile << "Case, Volume, Specter, Model, Mode, TG195Result_min, TG195Result, TG195Result_max, Result, Stddev, nEvents, SimulationTime\n";
    }

    void operator()(const ResultKeys& r, bool terminal = true, const std::string& units = "ev/hist")
    {
        print(r, terminal, units);
    }

    void print(const ResultKeys& r, bool terminal = true, const std::string& units = "ev/hist")
    {
        m_myfile << r.rCase << ", ";
        m_myfile << r.volume << ", ";
        m_myfile << r.specter << ", ";
        m_myfile << r.model << ", ";
        m_myfile << r.modus << ", ";
        m_myfile << r.TG195Result_min << ", ";
        m_myfile << r.TG195Result << ", ";
        m_myfile << r.TG195Result_max << ", ";
        m_myfile << r.result << ", ";
        m_myfile << r.result_std << ", ";
        m_myfile << r.nEvents << ", ";
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
    return std::make_pair(1, w);
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
    std::array<double, 3> TG195result_min = { 0, 0, 0 };
    std::array<double, 3> TG195result_max = { 0, 0, 0 };
    if constexpr (std::same_as<Beam, IsotropicBeamCircle<>>) { // we have a specter
        res.modus = "polyenergetic";
        if (mammo) {
            const auto specter = TG195_30KV();
            beam.setEnergySpecter(specter);
            res.specter = "30kVp";
            filter_thickness = { 0.03431, 0.07663 };
            TG195result = { 1.858E-01, 9.282E-02, 4.642E-02 };
            TG195result_min = { 1.857E-01, 9.282E-02, 4.64E-02 };
            TG195result_max = { 1.860E-01, 9.29E-02, 4.65E-02 };
        } else {
            const auto specter = TG195_100KV();
            beam.setEnergySpecter(specter);
            res.specter = "100kVp";
            filter_thickness = { 0.3950, 0.9840 };
            TG195result = { 3.452E-02, 1.725E-02, 8.621E-03 };
            TG195result_min = { 3.452E-02, 1.725E-02, 8.621E-03 };
            TG195result_max = { 3.452E-02, 1.725E-02, 8.621E-03 };
        }
    } else {
        res.modus = "monoenergetic";
        if (mammo) {
            beam.setEnergy(29.99);
            res.specter = "30keV";
            filter_thickness = { 0.2273, 0.4546 };
            TG195result = { 5.733E-02, 2.868E-02, 1.436E-02 };
            TG195result_min = { 5.733E-02, 2.868E-02, 1.436E-02 };
            TG195result_max = { 5.733E-02, 2.868E-02, 1.436E-02 };
        } else {
            beam.setEnergy(99.99);
            res.specter = "100keV";
            filter_thickness = { 1.511, 3.022 };
            TG195result = { 2.902E-02, 1.448E-02, 7.227E-03 };
            TG195result_min = { 2.902E-02, 1.448E-02, 7.227E-03 };
            TG195result_max = { 2.902E-02, 1.448E-02, 7.227E-03 };
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
    res.TG195Result_min = TG195result_min[0];
    res.TG195Result_max = TG195result_max[0];
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
    res.TG195Result_min = TG195result_min[1];
    res.TG195Result_max = TG195result_max[1];
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
    res.TG195Result_min = TG195result_min[2];
    res.TG195Result_max = TG195result_max[2];
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
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "120kVp" : "56.4keV";
    res.rCase = "Case 2";
    res.modus = tomo ? "tomosyntesis" : "radiography";
    if (LOWENERGYCORRECTION == 0)
        res.model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        res.model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        res.model = "IA";

    std::cout << "TG195 Case 2 for " << res.modus << " orientation and " << res.specter << " photons with low en model: " << res.model << std::endl;

    double TG195_value, TG195_value_min, TG195_value_max;
    std::array<double, 9> TG195_voi_values, TG195_voi_values_min, TG195_voi_values_max;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        if (tomo) {
            TG195_value = 30923.13;
            TG195_value_min = 30896.40;
            TG195_value_max = 30967.00;
            TG195_voi_values = { 30.35, 23.52, 31.64, 23.52, 8.90, 70.53, 47.74, 20.31, 12.51 };
            TG195_voi_values_min = { 30.27, 23.44, 31.31, 23.44, 8.85, 70.21, 46.87, 20.21, 12.37 };
            TG195_voi_values_max = { 30.43, 23.57, 31.82, 23.57, 8.96, 71.28, 48.18, 20.41, 12.79 };
        } else {
            TG195_value = 33125.98;
            TG195_value_min = 33078.00;
            TG195_value_max = 33189.00;
            TG195_voi_values = { 24.97, 24.95, 33.52, 24.96, 24.97, 72.70, 49.99, 21.73, 13.48 };
            TG195_voi_values_min = { 24.90, 24.88, 33.18, 24.89, 24.89, 72.33, 49.09, 21.64, 13.37 };
            TG195_voi_values_max = { 25.03, 25.03, 33.74, 25.03, 25.02, 73.55, 50.46, 21.84, 13.75 };
        }
    } else {
        if (tomo) {
            TG195_value = 30883.83;
            TG195_value_min = 30870.00;
            TG195_value_max = 30912.80;
            TG195_voi_values = { 33.08, 25.48, 34.63, 25.51, 9.79, 70.80, 51.06, 22.28, 13.54 };
            TG195_voi_values_min = { 33.02, 25.45, 34.48, 25.45, 9.76, 70.51, 50.36, 22.25, 13.44 };
            TG195_voi_values_max = { 33.14, 25.54, 34.71, 25.56, 9.83, 71.50, 51.33, 22.33, 13.76 };
        } else {
            TG195_value = 33171.40;
            TG195_value_min = 33134.90;
            TG195_value_max = 33205.40;
            TG195_voi_values = { 27.01, 27.00, 36.67, 27.01, 27.01, 72.86, 53.35, 23.83, 14.60 };
            TG195_voi_values_min = { 26.98, 26.94, 36.44, 26.97, 26.96, 72.44, 52.62, 23.80, 14.49 };
            TG195_voi_values_max = { 27.01, 27.08, 36.76, 27.08, 27.09, 73.68, 53.63, 23.89, 14.80 };
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
    res.TG195Result_min = TG195_value_min;
    res.TG195Result_max = TG195_value_max;
    res.result = total_ev / (total_hist / 1000);
    res.result_std = std::sqrt(total_ev_var) / (total_hist / 1000);
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
        res.TG195Result_min = TG195_voi_values_min[ind - 1];
        res.TG195Result_max = TG195_voi_values_max[ind - 1];
        res.result /= (total_hist / 1000);
        res.result_std = std::sqrt(res.result_std) / (total_hist / 1000);
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
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "30kVp" : "16.8keV";
    res.modus = tomo ? "tomosyntesis" : "radiography";
    std::string model = "NoneLC";
    if (LOWENERGYCORRECTION == 1)
        model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        model = "IA";

    std::cout << "TG195 Case 3 for " << res.modus << " orientation and " << res.specter << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 512;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 1000000 : 1000000;

    using Box = WorldBox<NShells, LOWENERGYCORRECTION>;
    using Breast = TG195World3Breast<NShells, LOWENERGYCORRECTION>;
    using World = World<Box, Breast>;
    using Material = Material<NShells>;

    auto [water_d, water_w] = TG195_water();
    auto [pmma_d, pmma_w] = TG195_pmma();
    auto [air_d, air_w] = TG195_air();
    auto [breast_d, breast_w] = TG195_breast_tissue();
    auto [skin_d, skin_w] = TG195_skin();

    auto water = Material::byWeight(water_w).value();
    auto pmma = Material::byWeight(pmma_w).value();
    auto air = Material::byWeight(air_w).value();
    auto breasttissue = Material::byWeight(breast_w).value();
    auto skin = Material::byWeight(skin_w).value();

    World world;
    world.setMaterialByWeight(air_w, air_d);
    world.reserveNumberOfItems(5);

    auto& body = world.template addItem<Box>({ { -17, -15, -15, 0, 15, 15 } });
    body.setMaterial(water, water_d);
    auto& uplate = world.template addItem<Box>({ { 0, -13, 2.5, 14, 13, 2.7 } });
    uplate.setMaterial(pmma, pmma_d);
    auto& lplate = world.template addItem<Box>({ { 0, -13, -2.7, 14, 13, -2.5 } });
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
    std::array<double, 7> sim_subvol, sim_subvol_min, sim_subvol_max;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        // specter
        if (tomo) {
            sim_ev = 4188.833;
            sim_ev_min = 4176.5;
            sim_ev_max = 4198.2;
            sim_subvol = { 14.390, 15.825, 15.972, 15.445, 17.171, 5.619, 49.022 };
            sim_subvol_min = { 14.346, 15.786, 15.944, 15.402, 17.124, 5.597, 48.883 };
            sim_subvol_max = { 14.432, 15.875, 16.007, 15.490, 17.219, 5.646, 49.154 };
        } else {
            sim_ev = 4293.433;
            sim_ev_min = 4287.4;
            sim_ev_max = 4296.8;
            sim_subvol = { 16.502, 16.658, 16.814, 16.249, 16.521, 6.041, 50.041 };
            sim_subvol_min = { 16.471, 16.635, 16.779, 16.227, 16.475, 6.023, 49.981 };
            sim_subvol_max = { 16.528, 16.691, 16.839, 16.271, 16.548, 6.051, 50.075 };
        }
    } else {
        if (tomo) {
            sim_ev = 4577.743;
            sim_ev_min = 4565.2;
            sim_ev_max = 4588.2;
            sim_subvol = { 15.217, 16.836, 16.943, 16.431, 18.370, 5.043, 54.974 };
            sim_subvol_min = { 14.965, 16.574, 16.651, 16.158, 18.031, 4.970, 53.898 };
            sim_subvol_max = { 15.351, 16.994, 17.074, 16.569, 18.595, 5.086, 55.480 };
        } else {
            sim_ev = 4697.333;
            sim_ev_min = 4691.0;
            sim_ev_max = 4700.5;
            sim_subvol = { 17.692, 18.070, 17.865, 17.262, 17.768, 5.417, 56.017 };
            sim_subvol_min = { 17.462, 17.931, 17.303, 16.472, 17.727, 4.965, 55.232 };
            sim_subvol_max = { 17.794, 18.329, 18.074, 17.550, 17.801, 5.587, 56.299 };
        }
    }

    res.model = model;

    constexpr auto evNormal = 1000.0 / (N_HISTORIES * N_EXPOSURES);

    res.volume = "Total body";
    res.TG195Result = sim_ev;
    res.TG195Result_max = sim_ev_max;
    res.TG195Result_min = sim_ev_min;
    res.result = breast.energyScored(8).energyImparted() * evNormal;
    res.result_std = breast.energyScored(8).standardDeviation() * evNormal;
    res.nEvents = breast.energyScored(8).numberOfEvents();
    res.nMilliseconds = time_elapsed.count();
    print(res, true);

    for (int i = 0; i < 7; ++i) {
        res.volume = "VOI " + std::to_string(i + 1);
        res.TG195Result = sim_subvol[i];
        res.TG195Result_min = sim_subvol_min[i];
        res.TG195Result_max = sim_subvol_max[i];
        res.result = breast.energyScored(i).energyImparted() * evNormal;
        res.result_std = breast.energyScored(i).standardDeviation() * evNormal;
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
        std::cout << "120kVp";
    else
        std::cout << "56.4keV";
    std::cout << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 512;
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
        ev_history_var[i] = std::sqrt(ev_history_var[i]) / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
    }

    ResultPrint print;
    ResultKeys res;
    res.specter = specter ? "120kVp" : "56.4keV";
    res.rCase = "Case 4.1";
    res.model = model;
    res.modus = large_collimation ? "80mm collimation" : "10mm collimation";
    res.nMilliseconds = time_elapsed.count();

    std::array<double, 4> tg195, tg195_min, tg195_max;
    if (specter && large_collimation) {
        tg195 = { 3586.59, 3537.84, 3378.99, 2672.21 };
        tg195_min = { 3576.06, 3518.80, 3359.52, 2640.62 };
        tg195_max = { 3596.06, 3554.34, 3394.30, 2687.50 };
    }
    if (specter && !large_collimation) {
        tg195 = { 13137.02, 2585.47, 1706.86, 1250.61 };
        tg195_min = { 13115.77, 2568.81, 1702.05, 1244.18 };
        tg195_max = { 13159.30, 2594.24, 1712.09, 1255.70 };
    }
    if (!specter && large_collimation) {
        tg195 = { 3380.39, 3332.64, 3176.44, 2559.58 };
        tg195_min = { 3370.58, 3326.68, 3159.20, 2531.15 };
        tg195_max = { 3388.13, 3335.66, 3186.53, 2573.84 };
    }
    if (!specter && !large_collimation) {
        tg195 = { 11592.27, 2576.72, 1766.85, 1330.53 };
        tg195_min = { 11581.99, 2564.74, 1759.77, 1329.80 };
        tg195_max = { 11618.30, 2581.85, 1774.79, 1331.27 };
    }

    for (std::size_t i = 0; i < voi_locations.size(); ++i) {
        res.volume = "VOI " + std::to_string(i + 1);
        res.TG195Result = tg195[i];
        res.TG195Result_max = tg195_max[i];
        res.TG195Result_min = tg195_min[i];
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
        std::cout << "120kVp";
    else
        std::cout << "56.4keV";
    std::cout << " photons with low en model: " << model << std::endl;

    const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? 24 : 1024;
    const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 10000 : 1000000;

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
    res.specter = std::same_as<Beam, IsotropicBeam<>> ? "120kVp" : "56.4keV";

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
    std::array<double, 36> sim_ev_center, sim_ev_center_max, sim_ev_center_min, sim_ev_pher, sim_ev_pher_min, sim_ev_pher_max;
    if constexpr (std::same_as<Beam, IsotropicBeam<>>) {
        if (large_collimation) {
            sim_ev_center = { 10.878025, 10.9243, 10.884625, 10.89795, 10.87265, 10.902675, 10.8994, 10.880875, 10.875475, 10.8862, 10.895975, 10.88105, 10.8996, 10.886225, 10.8934, 10.8942, 10.879025, 10.8855, 10.894125, 10.8898, 10.8916, 10.895875, 10.889525, 10.889775, 10.89365, 10.901875, 10.894475, 10.906975, 10.888025, 10.877475, 10.883325, 10.875925, 10.8881, 10.886775, 10.88975, 10.900075 };
            sim_ev_center_min = { 10.79, 10.87, 10.77, 10.85, 10.77, 10.87, 10.85, 10.80, 10.82, 10.80, 10.86, 10.82, 10.87, 10.83, 10.85, 10.83, 10.82, 10.84, 10.84, 10.83, 10.83, 10.84, 10.84, 10.83, 10.84, 10.86, 10.83, 10.88, 10.84, 10.80, 10.84, 10.82, 10.82, 10.78, 10.83, 10.83 };
            sim_ev_center_max = { 10.95, 10.96, 10.98, 10.96, 10.94, 10.95, 10.96, 10.94, 10.94, 10.95, 10.94, 10.94, 10.94, 10.97, 10.95, 10.96, 10.95, 10.92, 10.95, 10.97, 10.95, 10.98, 10.95, 10.96, 10.93, 10.97, 10.95, 10.95, 10.93, 10.96, 10.94, 10.94, 10.95, 10.98, 10.94, 10.96 };
            sim_ev_pher = { 115.34325, 113.76275, 109.16925, 101.706, 91.562975, 78.39105, 61.388325, 40.08625, 22.471075, 11.781725, 6.14551, 3.42218, 2.05605, 1.35319, 0.96088275, 0.743808, 0.61922025, 0.55457575, 0.5309405, 0.55428325, 0.6219885, 0.74445025, 0.96480125, 1.3481875, 2.0611025, 3.4154, 6.15532, 11.7854, 22.461525, 40.13715, 61.42595, 78.328975, 91.481375, 101.61325, 109.10425, 113.8365 };
            sim_ev_pher_min = { 115.22, 113.59, 109, 101.62, 91.4, 78.32, 61.34, 39.98, 22.43, 11.7, 6.12, 3.41, 2.04, 1.35, 0.95, 0.74, 0.61, 0.55, 0.53, 0.55, 0.61, 0.73, 0.96, 1.34, 2.04, 3.4, 6.14, 11.73, 22.36, 40.08, 61.36, 78.13, 91.4, 101.54, 108.92, 113.64 };
            sim_ev_pher_max = { 115.53, 113.99, 109.35, 101.89, 91.91, 78.49, 61.46, 40.15, 22.52, 11.87, 6.16, 3.43, 2.06, 1.36, 0.97, 0.75, 0.63, 0.56, 0.53, 0.56, 0.63, 0.75, 0.97, 1.37, 2.08, 3.44, 6.17, 11.86, 22.56, 40.18, 61.58, 78.41, 91.59, 101.66, 109.23, 114.13 };
        } else {
            sim_ev_center = { 11.35475, 11.399475, 11.38755, 11.402175, 11.400225, 11.370625, 11.402625, 11.37715, 11.385375, 11.4096, 11.399825, 11.376675, 11.3698, 11.3613, 11.38865, 11.377925, 11.3692, 11.371525, 11.3807, 11.36, 11.37645, 11.379075, 11.379975, 11.368975, 11.377675, 11.38375, 11.3871, 11.3844, 11.37395, 11.3821, 11.371825, 11.395575, 11.379075, 11.372125, 11.4001, 11.39895 };
            sim_ev_center_min = { 11.27, 11.35, 11.29, 11.34, 11.37, 11.28, 11.37, 11.29, 11.33, 11.38, 11.38, 11.30, 11.30, 11.26, 11.34, 11.31, 11.31, 11.32, 11.31, 11.30, 11.32, 11.34, 11.31, 11.29, 11.32, 11.31, 11.35, 11.28, 11.26, 11.28, 11.28, 11.35, 11.32, 11.26, 11.35, 11.33 };
            sim_ev_center_max = { 11.43, 11.45, 11.45, 11.45, 11.46, 11.43, 11.43, 11.41, 11.42, 11.44, 11.42, 11.45, 11.41, 11.41, 11.44, 11.42, 11.44, 11.46, 11.45, 11.41, 11.43, 11.44, 11.47, 11.44, 11.46, 11.44, 11.43, 11.44, 11.43, 11.44, 11.42, 11.42, 11.40, 11.42, 11.45, 11.45 };
            sim_ev_pher = { 116.7915, 115.31075, 110.532, 103.1535, 92.9892, 79.7158, 62.57135, 40.92065, 22.966025, 12.169, 6.4216075, 3.5657225, 2.1609475, 1.418235, 1.014725, 0.780394, 0.64618, 0.5763575, 0.559597, 0.5776625, 0.648646, 0.7762825, 1.008508, 1.412125, 2.1584, 3.577945, 6.4296425, 12.1714, 23.039025, 40.926825, 62.451325, 79.759575, 93.074725, 103.269, 110.595, 115.27 };
            sim_ev_pher_min = { 116.7, 115.25, 110.33, 102.99, 92.92, 79.49, 62.51, 40.67, 22.77, 12.09, 6.39, 3.55, 2.13, 1.41, 1.01, 0.77, 0.64, 0.57, 0.55, 0.57, 0.64, 0.77, 1, 1.39, 2.14, 3.56, 6.42, 12.07, 22.9, 40.69, 62.35, 79.65, 92.78, 103, 110.31, 115.14 };
            sim_ev_pher_max = { 116.94, 115.4, 110.64, 103.46, 93.14, 80.03, 62.65, 41.06, 23.14, 12.27, 6.45, 3.59, 2.18, 1.42, 1.02, 0.79, 0.65, 0.58, 0.57, 0.58, 0.66, 0.79, 1.02, 1.42, 2.17, 3.61, 6.46, 12.26, 23.13, 41.04, 62.55, 79.85, 93.45, 103.6, 110.98, 115.42 };
        }
    } else {
        if (large_collimation) {
            sim_ev_center = { 11.630625, 11.632925, 11.617175, 11.624825, 11.633075, 11.600225, 11.61655, 11.6235, 11.592875, 11.6258, 11.612, 11.612775, 11.608675, 11.623975, 11.611325, 11.6174, 11.6234, 11.627975, 11.60745, 11.632875, 11.628275, 11.6239, 11.61645, 11.617375, 11.621775, 11.6178, 11.6444, 11.61515, 11.626375, 11.64605, 11.63335, 11.628425, 11.622, 11.6198, 11.59835, 11.609925 };
            sim_ev_center_min = { 11.61, 11.60, 11.55, 11.60, 11.61, 11.55, 11.57, 11.59, 11.51, 11.58, 11.57, 11.60, 11.59, 11.62, 11.60, 11.59, 11.60, 11.60, 11.57, 11.61, 11.61, 11.61, 11.59, 11.60, 11.61, 11.59, 11.62, 11.58, 11.59, 11.60, 11.61, 11.62, 11.61, 11.61, 11.56, 11.54 };
            sim_ev_center_max = { 11.64, 11.68, 11.66, 11.65, 11.66, 11.63, 11.66, 11.65, 11.65, 11.68, 11.65, 11.63, 11.63, 11.64, 11.63, 11.64, 11.64, 11.68, 11.64, 11.66, 11.67, 11.65, 11.64, 11.63, 11.64, 11.63, 11.66, 11.64, 11.67, 11.68, 11.65, 11.65, 11.64, 11.64, 11.62, 11.65 };
            sim_ev_pher = { 99.665175, 98.300175, 94.509325, 88.48595, 80.1113, 69.261125, 55.124425, 37.2351, 21.754, 11.767, 6.269845, 3.460205, 2.073845, 1.3435025, 0.94542875, 0.71714775, 0.58643225, 0.52293525, 0.49996925, 0.5225535, 0.5875545, 0.719903, 0.94029425, 1.3461175, 2.07283, 3.4740625, 6.2371725, 11.79715, 21.73405, 37.28175, 55.1853, 69.2588, 80.036275, 88.397, 94.640025, 98.332825 };
            sim_ev_pher_min = { 99.57, 98.13, 94.39, 88.29, 79.95, 69.17, 55, 37.18, 21.71, 11.67, 6.23, 3.42, 2.06, 1.33, 0.93, 0.71, 0.58, 0.52, 0.5, 0.52, 0.58, 0.71, 0.93, 1.33, 2.06, 3.46, 6.21, 11.75, 21.66, 37.25, 55.14, 69.19, 79.95, 88.36, 94.51, 98.28 };
            sim_ev_pher_max = { 99.8, 98.44, 94.57, 88.74, 80.36, 69.32, 55.18, 37.29, 21.80, 11.83, 6.29, 3.49, 2.09, 1.35, 0.96, 0.72, 0.59, 0.53, 0.5, 0.53, 0.59, 0.73, 0.95, 1.36, 2.09, 3.49, 6.26, 11.83, 21.79, 37.3, 55.24, 69.33, 80.09, 88.47, 94.79, 98.41 };
        } else {
            sim_ev_center = { 12.168675, 12.11255, 12.163, 12.1358, 12.0881, 12.10275, 12.135925, 12.12425, 12.14915, 12.149575, 12.1482, 12.15955, 12.130775, 12.148475, 12.153925, 12.1383, 12.13475, 12.148675, 12.1576, 12.1452, 12.15665, 12.16225, 12.156625, 12.156875, 12.13905, 12.1561, 12.159325, 12.137025, 12.149725, 12.09575, 12.1471, 12.1231, 12.1211, 12.1368, 12.134825, 12.13955 };
            sim_ev_center_min = { 12.13, 12.08, 12.12, 12.12, 12.00, 12.02, 12.11, 12.11, 12.11, 12.12, 12.13, 12.13, 12.10, 12.10, 12.12, 12.12, 12.10, 12.12, 12.13, 12.11, 12.13, 12.13, 12.13, 12.12, 12.10, 12.12, 12.12, 12.11, 12.12, 12.02, 12.13, 12.10, 12.09, 12.10, 12.12, 12.12 };
            sim_ev_center_max = { 12.24, 12.15, 12.20, 12.15, 12.14, 12.16, 12.20, 12.14, 12.21, 12.18, 12.17, 12.21, 12.17, 12.20, 12.19, 12.18, 12.19, 12.18, 12.19, 12.19, 12.19, 12.19, 12.19, 12.20, 12.19, 12.19, 12.20, 12.18, 12.19, 12.14, 12.18, 12.14, 12.17, 12.19, 12.16, 12.18 };
            sim_ev_pher = { 101.29375, 99.802475, 96.108725, 89.966275, 81.369475, 70.5424, 56.2835, 38.1283, 22.348775, 12.166625, 6.515755, 3.643485, 2.180365, 1.41379, 0.984371, 0.7510235, 0.6165785, 0.54522075, 0.52163, 0.54515775, 0.61643625, 0.750186, 0.9914865, 1.4138375, 2.174915, 3.6299975, 6.5142425, 12.200725, 22.2865, 38.146275, 56.287325, 70.64165, 81.424625, 89.82015, 96.16995, 99.92435 };
            sim_ev_pher_min = { 101.17, 99.69, 95.98, 89.89, 81.24, 70.38, 56.08, 37.96, 22.30, 12.06, 6.49, 3.63, 2.17, 1.41, 0.98, 0.75, 0.61, 0.54, 0.52, 0.54, 0.61, 0.74, 0.98, 1.40, 2.16, 3.62, 6.49, 12.17, 22.09, 38.11, 56.11, 70.60, 81.08, 89.64, 96.05, 99.79 };
            sim_ev_pher_max = { 101.39, 99.95, 96.29, 90.05, 81.51, 70.68, 56.41, 38.23, 22.42, 12.28, 6.54, 3.66, 2.20, 1.43, 0.99, 0.76, 0.63, 0.55, 0.52, 0.55, 0.63, 0.76, 1.01, 1.43, 2.19, 3.66, 6.54, 12.25, 22.41, 38.17, 56.40, 70.70, 81.68, 89.97, 96.43, 99.99 };
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
        res.TG195Result_min = sim_ev_pher_min[i];
        res.TG195Result_max = sim_ev_pher_max[i];
        res.result = cylinder.energyScoredPeriferyCylinder().energyImparted() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        res.result_std = cylinder.energyScoredPeriferyCylinder().standardDeviation() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        res.nEvents = cylinder.energyScoredPeriferyCylinder().numberOfEvents();
        std::cout << " Pherifery: " << res.result << " sim/TG195: [" << (res.result / sim_ev_pher[i] - 1) * 100 << "%]";
        print(res, false);
        res.modus = large_collimation ? "Center 80mm collimation" : "Center 10mm collimation";
        res.volume = std::to_string(angInt);
        res.TG195Result = sim_ev_center[i];
        res.TG195Result_min = sim_ev_center_min[i];
        res.TG195Result_max = sim_ev_center_max[i];
        res.result = cylinder.energyScoredCenterCylinder().energyImparted() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
        res.result_std = cylinder.energyScoredCenterCylinder().standardDeviation() / ((N_HISTORIES * N_EXPOSURES) / 1000.0);
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
    std::vector<std::string> matFormula = {
        "C0.015019N78.443071O21.074800Ar0.467110",
        "H51.869709C36.108118N4.019974O8.002200",
        "C100.000000",
        "H63.000070C12.890598N1.165843O22.756479Na0.026307P0.039052S0.056594Cl0.034118K0.030937",
        "H63.688796C7.143744N1.278063O27.701991Na0.026851P0.039859S0.038509Cl0.034823K0.047365",
        "H63.731478C5.452396N1.380394O29.198156Na0.054259P0.040273S0.058363Cl0.052777K0.031904",
        "H63.217465C7.229913N1.338082O27.958043Na0.054349P0.060510S0.058460Cl0.035243K0.047936",
        "H63.000070C12.890598N1.165843O22.756479Na0.026307P0.039052S0.056594Cl0.034118K0.030937",
        "H63.655092C5.860784N1.423215O28.851671Na0.027097P0.060337S0.038862Cl0.035143K0.047799",
        "H64.343953C5.858428N0.961057O28.720940Na0.026615P0.019755S0.019085Cl0.034518K0.015650",
        "H64.343953C5.858428N0.961057O28.720940Na0.026615P0.019755S0.019085Cl0.034518K0.015650",
        "H63.939187C8.555183N0.955012O26.374094Na0.052895P0.039261S0.018965Cl0.034300K0.031102",
        "H63.000070C12.890598N1.165843O22.756479Na0.026307P0.039052S0.056594Cl0.034118K0.030937",
        "H63.845575C6.130922N1.060311O28.814466Na0.053834P0.019979S0.019302Cl0.034909K0.015827I0.004876",
        "H63.000070C12.890598N1.165843O22.756479Na0.026307P0.039052S0.056594Cl0.034118K0.030937",
        "H64.343953C5.858428N0.961057O28.720940Na0.026615P0.019755S0.019085Cl0.034518K0.015650",
        "H64.343953C5.858428N0.961057O28.720940Na0.026615P0.019755S0.019085Cl0.034518K0.015650",
        "H62.083429C10.628873N1.876505O25.228547Na0.054442P0.020204S0.039039Cl0.052955K0.016006",
        "H61.873627C28.698524N0.675867O8.736110P0.004495S0.004342K0.003561Ca0.003473",
        "H39.229963C15.009010N3.487490O31.621690Na0.050590Mg0.095705P3.867606S0.108832Ca6.529115",
    };

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

    const std::array<std::size_t, 3> dim = { 500, 320, 260 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    const std::array<double, 3> spacing = { 0.1f, 0.1f, 0.1f };

    AAVoxelGrid<NMATSHELLS, LOWENERGYCORRECTION, TRANSPARENTVOXEL> grid;

    auto matArray = readBinaryArray<std::uint8_t>("case5world.bin", size);

    std::vector<std::pair<double, std::string>> matInfo;
    matInfo.push_back(std::make_pair(.001205, "Air"));
    matInfo.push_back(std::make_pair(.075, "Cushion Foam"));
    matInfo.push_back(std::make_pair(1.20, "Carbon fiber"));
    matInfo.push_back(std::make_pair(1.03, "Soft Tissue"));
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
    matInfo.push_back(std::make_pair(1.92, "Cortial Bone"));

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

    using World = World<AAVoxelGrid<NShells, LOWENERGYCORRECTION, 255>>;
    auto [grid_object, matInf] = generateTG195World5<NShells, LOWENERGYCORRECTION, 255>();

    World world;
    auto& grid = world.addItem(grid_object);
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
        res.specter = "120kVp";
    } else {
        beam.setEnergy(56.4);
        res.specter = "56.4keV";
    }

    std::array<std::array<double, 17>, 8> tg195_doses;
    std::array<std::array<double, 17 * 2>, 8> tg195_doses_min_max;
    if constexpr (std::same_as<B, IsotropicBeam<>>) {
        tg195_doses[0] = { 12374.98, 2917.75, 1275.86, 612.31, 5.78, 16.68, 121.04, 15.16, 8.17, 0.15, 1.65, 40.66, 9.78, 33.37, 559.77, 21.49, 7727.77 };
        tg195_doses[1] = { 12594.50, 1801.82, 1007.28, 612.42, 5.74, 9.97, 76.98, 8.73, 5.86, 0.13, 1.39, 33.38, 7.12, 30.52, 538.17, 17.76, 6631.55 };
        tg195_doses[2] = { 10648.03, 737.54, 640.75, 447.05, 3.70, 6.58, 36.07, 3.29, 3.19, 0.10, 1.03, 15.51, 3.38, 21.38, 348.68, 6.82, 4814.57 };
        tg195_doses[3] = { 10137.70, 730.79, 716.21, 389.11, 2.62, 17.39, 56.73, 5.20, 4.38, 0.14, 1.24, 11.31, 3.72, 24.87, 223.37, 2.19, 7437.66 };
        tg195_doses[4] = { 10250.38, 1211.35, 1043.24, 385.78, 2.57, 31.09, 102.26, 9.74, 6.91, 0.18, 1.51, 10.66, 5.45, 30.25, 143.22, 2.96, 9718.55 };
        tg195_doses[5] = { 10069.33, 1121.22, 687.39, 243.21, 1.62, 33.91, 107.84, 11.19, 6.92, 0.15, 1.33, 8.57, 5.34, 29.45, 232.17, 2.14, 7265.83 };
        tg195_doses[6] = { 10666.10, 1503.58, 601.06, 164.33, 1.32, 30.15, 122.86, 16.10, 7.16, 0.12, 1.20, 10.51, 6.82, 22.47, 365.63, 7.18, 5072.79 };
        tg195_doses[7] = { 12488.68, 2558.90, 877.00, 375.60, 3.76, 24.31, 138.64, 18.52, 8.43, 0.13, 1.49, 26.49, 9.60, 27.68, 562.67, 18.29, 6436.26 };
        tg195_doses_min_max[0] = { 12363.00, 12397.00, 2902.87, 2942.20, 1271.08, 1283.75, 608.27, 618.20, 5.76, 5.81, 16.52, 16.88, 120.34, 122.25, 15.05, 15.35, 8.09, 8.31, 0.15, 0.15, 1.64, 1.67, 39.95, 41.06, 9.68, 9.94, 33.23, 33.54, 557.44, 561.22, 21.24, 21.97, 7709.95, 7745.70 };
        tg195_doses_min_max[1] = { 12523.10, 12623.00, 1796.56, 1807.00, 1004.47, 1010.00, 609.99, 615.15, 5.67, 5.80, 9.93, 10.03, 76.56, 77.33, 8.68, 8.78, 5.84, 5.89, 0.13, 0.14, 1.37, 1.41, 32.54, 33.74, 7.09, 7.16, 30.43, 30.67, 535.91, 539.85, 17.57, 18.03, 6599.30, 6664.10 };
        tg195_doses_min_max[2] = { 10637.70, 10666.00, 731.15, 745.57, 637.07, 645.89, 444.01, 451.57, 3.65, 3.72, 6.52, 6.64, 35.77, 36.66, 3.25, 3.36, 3.16, 3.24, 0.10, 0.11, 1.03, 1.04, 15.41, 15.62, 3.34, 3.42, 21.28, 21.50, 347.40, 349.40, 6.77, 6.92, 4789.72, 4852.34 };
        tg195_doses_min_max[3] = { 10081.40, 10159.60, 728.03, 734.45, 714.54, 718.82, 386.95, 391.16, 2.60, 2.63, 17.29, 17.48, 56.44, 57.14, 5.16, 5.22, 4.34, 4.45, 0.14, 0.15, 1.23, 1.25, 11.16, 11.41, 3.70, 3.74, 24.79, 24.97, 222.30, 224.14, 2.18, 2.20, 7394.90, 7474.50 };
        tg195_doses_min_max[4] = { 10237.80, 10270.80, 1200.94, 1228.77, 1038.39, 1051.39, 382.82, 389.95, 2.54, 2.62, 30.87, 31.41, 101.33, 103.57, 9.67, 9.84, 6.83, 7.00, 0.17, 0.18, 1.49, 1.52, 10.62, 10.73, 5.38, 5.57, 29.96, 30.49, 143.13, 143.34, 2.94, 3.00, 9702.03, 9750.60 };
        tg195_doses_min_max[5] = { 10021.20, 10091.30, 1116.27, 1124.20, 685.28, 688.73, 241.89, 244.32, 1.59, 1.65, 33.78, 34.07, 107.21, 108.30, 11.11, 11.26, 6.89, 6.98, 0.15, 0.16, 1.31, 1.37, 8.51, 8.64, 5.32, 5.37, 29.37, 29.54, 230.96, 233.04, 2.12, 2.16, 7225.89, 7302.00 };
        tg195_doses_min_max[6] = { 10651.90, 10698.50, 1492.70, 1520.50, 597.92, 605.55, 163.20, 166.28, 1.30, 1.36, 29.91, 30.49, 121.78, 124.54, 15.91, 16.42, 7.10, 7.22, 0.11, 0.12, 1.19, 1.21, 10.47, 10.60, 6.73, 6.93, 22.33, 22.61, 363.87, 367.24, 7.15, 7.24, 5047.17, 5112.62 };
        tg195_doses_min_max[7] = { 12432.60, 12514.10, 2550.60, 2564.40, 874.46, 879.88, 373.98, 377.21, 3.75, 3.78, 24.17, 24.46, 138.00, 139.15, 18.41, 18.65, 8.40, 8.48, 0.13, 0.14, 1.47, 1.50, 26.10, 26.70, 9.54, 9.66, 27.58, 27.78, 560.67, 564.63, 18.22, 18.46, 6407.19, 6466.50 };
    } else {
        tg195_doses[0] = { 11574.28, 3086.42, 1301.17, 679.47, 6.37, 17.57, 134.40, 16.73, 8.79, 0.15, 1.74, 44.22, 10.72, 36.90, 456.36, 21.68, 8761.23 };
        tg195_doses[1] = { 11761.25, 1932.04, 1045.65, 683.63, 6.33, 9.91, 82.98, 9.22, 6.18, 0.13, 1.47, 36.39, 7.64, 33.46, 437.27, 17.97, 7669.28 };
        tg195_doses[2] = { 9975.20, 786.73, 679.42, 495.34, 3.94, 6.27, 36.21, 3.15, 3.15, 0.10, 1.06, 16.90, 3.38, 23.03, 285.87, 6.80, 5611.52 };
        tg195_doses[3] = { 9581.67, 765.35, 755.19, 421.94, 2.62, 18.05, 58.05, 5.07, 4.38, 0.15, 1.27, 12.03, 3.61, 26.36, 189.37, 2.07, 8510.07 };
        tg195_doses[4] = { 9704.91, 1265.14, 1085.65, 411.34, 2.52, 33.79, 109.31, 10.18, 7.15, 0.18, 1.57, 11.22, 5.51, 32.11, 129.81, 2.82, 11003.75 };
        tg195_doses[5] = { 9487.98, 1202.94, 721.26, 251.05, 1.52, 37.44, 117.36, 12.02, 7.31, 0.16, 1.38, 8.95, 5.48, 31.53, 195.79, 2.00, 8367.44 };
        tg195_doses[6] = { 9970.40, 1625.79, 636.47, 168.17, 1.25, 33.40, 136.69, 17.88, 7.68, 0.12, 1.26, 11.26, 7.29, 24.39, 298.96, 7.18, 5876.97 };
        tg195_doses[7] = { 11649.90, 2725.15, 916.40, 409.55, 3.99, 26.54, 155.30, 20.68, 9.13, 0.13, 1.55, 28.93, 10.49, 30.43, 455.47, 18.39, 7391.31 };
        tg195_doses_min_max[0] = { 11557.00, 11596.10, 3074.90, 3110.85, 1296.90, 1309.36, 675.87, 686.53, 6.33, 6.46, 17.5, 17.74, 133.61, 136.08, 16.61, 17, 8.74, 8.9, 0.15, 0.16, 1.72, 1.75, 43.57, 44.49, 10.63, 10.9, 36.69, 37.23, 455.85, 456.99, 21.51, 22.06, 8747.74, 8785.40 };
        tg195_doses_min_max[1] = { 11745.00, 11783.50, 1924.80, 1946.92, 1042.10, 1051.22, 679.88, 690.8, 6.25, 6.47, 9.88, 10.01, 82.59, 83.83, 9.09, 9.39, 6.14, 6.25, 0.13, 0.13, 1.46, 1.48, 35.88, 36.66, 7.61, 7.73, 33.35, 33.63, 436.83, 437.9, 17.78, 18.4, 7652.67, 7699.95 };
        tg195_doses_min_max[2] = { 9962.60, 9990.23, 780.87, 797.63, 676.34, 685.09, 492.74, 501.17, 3.92, 3.99, 6.22, 6.37, 35.94, 36.84, 3.13, 3.18, 3.13, 3.2, 0.1, 0.11, 1.04, 1.07, 16.79, 16.96, 3.37, 3.4, 22.9, 23.31, 285.51, 286.13, 6.77, 6.85, 5589.64, 5661.26 };
        tg195_doses_min_max[3] = { 9526.84, 9612.64, 763.23, 768.82, 753.96, 756.73, 420.54, 423.97, 2.61, 2.64, 17.93, 18.16, 57.95, 58.3, 5.06, 5.09, 4.38, 4.39, 0.14, 0.15, 1.26, 1.27, 11.84, 12.1, 3.59, 3.63, 26.31, 26.4, 188.26, 189.99, 2.05, 2.08, 8466.28, 8528.96 };
        tg195_doses_min_max[4] = { 9687.10, 9723.43, 1256.84, 1284.20, 1080.90, 1094.77, 408.68, 416, 2.48, 2.59, 33.6, 34.2, 108.77, 110.52, 10.1, 10.35, 7.13, 7.18, 0.18, 0.18, 1.56, 1.58, 11.16, 11.32, 5.46, 5.59, 31.96, 32.47, 129.65, 130, 2.77, 2.91, 10999.30, 11010.00 };
        tg195_doses_min_max[5] = { 9443.60, 9513.53, 1199.17, 1206.59, 719.85, 722.48, 250, 252.77, 1.49, 1.54, 37.32, 37.62, 117.02, 117.74, 12, 12.05, 7.27, 7.34, 0.16, 0.16, 1.37, 1.38, 8.86, 9.04, 5.47, 5.51, 31.42, 31.74, 195.1, 196.14, 1.98, 2.03, 8331.94, 8383.47 };
        tg195_doses_min_max[6] = { 9952.40, 10000.60, 1615.96, 1647.62, 633.85, 641.47, 167.26, 170.77, 1.22, 1.3, 33.23, 33.83, 136.08, 138.46, 17.73, 18.25, 7.62, 7.78, 0.12, 0.13, 1.25, 1.29, 11.24, 11.28, 7.25, 7.35, 24.27, 24.64, 298.53, 299.7, 7.1, 7.37, 5854.81, 5928.06 };
        tg195_doses_min_max[7] = { 11596.20, 11678.40, 2719.50, 2733.11, 914.9, 918.56, 408.33, 411.42, 3.98, 4.02, 26.46, 26.78, 154.59, 156.35, 20.52, 20.85, 9.06, 9.19, 0.13, 0.14, 1.54, 1.57, 28.48, 29.11, 10.43, 10.56, 30.39, 30.45, 453.2, 456.72, 18.29, 18.61, 7363.20, 7403.80 };
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
                res.TG195Result_min = tg195_doses_min_max[angInt][tg195_idx * 2];
                res.TG195Result_max = tg195_doses_min_max[angInt][tg195_idx * 2 + 1];
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
                res.result_std = std::sqrt(ei_var) / ((N_HISTORIES * N_EXPOSURES) / 1000);
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
    /*
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
    */
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
                std::cout << "Binding correction type must be specified, see --help\n";
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
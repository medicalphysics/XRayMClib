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

Copyright 2023 Erlend Andersen
*/

#include "xraymc/beams/filters/bowtiefilter.hpp"
#include "xraymc/beams/filters/ctorganaecfilter.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <iostream>

bool testBowtie()
{
    std::vector<double> angle, weight;

    weight = {
        9084.336977,
        9056.683016,
        9010.555199,
        8976.621955,
        8937.858092,
        8847.304732,
        8754.532826,
        8656.860432,
        8527.312593,
        8422.551773,
        8306.212116,
        8125.62487,
        7931.430865,
        7786.041916,
        7639.243664,
        7424.480585,
        7197.561479,
        7024.472266,
        6851.343852,
        6621.742936,
        6388.017234,
        6152.921949,
        5974.012763,
        5794.802125,
        5568.133743,
        5353.997039,
        5134.814809,
        4912.471357,
        4695.75262,
        4498.066906,
        4310.238528,
        4119.691564,
        3935.618728,
        3777.102292,
        3629.564305,
        3476.787516,
        3289.33132,
        3128.526529,
        3013.761207,
        2897.170244,
        2762.216686,
        2654.233847,
        2548.605629,
        2441.276609,
        2178.949877,
        2116.337359,
        2015.555302,
        1964.260525,
        1900.254578,
        1857.235929,
        1790.316289,
        1756.816535,
        1691.007325,
        1654.787224,
        1621.951184,
        1579.044149,
        1547.725401,
        1520.275338,
        1492.313526,
        1473.686923,
        1446.729063,
        1427.908945,
        1399.920027,
        1248.556091
    };
    angle = {
        0.002678189,
        0.008034566,
        0.013390943,
        0.018747321,
        0.024103698,
        0.029460075,
        0.034816452,
        0.04017283,
        0.045529207,
        0.050885584,
        0.056241962,
        0.061598339,
        0.066954716,
        0.072311093,
        0.077667471,
        0.083023848,
        0.088380225,
        0.093736603,
        0.09909298,
        0.104449357,
        0.109805734,
        0.115162112,
        0.120518489,
        0.125874866,
        0.131231244,
        0.136587621,
        0.141943998,
        0.147300375,
        0.152656753,
        0.15801313,
        0.163369507,
        0.168725885,
        0.174082262,
        0.179438639,
        0.184795016,
        0.190151394,
        0.195507771,
        0.200864148,
        0.206220526,
        0.211576903,
        0.21693328,
        0.222289657,
        0.227646035,
        0.233002412,
        0.238358789,
        0.243715167,
        0.249071544,
        0.254427921,
        0.259784299,
        0.265140676,
        0.270497053,
        0.27585343,
        0.281209808,
        0.286566185,
        0.291922562,
        0.29727894,
        0.302635317,
        0.307991694,
        0.313348071,
        0.318704449,
        0.324060826,
        0.329417203,
        0.334773581,
        0.340129958
    };

    xraymc::BowtieFilter filter(angle, weight);
    /*
        // test values
        std::cout << "a, f1, f1s, f2, f2s, f3, f3s\n";
        for (std::size_t i = 0; i < angle.size(); ++i) {
            const auto a = angle[i];
            std::cout << a << ", ";
            std::cout << weight[i] << ",";
            std::cout << filter(a) << '\n';
        }
    */

    double w = 0;
    constexpr std::size_t N = 1E6;
    xraymc::RandomState state;
    for (std::size_t i = 0; i < N; ++i) {
        const auto ang = state.randomUniform(0.340129958);
        w += filter(ang);
    }
    w /= N;

    auto success = std::abs(w - 1) < 0.01;
    std::cout << "Test Bowtiefilter ";
    if (success)
        std::cout << "SUCCESS\n";
    else
        std::cout << "FAILURE\n";
    return success;
}

bool testOrganAECfilter()
{

    xraymc::CTOrganAECFilter f;

    const auto deg45 = std::numbers::pi_v<double> / 4;
    const auto deg10 = xraymc::DEG_TO_RAD<double>() * 10;
    f.setLowWeightFactor(0.25);
    f.setStartStopAngles(-deg45, deg45);
    f.setRampAngle(deg10);

    constexpr std::size_t N = 1E6;

    double acc = 0;
    xraymc::RandomState state;
    for (std::size_t i = 0; i < N; ++i) {
        const auto ang = state.randomUniform(-std::numbers::pi_v<double> * 1, std ::numbers::pi_v<double> * 1);
        acc += f(ang);
    }
    const auto mean = acc / N;
    const auto max = f.maxWeight();

    bool success = std::abs(mean - 1) < 1E-3;

    std::cout << "Test CTOrganFilter ";
    if (success)
        std::cout << "SUCCESS\n";
    else
        std::cout << "FAILURE\n";

    if (success)
        return true;

    return false;
}

int main()
{

    bool success = true;
    success = success && testOrganAECfilter();
    success = success && testBowtie();
    if (success) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}
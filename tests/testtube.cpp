

#include <algorithm>
#include <iostream>
#include <numeric>

#include "dxmc/beams/tube/tube.hpp"

using namespace dxmc;

constexpr double DOUBLEERRF = 1E-6;

bool testHalfLayerCalculation()
{

    Tube t;
    t.setVoltage(100);
    t.setAlFiltration(2.0);
    t.setAnodeAngleDeg(10);
    auto e = t.getEnergy();
    auto s = t.getSpecter(e);
    const auto& al = AtomHandler::Atom(13);

    auto p = interpolate(al.photoel, e);
    auto i = interpolate(al.incoherent, e);
    auto c = interpolate(al.coherent, e);
    auto att = addVectors(p, i, c);

    const auto al_dens = al.standardDensity;

    std::transform(att.cbegin(), att.cend(), att.begin(), [al_dens](auto a) {
        return a * al_dens;
    });

    const auto mmHVL = t.mmAlHalfValueLayer();

    auto air = Material<5>::byNistName("Air, Dry (near sea level)").value();
    std::vector<double> kerma(e.size());
    double I0 = 0;
    double I1 = 0;
    for (std::size_t i = 0; i < e.size(); ++i) {
        I0 += e[i] * s[i] * air.massEnergyTransferAttenuation(e[i]);
        I1 += e[i] * s[i] * air.massEnergyTransferAttenuation(e[i]) * std::exp(-att[i] * mmHVL * .1);
    }

    bool success = (std::abs(I1 / I0) - 0.5) < 0.01;
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    return success;
}

void printSpecter()
{
    Tube t;
    t.setAnodeAngleDeg(10);
    t.setVoltage(100);
    const auto s = t.getSpecter();
    std::cout << "Energy [keV], Intensity [A.U]\n";
    for (const auto [e, n] : s) {
        std::cout << e << ", " << n << '\n';
    }
}

int main(int argc, char* argv[])
{
    std::cout << "Testing tube ";
    bool success = testHalfLayerCalculation();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

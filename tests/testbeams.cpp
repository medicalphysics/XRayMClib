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

#include "dxmc/beams/beamtype.hpp"
#include "dxmc/beams/cbctbeam.hpp"
#include "dxmc/beams/ctsequentialbeam.hpp"
#include "dxmc/beams/ctspiralbeam.hpp"
#include "dxmc/beams/ctspiraldualenergybeam.hpp"
#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/isotropicbeamcircle.hpp"
#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/beams/isotropicmonoenergybeamcircle.hpp"
#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/vectormath.hpp"

#include <iostream>

template <dxmc::BeamType B>
bool initiateBeam(B& beam)
{
    auto e = beam.exposure(0);
    dxmc::RandomState state;
    auto p = e.sampleParticle(state);
    return true;
}

bool testIsotropicBeamCircle()
{
    dxmc::IsotropicBeamCircle<> beam;

    constexpr double angle = 10.0 * std::numbers::pi_v<double> / 180.0;

    beam.setCollimationHalfAngle(angle);

    dxmc::Tube tube;
    auto specter = tube.getSpecter();
    beam.setEnergySpecter(specter);
    return initiateBeam(beam);
}

bool testIsotropicMonoEnergyBeamCircle()
{
    dxmc::IsotropicMonoEnergyBeamCircle<> beam;

    constexpr double angle = 10.0 * std::numbers::pi_v<double> / 180.0;

    beam.setCollimationHalfAngle(angle);

    return initiateBeam(beam);
}

bool testDXBeam()
{
    dxmc::DXBeam<> beam;
    auto& tube = beam.tube();
    return initiateBeam(beam);
}

bool testpencilbeam()
{
    dxmc::PencilBeam<> beam;
    auto e = beam.exposure(0);
    dxmc::RandomState state;
    auto p = e.sampleParticle(state);
    return initiateBeam(beam);
}

bool testIsotropicMonoEnergyBeam()
{
    using Beam = dxmc::IsotropicMonoEnergyBeam<>;
    Beam beam;
    return initiateBeam(beam);
}

bool testCTSpiralBeam()
{
    using Beam = dxmc::CTSpiralBeam<>;
    Beam beam;
    beam.setStartStopPosition({ 0, 0, 0 }, { 0, 0, 1 });
    beam.setNumberOfParticlesPerExposure(1E2);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);

    return initiateBeam(beam);
}

bool testCTSeqBeam()
{
    using Beam = dxmc::CTSequentialBeam<>;
    Beam beam;
    beam.setNumberOfParticlesPerExposure(1E2);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);

    return initiateBeam(beam);
}

bool testCTSpiralDEBeam()
{
    using Beam = dxmc::CTSpiralDualEnergyBeam<>;
    Beam beam;
    beam.setTubeAVoltage(140);
    beam.setTubeBVoltage(80);

    return initiateBeam(beam);
}

bool testCBCTBeam()
{
    using Beam = dxmc::CBCTBeam<>;
    Beam beam;
    beam.setTubeVoltage(140);
    return initiateBeam(beam);
}

int main()
{
    std::cout << "Testing beams ";

    bool success = true;
    success = success && testIsotropicMonoEnergyBeamCircle();
    success = success && testIsotropicBeamCircle();
    success = success && testDXBeam();
    success = success && testIsotropicMonoEnergyBeam();
    success = success && testpencilbeam();
    success = success && testCTSpiralBeam();
    success = success && testCBCTBeam();
    success = success && testCTSpiralDEBeam();

    if (success)
        std::cout << "SUCCESS\n";
    else
        std::cout << "FAILURE\n";

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
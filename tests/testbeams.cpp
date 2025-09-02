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

#include "xraymc/beams/beamtype.hpp"
#include "xraymc/beams/cbctbeam.hpp"
#include "xraymc/beams/ctsequentialbeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/ctspiraldualenergybeam.hpp"
#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/isotropicbeamcircle.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/constants.hpp"
#include "xraymc/vectormath.hpp"

#include <fstream>
#include <iostream>

template <xraymc::BeamType B>
bool initiateBeam(B& beam)
{
    auto e = beam.exposure(0);
    xraymc::RandomState state;
    auto p = e.sampleParticle(state);
    return true;
}

bool testIsotropicBeamCircle()
{
    xraymc::IsotropicBeamCircle<> beam;

    constexpr double angle = 10.0 * std::numbers::pi_v<double> / 180.0;

    beam.setCollimationHalfAngle(angle);

    xraymc::Tube tube;
    auto specter = tube.getSpecter();
    beam.setEnergySpecter(specter);
    return initiateBeam(beam);
}

bool testIsotropicMonoEnergyBeamCircle()
{
    xraymc::IsotropicMonoEnergyBeamCircle<> beam;

    constexpr double angle = 10.0 * std::numbers::pi_v<double> / 180.0;

    beam.setCollimationHalfAngle(angle);

    return initiateBeam(beam);
}

bool testDXBeam()
{
    xraymc::DXBeam<> beam;
    auto& tube = beam.tube();
    return initiateBeam(beam);
}

bool testpencilbeam()
{
    xraymc::PencilBeam<> beam;
    auto e = beam.exposure(0);
    xraymc::RandomState state;
    auto p = e.sampleParticle(state);
    return initiateBeam(beam);
}

bool testIsotropicMonoEnergyBeam()
{
    using Beam = xraymc::IsotropicMonoEnergyBeam<>;
    Beam beam;
    return initiateBeam(beam);
}

bool testCTSpiralBeam()
{
    using Beam = xraymc::CTSpiralBeam<>;
    Beam beam;
    beam.setStartStopPosition({ 0, 0, 0 }, { 0, 0, 1 });
    beam.setNumberOfParticlesPerExposure(1E2);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);

    return initiateBeam(beam);
}

bool testCTSeqBeam()
{
    using Beam = xraymc::CTSequentialBeam<>;
    Beam beam;
    beam.setNumberOfParticlesPerExposure(1E2);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);

    return initiateBeam(beam);
}

bool testCTSpiralDEBeam()
{
    using Beam = xraymc::CTSpiralDualEnergyBeam<>;
    Beam beam;
    beam.setTubeAVoltage(140);
    beam.setTubeBVoltage(80);

    return initiateBeam(beam);
}

bool testCBCTBeam()
{
    using Beam = xraymc::CBCTBeam<>;
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
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

// CudaTransport example
//
// Transports 80 keV photons through a 30x30x30 cm water phantom using
// Woodcock delta-tracking on the GPU and prints the central-axis depth dose.
//
// Build requirements:
//   cmake -DXRAYMCLIB_ENABLE_CUDA=ON ...
//   Requires a CUDA-capable GPU (SM >= 6.0 recommended for native double atomicAdd).

#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp" // pulls in cuda_transport.hpp
#include "xraymc/xraymcrandom.hpp"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

void example()
{
    std::cout << "CudaTransport example\n";
    std::cout << "Woodcock transport of 80 keV photons in a 30x30x30 cm water phantom\n\n";

    // ── 1. Build a 50x50x50 water phantom, 0.6 cm voxels, centered at origin ──
    // AABB: [-15, -15, -15] to [15, 15, 15] cm
    constexpr std::size_t NX = 50, NY = 50, NZ = 50;
    constexpr double SPACING = 0.6; // cm

    // CudaTransport requires the Woodcock mode: AAVoxelGrid<16, 3, 255>
    using Grid = xraymc::CudaTransport::Grid;
    using Mat = xraymc::Material<16>;

    const std::size_t nVoxels = NX * NY * NZ;

    auto water = Mat::byChemicalFormula("H2O").value();
    std::vector<double> densities(nVoxels, 1.0); // water: 1.0 g/cm³
    std::vector<std::uint8_t> matIdx(nVoxels, 0); // all water (index 0)
    std::vector<Mat> materials = { water };

    Grid grid({ NX, NY, NZ }, { SPACING, SPACING, SPACING },
        densities, matIdx, materials);

    std::cout << "Phantom: " << NX << "x" << NY << "x" << NZ
              << " voxels, " << SPACING << " cm spacing\n"
              << "AABB: ["
              << grid.AABB()[0] << ", " << grid.AABB()[1] << ", " << grid.AABB()[2]
              << "] to ["
              << grid.AABB()[3] << ", " << grid.AABB()[4] << ", " << grid.AABB()[5]
              << "] cm\n\n";

    // ── 2. Upload geometry and cross-section tables to the GPU ────────────────
    xraymc::CudaTransport cuda(grid);
    cuda.uploadGrid();
    std::cout << "Grid uploaded to GPU\n";

    // ── 3. Sample a particle batch ────────────────────────────────────────────
    // Source: 60 cm upstream on -z axis, beam along +z.
    // Collimate to a 10x10 cm field at the phantom front face (SSD = 45 cm):
    //   half-angle = atan(5 / 45) ~ 6.3 degrees
    constexpr double ENERGY = 80.0; // keV
    constexpr double SSD = 45.0; // source-to-surface distance [cm]
    constexpr double FIELD_SIZE = 10.0; // full field size at surface [cm]
    const double hAngle = 2 * xraymc::DEG_TO_RAD();

    // Default beam direction is +z; default cosines are {1,0,0} and {0,1,0}.
    xraymc::IsotropicMonoEnergyBeam<> beam;
    beam.setPosition({ 0.0, 0.0, -60.0 });
    beam.setEnergy(ENERGY);

    beam.setCollimationHalfAngles(hAngle, hAngle);

    constexpr std::uint64_t N_PARTICLES = 1'000'000;
    beam.setNumberOfParticlesPerExposure(N_PARTICLES);

    std::vector<xraymc::CudaTransport::Particle> batch;
    batch.reserve(N_PARTICLES);

    xraymc::RandomState rng;
    const auto exp = beam.exposure(0);
    for (std::uint64_t k = 0; k < N_PARTICLES; ++k) {
        auto p = exp.sampleParticle(rng);
        batch.push_back(xraymc::CudaTransport::toCudaParticle(p));
    }
    std::cout << "Sampled " << N_PARTICLES << " photons at " << ENERGY << " keV, "
              << "field " << FIELD_SIZE << "x" << FIELD_SIZE << " cm at SSD=" << SSD << " cm\n";

    // ── 4. Run Woodcock transport on the GPU ──────────────────────────────────
    cuda.runBatch(batch);

    // Retrieve per-voxel energy scored and merge into the grid's EnergyScore
    // accumulators. Device buffer is zeroed so a subsequent runBatch starts fresh.
    cuda.downloadScores();
    std::cout << "Transport complete\n\n";

    // ── 5. Convert scored energy to absorbed dose ─────────────────────────────
    // calibration_factor = 1 gives dose in keV/g per history.
    // Divide by N_PARTICLES to get dose per photon history.
    const double calibration = 1.0 / static_cast<double>(N_PARTICLES);
    grid.addEnergyScoredToDoseScore(calibration);

    // ── 6. Print central-axis depth dose ─────────────────────────────────────
    // Central column: ix = NX/2, iy = NY/2, stepping through all NZ slices.
    const std::size_t cx = NX / 2;
    const std::size_t cy = NY / 2;

    std::cout << "Central-axis depth dose (voxel column ix=" << cx << ", iy=" << cy << "):\n";
    std::cout << std::left
              << std::setw(12) << "Depth[cm]"
              << std::setw(22) << "Dose[keV/g/photon]"
              << std::setw(16) << "Rel.unc.[%]"
              << "Events\n";

    for (std::size_t iz = 0; iz < NZ; ++iz) {
        const auto idx = grid.flatIndex(cx, cy, iz);
        const auto& ds = grid.doseScored(idx);
        const double depth = (iz + 0.5) * SPACING; // depth from phantom front face [cm]

        std::cout << std::left
                  << std::setw(12) << depth
                  << std::setw(22) << ds.dose()
                  << std::setw(16) << ds.relativeUncertainty() * 100.0
                  << ds.numberOfEvents() << '\n';
    }
}

int main()
{
    example();
    return 0;
}

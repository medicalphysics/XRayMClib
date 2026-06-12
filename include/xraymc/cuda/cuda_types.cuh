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

#pragma once

#include <cstddef>
#include <cstdint>

namespace xraymc::cuda_detail {

// Number of log-spaced energy points stored per material cross-section table.
// 128 points spans 1–151 keV with sufficient log-spacing for linear interpolation
// to match the CPU cubic-spline accuracy to better than 0.1 % across the range.
inline constexpr int CUDA_TABLE_N = 128;

// Number of log-spaced momentum-transfer points for form-factor F²(q) and
// incoherent scatter function S(q) tables.
inline constexpr int CUDA_FORM_N = 128;

// Fixed shell count — must match the NMaterialShells template parameter used
// when constructing the AAVoxelGrid that owns the materials.
inline constexpr int CUDA_N_SHELLS = 16;

// ─── Particle ─────────────────────────────────────────────────────────────────
// Plain-old-data mirror of xraymc::Particle for device transfer.
// Double precision throughout: position and direction errors accumulate across
// many transport steps and float precision is insufficient for large grids.
struct CudaParticle {
    double posX, posY, posZ; // world-space position [cm]
    double dirX, dirY, dirZ; // unit direction vector
    double energy;           // photon energy [keV]
    double weight;           // statistical weight
};

// ─── Per-thread PCG32 state ───────────────────────────────────────────────────
struct CudaRng {
    uint64_t state;
    uint64_t inc;
};

// ─── Voxel grid geometry ──────────────────────────────────────────────────────
// Passed by value to the kernel to avoid an extra pointer dereference per step.
struct CudaGridGeom {
    uint32_t nx, ny, nz;                            // voxel counts
    double   aabbXmin, aabbYmin, aabbZmin;          // AABB lower corner [cm]
    double   aabbXmax, aabbYmax, aabbZmax;          // AABB upper corner [cm]
    double   invSpacingX, invSpacingY, invSpacingZ; // 1/spacing [1/cm]
};

// ─── Material descriptor for LEC=3 physics ───────────────────────────────────
//
// LEC=3 uses:
//   Compton  — Klein-Nishina + Livermore incoherent scatter factor S(q)
//   Rayleigh — form-factor F²(q) with dipole rejection
//   Photo-el — full absorption (no fluorescence)
//
// All tables stored as float (single precision) to halve memory bandwidth;
// the interpolated result is cast back to double before use in kinematics.
//
// Energy table axis: log(E/keV), log-spaced over [log(1), log(151)].
// Momentum-transfer axis: log(q/Å⁻¹), log-spaced over [log(1e-3), log(qmax)].
struct CudaMaterial {
    // Total attenuation tables
    float logEnergy[CUDA_TABLE_N];        // log(E/keV) — shared x-axis
    float logPhotoelectric[CUDA_TABLE_N]; // log(μ_pe/ρ)  [cm²/g]
    float logIncoherent[CUDA_TABLE_N];    // log(μ_ic/ρ)  [cm²/g]
    float logCoherent[CUDA_TABLE_N];      // log(μ_coh/ρ) [cm²/g]
    int   nPoints;                        // active entries (≤ CUDA_TABLE_N)

    // Momentum-transfer tables for Rayleigh and Compton corrections
    float logMomTransfer[CUDA_FORM_N];    // log(q/Å⁻¹) — shared x-axis
    float logFormFactorSq[CUDA_FORM_N];  // log(F²(q)) for Rayleigh
    float logScatterFunc[CUDA_FORM_N];   // log(S(q))  for Compton Livermore
    int   nFormPoints;                    // active entries (≤ CUDA_FORM_N)

    // Scalar material properties
    float formFactorSqAt0; // F²(q→0) ≈ Z_eff² — maximum used in FF rejection
    float effectiveZ;      // Z_eff — normalisation for Livermore scatter factor
};

// ─── Woodcock majorant table ──────────────────────────────────────────────────
// max(μρ) over all materials and densities in the grid, as a function of energy.
struct CudaWoodcockTable {
    float logEnergy[CUDA_TABLE_N];  // log(E/keV)
    float logMaxAtt[CUDA_TABLE_N];  // log(max(μρ)) [cm⁻¹]
    int   nPoints;
};

// ─── Device buffer bookkeeping ────────────────────────────────────────────────
// Owned by CudaTransport; holds raw device pointers.
// All allocation and deallocation is done in aavoxelgrid_cuda.cu.
struct CudaDeviceBuffers {
    uint8_t*           d_matIndex     = nullptr; // per-voxel material index
    float*             d_density      = nullptr; // per-voxel density [g/cm³]
    double*            d_energyScored = nullptr; // per-voxel accumulated energy [keV]
    CudaMaterial*      d_materials    = nullptr; // material table array
    CudaWoodcockTable* d_woodcock     = nullptr; // majorant table
    CudaParticle*      d_particles    = nullptr; // particle batch (grown on demand)
    std::size_t        nVoxels        = 0;
    std::size_t        nParticleAlloc = 0;
};

// ─── Host API (implemented in src/cuda/aavoxelgrid_cuda.cu) ──────────────────

// Upload a complete voxel grid to device memory.
// Allocates all device buffers; call cudaFreeBuffers() before calling again.
void cudaUploadGrid(CudaDeviceBuffers&        bufs,
                    const uint8_t*            hMatIndex,
                    const float*              hDensity,
                    std::size_t               nVoxels,
                    const CudaMaterial*       hMaterials,
                    std::size_t               nMaterials,
                    const CudaWoodcockTable*  hWoodcock);

// Run a batch of particle histories on the device via the Woodcock kernel.
// Energy scored per voxel is accumulated atomically in d_energyScored.
void cudaLaunchWoodcock(CudaDeviceBuffers&       bufs,
                         const CudaParticle*     hParticles,
                         uint64_t                nParticles,
                         const CudaGridGeom&     geom,
                         uint64_t                rngSeedOffset);

// Copy per-voxel energy scores from device to hEnergy and zero the device buffer.
void cudaDownloadEnergyScored(CudaDeviceBuffers& bufs, double* hEnergy);

// Release all device allocations in bufs and zero the struct.
void cudaFreeBuffers(CudaDeviceBuffers& bufs);

} // namespace xraymc::cuda_detail

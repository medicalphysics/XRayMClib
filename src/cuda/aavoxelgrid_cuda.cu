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

// Compiled by nvcc as a CUDA translation unit.
// Implements the Woodcock delta-tracking kernel for AAVoxelGrid<16, 3, 255>
// and the host-side upload / launch / download functions declared in cuda_types.cuh.

#include "xraymc/cuda/cuda_types.cuh"

#include <cuda_runtime.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace xraymc::cuda_detail {

// ─── CUDA error checking ──────────────────────────────────────────────────────
#define CUDA_CHECK(expr)                                                        \
    do {                                                                        \
        cudaError_t _err = (expr);                                              \
        if (_err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d — %s\n",                      \
                    __FILE__, __LINE__, cudaGetErrorString(_err));              \
            ::std::abort();                                                     \
        }                                                                       \
    } while (0)

// ─── Physical and numerical constants ────────────────────────────────────────
// These mirror xraymc::constants.hpp values without pulling in host-only headers.
static constexpr double GEOM_EPS     = 1.0e-6;          // cm — GEOMETRIC_ERROR()
static constexpr double MIN_ENERGY_D = 1.0;              // keV — XRAYMCLIB_MINENERGY
static constexpr double M_E_KEV      = 510.9989461;      // electron rest mass [keV]
static constexpr double HC_KEV_AA    = 12.398520;        // hc [keV·Å]
static constexpr double TWO_PI       = 6.28318530717958647692;

// ─── PCG32 per-thread RNG ─────────────────────────────────────────────────────
// "Really Minimal PCG32" — M.E. O'Neill, 2014 (Apache-2.0).
// One CudaRng instance is held in registers for the lifetime of the thread.

__device__ __forceinline__ uint32_t pcg32Next(CudaRng& rng) noexcept
{
    uint64_t old = rng.state;
    rng.state    = old * 6364136223846793005ULL + rng.inc;
    uint32_t xsh = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = static_cast<uint32_t>(old >> 59u);
    return (xsh >> rot) | (xsh << ((-(int)rot) & 31));
}

// Returns a double in (0, 1] — the +0.5 offset avoids exact 0, which would
// produce -inf when passed to log().
__device__ __forceinline__ double pcg32Uniform(CudaRng& rng) noexcept
{
    return (static_cast<double>(pcg32Next(rng)) + 0.5) * (1.0 / 4294967296.0);
}

__device__ void pcg32Init(CudaRng& rng, uint64_t seed, uint64_t seq) noexcept
{
    rng.state = 0u;
    rng.inc   = (seq << 1u) | 1u;
    pcg32Next(rng);
    rng.state += seed;
    pcg32Next(rng);
}

// ─── Log-space binary search + linear interpolation ──────────────────────────
// logX must be sorted ascending. Returns the linearly interpolated logY value
// at logXq, clamped to the table endpoints.
__device__ __forceinline__ float logLerp(const float* __restrict__ logX,
                                          const float* __restrict__ logY,
                                          int n, float logXq) noexcept
{
    if (logXq <= logX[0])   return logY[0];
    if (logXq >= logX[n-1]) return logY[n-1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (logX[mid] <= logXq) lo = mid; else hi = mid;
    }
    float t = (logXq - logX[lo]) / (logX[hi] - logX[lo]);
    return logY[lo] + t * (logY[hi] - logY[lo]);
}

// ─── Attenuation lookup ───────────────────────────────────────────────────────
struct DevAtt { double photoelectric, incoherent, coherent; };

__device__ __forceinline__ DevAtt attenuationAt(const CudaMaterial& mat,
                                                  double energy) noexcept
{
    float lE = __double2float_rn(log(energy));
    DevAtt a;
    a.photoelectric = exp(static_cast<double>(
        logLerp(mat.logEnergy, mat.logPhotoelectric, mat.nPoints, lE)));
    a.incoherent    = exp(static_cast<double>(
        logLerp(mat.logEnergy, mat.logIncoherent,    mat.nPoints, lE)));
    a.coherent      = exp(static_cast<double>(
        logLerp(mat.logEnergy, mat.logCoherent,      mat.nPoints, lE)));
    return a;
}

__device__ __forceinline__ double attSum(const DevAtt& a) noexcept
{
    return a.photoelectric + a.incoherent + a.coherent;
}

// ─── Woodcock majorant lookup ─────────────────────────────────────────────────
__device__ __forceinline__ double woodcockMaxAtt(const CudaWoodcockTable& wt,
                                                   double energy) noexcept
{
    float lE = __double2float_rn(log(energy));
    return exp(static_cast<double>(
        logLerp(wt.logEnergy, wt.logMaxAtt, wt.nPoints, lE)));
}

// ─── Momentum transfer helpers ────────────────────────────────────────────────
__device__ __forceinline__ double momentumTransferMax(double energy) noexcept
{
    return energy / HC_KEV_AA; // Å⁻¹
}

// q for a given incident energy and cosine of scattering angle.
// q = E/(hc) · sqrt((1 − cosθ)/2)
__device__ __forceinline__ double momentumTransferFromCos(double energy,
                                                           double cosTheta) noexcept
{
    return momentumTransferMax(energy) * sqrt(fmax(0.0, (1.0 - cosTheta) * 0.5));
}

// ─── Form factor and scatter function lookups ─────────────────────────────────
__device__ __forceinline__ double scatterFuncAt(const CudaMaterial& mat,
                                                 double q) noexcept
{
    if (q <= 0.0) return 0.0;
    float lq = __double2float_rn(log(q));
    return exp(static_cast<double>(
        logLerp(mat.logMomTransfer, mat.logScatterFunc, mat.nFormPoints, lq)));
}

__device__ __forceinline__ double formFactorSqAt(const CudaMaterial& mat,
                                                  double q) noexcept
{
    if (q <= 0.0) return static_cast<double>(mat.formFactorSqAt0);
    float lq = __double2float_rn(log(q));
    return exp(static_cast<double>(
        logLerp(mat.logMomTransfer, mat.logFormFactorSq, mat.nFormPoints, lq)));
}

// ─── 3D vector helpers ────────────────────────────────────────────────────────
struct Vec3 { double x, y, z; };

__device__ __forceinline__ Vec3 cross3(Vec3 a, Vec3 b) noexcept
{
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}

__device__ __forceinline__ double dot3(Vec3 a, Vec3 b) noexcept
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ __forceinline__ void normalize3(Vec3& v) noexcept
{
    double inv = rsqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x *= inv; v.y *= inv; v.z *= inv;
}

// Rodrigues' rotation: rotate v around unit axis k by angle (sin, cos).
// result = (v·k)k + cos·(k×v)×k + sin·(k×v)
__device__ __forceinline__ Vec3 rodrigues(Vec3 v, Vec3 k,
                                           double sinA, double cosA) noexcept
{
    Vec3   kxv = cross3(k, v);
    double vk  = dot3(v, k);
    Vec3   v2  = cross3(kxv, k);
    return { vk * k.x + cosA * v2.x + sinA * kxv.x,
             vk * k.y + cosA * v2.y + sinA * kxv.y,
             vk * k.z + cosA * v2.z + sinA * kxv.z };
}

// Deflects unit direction dir by polar angle theta and azimuthal phi.
// Mirrors vectormath::peturb: constructs a frame-independent orthogonal axis,
// rotates it by phi around dir, then rotates dir by theta around that axis.
__device__ Vec3 peturb3(Vec3 dir, double theta, double phi) noexcept
{
    // Pick world axis with smallest |component| to avoid numerical cancellation.
    double ax = fabs(dir.x), ay = fabs(dir.y), az = fabs(dir.z);
    Vec3 k = { 0.0, 0.0, 0.0 };
    if (ax <= ay && ax <= az)      k.x = 1.0;
    else if (ay <= ax && ay <= az) k.y = 1.0;
    else                           k.z = 1.0;

    Vec3 perp = cross3(dir, k);
    normalize3(perp);

    // Rotate perp by phi around dir to obtain the scatter-plane normal.
    double sp, cp;
    sincos(phi, &sp, &cp);
    Vec3 axis = rodrigues(perp, dir, sp, cp);
    normalize3(axis);

    // Rotate dir by theta around the scatter-plane normal.
    double st, ct;
    sincos(theta, &st, &ct);
    Vec3 res = rodrigues(dir, axis, st, ct);
    normalize3(res);
    return res;
}

// ─── Geometry helpers ─────────────────────────────────────────────────────────

__device__ __forceinline__ bool insideAABB(double x, double y, double z,
                                            const CudaGridGeom& g) noexcept
{
    return x > g.aabbXmin && x < g.aabbXmax
        && y > g.aabbYmin && y < g.aabbYmax
        && z > g.aabbZmin && z < g.aabbZmax;
}

// Returns the distance along the ray (p.pos, p.dir) to the exit face of the AABB.
// Assumes the particle is already inside; returns 0 if all slabs are parallel.
__device__ double distToAABBExit(const CudaParticle& p,
                                  const CudaGridGeom& g) noexcept
{
    double tmax = 1.0e30;
    double t1, t2, tmp;

    if (fabs(p.dirX) > 1.0e-14) {
        t1 = (g.aabbXmin - p.posX) / p.dirX;
        t2 = (g.aabbXmax - p.posX) / p.dirX;
        if (t1 > t2) { tmp = t1; t1 = t2; t2 = tmp; }
        tmax = fmin(tmax, t2);
    }
    if (fabs(p.dirY) > 1.0e-14) {
        t1 = (g.aabbYmin - p.posY) / p.dirY;
        t2 = (g.aabbYmax - p.posY) / p.dirY;
        if (t1 > t2) { tmp = t1; t1 = t2; t2 = tmp; }
        tmax = fmin(tmax, t2);
    }
    if (fabs(p.dirZ) > 1.0e-14) {
        t1 = (g.aabbZmin - p.posZ) / p.dirZ;
        t2 = (g.aabbZmax - p.posZ) / p.dirZ;
        if (t1 > t2) { tmp = t1; t1 = t2; t2 = tmp; }
        tmax = fmin(tmax, t2);
    }
    return tmax;
}

// Convert a world-space position to a flat row-major voxel index (x-fastest).
// Clamps to valid range rather than checking bounds — the caller ensures we
// are inside the AABB.
__device__ __forceinline__ uint32_t worldToFlat(double x, double y, double z,
                                                  const CudaGridGeom& g) noexcept
{
    uint32_t ix = static_cast<uint32_t>((x - g.aabbXmin) * g.invSpacingX);
    uint32_t iy = static_cast<uint32_t>((y - g.aabbYmin) * g.invSpacingY);
    uint32_t iz = static_cast<uint32_t>((z - g.aabbZmin) * g.invSpacingZ);
    ix = min(ix, g.nx - 1u);
    iy = min(iy, g.ny - 1u);
    iz = min(iz, g.nz - 1u);
    return ix + g.nx * (iy + g.ny * iz);
}

// ─── Interaction physics — LEC = 3 ────────────────────────────────────────────
//
// LEC=3 maps to the following branch in interactions.hpp:
//   comptonScatter   → else (non-IA): Klein-Nishina + Livermore scatter factor
//   rayleightScatter → else (non-zero): form-factor corrected
//   photoelectricEffect → else (non-IA): full absorption, no fluorescence

struct DevInteract {
    double energyImparted; // keV · weight
    bool   alive;          // particle survives the interaction
    bool   energyChanged;  // true when energy changes (Compton); governs majorant update
};

// ── Compton: Klein-Nishina + Livermore incoherent scatter factor ──────────────
//
// Mirrors interactions::comptonScatter<N, 3>. Samples the energy ratio
// e = E'/E from the KN distribution with an additional Livermore S(q) weight.
__device__ DevInteract comptonLivermore(CudaParticle& p,
                                         const CudaMaterial& mat,
                                         CudaRng& rng) noexcept
{
    const double E       = p.energy;
    const double k       = E / M_E_KEV;         // E / m_e c²
    const double emin    = 1.0 / (1.0 + 2.0 * k);
    const double gmaxInv = emin / (1.0 + emin * emin); // 1 / g_max
    const double Zeff    = static_cast<double>(mat.effectiveZ);

    double e, cosTheta;
    for (;;) {
        const double r1   = pcg32Uniform(rng);
        e                 = r1 + (1.0 - r1) * emin; // uniform in [emin, 1]
        double t          = fmin((1.0 - e) / (k * e), 2.0); // guard rounding
        cosTheta          = 1.0 - t;
        double sinThetaSq = 1.0 - cosTheta * cosTheta;
        double g          = (1.0 / e + e - sinThetaSq) * gmaxInv; // ∈ (0, 1]

        const double q = momentumTransferFromCos(E, cosTheta);
        const double S = scatterFuncAt(mat, q); // S(q) ∈ [0, Z_eff]

        // Accept with probability g · S(q) / Z_eff
        if (pcg32Uniform(rng) * Zeff <= g * S) break;
    }

    const double Esc   = E * e;
    const double theta = acos(cosTheta);
    const double phi   = pcg32Uniform(rng) * TWO_PI;

    Vec3 dir = { p.dirX, p.dirY, p.dirZ };
    dir = peturb3(dir, theta, phi);
    p.dirX = dir.x; p.dirY = dir.y; p.dirZ = dir.z;
    p.energy = Esc;

    DevInteract r;
    r.energyImparted = (E - Esc) * p.weight;
    r.alive          = p.energy >= MIN_ENERGY_D;
    r.energyChanged  = true;
    return r;
}

// ── Rayleigh: form-factor corrected with dipole rejection ─────────────────────
//
// Mirrors interactions::rayleightScatter<N, 3>.
// Samples q² ∈ [0, qmax²] proportional to F²(q) via rejection on F²(q)/F²(0),
// then applies the dipole-distribution rejection (1 + cos²θ)/2.
__device__ DevInteract rayleighFF(CudaParticle& p,
                                   const CudaMaterial& mat,
                                   CudaRng& rng) noexcept
{
    const double qmax    = momentumTransferMax(p.energy);
    const double qmax2   = qmax * qmax;
    const double FF2max  = static_cast<double>(mat.formFactorSqAt0);

    double cosTheta;
    for (;;) {
        // Sample q² uniformly in [0, qmax²]; compute q and look up F²(q).
        const double q2  = pcg32Uniform(rng) * qmax2;
        const double q   = sqrt(q2);
        const double ff2 = formFactorSqAt(mat, q);

        // Form-factor rejection: accept with probability F²(q) / F²(0).
        if (pcg32Uniform(rng) * FF2max > ff2) continue;

        cosTheta = 1.0 - 2.0 * q2 / qmax2;

        // Dipole rejection: (1 + cos²θ)/2
        if (pcg32Uniform(rng) <= 0.5 * (1.0 + cosTheta * cosTheta)) break;
    }

    const double theta = acos(cosTheta);
    const double phi   = pcg32Uniform(rng) * TWO_PI;

    Vec3 dir = { p.dirX, p.dirY, p.dirZ };
    dir = peturb3(dir, theta, phi);
    p.dirX = dir.x; p.dirY = dir.y; p.dirZ = dir.z;

    DevInteract r;
    r.energyImparted = 0.0;
    r.alive          = true;
    r.energyChanged  = false;
    return r;
}

// ── Photoelectric: full absorption (LEC=3 — no fluorescence) ─────────────────
//
// Mirrors interactions::photoelectricEffect<N, 3> (the else-branch).
__device__ __forceinline__ DevInteract photoelectricAbsorb(CudaParticle& p) noexcept
{
    DevInteract r;
    r.energyImparted = p.energy * p.weight;
    r.alive          = false;
    r.energyChanged  = false;
    p.energy         = 0.0;
    return r;
}

// ── Interaction dispatcher ────────────────────────────────────────────────────
__device__ DevInteract deviceInteract(const DevAtt& att,
                                       CudaParticle& p,
                                       const CudaMaterial& mat,
                                       CudaRng& rng) noexcept
{
    const double r = pcg32Uniform(rng) * attSum(att);
    if (r < att.photoelectric)
        return photoelectricAbsorb(p);
    if (r < att.photoelectric + att.incoherent)
        return comptonLivermore(p, mat, rng);
    return rayleighFF(p, mat, rng);
}

// ─── Woodcock (delta-tracking) transport kernel ───────────────────────────────
//
// One CUDA thread per particle history. Each thread loops until the particle
// exits the grid AABB or is absorbed. The Woodcock majorant is refreshed after
// every Compton event (energy changes); Rayleigh and the initial step use the
// cached value.
//
// atomicAdd(double*) requires SM ≥ 6.0 for hardware support. Older devices
// fall back to a CAS-loop emulation automatically by the CUDA runtime.
__global__ void woodcockKernel(
    const CudaParticle*        __restrict__ inParticles,
    uint64_t                               nParticles,
    const uint8_t*             __restrict__ matIndex,
    const float*               __restrict__ density,
    double*                                energyScored,
    const CudaMaterial*        __restrict__ materials,
    const CudaWoodcockTable*   __restrict__ woodcock,
    CudaGridGeom                           geom,
    uint64_t                               rngOffset)
{
    const uint64_t tid = static_cast<uint64_t>(blockIdx.x) * blockDim.x
                         + static_cast<uint64_t>(threadIdx.x);
    if (tid >= nParticles) return;

    CudaParticle p = inParticles[tid];

    CudaRng rng;
    pcg32Init(rng, tid + rngOffset, /*seq=*/1u);

    bool   valid     = insideAABB(p.posX, p.posY, p.posZ, geom);
    bool   updateAtt = true;
    double attMaxInv = 0.0;

    while (valid) {
        if (updateAtt) {
            const double mx = woodcockMaxAtt(*woodcock, p.energy);
            attMaxInv  = mx > 0.0 ? 1.0 / mx : 0.0;
            updateAtt  = false;
        }

        const double step     = -log(pcg32Uniform(rng)) * attMaxInv;
        const double exitDist = distToAABBExit(p, geom);

        if (step < exitDist) {
            // Move to candidate interaction site.
            p.posX += p.dirX * step;
            p.posY += p.dirY * step;
            p.posZ += p.dirZ * step;

            const uint32_t fi     = worldToFlat(p.posX, p.posY, p.posZ, geom);
            const uint8_t  matIdx = matIndex[fi];

            // Material index 255 is TRANSPARENTVOXELS — treat as outside.
            if (matIdx == 255u) { valid = false; break; }

            const double  dens   = static_cast<double>(density[fi]);
            const DevAtt  att    = attenuationAt(materials[matIdx], p.energy);
            const double  attTot = attSum(att) * dens; // linear attenuation [1/cm]

            // Woodcock acceptance: real interaction with probability attTot·attMaxInv.
            if (pcg32Uniform(rng) < attTot * attMaxInv) {
                const DevInteract res = deviceInteract(att, p, materials[matIdx], rng);
                if (res.energyImparted > 0.0)
                    atomicAdd(energyScored + fi, res.energyImparted);
                valid     = res.alive;
                updateAtt = res.energyChanged;
            }
            // else: virtual interaction — no state change, continue loop.

        } else {
            // Particle exits the grid — nudge past the boundary.
            p.posX += p.dirX * (exitDist + GEOM_EPS);
            p.posY += p.dirY * (exitDist + GEOM_EPS);
            p.posZ += p.dirZ * (exitDist + GEOM_EPS);
            valid = false;
        }
    }
}

// ─── Host-side API implementations ───────────────────────────────────────────

void cudaUploadGrid(CudaDeviceBuffers&        bufs,
                    const uint8_t*            hMatIndex,
                    const float*              hDensity,
                    std::size_t               nVoxels,
                    const CudaMaterial*       hMaterials,
                    std::size_t               nMaterials,
                    const CudaWoodcockTable*  hWoodcock)
{
    bufs.nVoxels = nVoxels;

    CUDA_CHECK(cudaMalloc(&bufs.d_matIndex,     nVoxels                          ));
    CUDA_CHECK(cudaMalloc(&bufs.d_density,      nVoxels * sizeof(float)          ));
    CUDA_CHECK(cudaMalloc(&bufs.d_energyScored, nVoxels * sizeof(double)         ));
    CUDA_CHECK(cudaMalloc(&bufs.d_materials,    nMaterials * sizeof(CudaMaterial)));
    CUDA_CHECK(cudaMalloc(&bufs.d_woodcock,     sizeof(CudaWoodcockTable)        ));

    CUDA_CHECK(cudaMemcpy(bufs.d_matIndex,  hMatIndex,  nVoxels,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(bufs.d_density,   hDensity,   nVoxels * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(bufs.d_energyScored, 0, nVoxels * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(bufs.d_materials, hMaterials, nMaterials * sizeof(CudaMaterial),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(bufs.d_woodcock,  hWoodcock,  sizeof(CudaWoodcockTable),
                          cudaMemcpyHostToDevice));
}

void cudaLaunchWoodcock(CudaDeviceBuffers&    bufs,
                         const CudaParticle*  hParticles,
                         uint64_t             nParticles,
                         const CudaGridGeom&  geom,
                         uint64_t             rngOffset)
{
    // Grow the particle buffer only when needed; never shrink.
    if (nParticles > bufs.nParticleAlloc) {
        if (bufs.d_particles) CUDA_CHECK(cudaFree(bufs.d_particles));
        CUDA_CHECK(cudaMalloc(&bufs.d_particles, nParticles * sizeof(CudaParticle)));
        bufs.nParticleAlloc = nParticles;
    }
    CUDA_CHECK(cudaMemcpy(bufs.d_particles, hParticles,
                          nParticles * sizeof(CudaParticle), cudaMemcpyHostToDevice));

    constexpr int THREADS_PER_BLOCK = 128;
    const int blocks = static_cast<int>(
        (nParticles + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);

    woodcockKernel<<<blocks, THREADS_PER_BLOCK>>>(
        bufs.d_particles, nParticles,
        bufs.d_matIndex, bufs.d_density, bufs.d_energyScored,
        bufs.d_materials, bufs.d_woodcock,
        geom, rngOffset);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void cudaDownloadEnergyScored(CudaDeviceBuffers& bufs, double* hEnergy)
{
    CUDA_CHECK(cudaMemcpy(hEnergy, bufs.d_energyScored,
                          bufs.nVoxels * sizeof(double), cudaMemcpyDeviceToHost));
    // Clear device accumulator so the next batch starts from zero.
    CUDA_CHECK(cudaMemset(bufs.d_energyScored, 0, bufs.nVoxels * sizeof(double)));
}

void cudaFreeBuffers(CudaDeviceBuffers& bufs)
{
    if (bufs.d_matIndex)     CUDA_CHECK(cudaFree(bufs.d_matIndex));
    if (bufs.d_density)      CUDA_CHECK(cudaFree(bufs.d_density));
    if (bufs.d_energyScored) CUDA_CHECK(cudaFree(bufs.d_energyScored));
    if (bufs.d_materials)    CUDA_CHECK(cudaFree(bufs.d_materials));
    if (bufs.d_woodcock)     CUDA_CHECK(cudaFree(bufs.d_woodcock));
    if (bufs.d_particles)    CUDA_CHECK(cudaFree(bufs.d_particles));
    bufs = {};
}

} // namespace xraymc::cuda_detail

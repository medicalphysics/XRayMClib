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

// This header is included at the bottom of aavoxelgrid.hpp only when
// XRAYMCLIB_CUDA_ENABLED is defined. AAVoxelGrid<16,3,255> is therefore
// already visible here.

#pragma once

#include "xraymc/cuda/cuda_types.cuh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace xraymc {

// ─── CudaTransport ────────────────────────────────────────────────────────────
/**
 * @brief CUDA-accelerated Woodcock transport driver for AAVoxelGrid<16, 3, 255>.
 *
 * Uploads the voxel grid geometry and cross-section tables to the GPU once via
 * uploadGrid(), then executes batches of particle histories entirely on device
 * using the delta-tracking (Woodcock) algorithm.  Per-voxel energy scored is
 * accumulated with `atomicAdd` and merged back into the host grid's EnergyScore
 * accumulators by downloadScores().
 *
 * ### Physics (LOWENERGYCORRECTION = 3)
 * | Process     | Model                                        |
 * |-------------|----------------------------------------------|
 * | Compton     | Klein–Nishina + Livermore scatter factor S(q)|
 * | Rayleigh    | Form-factor corrected, dipole rejection      |
 * | Photoelectric | Full absorption (no fluorescence)          |
 *
 * ### Typical usage
 * @code
 *   AAVoxelGrid<16, 3> grid = ...;           // or <16, 3, 255> explicitly
 *   CudaTransport cuda(grid);
 *   cuda.uploadGrid();                        // once per geometry change
 *
 *   std::vector<CudaTransport::Particle> batch;
 *   for (auto& p : sampledParticles)
 *       batch.push_back(CudaTransport::toCudaParticle(p));
 *   cuda.runBatch(batch);
 *   cuda.downloadScores();                    // merges into grid EnergyScore
 *   grid.addEnergyScoredToDoseScore(cal);
 * @endcode
 *
 * ### Notes
 * - `atomicAdd(double*)` requires SM ≥ 6.0 (Pascal) for hardware support;
 *   older GPUs use a CAS-loop emulation which is slower but still correct.
 * - Only the Woodcock mode (`TRANSPARENTVOXELS = 255`) is supported; there is
 *   no CUDA path for voxel-by-voxel analog transport.
 * - The class owns device memory and is non-copyable.
 */
class CudaTransport {
public:
    using Grid    = AAVoxelGrid<16, 2, 255>;
    using Particle = cuda_detail::CudaParticle;

    // ── Construction / destruction ────────────────────────────────────────────

    explicit CudaTransport(Grid& grid) noexcept
        : m_grid(grid)
    {
    }

    ~CudaTransport()
    {
        if (m_uploaded)
            cuda_detail::cudaFreeBuffers(m_bufs);
    }

    CudaTransport(const CudaTransport&)            = delete;
    CudaTransport& operator=(const CudaTransport&) = delete;
    CudaTransport(CudaTransport&&)                 = delete;
    CudaTransport& operator=(CudaTransport&&)      = delete;

    // ── Public interface ──────────────────────────────────────────────────────

    /**
     * @brief Converts the host AAVoxelGrid to flat GPU-friendly arrays and
     *        uploads everything to device memory.
     *
     * Must be called once before the first runBatch() and again whenever the
     * grid geometry or material assignments change.
     *
     * @throws std::runtime_error if the grid has no materials or voxels.
     */
    void uploadGrid()
    {
        if (m_uploaded)
            cuda_detail::cudaFreeBuffers(m_bufs);

        const auto  nVoxels   = m_grid.size();
        const auto& hostMats  = m_grid.m_materials;

        if (nVoxels == 0 || hostMats.empty())
            throw std::runtime_error("CudaTransport::uploadGrid — grid is empty");

        // ── Flatten materials ─────────────────────────────────────────────────
        std::vector<cuda_detail::CudaMaterial> cudaMats(hostMats.size());
        for (std::size_t m = 0; m < hostMats.size(); ++m)
            flattenMaterial(hostMats[m], cudaMats[m]);

        // ── Flatten Woodcock majorant table ───────────────────────────────────
        cuda_detail::CudaWoodcockTable cudaWt;
        flattenWoodcockTable(m_grid.m_woodcockStepTableLin, cudaWt);

        // ── Pack voxel SoA arrays (better coalescing than AoS on GPU) ─────────
        std::vector<uint8_t> matIdx(nVoxels);
        std::vector<float>   density(nVoxels);
        for (std::size_t i = 0; i < nVoxels; ++i) {
            matIdx[i]  = m_grid.m_data[i].materialIndex;
            density[i] = static_cast<float>(m_grid.m_data[i].density);
        }

        cuda_detail::cudaUploadGrid(
            m_bufs,
            matIdx.data(), density.data(), nVoxels,
            cudaMats.data(), cudaMats.size(),
            &cudaWt);

        m_geom     = buildGeom();
        m_uploaded = true;
    }

    /**
     * @brief Transfers a pre-sampled particle batch to the device and runs the
     *        Woodcock kernel synchronously.
     *
     * Energy deposited per voxel accumulates on the device across multiple
     * runBatch() calls; call downloadScores() to retrieve the totals.
     *
     * @param particles Host-side particle batch; transferred to device memory.
     * @throws std::logic_error if uploadGrid() has not been called first.
     */
    void runBatch(const std::vector<Particle>& particles)
    {
        if (!m_uploaded)
            throw std::logic_error(
                "CudaTransport::runBatch — call uploadGrid() first");
        if (particles.empty()) return;

        cuda_detail::cudaLaunchWoodcock(
            m_bufs,
            particles.data(),
            static_cast<uint64_t>(particles.size()),
            m_geom,
            m_rngOffset);

        m_rngOffset += static_cast<uint64_t>(particles.size());
    }

    /**
     * @brief Downloads accumulated per-voxel energy from the device and merges
     *        it into the host grid's EnergyScore accumulators (additive).
     *        Clears the device buffer so the next batch starts from zero.
     *
     * @throws std::logic_error if uploadGrid() has not been called first.
     */
    void downloadScores()
    {
        if (!m_uploaded)
            throw std::logic_error(
                "CudaTransport::downloadScores — call uploadGrid() first");

        const auto nVoxels = m_grid.size();
        std::vector<double> hEnergy(nVoxels);
        cuda_detail::cudaDownloadEnergyScored(m_bufs, hEnergy.data());

        for (std::size_t i = 0; i < nVoxels; ++i) {
            if (hEnergy[i] > 0.0)
                m_grid.m_data[i].energyScored.scoreEnergy(hEnergy[i]);
        }
    }

    // ── Static helpers ────────────────────────────────────────────────────────

    /// @brief Converts an xraymc::Particle to the CudaParticle layout.
    static Particle toCudaParticle(const xraymc::Particle& p) noexcept
    {
        return { p.pos[0], p.pos[1], p.pos[2],
                 p.dir[0], p.dir[1], p.dir[2],
                 p.energy, p.weight };
    }

private:
    // ── Material flattening ───────────────────────────────────────────────────

    // Resamples all CPU-side cubic-spline tables into fixed log-spaced arrays
    // suitable for device binary-search interpolation.
    static void flattenMaterial(const Material<16>& src,
                                 cuda_detail::CudaMaterial& dst)
    {
        constexpr int N = cuda_detail::CUDA_TABLE_N;
        constexpr int M = cuda_detail::CUDA_FORM_N;

        // ── Attenuation tables — log-spaced energy grid 1–151 keV ─────────────
        const double logEmin = std::log(1.0);
        const double logEmax = std::log(151.0);
        const double eStep   = (logEmax - logEmin) / (N - 1);

        for (int i = 0; i < N; ++i) {
            const double logE = logEmin + i * eStep;
            const double E    = std::exp(logE);
            const auto   av   = src.attenuationValues(E);

            dst.logEnergy[i]        = static_cast<float>(logE);
            dst.logPhotoelectric[i] = static_cast<float>(
                std::log(av.photoelectric + 1.0e-30));
            dst.logIncoherent[i]    = static_cast<float>(
                std::log(av.incoherent    + 1.0e-30));
            dst.logCoherent[i]      = static_cast<float>(
                std::log(av.coherent      + 1.0e-30));
        }
        dst.nPoints = N;

        // ── Form-factor and scatter-function tables ───────────────────────────
        // Momentum-transfer range: 1e-3 to qmax(151 keV) ≈ 12.18 Å⁻¹.
        // S(q→0) = 0 and F(q→0) = Z_eff, so the low-q clamp is safe.
        const double logQmin = std::log(1.0e-3);
        const double logQmax = std::log(Material<16>::momentumTransferMax(151.0));
        const double qStep   = (logQmax - logQmin) / (M - 1);

        for (int i = 0; i < M; ++i) {
            const double logQ = logQmin + i * qStep;
            const double q    = std::exp(logQ);
            const double F    = src.formFactor(q);
            const double S    = src.scatterFactor(q);

            dst.logMomTransfer[i]  = static_cast<float>(logQ);
            dst.logFormFactorSq[i] = static_cast<float>(
                std::log(F * F + 1.0e-30));
            dst.logScatterFunc[i]  = static_cast<float>(
                std::log(S + 1.0e-30));
        }
        dst.nFormPoints = M;

        // F²(q→0): evaluate at a very small but non-zero q.
        const double F0        = src.formFactor(1.0e-5);
        dst.formFactorSqAt0    = static_cast<float>(F0 * F0);
        dst.effectiveZ         = static_cast<float>(src.effectiveZ());
    }

    // ── Woodcock majorant table flattening ────────────────────────────────────

    // m_woodcockStepTableLin stores {E_linear [keV], maxAtt_linear [1/cm]} pairs
    // sorted by energy. Resamples them onto a log-spaced GPU table.
    static void flattenWoodcockTable(
        const std::vector<std::pair<double, double>>& src,
        cuda_detail::CudaWoodcockTable& dst)
    {
        assert(!src.empty());
        constexpr int N = cuda_detail::CUDA_TABLE_N;

        const double logEmin = std::log(src.front().first);
        const double logEmax = std::log(src.back().first);
        const double step    = (logEmax - logEmin) / (N - 1);

        for (int i = 0; i < N; ++i) {
            const double logE = logEmin + i * step;
            const double E    = std::exp(logE);

            // Linear interpolation on the source table (values are in linear space).
            auto it = std::lower_bound(
                src.begin(), src.end(),
                std::make_pair(E, -1.0e30),
                [](const auto& a, const auto& b) { return a.first < b.first; });

            double att;
            if (it == src.end())        att = src.back().second;
            else if (it == src.begin()) att = src.front().second;
            else {
                const auto prev = std::prev(it);
                const double t  = (E - prev->first) / (it->first - prev->first);
                att = prev->second + t * (it->second - prev->second);
            }

            dst.logEnergy[i] = static_cast<float>(logE);
            dst.logMaxAtt[i] = static_cast<float>(std::log(att + 1.0e-30));
        }
        dst.nPoints = N;
    }

    // ── Geometry struct builder ───────────────────────────────────────────────

    cuda_detail::CudaGridGeom buildGeom() const noexcept
    {
        const auto& aabb = m_grid.AABB();
        const auto& sp   = m_grid.spacing();
        const auto& dim  = m_grid.m_dim; // private — accessible via friend

        cuda_detail::CudaGridGeom g;
        g.nx = static_cast<uint32_t>(dim[0]);
        g.ny = static_cast<uint32_t>(dim[1]);
        g.nz = static_cast<uint32_t>(dim[2]);
        g.aabbXmin = aabb[0]; g.aabbYmin = aabb[1]; g.aabbZmin = aabb[2];
        g.aabbXmax = aabb[3]; g.aabbYmax = aabb[4]; g.aabbZmax = aabb[5];
        g.invSpacingX = 1.0 / sp[0];
        g.invSpacingY = 1.0 / sp[1];
        g.invSpacingZ = 1.0 / sp[2];
        return g;
    }

    // ── Data members ──────────────────────────────────────────────────────────
    Grid&                          m_grid;
    cuda_detail::CudaDeviceBuffers m_bufs     {};
    cuda_detail::CudaGridGeom      m_geom     {};
    uint64_t                       m_rngOffset = 0;
    bool                           m_uploaded  = false;
};

} // namespace xraymc

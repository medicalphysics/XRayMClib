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

#pragma once

#include "xraymc/interactions.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/worldintersectionresult.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {

template <std::size_t NMaterialShells = 16, int LOWENERGYCORRECTION = 2, std::uint_fast8_t TRANSPARENTVOXELS = 255>
class AAVoxelGrid {
public:
    AAVoxelGrid()
    {
        // setting default constructor up with dummy data
        std::array<std::size_t, 3> dim = { 1, 1, 1 };
        std::vector<double> densities(1, 1);
        std::vector<std::uint8_t> mIdx(1, 0);
        std::vector<Material<NMaterialShells>> materials;

        auto air_cand = Material<NMaterialShells>::byNistName("Air, Dry (near sea level)");
        materials.push_back(air_cand.value());
        setData(dim, densities, mIdx, materials);
        setSpacing({ 1, 1, 1 });
    }

    AAVoxelGrid(const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing, const std::vector<double>& density, const std::vector<uint8_t>& materialIdx, const std::vector<Material<NMaterialShells>>& materials)
    {
        setData(dim, density, materialIdx, materials);
        setSpacing(spacing);
    }

    AAVoxelGrid(const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing, const std::vector<double>& ctNumbers)
    {
        setData(dim, ctNumbers);
        setSpacing(spacing);
    }

    bool setData(const std::array<std::size_t, 3>& dim, const std::vector<double>& density, const std::vector<uint8_t>& materialIdx, const std::vector<Material<NMaterialShells>>& materials)
    {
        const auto size = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
        if (density.size() != size || materialIdx.size() != size) {
            return false;
        }
        // finding max material idx;
        auto maxMatIdx = std::max_element(std::execution::par_unseq, materialIdx.cbegin(), materialIdx.cend());
        if (*maxMatIdx >= materials.size()) {
            return false;
        }
        m_dim = dim;
        m_data.resize(size);
        m_dose.clear();
        m_dose.resize(size);

        std::transform(std::execution::par_unseq, density.cbegin(), density.cend(), materialIdx.cbegin(), m_data.begin(), [](const auto d, const auto mIdx) -> DataElement {
            return { .energyScored = EnergyScore {}, .density = d, .materialIndex = mIdx };
        });

        m_materials = materials;
        generateWoodcockStepTable();

        updateAABB();
        return true;
    }

    bool setData(const std::array<std::size_t, 3>& dim, const std::vector<double>& ctNumbers)
    {
        const auto size = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
        m_dim = dim;
        // allocating
        m_data.clear();
        m_dose.clear();
        m_data.resize(size);
        m_dose.resize(size);

        const auto dens = shcneiderMaterialDensities(ctNumbers);
        const auto matIdx = shcneiderMaterialIndices(ctNumbers);
        for (std::size_t i = 0; i < dens.size(); ++i) {
            m_data[i].density = dens[i];
            m_data[i].materialIndex = matIdx[i];
        }

        // generate materials
        m_materials.clear();
        const auto shcneiderWeights = shcneiderMaterialWeights();
        m_materials.reserve(shcneiderWeights.size());
        for (const auto& w : shcneiderWeights)
            m_materials.push_back(Material<NMaterialShells>::byWeight(w).value());

        generateWoodcockStepTable();

        updateAABB();
        return true;
    }

    static std::array<std::string, 24> shcneiderMaterialNames()
    {
        std::array<std::string, 24> names;
        names[0] = "Air";
        names[1] = "Lung";
        for (std::size_t i = 2; i < 8; ++i)
            names[i] = "Soft " + std::to_string(i - 2);
        for (std::size_t i = 8; i < names.size(); ++i)
            names[i] = "Skeletal " + std::to_string(i - 8);
        return names;
    }

    static std::array<std::map<std::size_t, double>, 24> shcneiderMaterialWeights()
    {
        std::array<std::map<std::size_t, double>, 24> sw;
        sw[0] = { { 1, 75.5 }, { 6, 23.2 }, { 7, 1.3 } };
        sw[1] = { { 1, 10.3 }, { 6, 10.5 }, { 7, 3.1 }, { 8, 74.9 }, { 11, 0.2 }, { 15, 0.2 }, { 16, 0.3 }, { 17, 0.3 }, { 18, 0.2 } };
        sw[2] = { { 1, 11.6 }, { 6, 68.1 }, { 7, 0.2 }, { 8, 19.8 }, { 11, 0.1 }, { 15, 0.1 }, { 16, 0.1 } };
        sw[3] = { { 1, 11.3 }, { 6, 56.7 }, { 7, 0.9 }, { 8, 30.8 }, { 11, 0.1 }, { 15, 0.1 }, { 16, 0.1 } };
        sw[4] = { { 1, 11 }, { 6, 45.8 }, { 7, 1.5 }, { 8, 41.1 }, { 11, 0.1 }, { 15, 0.1 }, { 16, 0.2 }, { 17, 0.2 } };
        sw[5] = { { 1, 10.8 }, { 6, 35.6 }, { 7, 2.2 }, { 8, 50.9 }, { 11, 0.1 }, { 15, 0.2 }, { 16, 0.2 } };
        sw[6] = { { 1, 10.6 }, { 6, 28.4 }, { 7, 2.6 }, { 8, 57.8 }, { 11, 0.1 }, { 15, 0.2 }, { 16, 0.2 }, { 17, 0.1 } };
        sw[7] = { { 1, 10.3 }, { 6, 13.4 }, { 7, 3 }, { 8, 72.3 }, { 11, 0.2 }, { 15, 0.2 }, { 16, 0.2 }, { 17, 0.2 } };
        sw[8] = { { 1, 9.4 }, { 6, 20.7 }, { 7, 6.2 }, { 8, 62.2 }, { 12, 0.6 }, { 16, 0.6 }, { 17, 0.3 } };
        sw[9] = { { 1, 9.5 }, { 6, 45.5 }, { 7, 2.5 }, { 8, 35.5 }, { 11, 0.1 }, { 15, 2.1 }, { 16, 0.1 }, { 17, 0.1 }, { 18, 0.1 }, { 20, 4.5 } };
        sw[10] = { { 1, 8.9 }, { 6, 42.3 }, { 7, 2.7 }, { 8, 36.3 }, { 11, 0.1 }, { 15, 3 }, { 16, 0.1 }, { 17, 0.1 }, { 18, 0.1 }, { 20, 6.4 } };
        sw[11] = { { 1, 8.2 }, { 6, 39.1 }, { 7, 2.9 }, { 8, 37.2 }, { 11, 0.1 }, { 15, 3.9 }, { 16, 0.1 }, { 17, 0.1 }, { 18, 0.1 }, { 20, 8.3 } };
        sw[12] = { { 1, 7.6 }, { 6, 36.1 }, { 7, 3 }, { 8, 38 }, { 11, 0.1 }, { 15, 4.7 }, { 16, 0.2 }, { 17, 0.1 }, { 20, 10.1 } };
        sw[13] = { { 1, 7.1 }, { 6, 33.5 }, { 7, 3.2 }, { 8, 38.7 }, { 11, 0.1 }, { 15, 5.4 }, { 16, 0.2 }, { 20, 11.7 } };
        sw[14] = { { 1, 6.6 }, { 6, 31 }, { 7, 3.3 }, { 8, 39.4 }, { 11, 0.1 }, { 15, 6.1 }, { 16, 0.2 }, { 20, 13.2 } };
        sw[15] = { { 1, 6.1 }, { 6, 28.7 }, { 7, 3.5 }, { 8, 40 }, { 11, 0.1 }, { 15, 6.7 }, { 16, 0.2 }, { 20, 14.6 } };
        sw[16] = { { 1, 5.6 }, { 6, 26.5 }, { 7, 3.6 }, { 8, 40.5 }, { 11, 0.1 }, { 15, 7.3 }, { 16, 0.3 }, { 20, 15.9 } };
        sw[17] = { { 1, 5.2 }, { 6, 24.6 }, { 7, 3.7 }, { 8, 41.1 }, { 11, 0.1 }, { 15, 7.8 }, { 16, 0.3 }, { 20, 17 } };
        sw[18] = { { 1, 4.9 }, { 6, 22.7 }, { 7, 3.8 }, { 8, 41.6 }, { 11, 0.1 }, { 15, 8.3 }, { 16, 0.3 }, { 20, 18.1 } };
        sw[19] = { { 1, 4.5 }, { 6, 21 }, { 7, 3.9 }, { 8, 42 }, { 11, 0.1 }, { 15, 8.8 }, { 16, 0.3 }, { 20, 19.2 } };
        sw[20] = { { 1, 4.2 }, { 6, 19.4 }, { 7, 4 }, { 8, 42.5 }, { 11, 0.1 }, { 15, 9.2 }, { 16, 0.3 }, { 20, 20.1 } };
        sw[21] = { { 1, 3.9 }, { 6, 17.9 }, { 7, 4.1 }, { 8, 42.9 }, { 11, 0.1 }, { 15, 9.6 }, { 16, 0.3 }, { 20, 21 } };
        sw[22] = { { 1, 3.6 }, { 6, 16.5 }, { 7, 4.2 }, { 8, 43.2 }, { 11, 0.1 }, { 15, 10 }, { 16, 0.3 }, { 20, 21.9 } };
        sw[23] = { { 1, 3.4 }, { 6, 15.5 }, { 7, 4.2 }, { 8, 43.5 }, { 11, 0.1 }, { 15, 10.3 }, { 16, 0.3 }, { 20, 22.5 } };
        return sw;
    }
    static std::vector<double> shcneiderMaterialDensities(const std::vector<double>& ctNumbers)
    {
        std::vector<double> dens(ctNumbers.size());
        std::transform(std::execution::par_unseq, ctNumbers.cbegin(), ctNumbers.cend(), dens.begin(), [](double H) {
            double d;
            if (H < -98.0) {
                d = H / 1000.0 + 1; // linear
            } else if (H < 14.0) {
                d = (1.018 + 0.000893 * H);
            } else if (H < 23) {
                d = 1.03;
            } else if (H < 100) {
                d = 1.003 + 0.001169 * H;
            } else {
                d = 1.017 + 0.000592 * H;
            }
            return d;
        });
        return dens;
    }
    static std::vector<std::uint8_t> shcneiderMaterialIndices(const std::vector<double>& ctNumbers)
    {
        constexpr std::array<double, 24> schneiderUlim = { -950.0, -120, -83, -53, -23, 7, 18, 80, 120, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600 };
        std::vector<std::uint8_t> m(ctNumbers.size());
        std::transform(std::execution::par_unseq, ctNumbers.cbegin(), ctNumbers.cend(), m.begin(), [&](double H) {
            const auto it = std::upper_bound(schneiderUlim.cbegin(), schneiderUlim.cend(), H);
            const auto idx = it != schneiderUlim.cend() ? std::distance(schneiderUlim.cbegin(), it) : schneiderUlim.size() - 1;
            return static_cast<std::uint8_t>(idx);
        });
        return m;
    }

    void flipAxis(std::size_t axis)
    {
        auto data = m_data;
        auto dose = m_dose;
        for (std::size_t z = 0; z < m_dim[2]; ++z) {
            const auto zz = axis != 2 ? z : m_dim[2] - z - 1;
            for (std::size_t y = 0; y < m_dim[1]; ++y) {
                const auto yy = axis != 1 ? y : m_dim[1] - y - 1;
                for (std::size_t x = 0; x < m_dim[0]; ++x) {
                    const auto fIdx = x + m_dim[0] * (y + m_dim[1] * z);
                    const auto xx = axis != 0 ? x : m_dim[0] - x - 1;
                    const auto tIdx = xx + m_dim[0] * (yy + m_dim[1] * zz);
                    data[tIdx] = m_data[fIdx];
                    dose[tIdx] = m_dose[fIdx];
                }
            }
        }
        m_data = data;
        m_dose = dose;
    }

    void rollAxis(std::size_t from, std::size_t to)
    {
        if (from > 2 || to > 2 || from == to)
            return;

        const auto dim = m_dim;

        std::swap(m_dim[from], m_dim[to]);
        std::swap(m_invSpacing[from], m_invSpacing[to]);
        std::swap(m_spacing[from], m_spacing[to]);
        std::swap(m_aabb[from], m_aabb[to]);
        std::swap(m_aabb[from + 3], m_aabb[to + 3]);

        std::array<std::size_t, 3> swapped { 0, 1, 2 };
        std::swap(swapped[from], swapped[to]);

        auto data = m_data; // making copy
        auto dose = m_dose;

        for (std::size_t z = 0; z < dim[2]; ++z)
            for (std::size_t y = 0; y < dim[1]; ++y)
                for (std::size_t x = 0; x < dim[0]; ++x) {
                    const auto fIdx = x + dim[0] * (y + dim[1] * z);
                    const std::array v = { x, y, z };
                    const std::array w = { v[swapped[0]], v[swapped[1]], v[swapped[2]] };
                    const auto tIdx = w[0] + m_dim[0] * (w[1] + m_dim[1] * w[2]);
                    data[tIdx] = m_data[fIdx];
                    dose[tIdx] = m_dose[fIdx];
                }
        m_data = data;
        m_dose = dose;
    }

    std::size_t size() const
    {
        return m_dim[0] * m_dim[1] * m_dim[2];
    }

    void setSpacing(const std::array<double, 3>& spacing)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_spacing[i] = std::abs(spacing[i]);
            m_invSpacing[i] = 1 / m_spacing[i];
        }
        updateAABB();
    }

    const std::array<double, 3>& spacing() const
    {
        return m_spacing;
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<double, 3> center() const
    {
        std::array<double, 3> center;
        for (std::size_t i = 0; i < 3; ++i) {
            center[i] = (m_aabb[i] + m_aabb[i + 3]) / 2;
        }
        return center;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    inline std::size_t flatIndex(const std::array<std::size_t, 3>& index) const
    {
        return flatIndex(index[0], index[1], index[2]);
    }

    inline std::size_t flatIndex(const std::size_t x, const std::size_t y, const std::size_t z) const
    {
        return x + m_dim[0] * (y + m_dim[1] * z);
    }

    template <bool BOUNDSCHECK = true>
    inline std::size_t flatIndex(const std::array<double, 3>& pos) const
    {
        return flatIndex(index<BOUNDSCHECK>(pos));
    }

    template <bool BOUNDSCHECK = true>
    inline std::array<std::size_t, 3> index(const std::array<double, 3>& pos) const
    {
        if constexpr (BOUNDSCHECK) {
            const std::array<std::size_t, 3> idx = {
                static_cast<std::size_t>(std::clamp((pos[0] - m_aabb[0]) * m_invSpacing[0], 0.0, static_cast<double>(m_dim[0] - 1))),
                static_cast<std::size_t>(std::clamp((pos[1] - m_aabb[1]) * m_invSpacing[1], 0.0, static_cast<double>(m_dim[1] - 1))),
                static_cast<std::size_t>(std::clamp((pos[2] - m_aabb[2]) * m_invSpacing[2], 0.0, static_cast<double>(m_dim[2] - 1)))
            };
            return idx;
        } else {
            const std::array<std::size_t, 3> idx = {
                static_cast<std::size_t>((pos[0] - m_aabb[0]) * m_invSpacing[0]),
                static_cast<std::size_t>((pos[1] - m_aabb[1]) * m_invSpacing[1]),
                static_cast<std::size_t>((pos[2] - m_aabb[2]) * m_invSpacing[2])
            };
            return idx;
        }
    }

    inline std::array<std::size_t, 3> index(const std::size_t flatIndex) const
    {
        const auto z = flatIndex / (m_dim[0] * m_dim[1]);
        const auto y = (flatIndex - z * m_dim[0] * m_dim[1]) / m_dim[0];
        const auto x = flatIndex - z * m_dim[0] * m_dim[1] - y * m_dim[0];
        std::array<std::size_t, 3> arr = { x, y, z };
        return arr;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        auto inter = basicshape::AABB::intersect(p, m_aabb);
        if constexpr (TRANSPARENTVOXELS != 255) {
            if (inter.valid()) {
                voxelIntersect<WorldIntersectionResult, TRANSPARENTVOXELS>(p, inter);
            }
        }
        return inter;
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        auto res = basicshape::AABB::template intersectVisualization<U>(p, m_aabb);
        if constexpr (TRANSPARENTVOXELS != 255) {
            if (res.valid()) {
                res.normal = { 0, 0, 0 };
                voxelIntersect<VisualizationIntersectionResult<U>, TRANSPARENTVOXELS>(p, res);
            }
        } else {
            if (res.valid()) {
                res.normal = { 0, 0, 0 };
                voxelIntersect<VisualizationIntersectionResult<U>, 0>(p, res);
            }
        }
        return res;
    }

    std::vector<std::uint8_t> getMaterialIndex() const
    {
        const auto size = m_dim[0] * m_dim[1] * m_dim[2];
        std::vector<std::uint8_t> i(size);
        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), i.begin(), [](const auto& d) { return d.materialIndex; });
        return i;
    }

    std::vector<double> getDensity() const
    {
        const auto size = m_dim[0] * m_dim[1] * m_dim[2];
        std::vector<double> i(size);
        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), i.begin(), [](const auto& d) { return d.density; });
        return i;
    }

    std::vector<EnergyScore> getEnergyScores() const
    {
        const auto size = m_dim[0] * m_dim[1] * m_dim[2];
        std::vector<EnergyScore> i(size);
        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), i.begin(), [](const auto& d) { return d.energyScored; });
        return i;
    }

    const std::vector<DoseScore>& getDoseScores() const
    {
        return m_dose;
    }

    const EnergyScore& energyScored(std::size_t flatIndex = 0) const
    {
        return m_data.at(flatIndex).energyScored;
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto size = m_dim[0] * m_dim[1] * m_dim[2];
        const auto voxel_volume = m_spacing[0] * m_spacing[1] * m_spacing[2];

        for (std::size_t i = 0; i < size; ++i) {
            const auto& ei = m_data[i].energyScored;
            m_dose[i].addScoredEnergy(ei, voxel_volume, m_data[i].density, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t flatIndex = 0) const
    {
        return m_dose.at(flatIndex);
    }

    void clearEnergyScored()
    {
        std::for_each(std::execution::par_unseq, m_data.begin(), m_data.end(), [](auto& d) { d.energyScored.clear(); });
    }

    void clearDoseScored()
    {
        std::for_each(std::execution::par_unseq, m_dose.begin(), m_dose.end(), [](auto& d) { d.clear(); });
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        if constexpr (TRANSPARENTVOXELS != 255) {
            voxelTransport<TRANSPARENTVOXELS>(p, state);
        } else {
            woodcockTransport(p, state);
        }
    }

    double maxAttenuationValue(const double energy) const
    {
        return interpolate(m_woodcockStepTableLin, energy);
    }

protected:
    void updateAABB()
    {
        const auto c = center();
        for (std::size_t i = 0; i < 3; ++i) {
            const auto half_dist = m_dim[i] * 0.5 * m_spacing[i];
            m_aabb[i] = -half_dist;
            m_aabb[i + 3] = half_dist;
        }
        translate(c);
    }

    static inline std::uint_fast8_t argmin3(const std::array<double, 3>& a)
    {
        return a[0] < a[1] ? a[0] < a[2] ? 0 : 2 : a[1] < a[2] ? 1
                                                               : 2;
    }
    static inline std::uint_fast8_t argmax3(const std::array<double, 3>& a)
    {
        return a[0] > a[1] ? a[0] > a[2] ? 0 : 2 : a[1] > a[2] ? 1
                                                               : 2;
    }

    template <typename Intersection = WorldIntersectionResult, std::uint_fast8_t IGNOREIDX = 255>
    void voxelIntersect(const ParticleType auto& p, Intersection& intersection) const
    {
        std::array<std::size_t, 3> xyz;

        if (intersection.rayOriginIsInsideItem) {
            xyz = index<true>(p.pos);
        } else {
            // make sure we are well inside a voxel
            const auto t = intersection.intersection + GEOMETRIC_ERROR();
            const std::array<double, 3> pos = {
                p.pos[0] + p.dir[0] * t,
                p.pos[1] + p.dir[1] * t,
                p.pos[2] + p.dir[2] * t
            };
            xyz = index<true>(pos);
        }

        auto index_flat = flatIndex(xyz);
        if (m_data[index_flat].materialIndex != IGNOREIDX) {
            if constexpr (!std::is_same_v<Intersection, WorldIntersectionResult>) {
                intersection.value = m_dose[index_flat].dose();
                std::array<double, 3> tMax;
                for (std::size_t i = 0; i < 3; ++i) {
                    if (std::abs(p.dir[i]) > GEOMETRIC_ERROR()) {
                        const auto plane = p.dir[i] > 0 ? m_aabb[i] + xyz[i] * m_spacing[i] : m_aabb[i] + (xyz[i] + 1) * m_spacing[i];
                        tMax[i] = (plane - (p.pos[i])) / p.dir[i];
                    } else {
                        tMax[i] = std::numeric_limits<double>::max();
                    }
                }
                const auto planeIdx = argmax3(tMax);
                intersection.normal = { 0, 0, 0 };
                intersection.normal[planeIdx] = p.dir[planeIdx] < 0 ? 1 : -1;
            }
            return;
        }

        std::array<double, 3> tMax;
        for (std::size_t i = 0; i < 3; ++i) {
            if (std::abs(p.dir[i]) > GEOMETRIC_ERROR()) {
                const auto plane = p.dir[i] < 0 ? m_aabb[i] + xyz[i] * m_spacing[i] : m_aabb[i] + (xyz[i] + 1) * m_spacing[i];
                tMax[i] = (plane - (p.pos[i])) / p.dir[i];
            } else {
                tMax[i] = std::numeric_limits<double>::max();
            }
        }

        intersection.intersectionValid = false;
        bool still_inside;
        do {
            const auto planeIdx = argmin3(tMax);
            xyz[planeIdx] += p.dir[planeIdx] < 0 ? -1 : 1;
            if (xyz[planeIdx] < m_dim[planeIdx]) {
                still_inside = m_data[flatIndex(xyz)].materialIndex == IGNOREIDX;
                if (still_inside) {
                    tMax[planeIdx] += m_spacing[planeIdx] / std::abs(p.dir[planeIdx]);
                } else {
                    intersection.intersection = std::min({ tMax[0], tMax[1], tMax[2] });
                    intersection.rayOriginIsInsideItem = false;
                    intersection.intersectionValid = true;
                    if constexpr (!std::is_same_v<Intersection, WorldIntersectionResult>) {
                        intersection.value = m_dose[flatIndex(xyz)].dose();
                        intersection.normal = { 0, 0, 0 };
                        intersection.normal[planeIdx] = p.dir[planeIdx] < 0 ? 1 : -1;
                    }
                }
            } else {
                still_inside = false;
            }
        } while (still_inside);
    }

    template <std::uint_fast8_t IGNOREIDX = 255>
    void voxelTransport(ParticleType auto& p, RandomState& state)
    {

        // Assume valid voxel when we start

        std::array<std::size_t, 3> xyz = index<true>(p.pos);

        do {
            std::array<double, 3> tstep;
            for (std::size_t i = 0; i < 3; ++i) {
                if (std::abs(p.dir[i] > GEOMETRIC_ERROR())) {
                    const auto plane = p.dir[i] < 0 ? m_aabb[i] + xyz[i] * m_spacing[i] : m_aabb[i] + (xyz[i] + 1) * m_spacing[i];
                    tstep[i] = (plane - p.pos[i]) / p.dir[i];
                } else {
                    tstep[i] = std::numeric_limits<double>::max();
                }
            }
            auto nextPlane = argmin3(tstep);
            // probability
            const auto index_flat = flatIndex(xyz);
            const auto matInd = m_data[index_flat].materialIndex;
            if (matInd == IGNOREIDX) {
                return;
            }
            const auto density = m_data[index_flat].density;
            const auto att = m_materials[matInd].attenuationValues(p.energy);

            const auto prob = std::exp(-att.sum() * density * tstep[nextPlane]);
            const auto rand = state.randomUniform();
            if (rand < prob) {
                // no interaction, we translate to next plane
                xyz[nextPlane] = xyz[nextPlane] + (p.dir[nextPlane] < 0 ? -1 : 1);
                if (xyz[nextPlane] >= m_dim[nextPlane]) {
                    // We are outside of the volume, we translate to the border and exit
                    p.border_translate(tstep[nextPlane]);
                    return;
                } else {
                    p.translate(tstep[nextPlane]);
                }
            } else {
                // interaction
                const auto step_correction = -std::log(rand) / (att.sum() * density);
                p.translate(step_correction);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_materials[matInd], state);
                auto& energyScored = m_data[index_flat].energyScored;
                energyScored.scoreEnergy(intRes.energyImparted);
                if (!intRes.particleAlive) {
                    return;
                }
            }
        } while (true);
    }

    void generateWoodcockStepTable()
    {
        std::vector<double> energy;
        {
            auto e = std::log(MIN_ENERGY());
            const auto emax = std::log(MAX_ENERGY());
            const auto estep = (emax - e) / 10;
            while (e < emax) {
                energy.push_back(e);
                e += estep;
            }
            energy.push_back(emax);
        }
        // adding edges;
        for (const auto& mat : m_materials) {
            for (std::size_t i = 0; i < mat.numberOfShells(); ++i) {
                const auto& shell = mat.shell(i);
                const auto e = shell.bindingEnergy + 0.01;
                if (e > MIN_ENERGY()) {
                    energy.push_back(std::log(e));
                }
            }
        }
        std::sort(energy.begin(), energy.end());
        auto remove = std::unique(energy.begin(), energy.end());
        energy.erase(remove, energy.end());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), energy.begin(), [](const auto e) { return std::exp(e); });

        // finding max density for each material;
        std::vector<double> dens(m_materials.size(), 0.0);
        for (const auto& d : m_data) {
            const auto i = d.materialIndex;
            dens[i] = std::max(d.density, dens[i]);
        }

        // finding max attenuation for each energy
        std::vector<double> att(energy.size(), 0.0);
        for (std::size_t mIdx = 0; mIdx < m_materials.size(); ++mIdx) {
            const auto& mat = m_materials[mIdx];
            const auto d = dens[mIdx];
            for (std::size_t i = 0; i < energy.size(); ++i) {
                const auto aval = mat.attenuationValues(energy[i]);
                att[i] = std::max(aval.sum() * d, att[i]);
            }
        }
        std::vector<std::pair<double, double>> data(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.cbegin(), data.begin(), [](const auto e, const auto a) {
            return std::make_pair(e, a);
        });

        m_woodcockStepTableLin = data;
    }

    void woodcockTransport(ParticleType auto& p, RandomState& state)
    {
        bool valid = basicshape::AABB::pointInside(p.pos, m_aabb);
        bool updateAtt = true;

        double attMaxInv;
        while (valid) {
            if (updateAtt) {
                attMaxInv = 1 / maxAttenuationValue(p.energy);
                updateAtt = false;
            }
            const auto steplen = -std::log(state.randomUniform()) * attMaxInv;

            const auto intersection = basicshape::AABB::intersect(p, m_aabb);

            if (steplen < intersection.intersection) {
                p.translate(steplen);
                const auto flat_index = flatIndex<false>(p.pos);
                const auto matIdx = m_data[flat_index].materialIndex;
                const auto dens = m_data[flat_index].density;
                const auto att = m_materials[matIdx].attenuationValues(p.energy);
                const auto attTot = att.sum() * dens;
                // check if real or virtual interaction
                if (state.randomUniform() < attTot * attMaxInv) {
                    const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_materials[matIdx], state);
                    m_data[flat_index].energyScored.scoreEnergy(intRes.energyImparted);
                    valid = intRes.particleAlive;
                    updateAtt = intRes.particleEnergyChanged;
                }
            } else {
                p.border_translate(intersection.intersection);
                valid = false;
            }
        }
    }

    struct DataElement {
        EnergyScore energyScored;
        double density = 0;
        std::uint8_t materialIndex = 0;
    };

private:
    std::array<std::size_t, 3> m_dim = { 1, 1, 1 };
    std::array<double, 3> m_invSpacing = { 1, 1, 1 };
    std::array<double, 3> m_spacing = { 1, 1, 1 };
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<std::pair<double, double>> m_woodcockStepTableLin;
    std::vector<DataElement> m_data;
    std::vector<DoseScore> m_dose;
    std::vector<Material<NMaterialShells>> m_materials;
};
}
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

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <array>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, std::uint_fast8_t TRANSPARENTVOXELS = 255>
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
        m_dose.resize(size);
        EnergyScore dummy_dose;
        std::transform(std::execution::par_unseq, density.cbegin(), density.cend(), materialIdx.cbegin(), m_data.begin(), [=](const auto d, const auto mIdx) -> DataElement {
            return { .energyScored = dummy_dose, .density = d, .materialIndex = mIdx };
        });

        m_materials = materials;
        generateWoodcockStepTable();

        updateAABB();
        return true;
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
            const auto half_dist = (m_dim[i] * 0.5) * m_spacing[i];
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

    template <typename Intersection = WorldIntersectionResult, std::uint_fast8_t IGNOREIDX = 255>
    void voxelIntersect(const ParticleType auto& p, Intersection& intersection) const
    {
        std::array<std::size_t, 3> xyz;

        if (intersection.rayOriginIsInsideItem) {
            xyz = index<false>(p.pos);
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
            return;
        }

        const std::array<int, 3> xyz_step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };

        // by using int, max x*y size is about 2e9
        const std::array<int, 3> xyz_step_flat = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -static_cast<int>(m_dim[0]) : static_cast<int>(m_dim[0]),
            p.dir[2] < 0 ? -static_cast<int>(m_dim[0] * m_dim[1]) : static_cast<int>(m_dim[0] * m_dim[1])
        };

        const std::array<double, 3> delta = {
            m_spacing[0] / std::abs(p.dir[0]),
            m_spacing[1] / std::abs(p.dir[1]),
            m_spacing[2] / std::abs(p.dir[2])
        };

        std::array<double, 3> tMax;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto plane = p.dir[i] < 0 ? m_aabb[i] + xyz[i] * m_spacing[i] : m_aabb[i] + (xyz[i] + 1) * m_spacing[i];
            tMax[i] = (plane - (p.pos[i] + p.dir[i] * intersection.intersection)) / p.dir[i];
        }

        bool still_inside = true;
        std::uint_fast8_t dIdx;
        while (still_inside && m_data[index_flat].materialIndex == IGNOREIDX) {
            dIdx = argmin3(tMax);
            xyz[dIdx] += xyz_step[dIdx];
            still_inside = xyz[dIdx] < m_dim[dIdx];
            index_flat += xyz_step_flat[dIdx];
            tMax[dIdx] += delta[dIdx];
        }

        if (still_inside) {
            intersection.intersection += tMax[dIdx] - delta[dIdx];
            intersection.rayOriginIsInsideItem = false;
            if constexpr (!std::is_same<Intersection, WorldIntersectionResult>::value) {
                intersection.normal[dIdx] = p.dir[dIdx] < 0 ? -1 : 1;
                intersection.value = m_dose[index_flat].dose();
            }
            intersection.intersectionValid = true;
        } else {
            intersection.intersectionValid = false;
        }
    }

    template <std::uint_fast8_t IGNOREIDX = 255>
    void voxelTransport(ParticleType auto& p, RandomState& state)
    {
        bool still_inside;
        do {
            std::array<std::size_t, 3> xyz = index<false>(p.pos);

            auto index_flat = flatIndex(xyz);

            const std::array<int, 3> xyz_step = {
                p.dir[0] < 0 ? -1 : 1,
                p.dir[1] < 0 ? -1 : 1,
                p.dir[2] < 0 ? -1 : 1
            };

            const std::array<int, 3> xyz_step_flat = {
                p.dir[0] < 0 ? -1 : 1,
                p.dir[1] < 0 ? -static_cast<int>(m_dim[0]) : static_cast<int>(m_dim[0]),
                p.dir[2] < 0 ? -static_cast<int>(m_dim[0] * m_dim[1]) : static_cast<int>(m_dim[0] * m_dim[1])
            };

            const std::array<double, 3> delta = {
                m_spacing[0] / std::abs(p.dir[0]),
                m_spacing[1] / std::abs(p.dir[1]),
                m_spacing[2] / std::abs(p.dir[2])
            };

            std::array<double, 3> tMax;
            for (std::size_t i = 0; i < 3; ++i) {
                const auto plane = p.dir[i] < 0 ? m_aabb[i] + xyz[i] * m_spacing[i] : m_aabb[i] + (xyz[i] + 1) * m_spacing[i];
                tMax[i] = (plane - p.pos[i]) / p.dir[i];
            }
            // Endre denne??? sjekk voxel intersect
            double tCurrent = 0;

            double interaction_accum = 1;
            const auto interaction_thres = state.randomUniform();

            const auto tLimit = (*basicshape::AABB::intersectForwardInterval(p, m_aabb))[1];

            bool cont = true;
            std::uint_fast8_t dIdx;

            do {
                // current material
                const auto matInd = m_data[index_flat].materialIndex;
                const auto density = m_data[index_flat].density;
                auto& energyScored = m_data[index_flat].energyScored;

                // updating stepping
                dIdx = argmin3(tMax);

                // updating radiological path
                const auto att = m_materials[matInd].attenuationValues(p.energy);
                const auto att_tot = att.sum() * density;
                const auto step = tMax[dIdx] - tCurrent;
                const auto radio_step = std::exp(-att_tot * step);
                interaction_accum *= radio_step;
                if (interaction_accum < interaction_thres) {
                    // interaction happends
                    // Correcting for not traversing the whole voxel (step_correction is negative)
                    const auto step_correction = std::log(interaction_accum / interaction_thres) / att_tot;
                    tMax[dIdx] += step_correction;

                    if (tMax[dIdx] > tLimit) { // interaction happends outside volume, we exits
                        p.border_translate(tLimit);
                        still_inside = false;
                    } else {
                        // translate particle before interaction
                        p.translate(tMax[dIdx]);

                        const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_materials[matInd], state);
                        energyScored.scoreEnergy(intRes.energyImparted);
                        still_inside = intRes.particleAlive;
                        // particle energy or direction has changed, we restart stepping
                        cont = false;
                    }

                } else {
                    xyz[dIdx] += xyz_step[dIdx];
                    index_flat += xyz_step_flat[dIdx];
                    still_inside = xyz[dIdx] < m_dim[dIdx] && m_data[index_flat].materialIndex != IGNOREIDX;
                    if (!still_inside) {
                        // next step is outside, we exits
                        p.border_translate(tMax[dIdx]);
                    } else {
                        tCurrent = tMax[dIdx];
                        tMax[dIdx] += delta[dIdx];
                    }
                }
            } while (cont && still_inside);

        } while (still_inside);
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
                const auto flat_index = flatIndex<true>(p.pos);
                const auto matIdx = m_data[flat_index].materialIndex;
                const auto& dens = m_data[flat_index].density;
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
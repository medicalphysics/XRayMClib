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
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/particletracker.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/basicshapes/cylinder.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <limits>
#include <optional>
#include <vector>

namespace xraymc {

class FlatDetector {
public:
    FlatDetector()
    {
        const auto nPix = m_detector_dimensions[0] * m_detector_dimensions[1];
        m_energyScore.resize(nPix);
        m_doseScore.resize(nPix);
        m_aabb = calculateAABB();
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
        m_aabb = calculateAABB();
    }

    const std::array<double, 3>& center() const
    {
        return m_center;
    }
    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        m_aabb = calculateAABB();
    }

    void setDirectionCosines(const std::array<double, 3>& dir_cosX, const std::array<double, 3>& dir_cosY)
    {
        if (vectormath::dot(dir_cosX, dir_cosY) < 2 * GEOMETRIC_ERROR<>()) {
            m_direction_cosines[0] = vectormath::normalized(dir_cosX);
            m_direction_cosines[1] = vectormath::normalized(dir_cosY);
        }
        m_aabb = calculateAABB();
    }

    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_direction_cosines;
    }

    void setPixelSpacing(const std::array<double, 2>& spacing)
    {
        setPixelSpacing(spacing[0], spacing[1]);
    }
    void setPixelSpacing(double spacingX, double spacingY)
    {
        m_pixel_spacing[0] = std::max(spacingX, 0.00001);
        m_pixel_spacing[1] = std::max(spacingY, 0.00001);
        m_aabb = calculateAABB();
    }
    std::array<double, 2> pixelSpacing() const
    {
        return m_pixel_spacing;
    }

    void setDetectorDimensions(const std::array<std::size_t, 2>& dimensions)
    {
        setDetectorDimensions(dimensions[0], dimensions[1]);
    }
    void setDetectorDimensions(std::size_t dimX, std::size_t dimY)
    {
        constexpr std::size_t min = 1;
        m_detector_dimensions[0] = std::max(dimX, min);
        m_detector_dimensions[1] = std::max(dimY, min);
        m_energyScore.resize(dimX * dimY);
        m_doseScore.resize(dimX * dimY);
        m_aabb = calculateAABB();
    }

    const std::array<std::size_t, 2>& detectorDimensions() const
    {
        return m_detector_dimensions;
    }

    std::array<double, 3> normal() const
    {
        return vectormath::cross(m_direction_cosines[0], m_direction_cosines[1]);
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        const auto dir = vectormath::changeBasis(m_direction_cosines[0], m_direction_cosines[1], normal(), p.dir);
        const auto pos = vectormath::changeBasis(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto denom = dir[2];
        if (std::abs(denom) > GEOMETRIC_ERROR<>()) {
            const auto t = -pos[2] / denom;
            if (t > 0) {
                const auto x = pos[0] + dir[0] * t;
                const auto y = pos[1] + dir[1] * t;
                if (std::abs(x) <= m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5 && std::abs(y) <= m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) {
                    res.intersection = t;
                    res.rayOriginIsInsideItem = false;
                    res.intersectionValid = true;
                }
            }
        }
        return res;
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        VisualizationIntersectionResult<U> res;
        const auto dir = vectormath::changeBasis(m_direction_cosines[0], m_direction_cosines[1], normal(), p.dir);
        const auto pos = vectormath::changeBasis(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto denom = dir[2];
        if (std::abs(denom) > GEOMETRIC_ERROR<>()) {
            const auto t = -pos[2] / denom;
            if (t > 0) {
                const auto x = pos[0] + dir[0] * t;
                const auto y = pos[1] + dir[1] * t;
                if (std::abs(x) <= m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5 && std::abs(y) <= m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) {
                    res.intersection = t;
                    res.rayOriginIsInsideItem = false;
                    res.intersectionValid = true;
                    res.normal = denom < 0 ? normal() : vectormath::scale(normal(), -1.0);
                    const auto indx = (x + m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5) / m_pixel_spacing[0];
                    const auto indy = (y + m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) / m_pixel_spacing[1];
                    const auto flatIdx = static_cast<std::size_t>(std::floor(indx)) * m_detector_dimensions[1] + static_cast<std::size_t>(std::floor(indy));
                    res.value = m_doseScore[flatIdx].dose();
                }
            }
        }
        return res;
    }

    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }

        const auto pos = vectormath::changeBasis(m_direction_cosines[0], m_direction_cosines[1], normal(), vectormath::subtract(p.pos, m_center));
        const auto x = pos[0];
        const auto y = pos[1];
        const auto indx = (x + m_pixel_spacing[0] * m_detector_dimensions[0] * 0.5) / m_pixel_spacing[0];
        const auto indy = (y + m_pixel_spacing[1] * m_detector_dimensions[1] * 0.5) / m_pixel_spacing[1];
        const auto flatIdx = static_cast<std::size_t>(std::floor(indx)) * m_detector_dimensions[1] + static_cast<std::size_t>(std::floor(indy));
        m_energyScore[flatIdx].scoreEnergy(p.energy * p.weight);
        p.energy = 0;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScore[index];
    }

    void clearEnergyScored()
    {
        for (auto& d : m_energyScore) {
            d.clear();
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = m_pixel_spacing[0] * m_pixel_spacing[1];
        constexpr double density = 1.0;
        for (std::size_t i = 0; i < m_doseScore.size(); ++i) {
            m_doseScore[i].addScoredEnergy(m_energyScore[i], volume, density, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_doseScore[index];
    }

    void clearDoseScored()
    {
        for (auto& d : m_doseScore) {
            d.clear();
        }
    }
    const ParticleTracker& particleTracker() const
    {
        return m_tracker;
    }

    ParticleTracker& particleTracker()
    {
        return m_tracker;
    }

protected:
    static std::pair<double, double> minmax(const double a, const double b)
    {
        return a < b ? std::pair { a, b } : std::pair { b, a };
    }

    std::array<double, 6> calculateAABB() const
    {
        const auto normal_vec = normal();
        const auto half_span_x = vectormath::scale(m_direction_cosines[0], m_pixel_spacing[0] * m_detector_dimensions[0] / 2);
        const auto half_span_y = vectormath::scale(m_direction_cosines[1], m_pixel_spacing[1] * m_detector_dimensions[1] / 2);
        const auto span_vecs = std::array { half_span_x, half_span_y, vectormath::scale(normal_vec, GEOMETRIC_ERROR<>()) };
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (const auto& s : span_vecs)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                aabb[i] = std::min(aabb[i], mi);
                aabb[i + 3] = std::max(aabb[i + 3], ma);
            }
        return aabb;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 };
    std::array<double, 2> m_pixel_spacing = { 0.1, 0.1 };
    std::array<std::size_t, 2> m_detector_dimensions = { 128, 128 };
    std::array<std::array<double, 3>, 2> m_direction_cosines = { 1, 0, 0, 0, 1, 0 };
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<EnergyScore> m_energyScore;
    std::vector<DoseScore> m_doseScore;
    ParticleTracker m_tracker;
};
};

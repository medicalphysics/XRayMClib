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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace dxmc {

class FluenceScore {
public:
    FluenceScore(double radius = 16, const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 3>& normal = { 0, 0, 1 })
        : m_center(center)
        , m_radius(radius)
    {
        setPlaneNormal(normal);
        setEnergyStep(1);
    }

    void setEnergyStep(double step)
    {
        m_energy_step = std::max(step, 0.1);
        const auto N = static_cast<std::uint64_t>(MAX_ENERGY() / m_energy_step) + 1;
        m_intensity.resize(N);
        std::fill(m_intensity.begin(), m_intensity.end(), 0);
    }

    void setCenter(const std::array<double, 3>& c)
    {
        m_center = c;
        calculateAABB();
    }

    void setRadius(double r)
    {
        m_radius = std::max(std::abs(r), 0.00001);
        calculateAABB();
    }

    double area() const
    {
        return std::numbers::pi_v<double> * m_radius * m_radius;
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_center[i] += dist[i];
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    const std::array<double, 3>& center() const
    {
        return m_center;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    void setPlaneNormal(const std::array<double, 3>& normal)
    {
        m_normal = normal;
        vectormath::normalize(m_normal);
        calculateAABB();
    }

    std::vector<std::pair<double, std::uint64_t>> getSpecter() const
    {
        std::vector<std::pair<double, std::uint64_t>> spec(m_intensity.size());
        for (std::size_t i = 0; i < m_intensity.size(); ++i) {
            spec[i] = std::make_pair(m_energy_step * i + m_energy_step / 2, m_intensity[i]);
        }
        return spec;
    }

    std::vector<std::pair<double, double>> getFluenceSpecter(std::uint64_t N_primaries = 1) const
    {
        const auto specter = getSpecter();
        std::vector<std::pair<double, double>> fluence(specter.size());
        const auto a = area();
        const auto norm = N_primaries > 0 ? static_cast<double>(N_primaries) : 0.0;
        const auto norm_inv = 1.0 / norm;
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), fluence.begin(), [a, norm_inv](const auto& s) {
            return std::make_pair(s.first, s.second * norm_inv / a);
        });
        return fluence;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const auto aabb_inter = basicshape::AABB::intersectForwardInterval(p, m_aabb);
        return aabb_inter ? intersectDisc(p) : WorldIntersectionResult {};
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        const auto res = intersect(p);
        VisualizationIntersectionResult<U> w;
        if (res.valid()) {
            w.rayOriginIsInsideItem = false;
            w.intersection = res.intersection;
            w.intersectionValid = true;
            w.normal = vectormath::dot(p.dir, m_normal) <= 0 ? m_normal : vectormath::scale(m_normal, -1.0);
        }
        return w;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    void clearEnergyScored()
    {
        return;
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        // not defined for fluence counter
        return;
    }

    const DoseScore doseScored(std::size_t index = 0) const
    {
        return DoseScore {};
    }

    void clearDoseScored()
    {
        m_energyScored.clear();
        std::fill(m_intensity.begin(), m_intensity.end(), std::uint64_t { 0 });
        return;
    }

    void transport(ParticleType auto& particle, RandomState& state)
    {
        // Assuming particle is on the disc
        const auto eIdx = static_cast<std::size_t>(particle.energy / m_energy_step);
        auto counter = std::atomic_ref(m_intensity[eIdx]);
        counter++;
        m_energyScored.scoreEnergy(particle.energy);
        particle.border_translate(0);
    }

protected:
    static inline std::pair<double, double> minmax(double v1, double v2)
    {
        // use own minmax instead of std::minamx due to bug or weird feature of MSVC compiler with /O2
        return v1 <= v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
    }

    void calculateAABB()
    {
        const std::array<std::array<double, 3>, 3> span = {
            vectormath::scale(vectormath::cross(m_normal, { 1, 0, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 1, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 0, 1 }), m_radius)
        };
        m_aabb = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };

        for (const auto& s : span)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                m_aabb[i] = std::min(m_aabb[i], mi);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], ma);
            }
        // ensure a min size;
        for (std::size_t i = 0; i < 3; ++i) {
            if (m_aabb[i + 3] - m_aabb[i] < 2 * GEOMETRIC_ERROR()) {
                m_aabb[i] -= GEOMETRIC_ERROR();
                m_aabb[i + 3] += GEOMETRIC_ERROR();
            }
        }
    }

    WorldIntersectionResult intersectDisc(const ParticleType auto& p) const
    {
        WorldIntersectionResult res;
        const auto D = vectormath::dot(p.dir, m_normal);
        constexpr double minOrt = 1E-6;
        if (D < minOrt && D > -minOrt)
            return res; // dir and normal is orthogonal, we exits
        res.intersection = vectormath::dot(vectormath::subtract(m_center, p.pos), m_normal) / D;
        if (res.intersection <= 0)
            return res;

        // intersection point
        const auto p_int = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));

        // distance from center
        const auto c_dist = vectormath::subtract(m_center, p_int);
        // check if distance from center is less than radius
        if (vectormath::dot(c_dist, c_dist) <= m_radius * m_radius) {
            res.intersectionValid = true;
        }
        return res;
    }

private:
    std::array<double, 3> m_center = { 0, 0, 0 };
    std::array<double, 3> m_normal = { 0, 0, 1 };
    double m_radius = 16;
    double m_energy_step = 1;
    std::array<double, 6> m_aabb;
    std::vector<std::uint64_t> m_intensity;
    EnergyScore m_energyScored;
};
}
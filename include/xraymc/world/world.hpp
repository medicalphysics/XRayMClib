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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "xraymc/interactions.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/particle.hpp"
#include "xraymc/vectormath.hpp"
#include "xraymc/world/basicshapes/aabb.hpp"
#include "xraymc/world/kdtreeflat.hpp"
#include "xraymc/world/worlditems/worlditemtype.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <concepts>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace xraymc {

template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

// Template for at least one type of items
template <WorldItemType F, WorldItemType... Us>
class World {
    static constexpr std::size_t WorldShells() { return 12; }
    using MaterialType = Material<WorldShells()>;

public:
    World()
        : m_fillMaterial(MaterialType::byNistName("Air, Dry (near sea level)").value())
    {
    }

    World(std::size_t reserveNumberOfWorldItems)
        : m_fillMaterial(MaterialType::byNistName("Air, Dry (near sea level)").value())
    {
        m_items.reserve(reserveNumberOfWorldItems);
    }

    void setMaterial(const MaterialType& mat)
    {
        m_fillMaterial = mat;
    }
    void setMaterial(const MaterialType& mat, double dens)
    {
        m_fillMaterial = mat;
        m_fillMaterialDensity = std::abs(dens);
    }
    void setMaterialDensity(double dens)
    {
        m_fillMaterialDensity = std::abs(dens);
    }
    void setMaterialByWeight(const std::map<std::size_t, double>& composition)
    {
        m_fillMaterial = MaterialType::byWeight(composition).value();
    }
    void setMaterialByWeight(const std::map<std::size_t, double>& composition, double dens)
    {
        m_fillMaterial = MaterialType::byWeight(composition).value();
        m_fillMaterialDensity = std::abs(dens);
    }

    const MaterialType& fillMaterial() const { return m_fillMaterial; }
    double fillMaterialDensity() const { return m_fillMaterialDensity; }

    void reserveNumberOfItems(std::size_t size)
    {
        m_items.reserve(size);
    }

    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U& item)
    {
        m_items.push_back(item);
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U& item, std::string_view name)
    {
        m_items.push_back(item);
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U&& item)
    {
        m_items.push_back(std::move(item));
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U&& item, std::string_view name)
    {
        m_items.push_back(std::move(item));
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    template <AnyWorldItemType<F, Us...> U>
    auto& addItem()
    {
        U item;
        m_items.push_back(std::move(item));
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(std::string_view name)
    {
        U item;
        m_items.push_back(std::move(item));
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    const auto& items() const
    {
        return m_items;
    }

    auto& items()
    {
        return m_items;
    }

    const auto& itemNames() const
    {
        return m_item_names;
    }

    auto& itemNames()
    {
        return m_item_names;
    }

    std::vector<std::variant<F, Us...>*> getItemPointers()
    {
        std::vector<std::variant<F, Us...>*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return &v;
        });
        return ptrs;
    }

    std::vector<const std::variant<F, Us...>*> getItemPointers() const
    {
        std::vector<const std::variant<F, Us...>*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return &v;
        });
        return ptrs;
    }

    std::variant<F, Us...>* getItemPointerFromName(std::string_view name)
    {
        for (std::size_t i = 0; i < m_item_names.size(); ++i) {
            if (m_item_names[i].compare(name) == 0) {
                return &m_items[i];
            }
        }
        return nullptr;
    }

    const std::variant<F, Us...>* getItemPointerFromName(std::string_view name) const
    {
        for (std::size_t i = 0; i < m_item_names.size(); ++i) {
            if (m_item_names[i].compare(name) == 0) {
                return &m_items[i];
            }
        }
        return nullptr;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearEnergyScored(); }, v);
        }
    }

    void clearDoseScored()
    {
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearDoseScored(); }, v);
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        for (auto& v : m_items) {
            std::visit([calibration_factor](auto&& arg) { arg.addEnergyScoredToDoseScore(calibration_factor); }, v);
        }
    }

    void build(double AABB_padding = 10)
    {
        auto ptrs = getItemPointers();
        m_kdtree.setData(ptrs);
        m_aabb = m_kdtree.AABB();

        // adding padding
        const auto padding = std::max(AABB_padding, 0.1); // always at least 1 mm padding
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = m_aabb[i] - padding;
            m_aabb[i + 3] = m_aabb[i + 3] + padding;
        }
    }

    void build(const std::array<double, 6>& aabb)
    {
        build();
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = std::min(m_aabb[i], aabb[i]);
            m_aabb[i + 3] = std::max(m_aabb[i + 3], aabb[i + 3]);
        }
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    std::array<double, 3> center() const
    {
        const auto [l, r] = vectormath::splice(m_aabb);
        return vectormath::scale(0.5, vectormath::add(l, r));
    }

    void translate(const std::array<double, 3> dist)
    {
        for (auto& v : m_items)
            std::visit([&dist](auto&& arg) { arg.translate(dist); }, v);

        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
        m_kdtree.translate(dist);
    }

    inline auto intersect(const ParticleType auto& p)
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    inline auto intersectVisualization(const ParticleType auto& p) const
    {
        return m_kdtree.intersectVisualization(p, m_aabb);
    }

    inline bool transportParticleToWorld(ParticleType auto& p) const
    {
        if (!basicshape::AABB::pointInside(p.pos, m_aabb)) {
            const auto t = basicshape::AABB::intersect(p, m_aabb);
            if (t.valid()) {
                p.border_translate(t.intersection);
                return true;
            } else {
                return false;
            }
        }
        return true;
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        bool continueSampling = transportParticleToWorld(p);
        bool updateAttenuation = true;

        double attenuationTotalInv;
        AttenuationValues att;
        while (continueSampling) {
            if (updateAttenuation) {
                att = m_fillMaterial.attenuationValues(p.energy);
                attenuationTotalInv = 1 / (att.sum() * m_fillMaterialDensity);
                updateAttenuation = false;
            }

            const auto r1 = state.randomUniform();
            const auto stepLength = -std::log(r1) * attenuationTotalInv; // cm

            // where do we hit an object
            const auto intersection = m_kdtree.intersect(p, m_aabb);

            if (intersection.valid()) { // Do we intersect anything?
                if (intersection.intersection < stepLength) {
                    // Object is closer than free path.
                    if (!intersection.rayOriginIsInsideItem) { // if we are not already inside the object (we seldom are)
                        p.border_translate(intersection.intersection);
                    }

                    std::visit([&p, &state](auto& it) { it.transport(p, state); }, *intersection.item);
                    continueSampling = p.energy > 0;
                } else { // Free path is closer than object, we interact in the world empty space
                    p.translate(stepLength);
                    const auto interactionResult = interactions::template interact<WorldShells(), 2>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_energyScored.scoreEnergy(interactionResult.energyImparted);
                }
            } else { // We do not intersect any object
                p.translate(stepLength);
                if (basicshape::AABB::pointInside(p.pos, m_aabb)) { // Are we still inside world?
                    const auto interactionResult = interactions::template interact<WorldShells(), 2>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_energyScored.scoreEnergy(interactionResult.energyImparted);
                } else {
                    continueSampling = false;
                }
            }
        }
    }

private:
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<std::variant<F, Us...>> m_items;
    KDTreeFlat<F, Us...> m_kdtree;
    MaterialType m_fillMaterial;
    double m_fillMaterialDensity = 0.001225;
    EnergyScore m_energyScored;
    std::vector<std::string> m_item_names;
};
}

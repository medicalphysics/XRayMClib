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
#include "xraymc/serializer.hpp"
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

/**
 * @brief Concept satisfied when @p U is one of the types in the pack @p Us.
 */
template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

/**
 * @brief Container and transport engine for a heterogeneous collection of world items.
 *
 * Holds a flat list of world items as `std::variant<F, Us...>`, accelerated by a
 * KDTreeFlat for intersection queries. The space between items is filled with a
 * homogeneous material (default: dry air at 1.225×10⁻³ g/cm³). Particle transport
 * in the fill material uses analog Monte Carlo sampling with 16-shell cross-sections.
 *
 * Typical workflow:
 * 1. Add items with addItem().
 * 2. Call build() to construct the KD-tree and compute the AABB.
 * 3. Call transport() for each simulated particle.
 * 4. Call addEnergyScoredToDoseScore() and read dose from individual items.
 *
 * The world can be serialized to a byte buffer and reconstructed via deserialize().
 *
 * @tparam F     First (mandatory) world item type, satisfying WorldItemType.
 * @tparam Us... Additional world item types, each satisfying WorldItemType.
 */
template <WorldItemType F, WorldItemType... Us>
class World {
    static constexpr std::size_t WorldShells() { return 16; }
    using MaterialType = Material<WorldShells()>;

public:
    /// @brief Constructs an empty world with dry air as the fill material.
    World()
        : m_fillMaterial(MaterialType::byNistName("Air, Dry (near sea level)").value())
    {
    }

    /**
     * @brief Constructs an empty world and pre-reserves storage for @p reserveNumberOfWorldItems items.
     * @param reserveNumberOfWorldItems Hint for the initial capacity of the item vector.
     */
    World(std::size_t reserveNumberOfWorldItems)
        : m_fillMaterial(MaterialType::byNistName("Air, Dry (near sea level)").value())
    {
        m_items.reserve(reserveNumberOfWorldItems);
    }

    /**
     * @brief Sets the fill material between world items.
     * @param mat Material cross-section data.
     */
    void setMaterial(const MaterialType& mat)
    {
        m_fillMaterial = mat;
    }

    /**
     * @brief Sets both the fill material and its mass density.
     * @param mat  Material cross-section data.
     * @param dens Density in g/cm³; absolute value is used.
     */
    void setMaterial(const MaterialType& mat, double dens)
    {
        m_fillMaterial = mat;
        m_fillMaterialDensity = std::abs(dens);
    }

    /**
     * @brief Sets the mass density of the fill material.
     * @param dens Density in g/cm³; absolute value is used.
     */
    void setMaterialDensity(double dens)
    {
        m_fillMaterialDensity = std::abs(dens);
    }

    /**
     * @brief Sets the fill material from elemental weight fractions.
     * @param composition Map from atomic number to weight fraction.
     */
    void setMaterialByWeight(const std::map<std::size_t, double>& composition)
    {
        m_fillMaterial = MaterialType::byWeight(composition).value();
    }

    /**
     * @brief Sets the fill material from elemental weight fractions and a density.
     * @param composition Map from atomic number to weight fraction.
     * @param dens        Density in g/cm³; absolute value is used.
     */
    void setMaterialByWeight(const std::map<std::size_t, double>& composition, double dens)
    {
        m_fillMaterial = MaterialType::byWeight(composition).value();
        m_fillMaterialDensity = std::abs(dens);
    }

    /// @brief Returns the current fill material.
    const MaterialType& fillMaterial() const { return m_fillMaterial; }

    /// @brief Returns the fill material density in g/cm³.
    double fillMaterialDensity() const { return m_fillMaterialDensity; }

    /// @brief Pre-allocates storage for @p size items to avoid repeated reallocations.
    void reserveNumberOfItems(std::size_t size)
    {
        m_items.reserve(size);
    }

    /**
     * @brief Copies @p item into the world and assigns a default name ("Item N").
     * @tparam U One of the world item types (F or Us...).
     * @param item Item to copy.
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U& item)
    {
        m_items.push_back(item);
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }

    /**
     * @brief Copies @p item into the world with the given @p name.
     * @tparam U One of the world item types (F or Us...).
     * @param item Item to copy.
     * @param name Human-readable name for the item.
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U& item, std::string_view name)
    {
        m_items.push_back(item);
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    /**
     * @brief Moves @p item into the world and assigns a default name ("Item N").
     * @tparam U One of the world item types (F or Us...).
     * @param item Item to move.
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U&& item)
    {
        m_items.push_back(std::move(item));
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }

    /**
     * @brief Moves @p item into the world with the given @p name.
     * @tparam U One of the world item types (F or Us...).
     * @param item Item to move.
     * @param name Human-readable name for the item.
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(U&& item, std::string_view name)
    {
        m_items.push_back(std::move(item));
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    /**
     * @brief Default-constructs a new item of type @p U and adds it with a default name.
     * @tparam U One of the world item types (F or Us...).
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem()
    {
        U item;
        m_items.push_back(std::move(item));
        std::string name = "Item " + std::to_string(m_items.size());
        m_item_names.push_back(name);
        return std::get<U>(m_items.back());
    }

    /**
     * @brief Default-constructs a new item of type @p U and adds it with the given @p name.
     * @tparam U One of the world item types (F or Us...).
     * @param name Human-readable name for the item.
     * @return Reference to the stored item.
     */
    template <AnyWorldItemType<F, Us...> U>
    auto& addItem(std::string_view name)
    {
        U item;
        m_items.push_back(std::move(item));
        m_item_names.push_back(std::string(name));
        return std::get<U>(m_items.back());
    }

    /// @brief Returns the item vector (const overload).
    const auto& items() const
    {
        return m_items;
    }

    /// @brief Returns the item vector (mutable overload).
    auto& items()
    {
        return m_items;
    }

    /// @brief Returns the item name vector (const overload).
    const auto& itemNames() const
    {
        return m_item_names;
    }

    /// @brief Returns the item name vector (mutable overload).
    auto& itemNames()
    {
        return m_item_names;
    }

    /// @brief Returns a vector of mutable pointers to all stored item variants.
    std::vector<std::variant<F, Us...>*> getItemPointers()
    {
        std::vector<std::variant<F, Us...>*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return &v;
        });
        return ptrs;
    }

    /// @brief Returns a vector of const pointers to all stored item variants.
    std::vector<const std::variant<F, Us...>*> getItemPointers() const
    {
        std::vector<const std::variant<F, Us...>*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return &v;
        });
        return ptrs;
    }

    /**
     * @brief Returns a mutable pointer to the item with the given @p name, or nullptr if not found.
     * @param name Name to search for (exact match).
     */
    std::variant<F, Us...>* getItemPointerFromName(std::string_view name)
    {
        for (std::size_t i = 0; i < m_item_names.size(); ++i) {
            if (m_item_names[i].compare(name) == 0) {
                return &m_items[i];
            }
        }
        return nullptr;
    }

    /**
     * @brief Returns a const pointer to the item with the given @p name, or nullptr if not found.
     * @param name Name to search for (exact match).
     */
    const std::variant<F, Us...>* getItemPointerFromName(std::string_view name) const
    {
        for (std::size_t i = 0; i < m_item_names.size(); ++i) {
            if (m_item_names[i].compare(name) == 0) {
                return &m_items[i];
            }
        }
        return nullptr;
    }

    /// @brief Resets the energy-score accumulators in the world fill material and all items.
    void clearEnergyScored()
    {
        m_energyScored.clear();
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearEnergyScored(); }, v);
        }
    }

    /// @brief Resets the dose-score accumulators in all items.
    void clearDoseScored()
    {
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearDoseScored(); }, v);
        }
    }

    /**
     * @brief Converts each item's accumulated energy score to a dose score.
     * @param calibration_factor Optional scaling factor applied during the conversion (default 1).
     */
    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        for (auto& v : m_items) {
            std::visit([calibration_factor](auto&& arg) { arg.addEnergyScoredToDoseScore(calibration_factor); }, v);
        }
    }

    /**
     * @brief Builds the KD-tree from the current item list and computes the padded AABB.
     *
     * Must be called after all items are added and before transport() or intersect().
     * @param AABB_padding Extra margin added to each side of the item AABB in cm
     *                     (minimum 0.1 cm is always applied; default 10 cm).
     */
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

    /**
     * @brief Builds the KD-tree and expands the AABB to encompass @p aabb.
     *
     * Calls build() with default padding first, then unions the resulting AABB
     * with @p aabb so that the world always covers the supplied region.
     * @param aabb Minimum required bounding box as {xmin,ymin,zmin,xmax,ymax,zmax} in cm.
     */
    void build(const std::array<double, 6>& aabb)
    {
        build();
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = std::min(m_aabb[i], aabb[i]);
            m_aabb[i + 3] = std::max(m_aabb[i + 3], aabb[i + 3]);
        }
    }

    /// @brief Returns the padded world AABB computed by the last call to build().
    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    /// @brief Returns the center of the world AABB in cm.
    std::array<double, 3> center() const
    {
        const auto [l, r] = vectormath::splice(m_aabb);
        return vectormath::scale(0.5, vectormath::add(l, r));
    }

    /**
     * @brief Translates all items, the world AABB, and the KD-tree split planes by @p dist.
     * @param dist Displacement vector in cm.
     */
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

    /**
     * @brief Tests a particle ray against all world items via the KD-tree.
     * @param p Particle whose position and direction define the ray.
     * @return Closest transport intersection result within the world AABB.
     */
    inline auto intersect(const ParticleType auto& p)
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    /**
     * @brief Tests a particle ray for visualization intersection via the KD-tree.
     * @param p Particle whose position and direction define the ray.
     * @return Closest visualization intersection result within the world AABB.
     */
    inline auto intersectVisualization(const ParticleType auto& p) const
    {
        return m_kdtree.intersectVisualization(p, m_aabb);
    }

    /**
     * @brief Advances a particle to the world boundary if it is currently outside.
     *
     * If the particle is already inside the AABB this is a no-op and returns true.
     * If it is outside, the ray is tested against the AABB; if it hits, the particle
     * is advanced to the surface and true is returned. Returns false if the particle
     * will never enter the world.
     * @param p Particle to transport; position is updated in place.
     * @return true if the particle is inside (or has been moved inside) the world AABB.
     */
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

    /**
     * @brief Transports a particle through the world until it exits or is absorbed.
     *
     * First advances the particle to the world boundary via transportParticleToWorld().
     * Then alternates between sampling a free path in the fill material and testing
     * whether a world item is closer:
     * - If an item boundary is closer, the particle is advanced to it and the item's
     *   own transport() is called.
     * - Otherwise the particle takes a step in the fill material and analog interaction
     *   sampling is performed (LOWENERGYCORRECTION = 2, 16 shells).
     * Transport stops when the particle exits the AABB or its energy reaches zero.
     * @param p     Particle to transport; modified in place.
     * @param state Random number generator state.
     */
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
                    // most likely the eergy of the particle has changed when interacting with world items
                    updateAttenuation = true;
                } else { // Free path is closer than object, we interact in the world empty space
                    p.translate(stepLength);
                    const auto interactionResult = interactions::template interact<WorldShells(), 2>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    if (interactionResult.energyImparted > 0)
                        m_energyScored.scoreEnergy(interactionResult.energyImparted);
                }
            } else { // We do not intersect any object
                p.translate(stepLength);
                if (basicshape::AABB::pointInside(p.pos, m_aabb)) { // Are we still inside world?
                    const auto interactionResult = interactions::template interact<WorldShells(), 2>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    if (interactionResult.energyImparted > 0)
                        m_energyScored.scoreEnergy(interactionResult.energyImparted);
                } else {
                    continueSampling = false;
                }
            }
        }
    }

    /**
     * @brief Returns the 32-byte magic identifier for this World type.
     *
     * The magic ID is a fixed-width string ("World1" padded with spaces) used as a
     * header tag during serialization to identify the buffer as a World object.
     *
     * @return A 32-byte character array containing the magic identifier.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "World1";
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether a byte buffer begins with the World magic identifier.
     *
     * Reads the first 32 bytes of @p data and compares them against `magicID()`.
     * Used to validate a serialized buffer before attempting deserialization.
     *
     * @param data View of the byte buffer to inspect; must be at least 32 bytes.
     * @return True if the buffer starts with the expected magic ID, false otherwise.
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the world to a flat byte buffer.
     *
     * The buffer layout is:
     *   1. Item count (uint64).
     *   2. Each item's magic ID and serialized payload (via `Serializer::serializeItem`).
     *   3. Fill material composition weights.
     *   4. Fill material density (double).
     *   5. Item name strings.
     *
     * The resulting buffer can be reconstructed with `deserialize()`.
     *
     * @return A `std::vector<char>` containing the complete serialized world.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();

        Serializer::serialize(static_cast<std::uint64_t>(m_items.size()), buffer);
        for (const auto& item : m_items) {
            auto item_buffer = std::visit([](const auto& arg) { return arg.serialize(); }, item);
            const auto item_name = std::visit([](const auto& arg) { return arg.magicID(); }, item);
            Serializer::serializeItem(item_name, item_buffer, buffer);
        }

        Serializer::serializeMaterialWeights(m_fillMaterial.composition(), buffer);
        Serializer::serialize(m_fillMaterialDensity, buffer);
        Serializer::serialize(m_item_names, buffer);

        return buffer;
    }

    /**
     * @brief Returns the magic IDs of all world item types supported by this World.
     *
     * Collects `magicID()` from the first type `F` and each type in `Us...` via a
     * fold expression. Used during deserialization to match each stored item's tag
     * against the known item types.
     *
     * @return A vector of 32-byte arrays, one per supported item type, in the order
     *         `F, Us...`.
     */
    static std::vector<std::array<char, 32>> itemMagicIDs()
    {
        std::vector<std::array<char, 32>> names;
        names.push_back(F::magicID());
        ((names.push_back(Us::magicID())), ...); // Fold expression magic
        return names;
    }

    /**
     * @brief Reconstructs a World from a serialized byte buffer.
     *
     * Reads the item count, then deserializes each item by matching its stored magic
     * ID against the supported types (`F, Us...`) using their `validMagicID()` methods.
     * After all items are restored, the fill material and density are read, and
     * `build()` is called to rebuild the KD-tree acceleration structure.
     *
     * Returns `std::nullopt` if any item cannot be matched to a known type or if the
     * fill material cannot be reconstructed from the stored composition weights.
     *
     * @param buffer A read-only view of the byte buffer produced by `serialize()`.
     * @return An engaged `std::optional<World<F, Us...>>` on success, or
     *         `std::nullopt` on failure.
     */
    static std::optional<World<F, Us...>> deserialize(std::span<const char> buffer)
    {
        std::uint64_t n_elements;
        buffer = Serializer::deserialize(n_elements, buffer);

        World<F, Us...> world;

        std::vector<std::array<char, 32>> IDs = world.itemMagicIDs();

        for (std::size_t i = 0; i < n_elements; ++i) {
            std::vector<char> item_buffer;
            std::array<char, 32> item_name;
            buffer = Serializer::deserializeItem(item_name, item_buffer, buffer);

            std::optional<std::variant<F, Us...>> item_opt;

            // Match iten_name to type Us by Us::validMagicID
            if (F::validMagicID(item_name)) {
                item_opt = F::deserialize(item_buffer); // First type
            } else {
                if (!item_opt) {
                    ((Us::validMagicID(item_name) && (void(item_opt = Us::deserialize(item_buffer)), true)) || ...); // Fold expression for the rest of Us types
                }
            }

            if (item_opt)
                world.m_items.push_back(item_opt.value());
            else {
                return std::nullopt;
            }
        }

        std::map<std::uint8_t, double> mat_weights;
        buffer = Serializer::deserializeMaterialWeights(mat_weights, buffer);
        auto material_opt = Material<WorldShells()>::byWeight(mat_weights);
        if (material_opt) {
            world.m_fillMaterial = material_opt.value();
        } else {
            return std::nullopt;
        }

        buffer = Serializer::deserialize(world.m_fillMaterialDensity, buffer);
        buffer = Serializer::deserialize(world.m_item_names, buffer);

        world.build();
        return world;

        // return buffer;
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

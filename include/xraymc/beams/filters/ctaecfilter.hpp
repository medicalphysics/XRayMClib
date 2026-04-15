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

#include "xraymc/constants.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"

#include <array>
#include <vector>

namespace xraymc {

/**
 * @brief Automatic Exposure Control (AEC) filter for CT along a patient axis.
 *
 * Models a CT AEC system by storing a 1-D weight profile along a directed line
 * segment from `start` to `stop`. Given a 3-D world position, the filter projects
 * that position onto the axis and returns an interpolated weight via `operator()`.
 *
 * Weights are normalised so that their integral over the full axis length equals
 * the axis length, giving an expected weight of 1.0 for a uniform distribution.
 * An alternative `normalizeBetween()` normalises over an arbitrary sub-segment.
 *
 * The default-constructed filter is a flat (uniform) profile and is considered
 * "empty" (`isEmpty()` returns true). A non-trivial profile is set via `setData()`.
 *
 * Satisfies the `SerializeItemType` concept via `magicID()`, `validMagicID()`,
 * `serialize()`, and `deserialize()`.
 */
class CTAECFilter {
public:
    /**
     * @brief Default constructor — creates a flat (uniform weight = 1) profile.
     *
     * The filter axis runs from 0 to 1 in internal distance units.
     * `isEmpty()` returns true for a default-constructed filter.
     */
    CTAECFilter()
    {
        m_data.resize(2);
        m_data[0] = std::make_pair(0.0, 1.0);
        m_data[1] = std::make_pair(1.0, 1.0);
    }

    /**
     * @brief Constructs a filter with an explicit axis and weight profile.
     *
     * Delegates to `setData(start, stop, values)`.
     *
     * @param start   World-space start point of the AEC axis [cm].
     * @param stop    World-space end point of the AEC axis [cm].
     * @param values  Weight samples uniformly spaced along the axis. Must contain
     *                at least 2 elements; fewer falls back to a flat profile.
     */
    CTAECFilter(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& values)
    {
        setData(start, stop, values);
    }

    /**
     * @brief Re-normalises the weight profile over a world-space sub-segment.
     *
     * Scales all weights so that the integral over [start, stop] equals the
     * projected length of that sub-segment, giving a mean weight of 1.0 within
     * the specified range. Has no effect on a default (empty) filter.
     *
     * @param start  World-space start of the normalisation region [cm].
     * @param stop   World-space end of the normalisation region [cm].
     */
    void normalizeBetween(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        if (!isEmpty()) {
            normalize(start, stop);
        }
    }

    /// @brief Returns the number of sample points in the weight profile.
    std::size_t size() const
    {
        return m_data.size();
    }

    /**
     * @brief Returns true if the filter holds only the default flat profile.
     *
     * A filter is considered empty (flat) when it was default-constructed or
     * when `setData` was called with fewer than 2 weight samples.
     *
     * @return True if the profile has exactly 2 samples (uniform weight = 1).
     */
    bool isEmpty() const
    {
        return m_data.size() == 2;
    }

    /**
     * @brief Returns a copy of all weight values in the profile.
     * @return Vector of weight values, one per sample point along the axis.
     */
    std::vector<double> weights() const
    {
        std::vector<double> w(m_data.size());
        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), w.begin(), [](const auto& v) { return v.second; });
        return w;
    }

    /// @brief Returns the world-space start point of the AEC axis [cm].
    const std::array<double, 3>& start() const
    {
        return m_start;
    }

    /// @brief Returns the world-space end point of the AEC axis [cm].
    std::array<double, 3> stop() const
    {
        return xraymc::vectormath::add(m_start, vectormath::scale(m_dir, m_length));
    }

    /// @brief Returns the total length of the AEC axis [cm].
    double length() const
    {
        return m_length;
    }

    /**
     * @brief Sets the filter axis and weight profile.
     *
     * Computes the unit direction vector and length from @p start to @p stop.
     * The weight samples in @p data are mapped to evenly spaced positions along
     * the axis, then normalised so that the integral over the full length equals
     * the axis length (mean weight = 1.0). Falls back to a flat profile when
     * @p data has fewer than 2 elements or when start == stop.
     *
     * @param start  World-space start point of the AEC axis [cm].
     * @param stop   World-space end point of the AEC axis [cm].
     * @param data   Weight samples, uniformly distributed from start to stop.
     */
    void setData(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& data)
    {
        m_start = start;
        const auto dir = vectormath::subtract(stop, start);
        if (vectormath::length_sqr(dir) < GEOMETRIC_ERROR()) {
            m_dir = { 0, 0, 1 };
            m_length = 1;
        } else {
            m_dir = vectormath::normalized(dir);
            m_length = vectormath::length(dir);
        }

        if (data.size() < 2) {
            m_data.resize(2);
            m_data[0] = std::make_pair(0.0, 1.0);
            m_data[1] = std::make_pair(1.0, 1.0);
            return;
        }

        m_data.resize(data.size());
        const auto step = m_length / data.size();
        for (std::size_t i = 0; i < data.size(); ++i) {
            m_data[i].first = step * i;
            m_data[i].second = data[i];
        }
        normalize();
    }

    /**
     * @brief Returns the interpolated weight at a world-space position.
     *
     * Projects @p pos onto the filter axis and delegates to the scalar overload.
     *
     * @param pos  3-D world position [cm].
     * @return Interpolated weight at the projected axis coordinate.
     */
    double operator()(const std::array<double, 3>& pos) const
    {
        return this->operator()(positionToIndex(pos));
    }

    /**
     * @brief Returns the interpolated weight at an axis distance @p d [cm].
     *
     * @p d is measured from the axis start point. Values outside [0, length] are
     * handled by the underlying interpolant's extrapolation behaviour.
     *
     * @param d  Distance from the start of the axis [cm].
     * @return Interpolated weight at position @p d.
     */
    double operator()(double d) const
    {
        return interpolate(m_data, d);
    }

    /**
     * @brief Integrates the weight profile over the full axis length.
     *
     * Uses the trapezoidal rule over all stored sample points.
     *
     * @return Integral of the weight profile from axis start to axis end.
     */
    double integrate() const
    {
        return trapz(m_data);
    }

    /**
     * @brief Integrates the weight profile between two world-space points.
     *
     * Projects @p start_p and @p stop_p onto the axis, then applies the
     * trapezoidal rule over the corresponding sub-range of the profile.
     *
     * @param start_p  World-space start of the integration range [cm].
     * @param stop_p   World-space end of the integration range [cm].
     * @return Integral of the weight profile between the two projected positions.
     */
    double integrate(const std::array<double, 3>& start_p, const std::array<double, 3>& stop_p) const
    {
        const auto start = positionToIndex(start_p);
        const auto stop = positionToIndex(stop_p);
        return trapz(m_data, start, stop);
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "CTAECFilter" padded with spaces, used by the
     *         serializer to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "CTAECFilter";
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    /**
     * @brief Checks whether @p data begins with the expected magic identifier.
     * @param data  Byte span to inspect; must be at least 32 bytes.
     * @return True if the first 32 bytes match `magicID()`, false otherwise.
     */
    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    /**
     * @brief Serializes the filter to a byte buffer.
     *
     * Writes the axis start point, unit direction, length, and all
     * (distance, weight) sample pairs using the `Serializer` format.
     *
     * @return Byte buffer containing the complete filter state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_start, buffer);
        Serializer::serialize(m_dir, buffer);
        Serializer::serialize(m_length, buffer);

        std::vector<double> data;
        data.reserve(m_data.size() * 2);
        for (const auto& d : m_data) {
            data.push_back(d.first);
            data.push_back(d.second);
        }
        Serializer::serialize(data, buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `CTAECFilter` from a serialized byte buffer.
     *
     * Reads the axis geometry and weight profile written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CTAECFilter>` on success. Currently always
     *         returns a valid object; returns `std::nullopt` only if future
     *         validation logic is added.
     */
    static std::optional<CTAECFilter> deserialize(std::span<const char> buffer)
    {
        CTAECFilter item;
        buffer = Serializer::deserialize(item.m_start, buffer);
        buffer = Serializer::deserialize(item.m_dir, buffer);
        buffer = Serializer::deserialize(item.m_length, buffer);
        std::vector<double> data;
        buffer = Serializer::deserialize(data, buffer);
        std::size_t idx = 0;
        item.m_data.resize(data.size() / 2);
        for (auto& d : item.m_data) {
            d.first = data.at(idx++);
            d.second = data.at(idx++);
        }
        return item;
    }

protected:
    /**
     * @brief Projects a world-space position onto the AEC axis and returns the axis distance.
     *
     * Computes the dot product of (pos − m_start) with the unit axis direction,
     * clamped to [0, m_length].
     *
     * @param pos  3-D world position [cm].
     * @return Distance along the axis from the start point [cm], in [0, length].
     */
    double positionToIndex(const std::array<double, 3>& start) const
    {
        const auto dist_start = vectormath::subtract(start, m_start);
        return std::clamp(vectormath::dot(dist_start, m_dir), 0.0, m_length);
    }

    /**
     * @brief Normalises the weight profile over the full axis length.
     *
     * Scales all weights so that the trapezoidal integral over [0, length]
     * equals @p m_length, giving a mean weight of 1.0.
     */
    void normalize()
    {
        const auto area = integrate();
        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto k = m_length / area;
        for (auto& d : m_data)
            d.second *= k;
    }

    /**
     * @brief Normalises the weight profile over a world-space sub-segment.
     *
     * Projects @p start and @p stop onto the axis, then scales all weights so
     * that the integral over that sub-range equals the projected sub-length,
     * giving a mean weight of 1.0 within the specified region.
     *
     * @param start  World-space start of the normalisation region [cm].
     * @param stop   World-space end of the normalisation region [cm].
     */
    void normalize(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        const auto p_start = positionToIndex(start);
        const auto p_stop = positionToIndex(stop);

        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto area = integrate(start, stop);
        const auto k = std::abs(p_stop - p_start) / area;
        for (auto& d : m_data)
            d.second *= k;
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 }; ///< World-space start point of the AEC axis [cm].
    std::array<double, 3> m_dir = { 0, 0, 1 }; ///< Unit direction vector along the AEC axis.
    double m_length = 1; ///< Total length of the AEC axis [cm].
    std::vector<std::pair<double, double>> m_data; ///< Sampled (axis distance [cm], weight) pairs.
};
}
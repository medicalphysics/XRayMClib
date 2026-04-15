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

#include "xraymc/interpolation.hpp"
#include "xraymc/serializer.hpp"

#include <array>

namespace xraymc {

/**
 * @brief Angle-dependent intensity weighting filter for CT beam shaping.
 *
 * Models a bowtie (bow-tie) filter by storing a spline interpolant that maps
 * fan angle (in radians) to a relative intensity weight. The filter is
 * normalised so that its integral over the angular range equals the angular
 * span, preserving the mean intensity.
 *
 * The `ONESIDED` template parameter controls whether the input angle data is
 * treated as symmetric (absolute value taken before lookup) or as a full
 * signed angular range.
 *
 * Internally uses a 6-knot `AkimaSplineStatic` for compact, serialization-
 * friendly storage. Satisfies the `SerializeItemType` concept via `magicID()`,
 * `validMagicID()`, `serialize()`, and `deserialize()`.
 */
class BowtieFilter {
public:
    /**
     * @brief Constructs a bowtie filter from separate angle and intensity vectors.
     *
     * @tparam ONESIDED If true (default), angles are folded to their absolute value
     *                  so that a single half-profile describes both sides of the fan.
     * @param angles_r     Fan angles in radians, one per sample point.
     * @param intensity_r  Relative intensity at each corresponding angle. Negative
     *                     values are taken as absolute. Must have the same size as
     *                     @p angles_r (excess elements in the longer vector are ignored).
     */
    template <bool ONESIDED = true>
    BowtieFilter(const std::vector<double>& angles_r, const std::vector<double>& intensity_r)
    {
        setData<ONESIDED>(angles_r, intensity_r);
    }

    /**
     * @brief Constructs a bowtie filter from a vector of (angle, intensity) pairs.
     *
     * @tparam ONESIDED If true (default), angles are folded to their absolute value.
     * @param data  Vector of `{angle [rad], relative intensity}` pairs.
     */
    template <bool ONESIDED = true>
    BowtieFilter(const std::vector<std::pair<double, double>>& data)
    {
        setData<ONESIDED>(data);
    }

    /**
     * @brief Default constructor — loads a built-in profile from a Siemens Definition Flash CT.
     *
     * Provides a reasonable generic bowtie profile when no measured data is available.
     */
    BowtieFilter()
    {
        // generic filter from a Siemens Definition Flash
        static std::vector<std::pair<double, double>> data = {
            { 0.166511074, 3.53208 },
            { 0.000000000, 13.9167 },
            { 0.041992107, 12.5868 },
            { 0.083836642, 9.41943 },
            { 0.246954945, 1.96665 },
            { 0.324269441, 1.27605 },
            { 0.390607044, 0.947716 }
        };
        setData<true>(data);
    }

    /**
     * @brief Returns the filter intensity weight for the given fan angle.
     *
     * @tparam ONESIDED If true (default), the absolute value of @p angle is used,
     *                  assuming a symmetric filter profile.
     * @param angle  Fan angle in radians.
     * @return Normalised intensity weight at that angle (spline-interpolated).
     */
    template <bool ONESIDED = true>
    double operator()(double angle) const
    {
        if constexpr (ONESIDED)
            return m_inter(std::abs(angle));
        else
            return m_inter(angle);
    }

    /**
     * @brief Sets the filter profile from a vector of (angle, intensity) pairs.
     *
     * Sorts the data by angle, fits an Akima spline, then scales the spline so
     * that its integral over the angular range equals the angular span (normalisation).
     *
     * @tparam ONESIDED If true (default), all angles are converted to their absolute
     *                  value before sorting and fitting.
     * @param data  Vector of `{angle [rad], relative intensity}` pairs. Intensity
     *              values are taken as absolute.
     */
    template <bool ONESIDED = true>
    void setData(std::vector<std::pair<double, double>> data)
    {
        for (auto& d : data) {
            if constexpr (ONESIDED)
                d.first = std::abs(d.first);
            d.second = std::abs(d.second);
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        m_inter.setup(data);

        const auto start = data.front().first;
        const auto stop = data.back().first;
        const auto area = m_inter.integral(start, stop);

        m_inter.scale((stop - start) / area);
    }

    /**
     * @brief Sets the filter profile from separate angle and intensity vectors.
     *
     * Zips the two vectors (up to the shorter length) into pairs and delegates to
     * the pair-vector overload of `setData`.
     *
     * @tparam ONESIDED If true (default), angles are folded to their absolute value.
     * @param angles_r     Fan angles in radians.
     * @param intensity_r  Relative intensity at each angle; negative values are made positive.
     */
    template <bool ONESIDED = true>
    void setData(const std::vector<double>& angles_r, const std::vector<double>& intensity_r)
    {
        const auto N = std::min(angles_r.size(), intensity_r.size());
        std::vector<std::pair<double, double>> data(N);

        for (std::size_t i = 0; i < N; ++i) {
            data[i].first = angles_r[i];
            data[i].second = std::abs(intensity_r[i]);
        }

        setData<ONESIDED>(data);
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "BowtieFilter" padded with spaces, used by the serializer
     *         to identify and validate stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "BowtieFilter";
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
     * Stores the internal Akima spline knot data using the `Serializer` format.
     *
     * @return Byte buffer containing the serialized filter state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_inter.copyInternalData(), buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `BowtieFilter` from a serialized byte buffer.
     *
     * Reads the Akima spline knot data written by `serialize()` and reconstructs
     * the internal interpolant.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<BowtieFilter>` on success, or `std::nullopt` if
     *         the spline data is invalid or incomplete.
     */
    static std::optional<BowtieFilter> deserialize(std::span<const char> buffer)
    {
        BowtieFilter item;
        std::vector<double> filterdata;
        buffer = Serializer::deserialize(filterdata, buffer);
        auto inter_opt = item.m_inter.fromInternalData(filterdata);
        if (inter_opt)
            item.m_inter = inter_opt.value();
        else
            return std::nullopt;
        return std::make_optional(item);
    }

private:
    AkimaSplineStatic<double, 6> m_inter; ///< 6-knot Akima spline mapping angle → intensity weight.
};
}
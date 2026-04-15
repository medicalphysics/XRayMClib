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

#pragma once

#include "xraymc/constants.hpp"
#include "xraymc/interpolation.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/vectormath.hpp"

#include <array>
#include <numbers>
#include <vector>

namespace xraymc {

/**
 * @brief Organ-based AEC filter providing angle-dependent tube-current modulation.
 *
 * Models CT organ dose-control (ODC / organ AEC) by dividing the gantry rotation
 * into three angular zones:
 *
 * - **Low-weight zone** [start_angle, stop_angle]: the tube current is reduced by
 *   `lowWeight()` relative to the high-weight reference (e.g. anterior view over
 *   the eye lens or breast).
 * - **Ramp zone** (width `rampAngle()` on each side): a linear transition from the
 *   low weight back to the high weight, avoiding abrupt steps.
 * - **High-weight zone**: the remainder of the rotation. Weight equals `maxWeight()`,
 *   which is 1.0 when `compensateOutside` is false, or a compensated value > 1.0
 *   when `compensateOutside` is true to preserve the total integral dose.
 *
 * All angles are stored in radians in (−π, π]. Degree convenience setters/getters
 * are provided for every angle parameter.
 *
 * When `useFilter()` is false the filter is disabled; callers should check this
 * flag before applying the weight.
 *
 * Satisfies the `SerializeItemType` concept via `magicID()`, `validMagicID()`,
 * `serialize()`, and `deserialize()`.
 */
class CTOrganAECFilter {
public:
    /**
     * @brief Default constructor — creates a disabled filter with default parameters.
     *
     * Low weight: 0.6, start angle: 0 rad, stop angle: π/2 rad, ramp: 0,
     * compensate outside: true, use filter: false.
     */
    CTOrganAECFilter()
    {
    }

    /**
     * @brief Constructs and configures the organ AEC filter in one step.
     *
     * @param low_weight   Relative tube-current weight in the low zone (0.001–1.0).
     * @param start_angle  Start of the low-weight angular zone [rad].
     * @param stop_angle   End of the low-weight angular zone [rad].
     * @param ramp_angle   Half-width of the linear ramp transition [rad] (default: 0,
     *                     meaning an abrupt step). Clamped to [0, π/4].
     */
    CTOrganAECFilter(double low_weight, double start_angle, double stop_angle, double ramp_angle = 0)
    {
        setLowWeightFactor(low_weight);
        setStartStopAngles(start_angle, stop_angle);
        setRampAngle(ramp_angle);
    }

    /**
     * @brief Returns the angular weight for the given gantry angle.
     *
     * Normalises @p angle to (−π, π], then returns:
     * - `lowWeight()` in the low zone [start, stop],
     * - a linearly interpolated value in each ramp region, or
     * - `maxWeight()` in the high zone.
     *
     * @param angle  Gantry angle [rad].
     * @return Relative tube-current weight at that angle.
     */
    double operator()(double angle) const
    {
        const auto nang = normalize_angle(angle);
        if (m_start_angle < nang && nang < m_stop_angle) {
            // low region
            return m_weight_factor;
        } else if (m_start_angle - m_ramp_angle < nang && nang < m_stop_angle + m_ramp_angle) {
            // ramp region
            const auto delta = (maxWeight() - m_weight_factor) / m_ramp_angle;
            const auto dang = nang > m_stop_angle ? nang - m_stop_angle : m_start_angle - nang;
            return dang * delta + m_weight_factor;
        } else {
            // high region
            return maxWeight();
        }
    }

    /**
     * @brief Enables or disables the organ AEC filter.
     *
     * When disabled, callers should treat the weight as 1.0 for all angles.
     *
     * @param on  True to enable, false to disable.
     */
    void setUseFilter(bool on)
    {
        m_useFilter = on;
    }

    /// @brief Returns true if the organ AEC filter is active.
    bool useFilter() const { return m_useFilter; }

    /**
     * @brief Sets whether the high-weight zone compensates for the dose reduction.
     *
     * When true, `maxWeight()` is computed so that the integral of the weight
     * profile over a full 2π rotation equals 2π (mean weight = 1.0), preserving
     * the total dose despite the low-weight zone. When false, `maxWeight()` is 1.0.
     *
     * @param on  True to enable compensation (default), false to disable.
     */
    void setCompensateOutside(bool on)
    {
        m_compensate_outside = on;
    }

    /// @brief Returns true if the high-weight zone compensates for the low-weight zone.
    bool compensateOutside() const
    {
        return m_compensate_outside;
    }

    /**
     * @brief Sets the start angle of the low-weight zone [rad].
     *
     * The angle is normalised to (−π, π].
     *
     * @param ang  Start angle in radians.
     */
    void setStartAngle(double ang)
    {
        m_start_angle = normalize_angle(ang);
    }

    /**
     * @brief Sets the start angle of the low-weight zone [degrees].
     * @param ang  Start angle in degrees; converted to radians internally.
     */
    void setStartAngleDeg(double ang)
    {
        m_start_angle = normalize_angle(ang * DEG_TO_RAD<double>());
    }

    /**
     * @brief Sets the stop angle of the low-weight zone [rad].
     *
     * The angle is normalised to (−π, π].
     *
     * @param ang  Stop angle in radians.
     */
    void setStopAngle(double ang)
    {
        m_stop_angle = normalize_angle(ang);
    }

    /**
     * @brief Sets the stop angle of the low-weight zone [degrees].
     * @param ang  Stop angle in degrees; converted to radians internally.
     */
    void setStopAngleDeg(double ang)
    {
        m_stop_angle = normalize_angle(ang * DEG_TO_RAD<double>());
    }

    /// @brief Returns the stop angle of the low-weight zone [rad].
    double stopAngle() const { return m_stop_angle; }
    /// @brief Returns the start angle of the low-weight zone [rad].
    double startAngle() const { return m_start_angle; }
    /// @brief Returns the stop angle of the low-weight zone [degrees].
    double stopAngleDeg() const { return m_stop_angle * RAD_TO_DEG<double>(); }
    /// @brief Returns the start angle of the low-weight zone [degrees].
    double startAngleDeg() const { return m_start_angle * RAD_TO_DEG<double>(); }

    /**
     * @brief Sets both the start and stop angles of the low-weight zone [rad].
     * @param min  Start angle in radians.
     * @param max  Stop angle in radians.
     */
    void setStartStopAngles(double min, double max)
    {
        setStartAngle(min);
        setStopAngle(max);
    }

    /**
     * @brief Sets the half-width of the linear ramp transition [rad].
     *
     * The ramp smoothly interpolates from `lowWeight()` to `maxWeight()` over
     * this angular width on each side of the low-weight zone. Clamped to [0, π/4].
     *
     * @param ang  Ramp half-width in radians.
     */
    void setRampAngle(double ang)
    {
        m_ramp_angle = std::clamp(std::abs(ang), 0.0, std::numbers::pi_v<double> / 4.0);
    }

    /**
     * @brief Sets the ramp half-width [degrees].
     * @param ang  Ramp half-width in degrees; converted to radians internally.
     */
    void setRampAngleDeg(double ang) { setRampAngle(ang * DEG_TO_RAD<double>()); }
    /// @brief Returns the ramp half-width [rad].
    double rampAngle() const { return m_ramp_angle; }
    /// @brief Returns the ramp half-width [degrees].
    double rampAngleDeg() const { return m_ramp_angle * RAD_TO_DEG<double>(); }

    /**
     * @brief Sets the relative tube-current weight in the low-weight zone.
     *
     * Clamped to [0.001, 1.0] to avoid zero-weight (which would imply the tube
     * is completely off in that zone).
     *
     * @param w  Low-zone weight factor (0.001–1.0).
     */
    void setLowWeightFactor(double w)
    {
        m_weight_factor = std::clamp(std::abs(w), 0.001, 1.0);
    }

    /// @brief Returns the relative tube-current weight in the low-weight zone.
    double lowWeight() const { return m_weight_factor; }

    /**
     * @brief Returns the weight applied in the high-weight zone.
     *
     * When `compensateOutside()` is true, the returned value is greater than 1.0
     * and is calculated so that integrating the full angular weight profile over
     * 2π gives 2π (mean weight = 1.0). When false, returns 1.0.
     *
     * @return High-zone tube-current weight.
     */
    double maxWeight() const
    {
        if (m_compensate_outside) {
            const auto al = std::abs(m_start_angle - m_stop_angle);
            const auto ar = m_ramp_angle;
            const auto ah = 2 * std::numbers::pi_v<double> - 2 * ar - al;
            return (std::numbers::pi_v<double> * 2 - m_weight_factor * (al + ar)) / (ah + ar);
        } else {
            return 1.0;
        }
    }

    /**
     * @brief Returns the 32-byte magic identifier for this type.
     * @return Fixed-length tag "CTOrganFilter" padded with spaces, used by the
     *         serializer to identify stored data blocks.
     */
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "CTOrganFilter";
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
     * Writes start angle, stop angle, ramp angle, low-weight factor,
     * compensate-outside flag, and use-filter flag using the `Serializer` format.
     *
     * @return Byte buffer containing the complete filter state.
     */
    std::vector<char> serialize() const
    {
        auto buffer = Serializer::getEmptyBuffer();
        Serializer::serialize(m_start_angle, buffer);
        Serializer::serialize(m_stop_angle, buffer);
        Serializer::serialize(m_ramp_angle, buffer);
        Serializer::serialize(m_weight_factor, buffer);
        Serializer::serialize(static_cast<std::uint64_t>(m_compensate_outside), buffer);
        Serializer::serialize(static_cast<std::uint64_t>(m_useFilter), buffer);
        return buffer;
    }

    /**
     * @brief Reconstructs a `CTOrganAECFilter` from a serialized byte buffer.
     *
     * Reads the six fields written by `serialize()`.
     *
     * @param buffer  Byte span produced by a prior `serialize()` call.
     * @return An engaged `optional<CTOrganAECFilter>` containing the restored filter.
     */
    static std::optional<CTOrganAECFilter> deserialize(std::span<const char> buffer)
    {
        CTOrganAECFilter item;
        buffer = Serializer::deserialize(item.m_start_angle, buffer);
        buffer = Serializer::deserialize(item.m_stop_angle, buffer);
        buffer = Serializer::deserialize(item.m_ramp_angle, buffer);
        buffer = Serializer::deserialize(item.m_weight_factor, buffer);
        std::uint64_t outside, useFilter;
        buffer = Serializer::deserialize(outside, buffer);
        buffer = Serializer::deserialize(useFilter, buffer);
        item.m_compensate_outside = outside > 0;
        item.m_useFilter = useFilter > 0;
        return item;
    }

protected:
    /**
     * @brief Normalises an angle to the half-open interval (−π, π].
     *
     * Repeatedly adds or subtracts 2π until the value falls within range.
     *
     * @param angle  Input angle in radians (any value).
     * @return Equivalent angle in (−π, π].
     */
    static double normalize_angle(double angle)
    {
        while (angle < std::numbers::pi_v<double>) {
            angle += std::numbers::pi_v<double> * 2;
        }
        while (angle > std::numbers::pi_v<double>) {
            angle -= std::numbers::pi_v<double> * 2;
        }
        return angle;
    }

private:
    double m_start_angle = 0; ///< Start of the low-weight zone [rad].
    double m_stop_angle = std::numbers::pi_v<double> / 2; ///< End of the low-weight zone [rad].
    double m_ramp_angle = 0; ///< Half-width of the linear ramp transition [rad].
    double m_weight_factor = 0.6; ///< Relative tube-current weight in the low zone.
    bool m_compensate_outside = true; ///< If true, high-zone weight is boosted to preserve mean dose.
    bool m_useFilter = false; ///< If false, the filter is disabled.
};
}
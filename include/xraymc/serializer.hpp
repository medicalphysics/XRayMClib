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

Copyright 2026 Erlend Andersen
*/

#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <expected>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "xraymc/world/dosescore.hpp"

namespace xraymc {

/*
FILEDESCRIPTION the saved files contains a version header of 16 chars and binary data

Example:
"xraymc1        "


*/
/**
 * @brief Concept satisfied by types that can be round-tripped through a byte buffer.
 *
 * A conforming type must provide:
 * - `serialize()` — returns a `std::vector<char>` containing the serialized payload.
 * - `deserialize(buffer)` — static or member function returning `std::optional<U>`.
 * - `magicID()` — returns a `std::array<char, 32>` identifying the type.
 * - `validMagicID(buffer)` — returns `bool` indicating whether @p buffer starts
 *   with the expected magic tag.
 *
 * @tparam U The candidate type.
 */
template <typename U>
concept SerializeItemType = requires(U u, std::span<const char> buffer) {
    {
        u.deserialize(buffer)
    } -> std::same_as<std::optional<U>>;
    {
        u.serialize()
    } -> std::same_as<std::vector<char>>;
    {
        u.magicID()
    } -> std::same_as<std::array<char, 32>>;
    {
        u.validMagicID(buffer)
    } -> std::same_as<bool>;
};

/**
 * @brief Utility class for binary serialization and deserialization of XRayMClib objects.
 *
 * Provides a collection of static helper methods for reading and writing primitive
 * scalars, vectors, arrays, strings, material composition maps, and `DoseScore`
 * objects to/from flat `std::vector<char>` byte buffers.
 *
 * File I/O methods prepend/strip a 16-byte version header ("xraymc1        ").
 * Each named item in a buffer is preceded by a 32-byte magic tag and a uint64
 * payload length, allowing heterogeneous objects to be stored and recovered in order.
 *
 * All `deserialize` overloads consume bytes from a `std::span<const char>` and
 * return the remaining span, enabling sequential deserialization by chaining calls.
 * They throw `std::length_error` if the buffer is too short or a magic tag mismatches.
 */
class Serializer {
public:
    /**
     * @brief Error codes returned by the file `read()` method.
     */
    enum class parse_error {
        buffer_size_short,        ///< File is too small to contain a valid header.
        buffer_heading_mismatch,  ///< File header does not match the expected tag.
        buffer_version_mismatch,  ///< File version tag does not match `version()`.
        buffer_file_error         ///< File could not be opened or read.
    };

    /**
     * @brief Constructs a Serializer optionally bound to a file path.
     * @param filename Path used by the non-static `write()` overload (default: empty).
     */
    Serializer(const std::string& filename = "")
        : m_filename(filename)
    {
    }

    /**
     * @brief Returns an empty byte buffer suitable as a starting point for serialization.
     * @return An empty `std::vector<char>`.
     */
    static std::vector<char> getEmptyBuffer()
    {
        return std::vector<char> { };
    }

    /**
     * @brief Returns a 32-byte name-ID template filled with spaces.
     *
     * Used as a base when constructing magic ID arrays for named items.
     *
     * @return 32-byte array of space characters.
     */
    static constexpr std::array<char, 32> getNameIDTemplate()
    {
        std::array<char, 32> t;
        t.fill(' ');
        return t;
    }

    /**
     * @brief Returns the 16-byte file version tag written at the start of every file.
     * @return Span over the literal string "xraymc1        " (16 chars).
     */
    static std::span<const char, 16> version()
    {
        return std::span { "xraymc1        " };
    }

    /**
     * @brief Writes @p buffer to the file path set at construction time.
     *
     * Prepends the 16-byte version tag before the payload.
     *
     * @param buffer Serialized payload to write.
     * @return True on success, false if no filename was set or the file cannot be opened.
     */
    bool write(const std::vector<char>& buffer) const
    {
        if (m_filename.size() == 0)
            return false;
        else
            return write(m_filename, buffer);
    }

    /**
     * @brief Writes @p buffer to @p filename, prefixed by the version tag.
     * @param filename Destination file path.
     * @param buffer   Serialized payload to write.
     * @return True on success, false if the file cannot be opened.
     */
    static bool write(const std::string& filename, const std::vector<char>& buffer)
    {
        std::ofstream of;
        of.open(filename, std::ios::binary);
        if (of.is_open()) {
            auto ver = version();
            of.write(ver.data(), ver.size());
            of.write(buffer.data(), buffer.size());
            return true;
        }
        return false;
    }

    /**
     * @brief Reads a serialized file and returns the payload without the version header.
     *
     * Opens @p filename in binary mode, verifies the 16-byte version tag, strips it,
     * and returns the remaining bytes. Returns an unexpected value on any failure.
     *
     * @param filename Source file path.
     * @return The payload bytes on success, or a `parse_error` on failure.
     */
    static std::expected<std::vector<char>, parse_error> read(const std::string& filename)
    {
        // reading buffer
        std::ifstream buffer_file(filename, std::ios::binary);
        if (buffer_file.good()) {
            std::vector<char> data(std::istreambuf_iterator<char>(buffer_file), { });
            const auto ver = version();
            if (data.size() > ver.size()) {
                if (std::search(data.cbegin(), data.cbegin() + ver.size(), ver.cbegin(), ver.cend()) < data.cbegin() + ver.size()) {
                    data.erase(data.cbegin(), data.cbegin() + ver.size());
                    return data;
                } else {
                    return std::unexpected(parse_error::buffer_version_mismatch);
                }
            } else {
                return std::unexpected(parse_error::buffer_size_short);
            }
        }
        return std::unexpected(parse_error::buffer_file_error);
    }

    /**
     * @brief Reads the 32-byte magic tag from the front of @p buffer.
     * @param buffer Buffer whose first 32 bytes contain the item name.
     * @return The 32-byte name array.
     */
    static std::array<char, 32> getCurrentItemName(std::span<const char> buffer)
    {
        std::array<char, 32> name;
        std::copy(buffer.cbegin(), buffer.cbegin() + name.size(), name.begin());
        return name;
    }

    /**
     * @brief Appends a named item to @p buffer as [32-byte name][uint64 size][payload].
     * @param name   32-byte magic tag identifying the item type.
     * @param in     Serialized item payload.
     * @param buffer Destination buffer; the tag, size, and payload are appended.
     */
    static void serializeItem(const std::array<char, 32>& name, std::span<const char> in, std::vector<char>& buffer)
    {

        buffer.reserve(in.size() + name.size() + sizeof(std::uint64_t));
        std::copy(name.cbegin(), name.cend(), std::back_inserter(buffer));

        const std::uint64_t size = in.size();
        const auto char_size = reinterpret_cast<const char*>(&size);
        std::copy(char_size, char_size + sizeof(std::uint64_t), std::back_inserter(buffer));
        std::copy(in.cbegin(), in.cend(), std::back_inserter(buffer));
    }

    /**
     * @brief Reads a named item from the front of @p buffer.
     *
     * Extracts the 32-byte magic tag into @p name, reads the uint64 payload size,
     * copies the payload into @p out, and returns the remaining buffer span.
     * Throws `std::length_error` if the buffer is too short.
     *
     * @param name   Output: receives the 32-byte item tag.
     * @param out    Output: receives the item payload bytes.
     * @param buffer Input buffer positioned at the start of a serialized item.
     * @return Remaining buffer span after the consumed item.
     */
    static std::span<const char> deserializeItem(std::array<char, 32>& name, std::vector<char>& out, std::span<const char> buffer)
    {
        if (buffer.size() < name.size()) {
            throw std::length_error("Buffer lenght is too short to contain item requested.");
        }

        /*if (std::search(buffer.cbegin(), buffer.cbegin() + name.size(), name.cbegin(), name.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer do not contain item requested.");
        }*/
        std::copy(buffer.cbegin(), buffer.cbegin() + name.size(), name.begin());
        buffer = buffer.subspan(name.size());
        std::uint64_t size;
        buffer = deserialize(size, buffer);
        if (buffer.size() < size)
            throw std::length_error("Buffer lenght is too short to contain item requested.");
        out.resize(size);
        std::copy(buffer.cbegin(), buffer.cbegin() + size, out.begin());
        return buffer.subspan(size);
    }

    /**
     * @brief Serializes a `SerializeItemType` object and appends it to @p buffer.
     *
     * Calls `item.magicID()` and `item.serialize()`, then delegates to the
     * low-level `serializeItem(name, payload, buffer)` overload.
     *
     * @tparam U Type satisfying `SerializeItemType`.
     * @param item   Object to serialize.
     * @param buffer Destination buffer.
     */
    template <SerializeItemType U>
    static void serializeItem(const U& item, std::vector<char>& buffer)
    {
        const auto name = item.magicID();
        const auto ser = item.serialize();
        serializeItem(name, ser, buffer);
    }

    /**
     * @brief Appends a scalar value to @p buffer as raw bytes.
     *
     * Constrained to `double`, `uint64_t`, and `uint8_t`.
     *
     * @tparam T Scalar type.
     * @param in     Value to serialize.
     * @param buffer Destination buffer; sizeof(T) bytes are appended.
     */
    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value)
    static void serialize(T in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        auto in_c = reinterpret_cast<char*>(&in);
        std::copy(in_c, in_c + sizeof(T), dest);
    }

    /**
     * @brief Reads a scalar value from the front of @p begin.
     *
     * Constrained to `double`, `uint64_t`, and `uint8_t`.
     * Throws `std::length_error` if fewer than sizeof(T) bytes remain.
     *
     * @tparam T Scalar type.
     * @param value  Output: receives the deserialized value.
     * @param begin  Input span positioned at the scalar's first byte.
     * @return Remaining span after the consumed scalar.
     */
    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value)
    static std::span<const char> deserialize(T& value, std::span<const char> begin)
    {
        if (begin.size() < sizeof(T))
            throw std::length_error("Buffer lenght do not contain data requested.");

        auto val_ptr = reinterpret_cast<const T*>(begin.data());
        value = *val_ptr;
        return begin.subspan(sizeof(T));
    }

    /**
     * @brief Appends a vector of scalars to @p buffer as [uint64 count][raw elements].
     *
     * Constrained to `double`, `uint64_t`, `uint8_t`, and `uint32_t`.
     *
     * @tparam T Element type.
     * @param in     Vector to serialize.
     * @param buffer Destination buffer.
     */
    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value || std::is_same<T, std::uint32_t>::value)
    static void serialize(const std::vector<T>& in, std::vector<char>& buffer)
    {
        const std::uint64_t n_elements = in.size();
        const std::uint64_t size = n_elements * sizeof(T);
        serialize(n_elements, buffer);
        auto in_c = reinterpret_cast<const char*>(in.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }
    /**
     * @brief Appends a fixed-size array of scalars to @p buffer as [uint64 count][raw elements].
     *
     * Constrained to `double`, `uint64_t`, and `uint8_t`.
     *
     * @tparam T Element type.
     * @tparam N Array length.
     * @param in     Array to serialize.
     * @param buffer Destination buffer.
     */
    template <typename T, std::size_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value)
    static void serialize(const std::array<T, N>& in, std::vector<char>& buffer)
    {
        const std::uint64_t n_elements = in.size();
        const std::uint64_t size = n_elements * sizeof(T);
        serialize(n_elements, buffer);
        auto in_c = reinterpret_cast<const char*>(in.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }

    /**
     * @brief Reads a vector of scalars from @p begin.
     *
     * Reads the uint64 element count, then copies that many raw elements into @p out.
     * Throws `std::length_error` if the buffer is too short.
     * Constrained to `double`, `uint64_t`, `uint8_t`, and `uint32_t`.
     *
     * @tparam T Element type.
     * @param out   Output vector; cleared and populated.
     * @param begin Input span positioned at the serialized vector's first byte.
     * @return Remaining span after the consumed data.
     */
    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value || std::is_same<T, std::uint32_t>::value)
    static std::span<const char> deserialize(std::vector<T>& out, std::span<const char> begin)
    {
        std::uint64_t n_elements;
        auto data_start = deserialize(n_elements, begin);
        if (n_elements * sizeof(T) > data_start.size()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        out.clear();
        out.reserve(n_elements);
        auto start = reinterpret_cast<const T*>(data_start.data());
        std::copy(start, start + n_elements, std::back_inserter(out));
        return data_start.subspan(n_elements * sizeof(T));
    }

    /**
     * @brief Reads a fixed-size array of scalars from @p begin.
     *
     * Reads the uint64 element count and verifies it equals `N` before copying.
     * Throws `std::length_error` if the buffer is too short or the count mismatches.
     * Constrained to `double`, `uint64_t`, and `uint8_t`.
     *
     * @tparam T Element type.
     * @tparam N Expected array length.
     * @param out   Output array; overwritten in place.
     * @param begin Input span positioned at the serialized array's first byte.
     * @return Remaining span after the consumed data.
     */
    template <typename T, std::size_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value)
    static std::span<const char> deserialize(std::array<T, N>& out, std::span<const char> begin)
    {
        std::uint64_t n_elements;
        auto data_start = deserialize(n_elements, begin);
        if (n_elements * sizeof(T) > data_start.size()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        if (n_elements != out.size()) {
            throw std::length_error("Buffer data lenght do not match in_array lenght.");
        }

        auto start = reinterpret_cast<const T*>(data_start.data());
        std::copy(start, start + n_elements, out.data());
        return data_start.subspan(n_elements * sizeof(T));
    }

    /**
     * @brief Appends a vector of strings to @p buffer.
     *
     * Layout: [uint64 count] followed by [uint64 length][chars...] for each string.
     *
     * @param in     Strings to serialize.
     * @param buffer Destination buffer.
     */
    static void serialize(const std::vector<std::string>& in, std::vector<char>& buffer)
    {
        const std::uint64_t n_elements = in.size();
        serialize(n_elements, buffer);
        for (const auto& s : in) {
            const std::uint64_t len = s.size();
            serialize(len, buffer);
            std::copy(s.c_str(), s.c_str() + len, std::back_inserter(buffer));
        }
    }

    /**
     * @brief Reads a vector of strings from @p begin.
     *
     * Reads the uint64 count, then for each string reads its uint64 length and
     * the corresponding character bytes.
     *
     * @param out   Output vector; resized and populated.
     * @param begin Input span positioned at the serialized string vector's first byte.
     * @return Remaining span after the consumed data.
     */
    static std::span<const char> deserialize(std::vector<std::string>& out, std::span<const char> begin)
    {
        std::uint64_t n_elements;
        begin = deserialize(n_elements, begin);
        out.resize(n_elements);
        for (std::size_t i = 0; i < n_elements; ++i) {
            std::uint64_t len;
            begin = deserialize(len, begin);
            out[i] = std::string(begin.data(), len);
            begin = begin.subspan(len);
        }
        return begin;
    }

    /**
     * @brief Serializes a single material's elemental composition to @p buffer.
     *
     * Layout: 8-byte "Material" tag, uint64 material count (always 1),
     * uint64 element count, then [uint8 Z][double weight] pairs.
     *
     * @param map    Map of atomic number Z → mass fraction weight.
     * @param buffer Destination buffer.
     */
    static void serializeMaterialWeights(const std::map<std::uint8_t, double>& map, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };
        std::copy(mat.cbegin(), mat.cend(), std::back_inserter(buffer));
        serialize(std::uint64_t { 1 }, buffer);
        const auto size = static_cast<std::uint64_t>(map.size());
        serialize(size, buffer);
        if (size > 0) {
            for (const auto& [Z, w] : map) {
                serialize(Z, buffer);
                serialize(w, buffer);
            }
        }
        return;
    }

    /**
     * @brief Deserializes a single material's elemental composition from @p buffer.
     *
     * Verifies the "Material" tag and expects exactly one material entry.
     * Throws `std::length_error` if the buffer is malformed or contains multiple materials.
     *
     * @param out    Output map populated with Z → weight pairs.
     * @param buffer Input span positioned at the "Material" tag.
     * @return Remaining span after the consumed data.
     */
    static std::span<const char> deserializeMaterialWeights(std::map<std::uint8_t, double>& out, std::span<const char> buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };
        if (buffer.size() < 8) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        if (std::search(buffer.cbegin(), buffer.cbegin() + mat.size(), mat.cbegin(), mat.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        buffer = buffer.subspan(mat.size());

        std::uint64_t n_materials;
        buffer = deserialize(n_materials, buffer);
        if (n_materials != 1)
            throw std::length_error("Requested to deserialize ONE material but the buffer contains a material list");

        std::uint64_t size;
        buffer = deserialize(size, buffer);

        out.clear();
        for (std::uint64_t i = 0; i < size; i++) {
            std::uint8_t Z;
            double w;
            buffer = deserialize(Z, buffer);
            buffer = deserialize(w, buffer);
            out[Z] = w;
        }
        return buffer;
    }

    /**
     * @brief Serializes a list of materials' elemental compositions to @p buffer.
     *
     * Layout: 8-byte "Material" tag, uint64 material count, then for each material
     * uint64 element count followed by [uint8 Z][double weight] pairs.
     * Does nothing if @p maps is empty.
     *
     * @param maps   Vector of Z → weight maps, one per material.
     * @param buffer Destination buffer.
     */
    static void serializeMaterialWeights(const std::vector<std::map<std::uint8_t, double>>& maps, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };

        if (maps.size() == 0)
            return;
        std::copy(mat.cbegin(), mat.cend(), std::back_inserter(buffer));
        serialize(static_cast<std::uint64_t>(maps.size()), buffer);

        for (const auto& map : maps) {
            serialize(static_cast<std::uint64_t>(map.size()), buffer);
            for (const auto& [Z, w] : map) {
                serialize(Z, buffer);
                serialize(w, buffer);
            }
        }
        return;
    }

    /**
     * @brief Deserializes a list of materials' elemental compositions from @p buffer.
     *
     * Verifies the "Material" tag and reads all material entries.
     * Throws `std::length_error` if the buffer is malformed or contains zero materials.
     *
     * @param out    Output vector resized to the number of materials in the buffer.
     * @param buffer Input span positioned at the "Material" tag.
     * @return Remaining span after the consumed data.
     */
    static std::span<const char> deserializeMaterialWeights(std::vector<std::map<std::uint8_t, double>>& out, std::span<const char> buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };
        if (buffer.size() < 8) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        if (std::search(buffer.cbegin(), buffer.cbegin() + mat.size(), mat.cbegin(), mat.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        buffer = buffer.subspan(mat.size());

        std::uint64_t n_materials;
        buffer = deserialize(n_materials, buffer);
        if (n_materials == 0)
            throw std::length_error("Requested to deserialize material but the buffer contains zero materials");

        out.clear();
        out.resize(n_materials);
        for (std::uint64_t j = 0; j < n_materials; ++j) {
            std::uint64_t size;
            buffer = deserialize(size, buffer);
            for (std::uint64_t i = 0; i < size; i++) {
                std::uint8_t Z;
                double w;
                buffer = deserialize(Z, buffer);
                buffer = deserialize(w, buffer);
                out[j][Z] = w;
            }
        }
        return buffer;
    }

    /**
     * @brief Serializes a span of `DoseScore` objects to @p buffer.
     *
     * Layout: 8-byte "Dose    " tag, uint64 count, then for each score
     * [double dose][double variance][uint64 numberOfEvents].
     * Does nothing if @p in is empty.
     *
     * @param in     Span of dose scores to serialize.
     * @param buffer Destination buffer.
     */
    static void serializeDoseScore(std::span<const DoseScore> in, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> dose = { 'D', 'o', 's', 'e', ' ', ' ', ' ', ' ' };

        if (in.size() == 0)
            return;
        std::copy(dose.cbegin(), dose.cend(), std::back_inserter(buffer));

        const std::uint64_t size = in.size();
        serialize(size, buffer);

        for (const auto& d : in) {
            serialize(d.dose(), buffer);
            serialize(d.variance(), buffer);
            serialize(d.numberOfEvents(), buffer);
        }
        return;
    }

    /**
     * @brief Deserializes a vector of `DoseScore` objects from @p buffer.
     *
     * Verifies the "Dose    " tag and reads all entries.
     * Throws `std::length_error` if the buffer is malformed or empty.
     *
     * @param out    Output vector; cleared and populated.
     * @param buffer Input span positioned at the "Dose    " tag.
     * @return Remaining span after the consumed data.
     */
    static std::span<const char> deserializeDoseScore(std::vector<DoseScore>& out, std::span<const char> buffer)
    {
        constexpr std::array<char, 8> dose = { 'D', 'o', 's', 'e', ' ', ' ', ' ', ' ' };
        if (buffer.size() < 8) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        if (std::search(buffer.cbegin(), buffer.cbegin() + dose.size(), dose.cbegin(), dose.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        buffer = buffer.subspan(dose.size());

        std::uint64_t n_doses;
        buffer = deserialize(n_doses, buffer);
        if (n_doses == 0)
            throw std::length_error("Requested to deserialize doses but the buffer contains zero doses");

        out.clear();
        out.resize(n_doses);
        for (std::uint64_t i = 0; i < n_doses; ++i) {
            double d, v;
            std::uint64_t n;
            buffer = deserialize(d, buffer);
            buffer = deserialize(v, buffer);
            buffer = deserialize(n, buffer);
            out[i].set(d, v, n);
        }
        return buffer;
    }

    /**
     * @brief Deserializes exactly N `DoseScore` objects from @p buffer into a fixed-size array.
     *
     * Verifies the "Dose    " tag and that the stored count equals `N`.
     * Throws `std::length_error` if the buffer is malformed, empty, or the count mismatches.
     *
     * @tparam N Expected number of dose scores.
     * @param out    Output array of length N; overwritten in place.
     * @param buffer Input span positioned at the "Dose    " tag.
     * @return Remaining span after the consumed data.
     */
    template <std::size_t N>
    static std::span<const char> deserializeDoseScore(std::array<DoseScore, N>& out, std::span<const char> buffer)
    {
        constexpr std::array<char, 8> dose = { 'D', 'o', 's', 'e', ' ', ' ', ' ', ' ' };
        if (buffer.size() < 8) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        if (std::search(buffer.cbegin(), buffer.cbegin() + dose.size(), dose.cbegin(), dose.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        buffer = buffer.subspan(dose.size());

        std::uint64_t n_doses;
        buffer = deserialize(n_doses, buffer);
        if (n_doses == 0)
            throw std::length_error("Requested to deserialize doses but the buffer contains zero doses");

        if (n_doses != N)
            throw std::length_error("Requested to deserialize doses but the buffer contains different number of items");

        for (std::uint64_t i = 0; i < n_doses; ++i) {
            double d, v;
            std::uint64_t n;
            buffer = deserialize(d, buffer);
            buffer = deserialize(v, buffer);
            buffer = deserialize(n, buffer);
            out[i].set(d, v, n);
        }
        return buffer;
    }

    /**
     * @brief Serializes a single `DoseScore` to @p buffer.
     *
     * Writes the "Dose    " tag, count = 1, then [double dose][double variance][uint64 events].
     *
     * @param in     The dose score to serialize.
     * @param buffer Destination buffer.
     */
    static void serializeDoseScore(const DoseScore& in, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> dose = { 'D', 'o', 's', 'e', ' ', ' ', ' ', ' ' };
        std::copy(dose.cbegin(), dose.cend(), std::back_inserter(buffer));
        constexpr std::uint64_t size = 1;
        serialize(size, buffer);
        serialize(in.dose(), buffer);
        serialize(in.variance(), buffer);
        serialize(in.numberOfEvents(), buffer);
        return;
    }
    /**
     * @brief Deserializes a single `DoseScore` from @p buffer.
     *
     * Verifies the "Dose    " tag and that the stored count is exactly 1.
     * Throws `std::length_error` if the buffer is malformed or the count is not 1.
     *
     * @param out    Output dose score; overwritten.
     * @param buffer Input span positioned at the "Dose    " tag.
     * @return Remaining span after the consumed data.
     */
    static std::span<const char> deserializeDoseScore(DoseScore& out, std::span<const char> buffer)
    {
        constexpr std::array<char, 8> dose = { 'D', 'o', 's', 'e', ' ', ' ', ' ', ' ' };
        if (buffer.size() < 8) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }

        if (std::search(buffer.cbegin(), buffer.cbegin() + dose.size(), dose.cbegin(), dose.cend()) != buffer.cbegin()) {
            throw std::length_error("Buffer lenght do not contain data requested.");
        }
        buffer = buffer.subspan(dose.size());

        std::uint64_t n_doses;
        buffer = deserialize(n_doses, buffer);
        if (n_doses != 1)
            throw std::length_error("Requested to deserialize ONE doses but the buffer contains not one dose");

        double d, v;
        std::uint64_t n;
        buffer = deserialize(d, buffer);
        buffer = deserialize(v, buffer);
        buffer = deserialize(n, buffer);
        out.set(d, v, n);

        return buffer;
    }

private:
    std::string m_filename;
};

} // namespace xraymc

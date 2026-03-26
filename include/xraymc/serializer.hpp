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

class Serializer {
public:
    enum class parse_error {
        buffer_size_short,
        buffer_heading_mismatch,
        buffer_version_mismatch,
        buffer_file_error
    };

    Serializer(const std::string& filename = "")
        : m_filename(filename)
    {
    }

    static std::vector<char> getEmptyBuffer()
    {
        return std::vector<char> { };
    }

    static std::span<const char, 16> version()
    {
        return std::span { "xraymc1        " };
    }

    bool write(const std::vector<char>& buffer) const
    {
        if (m_filename.size() == 0)
            return false;
        else
            return write(m_filename, buffer);
    }

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

    static std::array<char, 32> getCurrentItemName(std::span<const char> buffer)
    {
        std::array<char, 32> name;
        std::copy(buffer.cbegin(), buffer.cbegin() + name.size(), name.begin());
        return name;
    }

    static void serializeItem(const std::array<char, 32>& name, std::span<const char> in, std::vector<char>& buffer)
    {

        buffer.reserve(in.size() + name.size() + sizeof(std::uint64_t));
        std::copy(name.cbegin(), name.cend(), std::back_inserter(buffer));

        const std::uint64_t size = in.size();
        const auto char_size = reinterpret_cast<const char*>(&size);
        std::copy(char_size, char_size + sizeof(std::uint64_t), std::back_inserter(buffer));
        std::copy(in.cbegin(), in.cend(), std::back_inserter(buffer));
    }

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

    template <SerializeItemType U>
    static void serializeItem(const U& item, std::vector<char>& buffer)
    {
        const auto name = item.magicID();
        const auto ser = item.serialize();
        serializeItem(name, ser, buffer);
    }

    // Serialize a double or uint64 or uint8_t value
    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value || std::is_same<T, std::uint8_t>::value)
    static void serialize(T in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        auto in_c = reinterpret_cast<char*>(&in);
        std::copy(in_c, in_c + sizeof(T), dest);
    }

    // Deserialize a double or uint64 or uint8_t value
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

    // Serialize a vector of doubles or uint64 values or uint8
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

    static void serializeMaterialWeights(const std::map<std::uint8_t, double>& map, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };
        const auto size = static_cast<std::uint64_t>(map.size());
        if (size == 0)
            return;
        std::copy(mat.cbegin(), mat.cend(), std::back_inserter(buffer));
        serialize(std::uint64_t { 1 }, buffer);
        serialize(size, buffer);
        for (const auto& [Z, w] : map) {
            serialize(Z, buffer);
            serialize(w, buffer);
        }
        return;
    }

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

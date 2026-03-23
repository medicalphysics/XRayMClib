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
#include <concepts>
#include <expected>
#include <fstream>
#include <iterator>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace xraymc {

/*
FILEDESCRIPTION the saved files contains a version header of 16 chars and binary data

Example:
"xraymc1        "


*/

class Serializer {
public:
    enum class parse_error {
        buffer_size_short,
        buffer_heading_mismatch,
        buffer_version_mismatch,
        buffer_file_error
    };

    Serializer() { }
    Serializer(const std::string& filename = "")
        : m_filename(filename)
    {
    }

    static std::vector<char> getEmptyBuffer()
    {
        return std::vector<char> {};
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
            std::vector<char> data(std::istreambuf_iterator<char>(buffer_file), {});
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

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(T in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        auto in_c = reinterpret_cast<char*>(&in);
        std::copy(in_c, in_c + sizeof(T), dest);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static std::span<const char> deserialize(T& value, std::span<const char> begin)
    {
        if (begin.size() < sizeof(T))
            throw std::length_error("Buffer lenght do not contain data requested.");

        auto val_ptr = reinterpret_cast<const T*>(begin.data());
        value = *val_ptr;
        return begin.subspan(sizeof(T));
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(std::span<const T> in, std::vector<char>& buffer)
    {
        const std::uint64_t n_elements = in.size();
        const std::uint64_t size = n_elements * sizeof(T);
        serialize(n_elements, buffer);
        auto in_c = reinterpret_cast<const char*>(in.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(const std::vector<T>& in, std::vector<char>& buffer)
    {
        serialize(std::span<const T> { in }, buffer);
    }
    template <typename T, std::uint64_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(const std::array<T, N>& in, std::vector<char>& buffer)
    {
        serialize(std::span<const T> { in }, buffer);
    }

    template <std::uint64_t N>
    static void serialize(const std::array<char, N>& in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        std::copy(in.cbegin(), in.cend(), dest);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
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

    template <typename T, std::uint64_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
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

    static void serializeMaterialWeights(const std::map<std::uint64_t, double>& map, std::vector<char>& buffer)
    {
        constexpr std::array<char, 8> mat = { 'M', 'a', 't', 'e', 'r', 'i', 'a', 'l' };
        const auto size = static_cast<std::uint64_t>(map.size());
        if (size == 0)
            return;
        serialize(mat, buffer);
        serialize(std::uint64_t { 1 }, buffer);
        serialize(size, buffer);
        for (const auto& [Z, w] : map) {
            serialize(Z, buffer);
            serialize(w, buffer);
        }
        return;
    }
    static std::span<const char> deserializeMaterialWeights(std::map<std::uint64_t, double>& out, std::span<const char>& buffer)
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

        std::map<std::uint64_t, double> map;
        for (std::uint64_t i = 0; i < size; i++) {
            std::uint64_t Z;
            double w;
            buffer = deserialize(Z, buffer);
            buffer = deserialize(w, buffer);
            map[Z] = w;
        }
        buffer;
    }

private:
    std::string m_filename;
};

} // namespace xraymc

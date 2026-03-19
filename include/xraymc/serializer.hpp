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
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
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
    Serializer() { }
    Serializer(const std::string& filename)
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

    static std::optional<std::vector<char>> read(const std::string& filename)
    {
        // reading buffer
        std::ifstream buffer_file(filename, std::ios::binary);
        if (buffer_file.good()) {
            std::vector<char> data(std::istreambuf_iterator<char>(buffer_file), {});
            auto ver = version();
            if (data.size() > ver.size()) {
                if (std::search(data.cbegin(), data.cbegin() + ver.size(), ver.cbegin(), ver.cend()) < data.cbegin() + ver.size()) {
                    data.erase(data.cbegin(), data.cbegin() + ver.size());
                    return data;
                }
            }
        }
        return std::nullopt;
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
    static char* deserialize(T& value, char* begin)
    {
        auto val_ptr = reinterpret_cast<T*>(begin);
        value = *val_ptr;
        return begin + sizeof(T);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(std::span<T> in, std::vector<char>& buffer)
    {
        const std::uint64_t n_elements = in.size();
        const std::uint64_t size = n_elements * sizeof(T);
        serialize(size, buffer);
        auto in_c = reinterpret_cast<const char*>(in.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(std::vector<T>& in, std::vector<char>& buffer)
    {
        serialize(std::span { in }, buffer);
    }
    template <typename T, std::uint64_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static void serialize(std::array<T, N>& in, std::vector<char>& buffer)
    {
        serialize(std::span { in }, buffer);
    }

    template <typename T>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static char* deserialize(std::vector<T>& out, char* begin)
    {
        std::uint64_t size;
        auto data_start = deserialize(size, begin);
        const auto n_elements = size / sizeof(T);

        out.clear();
        out.reserve(n_elements);
        auto start = reinterpret_cast<T*>(data_start);
        std::copy(start, start + n_elements, std::back_inserter(out));
        return data_start + n_elements * sizeof(T);
    }

    template <typename T, std::uint64_t N>
        requires(std::is_same<T, double>::value || std::is_same<T, std::uint64_t>::value)
    static char* deserialize(std::array<T, N>& out, char* begin)
    {
        std::uint64_t size;
        auto data_start = deserialize(size, begin);
        const auto n_elements = size / sizeof(T);

        if (n_elements != out.size()) {
            throw std::length_error("Buffer data lenghtr do not match in array lenght.");
        }

        out.clear();
        out.reserve(n_elements);
        auto start = reinterpret_cast<T*>(data_start);
        std::copy(start, start + n_elements, out.data());
        return data_start + n_elements * sizeof(T);
    }

private:
    std::string m_filename;
};

} // namespace xraymc

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

#include "xraymc/world/worlditems/triangulatedmesh/triangle.hpp"

#include <cstring>
#include <fstream>
#include <string>

namespace xraymc {

/**
 * @brief Reads a triangulated surface mesh from an STL file (ASCII or binary).
 *
 * The format is detected automatically by inspecting the first line of the file:
 * if it contains the keyword "solid" the file is treated as ASCII; otherwise it is
 * read as binary. Both overloads of `operator()` return a flat vector of Triangle
 * objects. Any error is recorded in an internal message string accessible via
 * message().
 *
 * Typical use:
 * @code
 * STLReader reader("mesh.stl");
 * auto triangles = reader();
 * if (triangles.empty())
 *     std::cerr << reader.message();
 * @endcode
 */
class STLReader {
public:
    /**
     * @brief Constructs a reader with a pre-set file path.
     * @param path Path to the STL file to read.
     */
    STLReader(const std::string& path)
        : m_filePath(path)
    {
    }

    /// @brief Constructs a reader with no file path set.
    STLReader() { }

    /**
     * @brief Sets the file path to use when operator()() is called without an argument.
     * @param path Path to the STL file.
     */
    void setFilePath(const std::string& path)
    {
        m_filePath = path;
    }

    /// @brief Returns the last error or status message produced by a read operation.
    const std::string& message() const { return m_error; }

    /**
     * @brief Reads the STL file set by the constructor or setFilePath().
     * @return Vector of triangles parsed from the file, or an empty vector on failure.
     */
    std::vector<Triangle> operator()()
    {
        return operator()(m_filePath);
    }

    /**
     * @brief Reads an STL file from @p path, auto-detecting ASCII vs binary format.
     *
     * Opens the file, checks whether the first line contains "solid" to select the
     * parser, then delegates to readSTLfileASCII() or readSTLfileBinary(). Sets the
     * internal file path to @p path.
     * @param path Path to the STL file.
     * @return Vector of triangles, or an empty vector if the file could not be opened.
     */
    std::vector<Triangle> operator()(const std::string& path)
    {
        m_filePath = path;
        std::vector<Triangle> data;
        // First we chech for binary or ascii file
        std::ifstream f(path);
        if (f.is_open()) {
            std::string line;
            std::getline(f, line);
            const auto pos = line.find("solid", 0);
            f.close();
            if (pos != std::string::npos) {
                data = readSTLfileASCII(path);
            } else {
                data = readSTLfileBinary(path);
            }
        } else {
            m_error = "Could not open file: " + path;
            f.close();
        }
        return data;
    }

protected:
    /**
     * @brief Reads a value of type @p S from a byte buffer using memcpy.
     * @tparam S Trivially copyable target type.
     * @param bytes Pointer to at least sizeof(S) contiguous bytes.
     * @return The value reconstructed from the raw bytes.
     */
    template <typename S>
    static S readFromBytes(const std::uint8_t* bytes)
    {
        S value;
        std::memcpy(&value, bytes, sizeof(S));
        return value;
    }

    /**
     * @brief Extracts per-triangle float data from a binary STL byte buffer.
     *
     * Each triangle record in the binary STL format is 50 bytes: 12 bytes for the
     * normal (3 floats) followed by 36 bytes for the three vertices (9 floats), plus
     * a 2-byte attribute count. This method copies all 12 floats (normal + vertices)
     * per triangle into a flat vector, skipping the attribute bytes.
     * @param buffer Raw file bytes; must start with the 84-byte header.
     * @param header Byte offset to the first triangle record (default 84).
     * @param offset Byte stride per triangle record (default 50).
     * @return Flat vector of floats: 12 values per triangle (3 normal + 9 vertex coords).
     */
    std::vector<float> readTrianglesFromBytes(const std::vector<std::uint8_t>& buffer, const std::size_t header = 84, const std::size_t offset = 50) const
    {
        const auto n_elements = ((buffer.size() - 84) / 50);
        std::vector<float> vec(n_elements * 12);

        std::size_t vec_pos = 0;
        for (std::size_t i = 0; i < n_elements; i++) {
            const auto b = &buffer[84 + 50 * i];
            std::memcpy(&vec[vec_pos], b, sizeof(float) * 12);
            vec_pos += 12;
        }
        return vec;
    }

    /**
     * @brief Reads triangles from a binary STL file.
     *
     * Binary STL layout: 80-byte header, 4-byte triangle count, then N × 50-byte
     * records (12-byte normal + 36-byte vertices + 2-byte attribute). The triangle
     * count is validated against the file size before parsing. Vertex coordinates
     * are stored as single-precision floats in the file and promoted to double on
     * output.
     * @param path Path to the binary STL file.
     * @return Vector of triangles, or an empty vector on error (message() is set).
     */
    std::vector<Triangle> readSTLfileBinary(const std::string& path)
    {
        constexpr auto MIN_STL_SIZE = 80 + 4 + 50;

        std::vector<Triangle> data;
        std::vector<std::uint8_t> buffer;

        std::ifstream f(path, std::ios::binary | std::ios::in | std::ios::ate);
        if (f.is_open()) {
            const auto filesize = f.tellg();
            f.seekg(0, std::ios::beg);
            buffer.resize(filesize);
            f.read((char*)&buffer[0], filesize);
        } else {
            m_error = "Could not open file: " + path;
        }
        f.close();

        if (buffer.size() < MIN_STL_SIZE) {
            m_error = "File do not contains any triangles";
            return data;
        }

        const std::uint32_t n_triangles = readFromBytes<std::uint32_t>(&buffer[80]);
        const bool valid_file = MIN_STL_SIZE + 50 * (n_triangles - 1) == buffer.size();
        if (!valid_file) {
            m_error = "Error reading triangles";
            return data;
        } else {
            const std::vector<float> data_buffer = readTrianglesFromBytes(buffer);

            std::vector<double> vertices;
            std::vector<double> normals;
            vertices.reserve(data_buffer.size() * 9);
            normals.reserve(data_buffer.size() * 3);
            for (std::size_t i = 0; i < data_buffer.size(); i = i + 12) {
                for (std::size_t j = i; j < i + 3; j++) {
                    normals.push_back(data_buffer[j]);
                }
                for (std::size_t j = i + 3; j < i + 12; j++) {
                    vertices.push_back(data_buffer[j]);
                }
            }

            const bool valid_n_points = vertices.size() % 9 == 0;

            if (valid_n_points) {
                const auto n_triangles = vertices.size() / 9;
                data.clear();
                data.reserve(n_triangles);
                for (std::size_t i = 0; i < vertices.size(); i += 9) {
                    data.emplace_back(&vertices[i]);
                }
            }
        }

        return data;
    }

    /**
     * @brief Splits @p text on the separator character @p sep, skipping empty tokens.
     * @param text Input string to split.
     * @param sep  Delimiter character.
     * @return Vector of non-empty substrings between occurrences of @p sep.
     */
    static std::vector<std::string> stringSplit(const std::string& text, char sep)
    {
        // this function splits a string into a vector based on a sep character
        // it will skip empty tokens
        std::vector<std::string> tokens;
        std::size_t start = 0, end = 0;
        while ((end = text.find(sep, start)) != std::string::npos) {
            if (end != start) {
                tokens.push_back(text.substr(start, end - start));
            }
            start = end + 1;
        }
        if (end != start) {
            tokens.push_back(text.substr(start));
        }
        return tokens;
    }

    /**
     * @brief Reads triangles from an ASCII STL file.
     *
     * Scans the file line by line. Lines whose first whitespace-delimited token is
     * "vertex" (case-insensitive) contribute the next three tokens as x, y, z
     * coordinates. Every nine collected coordinates are assembled into one Triangle.
     * @param path Path to the ASCII STL file.
     * @return Vector of triangles, or an empty vector on error (message() is set).
     */
    std::vector<Triangle> readSTLfileASCII(const std::string& path)
    {
        auto processLine = [](const std::string& line) -> std::optional<std::array<double, 3>> {
            auto words = stringSplit(line, ' ');
            if (words.size() < 4)
                return std::nullopt;
            const std::string cmp("vertex");
            std::transform(words[0].begin(), words[0].end(), words[0].begin(), [](unsigned char u) { return std::tolower(u); });
            if (cmp.compare(words[0]) == 0) {
                try {
                    std::array<double, 3> res {
                        std::stod(words[1]),
                        std::stod(words[2]),
                        std::stod(words[3])
                    };
                    return std::optional(res);
                } catch (std::invalid_argument) {
                    return std::nullopt;
                }
            }
            return std::nullopt;
        };

        std::vector<double> vertices;

        std::ifstream f(path, std::ios::in);
        if (f.is_open()) {
            for (std::string line; std::getline(f, line);) {
                auto res = processLine(line);
                if (res) {
                    for (auto val : *res)
                        vertices.push_back(val);
                }
            }

        } else {
            m_error = "Could not open file: " + path;
        }
        f.close();

        std::vector<Triangle> mesh;
        if (vertices.size() < 9) {
            m_error = "File do not contains any triangles";
            return mesh;
        } else {
            const auto n_triangles = vertices.size() / 9;
            mesh.clear();
            mesh.reserve(n_triangles);
            for (std::size_t i = 0; i < vertices.size(); i += 9) {
                mesh.emplace_back(&vertices[i]);
            }
        }

        return mesh;
    }

private:
    std::string m_error;
    std::string m_filePath;
};

}
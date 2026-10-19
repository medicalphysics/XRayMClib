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

#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangle.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

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
     * @brief Returns the number of triangles skipped by the last read operation.
     *
     * Triangles are skipped if they are degenerate (zero area), have non-finite
     * coordinates, or (ASCII files) belong to a facet without exactly three vertices.
     */
    std::size_t numberOfDroppedTriangles() const { return m_dropped; }

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
     * A file is read as binary if its size equals 84 + 50 × (triangle count in the
     * header), since many exporters also write "solid" in the header of binary files.
     * Otherwise the file is read as ASCII if the first line contains "solid", and as
     * binary if not. Delegates to readSTLfileASCII() or readSTLfileBinary(). Sets the
     * internal file path to @p path.
     * @param path Path to the STL file.
     * @return Vector of triangles, or an empty vector if the file could not be opened.
     */
    std::vector<Triangle> operator()(const std::string& path)
    {
        m_filePath = path;
        m_error.clear();
        m_dropped = 0;
        std::vector<Triangle> data;
        // First we chech for binary or ascii file
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (f.is_open()) {
            const auto filesize = static_cast<std::uint64_t>(f.tellg());
            f.seekg(0, std::ios::beg);
            std::array<char, 84> header {};
            f.read(header.data(), header.size());
            f.close();
            bool binary = false;
            if (filesize >= header.size()) {
                std::uint32_t n_triangles;
                std::memcpy(&n_triangles, header.data() + 80, sizeof(n_triangles));
                binary = filesize == 84 + 50 * static_cast<std::uint64_t>(n_triangles);
            }
            const std::string first_line(header.cbegin(), std::find(header.cbegin(), header.cend(), '\n'));
            if (!binary && first_line.find("solid", 0) != std::string::npos) {
                data = readSTLfileASCII(path);
            } else {
                data = readSTLfileBinary(path);
            }
        } else {
            m_error = "Could not open file: " + path;
            f.close();
        }
        return removeInvalidTriangles(data);
    }

protected:
    /**
     * @brief Appends a line to the status message.
     * @param msg Text to append.
     */
    void addMessage(const std::string& msg)
    {
        if (!m_error.empty())
            m_error += "\n";
        m_error += msg;
    }

    /**
     * @brief Returns true if a triangle has non-finite coordinates or (numerically) zero area.
     *
     * The area test is scale invariant: the edge cross product is compared to the square of
     * the longest edge.
     */
    static bool invalidTriangle(const Triangle& tri)
    {
        const auto& v = tri.vertices();
        for (const auto& p : v)
            for (const auto c : p)
                if (!std::isfinite(c))
                    return true;
        const auto e1 = vectormath::subtract(v[1], v[0]);
        const auto e2 = vectormath::subtract(v[2], v[0]);
        const auto e3 = vectormath::subtract(v[2], v[1]);
        const auto maxEdge2 = std::max({ vectormath::dot(e1, e1), vectormath::dot(e2, e2), vectormath::dot(e3, e3) });
        const auto crossLength = vectormath::length(vectormath::cross(e1, e2));
        return !(crossLength > 64 * std::numeric_limits<double>::epsilon() * maxEdge2);
    }

    /**
     * @brief Removes degenerate or non-finite triangles, recording how many were dropped.
     * @param triangles Triangles read from file.
     * @return The remaining valid triangles.
     */
    std::vector<Triangle> removeInvalidTriangles(std::vector<Triangle>& triangles)
    {
        const auto n = std::erase_if(triangles, [](const auto& tri) { return invalidTriangle(tri); });
        if (n > 0) {
            m_dropped += n;
            addMessage("Skipped " + std::to_string(n) + " degenerate or non-finite triangle(s)");
        }
        return std::move(triangles);
    }

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
     * @brief Reads triangles from a binary STL file.
     *
     * Binary STL layout: 80-byte header, 4-byte triangle count, then N × 50-byte
     * records (12-byte normal + 36-byte vertices + 2-byte attribute). A file shorter
     * than the triangle count requires is rejected as truncated; bytes trailing the last
     * record are ignored (noted in message()). Vertex coordinates are stored as
     * single-precision floats in the file and promoted to double on output; the stored
     * normals are ignored since normals are computed from the vertex order.
     * @param path Path to the binary STL file.
     * @return Vector of triangles, or an empty vector on error (message() is set).
     */
    std::vector<Triangle> readSTLfileBinary(const std::string& path)
    {
        constexpr std::size_t HEADER_SIZE = 80 + 4;
        constexpr std::size_t RECORD_SIZE = 50;

        std::vector<Triangle> data;
        std::vector<std::uint8_t> buffer;

        std::ifstream f(path, std::ios::binary | std::ios::in | std::ios::ate);
        if (f.is_open()) {
            const auto filesize = f.tellg();
            f.seekg(0, std::ios::beg);
            buffer.resize(filesize);
            f.read((char*)buffer.data(), filesize);
        } else {
            m_error = "Could not open file: " + path;
            return data;
        }
        f.close();

        if (buffer.size() < HEADER_SIZE) {
            m_error = "File is too small to be a binary STL file: " + path;
            return data;
        }

        const std::uint32_t n_triangles = readFromBytes<std::uint32_t>(&buffer[80]);
        if (n_triangles == 0) {
            m_error = "File do not contains any triangles";
            return data;
        }
        const std::uint64_t required = HEADER_SIZE + RECORD_SIZE * static_cast<std::uint64_t>(n_triangles);
        if (buffer.size() < required) {
            m_error = "File is truncated: header states " + std::to_string(n_triangles) + " triangles, but the file only holds " + std::to_string((buffer.size() - HEADER_SIZE) / RECORD_SIZE);
            return data;
        }
        if (buffer.size() > required)
            addMessage("Ignored " + std::to_string(buffer.size() - required) + " bytes after the last triangle");

        data.reserve(n_triangles);
        for (std::size_t i = 0; i < n_triangles; ++i) {
            // skip the normal vector (3 floats)
            const auto record = &buffer[HEADER_SIZE + RECORD_SIZE * i + 3 * sizeof(float)];
            std::array<float, 9> v;
            std::memcpy(v.data(), record, sizeof(float) * v.size());
            std::array<double, 9> vd;
            std::copy(v.cbegin(), v.cend(), vd.begin());
            data.emplace_back(vd.data());
        }
        return data;
    }

    /**
     * @brief Reads triangles from an ASCII STL file.
     *
     * The file is read as whitespace separated tokens (spaces, tabs and line endings are
     * all accepted). Each `facet ... endfacet` block must contain exactly three `vertex x y z`
     * entries; facets with missing, extra or unparsable vertices are skipped as a whole (and
     * counted in numberOfDroppedTriangles()) so they do not shift the remaining triangles.
     * Keywords are case-insensitive. Files without facet keywords are read as consecutive
     * groups of three vertices.
     * @param path Path to the ASCII STL file.
     * @return Vector of triangles, or an empty vector on error (message() is set).
     */
    std::vector<Triangle> readSTLfileASCII(const std::string& path)
    {
        std::vector<Triangle> mesh;
        std::ifstream f(path, std::ios::in);
        if (!f.is_open()) {
            m_error = "Could not open file: " + path;
            return mesh;
        }

        const auto parseDouble = [](const std::string& token, double& val) {
            const auto end = token.data() + token.size();
            auto [ptr, ec] = std::from_chars(token.data(), end, val);
            return ec == std::errc() && ptr == end;
        };

        std::array<std::array<double, 3>, 3> facet;
        std::size_t nVertices = 0;
        bool facetValid = true;
        bool inFacet = false;
        std::size_t badFacets = 0;

        const auto finishFacet = [&]() {
            if (nVertices == 3 && facetValid)
                mesh.emplace_back(facet);
            else if (nVertices > 0 || inFacet)
                ++badFacets;
            nVertices = 0;
            facetValid = true;
        };

        std::string token;
        while (f >> token) {
            std::transform(token.begin(), token.end(), token.begin(), [](unsigned char u) { return std::tolower(u); });
            if (token == "facet") {
                if (inFacet)
                    finishFacet(); // previous facet lacks endfacet
                inFacet = true;
            } else if (token == "endfacet") {
                finishFacet();
                inFacet = false;
            } else if (token == "vertex") {
                std::array<double, 3> v;
                bool ok = true;
                for (auto& c : v) {
                    std::string num;
                    ok = ok && (f >> num) && parseDouble(num, c);
                }
                if (ok && nVertices < 3)
                    facet[nVertices] = v;
                facetValid = facetValid && ok;
                ++nVertices;
                if (!inFacet && nVertices == 3)
                    finishFacet();
            }
        }
        if (inFacet || nVertices > 0)
            finishFacet();
        f.close();

        if (badFacets > 0) {
            m_dropped += badFacets;
            addMessage("Skipped " + std::to_string(badFacets) + " facet(s) without exactly three valid vertices");
        }
        if (mesh.empty())
            addMessage("File do not contains any triangles");
        return mesh;
    }

private:
    std::string m_error;
    std::size_t m_dropped = 0;
    std::string m_filePath;
};

}
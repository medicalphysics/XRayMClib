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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include "xraymc/vectormath.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace xraymc {

/**
 * @brief Reader for tetrahedral mesh files used with ICRP and IETR computational phantoms.
 *
 * Parses TetGen-style node and element files and one of several supported material/organ
 * file formats (ICRP 145 separate material + organ files, ICRP pregnant-woman combined
 * format, or IETR combined CSV format) into a TetrahedalMeshData struct. The resulting
 * data can be passed directly to TetrahedalMesh::setData(). After reading, call valid()
 * to confirm that all collections were matched with density and material data.
 */
class TetrahedalMeshReader {
public:
    /// @brief Default constructor. Creates an empty reader; call readData() before use.
    TetrahedalMeshReader() { }

    /**
     * @brief Reads node and element files, assigning default water material to all collections.
     * @param nodeFile    Path to the TetGen .node file (vertex positions).
     * @param elementFile Path to the TetGen .ele file (tetrahedron connectivity).
     */
    TetrahedalMeshReader(const std::string& nodeFile, const std::string& elementFile)
    {
        readData(nodeFile, elementFile);
    }

    /**
     * @brief Reads node and element files with separate ICRP 145 material and organ files.
     * @param nodeFile      Path to the TetGen .node file.
     * @param elementFile   Path to the TetGen .ele file.
     * @param matfilePath   Path to the ICRP 145 material definition file.
     * @param organFilePath Path to the ICRP 145 organ index file.
     */
    TetrahedalMeshReader(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matfilePath,
        const std::string& organFilePath)
    {
        readData(nodeFile, elementFile, matfilePath, organFilePath);
    }

    /**
     * @brief Reads node and element files with a combined material+organ file.
     *
     * @param nodeFile         Path to the TetGen .node file.
     * @param elementFile      Path to the TetGen .ele file.
     * @param matorganfilePath Path to the combined material/organ file.
     * @param IETR_format      When true, parses @p matorganfilePath as IETR CSV format;
     *                         otherwise as ICRP pregnant-woman combined format. Default: false.
     */
    TetrahedalMeshReader(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matorganfilePath, bool IETR_format=false)
    {
        readData(nodeFile, elementFile, matorganfilePath, IETR_format);
    }

    /**
     * @brief Reads node and element files, assigning default water material to all collections.
     *
     * Each unique collection index found in the element file is assigned a water material
     * (H 11.2%, O 88.8%) at 1 g/cm³.
     * @param nodeFile    Path to the TetGen .node file.
     * @param elementFile Path to the TetGen .ele file.
     */
    void readData(const std::string& nodeFile, const std::string& elementFile)
    {
        if (!readGeometry(nodeFile, elementFile))
            return;

        // generating collections
        auto coll = m_data.collectionIndices;
        std::sort(std::execution::par_unseq, coll.begin(), coll.end());
        coll.erase(std::unique(coll.begin(), coll.end()), coll.end());
        std::vector<ICRPCollection> collection(coll.size());
        for (std::uint32_t i = 0; i < coll.size(); ++i) {
            auto& c = collection[i];
            c.density = 1;
            c.index = coll[i];
            c.material_weights[1] = 11.2;
            c.material_weights[8] = 88.8;
            c.name = "Default water " + std::to_string(i + 1);
        }

        finalize(collection);
    }

    /**
     * @brief Reads node and element files with a combined material+organ file.
     * @param nodeFile         Path to the TetGen .node file.
     * @param elementFile      Path to the TetGen .ele file.
     * @param matorganfilePath Path to the combined material/organ file.
     * @param ietr_format      When true, parses the file as an IETR CSV format;
     *                         when false (default), uses ICRP pregnant-woman format.
     */
    void readData(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matorganfilePath, bool ietr_format = false)
    {
        auto coll = ietr_format ? readIETROrganAndMaterial(matorganfilePath) : readICRPPregnantOrganAndMaterial(matorganfilePath);
        if (readGeometry(nodeFile, elementFile))
            finalize(coll);
    }

    /**
     * @brief Reads node and element files with separate ICRP 145 material and organ files.
     * @param nodeFile      Path to the TetGen .node file.
     * @param elementFile   Path to the TetGen .ele file.
     * @param matfilePath   Path to the ICRP 145 material definition file.
     * @param organFilePath Path to the ICRP 145 organ index file.
     */
    void readData(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matfilePath,
        const std::string& organFilePath)
    {
        auto coll = readICRP145PhantomMaterialAndOrgans(matfilePath, organFilePath);
        if (readGeometry(nodeFile, elementFile))
            finalize(coll);
    }

    /// @brief Returns a mutable reference to the parsed mesh data.
    TetrahedalMeshData& data()
    {
        return m_data;
    }

    /**
     * @brief Returns true if the last readData() call succeeded and the mesh has nodes and elements.
     *
     * A false result indicates that the node/element files were missing, empty or malformed
     * (e.g. missing or duplicated node/element indices, element counts not matching the
     * file header), that the mesh is geometrically invalid (see
     * TetrahedalMeshData::validationError()), or that one or more collections could not be
     * matched with density or material data. message() describes the problem.
     */
    bool valid() const
    {
        return m_valid && m_data.nodes.size() > 0 && m_data.elements.size() > 0;
    }

    /**
     * @brief Returns the error or status message from the last readData() call.
     *
     * Empty if reading succeeded without remarks. On success it may note corrections that
     * were applied, such as inverted tetrahedra that were reoriented.
     */
    const std::string& message() const
    {
        return m_message;
    }

    /// @brief Resets all parsed data, leaving the reader in the same state as after default construction.
    void clear()
    {
        m_data = TetrahedalMeshData { };
        m_valid = false;
        m_message.clear();
        m_nodeIndexBase = 0;
    }

    /**
     * @brief Rotates all vertex node positions around @p axis by @p angle radians about the origin.
     * @param axis  Rotation axis (need not be normalized).
     * @param angle Rotation angle in radians.
     */
    void rotate(const std::array<double, 3>& axis, double angle)
    {
        std::transform(std::execution::par_unseq, m_data.nodes.cbegin(), m_data.nodes.cend(), m_data.nodes.begin(), [angle, &axis](const auto& v) {
            return vectormath::rotate(v, axis, angle);
        });
    }

    /**
     * @brief Uniformly scales all vertex node positions by factor @p s.
     * @param s Scale factor applied to every coordinate of every node.
     */
    void scale(double s)
    {
        std::transform(std::execution::par_unseq, m_data.nodes.cbegin(), m_data.nodes.cend(), m_data.nodes.begin(), [s](const auto& v) {
            return vectormath::scale(v, s);
        });
    }

protected:
    /**
     * @brief Intermediate representation of one organ/material collection during parsing.
     *
     * Holds all the per-collection data extracted from the various file formats before
     * it is merged into TetrahedalMeshData. The @p index field matches the collection
     * index stored in the element file; @p material_index is used internally to cross-
     * reference ICRP 145 material entries with organ entries.
     */
    struct ICRPCollection {
        std::uint32_t index = 1; ///< Collection index as read from the element file.
        std::uint32_t material_index = 0; ///< Internal material cross-reference index (ICRP 145).
        double density = -1; ///< Mass density in g/cm³; -1 means not yet assigned.
        std::map<std::uint8_t, double> material_weights; ///< Elemental weight fractions keyed by atomic number.
        std::string name; ///< Human-readable collection name.
    };

    /**
     * @brief Reads ICRP 145 material and organ files and merges them into a collection list.
     * @param matfilePath   Path to the ICRP 145 material file.
     * @param organFilePath Path to the ICRP 145 organ file.
     * @return Merged collection list with density, composition, and name per organ.
     */
    std::vector<ICRPCollection> readICRP145PhantomMaterialAndOrgans(
        const std::string& matfilePath,
        const std::string& organFilePath)
    {
        auto orgs = readICRP145PhantomOrgans(organFilePath);
        auto mats = readICRP145PhantomMaterials(matfilePath);
        return mergeOrganMaterialData(orgs, mats);
    }
    /**
     * @brief Combines organ index entries with material property entries.
     *
     * For each organ, finds the matching material by material_index and fills in any
     * missing density, weight fractions, or name from the material entry.
     * @param organs    Organ list (from the organ index file).
     * @param materials Material list (from the material definition file).
     * @return Merged collection list ready for mergeTetAndCollData().
     */
    static std::vector<ICRPCollection> mergeOrganMaterialData(const std::vector<ICRPCollection>& organs, const std::vector<ICRPCollection>& materials)
    {
        std::vector<ICRPCollection> res;
        res.reserve(organs.size());

        for (const auto& o : organs) {
            auto organ = o;
            for (const auto& m : materials) {
                if (organ.material_index == m.material_index) {
                    if (organ.density <= 0)
                        organ.density = m.density;
                    if (o.material_weights.size() == 0)
                        organ.material_weights = m.material_weights;
                    if (o.name.size() == 0)
                        organ.name = m.name;
                    break;
                }
            }
            res.push_back(organ);
        }
        return res;
    }

    /**
     * @brief Matches the collection list against element collection indices and populates
     *        the density, material composition, and name arrays in m_data.
     *
     * Remaps the raw collection indices from the element file to a contiguous range
     * starting at 0, then copies density, weight fractions, and names from @p collection.
     * @param collection Parsed collection list from one of the readXxx() helpers.
     * @return true if every remapped collection index was matched; false if any density
     *         remains at the sentinel value -1 (unmatched collection).
     */
    bool mergeTetAndCollData(const std::vector<ICRPCollection>& collection)
    {
        // matching collection with collection indices;
        // unique collection from elements

        const auto lut = [](const std::vector<std::uint32_t>& collIdx) {
            auto coll = collIdx;
            std::sort(std::execution::par_unseq, coll.begin(), coll.end());
            coll.erase(std::unique(coll.begin(), coll.end()), coll.end());
            std::map<std::uint32_t, std::uint32_t> lut;
            for (std::uint32_t i = 0; i < coll.size(); ++i) {
                lut[coll[i]] = i;
            }
            return lut;
        }(m_data.collectionIndices);

        std::transform(std::execution::par_unseq, m_data.collectionIndices.cbegin(), m_data.collectionIndices.cend(), m_data.collectionIndices.begin(), [&lut](auto idx) { return lut.at(idx); });

        m_data.collectionDensities.resize(lut.size());
        std::fill(m_data.collectionDensities.begin(), m_data.collectionDensities.end(), -1.0); // for completeness check
        m_data.collectionMaterialComposition.resize(lut.size());
        m_data.collectionNames.resize(lut.size());

        for (const auto& c : collection) {
            if (lut.contains(c.index)) {
                auto index = lut.at(c.index);
                m_data.collectionDensities[index] = c.density;
                m_data.collectionMaterialComposition[index] = c.material_weights;
                m_data.collectionNames[index] = c.name;
            }
        }
        // test if we are missing some data by checking density = -1
        std::string missing;
        std::size_t nMissing = 0;
        for (const auto& [fileIndex, index] : lut) {
            if (m_data.collectionDensities[index] < 0.0) {
                if (nMissing < 10)
                    missing += (nMissing == 0 ? "" : ", ") + std::to_string(fileIndex);
                ++nMissing;
            }
        }
        if (nMissing > 0) {
            m_message = "No density/material data for " + std::to_string(nMissing) + " collection(s) used by elements: " + missing + (nMissing > 10 ? ", ..." : "");
            return false;
        }
        return true;
    }

    /**
     * @brief Clears previous data and reads the node and element files.
     * @return true if both files were read without errors; message() describes any error.
     */
    bool readGeometry(const std::string& nodeFile, const std::string& elementFile)
    {
        clear();
        return readNodes(nodeFile) && readElements(elementFile);
    }

    /**
     * @brief Merges collection data into the mesh, reorients inverted tetrahedra and
     *        validates the result, setting m_valid and m_message.
     * @param collection Parsed collection list from one of the readXxx() helpers.
     */
    void finalize(const std::vector<ICRPCollection>& collection)
    {
        m_valid = mergeTetAndCollData(collection);
        if (!m_valid)
            return;
        const auto nFlipped = m_data.fixTetrahedronOrientation();
        const auto error = m_data.validationError();
        m_valid = error.empty();
        if (!m_valid)
            m_message = error;
        else if (nFlipped > 0)
            m_message = "Reoriented " + std::to_string(nFlipped) + " inverted element(s)";
    }

    /**
     * @brief Parses a single field from a character buffer, advancing past @p sep delimiters.
     *
     * Handles arithmetic types via std::from_chars and std::string by scanning to the
     * next separator (trimming trailing whitespace for strings).
     * @tparam U     Field type; must be arithmetic or std::string.
     * @param start  Pointer to the start of the remaining input.
     * @param end    Pointer one past the end of the input.
     * @param sep    Delimiter character separating fields.
     * @param val    Output variable that receives the parsed value.
     * @return Pointer to the position immediately after the parsed field.
     */
    template <typename U>
    static const char* parseLine(const char* start, const char* end, char sep, U& val)
    {
        while ((std::isspace(*start) || *start == sep) && start != end)
            ++start;

        if constexpr (std::is_arithmetic<U>::value) {
            while (*start != sep && start != end) {
                auto [ptr, ec] = std::from_chars(start, end, val);
                if (ec == std::errc())
                    return ptr;
                else if (ec == std::errc::invalid_argument)
                    ++start;
                else if (ec == std::errc::result_out_of_range)
                    return ptr;
            }
            return start;
        } else if constexpr (std::is_same<U, std::string>::value) {
            auto wstop = start;
            while (*wstop != sep && wstop != end) {
                ++wstop;
            }
            // trimming spaces from back
            if (std::distance(start, wstop) > 1) {
                while (std::isspace(*(wstop - 1)) && wstop - 1 != start) {
                    --wstop;
                }
            }
            val = std::string(start, wstop);
            return wstop;
        }
    }

    /**
     * @brief Variadic overload of parseLine() that parses multiple fields in sequence.
     * @tparam U    Type of the first field.
     * @tparam Args Types of the remaining fields.
     * @param start Pointer to the start of the remaining input.
     * @param end   Pointer one past the end of the input.
     * @param sep   Delimiter character separating fields.
     * @param val   Output variable for the first field.
     * @param args  Output variables for the remaining fields.
     * @return Pointer to the position after all fields have been parsed.
     */
    template <typename U, typename... Args>
    static const char* parseLine(const char* start, const char* end, char sep, U& val, Args&... args)
    {
        if (start == end)
            return end;

        auto next_start = parseLine(start, end, sep, val);
        if (next_start != end)
            return parseLine(next_start, end, sep, args...);
        return end;
    }

    /**
     * @brief Parses an IETR combined organ+material CSV file.
     *
     * Expects a header line followed by one row per collection. Each row contains
     * collection index, name, density, and weight fractions for H, C, N, O, Na, Mg,
     * P, S, Cl, K, Ca, Fe, I (in that order). Zero-weight elements are discarded.
     * @param matfilePath Path to the IETR CSV file.
     * @return Collection list; empty if the file could not be read or has fewer than 2 lines.
     */
    static std::vector<ICRPCollection> readIETROrganAndMaterial(const std::string& matfilePath)
    {
        std::vector<ICRPCollection> res;
        const auto data = readBufferFromFile(matfilePath);
        if (data.size() == 0)
            return res;

        // find file lines
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        std::size_t start = 0;
        std::size_t stop = 0;
        while (stop != data.size()) {
            const auto c = data[stop];
            if (c == '\n') {
                lineIdx.push_back(std::make_pair(start, stop));
                start = stop + 1;
            }
            ++stop;
        }
        if (start != data.size()) {
            lineIdx.push_back(std::make_pair(start, data.size()));
        }

        // Assume first line is header
        if (lineIdx.size() < 2)
            return res;
        res.resize(lineIdx.size() - 1);
        for (std::size_t i = 1; i < lineIdx.size(); ++i) {
            auto& organ = res[i - 1];
            parseLine(
                data.data() + lineIdx[i].first, data.data() + lineIdx[i].second, ',', organ.index, organ.name, organ.density,
                organ.material_weights[1],
                organ.material_weights[6],
                organ.material_weights[7],
                organ.material_weights[8],
                organ.material_weights[11],
                organ.material_weights[12],
                organ.material_weights[15],
                organ.material_weights[16],
                organ.material_weights[17],
                organ.material_weights[19],
                organ.material_weights[20],
                organ.material_weights[26],
                organ.material_weights[53]);
        }

        // cleanup unneeded atoms
        for (auto& c : res) {
            std::erase_if(c.material_weights, [](const auto& item) {
                const auto& [key, value] = item;
                return value <= 0;
            });
        }

        return res;
    }
    /**
     * @brief Parses an ICRP 145 phantom material definition file.
     *
     * Data begins at line 3 (0-indexed). Each line contains a material index, a
     * name, weight fractions for Z = {1,6,7,8,11,12,15,16,17,19,20,26,53}, and a
     * density value. Returns one ICRPCollection per material entry.
     * @param matfilePath Path to the ICRP 145 material file.
     * @return Material list; empty if the file could not be read.
     */
    static std::vector<ICRPCollection> readICRP145PhantomMaterials(const std::string& matfilePath)
    {
        std::vector<ICRPCollection> res;
        const auto data = readBufferFromFile(matfilePath);
        if (data.size() == 0)
            return res;

        // find file lines
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        std::size_t start = 0;
        std::size_t stop = 0;
        while (stop != data.size()) {
            const auto c = data[stop];
            if (c == '\n') {
                lineIdx.push_back(std::make_pair(start, stop));
                start = stop + 1;
            }
            ++stop;
        }
        if (start != data.size()) {
            lineIdx.push_back(std::make_pair(start, data.size()));
        }

        if (lineIdx.size() < 5)
            return res;

        // data start at line 3
        for (std::size_t i = 3; i < lineIdx.size(); ++i) {

            auto start = lineIdx[i].first;
            const auto end = lineIdx[i].second;

            auto d = data.data();
            std::uint16_t organIndex = 0;
            auto [ptr, ec] = std::from_chars(d + start, d + end, organIndex);
            if (ec == std::errc())
                start = std::distance(d, ptr);
            else if (ec == std::errc::invalid_argument)
                ++start;
            else if (ec == std::errc::result_out_of_range)
                start = std::distance(d, ptr);

            // reading to a alpha
            while (start != end && !std::isalpha(*(d + start))) {
                ++start;
            }
            const auto end_word = data.find("  ", start);

            std::string organName;
            if (end_word != std::string::npos) {
                organName = data.substr(start, end_word - start);
                start = end_word + 1;
            }

            std::array<std::size_t, 13> Z = { 1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20, 26, 53 };
            std::size_t zIdx = 0;
            std::map<std::uint8_t, double> frac;
            while (start < end) {
                double w = -1;
                auto [ptr, ec] = std::from_chars(d + start, d + end, w);
                if (ec == std::errc())
                    start = std::distance(d, ptr);
                else if (ec == std::errc::invalid_argument)
                    ++start;
                else if (ec == std::errc::result_out_of_range)
                    start = std::distance(d, ptr);

                if (w > 0) {
                    if (zIdx < Z.size()) {
                        frac[Z[zIdx]] = w;
                    } else {
                        frac[0] = w;
                    }
                }
                if (w > -1) {
                    ++zIdx;
                }
            }

            double organDensity = -1;
            if (frac.contains(0)) {
                organDensity = frac.at(0);
                frac.erase(0);
            }

            ICRPCollection organ;
            organ.material_weights = frac;
            organ.material_index = organIndex;
            organ.name = organName;
            if (organDensity > 0)
                organ.density = organDensity;
            res.push_back(organ);
        }
        return res;
    }

    /**
     * @brief Parses an ICRP 145 phantom organ index file.
     *
     * The first line is treated as a header and skipped. Each subsequent CSV line
     * contains organ index, name, material_index, and density.
     * @param organfilePath Path to the ICRP 145 organ file.
     * @return Organ list; empty if the file could not be read.
     */
    static std::vector<ICRPCollection> readICRP145PhantomOrgans(const std::string& organfilePath)
    {
        std::vector<ICRPCollection> res;
        const std::string data = readBufferFromFile(organfilePath);
        if (data.size() == 0)
            return res;

        // find file lines
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        std::size_t start = 0;
        std::size_t stop = 0;
        while (stop != data.size()) {
            const auto c = data[stop];
            if (c == '\n') {
                lineIdx.push_back(std::make_pair(start, stop));
                start = stop + 1;
            }
            ++stop;
        }
        if (start != data.size()) {
            lineIdx.push_back(std::make_pair(start, data.size()));
        }
        res.resize(lineIdx.size() - 1);
        for (std::size_t i = 1; i < lineIdx.size(); ++i) {
            auto& organ = res[i - 1];
            parseLine(data.data() + lineIdx[i].first, data.data() + lineIdx[i].second, ',', organ.index, organ.name, organ.material_index, organ.density);
        }
        return res;
    }

    /**
     * @brief Parses an ICRP pregnant-woman phantom combined material+organ file.
     *
     * The file is segmented by 'C' marker lines. Each segment encodes one collection:
     * the first line (after 'C') gives name and density; subsequent lines give the
     * collection index (prefixed 'm') and elemental weight fractions as Z*1000 + w
     * pairs.
     * @param matfile Path to the combined ICRP pregnant-woman file.
     * @return Collection list parsed from all segments in the file.
     */
    static std::vector<ICRPCollection> readICRPPregnantOrganAndMaterial(const std::string& matfile)
    {
        std::vector<ICRPCollection> organs;
        const std::string data = readBufferFromFile(matfile);

        std::vector<std::pair<std::size_t, std::size_t>> segments;
        bool alternator = true;
        for (std::size_t teller = 0; teller < data.size(); teller++) {
            if (data[teller] == 'C') {
                if (alternator) {
                    segments.push_back(std::make_pair(teller, teller));
                    alternator = false;
                } else {
                    if (data[teller - 1] == '\n') {
                        segments.back().second = teller;
                        alternator = true;
                    }
                }
            }
        }
        organs.resize(segments.size());
        std::transform(segments.cbegin(), segments.cend(), organs.begin(), [&data](const auto& seg) {
            ICRPCollection organ;
            auto begin = data.data() + seg.first + 1;
            const auto end = data.data() + seg.second;

            std::string units_number;
            begin = parseLine(begin, end, ' ', organ.name, organ.density);
            begin = std::find(begin, end, '\n');
            if (begin == end) {
                return organ; // this will result in error
            } else {
                begin++;
            }
            begin = parseLine(begin, end, 'm', organ.index);

            while (begin < end) {
                std::size_t Z = 0;
                double w = 0;
                begin = parseLine(begin, end, ' ', Z, w) + 1;
                if (begin < end) {
                    const auto Za = static_cast<std::uint8_t>(Z / 1000);
                    organ.material_weights[Za] = std::abs(w);
                }
            }
            return organ;
        });

        return organs;
    }

    /**
     * @brief Splits a buffer into lines, returning [begin, end) offsets for each line.
     * @param data      Buffer to split.
     * @param collength Estimated average line length used for pre-allocating the result.
     */
    static std::vector<std::pair<std::size_t, std::size_t>> splitLines(const std::string& data, std::size_t collength = 80)
    {
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        lineIdx.reserve(data.size() / std::max(collength, std::size_t { 1 }) + 1);
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                lineIdx.push_back({ start, i });
                start = i + 1;
            }
        }
        // last line without line ending
        if (start < data.size())
            lineIdx.push_back({ start, data.size() });
        return lineIdx;
    }

    /**
     * @brief Splits a line into whitespace (or comma) separated tokens, ignoring anything after '#'.
     * @tparam N     Maximum number of tokens stored.
     * @param start  Start of the line.
     * @param stop   One past the end of the line.
     * @param tokens Output array of tokens.
     * @return Number of tokens found; may exceed N, in which case only the first N are stored.
     */
    template <std::size_t N>
    static std::size_t tokenize(const char* start, const char* stop, std::array<std::string_view, N>& tokens)
    {
        const auto isSep = [](char c) { return std::isspace(static_cast<unsigned char>(c)) || c == ','; };
        std::size_t n = 0;
        while (start != stop && *start != '#') {
            while (start != stop && isSep(*start))
                ++start;
            if (start == stop || *start == '#')
                break;
            auto wend = start;
            while (wend != stop && !isSep(*wend) && *wend != '#')
                ++wend;
            if (n < N)
                tokens[n] = std::string_view(start, wend - start);
            ++n;
            start = wend;
        }
        return n;
    }

    /**
     * @brief Parses a complete token as a number.
     * @return true if the whole token was consumed and the value fits in @p val.
     */
    template <typename U>
    static bool parseToken(std::string_view token, U& val)
    {
        const auto end = token.data() + token.size();
        auto [ptr, ec] = std::from_chars(token.data(), end, val);
        return ec == std::errc() && ptr == end;
    }

    /**
     * @brief Parses an unsigned integer token, also accepting integral floating point
     *        notation such as "3.0" (TetGen may write region attributes as floats).
     */
    static bool parseIndexToken(std::string_view token, std::uint32_t& val)
    {
        if (parseToken(token, val))
            return true;
        double d;
        if (parseToken(token, d) && d >= 0 && d <= std::numeric_limits<std::uint32_t>::max() && std::floor(d) == d) {
            val = static_cast<std::uint32_t>(d);
            return true;
        }
        return false;
    }

    /**
     * @brief Reads the counts on the header line of a TetGen file.
     * @return The leading integer tokens of the header line (empty if none).
     */
    static std::vector<std::uint64_t> readHeaderCounts(const std::string& data, const std::pair<std::size_t, std::size_t>& line)
    {
        std::array<std::string_view, 4> tokens;
        const auto n = std::min(tokenize(data.data() + line.first, data.data() + line.second, tokens), tokens.size());
        std::vector<std::uint64_t> counts;
        for (std::size_t i = 0; i < n; ++i) {
            std::uint64_t v;
            if (!parseToken(tokens[i], v))
                break;
            counts.push_back(v);
        }
        return counts;
    }

    /// @brief Parse status of one line in a node or element file.
    enum class LineStatus : std::uint8_t {
        Empty, ///< Blank or comment-only line.
        Ok, ///< Line parsed successfully.
        Malformed ///< Line has content that could not be parsed.
    };

    /**
     * @brief Checks that sorted file indices form a contiguous sequence starting at 0 or 1.
     * @param indices Sorted index sequence.
     * @param what    Name used in the error message ("node" or "element").
     * @param path    File path used in the error message.
     * @return Empty string if valid, otherwise an error message.
     */
    static std::string checkIndexSequence(const std::vector<std::uint32_t>& indices, const std::string& what, const std::string& path)
    {
        if (indices.empty())
            return "No " + what + "s found in " + path;
        if (indices[0] > 1)
            return what + " indices in " + path + " must start at 0 or 1, first index is " + std::to_string(indices[0]);
        for (std::size_t i = 1; i < indices.size(); ++i) {
            if (indices[i] == indices[i - 1])
                return "Duplicate " + what + " index " + std::to_string(indices[i]) + " in " + path;
            if (indices[i] != indices[i - 1] + 1)
                return "Missing " + what + " index " + std::to_string(indices[i - 1] + 1) + " in " + path;
        }
        return {};
    }

    /**
     * @brief Reads a TetGen .ele element file into m_data.elements and m_data.collectionIndices.
     *
     * Expected line format: `<index> <n0> <n1> <n2> <n3> [collection_index]`. Text after '#'
     * is ignored. Element indices may start at 0 or 1 and must be contiguous without
     * duplicates; if the header line starts with a count it must match the number of
     * elements read. Node references are shifted by the index base of the node file, so
     * readNodes() must be called first. Either all or no lines may carry a collection index.
     * @param path         Path to the .ele file.
     * @param nHeaderLines Number of header lines to skip; default 1.
     * @param collength    Estimated average line length used for pre-allocating the line index.
     * @return true on success; on failure m_message describes the problem and no elements are stored.
     */
    bool readElements(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index n0 n1 n2 n3 matIdx000>, i.e "512 51 80 90 101"
        const auto data = readBufferFromFile(path);
        if (data.size() == 0) {
            m_message = "Could not read element file: " + path;
            return false;
        }

        const auto lineIdx = splitLines(data, collength);
        const auto nHeader = static_cast<std::size_t>(std::max(nHeaderLines, 0));
        if (lineIdx.size() <= nHeader) {
            m_message = "Element file contains no elements: " + path;
            return false;
        }
        const auto header = nHeader > 0 ? readHeaderCounts(data, lineIdx[0]) : std::vector<std::uint64_t> {};
        // TetGen header: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
        const std::size_t nodesPerTet = header.size() > 1 && header[1] >= 4 ? static_cast<std::size_t>(header[1]) : 4;
        const std::uint32_t nodeBase = m_nodeIndexBase;

        struct ElementLine {
            std::uint32_t index = 0;
            std::array<std::uint32_t, 4> nodes = { 0, 0, 0, 0 };
            std::uint32_t collection = 0;
            bool hasCollection = false;
            LineStatus status = LineStatus::Empty;
        };
        std::vector<ElementLine> lines(lineIdx.size() - nHeader);

        std::transform(std::execution::par_unseq, lineIdx.cbegin() + nHeader, lineIdx.cend(), lines.begin(), [&data, nodesPerTet, nodeBase](const auto& idx) {
            ElementLine e;
            std::array<std::string_view, 16> tokens;
            const auto n = tokenize(data.data() + idx.first, data.data() + idx.second, tokens);
            if (n == 0)
                return e;
            e.status = LineStatus::Malformed;
            if (n < nodesPerTet + 1 || nodesPerTet + 1 > tokens.size())
                return e;
            if (!parseToken(tokens[0], e.index))
                return e;
            for (std::size_t i = 0; i < 4; ++i) {
                if (!parseToken(tokens[i + 1], e.nodes[i]) || e.nodes[i] < nodeBase)
                    return e;
                e.nodes[i] -= nodeBase;
            }
            if (n > nodesPerTet + 1) {
                if (nodesPerTet + 1 >= tokens.size() || !parseIndexToken(tokens[nodesPerTet + 1], e.collection))
                    return e;
                e.hasCollection = true;
            }
            e.status = LineStatus::Ok;
            return e;
        });

        const auto malformed = std::find_if(lines.cbegin(), lines.cend(), [](const auto& e) { return e.status == LineStatus::Malformed; });
        if (malformed != lines.cend()) {
            const auto nMalformed = std::count_if(lines.cbegin(), lines.cend(), [](const auto& e) { return e.status == LineStatus::Malformed; });
            const auto lineNumber = std::distance(lines.cbegin(), malformed) + nHeader + 1;
            m_message = std::to_string(nMalformed) + " malformed line(s) in element file " + path + ", first at line " + std::to_string(lineNumber);
            return false;
        }
        std::erase_if(lines, [](const auto& e) { return e.status == LineStatus::Empty; });

        const auto nWithCollection = std::count_if(lines.cbegin(), lines.cend(), [](const auto& e) { return e.hasCollection; });
        if (nWithCollection != 0 && static_cast<std::size_t>(nWithCollection) != lines.size()) {
            m_message = std::to_string(lines.size() - nWithCollection) + " element(s) in " + path + " lack a collection index";
            return false;
        }

        // sort elements
        std::sort(std::execution::par_unseq, lines.begin(), lines.end(), [](const auto& lh, const auto& rh) { return lh.index < rh.index; });
        std::vector<std::uint32_t> indices(lines.size());
        std::transform(lines.cbegin(), lines.cend(), indices.begin(), [](const auto& e) { return e.index; });
        if (auto error = checkIndexSequence(indices, "element", path); !error.empty()) {
            m_message = error;
            return false;
        }
        if (!header.empty() && header[0] != lines.size()) {
            m_message = "Element file " + path + " header states " + std::to_string(header[0]) + " elements, but " + std::to_string(lines.size()) + " were read";
            return false;
        }

        // setting data to mesh data
        m_data.elements.resize(lines.size());
        std::transform(std::execution::par_unseq, lines.cbegin(), lines.cend(), m_data.elements.begin(), [](const auto& e) {
            return e.nodes;
        });
        m_data.collectionIndices.resize(lines.size());
        std::transform(std::execution::par_unseq, lines.cbegin(), lines.cend(), m_data.collectionIndices.begin(), [](const auto& e) {
            return e.collection;
        });
        return true;
    }

    /**
     * @brief Reads a TetGen .node vertex file into m_data.nodes.
     *
     * Expected line format: `<index> <x> <y> <z> [attributes...]`. Text after '#' is
     * ignored. Node indices may start at 0 or 1 (the base is remembered and applied to
     * the element file) and must be contiguous without duplicates; if the header line
     * starts with a count it must match the number of nodes read.
     * @param path         Path to the .node file.
     * @param nHeaderLines Number of header lines to skip; default 1.
     * @param collength    Estimated average line length used for pre-allocating the line index.
     * @return true on success; on failure m_message describes the problem and no nodes are stored.
     */
    bool readNodes(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index v0 v1 v2>, i.e "512 0.2 0.4523 -0.974"
        // # is treated as comment start
        const auto data = readBufferFromFile(path);
        if (data.size() == 0) {
            m_message = "Could not read node file: " + path;
            return false;
        }

        const auto lineIdx = splitLines(data, collength);
        const auto nHeader = static_cast<std::size_t>(std::max(nHeaderLines, 0));
        if (lineIdx.size() <= nHeader) {
            m_message = "Node file contains no nodes: " + path;
            return false;
        }
        const auto header = nHeader > 0 ? readHeaderCounts(data, lineIdx[0]) : std::vector<std::uint64_t> {};

        struct NodeLine {
            std::uint32_t index = 0;
            std::array<double, 3> pos = { 0, 0, 0 };
            LineStatus status = LineStatus::Empty;
        };
        std::vector<NodeLine> lines(lineIdx.size() - nHeader);

        std::transform(std::execution::par_unseq, lineIdx.cbegin() + nHeader, lineIdx.cend(), lines.begin(), [&data](const auto& idx) {
            NodeLine v;
            std::array<std::string_view, 4> tokens;
            const auto n = tokenize(data.data() + idx.first, data.data() + idx.second, tokens);
            if (n == 0)
                return v;
            v.status = LineStatus::Malformed;
            if (n < 4 || !parseToken(tokens[0], v.index))
                return v;
            for (std::size_t i = 0; i < 3; ++i)
                if (!parseToken(tokens[i + 1], v.pos[i]) || !std::isfinite(v.pos[i]))
                    return v;
            v.status = LineStatus::Ok;
            return v;
        });

        const auto malformed = std::find_if(lines.cbegin(), lines.cend(), [](const auto& v) { return v.status == LineStatus::Malformed; });
        if (malformed != lines.cend()) {
            const auto nMalformed = std::count_if(lines.cbegin(), lines.cend(), [](const auto& v) { return v.status == LineStatus::Malformed; });
            const auto lineNumber = std::distance(lines.cbegin(), malformed) + nHeader + 1;
            m_message = std::to_string(nMalformed) + " malformed line(s) in node file " + path + ", first at line " + std::to_string(lineNumber);
            return false;
        }
        std::erase_if(lines, [](const auto& v) { return v.status == LineStatus::Empty; });

        // sort vertices
        std::sort(std::execution::par_unseq, lines.begin(), lines.end(), [](const auto& lh, const auto& rh) { return lh.index < rh.index; });
        std::vector<std::uint32_t> indices(lines.size());
        std::transform(lines.cbegin(), lines.cend(), indices.begin(), [](const auto& v) { return v.index; });
        if (auto error = checkIndexSequence(indices, "node", path); !error.empty()) {
            m_message = error;
            return false;
        }
        if (!header.empty() && header[0] != lines.size()) {
            m_message = "Node file " + path + " header states " + std::to_string(header[0]) + " nodes, but " + std::to_string(lines.size()) + " were read";
            return false;
        }
        m_nodeIndexBase = indices[0];

        // populating mesh data
        m_data.nodes.resize(lines.size());
        std::transform(std::execution::par_unseq, lines.cbegin(), lines.cend(), m_data.nodes.begin(), [](const auto& v) {
            return v.pos;
        });
        return true;
    }

    /**
     * @brief Reads the entire contents of a file into a string.
     * @param path Path to the file to read.
     * @return File contents as a string, or an empty string if the file could not be opened.
     */
    static std::string readBufferFromFile(const std::string& path)
    {
        std::string buffer_str;

        // Open the file for reading
        std::ifstream file(path);
        if (file.is_open()) {
            std::stringstream buffer;
            buffer << file.rdbuf();
            buffer_str = buffer.str();
            file.close();
        }
        return buffer_str;
    }

private:
    TetrahedalMeshData m_data;
    bool m_valid = false;
    std::string m_message;
    std::uint32_t m_nodeIndexBase = 0; ///< Index of the first node in the node file (0 or 1).
};
}
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
     * Uses ICRP pregnant-woman format by default; pass the three-argument readData()
     * overload with @p ietr_format = true for IETR CSV format.
     * @param nodeFile         Path to the TetGen .node file.
     * @param elementFile      Path to the TetGen .ele file.
     * @param matorganfilePath Path to the combined material/organ file.
     */
    TetrahedalMeshReader(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matorganfilePath)
    {
        readData(nodeFile, elementFile, matorganfilePath);
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
        readNodes(nodeFile);
        readElements(elementFile);

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

        m_valid = mergeTetAndCollData(collection);
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
        if (ietr_format) {
            auto coll = readIETROrganAndMaterial(matorganfilePath);
            readNodes(nodeFile);
            readElements(elementFile);
            m_valid = mergeTetAndCollData(coll);
        } else {
            auto coll = readICRPPregnantOrganAndMaterial(matorganfilePath);
            readNodes(nodeFile);
            readElements(elementFile);
            m_valid = mergeTetAndCollData(coll);
        }
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
        readNodes(nodeFile);
        readElements(elementFile);
        m_valid = mergeTetAndCollData(coll);
    }

    /// @brief Returns a mutable reference to the parsed mesh data.
    TetrahedalMeshData& data()
    {
        return m_data;
    }

    /**
     * @brief Returns true if the last readData() call succeeded and the mesh has nodes and elements.
     *
     * A false result indicates that one or more collections could not be matched with
     * density or material data, or that the node/element files were empty or malformed.
     */
    bool valid() const
    {
        return m_valid && m_data.nodes.size() > 0 && m_data.elements.size() > 0;
    }

    /// @brief Resets all parsed data, leaving the reader in the same state as after default construction.
    void clear()
    {
        m_data = TetrahedalMeshData { };
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
        std::uint32_t index = 1;                         ///< Collection index as read from the element file.
        std::uint32_t material_index = 0;                ///< Internal material cross-reference index (ICRP 145).
        double density = -1;                             ///< Mass density in g/cm³; -1 means not yet assigned.
        std::map<std::uint8_t, double> material_weights; ///< Elemental weight fractions keyed by atomic number.
        std::string name;                                ///< Human-readable collection name.
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

        const auto lut = [&collection](const std::vector<std::uint32_t>& collIdx) {
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
        auto loc = std::find_if(m_data.collectionDensities.cbegin(), m_data.collectionDensities.cend(), [](const auto& dens) { return dens < 0.0; });
        return loc == m_data.collectionDensities.cend();
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
     * @brief Reads a TetGen .ele element file into m_data.elements and m_data.collectionIndices.
     *
     * Expected line format: `<index> <n0> <n1> <n2> <n3> [collection_index]`.
     * Lines are parsed in parallel; results are sorted by index and any sentinel-valued
     * (invalid) entries are removed. Indexing must start at 0.
     * @param path         Path to the .ele file.
     * @param nHeaderLines Number of header lines to skip; default 1.
     * @param collength    Estimated average line length used for pre-allocating the line index.
     */
    void readElements(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index n0 n1 n2 n3 matIdx000>, i.e "512 51 80 90 101"

        std::vector<std::tuple<std::uint32_t, std::array<std::uint32_t, 4>, std::uint32_t>> nodes;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return;

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        lineIdx.reserve(data.size() / collength);
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                lineIdx.push_back({ start, i });
                start = i;
            }
        }
        auto begin = lineIdx.cbegin() + nHeaderLines;
        auto end = lineIdx.cend();
        if (end - begin < 4)
            return;

        nodes.resize(end - begin);

        std::transform(std::execution::par_unseq, begin, end, nodes.begin(), [&data](const auto& idx) {
            auto start = data.data() + idx.first;
            auto stop = data.data() + idx.second;

            std::uint32_t index = std::numeric_limits<uint32_t>::max();
            std::array<std::uint32_t, 4> v {
                std::numeric_limits<std::uint32_t>::max(),
                std::numeric_limits<std::uint32_t>::max(),
                std::numeric_limits<std::uint32_t>::max(),
                std::numeric_limits<std::uint32_t>::max()
            };
            std::uint32_t collection = 0;

            parseLine(start, stop, ' ', index, v[0], v[1], v[2], v[3], collection);
            return std::make_tuple(index, v, collection);
        });

        // sort nodes
        std::sort(std::execution::par_unseq, nodes.begin(), nodes.end(), [](const auto& lh, const auto& rh) { return std::get<0>(lh) < std::get<0>(rh); });

        // removing nan values and validating;
        auto delete_from = nodes.cbegin();
        while (delete_from != nodes.cend() && std::get<0>(*delete_from) != std::numeric_limits<uint32_t>::max()) {
            ++delete_from;
        }

        if (delete_from != nodes.cend())
            nodes.erase(delete_from, nodes.cend());

        // We dont allow index start from other than 0
        if (nodes.size() > 0)
            if (std::get<0>(nodes[0]) != 0)
                return;

        // setting data to mesh data
        m_data.elements.resize(nodes.size());
        std::transform(std::execution::par_unseq, nodes.cbegin(), nodes.cend(), m_data.elements.begin(), [](const auto& e) {
            return std::get<1>(e);
        });
        m_data.collectionIndices.resize(nodes.size());
        std::transform(std::execution::par_unseq, nodes.cbegin(), nodes.cend(), m_data.collectionIndices.begin(), [](const auto& e) {
            return std::get<2>(e);
        });
    }

    /**
     * @brief Reads a TetGen .node vertex file into m_data.nodes.
     *
     * Expected line format: `<index> <x> <y> <z>`. Lines starting with '#' are treated
     * as comments. Lines are parsed in parallel; results are sorted by index and any
     * NaN-valued entries are removed. Indexing must start at 0.
     * @param path         Path to the .node file.
     * @param nHeaderLines Number of header lines to skip; default 1.
     * @param collength    Estimated average line length used for pre-allocating the line index.
     */
    void readNodes(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index v0 v1 v2>, i.e "512 0.2 0.4523 -0.974"
        // # is treated as comment start

        std::vector<std::pair<std::uint32_t, std::array<double, 3>>> vertices;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return;

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        lineIdx.reserve(data.size() / collength);
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                lineIdx.push_back({ start, i });
                start = i + 1;
            }
        }
        auto begin = lineIdx.cbegin() + nHeaderLines;
        auto end = lineIdx.cend();
        if (end - begin < 4)
            return;

        vertices.resize(end - begin);

        std::transform(std::execution::par_unseq, begin, end, vertices.begin(), [&data](const auto& idx) {
            auto start = data.data() + idx.first;
            auto stop = data.data() + idx.second;
            std::uint32_t index = std::numeric_limits<std::uint32_t>::max();
            std::array<double, 3> v = {
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()
            };
            parseLine(start, stop, ' ', index, v[0], v[1], v[2]);
            return std::make_pair(index, v);
        });

        // sort vertices
        std::sort(std::execution::par_unseq, vertices.begin(), vertices.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        // removing nan values and validating;
        auto delete_from = vertices.cbegin();
        while (delete_from != vertices.cend() && delete_from->first != std::numeric_limits<uint32_t>::max()) {
            ++delete_from;
        }

        if (delete_from != vertices.cend())
            vertices.erase(delete_from, vertices.cend());

        // We dont allow index start from other than 0
        if (vertices.size() > 0)
            if (vertices[0].first != 0)
                return;

        // populating mesh data
        m_data.nodes.resize(vertices.size());
        std::transform(std::execution::par_unseq, vertices.cbegin(), vertices.cend(), m_data.nodes.begin(), [](const auto& e) {
            return e.second;
        });
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
};
}
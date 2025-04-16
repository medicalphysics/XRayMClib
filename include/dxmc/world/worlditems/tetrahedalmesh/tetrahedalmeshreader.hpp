/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/material/material.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <execution>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace dxmc {

template <std::size_t Nshells = 5, int LOWENERGYCORRECTION = 2, bool FLUENCESCORING = true>
class TetrahedalmeshReader {
public:
    TetrahedalmeshReader() { }

    TetrahedalmeshReader(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matfilePath,
        const std::string& organFilePath)
    {
        readICRP145Phantom(nodeFile, elementFile, matfilePath, organFilePath);
    }

    TetrahedalmeshReader(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matorganfilePath)
    {
        readICRPPregnantPhantom(nodeFile, elementFile, matorganfilePath);
    }

    void readICRP145Phantom(
        const std::string& nodeFile,
        const std::string& elementFile,
        const std::string& matfilePath,
        const std::string& organFilePath)
    {
        auto orgs = readICRP145PhantomOrgans(organFilePath);
        auto mats = readICRP145PhantomMaterials(matfilePath);

        // harmonizing organ and material indices
        std::unordered_map<std::uint16_t, std::uint16_t> mat_lut;
        for (std::uint16_t i = 0; i < mats.size(); ++i) {
            const auto& m = mats[i];
            const auto key = m.index;
            mat_lut[key] = i;
        }
        for (auto& o : orgs) {
            o.materialIdx = mat_lut[o.materialIdx];
        }
        std::unordered_map<std::uint16_t, std::uint16_t> org_lut;
        for (std::uint16_t i = 0; i < orgs.size(); ++i) {
            const auto& o = orgs[i];
            const auto key = o.index;
            org_lut[key] = i;
        }

        m_tets = readICRPPhantomGeometry(nodeFile, elementFile);

        // updating material ond organ indices in tets
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&orgs, &org_lut](auto& t) {
            const auto oldOidx = t.collection();
            const auto newOidx = org_lut[oldOidx];
            t.setCollection(newOidx);
            t.setMaterialIndex(orgs[newOidx].materialIdx);
        });

        m_densities.resize(orgs.size());
        m_organNames.resize(orgs.size());
        for (std::size_t i = 0; i < orgs.size(); ++i) {
            const auto matIdx = orgs[i].materialIdx;
            m_densities[i] = mats[matIdx].density;
            m_organNames[i] = orgs[i].name;
        }

        m_materials.reserve(mats.size());
        for (const auto& m : mats)
            m_materials.emplace_back(m.material);
    }

    void readICRPPregnantPhantom(const std::string& nodeFile, const std::string& elementFile, const std::string& matorganfilePath)
    {
        auto organs = readICRPPregnantOrganMaterial(matorganfilePath);

        // splitting materials and setting material index
        {
            std::vector<std::map<std::size_t, double>> material_weights;
            material_weights.reserve(organs.size());
            for (std::size_t i = 0; i < organs.size(); ++i) {
                const auto& organ = organs[i];
                auto idx_pos = std::find(material_weights.cbegin(), material_weights.cend(), organ.material_weights);
                if (idx_pos == material_weights.cend()) {
                    organs[i].material_index = material_weights.size();
                    material_weights.push_back(organ.material_weights);
                } else {
                    organs[i].material_index = std::distance(material_weights.cbegin(), idx_pos);
                }
            }

            m_materials.clear();
            for (const auto& w : material_weights)
                m_materials.emplace_back(Material<Nshells>::byWeight(w).value());
        }

        m_densities.clear();
        m_densities.reserve(organs.size());
        m_organNames.clear();
        m_organNames.reserve(organs.size());
        for (const auto& o : organs) {
            m_densities.push_back(o.density);
            m_organNames.push_back(o.name);
        }

        m_tets = readICRPPhantomGeometry(nodeFile, elementFile);

        // map organ to collection index
        std::map<std::size_t, std::uint16_t> organIdx_map, materialIdx_map;
        for (std::uint16_t i = 0; i < organs.size(); ++i) {
            organIdx_map[organs[i].index] = i;
            materialIdx_map[organs[i].index] = organs[i].material_index;
        }
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&organIdx_map, &materialIdx_map](auto& tet) {
            tet.setMaterialIndex(materialIdx_map.at(tet.collection()));
            tet.setCollection(organIdx_map.at(tet.collection()));
        });
    }

    TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> getMesh(int depth = 8)
    {
        TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> mesh(m_tets, m_densities, m_materials, m_organNames, depth);
        return mesh;
    }

    TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> getMesh(int x, int y, int z)
    {
        std::array<int, 3> d = { x, y, z };
        TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> mesh(m_tets, m_densities, m_materials, m_organNames, d);
        return mesh;
    }

    TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> getMesh(const std::array<int, 3>& depth)
    {
        TetrahedalMesh<Nshells, LOWENERGYCORRECTION, FLUENCESCORING> mesh(m_tets, m_densities, m_materials, m_organNames, depth);
        return mesh;
    }

    void rotate(const std::array<double, 3>& axis, double angle)
    {
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&axis, angle](auto& v) {
            v.rotate(axis, angle);
        });
    }

protected:
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

    struct ICRP145Materials {
        std::uint16_t index = 0;
        Material<Nshells> material;
        double density = 1;
        std::string name;

        ICRP145Materials(Material<Nshells>& mat)
            : material(mat)
        {
        }
    };

    static std::vector<ICRP145Materials> readICRP145PhantomMaterials(const std::string& matfilePath)
    {
        std::vector<ICRP145Materials> res;
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
            std::map<std::size_t, double> frac;
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

            double organDensity = 1;
            if (frac.contains(0)) {
                organDensity = frac.at(0);
                frac.erase(0);
            }
            auto mat = Material<Nshells>::byWeight(frac);
            if (mat) {
                ICRP145Materials organ(mat.value());
                organ.density = organDensity;
                organ.index = organIndex;
                organ.name = organName;
                res.push_back(organ);
            }
        }
        return res;
    }

    struct ICRP145Organs {
        std::uint16_t index = 0;
        std::uint16_t materialIdx = 0;
        std::string name;
    };

    static std::vector<ICRP145Organs> readICRP145PhantomOrgans(const std::string& organfilePath)
    {
        std::vector<ICRP145Organs> res;
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
            parseLine(data.data() + lineIdx[i].first, data.data() + lineIdx[i].second, ',', organ.index, organ.name, organ.materialIdx);
        }
        return res;
    }

    struct ICRPPregnantOrgan {
        std::size_t index = 0;
        std::size_t material_index = 0;
        double density = 0;
        std::map<std::size_t, double> material_weights;
        std::string name;
    };

    static std::vector<ICRPPregnantOrgan> readICRPPregnantOrganMaterial(const std::string& matfile)
    {
        std::vector<ICRPPregnantOrgan> organs;
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
            ICRPPregnantOrgan organ;
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
                if (begin < end)
                    organ.material_weights[Z / 1000] = std::abs(w);
            }
            return organ;
        });

        return organs;
    }

    static std::vector<Tetrahedron> readICRPPhantomGeometry(const std::string& nodeFile, const std::string& elementFile)
    {
        auto vertices = readVertices(nodeFile, 1, 80);
        auto nodes = readTetrahedalIndices(elementFile, 1, 80);

        const auto max_ind = std::transform_reduce(
            std::execution::par_unseq, nodes.cbegin(), nodes.cend(), std::size_t { 0 },
            [](const auto lh, const auto rh) { return std::max(rh, lh); }, [](const auto& node) {
                const auto& v = std::get<1>(node);
                return std::max(v[0], std::max(v[1], v[2])); });

        std::vector<Tetrahedron> tets;
        if (max_ind < vertices.size()) {
            tets.resize(nodes.size());
            std::transform(std::execution::par_unseq, nodes.cbegin(), nodes.cend(), tets.begin(), [&vertices](const auto& n) {
                const auto& vIdx = std::get<1>(n);
                const auto collection = static_cast<std::uint16_t>(std::get<2>(n));

                const auto& v0 = vertices[vIdx[0]].second;
                const auto& v1 = vertices[vIdx[1]].second;
                const auto& v2 = vertices[vIdx[2]].second;
                const auto& v3 = vertices[vIdx[3]].second;

                return Tetrahedron { v0, v1, v2, v3, collection };
            });
        }

        return tets;
    }

    static std::vector<std::tuple<std::size_t, std::array<std::size_t, 4>, std::size_t>> readTetrahedalIndices(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index n0 n1 n2 n3 matIdx000>, i.e "512 51 80 90 101"

        std::vector<std::tuple<std::size_t, std::array<std::size_t, 4>, std::size_t>> nodes;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return nodes;

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
            return nodes;

        nodes.resize(end - begin);

        std::transform(std::execution::par_unseq, begin, end, nodes.begin(), [&data](const auto& idx) {
            auto start = data.data() + idx.first;
            auto stop = data.data() + idx.second;

            std::size_t index = std::numeric_limits<std::size_t>::max();
            std::array<std::size_t, 4> v {
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max()
            };
            std::size_t collection = std::numeric_limits<std::size_t>::max();

            parseLine(start, stop, ' ', index, v[0], v[1], v[2], v[3], collection);
            return std::make_tuple(index, v, collection);
        });

        // sort nodes
        std::sort(std::execution::par_unseq, nodes.begin(), nodes.end(), [](const auto& lh, const auto& rh) { return std::get<0>(lh) < std::get<0>(rh); });

        // removing nan values and validating;
        auto delete_from = nodes.cbegin();
        std::size_t index = 0;
        while (delete_from != nodes.cend() && index == std::get<0>(*delete_from)) {
            ++index;
            ++delete_from;
        }

        if (delete_from != nodes.cend())
            nodes.erase(delete_from, nodes.cend());

        return nodes;
    }

    static std::vector<std::pair<std::size_t, std::array<double, 3>>> readVertices(const std::string& path, int nHeaderLines = 1, std::size_t collength = 80)
    {
        // reads a file formatted as <index v0 v1 v2>, i.e "512 0.2 0.4523 -0.974"
        // # is treated as comment start

        std::vector<std::pair<std::size_t, std::array<double, 3>>> vertices;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return vertices;

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
            return vertices;

        vertices.resize(end - begin);

        std::transform(std::execution::par_unseq, begin, end, vertices.begin(), [&data](const auto& idx) {
            auto start = data.data() + idx.first;
            auto stop = data.data() + idx.second;
            std::size_t index = std::numeric_limits<std::size_t>::max();
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
        std::size_t index = 0;
        while (delete_from != vertices.cend() && index == delete_from->first) {
            ++index;
            ++delete_from;
        }

        if (delete_from != vertices.cend())
            vertices.erase(delete_from, vertices.cend());

        return vertices;
    }

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
    std::vector<Tetrahedron> m_tets;
    std::vector<double> m_densities;
    std::vector<Material<Nshells>> m_materials;
    std::vector<std::string> m_organNames;
};
}
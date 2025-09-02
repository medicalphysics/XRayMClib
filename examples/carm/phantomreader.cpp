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

#include "phantomreader.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <execution>
#include <fstream>
#include <numeric>
#include <optional>
#include <sstream>

// trim from start (in place)
inline void ltrim(std::string& s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(),
        s.end());
}

// trim from both ends (in place)
inline void trim(std::string& s)
{
    rtrim(s);
    ltrim(s);
}

auto to_uchar(std::string_view s) -> std::optional<std::uint8_t>
{
    std::uint8_t value {};
    if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
        return value;
    else
        return std::nullopt;
};
auto to_double(std::string_view s) -> std::optional<double>
{
    double value {};
    if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
        return value;
    else
        return std::nullopt;
};

std::vector<std::uint8_t> readASCIIData(const std::string& path)
{
    std::ifstream t;
    std::vector<std::uint8_t> buffer;
    t.open(path); // open input file
    if (t.is_open()) {
        t.seekg(0, std::ios::end); // go to the end
        auto length = t.tellg(); // report location (this is the length)

        t.seekg(0, std::ios::beg); // go back to the beginning
        std::string test(length, '0');
        t.read(&test[0], length);
        t.close();

        test.erase(std::remove_if(std::execution::par_unseq, test.begin(), test.end(), [](unsigned char x) { return x == '\t' || x == '\n'; }), test.end());

        const auto buffer_size = test.size() / 4;

        buffer.resize(buffer_size, 255);
        for (std::size_t i = 0; i < buffer.size(); ++i) {
            char* start = test.data() + 4 * i;
            char* end = test.data() + 4 * (i + 1);

            while (std::isspace(*start))
                ++start;
            std::from_chars(start, end, buffer[i]);
        }
    }
    return buffer;
}

struct Organ {
    std::string name;
    double density = 0;
    std::uint8_t ID = 0;
    std::uint8_t mediaID = 0;
};

std::vector<Organ> readASCIIOrgans(const std::string& organ_path)
{
    std::vector<Organ> organs;
    std::ifstream t(organ_path);
    std::stringstream buffer;
    if (t.is_open()) {
        buffer << t.rdbuf();
        t.close();
    } else {
        return organs;
    }

    std::size_t lineNr = 0;

    for (std::string line; std::getline(buffer, line);) {
        lineNr++;
        if (line.size() >= 66 && lineNr > 4) {
            Organ organ;
            // ID
            auto IDs = line.substr(0, 6);
            auto ID = to_uchar(IDs);
            if (ID)
                organ.ID = ID.value();
            auto name = line.substr(6, 49);
            trim(name);
            organ.name = name;
            auto TNRs = line.substr(55, 3);
            ltrim(TNRs);
            auto TNR = to_uchar(TNRs);
            if (TNR)
                organ.mediaID = TNR.value();
            auto denss = line.substr(61, 5);
            auto dens = to_double(denss);
            if (dens)
                organ.density = dens.value();

            organs.push_back(organ);
        }
    }
    organs.push_back({ .name = "Air", .density = 0.001, .ID = 0, .mediaID = 0 });

    return organs;
}

struct Media {
    std::uint8_t ID = 0;
    std::map<std::size_t, double> composition;
    std::string name;
};

std::vector<Media> readASCIIMedia(const std::string& media_path)
{
    std::vector<Media> medias;
    std::ifstream t(media_path);
    std::stringstream buffer;
    if (t.is_open()) {
        buffer << t.rdbuf();
        t.close();
    } else {
        return medias;
    }

    std::size_t lineNr = 0;

    for (std::string line; std::getline(buffer, line);) {
        lineNr++;
        if (line.size() >= 154 && lineNr > 3) {
            Media media;
            // ID
            auto IDs = line.substr(0, 6);
            auto ID = to_uchar(IDs);
            if (ID)
                media.ID = ID.value();
            auto name = line.substr(6, 72);
            trim(name);
            media.name = name;

            // composition
            const std::array Zs = { 1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20, 26, 53 };
            for (std::size_t i = 0; i < Zs.size(); ++i) {
                constexpr std::size_t offset = 76;
                constexpr std::size_t step = 6;
                auto s = line.substr(offset + step * i, step);
                trim(s);
                auto w = to_double(s);
                if (w)
                    if (w.value() > 0.0)
                        media.composition[Zs[i]] = w.value();
            }
            medias.push_back(media);
        }
    }

    medias.push_back({ .ID = 0, .composition = { { 7, 80.0 }, { 8, 20.0 } }, .name = "Air" });
    return medias;
}

void sanitizeIDs(std::vector<std::uint8_t>& organ_data, std::vector<Organ>& organs, std::vector<Media>& media)
{

    // do we have all organs
    {
        // finding organ IDs in organ data;
        auto organIDs = organ_data; // copy of vector
        std::sort(std::execution::par_unseq, organIDs.begin(), organIDs.end());
        organIDs.erase(std::unique(std::execution::par_unseq, organIDs.begin(), organIDs.end()), organIDs.end());

        std::sort(organs.begin(), organs.end(), [](const auto& lh, const auto& rh) { return lh.ID < rh.ID; });
        std::uint8_t teller = 0;
        while (teller < organIDs.size()) {
            const auto oID = organIDs[teller];
            auto pos = std::find_if(organs.cbegin(), organs.cend(), [oID](const auto& o) { return o.ID == oID; });
            if (pos == organs.cend()) { // could not find organ id in organ list set organ to 0
                organIDs[teller] = 0;
                std::replace(organ_data.begin(), organ_data.end(), oID, std::uint8_t { 0 });
                std::sort(organIDs.begin(), organIDs.end());
                organIDs.erase(std::unique(organIDs.begin(), organIDs.end()), organIDs.end());
            } else {
                ++teller;
            }
        }

        // remove organs not in data
        teller = 0;
        while (teller < organs.size()) {
            const auto oID = organs[teller].ID;
            auto pos = std::find(organIDs.cbegin(), organIDs.cend(), oID);
            if (pos == organIDs.cend()) { // extra organ, we delete it
                organs.erase(organs.begin() + teller);
            } else {
                ++teller;
            }
        }

        // making organ IDs consecutive
        std::sort(organs.begin(), organs.end(), [](const auto& lh, const auto& rh) { return lh.ID < rh.ID; });
        for (std::uint8_t i = 0; i < organs.size(); ++i) {
            if (organs[i].ID != i) {
                std::replace(std::execution::par_unseq, organ_data.begin(), organ_data.end(), organs[i].ID, i);
                organs[i].ID = i;
            }
        }
    }

    // making sure we have all media
    {
        std::vector<std::uint8_t> mediaIDs;
        for (const auto& o : organs)
            mediaIDs.push_back(o.mediaID);
        std::sort(mediaIDs.begin(), mediaIDs.end());
        mediaIDs.erase(std::unique(mediaIDs.begin(), mediaIDs.end()), mediaIDs.end());
        // do we miss any media, in case set media to air
        for (const auto id : mediaIDs) {
            auto pos = std::find_if(media.cbegin(), media.cend(), [id](const auto& m) { return m.ID == id; });
            if (pos == media.cend()) { // missing media, setting to air
                for (auto& o : organs) {
                    if (o.mediaID == id) {
                        o.mediaID = 0;
                    }
                }
            }
        }

        // removing unneeded media
        std::uint8_t teller = 0;
        while (teller < media.size()) {
            const auto id = media[teller].ID;
            auto pos = std::find_if(organs.cbegin(), organs.cend(), [id](const auto& o) { return o.mediaID == id; });
            if (pos == organs.cend()) { // unneeded media, remove
                media.erase(media.cbegin() + teller);
            } else {
                ++teller;
            }
        }
        std::sort(media.begin(), media.end(), [](const auto& lh, const auto& rh) { return lh.ID < rh.ID; });
        // making media consecutive
        for (std::uint8_t i = 0; i < media.size(); ++i) {
            if (media[i].ID != i) {
                for (auto& o : organs) {
                    if (o.mediaID == media[i].ID) {
                        o.mediaID = i;
                    }
                }
                media[i].ID = i;
            }
        }
    }
}

template <bool FEMALE>
ICRP110PhantomReader ICRP110PhantomReader::readPhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path)
{

    ICRP110PhantomReader data;

    if constexpr (FEMALE) {
        data.m_dim = { 299, 137, 348 };
        data.m_spacing = { .1775, .1775, .484 };
    } else {
        data.m_dim = { 254, 127, 222 };
        data.m_spacing = { .2137, .2137, 0.8 };
    }
    const auto size = std::reduce(data.m_dim.cbegin(), data.m_dim.cend(), std::size_t { 1 }, std::multiplies<>());
    data.m_organ_data = readASCIIData(phantom_path);

    if (size != data.m_organ_data.size())
        return data;

    auto organs = readASCIIOrgans(organ_path);
    auto media = readASCIIMedia(media_path);
    sanitizeIDs(data.m_organ_data, organs, media);

    std::vector<std::uint8_t> organ_to_media_map;
    for (const auto& o : organs)
        organ_to_media_map.push_back(o.mediaID);
    data.m_media_data.resize(data.m_organ_data.size());
    std::transform(std::execution::par_unseq, data.m_organ_data.cbegin(), data.m_organ_data.cend(), data.m_media_data.begin(), [organ_to_media_map](const auto oID) { return organ_to_media_map[oID]; });

    std::vector<double> organ_to_density_map;
    for (const auto& o : organs)
        organ_to_density_map.push_back(o.density);
    data.m_density_data.resize(data.m_organ_data.size());
    std::transform(std::execution::par_unseq, data.m_organ_data.cbegin(), data.m_organ_data.cend(), data.m_density_data.begin(), [organ_to_density_map](const auto oID) { return organ_to_density_map[oID]; });

    for (const auto& o : organs)
        data.m_organ_name.push_back(o.name);

    for (const auto& m : media)
        data.m_media_composition.push_back(m.composition);

    return data;
}

ICRP110PhantomReader ICRP110PhantomReader::readFemalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path)
{
    return readPhantom<true>(phantom_path, media_path, organ_path);
}

ICRP110PhantomReader ICRP110PhantomReader::readMalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path)
{
    return readPhantom<false>(phantom_path, media_path, organ_path);
}
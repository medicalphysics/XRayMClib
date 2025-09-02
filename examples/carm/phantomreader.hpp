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

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#pragma once

class ICRP110PhantomReader {
public:
    const std::vector<std::uint8_t>& mediaData() const { return m_media_data; }
    const std::vector<std::uint8_t>& organData() const { return m_organ_data; }
    const std::vector<double>& densityData() const { return m_density_data; }
    const std::vector<std::map<std::size_t, double>>& mediaComposition() const { return m_media_composition; }
    const std::vector<std::string>& organNames() const { return m_organ_name; }
    std::array<double, 3> spacing() const { return m_spacing; }
    std::array<std::size_t, 3> dimensions() const { return m_dim; }

    static ICRP110PhantomReader readFemalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path);
    static ICRP110PhantomReader readMalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path);

protected:
    template <bool FEMALE = true>
    static ICRP110PhantomReader readPhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path);
    ICRP110PhantomReader() = default;

private:
    std::vector<std::uint8_t> m_organ_data;
    std::vector<std::uint8_t> m_media_data;
    std::vector<double> m_density_data;
    std::array<double, 3> m_spacing;
    std::array<std::size_t, 3> m_dim;
    std::vector<std::string> m_organ_name;
    std::vector<std::map<std::size_t, double>> m_media_composition;
};

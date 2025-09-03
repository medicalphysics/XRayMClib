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

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace dxmc {
namespace dxmclodepng {

    std::vector<std::uint8_t> encodePNG(const std::vector<std::uint8_t>& image, std::size_t width, size_t height);

    template <typename T>
        requires std::is_same_v<T, float> || std::is_same_v<T, double>
    std::vector<std::uint8_t> encodePNG(const std::vector<T>& image, std::size_t width, size_t height)
    {
        std::vector<std::uint8_t> bytes(image.size());
        std::transform(image.cbegin(), image.cend(), bytes.begin(), [](const T& v) {
            return static_cast<std::uint8_t>(std::clamp(v, T { 0 }, T { 1 }) * 255);
        });
        return encodePNG(bytes, width, height);
    }

    bool savePNG(const std::string& filename, const std::vector<std::uint8_t>& encodedImage);

    template <typename T>
        requires std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, std::uint8_t>
    bool
    savePNG(const std::string& filename, const std::vector<T>& image, std::size_t width, size_t height)
    {
        return savePNG(filename, encodePNG(image, width, height));
    }

};
}
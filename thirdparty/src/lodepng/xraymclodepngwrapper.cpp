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

#include "lodepng.h"
#include "lodepng/xraymclodepngwrapper.hpp"

#include <fstream>

namespace xraymc {
namespace xraymclodepng {

    std::vector<std::uint8_t> encodePNG(const std::vector<std::uint8_t>& image, std::size_t width, size_t height)
    {
        std::vector<unsigned char> png;
        if (image.size() == width * height * 4) {
            lodepng::encode(png, image, width, height);
        }
        return png;
    }

    bool writeImage(const std::string& filename, const std::vector<std::uint8_t>& buffer)
    {
        std::ofstream file(filename, std::ios::out | std::ios::binary);
        if (file.is_open()) {
            file.write((char*)buffer.data(), buffer.size() * sizeof(std::uint8_t));
            file.close();
        } else {
            return false;
        }
        return true;
    }

    bool savePNG(const std::string& filename, const std::vector<std::uint8_t>& encodedImage)
    {
        if (encodedImage.size() == 0)
            return false;
        return writeImage(filename, encodedImage);
    }

};
}
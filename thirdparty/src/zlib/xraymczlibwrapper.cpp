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

#include "zlib/xraymczlibwrapper.hpp"

#include <zlib.h>

namespace xraymc {
namespace xraymczlib {

    std::expected<std::vector<char>, error> compress(std::span<const char> in, int level)
    {
        if (in.empty())
            return std::vector<char> {};

        uLongf destLen = compressBound(static_cast<uLong>(in.size()));
        std::vector<char> out(destLen);

        const int ret = compress2(
            reinterpret_cast<Bytef*>(out.data()), &destLen,
            reinterpret_cast<const Bytef*>(in.data()), static_cast<uLong>(in.size()),
            level);

        switch (ret) {
        case Z_OK:
            out.resize(destLen);
            return out;
        case Z_MEM_ERROR:
            return std::unexpected(error::memory_error);
        case Z_BUF_ERROR:
            return std::unexpected(error::buffer_error);
        default:
            return std::unexpected(error::unknown_error);
        }
    }

    std::expected<std::vector<char>, error> decompress(std::span<const char> in, std::uint64_t uncompressedSize)
    {
        if (uncompressedSize == 0)
            return std::vector<char> {};

        std::vector<char> out(uncompressedSize);
        uLongf destLen = static_cast<uLongf>(uncompressedSize);

        const int ret = uncompress(
            reinterpret_cast<Bytef*>(out.data()), &destLen,
            reinterpret_cast<const Bytef*>(in.data()), static_cast<uLong>(in.size()));

        switch (ret) {
        case Z_OK:
            out.resize(destLen);
            return out;
        case Z_MEM_ERROR:
            return std::unexpected(error::memory_error);
        case Z_BUF_ERROR:
            return std::unexpected(error::buffer_error);
        case Z_DATA_ERROR:
            return std::unexpected(error::data_error);
        default:
            return std::unexpected(error::unknown_error);
        }
    }

}
}

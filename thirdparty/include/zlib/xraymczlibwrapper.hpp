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

#pragma once

#include <cstdint>
#include <expected>
#include <span>
#include <vector>

namespace xraymc {
namespace xraymczlib {

    enum class error {
        memory_error, ///< zlib ran out of memory.
        buffer_error, ///< Output buffer was not large enough.
        data_error, ///< Compressed stream is corrupt or incomplete.
        unknown_error ///< Any other zlib failure.
    };

    /**
     * @brief Compresses @p in with zlib.
     * @param level 0 (no compression) .. 9 (best compression); -1 selects the zlib default.
     * @return Compressed bytes, or an error on failure.
     */
    std::expected<std::vector<char>, error> compress(std::span<const char> in, int level = -1);

    /**
     * @brief Decompresses @p in, which must inflate to exactly @p uncompressedSize bytes.
     * @param in               Compressed bytes, as produced by compress().
     * @param uncompressedSize Exact size of the original, uncompressed data.
     * @return Decompressed bytes, or an error on failure.
     */
    std::expected<std::vector<char>, error> decompress(std::span<const char> in, std::uint64_t uncompressedSize);

}
}

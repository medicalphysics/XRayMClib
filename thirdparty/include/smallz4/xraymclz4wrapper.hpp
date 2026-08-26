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

#include "smallz4/smallz4.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace xraymc {
namespace xraymclz4 {

    namespace detail {

        struct CompressContext {
            const std::vector<char>* input;
            std::size_t pos = 0;
            std::vector<char>* output;
        };

        inline std::size_t readBytes(void* data, std::size_t numBytes, void* userPtr)
        {
            auto* ctx = static_cast<CompressContext*>(userPtr);
            const std::size_t remaining = ctx->input->size() - ctx->pos;
            const std::size_t toCopy = std::min(numBytes, remaining);
            if (toCopy > 0) {
                std::memcpy(data, ctx->input->data() + ctx->pos, toCopy);
                ctx->pos += toCopy;
            }
            return toCopy;
        }

        inline void writeBytes(const void* data, std::size_t numBytes, void* userPtr)
        {
            auto* ctx = static_cast<CompressContext*>(userPtr);
            const auto* bytes = static_cast<const char*>(data);
            ctx->output->insert(ctx->output->end(), bytes, bytes + numBytes);
        }

        // decode a single raw LZ4 block, appending the result to out (out also
        // serves as the match window since smallz4 emits block-dependent frames)
        inline void decodeBlock(const char* src, std::size_t srcSize, std::vector<char>& out)
        {
            const auto* ip = reinterpret_cast<const std::uint8_t*>(src);
            const auto* ipEnd = ip + srcSize;

            while (ip < ipEnd) {
                const std::uint8_t token = *ip++;

                std::size_t literalLength = token >> 4;
                if (literalLength == 15) {
                    std::uint8_t s;
                    do {
                        if (ip >= ipEnd)
                            throw std::runtime_error("xraymclz4::decompress: corrupt block (literal length)");
                        s = *ip++;
                        literalLength += s;
                    } while (s == 255);
                }
                if (static_cast<std::size_t>(ipEnd - ip) < literalLength)
                    throw std::runtime_error("xraymclz4::decompress: corrupt block (literal overrun)");
                out.insert(out.end(), reinterpret_cast<const char*>(ip), reinterpret_cast<const char*>(ip) + literalLength);
                ip += literalLength;

                if (ip >= ipEnd)
                    break; // final sequence of a block is literals-only

                if (ipEnd - ip < 2)
                    throw std::runtime_error("xraymclz4::decompress: corrupt block (missing offset)");
                std::uint16_t offset;
                std::memcpy(&offset, ip, 2);
                ip += 2;
                if (offset == 0 || offset > out.size())
                    throw std::runtime_error("xraymclz4::decompress: corrupt block (bad offset)");

                std::size_t matchLength = token & 0x0F;
                if (matchLength == 15) {
                    std::uint8_t s;
                    do {
                        if (ip >= ipEnd)
                            throw std::runtime_error("xraymclz4::decompress: corrupt block (match length)");
                        s = *ip++;
                        matchLength += s;
                    } while (s == 255);
                }
                matchLength += 4; // minimum match length

                std::size_t matchPos = out.size() - offset;
                out.reserve(out.size() + matchLength);
                for (std::size_t i = 0; i < matchLength; ++i)
                    out.push_back(out[matchPos + i]);
            }
        }

    } // namespace detail

    /// Compress a buffer with smallz4 (LZ4 frame format).
    /// maxChainLength: 0 = store uncompressed, <=3 greedy, <=6 lazy, higher = optimal parsing (slower, smaller).
    inline std::vector<char> compress(const std::vector<char>& input, unsigned short maxChainLength = 65535)
    {
        std::vector<char> output;
        detail::CompressContext ctx { &input, 0, &output };
        smallz4::lz4(detail::readBytes, detail::writeBytes, maxChainLength, false, &ctx);
        return output;
    }

    /// Decompress a buffer produced by compress().
    inline std::vector<char> decompress(const std::vector<char>& compressed)
    {
        if (compressed.size() < 7)
            throw std::length_error("xraymclz4::decompress: buffer too small to contain a frame header");

        const auto* p = reinterpret_cast<const std::uint8_t*>(compressed.data());
        const auto* end = p + compressed.size();

        static constexpr std::uint8_t magic[4] = { 0x04, 0x22, 0x4D, 0x18 };
        if (std::memcmp(p, magic, 4) != 0)
            throw std::runtime_error("xraymclz4::decompress: bad magic number, not an LZ4 frame");
        p += 4;

        const std::uint8_t flags = *p++;
        p += 1; // block descriptor byte, not needed since each block is size-prefixed

        const bool blockChecksum = (flags >> 4) & 1;
        const bool contentSize = (flags >> 3) & 1;
        const bool contentChecksum = (flags >> 2) & 1;
        const bool dictID = flags & 1;

        if (contentSize)
            p += 8;
        if (dictID)
            p += 4;
        p += 1; // header checksum, not verified here

        std::vector<char> out;
        while (true) {
            if (end - p < 4)
                throw std::length_error("xraymclz4::decompress: truncated block size");
            std::uint32_t blockSizeTagged;
            std::memcpy(&blockSizeTagged, p, 4);
            p += 4;

            if (blockSizeTagged == 0)
                break; // end mark

            const bool uncompressed = (blockSizeTagged & 0x80000000u) != 0;
            const std::uint32_t blockSize = blockSizeTagged & 0x7FFFFFFFu;

            if (static_cast<std::uint32_t>(end - p) < blockSize)
                throw std::length_error("xraymclz4::decompress: truncated block data");

            if (uncompressed)
                out.insert(out.end(), reinterpret_cast<const char*>(p), reinterpret_cast<const char*>(p) + blockSize);
            else
                detail::decodeBlock(reinterpret_cast<const char*>(p), blockSize, out);

            p += blockSize;
            if (blockChecksum)
                p += 4;
        }

        if (contentChecksum)
            p += 4;

        return out;
    }

} // namespace xraymclz4
} // namespace xraymc

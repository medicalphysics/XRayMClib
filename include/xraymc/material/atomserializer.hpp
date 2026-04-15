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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "xraymc/material/atomicelement.hpp"
#include "xraymc/material/atomicshell.hpp"

#include <algorithm>
#include <concepts>
#include <iterator>
#include <map>
#include <vector>

namespace xraymc {
/**
 * @brief Concept satisfied by any integral or floating-point type.
 *
 * Used to constrain the numeric overloads of `AtomSerializer::serialize` and
 * `AtomSerializer::deserialize` to primitive arithmetic types only.
 */
template <typename T>
concept Number = std::is_integral<T>::value || std::is_floating_point<T>::value;

/**
 * @brief Binary serializer and deserializer for `AtomicElement` and `AtomicShell` data.
 *
 * Provides static methods to pack a map of `AtomicElement` objects (including all
 * their cross-section tables, orbital parameters, and fluorescence transition lists)
 * into a flat byte buffer, and to reconstruct them from such a buffer.
 *
 * The binary format is little-endian, length-prefixed: each variable-length block
 * (vector or sub-struct) is preceded by a `uint64_t` byte count so that the reader
 * can skip unknown fields or validate sizes.
 *
 * Public interface:
 * - `serializeAtoms`   — encodes the full element map to `vector<char>`.
 * - `deserializeAtoms` — decodes a buffer produced by `serializeAtoms`.
 *
 * Protected helpers handle individual scalars, vectors, shells, and elements and
 * are available to derived serializers.
 */
class AtomSerializer {
public:
    /**
     * @brief Serializes a map of `AtomicElement` objects to a byte buffer.
     *
     * Writes the element count followed by each element's data (via
     * `serializeAtomicElement`).
     *
     * @param elements  Map of `{Z → AtomicElement}` to serialize.
     * @return Byte buffer containing the encoded element data.
     */
    static std::vector<char> serializeAtoms(const std::map<std::uint64_t, AtomicElement>& elements)
    {
        std::vector<char> buffer;
        std::uint64_t n_elements = elements.size();
        serialize(n_elements, buffer);
        for (const auto& [Z, atom] : elements) {
            serializeAtomicElement(atom, buffer);
        }
        return buffer;
    }

    /**
     * @brief Deserializes a map of `AtomicElement` objects from a byte buffer.
     *
     * Reads the element count and reconstructs each `AtomicElement` (via
     * `deserializeAtomicElement`), inserting them into the returned map keyed
     * by atomic number cast to `T`.
     *
     * @tparam T  Key type for the returned map; must be an integral type.
     *            Defaults to `uint8_t` (suitable for Z ≤ 255).
     * @param buffer  Byte buffer produced by `serializeAtoms`.
     * @return Map of `{T(Z) → AtomicElement}`.
     */
    template <std::integral T = std::uint8_t>
    static std::map<T, AtomicElement> deserializeAtoms(std::vector<char>& buffer)
    {
        std::map<T, AtomicElement> elements;
        if (buffer.size() == 0)
            return elements;
        std::uint64_t number_elements;
        auto start = deserialize(number_elements, &(buffer[0]));

        for (std::uint64_t i = 0; i < number_elements; i++) {
            AtomicElement atom;
            start = deserializeAtomicElement(atom, start);
            if constexpr (std::is_same<T, std::uint64_t>::value) {
                elements[atom.Z] = atom;
            } else {
                elements[static_cast<T>(atom.Z)] = atom;
            }
        }
        return elements;
    }

protected:
    /**
     * @brief Appends a single numeric value to @p buffer as raw bytes.
     * @tparam T  Integral or floating-point type (`Number` concept).
     * @param in      Value to serialize.
     * @param buffer  Destination byte buffer.
     */
    template <Number T>
    static void serialize(T in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        auto in_c = reinterpret_cast<char*>(&in);
        std::copy(in_c, in_c + sizeof(T), dest);
    }

    /**
     * @brief Appends the contents of @p data to @p buffer.
     * @param data    Bytes to append.
     * @param buffer  Destination byte buffer.
     */
    static void appendToBuffer(const std::vector<char>& data, std::vector<char>& buffer)
    {
        std::copy(data.cbegin(), data.cend(), std::back_inserter(buffer));
    }

    /**
     * @brief Appends a length-prefixed vector of non-char elements to @p buffer.
     *
     * Writes a `uint64_t` byte count followed by the raw element data.
     *
     * @tparam T  Element type (must not be `char` to avoid ambiguity with raw byte buffers).
     * @param data    Vector to serialize.
     * @param buffer  Destination byte buffer.
     */
    template <typename T>
        requires(!std::is_same<T, char>::value)
    static void serialize(const std::vector<T>& data, std::vector<char>& buffer)
    {
        std::uint64_t size = data.size() * sizeof(T);
        serialize(size, buffer);
        auto in_c = reinterpret_cast<const char*>(data.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }

    /**
     * @brief Serializes an `AtomicShell` as a length-prefixed block into @p buffer.
     *
     * Writes all scalar fields, the photoelectric cross-section table, and each
     * `AtomicShellRadiativeEmission` entry. The entire payload is preceded by a
     * `uint64_t` byte count.
     *
     * @param shell   Shell to serialize.
     * @param buffer  Destination byte buffer.
     */
    static void serializeAtomicShell(const AtomicShell& shell, std::vector<char>& buffer)
    {
        std::vector<char> data;
        serialize(shell.shell, data);
        serialize(shell.numberOfElectrons, data);
        serialize(shell.bindingEnergy, data);
        serialize(shell.kineticEnergy, data);
        serialize(shell.HartreeFockOrbital_0, data);
        serialize(shell.numberOfPhotonsPerInitVacancy, data);
        serialize(shell.energyOfPhotonsPerInitVacancy, data);
        serialize(shell.photoel, data);

        // serialize atomic shell radiative transitions
        const std::uint64_t shellTransSize = shell.radiativeEmissions.size();
        serialize(shell.radiativeEmissions.size(), data);
        for (std::uint64_t i = 0; i < shellTransSize; ++i) {
            serialize(shell.radiativeEmissions[i].vacancy, data);
            serialize(shell.radiativeEmissions[i].probability, data);
            serialize(shell.radiativeEmissions[i].energy, data);
        }

        std::uint64_t data_size = data.size();
        serialize(data_size, buffer); // adding data size to buffer
        appendToBuffer(data, buffer); // adding data
    }

    /**
     * @brief Serializes an `AtomicElement` as a length-prefixed block into @p buffer.
     *
     * Writes Z, atomic weight, density, all cross-section tables, and then each
     * subshell (via `serializeAtomicShell`). The entire payload is preceded by a
     * `uint64_t` byte count.
     *
     * @param atom    Element to serialize.
     * @param buffer  Destination byte buffer.
     */
    static void serializeAtomicElement(const AtomicElement& atom, std::vector<char>& buffer)
    {
        std::vector<char> data;
        serialize(atom.Z, data);
        serialize(atom.atomicWeight, data);
        serialize(atom.standardDensity, data);
        serialize(atom.coherent, data);
        serialize(atom.incoherent, data);
        serialize(atom.photoel, data);
        serialize(atom.formFactor, data);
        serialize(atom.incoherentSF, data);
        serialize(atom.incoherentMeanScatterEnergy, data);

        // adding shells
        // adding number of shells
        std::uint64_t n_shells = atom.shells.size();
        serialize(n_shells, data);
        // adding each shell
        for (const auto& [id, shell] : atom.shells) {
            serializeAtomicShell(shell, data);
        }
        std::uint64_t data_size = data.size();
        serialize(data_size, buffer); // adding data size to buffer
        appendToBuffer(data, buffer); // adding data to buffer
    }

    /**
     * @brief Reads a single numeric value from @p begin and advances the pointer.
     * @tparam T  Integral or floating-point type (`Number` concept).
     * @param value  Output: the deserialized value.
     * @param begin  Pointer to the start of the encoded value.
     * @return Pointer to the byte immediately after the consumed value.
     */
    template <Number T>
    static char* deserialize(T& value, char* begin)
    {
        auto val_ptr = reinterpret_cast<T*>(begin);
        value = *val_ptr;
        return begin + sizeof(T);
    }
    /**
     * @brief Reads @p size bytes from @p begin into @p val and advances the pointer.
     *
     * The number of elements is `size / sizeof(T)`. The vector is cleared and
     * re-populated.
     *
     * @tparam T  Element type.
     * @param val    Output vector.
     * @param begin  Pointer to the start of the raw element data (no size prefix).
     * @param size   Number of bytes to read.
     * @return Pointer to the byte immediately after the consumed data.
     */
    template <typename T>
    static char* deserialize(std::vector<T>& val, char* begin, std::size_t size)
    {
        auto n_elements = size / (sizeof(T));
        val.clear();
        val.reserve(n_elements);
        auto start = reinterpret_cast<T*>(begin);
        std::copy(start, start + n_elements, std::back_inserter(val));
        return begin + n_elements * sizeof(T);
    }

    /**
     * @brief Reads a length-prefixed vector from @p begin and advances the pointer.
     *
     * Reads a leading `uint64_t` byte count, then delegates to the sized overload.
     *
     * @tparam T  Element type.
     * @param val    Output vector.
     * @param begin  Pointer to the `uint64_t` size prefix.
     * @return Pointer to the byte immediately after the consumed data.
     */
    template <typename T>
    static char* deserialize(std::vector<T>& val, char* begin)
    {
        std::uint64_t size;
        auto start = deserialize(size, begin);
        return deserialize(val, start, size);
    }
    /**
     * @brief Deserializes an `AtomicShell` from a length-prefixed block.
     *
     * Reads the block size, then reconstructs all scalar fields, the photoelectric
     * table, and the radiative emission list in the order written by
     * `serializeAtomicShell`.
     *
     * @param shell  Output: populated `AtomicShell`.
     * @param begin  Pointer to the `uint64_t` size prefix of the shell block.
     * @return Pointer to the byte immediately after the consumed block.
     */
    static char* deserializeAtomicShell(AtomicShell& shell, char* begin)
    {
        std::uint64_t size { 0 };
        auto start = deserialize(size, begin);
        start = deserialize(shell.shell, start);
        start = deserialize(shell.numberOfElectrons, start);
        start = deserialize(shell.bindingEnergy, start);
        start = deserialize(shell.kineticEnergy, start);
        start = deserialize(shell.HartreeFockOrbital_0, start);
        start = deserialize(shell.numberOfPhotonsPerInitVacancy, start);
        start = deserialize(shell.energyOfPhotonsPerInitVacancy, start);
        start = deserialize(shell.photoel, start);

        std::uint64_t shellTransSize = 0;
        start = deserialize(shellTransSize, start);
        shell.radiativeEmissions.resize(shellTransSize);
        for (std::uint64_t i = 0; i < shellTransSize; ++i) {
            start = deserialize(shell.radiativeEmissions[i].vacancy, start);
            start = deserialize(shell.radiativeEmissions[i].probability, start);
            start = deserialize(shell.radiativeEmissions[i].energy, start);
        }
        return start;
    }
    /**
     * @brief Deserializes an `AtomicElement` from a length-prefixed block.
     *
     * Reads the block size, then reconstructs Z, atomic weight, density, all
     * cross-section tables, and each subshell (via `deserializeAtomicShell`) in
     * the order written by `serializeAtomicElement`.
     *
     * @param atom   Output: populated `AtomicElement`.
     * @param begin  Pointer to the `uint64_t` size prefix of the element block.
     * @return Pointer to the byte immediately after the consumed block.
     */
    static char* deserializeAtomicElement(AtomicElement& atom, char* begin)
    {
        std::uint64_t size;
        auto start = deserialize(size, begin);
        start = deserialize(atom.Z, start);
        start = deserialize(atom.atomicWeight, start);
        start = deserialize(atom.standardDensity, start);
        start = deserialize(atom.coherent, start);
        start = deserialize(atom.incoherent, start);
        start = deserialize(atom.photoel, start);
        start = deserialize(atom.formFactor, start);
        start = deserialize(atom.incoherentSF, start);
        start = deserialize(atom.incoherentMeanScatterEnergy, start);

        std::uint64_t n_shells { 0 };
        atom.shells.clear();
        start = deserialize(n_shells, start);
        for (std::uint64_t i = 0; i < n_shells; ++i) {
            AtomicShell shell;
            start = deserializeAtomicShell(shell, start);
            atom.shells[shell.shell] = shell;
        }
        return start;
    }
};

}
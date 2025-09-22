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

#include "atomicelementhandler.hpp"

#include <map>
#include <string>
#include <vector>

class EPICSparser {
public:
    EPICSparser(const std::string& path);
    EPICSparser(std::vector<char>& data);
    void read(const std::string& path);
    void readHartreeFockProfiles(const std::string& path);
    void readComptonProfiles(const std::string& path);
    void readStandardDensities(const std::string& path);

    std::map<std::uint64_t, AtomicElementHandler>& getElements() { return m_elements; }

    std::vector<char> serializeElements() const;
    void deserializeElements(std::vector<char>& data);

protected:
    constexpr static std::size_t endIdx() { return 71; }

private:
    std::map<std::uint64_t, AtomicElementHandler> m_elements;
};
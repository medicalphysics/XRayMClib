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

#include "xraymc/serializer.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>

template <typename T>
bool testItem(const T& item)
{

    xraymc::Serializer s;

    auto buffer = xraymc::Serializer::getEmptyBuffer();
    xraymc::Serializer::serializeItem(item, buffer);

    xraymc::Serializer::write("testitem.xr", buffer);

    auto rbuffer = s.read("testitem.xr").value();

    auto name = T::magicID();
    std::vector<char> itemBuffer;
    s.deserializeItem(name, itemBuffer, rbuffer);
    auto item_opt = T::deserialize(itemBuffer);
    if (item_opt) {
        auto& item_read = item_opt.value();
        return item_read == item;
    }

    return false;
}

bool testString()
{
    xraymc::Serializer s;

    std::vector<std::string> ss;
    ss.push_back("Test 1");
    ss.push_back("Testing 2 ");

    auto buffer = s.getEmptyBuffer();
    s.serialize(ss, buffer);

    std::vector<std::string> ss_out;
    auto end = s.deserialize(ss_out, buffer);
    return ss_out == ss;
}

int main()
{

    bool success = testString();

    xraymc::WorldSphere<12, 2, true> sphere;
    sphere.setCenter({ 1, 2, 3 });
    sphere.setMaterialDensity(1.9);
    sphere.setRadius(5);

    success = success && testItem(sphere);

    xraymc::WorldCylinder<16, 2> cylinder;
    success = success && testItem(cylinder);

    xraymc::AAVoxelGrid<12, 2, 255> vgrid;
    success = success && testItem(vgrid);

    xraymc::WorldBoxGrid<15, 2> boxgrid;
    success = success && testItem(boxgrid);

    xraymc::WorldBox<12, 2> box;
    success = success && testItem(box);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
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
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <iostream>
int main()
{

    xraymc::Serializer s;
    auto buffer = s.getEmptyBuffer();

    xraymc::WorldSphere<12, 2, true> sphere;
    sphere.setCenter({ 1, 2, 3 });
    sphere.setMaterialDensity(1.9);
    sphere.setRadius(5);

    sphere.serialize(buffer);

    s.write("testsphere.xr", buffer);

    auto rbuffer = s.read("testsphere.xr").value();
    std::span<const char> span(rbuffer);

    auto sphere_opt = xraymc::WorldSphere<12, 2, true>::deserialize(span);

    std::cout << span.size() << std::endl;


    return EXIT_SUCCESS;
}
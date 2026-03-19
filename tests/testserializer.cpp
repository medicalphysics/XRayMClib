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

int main()
{

    xraymc::Serializer s;
    auto buffer = s.getEmptyBuffer();

    double in = 1.0;
    s.serialize(in, buffer);

    std::vector<double> vec_in(5, 1.0);
    s.serialize(vec_in, buffer);

    double in2 = 2;
    s.serialize(in2, buffer);

    s.write("test.xray", buffer);

    auto t = s.read("test.xray");

    double out;
    auto start = s.deserialize(out, buffer.data());
    std::vector<double> vec_out;
    start = s.deserialize(vec_out, start);
    double out2;
    start = s.deserialize(out2, start);

    return EXIT_SUCCESS;
}
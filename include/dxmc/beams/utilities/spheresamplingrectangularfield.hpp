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

Copyright 2025 Erlend Andersen
*/

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {
struct SphereSamplingRectangularField {
    // Sampling of uniform points on a sphere constrained by a rectangular field
    // Based on random sampling of points on a whole sphere by sampling x, y, z in [-1, 1] and normalizing the vector
    // Uses simple random sampling of x, y, z inside a constrained box and rejecting points outside field
    // perhaps https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone is better?

    SphereSamplingRectangularField(const std::array<double, 2>& angles)
    {
        setData(angles[0], angles[1]);
    }
    SphereSamplingRectangularField(double collimationHalfAngle_x, double collimationHalfAngle_y)
    {
        setData(collimationHalfAngle_x, collimationHalfAngle_y);
    }

    void setData(double collimationHalfAngle_x, double collimationHalfAngle_y)
    {
        borderx = std::sin(collimationHalfAngle_x);
        bordery = std::sin(collimationHalfAngle_y);
        limx = std::tan(collimationHalfAngle_x);
        limy = std::tan(collimationHalfAngle_y);
        const auto sinz = std::sqrt(limx * limx + limy * limy);
        limz = std::sqrt(1 - sinz * sinz);
    }

    std::array<double, 3> operator()(RandomState& state) const
    {
        std::array<double, 3> dir;
        do {
            dir = {
                state.randomUniform(-limx, limx),
                state.randomUniform(-limy, limy),
                state.randomUniform(limz, 1.0)
            };
            vectormath::normalize(dir);
        } while (std::abs(dir[0]) > borderx || std::abs(dir[1]) > bordery);
        return dir;

        // safe naive version
        // const auto sinz = std::sqrt(borderx * borderx + bordery * bordery);
        // const auto cosz = std::sqrt(1 - sinz * sinz);
        // std::array<double, 3> dir;
        // do {
        //     dir = {
        //         state.randomUniform(-sinz, sinz),
        //         state.randomUniform(-sinz, sinz),
        //         state.randomUniform(cosz, 1.0)
        //     };
        //     vectormath::normalize(dir);
        // } while (std::abs(dir[0]) > borderx || std::abs(dir[1]) > bordery);
        // return dir;
    }

    double borderx = 1;
    double bordery = 1;
    double limx = 0;
    double limy = 0;
    double limz = 1;
};

struct SphereSamplingRectangularFieldOffset {
    // Sampling of uniform points on a sphere constrained by a rectangular field
    // Based on random sampling of points on a whole sphere by sampling x, y, z in [-1, 1] and normalizing the vector
    // Uses simple random sampling of x, y, z inside a constrained box and rejecting points outside field
    // perhaps https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone is better?
    SphereSamplingRectangularFieldOffset() { }
    SphereSamplingRectangularFieldOffset(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }
    SphereSamplingRectangularFieldOffset(double xmin, double ymin, double xmax, double ymax)
    {
        setData(xmin, ymin, xmax, ymax);
    }
    
    void setData(const std::array<double, 4>& angles)
    {
        setData(angles[0], angles[1], angles[2], angles[3]);
    }

    void setData(double xmin, double ymin, double xmax, double ymax)
    {
        borderx_min = std::sin(xmin);
        borderx_max = std::sin(xmax);
        bordery_min = std::sin(ymin);
        bordery_max = std::sin(ymax);

        limx_min = std::tan(xmin);
        limx_max = std::tan(xmax);
        limy_min = std::tan(ymin);
        limy_max = std::tan(ymax);
        const auto sinz = std::sqrt(std::max(limx_min * limx_min, limx_max * limx_max) + std::max(limy_min * limy_min, limy_max * limy_max));
        limz = std::sqrt(1 - sinz * sinz);
    }

    std::array<double, 3> operator()(RandomState& state) const
    {
        std::array<double, 3> dir;
        do {
            dir = {
                state.randomUniform(limx_min, limx_max),
                state.randomUniform(limy_min, limy_max),
                state.randomUniform(limz, 1.0)
            };
            vectormath::normalize(dir);
        } while (dir[0] < borderx_min || dir[0] > borderx_max || dir[1] < bordery_min || dir[1] > bordery_max);
        return dir;
    }

    double borderx_min = -1;
    double borderx_max = 1;
    double bordery_min = -1;
    double bordery_max = 1;
    double limx_min = 0;
    double limx_max = 0;
    double limy_min = 0;
    double limy_max = 0;
    double limz = 1;
};
}
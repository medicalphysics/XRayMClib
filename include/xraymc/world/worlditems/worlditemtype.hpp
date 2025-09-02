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
#include "xraymc/particle.hpp"
#include "xraymc/world/dosescore.hpp"
#include "xraymc/world/energyscore.hpp"
#include "xraymc/world/visualizationintersectionresult.hpp"
#include "xraymc/world/worldintersectionresult.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <optional>

namespace xraymc {

template <typename U>
concept WorldItemType = requires(U u, Particle p, std::array<double, 3> vec, std::size_t index, double factor, RandomState state) {
    u.translate(vec);
    {
        u.center()
    } -> std::convertible_to<std::array<double, 3>>;
    {
        u.AABB()
    } -> std::convertible_to<std::array<double, 6>>;
    {
        u.intersect(p)
    } -> std::same_as<WorldIntersectionResult>;
    {
        u.energyScored(index)
    } -> std::convertible_to<EnergyScore>;
    {
        u.doseScored(index)
    } -> std::convertible_to<DoseScore>;
    u.clearDoseScored();
    u.clearEnergyScored();
    u.addEnergyScoredToDoseScore(factor);
    u.transport(p, state);
};
}
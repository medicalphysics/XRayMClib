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
#include "xraymc/transportprogress.hpp"
#include "xraymc/xraymcrandom.hpp"

#include <array>

namespace xraymc {

template <typename B>
concept BeamType = requires(B beam, std::array<double, 3> vec, RandomState state, std::uint64_t index, TransportProgress* progress) {
    {
        beam.exposure(index).sampleParticle(state)
    }; // -> std::same_as<Particle>;

    {
        beam.exposure(index).numberOfParticles()
    } -> std::same_as<std::uint64_t>;

    {
        beam.exposure(index).position()
    } -> std::convertible_to<std::array<double, 3>>;

    {
        beam.numberOfExposures()
    } -> std::same_as<std::uint64_t>;

    {
        beam.numberOfParticles()
    } -> std::same_as<std::uint64_t>;

    {
        beam.calibrationFactor(progress)
    } -> std::convertible_to<double>;
};
}
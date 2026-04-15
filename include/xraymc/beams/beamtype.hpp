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

/**
 * @brief Concept constraining a type to the beam interface required by the transport engine.
 *
 * A type `B` satisfies `BeamType` when it provides the following expressions:
 *
 * | Expression | Return type | Description |
 * |---|---|---|
 * | `beam.exposure(index).sampleParticle(state)` | (any) | Samples a photon from exposure `index` using PRNG `state`. |
 * | `beam.exposure(index).numberOfParticles()` | `std::uint64_t` | Number of particles in exposure `index`. |
 * | `beam.exposure(index).position()` | convertible to `std::array<double,3>` | Source position for exposure `index` [cm]. |
 * | `beam.numberOfExposures()` | `std::uint64_t` | Total number of exposures (i.e. gantry positions / tube firings). |
 * | `beam.numberOfParticles()` | `std::uint64_t` | Total number of particles across all exposures. |
 * | `beam.calibrationFactor(progress)` | convertible to `double` | Dose calibration factor applied after transport. |
 *
 * The `Transport` engine distributes exposures across worker threads and calls
 * `beam.calibrationFactor(progress)` once after all threads complete to convert
 * accumulated energy scores to dose.
 *
 * Conforming beam types include (non-exhaustively): `CTSpiralBeam`, `CTAxialBeam`,
 * `DXBeam`, and `IsotropicBeam`.
 *
 * @tparam B  The beam type to check.
 */
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
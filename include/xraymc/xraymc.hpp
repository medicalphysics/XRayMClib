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

/**
 * @file xraymc.hpp
 * @brief Convenience umbrella header for XRayMClib.
 *
 * Including this single header pulls in the entire public API of XRayMClib,
 * an X-ray Monte Carlo photon transport library for medical physics applications.
 *
 * ## Library overview
 *
 * XRayMClib simulates the transport of diagnostic X-ray photons through geometric
 * worlds composed of heterogeneous materials. The library is structured around four
 * main concepts:
 *
 * ### Beams
 * Source definitions that sample photon energies, positions, and directions:
 * - **CT beams**: `CTSpiralBeam`, `CTSequentialBeam`, `CTSpiralDualEnergyBeam`, `CBCTBeam`
 * - **Radiography**: `DXBeam`
 * - **Generic**: `IsotropicBeam`, `IsotropicBeamCircle`, `IsotropicMonoEnergyBeam`,
 *                `IsotropicMonoEnergyBeamCircle`, `PencilBeam`
 *
 * ### World and world items
 * The `World` class holds a heterogeneous collection of geometric items accelerated
 * by an internal KD-tree. Available item types include:
 * - **Primitive shapes**: `WorldBox`, `WorldSphere`, `WorldCylinder`, `WorldBoxGrid`
 * - **Mesh-based geometry**: `TetrahedralMesh` (with file reader), `TriangulatedMesh`
 *   (STL reader), `TriangulatedOpenSurface`
 * - **Scoring volumes**: `FlatDetector`, `DepthDose`, `FluenceScore`, `AAVoxelGrid`
 * - **Phantoms**: `CTDIPhantom`, `EnclosedRoom`
 *
 * ### Materials
 * `Material` describes a medium by its elemental composition and density.
 * Pre-built NIST reference compositions are available via `NistMaterials`.
 *
 * ### Transport
 * `Transport` drives multi-threaded Monte Carlo photon transport, distributing
 * exposures across worker threads and applying beam calibration factors to convert
 * accumulated energy scores to dose.
 *
 * ## Typical usage
 * @code
 * #include <xraymc/xraymc.hpp>
 *
 * xraymc::World world;
 * world.addItem(xraymc::WorldBox{ ... });
 * world.build();
 *
 * xraymc::CTSpiralBeam beam{ ... };
 * xraymc::Transport::runConsole(world, beam);
 * @endcode
 */
#pragma once

#include "xraymc/beams/cbctbeam.hpp"
#include "xraymc/beams/ctsequentialbeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/ctspiraldualenergybeam.hpp"
#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/isotropicbeamcircle.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/material/material.hpp"
#include "xraymc/material/nistmaterials.hpp"
#include "xraymc/transport.hpp"
#include "xraymc/world/visualization/visualizeworld.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/aavoxelgrid.hpp"
#include "xraymc/world/worlditems/ctdiphantom.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/flatdetector.hpp"
#include "xraymc/world/worlditems/fluencescore.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh.hpp"
#include "xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "xraymc/world/worlditems/triangulatedmesh.hpp"
#include "xraymc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"
#include "xraymc/world/worlditems/triangulatedopensurface.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldboxgrid.hpp"
#include "xraymc/world/worlditems/worldcylinder.hpp"
#include "xraymc/world/worlditems/worlditemtype.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

/**
 * @namespace xraymc
 * @brief Primary namespace for all XRayMClib functionality.
 *
 * Every public class, function, and type in XRayMClib lives in this namespace.
 * Internal implementation details are nested under sub-namespaces such as
 * `xraymc::basicshape`, `xraymc::vectormath`, and `xraymc::interactions`.
 *
 * The top-level namespace exposes:
 * - Beam types (e.g. `CTSpiralBeam`, `DXBeam`, `PencilBeam`)
 * - `World` and the full set of world item types
 * - `Material` and `NistMaterials`
 * - `Transport` for running simulations
 * - Visualization utilities via `visualizeWorld()`
 */
namespace xraymc {
}

# XrayMClib

A C++ library for Monte Carlo simulation of x-ray dose scoring, supporting voxel geometries, triangulated and tetrahedral mesh geometries, and basic shapes (spheres, boxes, cylinders). Designed for energy levels used in diagnostic radiology. 

XrayMC is benchmarked against EGSnrc, Geant4, MCNP, and PENELOPE in the research paper [Development and validation of XrayMC: a dedicated Monte Carlo tool for X-ray imaging and radiation protection](https://iopscience.iop.org/article/10.1088/1361-6560/ae8353).

XrayMClib is the primary simulation engine of [OpenXrayMC](https://github.com/medicalphysics/Openxraymc), a GUI application for Monte Carlo patient dose simulation of CT scans, conventional x-rays, and CBCT scans.

---

## Features

- **Dose scoring** in voxel, mesh, and basic shape geometries
- **X-ray spectrum generation** — semi-analytical tungsten-anode model (50–150 kV) ([Poludniowski & Evans](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725))
- **Multiple beam types**: radiography, CT sequential/spiral/dual-energy/cone-beam and pencil beam
- **Three binding-energy correction levels**: None, Livermore, or Impulse Approximation
- **Variance reduction**: forced photoelectric interactions for low-density scoring volumes (e.g. air cavities in CTDI phantoms)
- **Material library**: define materials by atomic number, chemical formula (e.g. `H2O`), or NIST name; 180 NIST materials included
- **Multithreaded transport** via exposure-level parallelism using `std::atomic_ref` (C++20)
- **Built-in visualization**: render geometry and dose distributions to PNG

## Physics overview

XrayMC uses Woodcock (delta) tracking for photon transport in voxelized volumes and accelerated Siddons path tracing for tetrahedral mesh traversal. Three binding energy corrections can be set at compile time, they are:

| Interaction | None | Livermore | Impulse Approximation (IA) |
|---|---|---|---|
| Photoelectric | Full energy deposit | Same as None | Characteristic x-ray emission from selected shell |
| Compton | Klein–Nishina (free electron) | + Livermore scatter function correction | Relativistic impulse approximation sampling|
| Rayleigh | Thomson free-electron | + Hubbell atomic form factor | Same as Livermore |

Secondary electrons are assumed to deposit their energy locally (valid for photon energies ≤ 150 keV).


## Dependencies

| Dependency | Purpose |
|---|---|
| [EPICS 2025](https://www-nds.iaea.org/epics/) | Cross sections, atomic form factors, electron shell binding energies, atomic data |
| CMake ≥ 3.20 | Build system |

**Compiler requirements:** MSVC ≥ 16.8 or Clang ≥ 13.0 or GCC ≥ 14 (C++20 required, C++23 recommended).

## Installation

XrayMClib uses CMake. The recommended approach is `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
    libxraymc
    GIT_REPOSITORY https://github.com/medicalphysics/xraymclib.git
)
FetchContent_MakeAvailable(libxraymc)

add_executable(your_executable your_example.cpp)
target_include_directories(your_executable PRIVATE ${libxraymc_SOURCE_DIR}/include)
target_link_libraries(your_executable PRIVATE libxraymc)

# Copies physicslist.bin to your executable folder at build time
xraymclib_add_physics_list(your_executable)
```

### EPICS data

Set one of these CMake variables during configure:

- `XRAYMCLIB_EPICS_DATA_DIRPATH` — path to a folder containing `EPDL2025.ALL` and `EADL2025.ALL` in ENDL format
- `XRAYMCLIB_EPICS_DOWNLOAD=ON` — download EPICS data from IAEA automatically during configure step

## CMake build options

| Option | Default | Description |
|---|---|---|
| `XRAYMCLIB_EPICS_DOWNLOAD` | `OFF` | Download EPICS data during configure |
| `XRAYMCLIB_EPICS_DATA_DIRPATH` | _(empty)_ | Path to local EPICS data files |
| `XRAYMCLIB_USE_LOADPNG` | `ON` | Use the bundled lodepng library for writing PNG files |
| `XRAYMCLIB_USE_ZLIB` | `ON` | Use zlib, if found, for buffer compression |
| `XRAYMCLIB_BUILD_TESTS` | `ON` when top level | Build the test suite |
| `XRAYMCLIB_BUILD_EXAMPLES` | `ON` | Build the examples (top level only) |

### Compiler and optimization options

| Option | Default | Description |
|---|---|---|
| `XRAYMCLIB_USE_FAST_MATH` | `ON` | Fast floating point math, see below |
| `XRAYMCLIB_OPTIMIZE_FOR_NATIVE` | `ON` | Build for the host CPU: `-march=native` (or `-mcpu=native`) on GCC/Clang, `/arch:<XRAYMCLIB_MSVC_NATIVE_OPTION>` on MSVC |
| `XRAYMCLIB_MSVC_NATIVE_OPTION` | `AVX2` | MSVC instruction set: `SSE2`, `AVX`, `AVX2` or `AVX512` |
| `XRAYMCLIB_WHOLE_PROGRAM_OPTIMIZATION` | `ON` | Link time optimization (`-flto` on GCC/Clang, `/GL` + `/LTCG` on MSVC), if the toolchain supports it |
| `XRAYMCLIB_PGO` | `OFF` | Profile guided optimization: `OFF`, `GENERATE` or `USE`, see below |
| `XRAYMCLIB_PGO_DIR` | `<build>/pgo` | Directory for profile data |

The build type defaults to `Release` when none is given. Fast math, link time optimization, profile guided optimization and the MSVC `/Ob3 /Gw` flags apply only to `Release` and `RelWithDebInfo`. The native architecture flag applies to all configurations.

These flags are also passed on to your targets when you link `libxraymc` through `FetchContent` or `add_subdirectory`, since the header only code is compiled as part of your program. They are not part of an installed package.

**Fast math.** On GCC and Clang (including clang-cl), the individual parts of `-ffast-math` are enabled: `-ffinite-math-only -fno-math-errno -fno-trapping-math -fno-signed-zeros -fassociative-math -freciprocal-math -ffp-contract=fast`. On MSVC, `/fp:fast` is used. With these flags the compiler assumes that inf and NaN do not occur.

**Binaries built for the host CPU** may not run on other machines. Turn `XRAYMCLIB_OPTIMIZE_FOR_NATIVE` off, or pick a lower `XRAYMCLIB_MSVC_NATIVE_OPTION`, when distributing binaries.

**Profile guided optimization** is a two step build in the same build directory. It is supported on GCC, Clang and MSVC, but not on clang-cl:

```bash
cmake -B build -DXRAYMCLIB_PGO=GENERATE
cmake --build build
cmake --build build --target pgo-train   # runs the test suite to collect profiles
cmake -B build -DXRAYMCLIB_PGO=USE
cmake --build build
```

`pgo-train` needs `llvm-profdata` with Clang to merge the profiles. `pgo-clean` removes collected profiles, which otherwise add up over training runs. Each executable is optimized from its own profile, so only programs that ran during training benefit. When your own project consumes an instrumented build, it is instrumented as well.

**MSVC C runtime.** To link the C runtime statically (`/MT`), configure with `-DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>"`. Every library linked into the same executable must use the same runtime.

**Parallel algorithms.** With libstdc++ (GCC, and usually Clang on Linux), `std::execution::par` runs in parallel only when TBB is available. TBB is linked automatically when CMake finds it, otherwise these algorithms run serially.

## Testing

Tests are built when XrayMClib is the top level project (or with `XRAYMCLIB_BUILD_TESTS=ON`). Build all test executables with the `xraymclib_tests` target and run them from the build directory with CTest:

```bash
cmake --build build --config Release --target xraymclib_tests
cd build
ctest -C Release            # all tests
ctest -C Release -L unit    # fast tests of math, geometry, materials, beams and world items
ctest -C Release -L transport  # Monte Carlo transport tests compared against analytical results
```

Monte Carlo tests compare results statistically, build in Release mode since Debug builds are slow.

## Quick start examples

Examples can be explored at 

- [XrayMC example 1](https://github.com/medicalphysics/xraymcExample1): Depth dose simulation of a pencil beam onto a slab of aluminum
- [XrayMC example 2](https://github.com/medicalphysics/xraymcExample2): Simple beam onto a tetrahedral mesh water sphere with an internal aluminum tourus.
- [XrayMC example 3](https://github.com/medicalphysics/xraymcExample3): Simulation of scattered radiation in an interventional X-ray image guided laboratory. 
- [XrayMC example 4](https://github.com/medicalphysics/xraymcExample4): Minimal example on how to implement a spherical shell as a user defined world object in XrayMC. 



## Documentation

Full API reference and physics model description: [xraymclib.readthedocs.io](https://xraymclib.readthedocs.io/en/latest/)

## References

- G. Poludniowski & P. Evans, [Calculation of x-ray spectra emerging from an x-ray tube. Part I](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725), *Med. Phys.* 34, 2164–2174 (2007)
- G. Poludniowski, [Calculation of x-ray spectra emerging from an x-ray tube. Part II](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734726), *Med. Phys.* 34, 2175–2186 (2007)
- T. Schoonjans et al., [xraylib](https://github.com/tschoonj/xraylib)
- J.H. Hubbell et al., *Atomic form factors, incoherent scattering functions, and photon scattering cross sections*, J. Phys. Chem. Ref. Data 4(3) (1975)
- I. Kawrakow et al., [The EGSnrc Code System](https://nrc-cnrc.github.io/EGSnrc/), Technical Report PIRS-701, NRC Canada (2020)
- [EPICS 2025 database](https://www-nds.iaea.org/epics/), IAEA Nuclear Data Section. In the future the EPICS library will be available from the [Lawrence Livermore National Laboratory](https://nuclear.llnl.gov/EPICS/)

## How to cite

If you use XRayMClib in your research, please cite:

> Andersen E K F, Østerås B H, Samnøy A T and Ytre-Hauge K S 2026 Development and validation of XrayMC: a dedicated Monte Carlo tool for X-ray imaging and radiation protection *Phys. Med. Biol.* https://doi.org/10.1088/1361-6560/ae8353



## License

GnuGPLv3, see [LICENSE](LICENSE).

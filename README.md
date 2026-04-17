# XRayMClib

A C++ library for Monte Carlo simulation of x-ray dose scoring, supporting voxel geometries, triangulated and tetrahedral mesh geometries, and basic shapes (spheres, boxes, cylinders). Designed for energy levels used in diagnostic radiology.

XRayMClib is the primary simulation engine of [OpenXRayMC](https://github.com/medicalphysics/Openxraymc), a GUI application for Monte Carlo patient dose simulation of CT scans, conventional x-rays, and CBCT scans.

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

XRayMClib uses Woodcock (delta) tracking for photon transport in voxelized volumes and accelerated Siddons path tracing for tetrahedral mesh traversal. The three binding energy corrections are:

| Interaction | None | Livermore | Impulse Approx. |
|---|---|---|---|
| Photoelectric | Full energy deposit | Same as None | Characteristic x-ray emission from selected shell |
| Compton | Klein–Nishina (free electron) | + Hartree–Fock scatter function correction | + Shell-specific momentum sampling (Doppler broadening) |
| Rayleigh | Thomson free-electron | + Hubbell atomic form factor | Same as Livermore |

Secondary electrons are assumed to deposit their energy locally (valid for diagnostic energies ≤ 150 keV).


## Dependencies

| Dependency | Purpose |
|---|---|
| [EPICS 2025](https://www-nds.iaea.org/epics/) | Cross sections, atomic form factors, electron shell binding energies, atomic data |
| [xraylib](https://github.com/tschoonj/xraylib) | Hartree–Fock Compton profiles from the [DABAX](https://www.esrf.fr/Instrumentation/software/data-analysis/Resources) library |
| CMake ≥ 3.20 | Build system |

**Compiler requirements:** MSVC ≥ 16.8 or Clang ≥ 13.0 (C++20 required, C++23 recommended).

## Installation

XRayMClib uses CMake. The recommended approach is `FetchContent`:

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
- `XRAYMCLIB_EPICS_DOWNLOAD=ON` — download EPICS data from IAEA automatically during configure

## CMake build options

| Option | Default | Description |
|---|---|---|
| `XRAYMCLIB_USE_FAST_MATH` | `ON` | Enable `-ffast-math` / `/fp:fast` |
| `XRAYMCLIB_OPTIMIZE_FOR_NATIVE` | `ON` | Enable `-march=native` / AVX on MSVC |
| `XRAYMCLIB_EPICS_DOWNLOAD` | `OFF` | Download EPICS data during configure |
| `XRAYMCLIB_EPICS_DATA_DIRPATH` | _(empty)_ | Path to local EPICS data files |

## Quick start examples

Examples can be explored at 

- [XrayMC example 1](github.com/medicalphysics/xraymcExample1): Depth dose simulation of a pencil beam onto a slab of aluminum
- [XrayMC example 2](github.com/medicalphysics/xraymcExample2): Simple beam onto a tetrahedral mesh water sphere with an internal aluminum tourus.
- [XrayMC example 3](github.com/medicalphysics/xraymcExample3): Simulation of scattered radiation in an interventional X-ray image guided laboratory. 



## Documentation

Full API reference and physics model description: [xraymclib.readthedocs.io](https://xraymclib.readthedocs.io/en/latest/)

## References

- G. Poludniowski & P. Evans, [Calculation of x-ray spectra emerging from an x-ray tube. Part I](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725), *Med. Phys.* 34, 2164–2174 (2007)
- G. Poludniowski, [Calculation of x-ray spectra emerging from an x-ray tube. Part II](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734726), *Med. Phys.* 34, 2175–2186 (2007)
- T. Schoonjans et al., [xraylib](https://github.com/tschoonj/xraylib)
- European Synchrotron Radiation Facility, [DABAX: A Database of Atomic X-ray Parameters](https://www.esrf.fr/Instrumentation/software/data-analysis/Resources), ESRF
- J.H. Hubbell et al., *Atomic form factors, incoherent scattering functions, and photon scattering cross sections*, J. Phys. Chem. Ref. Data 4(3) (1975)
- I. Kawrakow et al., [The EGSnrc Code System](https://nrc-cnrc.github.io/EGSnrc/), Technical Report PIRS-701, NRC Canada (2020)
- [EPICS 2025 database](https://www-nds.iaea.org/epics/), IAEA Nuclear Data Section

## License

See [LICENSE](LICENSE).

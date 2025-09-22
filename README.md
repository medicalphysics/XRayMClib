# xraymclib
C++ library for x-ray dose scoring in voxel and triangulated mesh geometries in addition to basic shapes like spheres and boxes. 

xraymclib aims to be an easy to use C++ dose scoring library for energy levels used in diagnostic radiology. It is the primary simulation engine of [OpenXRayMC](https://github.com/medicalphysics/Openxraymc), a GUI application for Monte Carlo dose simulation of CT scans, conventional x-rays and CBCT scans. XRayMClib also includes a x-ray specter generator based on the work by Gavin Poludniowski and Phil Evans; [Calculation of x‐ray spectra emerging from an x‐ray tube. Part I. Electron penetration characteristics in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725) and [Calculation of x‐ray spectra emerging from an x‐ray tube. Part II. X‐ray production and filtration in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734726).

XRayMClib is dependent on the [EPICS 2025](https://www-nds.iaea.org/epics/) dataset for interaction cross sections, atomic form factors and electron shell binding energies. Harthree-Fock orbital profiles is obtained from [xraylib library by Tom Schoonjans](https://github.com/tschoonj/xraylib). 

### Compilation
XRayMClib uses CMake as build generator, to include XRayMClib in a CMake project it is recommended to use CMakes 'FetchContent' module. Example to include XRayMClib in your CMakeLists.txt:

    include(FetchContent)
    ## Adding xraymclib package
    FetchContent_Declare(
        libxraymc
        GIT_REPOSITORY https://github.com/medicalphysics/xraymclib.git
        )
    FetchContent_MakeAvailable(libxraymc)

    # Example target you develop
    add_executable(your_executable your_example.cpp)
    # Adding xraymclib headers
    target_include_directories(your_executable PRIVATE ${libxraymc_SOURCE_DIR}/include)
    # Linking to xraymclib
    target_link_libraries(your_executable PRIVATE libxraymc)
    # Make atomic data and attenuation data tables available at runtime for 
    # your executable
    # This cmake function vil copy physicslist.bin to your executable folder
    xraymclib_add_physics_list(your_executable)

xraymclib takes advantage of concepts and std::atomic_ref introduced in C++20. Currently this library is tested on MSVC >= 16.8 and Clang >= 13.0 compilers. Set cmake variable XRAYMCLIB_EPICS_DATA_DIRPATH to folder containing EPDL2023.ALL and EADL2023.ALL in ENDL format or set XRAYMCLIB_EPICS_DOWNLOAD to true for download EPICS data from IAEA during the configure step.
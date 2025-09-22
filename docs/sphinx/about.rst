About
-----

xraymclib aims to be an easy to use C++ dose scoring library for voxelized geometry. It is the primary simulation engine of [Openxraymc](https://github.com/medicalphysics/Openxraymc), a GUI application for Monte Carlo dose simulation for diagnostic beam energy levels.

It is possible to simulate dose from conventional x-ray and CT examinations in arbitrary materials. xraymclib also includes a x-ray specter generator based on the work by Gavin Poludniowski and Phil Evans; [Calculation of x‐ray spectra emerging from an x‐ray tube. Part I. Electron penetration characteristics in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725) and [Calculation of x‐ray spectra emerging from an x‐ray tube. Part II. X‐ray production and filtration in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734726).


Compilation
_________
xraymclib uses CMake as build generator, to include xraymclib in a CMake project it is recommended to use CMakes 'FetchContent' module. Example to include xraymclib in your CMakeLists.txt:

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
    # CMake function to copy material data to executable folder
    xraymclib_add_physics_list(your_executable)


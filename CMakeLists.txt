# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.20)

project(xraymclib VERSION 0.7.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Single-config generators (Ninja, Makefiles) otherwise build without any optimization
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, defaulting to Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(CheckCXXCompilerFlag)
include(CheckIPOSupported)


# Optimization flags only apply to optimized configurations
set(XRAYMCLIB_OPT_CONFIG "$<CONFIG:Release,RelWithDebInfo>")

if(MSVC) # cl and clang-cl
    set(XRAYMCLIB_GNU_LIKE OFF)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(XRAYMCLIB_GNU_LIKE ON)
endif()

#MSVC fix for heavy templated code
if(MSVC)
    add_compile_options(/bigobj)
endif()

# /MP: compile source files in parallel within one cl invocation. Only has an effect with the
# Visual Studio generators, Ninja and Makefiles already run one compiler process per file.
if(MSVC AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(/MP)
endif()

# Static C runtime (/MT) instead of the DLL runtime (/MD): calls into the CRT (exp, log, sin, ...)
# become direct calls instead of going through DLL import stubs. Every library linked into the
# same executable must use the same runtime, including libraries built by this project (such as
# lodepng) when it is used as a subproject.
option(XRAYMCLIB_MSVC_STATIC_RUNTIME "Link the MSVC C runtime statically (/MT)" OFF)
if(MSVC AND XRAYMCLIB_MSVC_STATIC_RUNTIME)
    message(STATUS "Building with static MSVC runtime")
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

option(XRAYMCLIB_USE_FAST_MATH "Build with fast floating point math" ON)
if(XRAYMCLIB_USE_FAST_MATH)
    # The individual parts of -ffast-math, except -ffinite-math-only: some code relies on
    # inf/NaN propagating correctly (clamping of inf, NaN-aware comparisons).
    #   -fno-math-errno    : inline sqrt, allows vectorized log/exp/sin/cos
    #   -fno-trapping-math : FP operations are assumed not to trap
    #   -fno-signed-zeros  : -0.0 and +0.0 may be treated alike (required by -fassociative-math)
    #   -fassociative-math : reordering of sums, enables vectorized reductions
    #   -freciprocal-math  : x / y may become x * (1 / y)
    #   -ffp-contract=fast : fuse multiply-add into FMA instructions
    set(XRAYMCLIB_FAST_MATH_GNU_FLAGS
        -fno-math-errno
        -fno-trapping-math
        -fno-signed-zeros
        -fassociative-math
        -freciprocal-math
        -ffp-contract=fast
    )
    if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # clang-cl maps /fp:fast to full -ffast-math including -ffinite-math-only, pass the parts instead
        list(TRANSFORM XRAYMCLIB_FAST_MATH_GNU_FLAGS PREPEND "/clang:" OUTPUT_VARIABLE XRAYMCLIB_FAST_MATH_FLAGS)
    elseif(MSVC)
        set(XRAYMCLIB_FAST_MATH_FLAGS /fp:fast)
    elseif(XRAYMCLIB_GNU_LIKE)
        set(XRAYMCLIB_FAST_MATH_FLAGS ${XRAYMCLIB_FAST_MATH_GNU_FLAGS})
    endif()
    if(XRAYMCLIB_FAST_MATH_FLAGS)
        message(STATUS "Building with fast math flags: ${XRAYMCLIB_FAST_MATH_FLAGS}")
        add_compile_options("$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_FAST_MATH_FLAGS}>")
    endif()
endif()

option(XRAYMCLIB_OPTIMIZE_FOR_NATIVE "Build for the host CPU (-march=native or /arch:X on MSVC), binaries may not run on other machines" ON)
if(XRAYMCLIB_OPTIMIZE_FOR_NATIVE)
    if(MSVC)
        set(XRAYMCLIB_MSVC_NATIVE_OPTION "AVX2" CACHE STRING "Machine instructions support")
        set_property(CACHE XRAYMCLIB_MSVC_NATIVE_OPTION PROPERTY STRINGS SSE2 AVX AVX2 AVX512)
        message(STATUS "Building with /arch:${XRAYMCLIB_MSVC_NATIVE_OPTION} flag")
        add_compile_options("/arch:${XRAYMCLIB_MSVC_NATIVE_OPTION}")
    elseif(XRAYMCLIB_GNU_LIKE)
        check_cxx_compiler_flag("-march=native" XRAYMCLIB_HAS_MARCH_NATIVE)
        check_cxx_compiler_flag("-mcpu=native" XRAYMCLIB_HAS_MCPU_NATIVE) # older Apple clang on arm64
        if(XRAYMCLIB_HAS_MARCH_NATIVE)
            message(STATUS "Building with -march=native flag")
            add_compile_options(-march=native)
        elseif(XRAYMCLIB_HAS_MCPU_NATIVE)
            message(STATUS "Building with -mcpu=native flag")
            add_compile_options(-mcpu=native)
        endif()
    endif()
endif()

if(MSVC)
    # /Ob3: more aggressive inlining than the Release default /Ob2
    # /Gw: put global data in separate sections so the linker can drop unused data
    add_compile_options("$<${XRAYMCLIB_OPT_CONFIG}:/Ob3;/Gw>")
endif()

# Link time optimization: -flto on GCC/Clang, /GL + /LTCG on MSVC. 
option(XRAYMCLIB_WHOLE_PROGRAM_OPTIMALIZATION "Build with link time optimization" OFF)
if(XRAYMCLIB_WHOLE_PROGRAM_OPTIMALIZATION)
    check_ipo_supported(RESULT XRAYMCLIB_IPO_SUPPORTED OUTPUT XRAYMCLIB_IPO_ERROR)
    if(XRAYMCLIB_IPO_SUPPORTED)
        message(STATUS "Building with link time optimization")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
        if(MSVC)
            add_link_options("$<${XRAYMCLIB_OPT_CONFIG}:/CGTHREADS:8>")
        endif()
    else()
        message(WARNING "Link time optimization not supported: ${XRAYMCLIB_IPO_ERROR}")
    endif()
endif()

# Profile guided optimization, a two step build in the same build directory:
#   1. configure with -DXRAYMCLIB_PGO=GENERATE, build, then run the programs on representative
#      workloads (the pgo-train target runs the test suite)
#   2. reconfigure with -DXRAYMCLIB_PGO=USE and rebuild
# Each executable is optimized from its own profile, so only programs that were run benefit.
set(XRAYMCLIB_PGO "OFF" CACHE STRING "Profile guided optimization: OFF, GENERATE (instrumented build) or USE (optimize from profile)")
set_property(CACHE XRAYMCLIB_PGO PROPERTY STRINGS OFF GENERATE USE)
set(XRAYMCLIB_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "Directory for profile guided optimization data")
if(NOT XRAYMCLIB_PGO STREQUAL "OFF")
    if(NOT XRAYMCLIB_PGO MATCHES "^(GENERATE|USE)$")
        message(FATAL_ERROR "XRAYMCLIB_PGO must be OFF, GENERATE or USE, not '${XRAYMCLIB_PGO}'")
    endif()
    set(XRAYMCLIB_PGO_COMPILE_FLAGS "")
    set(XRAYMCLIB_PGO_LINK_FLAGS "")
    if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR "XRAYMCLIB_PGO is not supported with clang-cl")
    elseif(MSVC)
        # The .pgd file is written next to each executable, .pgc files from training runs are merged at link time.
        # /GENPROFILE:EXACT uses thread safe counters since transport is multithreaded.
        set(XRAYMCLIB_PGO_COMPILE_FLAGS /GL)
        if(XRAYMCLIB_PGO STREQUAL "GENERATE")
            set(XRAYMCLIB_PGO_LINK_FLAGS /LTCG /GENPROFILE:EXACT)
        else()
            set(XRAYMCLIB_PGO_LINK_FLAGS /LTCG /USEPROFILE)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(XRAYMCLIB_PGO STREQUAL "GENERATE")
            # atomic counter updates since transport is multithreaded
            set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-generate=${XRAYMCLIB_PGO_DIR} -fprofile-update=atomic)
            set(XRAYMCLIB_PGO_LINK_FLAGS ${XRAYMCLIB_PGO_COMPILE_FLAGS})
        else()
            # code not run in training is optimized as usual instead of for size,
            # and a profile from slightly different sources is a warning instead of an error
            set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-use=${XRAYMCLIB_PGO_DIR} -fprofile-partial-training
                -Wno-missing-profile -Wno-error=coverage-mismatch)
            set(XRAYMCLIB_PGO_LINK_FLAGS -fprofile-use=${XRAYMCLIB_PGO_DIR})
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # raw profiles must be merged with llvm-profdata before use, the pgo-train target does this
        get_filename_component(XRAYMCLIB_CXX_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
        string(REGEX MATCH "^[0-9]+" XRAYMCLIB_CLANG_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
        find_program(XRAYMCLIB_LLVM_PROFDATA NAMES llvm-profdata-${XRAYMCLIB_CLANG_MAJOR} llvm-profdata
            HINTS "${XRAYMCLIB_CXX_COMPILER_DIR}"
            DOC "llvm-profdata used to merge raw profiles")
        set(XRAYMCLIB_PGO_PROFDATA "${XRAYMCLIB_PGO_DIR}/default.profdata")
        if(XRAYMCLIB_PGO STREQUAL "GENERATE")
            set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-generate=${XRAYMCLIB_PGO_DIR} -fprofile-update=atomic)
            set(XRAYMCLIB_PGO_LINK_FLAGS ${XRAYMCLIB_PGO_COMPILE_FLAGS})
        else()
            if(NOT EXISTS "${XRAYMCLIB_PGO_PROFDATA}")
                message(WARNING "No merged profile at ${XRAYMCLIB_PGO_PROFDATA}, build and run pgo-train with XRAYMCLIB_PGO=GENERATE first")
            endif()
            set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-use=${XRAYMCLIB_PGO_PROFDATA}
                -Wno-profile-instr-unprofiled -Wno-profile-instr-out-of-date -Wno-backend-plugin)
            set(XRAYMCLIB_PGO_LINK_FLAGS -fprofile-use=${XRAYMCLIB_PGO_PROFDATA})
        endif()
    else()
        message(FATAL_ERROR "XRAYMCLIB_PGO is not supported for compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
    message(STATUS "Profile guided optimization: ${XRAYMCLIB_PGO} (profile data in ${XRAYMCLIB_PGO_DIR})")
    add_compile_options("$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_PGO_COMPILE_FLAGS}>")
    add_link_options("$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_PGO_LINK_FLAGS}>")
endif()

# Find packages go here.
find_package(Threads REQUIRED) # for threading support

## setting up material library since we use same headers for both xraymclib and EPICSparser
set(MATERIAL_HEADER_LIST 
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomicelement.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomicshell.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomhandler.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomserializer.hpp"
)

add_library(material INTERFACE ${MATERIAL_HEADER_LIST})
target_include_directories(material INTERFACE ${xraymclib_SOURCE_DIR}/include)

set(HEADER_LIST  
    "${xraymclib_SOURCE_DIR}/include/xraymc/xraymc.hpp"   
    "${xraymclib_SOURCE_DIR}/include/xraymc/xraymcrandom.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/interactions.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/transport.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/transportprogress.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/vectormath.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/utilities/spheresamplingrectangularfield.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/utilities/spheresamplingcircularfield.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/filters/bowtiefilter.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/filters/ctaecfilter.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/filters/ctorganaecfilter.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/tube/tube.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/tube/betheHeitlerCrossSection.hpp"  
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/beamtype.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/pencilbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropicmonoenergybeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropicmonoenergybeamcircle.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropicbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropicbeamcircle.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropiccircularmonoenergybeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/isotropiccircularbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/cbctbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/ctspiralbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/ctsequentialbeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/ctspiraldualenergybeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/ctdibeam.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/beams/flatmonoenergyfield.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/world.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/collision.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/dosescore.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/energyscore.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/statickdtree.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/kdtree.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/kdtreeflat.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/kdtreeintersectionresult.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/visualization/visualizeworld.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/visualization/vizualizationprops.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/collisions/collisiontest.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/collisions/meshcollision.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worldintersectionresult.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/personaldosimeter.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/ctdiphantom.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/flatdetector.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/depthdose.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worlditemtype.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/enclosedroom.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/fluencescore.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/aavoxelgrid.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/aafluencegrid.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worldbox.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worldboxgrid.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worldsphere.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worldcylinder.hpp"  
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedopensurface.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh/triangulatedmeshkdtree.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreeflat.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreetype.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/triangulatedmesh/triangle.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcollection.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"        
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcontourkdtree.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/basicshapes/aabb.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/basicshapes/cylinder.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/basicshapes/cylinderz.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/basicshapes/sphere.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/basicshapes/tetrahedron.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomicelement.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomicshell.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomhandler.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/atomserializer.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/nistmaterials.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/material/material.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/floating.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/constants.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/particle.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/particletracker.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/interpolation.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/serializer.hpp"
    )

add_library(libxraymc INTERFACE
    ${HEADER_LIST}
	)

target_include_directories(
    libxraymc
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)


# Linking to other libraries. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(libxraymc INTERFACE 
    Threads::Threads # Not needed on all platforms, but doesn't hurt
    )

# libstdc++ runs std::execution::par algorithms in parallel only with TBB
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    find_package(TBB QUIET)
    if(TBB_FOUND)
        target_link_libraries(libxraymc INTERFACE TBB::tbb)
    else()
        message(STATUS "TBB not found, std::execution parallel algorithms run serially")
    endif()
endif()

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

## Install
install(TARGETS libxraymc 
    EXPORT libxraymc    
    )

# EPICS data parser and physics data
add_subdirectory(data)

# Thirdparty libraries
# Including lodePNGwrapper for xraymc
add_subdirectory(thirdparty)

# Testing
option(BUILD_TESTING "Build tests" OFF)
if(BUILD_TESTING OR CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # enable_testing must be called from the top level directory for ctest to find the tests from the build root
    enable_testing()
    add_subdirectory(tests)

    if(XRAYMCLIB_PGO STREQUAL "GENERATE")
        # Runs the test suite with the instrumented build to collect profiles
        set(XRAYMCLIB_PGO_TRAIN_COMMANDS
            COMMAND "${CMAKE_CTEST_COMMAND}" -C $<CONFIG> --output-on-failure)
        if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT MSVC)
            if(NOT XRAYMCLIB_LLVM_PROFDATA)
                message(FATAL_ERROR "llvm-profdata is required for XRAYMCLIB_PGO with Clang, set XRAYMCLIB_LLVM_PROFDATA")
            endif()
            # merge every raw profile in the directory, a shell glob is not portable so let CMake expand it
            list(APPEND XRAYMCLIB_PGO_TRAIN_COMMANDS
                COMMAND "${CMAKE_COMMAND}" -DPROFDATA=${XRAYMCLIB_LLVM_PROFDATA} -DPGO_DIR=${XRAYMCLIB_PGO_DIR}
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/MergeClangProfiles.cmake")
        endif()
        add_custom_target(pgo-train
            ${XRAYMCLIB_PGO_TRAIN_COMMANDS}
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Collecting profiles for profile guided optimization"
            USES_TERMINAL
            VERBATIM)
    endif()
endif()

if(NOT XRAYMCLIB_PGO STREQUAL "OFF")
    # Removes collected profiles, counts from earlier training runs are otherwise accumulated
    add_custom_target(pgo-clean
        COMMAND "${CMAKE_COMMAND}" -E rm -rf "${XRAYMCLIB_PGO_DIR}"
        COMMENT "Removing profile data in ${XRAYMCLIB_PGO_DIR}"
        VERBATIM)
endif()

# Validation
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(validation)
endif()

# Docs by doxygen
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    find_package(Doxygen)
    find_package(Sphinx)
    if(Doxygen_FOUND AND Sphinx_FOUND)
	    add_subdirectory(docs)
    else()
	    message(STATUS "Doxygen or Sphinx not found, not building docs")
    endif()
endif()

# Examples
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(XRAYMCLIB_BUILD_EXAMPLES "Build examples" ON)
    if(XRAYMCLIB_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()
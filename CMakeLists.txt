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

option(XRAYMCLIB_USE_FAST_MATH "Build with fast floating point math" ON)
if(XRAYMCLIB_USE_FAST_MATH)
    if(MSVC)
        message(STATUS "Building with /fp:fast flag")
        add_compile_options("$<${XRAYMCLIB_OPT_CONFIG}:/fp:fast>")
    elseif(XRAYMCLIB_GNU_LIKE)
        # keep NaN/inf semantics so std::isfinite checks in mesh readers still work
        message(STATUS "Building with -ffast-math -fno-finite-math-only flags")
        add_compile_options("$<${XRAYMCLIB_OPT_CONFIG}:-ffast-math;-fno-finite-math-only>")
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

# Link time optimization: -flto on GCC/Clang, /GL + /LTCG on MSVC. Each test and example is a
# single translation unit of a header-only library, so the gain is usually small and link time grows a lot.
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
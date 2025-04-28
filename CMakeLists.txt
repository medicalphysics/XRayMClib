# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.20)

project(dxmclib VERSION 0.7.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

#Function for adding extra compiler flags if available
include(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()
# We want to add extra compiler warnings if supported
#enable_cxx_compiler_flag_if_supported("-Wall")
#enable_cxx_compiler_flag_if_supported("-Wextra")
#enable_cxx_compiler_flag_if_supported("-pedantic")

#MSVC fix for heavy templated code
if(MSVC)        
    MESSAGE("-- Building with /bigobj flag")
    enable_cxx_compiler_flag_if_supported("/bigobj")  
endif()

option(DXMCLIB_USE_FAST_MATH "Build with -fast-math or /fp:fast" ON)
if(DXMCLIB_USE_FAST_MATH)
    if(MSVC)        
        MESSAGE("-- Building with /fp:fast flag")
        enable_cxx_compiler_flag_if_supported("/fp:fast")  
    endif()
    if(UNIX)
        MESSAGE("-- Building with -ffast-math flag")
        enable_cxx_compiler_flag_if_supported("-ffast-math")  
    endif()
endif()

option(DXMCLIB_OPTIMIZE_FOR_NATIVE "Build with -march=native on gcc/clang or AVX_x on MSVC " ON)
if(DXMCLIB_OPTIMIZE_FOR_NATIVE)
    if(MSVC)
        set(DXMCLIB_MSVC_NATIVE_OPTION "AVX" CACHE STRING "Machine instructions support")
        set_property(CACHE DXMCLIB_MSVC_NATIVE_OPTION PROPERTY STRINGS SSE2 AVX AVX2 AVX512)
        MESSAGE("-- Trying to build for ${DXMCLIB_MSVC_NATIVE_OPTION} instruction set")
        enable_cxx_compiler_flag_if_supported("-arch:${DXMCLIB_MSVC_NATIVE_OPTION}")  
    endif()
    if(UNIX)
        MESSAGE("-- Trying to build for native (-march=native)")
        enable_cxx_compiler_flag_if_supported("-march=native")  
    endif()
endif()

# Find packages go here.
find_package(Threads REQUIRED) # for threading support

## setting up material library since we use same headers for both dxmclib and EPICSparser
set(MATERIAL_HEADER_LIST 
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomicelement.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomicshell.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomhandler.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomserializer.hpp"
)

add_library(material INTERFACE ${MATERIAL_HEADER_LIST})
target_include_directories(material INTERFACE ${dxmclib_SOURCE_DIR}/include)

set(HEADER_LIST         
    "${dxmclib_SOURCE_DIR}/include/dxmc/dxmcrandom.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/interactions.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/transport.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/transportprogress.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/vectormath.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/filters/bowtiefilter.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/filters/ctaecfilter.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/filters/ctorganaecfilter.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/tube/tube.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/tube/betheHeitlerCrossSection.hpp"  
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/beamtype.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/pencilbeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/isotropicmonoenergybeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/isotropicmonoenergybeamcircle.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/isotropicbeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/isotropicbeamcircle.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/cbctbeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/ctspiralbeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/ctsequentialbeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/ctspiraldualenergybeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/beams/ctdibeam.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/world.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/dosescore.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/energyscore.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/statickdtree.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/kdtree.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/kdtreeflat.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/kdtreeintersectionresult.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/visualization/visualizeworld.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/visualization/vizualizationprops.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/collisions/collisiontest.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/collisions/meshcollision.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worldintersectionresult.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/ctdiphantom.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/depthdose.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/worlditemtype.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/enclosedroom.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/fluencescore.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/aavoxelgrid.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/worldbox.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/worldboxgrid.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/worldsphere.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/worldcylinder.hpp"  
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedopensurface.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtree.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreeflat.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtreetype.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/triangulatedmesh/triangle.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh2.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshgrid.hpp"    
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshkdtree.hpp" 
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/basicshapes/aabb.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/basicshapes/cylinder.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/basicshapes/cylinderz.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/basicshapes/sphere.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/basicshapes/tetrahedron.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomicelement.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomicshell.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomhandler.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/atomserializer.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/nistmaterials.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material/material.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/floating.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/constants.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/particle.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/particletracker.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/interpolation.hpp"        
    )

add_library(libdxmc INTERFACE
    ${HEADER_LIST}
	)

target_include_directories(
    libdxmc
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# Linking to other libraries. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(libdxmc INTERFACE 
    Threads::Threads # Not needed on all platforms, but doesn't hurt
    )

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

## Install
install(TARGETS libdxmc 
    EXPORT libdxmc    
    )

# EPICS data parser and physics data
add_subdirectory(data)

# Thirdparty libraries
# Including lodePNGwrapper for dxmc
option(DXMCLIB_USE_LOADPNG "Use thirdparty loadpng library for writing png files" OFF)
add_subdirectory(thirdparty)

# Testing
option(BUILD_TESTING "Build tests" OFF)
if(BUILD_TESTING OR CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
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
    option(DXMCLIB_BUILD_EXAMPLES "Build examples" ON)
    if(DXMCLIB_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()
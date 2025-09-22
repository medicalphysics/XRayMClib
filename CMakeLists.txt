# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.20)

project(xraymclib VERSION 0.7.0 LANGUAGES CXX)

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

option(XRAYMCLIB_USE_FAST_MATH "Build with -fast-math or /fp:fast" ON)
if(XRAYMCLIB_USE_FAST_MATH)
    if(MSVC)        
        MESSAGE("-- Building with /fp:fast flag")
        enable_cxx_compiler_flag_if_supported("/fp:fast")  
    endif()
    if(UNIX)
        MESSAGE("-- Building with -ffast-math flag")
        enable_cxx_compiler_flag_if_supported("-ffast-math")  
    endif()
endif()

option(XRAYMCLIB_OPTIMIZE_FOR_NATIVE "Build with -march=native on gcc/clang or AVX_x on MSVC " ON)
if(XRAYMCLIB_OPTIMIZE_FOR_NATIVE)
    if(MSVC)
        set(XRAYMCLIB_MSVC_NATIVE_OPTION "AVX2" CACHE STRING "Machine instructions support")
        set_property(CACHE XRAYMCLIB_MSVC_NATIVE_OPTION PROPERTY STRINGS SSE2 AVX AVX2 AVX512)
        MESSAGE("-- Building with /arch:${XRAYMCLIB_MSVC_NATIVE_OPTION} flag")
        enable_cxx_compiler_flag_if_supported("/arch:${XRAYMCLIB_MSVC_NATIVE_OPTION}")          
    endif()
    if(UNIX)
        MESSAGE("-- Building with -march=native flag")
        enable_cxx_compiler_flag_if_supported("-march=native")  
    endif()
endif()

if(MSVC)
    option(XRAYMCLIB_WHOLE_PROGRAM_OPTIMALIZATION "Build with /GL flag" ON)
    if(XRAYMCLIB_WHOLE_PROGRAM_OPTIMALIZATION)
        MESSAGE("-- Building with /GL flag")
        enable_cxx_compiler_flag_if_supported("/GL")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LTCG")
        list(REMOVE_ITEM CMAKE_EXE_LINKER_FLAGS "/INCREMENTAL")
        enable_cxx_compiler_flag_if_supported("/Gw")  
        enable_cxx_compiler_flag_if_supported("/cgthreads8")  
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
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/world.hpp" 
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
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/ctdiphantom.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/depthdose.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/worlditemtype.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/enclosedroom.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/fluencescore.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/aavoxelgrid.hpp"
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
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh2.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh2/tetrahedalmeshgrid2.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshcollection.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshdata.hpp"
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshgrid.hpp"    
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedalmeshkdtree.hpp" 
    "${xraymclib_SOURCE_DIR}/include/xraymc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"
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

# Set type of material class
option(XRAYMCLIB_USE_LINEAR_INTERPOLATION "Use linear attenuation coeff interpolation" OFF)
if(XRAYMCLIB_USE_LINEAR_INTERPOLATION)
    target_compile_definitions(material INTERFACE -DMATERIAL_LINEAR_INTERPOLATION=1)
    target_compile_definitions(libxraymc INTERFACE -DMATERIAL_LINEAR_INTERPOLATION=1)
endif()

# Linking to other libraries. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(libxraymc INTERFACE 
    Threads::Threads # Not needed on all platforms, but doesn't hurt
    )

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
option(XRAYMCLIB_USE_LOADPNG "Use thirdparty loadpng library for writing png files" OFF)
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
    option(XRAYMCLIB_BUILD_EXAMPLES "Build examples" ON)
    if(XRAYMCLIB_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()
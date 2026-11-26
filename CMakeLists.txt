cmake_minimum_required(VERSION 3.22)

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

# Compiler families. MSVC is true for both cl and clang-cl.
set(XRAYMCLIB_CLANG_CL OFF)
set(XRAYMCLIB_GNU_LIKE OFF)
if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(XRAYMCLIB_CLANG_CL ON)
elseif(NOT MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(XRAYMCLIB_GNU_LIKE ON)
endif()

# Build and performance flags (fast math, native arch, LTO, PGO, ...). Every target in this project
# links xraymclib_flags, and consumers (FetchContent/add_subdirectory) get it through libxraymc.
# With PGO the consumer must also be instrumented and linked with the profiling runtime, otherwise
# linking against the instrumented helper libraries fails (undefined __gcov_* symbols).
add_library(xraymclib_flags INTERFACE)
# Optimization flags only apply to optimized configurations
set(XRAYMCLIB_OPT_CONFIG "$<CONFIG:Release,RelWithDebInfo>")

if(MSVC)
    # Heavy templated code
    target_compile_options(xraymclib_flags INTERFACE /bigobj)
    # /MP: compile source files in parallel within one cl invocation. Only has an effect with the
    # Visual Studio generators, Ninja and Makefiles already run one compiler process per file.
    if(NOT XRAYMCLIB_CLANG_CL)
        add_compile_options(/MP)
    endif()
    # For the static C runtime (/MT), configure with
    # -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>". Calls into the CRT then
    # become direct calls instead of going through DLL import stubs. Every library linked into the
    # same executable must use the same runtime.
endif()

option(XRAYMCLIB_USE_FAST_MATH "Build with fast floating point math" ON)
if(XRAYMCLIB_USE_FAST_MATH)
    # The individual parts of -ffast-math.
    #   -ffinite-math-only : assume no inf/NaN; code that may meet them must test explicitly
    #                        (xraymc::isFinite inspects the bits and is not optimized away)
    #   -fno-math-errno    : inline sqrt, allows vectorized log/exp/sin/cos
    #   -fno-trapping-math : FP operations are assumed not to trap
    #   -fno-signed-zeros  : -0.0 and +0.0 may be treated alike (required by -fassociative-math)
    #   -fassociative-math : reordering of sums, enables vectorized reductions
    #   -freciprocal-math  : x / y may become x * (1 / y)
    #   -ffp-contract=fast : fuse multiply-add into FMA instructions
    set(XRAYMCLIB_FAST_MATH_FLAGS
        -ffinite-math-only
        -fno-math-errno
        -fno-trapping-math
        -fno-signed-zeros
        -fassociative-math
        -freciprocal-math
        -ffp-contract=fast
    )
    if(XRAYMCLIB_CLANG_CL)
        # pass the same parts as GCC/Clang so all clang front ends build alike
        list(TRANSFORM XRAYMCLIB_FAST_MATH_FLAGS PREPEND "/clang:")
    elseif(MSVC)
        # MSVC has no separate finite-math switch, /fp:fast already lets the optimizer
        # assume that inf/NaN do not occur
        set(XRAYMCLIB_FAST_MATH_FLAGS /fp:fast)
    elseif(NOT XRAYMCLIB_GNU_LIKE)
        set(XRAYMCLIB_FAST_MATH_FLAGS "")
    endif()
    if(XRAYMCLIB_FAST_MATH_FLAGS)
        message(STATUS "Building with fast math flags: ${XRAYMCLIB_FAST_MATH_FLAGS}")
        target_compile_options(xraymclib_flags INTERFACE "$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_FAST_MATH_FLAGS}>")
    endif()
endif()

option(XRAYMCLIB_OPTIMIZE_FOR_NATIVE "Build for the host CPU (-march=native or /arch:X on MSVC), binaries may not run on other machines" ON)
if(XRAYMCLIB_OPTIMIZE_FOR_NATIVE)
    if(MSVC)
        set(XRAYMCLIB_MSVC_NATIVE_OPTION "AVX2" CACHE STRING "Machine instructions support")
        set_property(CACHE XRAYMCLIB_MSVC_NATIVE_OPTION PROPERTY STRINGS SSE2 AVX AVX2 AVX512)
        message(STATUS "Building with /arch:${XRAYMCLIB_MSVC_NATIVE_OPTION} flag")
        # /Ob3: more aggressive inlining than the Release default /Ob2
        # /Gw: put global data in separate sections so the linker can drop unused data
        target_compile_options(xraymclib_flags INTERFACE
            "/arch:${XRAYMCLIB_MSVC_NATIVE_OPTION}"
            "$<${XRAYMCLIB_OPT_CONFIG}:/Ob3;/Gw>")
    elseif(XRAYMCLIB_GNU_LIKE)
        check_cxx_compiler_flag("-march=native" XRAYMCLIB_HAS_MARCH_NATIVE)
        if(XRAYMCLIB_HAS_MARCH_NATIVE)
            set(XRAYMCLIB_NATIVE_FLAG -march=native)
        else()
            check_cxx_compiler_flag("-mcpu=native" XRAYMCLIB_HAS_MCPU_NATIVE) # older Apple clang on arm64
            if(XRAYMCLIB_HAS_MCPU_NATIVE)
                set(XRAYMCLIB_NATIVE_FLAG -mcpu=native)
            endif()
        endif()
        if(XRAYMCLIB_NATIVE_FLAG)
            message(STATUS "Building with ${XRAYMCLIB_NATIVE_FLAG} flag")
            target_compile_options(xraymclib_flags INTERFACE ${XRAYMCLIB_NATIVE_FLAG})
        endif()
    endif()
endif()

# Link time optimization: -flto on GCC/Clang, /GL + /LTCG on MSVC.
option(XRAYMCLIB_WHOLE_PROGRAM_OPTIMIZATION "Build with link time optimization" ON)
if(XRAYMCLIB_WHOLE_PROGRAM_OPTIMIZATION)
    check_ipo_supported(RESULT XRAYMCLIB_IPO_SUPPORTED OUTPUT XRAYMCLIB_IPO_ERROR)
    if(XRAYMCLIB_IPO_SUPPORTED)
        message(STATUS "Building with link time optimization")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
        # The IPO variables above do not reach consumers, pass the same flags CMake uses for IPO
        # (compile: e.g. -flto=auto -fno-fat-lto-objects or /GL, link: the same or /LTCG)
        set(XRAYMCLIB_IPO_LINK_FLAGS ${CMAKE_CXX_LINK_OPTIONS_IPO})
        if(XRAYMCLIB_CLANG_CL)
            # lld-link runs (Thin)LTO on bitcode objects by itself and does not understand
            # compiler flags such as -flto=thin or link.exe switches such as /CGTHREADS
            set(XRAYMCLIB_IPO_LINK_FLAGS "")
        elseif(MSVC)
            list(APPEND XRAYMCLIB_IPO_LINK_FLAGS /CGTHREADS:8)
        elseif(NOT XRAYMCLIB_IPO_LINK_FLAGS)
            set(XRAYMCLIB_IPO_LINK_FLAGS ${CMAKE_CXX_COMPILE_OPTIONS_IPO})
        endif()
        target_compile_options(xraymclib_flags INTERFACE "$<${XRAYMCLIB_OPT_CONFIG}:${CMAKE_CXX_COMPILE_OPTIONS_IPO}>")
        target_link_options(xraymclib_flags INTERFACE "$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_IPO_LINK_FLAGS}>")
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
    if(XRAYMCLIB_CLANG_CL)
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
    elseif(XRAYMCLIB_GNU_LIKE AND XRAYMCLIB_PGO STREQUAL "GENERATE")
        # atomic counter updates since transport is multithreaded
        set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-generate=${XRAYMCLIB_PGO_DIR} -fprofile-update=atomic)
        set(XRAYMCLIB_PGO_LINK_FLAGS ${XRAYMCLIB_PGO_COMPILE_FLAGS})
    elseif(XRAYMCLIB_GNU_LIKE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # code not run in training is optimized as usual instead of for size,
        # and a profile from slightly different sources is a warning instead of an error
        set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-use=${XRAYMCLIB_PGO_DIR} -fprofile-partial-training
            -Wno-missing-profile -Wno-error=coverage-mismatch)
        set(XRAYMCLIB_PGO_LINK_FLAGS -fprofile-use=${XRAYMCLIB_PGO_DIR})
    elseif(XRAYMCLIB_GNU_LIKE)
        # Clang: raw profiles are merged into one file by the pgo-train target
        set(XRAYMCLIB_PGO_PROFDATA "${XRAYMCLIB_PGO_DIR}/default.profdata")
        if(NOT EXISTS "${XRAYMCLIB_PGO_PROFDATA}")
            message(WARNING "No merged profile at ${XRAYMCLIB_PGO_PROFDATA}, build and run pgo-train with XRAYMCLIB_PGO=GENERATE first")
        endif()
        set(XRAYMCLIB_PGO_COMPILE_FLAGS -fprofile-use=${XRAYMCLIB_PGO_PROFDATA}
            -Wno-profile-instr-unprofiled -Wno-profile-instr-out-of-date -Wno-backend-plugin)
        set(XRAYMCLIB_PGO_LINK_FLAGS -fprofile-use=${XRAYMCLIB_PGO_PROFDATA})
    else()
        message(FATAL_ERROR "XRAYMCLIB_PGO is not supported for compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
    message(STATUS "Profile guided optimization: ${XRAYMCLIB_PGO} (profile data in ${XRAYMCLIB_PGO_DIR})")
    target_compile_options(xraymclib_flags INTERFACE "$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_PGO_COMPILE_FLAGS}>")
    target_link_options(xraymclib_flags INTERFACE "$<${XRAYMCLIB_OPT_CONFIG}:${XRAYMCLIB_PGO_LINK_FLAGS}>")
endif()

find_package(Threads REQUIRED)

# Material headers are shared by libxraymc and the EPICS parser
add_library(material INTERFACE
    include/xraymc/material/atomicelement.hpp
    include/xraymc/material/atomicshell.hpp
    include/xraymc/material/atomhandler.hpp
    include/xraymc/material/atomserializer.hpp
)
target_include_directories(material INTERFACE ${PROJECT_SOURCE_DIR}/include)

# Header only library, headers are listed for IDEs
file(GLOB_RECURSE XRAYMCLIB_HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/xraymc/*.hpp")
add_library(libxraymc INTERFACE ${XRAYMCLIB_HEADERS})
target_include_directories(libxraymc INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
# Code including the xraymc headers must be compiled as C++23
target_compile_features(libxraymc INTERFACE cxx_std_23)
# BUILD_INTERFACE: host specific flags and absolute profile paths are not exported on install
target_link_libraries(libxraymc INTERFACE Threads::Threads $<BUILD_INTERFACE:xraymclib_flags>)

# libstdc++ runs std::execution::par algorithms in parallel only with TBB
if(XRAYMCLIB_GNU_LIKE)
    find_package(TBB QUIET)
    if(TBB_FOUND)
        target_link_libraries(libxraymc INTERFACE TBB::tbb)
    else()
        message(STATUS "TBB not found, std::execution parallel algorithms run serially")
    endif()
endif()

include(GNUInstallDirs)
install(TARGETS libxraymc EXPORT libxraymc)
install(DIRECTORY include/xraymc TYPE INCLUDE FILES_MATCHING PATTERN "*.hpp")

# EPICS data parser and physics data
add_subdirectory(data)

# Thirdparty libraries
# Including lodePNGwrapper for xraymc
add_subdirectory(thirdparty)

option(XRAYMCLIB_BUILD_TESTS "Build tests" ${PROJECT_IS_TOP_LEVEL})
if(XRAYMCLIB_BUILD_TESTS)
    # enable_testing must be called from the top level directory for ctest to find the tests from the build root
    enable_testing()
    add_subdirectory(tests)

    if(XRAYMCLIB_PGO STREQUAL "GENERATE")
        # Runs the test suite with the instrumented build to collect profiles
        set(XRAYMCLIB_PGO_TRAIN_COMMANDS
            COMMAND "${CMAKE_CTEST_COMMAND}" -C $<CONFIG> --output-on-failure)
        if(XRAYMCLIB_GNU_LIKE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            get_filename_component(XRAYMCLIB_CXX_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
            string(REGEX MATCH "^[0-9]+" XRAYMCLIB_CLANG_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
            find_program(XRAYMCLIB_LLVM_PROFDATA NAMES llvm-profdata-${XRAYMCLIB_CLANG_MAJOR} llvm-profdata
                HINTS "${XRAYMCLIB_CXX_COMPILER_DIR}"
                DOC "llvm-profdata used to merge raw profiles")
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
        # build the test executables before running them
        add_dependencies(pgo-train xraymclib_tests)
    endif()
endif()

if(NOT XRAYMCLIB_PGO STREQUAL "OFF")
    # Removes collected profiles, counts from earlier training runs are otherwise accumulated
    add_custom_target(pgo-clean
        COMMAND "${CMAKE_COMMAND}" -E rm -rf "${XRAYMCLIB_PGO_DIR}"
        COMMENT "Removing profile data in ${XRAYMCLIB_PGO_DIR}"
        VERBATIM)
endif()

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(validation)

    # Docs by doxygen and sphinx
    find_package(Doxygen)
    find_package(Sphinx)
    if(Doxygen_FOUND AND Sphinx_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen or Sphinx not found, not building docs")
    endif()

    option(XRAYMCLIB_BUILD_EXAMPLES "Build examples" ON)
    if(XRAYMCLIB_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()

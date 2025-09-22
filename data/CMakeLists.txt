

# EPICS data parser
SET(XRAYMCLIB_EPICS_PHYSICSLISTS "${CMAKE_CURRENT_BINARY_DIR}/physicslists.bin" CACHE INTERNAL "Path to EPICS physicslists")

SET(XRAYMCLIB_EPICS_DATA_DIRPATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Path to directory holding EPDL2025.ALL and EADL2025.ALL in ENDL format")
if(XRAYMCLIB_EPICS_DATA_DIRPATH)
    SET(XRAYMCLIB_EPICS_DATA_EPDLPATH "${XRAYMCLIB_EPICS_DATA_DIRPATH}/EPDL2025.ALL" CACHE FILEPATH "Path to EPDL2025.ALL" FORCE)
    SET(XRAYMCLIB_EPICS_DATA_EADLPATH "${XRAYMCLIB_EPICS_DATA_DIRPATH}/EADL2025.ALL" CACHE FILEPATH "Path to EADL2025.ALL" FORCE)

    SET(XRAYMCLIB_EPICS_DOWNLOAD OFF CACHE BOOL "Download EPICS data from https://nuclear.llnl.gov/EPICS/ENDL2025 if necessary.")

    #CHECK IF EPICS DATAFILES EXISTS
    if(NOT EXISTS "${XRAYMCLIB_EPICS_DATA_EPDLPATH}" AND "${XRAYMCLIB_EPICS_DOWNLOAD}" )
        MESSAGE("Downloading EPDL")
        FILE(DOWNLOAD "https://nuclear.llnl.gov/EPICS/ENDL2025/EPDL2025.ALL" "${XRAYMCLIB_EPICS_DATA_EPDLPATH}" SHOW_PROGRESS)
    endif()
    if(NOT EXISTS "${XRAYMCLIB_EPICS_DATA_EADLPATH}" AND "${XRAYMCLIB_EPICS_DOWNLOAD}" )
        MESSAGE("Downloading EADL")
        FILE(DOWNLOAD "https://nuclear.llnl.gov/EPICS/ENDL2025/EADL2025.ALL" "${XRAYMCLIB_EPICS_DATA_EADLPATH}" SHOW_PROGRESS)
    endif()

    if(NOT EXISTS "${XRAYMCLIB_EPICS_DATA_EADLPATH}" OR NOT EXISTS "${XRAYMCLIB_EPICS_DATA_EPDLPATH}")
        message(FATAL_ERROR "Could not find data files EPDL2025.ALL or EPDL2025.ALL in ${XRAYMCLIB_EPICS_DATA_DIRPATH}")
    endif()

    message(STATUS "Using EPDL file: ${XRAYMCLIB_EPICS_DATA_EPDLPATH}")
    message(STATUS "Using EADL file: ${XRAYMCLIB_EPICS_DATA_EADLPATH}")
endif()

# Option to specify max energy range
set(XRAYMCLIB_MAX_ENERGY_KEV "500.0" CACHE STRING "Max energy in keV for photon energies.")
set(XRAYMCLIB_MIN_ENERGY_KEV "1.0")
set(XRAYMCLIB_MIN_TABLE_ENERGY_KEV "1.0")

add_executable(EPICSparser         
    EPICSParser/atomicelementhandler.cpp 
    EPICSParser/epicsparser.cpp 
    EPICSParser/epicsparsermain.cpp        
)

OPTION(XRAYMCLIB_USE_NIST_DATA "Use NIST attenuation data" OFF)
if(${XRAYMCLIB_USE_NIST_DATA})
IF(WIN32)
    set(XRAYLIB_INSTALL_PATH "C:/Program Files/xraylib 64-bit" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "${XRAYLIB_INSTALL_PATH}/Include")
	set(XRAYLIB_LIBRARY_SHARED_NAME "libxrl-11.dll")
	set(XRAYLIB_LIBRARY_SHARED_PATH "C:/Windows/sysnative/${XRAYLIB_LIBRARY_SHARED_NAME}")
    set(XRAYLIB_LIBRARY_SHARED_PATH_64 "C:/Windows/System32/${XRAYLIB_LIBRARY_SHARED_NAME}")
    find_library(XRAYLIB NAMES libxrl.lib libxrl libxrl-11 xrl HINTS "${XRAYLIB_INSTALL_PATH}/LIB")
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
    set(XRAYLIB_INSTALL_PATH "/usr/local/lib" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "/usr/local/include/xraylib")
    find_library(XRAYLIB NAMES libxrl libxrl-11 xrl HINTS "${XRAYLIB_INSTALL_PATH}")
ENDIF()

target_include_directories(EPICSparser PRIVATE ${XRAYLIB_INCLUDES})
# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(EPICSparser PRIVATE ${XRAYLIB})
target_compile_definitions(EPICSparser PRIVATE -DUSE_NIST)
endif()


OPTION(XRAYMCLIB_ALWAYS_REBUILD_PHYSICSLISTS "Rebuild physicslists on compile" OFF)
if(${XRAYMCLIB_ALWAYS_REBUILD_PHYSICSLISTS})
    target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_REBUILD_PHYSICSLISTS=1)
else()
    target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_REBUILD_PHYSICSLISTS=0)
endif()

target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_MAXENERGY="${XRAYMCLIB_MAX_ENERGY_KEV}")
target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_MINTABLEENERGY="${XRAYMCLIB_MIN_TABLE_ENERGY_KEV}")
target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_EADLPATH="${XRAYMCLIB_EPICS_DATA_EADLPATH}")
target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_EPDLPATH="${XRAYMCLIB_EPICS_DATA_EPDLPATH}")
target_compile_definitions(EPICSparser PUBLIC XRAYMCLIB_PHYSICSLISTSPATH="${XRAYMCLIB_EPICS_PHYSICSLISTS}")
target_include_directories(EPICSparser PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/include)     
target_link_libraries(EPICSparser PRIVATE material )

## Copying standard densities for all elements up to 92
add_custom_command(TARGET EPICSparser PRE_BUILD
    COMMENT "Copying standarddensities.csv to ${CMAKE_CURRENT_BINARY_DIR}/standarddensities.csv"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/standarddensities.csv
            "$<TARGET_FILE_DIR:EPICSparser>/standarddensities.csv"
)

## Copying Compton profiles since this data is lacking from EPICS data
add_custom_command(TARGET EPICSparser PRE_BUILD
    COMMENT "Copying ComptonProfiles.dat to ${CMAKE_CURRENT_BINARY_DIR}/ComptonProfiles.dat"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/ComptonProfiles.dat
            "$<TARGET_FILE_DIR:EPICSparser>/ComptonProfiles.dat"
)
# Copying standard densities for all elements up to 92
add_custom_command(TARGET EPICSparser PRE_BUILD
COMMENT "Copying standarddensities.csv to ${CMAKE_CURRENT_BINARY_DIR}/standarddensities.csv"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/standarddensities.csv
            ${CMAKE_CURRENT_BINARY_DIR}/standarddensities.csv
)

## Copying Compton profiles since this data is lacking from EPICS data
add_custom_command(TARGET EPICSparser PRE_BUILD
    COMMENT "Copying ComptonProfiles.dat to ${CMAKE_CURRENT_BINARY_DIR}/ComptonProfiles.dat"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/ComptonProfiles.dat
            ${CMAKE_CURRENT_BINARY_DIR}/ComptonProfiles.dat
)


## make parsed EPICS data available to all binaries consuming xraymc on build machine
target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_PHYSICSLISTSPATH="${XRAYMCLIB_EPICS_PHYSICSLISTS}")

## make max and min energy exposed to libxraymc
target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_MAXENERGY=${XRAYMCLIB_MAX_ENERGY_KEV})
target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_MINENERGY=${XRAYMCLIB_MIN_ENERGY_KEV})

## Generating physicslist
add_custom_target(libxraymc_physics BYPRODUCTS "${XRAYMCLIB_EPICS_PHYSICSLISTS}" COMMAND EPICSparser COMMENT "Generating physicslist.bin")

add_dependencies(libxraymc libxraymc_physics)
# Adding physics data to install folder
install(FILES ${XRAYMCLIB_EPICS_PHYSICSLISTS} DESTINATION ${CMAKE_INSTALL_BINDIR})

# Function to copy generated physicslist binary file to target binary directory
# Dependency on physicslist.bin file is handles by linking target to libxraymc
function(xraymclib_add_physics_list func_target)
        add_custom_command(
            TARGET ${func_target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${XRAYMCLIB_EPICS_PHYSICSLISTS}"
                "$<TARGET_FILE_DIR:${func_target}>/physicslist.bin"                         
        )        
endfunction()

# Function to copy icrp145 organ data to target binary directory
SET(XRAYMCLIB_ICRP145ORGANDATA "${CMAKE_CURRENT_SOURCE_DIR}/phantomdata/icrp145/icrp145organs.csv" CACHE INTERNAL "Path to ICRP145 organ data")
function(xraymclib_add_icrp145_organ_data func_target)
        add_custom_command(
            TARGET ${func_target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${XRAYMCLIB_ICRP145ORGANDATA}"
                "$<TARGET_FILE_DIR:${func_target}>/icrp145organs.csv"            
        )        
endfunction()

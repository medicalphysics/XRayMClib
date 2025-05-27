

# EPICS data parser
SET(DXMCLIB_EPICS_PHYSICSLISTS "${CMAKE_CURRENT_BINARY_DIR}/physicslists.bin" CACHE INTERNAL "Path to EPICS physicslists")

SET(DXMCLIB_EPICS_DATA_DIRPATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Path to directory holding EPDL2025.ALL and EADL2025.ALL in ENDL format")
if(DXMCLIB_EPICS_DATA_DIRPATH)
    SET(DXMCLIB_EPICS_DATA_EPDLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EPDL2025.ALL" CACHE FILEPATH "Path to EPDL2025.ALL" FORCE)
    SET(DXMCLIB_EPICS_DATA_EADLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EADL2025.ALL" CACHE FILEPATH "Path to EADL2025.ALL" FORCE)

    SET(DXMCLIB_EPICS_DOWNLOAD OFF CACHE BOOL "Download EPICS data from https://www-nds.iaea.org/epics/ if necessary.")

    #CHECK IF EPICS DATAFILES EXISTS
    if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EPDLPATH}" AND "${DXMCLIB_EPICS_DOWNLOAD}" )
        MESSAGE("Downloading EPDL")
        FILE(DOWNLOAD "https://www-nds.iaea.org/epics/ENDL2023/EPDL2025.ALL" "${DXMCLIB_EPICS_DATA_EPDLPATH}" SHOW_PROGRESS)
    endif()
    if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EADLPATH}" AND "${DXMCLIB_EPICS_DOWNLOAD}" )
        MESSAGE("Downloading EADL")
        FILE(DOWNLOAD "https://www-nds.iaea.org/epics/ENDL2023/EADL2025.ALL" "${DXMCLIB_EPICS_DATA_EADLPATH}" SHOW_PROGRESS)
    endif()

    if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EADLPATH}" OR NOT EXISTS "${DXMCLIB_EPICS_DATA_EPDLPATH}")
        message(FATAL_ERROR "Could not find data files EPDL2025.ALL or EPDL2025.ALL in ${DXMCLIB_EPICS_DATA_DIRPATH}")
    endif()

    message(STATUS "Using EPDL file: ${DXMCLIB_EPICS_DATA_EPDLPATH}")
    message(STATUS "Using EADL file: ${DXMCLIB_EPICS_DATA_EADLPATH}")
endif()

# Option to specify max energy range
set(DXMCLIB_MAX_ENERGY_KEV "500.0" CACHE STRING "Max energy in keV for photon energies.")
set(DXMCLIB_MIN_ENERGY_KEV "1.0")
set(DXMCLIB_MIN_TABLE_ENERGY_KEV "1.0")

add_executable(EPICSparser         
    EPICSParser/atomicelementhandler.cpp 
    EPICSParser/epicsparser.cpp 
    EPICSParser/epicsparsermain.cpp        
)

OPTION(DXMCLIB_ALWAYS_REBUILD_PHYSICSLISTS "Rebuild physicslists on compile" OFF)
if(${DXMCLIB_ALWAYS_REBUILD_PHYSICSLISTS})
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_REBUILD_PHYSICSLISTS=1)
else()
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_REBUILD_PHYSICSLISTS=0)
endif()

target_compile_definitions(EPICSparser PUBLIC DXMCLIB_MAXENERGY="${DXMCLIB_MAX_ENERGY_KEV}")
target_compile_definitions(EPICSparser PUBLIC DXMCLIB_MINTABLEENERGY="${DXMCLIB_MIN_TABLE_ENERGY_KEV}")
target_compile_definitions(EPICSparser PUBLIC DXMCLIB_EADLPATH="${DXMCLIB_EPICS_DATA_EADLPATH}")
target_compile_definitions(EPICSparser PUBLIC DXMCLIB_EPDLPATH="${DXMCLIB_EPICS_DATA_EPDLPATH}")
target_compile_definitions(EPICSparser PUBLIC DXMCLIB_PHYSICSLISTSPATH="${DXMCLIB_EPICS_PHYSICSLISTS}")
target_include_directories(EPICSparser PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/include)     
target_link_libraries(EPICSparser PRIVATE material )

## Copying standard densities for all elements up to 92
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

## make parsed EPICS data available to all binaries consuming dxmc on build machine
target_compile_definitions(libdxmc INTERFACE DXMCLIB_PHYSICSLISTSPATH="${DXMCLIB_EPICS_PHYSICSLISTS}")

## make max and min energy exposed to libdxmc
target_compile_definitions(libdxmc INTERFACE DXMCLIB_MAXENERGY=${DXMCLIB_MAX_ENERGY_KEV})
target_compile_definitions(libdxmc INTERFACE DXMCLIB_MINENERGY=${DXMCLIB_MIN_ENERGY_KEV})

## Generating physicslist
add_custom_target(libdxmc_physics BYPRODUCTS "${DXMCLIB_EPICS_PHYSICSLISTS}" COMMAND EPICSparser COMMENT "Generating physicslist.bin")

add_dependencies(libdxmc libdxmc_physics)
# Adding physics data to install folder
install(FILES ${DXMCLIB_EPICS_PHYSICSLISTS} DESTINATION ${CMAKE_INSTALL_BINDIR})

# Function to copy generated physicslist binary file to target binary directory
# Dependency on physicslist.bin file is handles by linking target to libdxmc
function(dxmclib_add_physics_list func_target)
        add_custom_command(
            TARGET ${func_target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${DXMCLIB_EPICS_PHYSICSLISTS}"
                "$<TARGET_FILE_DIR:${func_target}>/physicslist.bin"                         
        )        
endfunction()

# Function to copy icrp145 organ data to target binary directory
SET(DXMCLIB_ICRP145ORGANDATA "${CMAKE_CURRENT_SOURCE_DIR}/phantomdata/icrp145/icrp145organs.csv" CACHE INTERNAL "Path to ICRP145 organ data")
function(dxmclib_add_icrp145_organ_data func_target)
        add_custom_command(
            TARGET ${func_target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${DXMCLIB_ICRP145ORGANDATA}"
                "$<TARGET_FILE_DIR:${func_target}>/icrp145organs.csv"            
        )        
endfunction()

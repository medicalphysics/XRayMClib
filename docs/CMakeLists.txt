
option(XRAYMCLIB_BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen and Sphinx)" ${DOXYGEN_FOUND})

if(XRAYMCLIB_BUILD_DOCUMENTATION)
    
    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doxygen ALL
        DEPENDS libxraymc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION documentation)

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
    add_custom_target(sphinx ALL
        DEPENDS doxygen
        COMMAND ${SPHINX_EXECUTABLE} -b html
        # Tell Breathe where to find the Doxygen output
        -Dbreathe_projects.xraymclib=${CMAKE_CURRENT_BINARY_DIR}/xml
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx DESTINATION help)
endif()
# Tests are split in two groups:
#  - unit: fast, mostly deterministic tests of geometry, math, materials, beams and items
#  - transport: Monte Carlo transport tests that compare results statistically
#
# Run all tests from the build directory with
#   ctest -C Release
# or a single group with
#   ctest -C Release -L unit
#
# Build every test executable at once with the target xraymclib_tests.

set(XRAYMCLIB_UNIT_TESTS
    "testvectormath.cpp"
    "testrandom.cpp"
    "testinterpolation.cpp"
    "testbasicshapes.cpp"
    "testmaterials.cpp"
    "testinteractions.cpp"
    "testdirectionsampler.cpp"
    "testtube.cpp"
    "testbeamfilters.cpp"
    "testbeams.cpp"
    "testworlditems.cpp"
    "testserializer.cpp"
    "testscoring.cpp"
    "testflatdetector.cpp"
    "testfluencescore.cpp"
    "testtriangulatedmesh.cpp"
    "testtetrahedalmesh.cpp"
    "testworld.cpp"
    "testvisualization.cpp"
    "testcollision.cpp"
)

set(XRAYMCLIB_TRANSPORT_TESTS
    "testtransport.cpp"
    "testpersonaldosimeter.cpp"
    "testvoxelgrid.cpp"
    "testctdiphantom.cpp"
    "testworldsphere.cpp"
)

add_custom_target(xraymclib_tests)

function(xraymclib_add_test sourcefile label timeout)
    get_filename_component(testname ${sourcefile} NAME_WLE)
    add_executable(${testname} ${sourcefile} testutils.hpp testgeometry.hpp)
    target_link_libraries(${testname} libxraymc)
    add_test(NAME ${testname} COMMAND ${testname})
    set_tests_properties(${testname} PROPERTIES LABELS ${label} TIMEOUT ${timeout})
    add_dependencies(xraymclib_tests ${testname})
    # Adding interaction binary data to test
    xraymclib_add_physics_list(${testname})
endfunction()

foreach(sourcefile ${XRAYMCLIB_UNIT_TESTS})
    xraymclib_add_test(${sourcefile} unit 300)
endforeach()

foreach(sourcefile ${XRAYMCLIB_TRANSPORT_TESTS})
    xraymclib_add_test(${sourcefile} transport 900)
endforeach()

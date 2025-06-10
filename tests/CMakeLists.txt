include(CTest)
ENABLE_TESTING()

set(TEST_SOURCES 
    "testmaterials.cpp"
    "testbeams.cpp"
    "testctaecfilter.cpp"
    "testworlditems.cpp" 
    "testbasicshapes.cpp"
    "testtube.cpp"    
    "testbeamfilters.cpp"   
    "testrandom.cpp"
    "testvectormath.cpp"
    "testworldsphere.cpp"
    "testmassenergytransfer.cpp"
    "testtriangulatedmesh.cpp"
    "testcollision.cpp"
    "testtetrahedalmeshreader.cpp"
    "testtetrahedalmesh2.cpp"
)


foreach(sourcefile ${TEST_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WLE)
    add_executable(${filename} ${sourcefile})    
    target_link_libraries(${filename} libdxmc )
    add_test(NAME ${filename} COMMAND ${filename})    
    
    #Adding interaction binary data to example 
    dxmclib_add_physics_list(${filename})
endforeach()

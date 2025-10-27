

## Use loadpng thirdparty library
option(XRAYMCLIB_USE_LOADPNG "Use thirdparty loadpng library for writing png files" ON)
if(XRAYMCLIB_USE_LOADPNG)   
    ## setting up lodepng library files
    set(LODEPNG_FILES 
        "${xraymclib_SOURCE_DIR}/thirdparty/include/lodepng/xraymclodepngwrapper.hpp"
        "${xraymclib_SOURCE_DIR}/thirdparty/src/lodepng/xraymclodepngwrapper.cpp"
        "${xraymclib_SOURCE_DIR}/thirdparty/src/lodepng/lodepng.cpp"
        "${xraymclib_SOURCE_DIR}/thirdparty/src/lodepng/lodepng.h"
    )

    add_library(lodepngwrapper STATIC ${LODEPNG_FILES})
    target_include_directories(lodepngwrapper PUBLIC ${xraymclib_SOURCE_DIR}/thirdparty/include)
    target_link_libraries(libxraymc INTERFACE lodepngwrapper)
    target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_USE_LOADPNG=1)
    message(STATUS "Using loadpng thirdparty library for png files")
endif(XRAYMCLIB_USE_LOADPNG)

## Add other thirdparty libraries here...
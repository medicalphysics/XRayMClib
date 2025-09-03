

## Use loadpng thirdparty library
if(DXMCLIB_USE_LOADPNG)   
    ## setting up lodepng library files
    set(LODEPNG_FILES 
        "${dxmclib_SOURCE_DIR}/thirdparty/include/lodepng/dxmclodepngwrapper.hpp"
        "${dxmclib_SOURCE_DIR}/thirdparty/src/lodepng/dxmclodepngwrapper.cpp"
        "${dxmclib_SOURCE_DIR}/thirdparty/src/lodepng/lodepng.cpp"
        "${dxmclib_SOURCE_DIR}/thirdparty/src/lodepng/lodepng.h"
    )

    add_library(lodepngwrapper STATIC ${LODEPNG_FILES})
    target_include_directories(lodepngwrapper PUBLIC ${dxmclib_SOURCE_DIR}/thirdparty/include)
    target_link_libraries(libdxmc INTERFACE lodepngwrapper)
    target_compile_definitions(libdxmc INTERFACE DXMCLIB_USE_LOADPNG=1)
    message(STATUS "Using loadpng thirdparty library for png files")
endif(DXMCLIB_USE_LOADPNG)

## Add other thirdparty libraries here...
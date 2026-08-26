

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

## Use smallz4 thirdparty library for LZ4 buffer compression (used by Serializer)
set(SMALLZ4_FILES
    "${xraymclib_SOURCE_DIR}/thirdparty/include/smallz4/smallz4.h"
    "${xraymclib_SOURCE_DIR}/thirdparty/include/smallz4/xraymclz4wrapper.hpp"
)

add_library(smallz4wrapper INTERFACE ${SMALLZ4_FILES})
target_include_directories(smallz4wrapper INTERFACE ${xraymclib_SOURCE_DIR}/thirdparty/include)
target_link_libraries(libxraymc INTERFACE smallz4wrapper)
message(STATUS "Using smallz4 thirdparty library for buffer compression")

## Add other thirdparty libraries here...
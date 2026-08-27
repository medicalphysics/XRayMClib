

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


# Use zlib thirdparty library for buffer compression
option(XRAYMCLIB_USE_ZLIB "Use zlib thirdparty library, if found, for buffer compression" ON)
if(XRAYMCLIB_USE_ZLIB)
    ## Use zlib for buffer compression (used by Serializer)
    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
        message(STATUS "Found zlib library, can be used for buffer compression")            
        set(ZLIBWRAPPER_FILES
            "${xraymclib_SOURCE_DIR}/thirdparty/include/zlib/xraymczlibwrapper.hpp"
            "${xraymclib_SOURCE_DIR}/thirdparty/src/zlib/xraymczlibwrapper.cpp"
        )

        add_library(zlibwrapper STATIC ${ZLIBWRAPPER_FILES})
        target_include_directories(zlibwrapper PUBLIC ${xraymclib_SOURCE_DIR}/thirdparty/include)
        target_link_libraries(zlibwrapper PUBLIC ZLIB::ZLIB)
        target_link_libraries(libxraymc INTERFACE zlibwrapper)
        target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_USE_ZLIB=1)    

    else()
        message(STATUS "zlib NOT FOUND. Not using zlib library for buffer compression")
    endif(ZLIB_FOUND)
endif(XRAYMCLIB_USE_ZLIB)


 #   target_compile_definitions(libxraymc INTERFACE XRAYMCLIB_USE_LOADPNG=1)
 #   add_compile_definitions(XRAYMCLIB_USE_ZLIB=1)
 #   message(STATUS "Using zlib thirdparty library for buffer compression")



## Add other thirdparty libraries here...
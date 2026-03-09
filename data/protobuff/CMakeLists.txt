

SET(XRAYMC_PROTOC_PATH_HINT "C:/Users/ander/Downloads/protoc-34.0-win64/bin")
find_program(XRAYMC_PROTOBUF_PROTOC protoc HINTS ${XRAYMC_PROTOC_PATH_HINT} REQUIRED)

SET(XRAYMC_PROTOBUFF_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${XRAYMC_PROTOBUFF_DIR})

# set the arguments for the protobuf compiler
set(XRAYMC_PROTOBUFF_ARGS 
  "--cpp_out=${XRAYMC_PROTOBUFF_DIR}"
  "-I"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/dosecollection.proto"
)

set(XRAYMC_PROTOBUFF_OUTPUT
  "${XRAYMC_PROTOBUFF_DIR}/dosecollection.pb.h" 
  "${XRAYMC_PROTOBUFF_DIR}/dosecollection.pb.cc"
)

# add custom command to generate the protobuf sources
add_custom_command(
  OUTPUT ${XRAYMC_PROTOBUFF_OUTPUT}
  COMMAND ${XRAYMC_PROTOBUF_PROTOC}
  ARGS ${XRAYMC_PROTOBUFF_ARGS}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/dosecollection.proto
  COMMENT "Running C++ protocol buffer compiler on dosecollection.proto"
  VERBATIM
)


#set_source_files_properties(${XRAYMC_PROTOBUFF_OUTPUT} PROPERTIES GENERATED TRUE)
add_custom_target(libxraymc_protobuff DEPENDS ${XRAYMC_PROTOBUFF_OUTPUT})
#add_dependencies(libxraymc libxraymc_protobuff)

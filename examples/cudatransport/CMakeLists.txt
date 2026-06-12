
# CudaTransport example — only built when libxraymc_cuda is available.
if(NOT TARGET libxraymc_cuda)
    return()
endif()

add_executable(cudatransport cudatransport.cpp)
target_link_libraries(cudatransport PRIVATE libxraymc_cuda)
xraymclib_add_physics_list(cudatransport)

install(TARGETS cudatransport
    RUNTIME DESTINATION examples)

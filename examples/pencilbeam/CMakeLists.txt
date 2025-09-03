
# testing to see if dxmclib target is defined. If it's not we will download 
# dxmclib at configure time.
IF( NOT TARGET libdxmc)
	include(FetchContent)
	## Adding DXMClib package
	FetchContent_Declare(
		libdxmc
		GIT_REPOSITORY https://github.com/medicalphysics/DXMClib.git
		GIT_TAG master
		)
	FetchContent_MakeAvailable(libdxmc)
ENDIF()

#Define an executable
add_executable(pencilbeam pencilbeam.cpp)

#Linking to dxmclib
target_link_libraries(pencilbeam PRIVATE libdxmc)

#Adding interaction binary data to example 
#This cmake function vil copy physicslist.bin to executable folder
dxmclib_add_physics_list(pencilbeam)

install(TARGETS pencilbeam
	RUNTIME DESTINATION examples)

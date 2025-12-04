
# testing to see if xraymclib target is defined. If it's not we will download 
# xraymclib at configure time.
IF( NOT TARGET libxraymc)
	include(FetchContent)
	## Adding xraymclib package
	FetchContent_Declare(
		libxraymc
		GIT_REPOSITORY https://github.com/medicalphysics/xraymclib.git
		GIT_TAG master
		)
	FetchContent_MakeAvailable(libxraymc)
ENDIF()

#Define an executable
add_executable(donut donut.cpp)

#Linking to xraymclib
target_link_libraries(donut PRIVATE libxraymc)

#Adding interaction binary data to example 
#This cmake function vil copy physicslist.bin to executable folder
xraymclib_add_physics_list(donut)

install(TARGETS donut
	RUNTIME DESTINATION examples)

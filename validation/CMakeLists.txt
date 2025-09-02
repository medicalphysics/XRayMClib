

add_executable(validation validation.cpp tg195world42.hpp tg195world3breast.hpp)
target_link_libraries(validation libxraymc)

add_executable(scattervalidation scattervalidation.cpp)
target_link_libraries(scattervalidation libxraymc)

## Adding post build step to copy physicslist.bin to validation folder
xraymclib_add_physics_list(validation)
xraymclib_add_physics_list(scattervalidation)

set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(CASE5WORLD_BINARY "case5world.bin")
set(CASE5WORLD_COMPRESSED "case5world.tar.gz")

if(EXISTS "${DATA_PATH}/${CASE5WORLD_BINARY}")
	# If this binary file is in sources we compress it
	execute_process(COMMAND
	  ${CMAKE_COMMAND} -E tar -czvf ${CASE5WORLD_COMPRESSED} ${CASE5WORLD_BINARY}
	  WORKING_DIRECTORY ${DATA_PATH}
	  )
	# Then we delete it
	file(REMOVE "${DATA_PATH}/${CASE5WORLD_BINARY}")
endif()

add_custom_command(TARGET validation POST_BUILD
	COMMENT "Extracting case 5 world ${DATA_PATH}/${CASE5WORLD_COMPRESSED}."	
	COMMAND ${CMAKE_COMMAND} -E tar -xvzf "${DATA_PATH}/${CASE5WORLD_COMPRESSED}" ${CASE5WORLD_BINARY}
	WORKING_DIRECTORY $<TARGET_FILE_DIR:validation>
	)

set(PLOTSCRIPT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/validationPlots.py")
add_custom_command(TARGET validation POST_BUILD
	COMMENT "Copying plotting script ${PLOTSCRIPT_PATH}."	
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PLOTSCRIPT_PATH} validationPlots.py
	WORKING_DIRECTORY $<TARGET_FILE_DIR:validation>
	)

add_custom_command(TARGET validation POST_BUILD
	COMMENT "Copying plotting TG195 results ${CMAKE_CURRENT_SOURCE_DIR}/data/TG195results.txt."	
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/TG195results.txt" TG195results.txt
	WORKING_DIRECTORY $<TARGET_FILE_DIR:validation>
	)

set(PLOTSCRIPTSCATTER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/scatterPlots.py")
add_custom_command(TARGET scattervalidation POST_BUILD
	COMMENT "Copying plotting script ${PLOTSCRIPTSCATTER_PATH}."	
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PLOTSCRIPTSCATTER_PATH} scatterPlots.py
	WORKING_DIRECTORY $<TARGET_FILE_DIR:scattervalidation>
	)



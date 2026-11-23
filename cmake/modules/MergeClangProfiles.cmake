# Merges the raw Clang profiles (*.profraw) in PGO_DIR into PGO_DIR/default.profdata.
# Usage: cmake -DPROFDATA=<llvm-profdata> -DPGO_DIR=<dir> -P MergeClangProfiles.cmake
file(GLOB_RECURSE RAW_PROFILES "${PGO_DIR}/*.profraw")
if(NOT RAW_PROFILES)
    message(FATAL_ERROR "No raw profiles (*.profraw) found in ${PGO_DIR}, did the training run produce any?")
endif()
execute_process(
    COMMAND "${PROFDATA}" merge "--output=${PGO_DIR}/default.profdata" ${RAW_PROFILES}
    RESULT_VARIABLE MERGE_RESULT)
if(NOT MERGE_RESULT EQUAL 0)
    message(FATAL_ERROR "llvm-profdata merge failed")
endif()
list(LENGTH RAW_PROFILES N_PROFILES)
message(STATUS "Merged ${N_PROFILES} raw profile(s) into ${PGO_DIR}/default.profdata")

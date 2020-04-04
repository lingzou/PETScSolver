##########################################################################
# Step 1. Run the input file and check if it executes successfully
#
EXECUTE_PROCESS(
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMAND ${TEST_PROG} input/${TEST_NAME}.i
  RESULT_VARIABLE EXE_HAD_ERROR
)

if(EXE_HAD_ERROR)
	message(FATAL_ERROR "Test failed to run, i.e., exit value != 0")
endif()

##########################################################################
# Step 2. Check if the output files are (numerically) the same as expected
#         Will check two time stamps
# 2.1 First time stamp
message("Now check the first time stamp at N = ${STEP_1}")
EXECUTE_PROCESS(
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_SOURCE_DIR}/vtk_diff.py output/${TEST_NAME}_step_${STEP_1}.vtk output/expected/${TEST_NAME}_step_${STEP_1}.vtk
  RESULT_VARIABLE FILES_ARE_DIFF
)

if(FILES_ARE_DIFF)
	message(FATAL_ERROR "TEST FAILED - files are different.")
endif()

# 2.2 Second time stamp
message("") # an empty line
message("Now check the second time stamp at N = ${STEP_2}")
EXECUTE_PROCESS(
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_SOURCE_DIR}/vtk_diff.py output/${TEST_NAME}_step_${STEP_2}.vtk output/expected/${TEST_NAME}_step_${STEP_2}.vtk
  RESULT_VARIABLE FILES_ARE_DIFF
)

if(FILES_ARE_DIFF)
	message(FATAL_ERROR "TEST FAILED - files are different.")
endif()

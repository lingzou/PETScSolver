#######################################################################################
# The original CMakeList.txt file is downloaded from:
# https://github.com/jedbrown/dohp/blob/master/CMakeLists.txt
#######################################################################################

SET (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")

# Change this line to your gcc location
SET (CMAKE_C_COMPILER             "/usr/bin/gcc")
#SET (CMAKE_C_COMPILER             "/usr/local/bin/gcc")
#SET (CMAKE_C_FLAGS                "-Wall -std=c99")
#SET (CMAKE_C_FLAGS_DEBUG          "-g")
#SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

# Change this line to your g++ location
SET (CMAKE_CXX_COMPILER             "/usr/bin/g++")
#SET (CMAKE_CXX_COMPILER             "/usr/local/bin/g++")
#SET (CMAKE_CXX_FLAGS                "-Wall")
#SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
#SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
#SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")


cmake_minimum_required (VERSION 3.10)
project (PETScSolver)

# Get root directory path, show it
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
#message ("The root dir is: " ${CMAKE_SOURCE_DIR})

# Normally PETSc is built with MPI, if not, use CC=mpicc, etc
find_package (PETSc REQUIRED)

include (CheckCSourceCompiles) # The name is misleading, this also tries to link
check_c_source_compiles ("
#define PragmaQuote(a) _Pragma(#a)
#define PragmaGCC(a) PragmaQuote(GCC a)
int main(int argc,char *argv[]) {
  PragmaGCC(diagnostic ignored \"-Wconversion\")
  char c = (int)argv[0][0] + argv[argc-1][0];
  return c;
}" dHAVE_PRAGMA_GCC)

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories (${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})
# additional include(s)
include_directories (${CMAKE_SOURCE_DIR})
include_directories (${CMAKE_SOURCE_DIR}/problems)
include_directories (${CMAKE_SOURCE_DIR}/utils)
include_directories (${CMAKE_SOURCE_DIR}/GetPot)

set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Lists of all source file subdirectories
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR} SRC_ROOT)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/problems SRC_PROBLEMS)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/utils SRC_UTILS)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/GetPot SRC_GETPOT)
#message ("SRC_ROOT" ${SRC_ROOT})
#message ("SRC_PROBLEMS" ${SRC_PROBLEMS})

ADD_EXECUTABLE (PETScSolver ${SRC_ROOT} ${SRC_PROBLEMS} ${SRC_UTILS} ${SRC_GETPOT})
TARGET_LINK_LIBRARIES(PETScSolver ${PETSC_LIBRARIES})

# CTEST
ENABLE_TESTING()

set (ALL_TESTS HeatCond1D SodShockTube LaxRiemannProblem WaterFaucet_5E2P_SG_BDF1_1st)
# The following are for each test in the format of a LIST: (test_name, 1st_step_to_test, 2nd_step_to_test)
#   (Name                           TestName;           step; step)
set (HeatCond1D                     HeatCond1D;         20;   40)
set (SodShockTube                   SodShockTube;       20;   40)
set (LaxRiemannProblem              LaxRiemannProblem;  80;   160)
set (WaterFaucet_5E2P_SG_BDF1_1st   WaterFaucet_5E2P_SG_BDF1_1st;       20;   40)

foreach(TEST IN LISTS ALL_TESTS)
  LIST(GET ${TEST} 0 TEST_NAME)
  LIST(GET ${TEST} 1 STEP_VAL_1)
  LIST(GET ${TEST} 2 STEP_VAL_2)
  ADD_TEST( NAME ${TEST_NAME}
            COMMAND ${CMAKE_COMMAND}
            -DTEST_PROG=$<TARGET_FILE:PETScSolver>
            -DSOURCEDIR=${CMAKE_SOURCE_DIR}
            -DTEST_NAME=${TEST}
            -DSTEP_1=${STEP_VAL_1}
            -DSTEP_2=${STEP_VAL_2}
            -P runtest.cmake)
endforeach()

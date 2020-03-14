#######################################################################################
# The original CMakeList.txt file is downloaded from:
# https://github.com/jedbrown/dohp/blob/master/CMakeLists.txt
#######################################################################################

SET (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")

# Change this line to your gcc location
SET (CMAKE_C_COMPILER             "/opt/moose/gcc-9.2.0/bin/gcc")
#SET (CMAKE_C_COMPILER             "/usr/local/bin/gcc")
#SET (CMAKE_C_FLAGS                "-Wall -std=c99")
#SET (CMAKE_C_FLAGS_DEBUG          "-g")
#SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

# Change this line to your g++ location
SET (CMAKE_CXX_COMPILER             "/opt/moose/gcc-9.2.0/bin/g++")
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

set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Lists of all source file subdirectories
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR} SRC_ROOT)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/problems SRC_PROBLEMS)
#message ("SRC_ROOT" ${SRC_ROOT})
#message ("SRC_PROBLEMS" ${SRC_PROBLEMS})

ADD_EXECUTABLE (PETScSolver ${SRC_ROOT} ${SRC_PROBLEMS})
TARGET_LINK_LIBRARIES(PETScSolver ${PETSC_LIBRARIES})

# PETSc for PDEs
*PETSc as a solver for non-linear PDEs*

"PETSc is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations" [See: https://www.mcs.anl.gov/petsc/]

## Compiling
This assumes that you have PETSc correctly installed with your specified gcc and g++ compilers.

  * Goto CMakeLists.txt and change line 9 and 17 for your gcc and g++ compilers.
  * Goto FindPETSc.cmake and change line 26/27 for PETSC_DIR and PETSC_ARCH.
  * At root directory, type `cmake .`
  * Then, type `make` or  `make -j <n>`

## Running an input file
PETSc option style is used for input file, e.g., `-option_name option_value`. See, for example, `input/HeatCond1D.i`

  * Type `./PETScSolver <input_file_name>` to run an input file. For example, `./PETScSolver input/HeatCond1D.i`

## Regression test
CMake-style testing (ctest) is used for regression test

  * Type `ctest` or `ctest -j <n>`

#!/bin/bash

# run
./PETScSolver input/HeatCond1D.i > run_test_output

# do diff for each case
temp=$(diff output/Solution_step_20.vtk output/expected/Solution_step_20.vtk)
error=$?
if [ $error -eq 0 ]
then
  printf "................%24b" "\e[0;32m[Step 20; OK]\e[m" "\n"
else
  printf "%24b........" "\e[1;31m[Step 20; DIFF ERROR]\e[m" "\n"
fi

temp=$(diff output/Solution_step_40.vtk output/expected/Solution_step_40.vtk)
error=$?
if [ $error -eq 0 ]
then
  printf "................%24b" "\e[0;32m[Step 40; OK]\e[m" "\n"
else
  printf "%24b........" "\e[1;31m[Step 40; DIFF ERROR]\e[m" "\n"
fi

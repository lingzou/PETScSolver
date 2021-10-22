#include <iostream>
#include <petscsnes.h>

#include "PETScProblemInterface.h"
#include "InputParser.h"
#include <GetPot>

int main(int argc, char **argv)
{
  /*
   * Input file parser
   */
  if (argc < 2)
    sysError("Invalid input command line.\n"
             "Please use: './PETScSolver <input_file_name> [PETSc options]'\n"
             "<required> [optional]");
  InputParser input_parser(argv[1]);

  // PETSc application starts with PetscInitialize
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  /*
   * Initialize PETSc App
   */
  ApplicationCtx AppCtx;
  AppCtx.initializePETScApp(&input_parser);

  /*
   *  Setup PETSc work space
   */
  AppCtx.setupPETScWorkSpace();

  /*
   *  Setup initial conditions for the application
   */
  AppCtx.setupInitialConditions();

  /*
   *  Solving
   */
  InputParameterList& globalParamList = input_parser.getGlobalParamList();
  TimeScheme ts = globalParamList.getParameterValue<TimeScheme>("ts");
  int N_Steps = globalParamList.getParameterValue<int>("n_steps");
  for (unsigned int step = 1; step <= N_Steps; step++)
  {
    // 1. Before PETSc Solving
    // 1.1 (If applicable) preparing CN old time step RHS (only for the very first time step)
    if ((ts == CN) && (step == 1))
    {
      PetscReal * res_RHS_old;
      VecGetArray(AppCtx.res_RHS_old, &res_RHS_old);
      AppCtx.myProblemSystem->RHS(res_RHS_old);
      VecRestoreArray(AppCtx.res_RHS_old, &res_RHS_old);
    }

    // 2. PETSc solving: solve the system by iterating the solution vector
    PetscPrintf(PETSC_COMM_WORLD, "Time step = %d, dt = %g\n", step, AppCtx.myProblemSystem->getDt());
    SNESSolve(AppCtx.snes, NULL, AppCtx.u);

    // 3. After PETSc solving: March time forward; write solutions; update old solutions etc.
    AppCtx.myProblemSystem->onTimestepEnd();

    // 3.1. (If applicable) save the NEW RHS to OLD RHS for the next time step CN
    if (ts == CN)   VecCopy(AppCtx.res_RHS, AppCtx.res_RHS_old);

    // 3.2. Print some additional info.
    double current_t = AppCtx.myProblemSystem->getCurrentTime();
    PetscPrintf(PETSC_COMM_WORLD, "End of Time step = %d, current time = %g\n\n", step, current_t);
  }

  /*
   *  Free work space
   */
  AppCtx.FreePETScWorkSpace();

  // All PETSc applications have to end with calling PetscFinalize()
  PetscFinalize();

  return 0;
}

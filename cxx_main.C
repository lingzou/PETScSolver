#include <iostream>
#include <petscsnes.h>
#include <petsc.h>

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
  bool converged = false;
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
    VecCopy(AppCtx.u, AppCtx.u_backup); // backup the current u, in case to retry a failed time step

    // 2. PETSc solving: solve the system by iterating the solution vector
    PetscPrintf(PETSC_COMM_WORLD, "Time step = %d, dt = %g\n", step, AppCtx.myProblemSystem->getDt());
    converged = false; // reset converged to be false at the beginning of every time step
    while(!converged)
    {
      SNESSolve(AppCtx.snes, NULL, AppCtx.u);
      SNESConvergedReason snes_converged_reason;
      SNESGetConvergedReason(AppCtx.snes, &snes_converged_reason);  // https://petsc.org/main/docs/manualpages/SNES/SNESGetConvergedReason.html
      converged = (snes_converged_reason > 0); // https://petsc.org/main/docs/manualpages/SNES/SNESConvergedReason.html#SNESConvergedReason

      if (!converged) // PETSc solving not converged
      {
        AppCtx.myProblemSystem->adjustTimeStepSize(0.8);
        PetscPrintf(PETSC_COMM_WORLD, "  Did not converge, try a smaller dt = %g.\n", AppCtx.myProblemSystem->getDt());
        if (AppCtx.myProblemSystem->getDt() < globalParamList.getParameterValue<double>("dt_min"))
        {
          PetscPrintf(PETSC_COMM_WORLD, "\nERROR: dt too small, < %g (dt_min).\n", globalParamList.getParameterValue<double>("dt_min"));
          break; // break the while loop
        }
        VecCopy(AppCtx.u_backup, AppCtx.u);
      }
    }
    if (!converged) break; // failed, break the for N_Steps loop

    // 3. After PETSc solving: March time forward; write solutions; update old solutions etc.
    AppCtx.myProblemSystem->onTimestepEnd();

    // 3.1. (If applicable) save the NEW RHS to OLD RHS for the next time step CN
    if (ts == CN)   VecCopy(AppCtx.res_RHS, AppCtx.res_RHS_old);

    // 3.2. Print some additional info.
    double current_t = AppCtx.myProblemSystem->getCurrentTime();
    PetscPrintf(PETSC_COMM_WORLD, "End of Time step = %d, current time = %g\n\n", step, current_t);

    // 3.3 Now it is the right time to adjust time step size back
    AppCtx.myProblemSystem->adjustTimeStepSize(1.25);
  }

  /*
   *  Free work space
   */
  AppCtx.FreePETScWorkSpace();

  // All PETSc applications have to end with calling PetscFinalize()
  PetscFinalize();

  if (!converged) { std::cout << "\n\n*** Simulation Failed. ***\n";  return 1; }
  return 0;
}

#include <iostream>
#include <petscsnes.h>

#include "PETScProblemInterface.h"
#include "HeatConduction1D.h"

int main(int argc, char **argv)
{
  // Require an input file
  if (argc != 2)
  {
    std::cerr << "Invalid input command line." << std::endl;
    std::cerr << "Please use: './PETScSolver <input_file_name>" << std::endl;
    exit(1);
  }

  ApplicationCtx AppCtx;

  // PETSc application starts with PetscInitialize
  PetscInitialize(&argc, &argv, argv[1], PETSC_NULL);
  //  MPI_Comm_size(PETSC_COMM_WORLD, &(AppCtx.size));
  //  if (AppCtx.size != 1)
  //    SETERRQ(PETSC_COMM_SELF, 1, "This is a uniprocessor example only!");

  double t_start = 0.0;

  double dt          = PetscOptionsGetRequiredReal("-dt");
  std::string ts_str = PetscOptionsGetRequiredString("-ts");
  TimeScheme ts = StringToEnum(ts_str);

  AppCtx.myPETScProblem = new HeatConduction1D(ts, t_start, dt);

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
  int N_Steps = 40;
  for (unsigned int step = 1; step <= N_Steps; step++)
  {
    // 1. Before PETSc Solving
    // 1.1 (If applicable) preparing CN old time step RHS (only for the very first time step)
    if ((ts == CN) && (step == 1))
    {
      PetscReal * res_RHS_old;
      VecGetArray(AppCtx.res_RHS_old, &res_RHS_old);
      AppCtx.myPETScProblem->RHS(res_RHS_old);
      VecRestoreArray(AppCtx.res_RHS_old, &res_RHS_old);
    }

    // 2. PETSc solving
    PetscPrintf(PETSC_COMM_WORLD, "Time step = %d, dt = %g\n", step, dt);
    SNESSolve(AppCtx.snes, NULL, AppCtx.u);

    // 3. After PETSc solving
    //    March time forward; write solutions; update old solutions etc.
    AppCtx.myPETScProblem->onTimestepEnd();

    // 3.1. (If applicable) save the NEW RHS to OLD RHS for the next time step CN
    if (ts == CN)
      VecCopy(AppCtx.res_RHS, AppCtx.res_RHS_old);

    // 3.2. Print some additional info.
    double current_t = AppCtx.myPETScProblem->getCurrentTime();
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

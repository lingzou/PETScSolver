#include <iostream>
#include <petscsnes.h>

#include "HeatConduction1D.h"

PetscErrorCode SNESFormFunction(SNES, Vec, Vec, void*);
PetscErrorCode SNESMonitor(SNES, PetscInt, PetscReal, void*);
PetscErrorCode KSPMonitor(KSP, PetscInt, PetscReal, void*);

struct ApplicationCtx
{
  HeatConduction1D * myPETScProblem;
  PetscInt        N_DOFs;     // Number of degrees of freedom

  SNES            snes;       // SNES
  KSP             ksp;        // KSP

  MatFDColoring   fdcoloring;   // MatFDColoring
  Mat             J_Mat;        // Jacobian matrix
  Mat             P_Mat;        // Preconditioning matrix
  Mat             J_MatrixFree; // Jacobian-free

  Vec             u;          // unknown vector
  Vec             u_old;      // unknown vector old
  Vec             u_oldold;   // unknown vector oldold

  Vec             r;               // total residual = res_transient - res_RHS
  Vec             res_transient;   /* residual from transient term */
  Vec             res_RHS;         /* residual from RHS */

  void setupPETScWorkSpace();
  void setupInitialConditions();
  void FreePETScWorkSpace();
  void setupMatrices();
};

void
ApplicationCtx::setupPETScWorkSpace()
{
  // Get total number of DOF
  N_DOFs = myPETScProblem->getNDOF();

  // Prepare NULL matrix
  J_Mat = NULL;
  J_MatrixFree = NULL;
  P_Mat = NULL;

  // Prepare PETSc vectors
  VecCreate(PETSC_COMM_SELF, &u);
  VecSetSizes(u, PETSC_DECIDE, N_DOFs);
  VecSetFromOptions(u);
  VecDuplicate(u, &u_old);
  VecDuplicate(u, &u_oldold);
  VecDuplicate(u, &r);
  VecDuplicate(u, &res_transient);
  VecDuplicate(u, &res_RHS);

  // Setup SNES
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetFunction(snes, r, SNESFormFunction, (void*)(this));
  SNESMonitorSet(snes, SNESMonitor, (void*)(this), NULL);

  // Setup Matrix
  setupMatrices();

  // Setup KSP
  PetscReal ksp_rtol = 1.0e-3;
  PetscInt ksp_maxits = 30;
  SNESGetKSP(snes, &ksp);
  KSPSetTolerances(ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, ksp_maxits);

  // Finalize SNES setup
  SNESSetFromOptions(snes);
}

void
ApplicationCtx::setupMatrices()
{
  bool fd_coloring = true;
  if(fd_coloring)
  {
    // Create Matrix-free context
    MatCreateSNESMF(snes, &J_MatrixFree);

    // Let the problem setup Jacobian matrix sparsity
    myPETScProblem->FillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc examples:
    // https://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex14.c.html
    // https://www.mcs.anl.gov/petsc/petsc-current/src/mat/examples/tutorials/ex16.c.html
    ISColoring      iscoloring;
    MatColoring     mc;
    MatColoringCreate(P_Mat, &mc);
    MatColoringSetType(mc, MATCOLORINGSL);
    MatColoringSetFromOptions(mc);
    MatColoringApply(mc, &iscoloring);
    MatColoringDestroy(&mc);
    MatFDColoringCreate(P_Mat, iscoloring, &fdcoloring);
    MatFDColoringSetFunction(fdcoloring, (PetscErrorCode (*)(void))SNESFormFunction, this);
    MatFDColoringSetFromOptions(fdcoloring);
    MatFDColoringSetUp(P_Mat, iscoloring, fdcoloring);
    ISColoringDestroy(&iscoloring);

    SNESSetJacobian(snes,           // snes
                    J_MatrixFree,   // Jacobian-free
                    P_Mat,          // Preconditioning matrix
                    SNESComputeJacobianDefaultColor,  // Use finite differencing and coloring
                    fdcoloring);    // fdcoloring
  }
  else
  {
    // See PETSc example:
    // https://www.mcs.anl.gov/petsc/petsc-current/src/ts/examples/tutorials/ex10.c.html
    MatCreateSeqAIJ(PETSC_COMM_SELF, N_DOFs, N_DOFs, PETSC_DEFAULT, PETSC_NULL, &J_Mat);
    SNESSetJacobian(snes,     // snes
                    J_Mat,    // Jacobian matrix
                    J_Mat,    // Preconditioning mat, use the same Jacobian mat
                    SNESComputeJacobianDefault,
                    PETSC_NULL);
  }
}

void
ApplicationCtx::setupInitialConditions()
{
  PetscScalar *uu;
  VecGetArray(u, &uu);
  myPETScProblem->SetupInitialCondition(uu);
  VecRestoreArray(u, &uu);
}

void
ApplicationCtx::FreePETScWorkSpace()
{
  // Destroy PETSc vectors
  VecDestroy(&u);
  VecDestroy(&u_old);
  VecDestroy(&u_oldold);
  VecDestroy(&r);
  VecDestroy(&res_transient);
  VecDestroy(&res_RHS);

  // Destroy PETSc matrix
  if (J_Mat != NULL)        MatDestroy(&J_Mat);
  if (J_MatrixFree != NULL) MatDestroy(&J_MatrixFree);
  if (P_Mat != NULL)        MatDestroy(&P_Mat);

  // Destroy SNES
  SNESDestroy(&snes);

  // Destroy MatFDColoring
  MatFDColoringDestroy(&fdcoloring);

  delete myPETScProblem;
}

PetscErrorCode SNESFormFunction(SNES snes, Vec u, Vec f, void * AppCtx)
{
  ApplicationCtx * appCtx = (ApplicationCtx *) AppCtx;
  HeatConduction1D * myProblem = appCtx->myPETScProblem;

  // get vectors
  PetscScalar *uu, *res_tran, *res_RHS;
  VecGetArray(u, &uu);
  VecGetArray(appCtx->res_transient, &res_tran);
  VecGetArray(appCtx->res_RHS, &res_RHS);

  // use the most updated solution vector to update solution, to compute RHS and transient residuals
  myProblem->updateSolution(uu, NEW);
  myProblem->RHS(res_RHS);
  myProblem->transientResidual(res_tran);

  // restore vectors
  VecRestoreArray(u, &uu);
  VecRestoreArray(appCtx->res_transient, &res_tran);
  VecRestoreArray(appCtx->res_RHS, &res_RHS);

  // assemble final residuals: f = transient - RHS
  VecWAXPY(f, -1.0, appCtx->res_RHS, appCtx->res_transient);

  return 0;
}

PetscErrorCode SNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, void* /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "        Non-linear Step = %2D, fnorm = %12.5E\n", its, fnorm);
  return 0;
}

PetscErrorCode KSPMonitor(KSP ksp, PetscInt its, PetscReal rnorm, void* /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "          Linear step = %2D, rnorm = %12.5E\n", its, rnorm);
  return 0;
}

int main(int argc, char **argv)
{
  ApplicationCtx AppCtx;

  // PETSc application starts with PetscInitialize
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  //  MPI_Comm_size(PETSC_COMM_WORLD, &(AppCtx.size));
  //  if (AppCtx.size != 1)
  //    SETERRQ(PETSC_COMM_SELF, 1, "This is a uniprocessor example only!");

  AppCtx.myPETScProblem = new HeatConduction1D();

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
  //std::string time_scheme = "CN";  // "BDF1", "BDF2"
  double t_start = 0.0;
  double dt = 1.0;
  int N_Steps = 10;

  for (int step = 1; step <= N_Steps; step++)
  {
    // PETSc solving
    PetscPrintf(PETSC_COMM_WORLD, "Solving time step %d, dt = %g. \n", step, dt);
    SNESSolve(AppCtx.snes, NULL, AppCtx.u);

    // After solving each time step, the NEW solution now becomes the OLD solution
    VecCopy(AppCtx.u, AppCtx.u_old);
    PetscScalar *uu;
    VecGetArray(AppCtx.u, &uu);
    AppCtx.myPETScProblem->updateSolution(uu, OLD);
    VecRestoreArray(AppCtx.u, &uu);
  }
  AppCtx.myPETScProblem->writeSolution();

  /*
   *  Free work space
   */
  AppCtx.FreePETScWorkSpace();

  // All PETSc applications have to end with calling PetscFinalize()
  PetscFinalize();

  return 0;
}

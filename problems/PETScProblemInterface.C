#include <iostream>
#include "PETScProblemInterface.h"
#include "HeatConduction1D.h"
#include "EulerEquation1D.h"

void
ApplicationCtx::initializePETScApp()
{
  std::string problem_name = PetscOptionsGetRequiredString("-problem");
  if (problem_name.compare("HeatConduction1D") == 0)
    myPETScProblem = new HeatConduction1D();
  else if (problem_name.compare("EulerEquation1D") == 0)
    myPETScProblem = new EulerEquation1D();
  else
  {
    std::cerr << "ERROR: UNKNOWN problem: '" << problem_name << "'." << std::endl;
    exit(1);
  }

  // Get total number of DOF
  N_DOFs = myPETScProblem->getNDOF();
}

void
ApplicationCtx::setupPETScWorkSpace()
{
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
  VecDuplicate(u, &res_RHS_old);

  // Setup SNES
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetFunction(snes, r, SNESFormFunction, (void*)(this));
  SNESMonitorSet(snes, SNESMonitor, (void*)(this), NULL);

  // Setup Matrix
  setupMatrices();

  // Setup KSP
  PetscReal ksp_rtol = 1.0e-3;  // Hard-coded value
  PetscInt ksp_maxits = 30;     // Hard-coded value
  SNESGetKSP(snes, &ksp);
  KSPSetTolerances(ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, ksp_maxits);

  // Finalize SNES setup
  SNESSetFromOptions(snes);
}

void
ApplicationCtx::setupMatrices()
{
  bool exact_jacobian = false;     // True if want to use hand-calculated Jacobian
  hasFDColoring = PETSC_TRUE;      // True if finite differencing Jacobian
  if(exact_jacobian)
  {
    // Create Matrix-free context
    // MatCreateSNESMF(snes, &J_MatrixFree); Why Jacobian-free not working here?

    // Let the problem setup Jacobian matrix sparsity
    myPETScProblem->FillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc example:
    // https://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex1.c.html
    // Use hand-calculated Jacobian
    SNESSetJacobian(snes, P_Mat, P_Mat, FormJacobian, this);
  }
  else if(hasFDColoring)
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
  VecDestroy(&res_RHS_old);

  // Destroy PETSc matrix
  if (J_Mat != NULL)        MatDestroy(&J_Mat);
  if (J_MatrixFree != NULL) MatDestroy(&J_MatrixFree);
  if (P_Mat != NULL)        MatDestroy(&P_Mat);

  // Destroy SNES
  SNESDestroy(&snes);

  // Destroy MatFDColoring
  if (hasFDColoring)
    MatFDColoringDestroy(&fdcoloring);

  delete myPETScProblem;
}

PetscErrorCode SNESFormFunction(SNES snes, Vec u, Vec f, void * AppCtx)
{
  ApplicationCtx * appCtx = (ApplicationCtx *) AppCtx;
  PETScProblem * myProblem = appCtx->myPETScProblem;
  TimeScheme ts = myProblem->getTimeScheme();

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

  switch (ts)
  {
    case BDF1:
    case BDF2:
      // assemble final residuals: f = transient - RHS
      VecWAXPY(f, -1.0, appCtx->res_RHS, appCtx->res_transient);
      break;

    case CN:
      // assemble final residuals: f = transient - 0.5(RHS + RHS_old)
      VecAXPBYPCZ(f, -0.5, -0.5, 0, appCtx->res_RHS, appCtx->res_RHS_old); // f = -0.5(RHS + RHS_old)
      VecAXPY(f, 1, appCtx->res_transient); // f = f + transient = transient - 0.5(RHS + RHS_old)
      break;

    defaut:
      std::cerr << "ERROR: not implemented." << std::endl;
      exit(1);
  }

  return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec u, Mat jac, Mat B, void * AppCtx)
{
  ApplicationCtx * appCtx = (ApplicationCtx *) AppCtx;
  PETScProblem * myProblem = appCtx->myPETScProblem;
  myProblem->computeJacobianMatrix(B);

  return 0;
}

PetscErrorCode SNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, void* /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "    NL Step = %2D, fnorm = %12.5E\n", its, fnorm);
  return 0;
}

PetscErrorCode KSPMonitor(KSP ksp, PetscInt its, PetscReal rnorm, void* /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "      Linear step = %2D, rnorm = %12.5E\n", its, rnorm);
  return 0;
}

double PetscOptionsGetRequiredReal(std::string name)
{
  PetscBool hasInput = PETSC_FALSE;
  double value = 0.0;

  PetscOptionsGetReal(NULL, NULL, name.c_str(), &value, &hasInput);

  if (!hasInput)
  {
    std::cerr << "Required PETSc <Real> input '" << name << "' is not found." << std::endl;
    exit(1);
  }

  return value;
}

std::string PetscOptionsGetRequiredString(std::string name)
{
  PetscBool hasInput = PETSC_FALSE;
  char value[4096];
  PetscOptionsGetString(NULL, NULL, name.c_str(), value, 4096, &hasInput);

  if (!hasInput)
  {
    std::cerr << "Required PETSc <String> input '" << name << "' is not found." << std::endl;
    exit(1);
  }

  return std::string(value);
}

#include <iostream>
#include "PETScProblemInterface.h"
#include "utils.h"
#include "ParameterList.h"

#include "HeatConduction1D.h"
#include "EulerEquation1D.h"
#include "FiveEqnTwoP_StagGrid.h"
#include "SinglePhaseFlow.h"

void
ApplicationCtx::initializePETScApp(InputParser* input_parser)
{
  _input_parser = input_parser;
  myProblemSystem = new ProblemSystem(input_parser);

  // Insert PETSc options provided from input file
  std::string option = _input_parser->getGlobalParamList().getParameterValue<std::string>("petsc_options");
  PetscOptionsInsertString(NULL, option.c_str());

  // Get total number of DOF
  N_DOFs = myProblemSystem->getNDOF();
}

void
ApplicationCtx::setupPETScWorkSpace()
{
  // Prepare NULL matrix
  J_Mat = NULL;
  J_MatrixFree = NULL;
  P_Mat = NULL;

  // Prepare NULL MatFDColoring
  fdcoloring = NULL;

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
  PetscReal ksp_rtol = _input_parser->getGlobalParamList().getParameterValue<double>("linear_rtol");
  PetscInt ksp_maxits = _input_parser->getGlobalParamList().getParameterValue<int>("linear_max_its");
  SNESGetKSP(snes, &ksp);
  KSPSetTolerances(ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, ksp_maxits);

  // Finalize SNES setup
  SNESSetFromOptions(snes);
}

void
ApplicationCtx::setupMatrices()
{
  int solver_option = _input_parser->getGlobalParamList().getParameterValue<int>("solver_option");

  if(solver_option == 0) // Newton's method + Hand-coded "(hopefully) exact" Jacobian
  {
    // Let the problem setup Jacobian matrix sparsity
    myProblemSystem->FillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc example:
    // https://petsc.org/release/src/snes/tutorials/ex1.c.html
    // Use hand-calculated Jacobian as both the Jacobian and Preconditioning Jacobian
    //SNESSetJacobian(snes, J_MatrixFree, P_Mat, FormJacobian, this);
    SNESSetJacobian(snes, P_Mat, P_Mat, FormJacobian, this);
  }
  else if(solver_option == 1) // Matrix-free + Hand-coded Jacobian as the Preconditioning Jacobian
  {
    // Create Matrix-free context
    MatCreateSNESMF(snes, &J_MatrixFree);

    // Let the problem setup Jacobian matrix sparsity
    myProblemSystem->FillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc example:
    // https://petsc.org/release/src/ts/tutorials/ex15.c.html
    // Use hand-calculated Jacobian as the Preconditioning Jacobian
    SNESSetJacobian(snes, J_MatrixFree, P_Mat, FormJacobian, this);
  }
  else if(solver_option == 2) // Matrix-free + Finite-differencing Preconditioning Jacobian (using coloring)
  {
    // Create Matrix-free context
    MatCreateSNESMF(snes, &J_MatrixFree);

    // Let the problem setup Jacobian matrix sparsity
    myProblemSystem->FillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc examples:
    // https://petsc.org/release/src/snes/tutorials/ex14.c.html
    // https://petsc.org/release/src/mat/tutorials/ex16.c.html
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
  else if (solver_option == 3) // Finite-differencing, no coloring, slowest
  {
    // See PETSc example:
    // https://petsc.org/release/src/ts/tutorials/ex10.c.html
    MatCreateSeqAIJ(PETSC_COMM_SELF, N_DOFs, N_DOFs, PETSC_DEFAULT, PETSC_NULL, &J_Mat);
    SNESSetJacobian(snes,     // snes
                    J_Mat,    // Jacobian matrix
                    J_Mat,    // Preconditioning mat, use the same Jacobian mat
                    SNESComputeJacobianDefault,
                    PETSC_NULL);
  }
  else  sysError("Unknown Jacobian option.");
}

void
ApplicationCtx::setupInitialConditions()
{
  PetscScalar *uu;
  VecGetArray(u, &uu);
  myProblemSystem->SetupInitialCondition(uu);
  VecRestoreArray(u, &uu);

  VecCopy(u, u_old);
  VecCopy(u, u_oldold);

  // Write output for t = 0
  myProblemSystem->writeOutput(0);
}

void
ApplicationCtx::FreePETScWorkSpace()
{
  _input_parser->checkUnusedVariables();
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
  if (fdcoloring != NULL)   MatFDColoringDestroy(&fdcoloring);

  delete myProblemSystem;
}

PetscErrorCode SNESFormFunction(SNES snes, Vec u, Vec f, void * AppCtx)
{
  ApplicationCtx * appCtx = (ApplicationCtx *) AppCtx;
  ProblemSystem * myProblemSystem = appCtx->myProblemSystem;
  TimeScheme ts = myProblemSystem->getTimeScheme();

  // get vectors
  PetscScalar *uu, *res_tran, *res_RHS;
  VecGetArray(u, &uu);
  VecGetArray(appCtx->res_transient, &res_tran);
  VecGetArray(appCtx->res_RHS, &res_RHS);

  // use the most updated solution vector to update solution, to compute RHS and transient residuals
  myProblemSystem->updateSolution(uu);
  myProblemSystem->RHS(res_RHS);
  myProblemSystem->transientResidual(res_tran);

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

    case INVALID:
    default:
      sysError("ERROR: not implemented.");
  }

  return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec u, Mat jac, Mat B, void * AppCtx)
{
  ApplicationCtx * appCtx = (ApplicationCtx *) AppCtx;
  ProblemSystem * myProblemSystem = appCtx->myProblemSystem;
  myProblemSystem->computeJacobianMatrix(B);

  if (jac != B)
  {
    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
  }

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

#ifndef PETSC_PROBLEM_INTERFACE_H
#define PETSC_PROBLEM_INTERFACE_H

#include "PETScProblem.h"

PetscErrorCode SNESFormFunction(SNES, Vec, Vec, void*);
PetscErrorCode SNESMonitor(SNES, PetscInt, PetscReal, void*);
PetscErrorCode KSPMonitor(KSP, PetscInt, PetscReal, void*);
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);

// Wrap up several PETSc built-in function with additional check;
//   e.g., PETSc built-in function 'PetscOptionsGetReal' -> 'PetscOptionsGetRequiredReal'
double      PetscOptionsGetRequiredReal(std::string name);
std::string PetscOptionsGetRequiredString(std::string name); /* max string length 4096 */

struct ApplicationCtx
{
  PETScProblem *  myPETScProblem;
  PetscInt        N_DOFs;     // Number of degrees of freedom

  SNES            snes;       // SNES
  KSP             ksp;        // KSP

  PetscBool       hasFDColoring;
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
  Vec             res_RHS_old;     /* residual from RHS (old time step)*/

  void initializePETScApp();
  void setupPETScWorkSpace();
  void setupInitialConditions();
  void FreePETScWorkSpace();
  void setupMatrices();
};

#endif /*PETSC_PROBLEM_INTERFACE_H*/

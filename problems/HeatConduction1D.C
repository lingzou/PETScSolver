#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "HeatConduction1D.h"

HeatConduction1D::HeatConduction1D()
{
  // Using Finite Difference Method for heat conduction problem
  // 1/alpha * dT/dt = Laplacian(T) + q/k

  // 0 <= x <= 1
  // At steady state, let T(x) = sin(pi * x)
  // such that: q(x) = k * pi * pi * sin(pi * x)

  _time_scheme = BDF1;
  _dt = 1.0;

  length = 1.0;
  n_Cell = 40;
  k = 1.0;
  alpha = 1.0;

  n_Node = n_Cell + 1;
  n_DOFs = n_Node;

  dx = length / n_Cell;
  dx2 = dx * dx;

  x.resize(n_DOFs);
  T.resize(n_DOFs);
  T_old.resize(n_DOFs);
  T_oldold.resize(n_DOFs);

  // locations of nodes
  for(unsigned int i = 0; i < x.size(); i++)
    x[i] = dx * i;
}

HeatConduction1D::~HeatConduction1D()
{
}

void
HeatConduction1D::SetupInitialCondition(double * u)
{
  for(unsigned int i = 0; i < n_DOFs; i++)
  {
    u[i] = 0.0;
    T[i] = 0.0; T_old[i] = 0.0; T_oldold[i] = 0.0;
  }
}

void
HeatConduction1D::updateSolution(double * u, TimeStepIndex index)
{
  switch (index)
  {
    case NEW:
      for(unsigned int i = 0; i < T.size(); i++)
        T[i] = u[i];
      break;

    case OLD:
      for(unsigned int i = 0; i < T_old.size(); i++)
        T_old[i] = u[i];
      break;

    case OLDOLD:
      for(unsigned int i = 0; i < T_oldold.size(); i++)
        T_oldold[i] = u[i];
      break;

    default:
      std::cerr << "ERROR.\n";
      exit(1);
      break;
  }
}

void
HeatConduction1D::transientResidual(double * res)
{
  // hard-coded BC; T_left = 0; T_right = 0
  res[0] = T[0] - 0;
  res[n_Node-1] = T[n_DOFs-1] - 0;

  // The remaining nodes
  for (unsigned int i = 1; i < n_DOFs-1; i++)
    res[i] = (T[i] - T_old[i]) / _dt / alpha;
}

void
HeatConduction1D::RHS(double * rhs)
{
  // no spatial terms for the first and last node (applied Direchlet BC)

  for (unsigned int i = 1; i < n_DOFs-1; i++)
  {
    double qx = k * PI * PI * sin(PI * x[i]);
    // RHS = Laplacian(T) + q/k
    rhs[i] = (T[i-1] - 2*T[i] + T[i+1]) / dx2 + qx / k;
  }
}

void
HeatConduction1D::writeSolution()
{
  std::cout << std::setw(20) << "x"
            << std::setw(20) << "T"
            << std::setw(20) << "T_exact" << std::endl;
  for (unsigned int i = 0; i < T.size(); i++)
    std::cout << std::setw(20) << x[i]
              << std::setw(20) << T[i]
              << std::setw(20) << sin(PI * x[i]) << std::endl;
}

void
HeatConduction1D::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOFs, n_DOFs, 3, NULL, &P_Mat);

  PetscInt     row, col;
  PetscScalar  v = 1.0;

  row = 0; col = 0;
  MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);
  row = n_DOFs-1; col = n_DOFs-1;
  MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);

  for (row = 1; row < n_DOFs-1; row++)
    for (col = row-1; col <= row+1; col++)
      MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);

  /*
  MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD); */
}

void
HeatConduction1D::computeJacobianMatrix(Mat & P_Mat)
{
  PetscInt     row, col;
  PetscScalar  v = 1.0;

  row = 0; col = 0;
  MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);
  row = n_DOFs-1; col = n_DOFs-1;
  MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);

  PetscInt cols[3];
  PetscReal jac[3];
  for (row = 1; row < n_DOFs-1; row++)
  {
    cols[0] = row - 1; cols[1] = row; cols[2] = row + 1;
    jac[0] = -1.0 / dx2;
    jac[1] = 1.0 / _dt / alpha + 2.0 / dx2;
    jac[2] = -1.0 / dx2;

    MatSetValues(P_Mat, 1, &row, 3, cols, jac, INSERT_VALUES);
  }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
}

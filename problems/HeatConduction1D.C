#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "HeatConduction1D.h"

HeatConduction1D::HeatConduction1D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  // Using Finite Difference Method for heat conduction problem
  // 1/alpha * dT/dt = Laplacian(T) + q/k

  // 0 <= x <= 1
  // At steady state, let T(x) = sin(pi * x)
  // such that: q(x) = k * pi * pi * sin(pi * x)

  // This problem also has time-dependent analytical solutions:
  // T(x, t) = (1-exp(-alpha*pi*pi*t))*sin(pi*x)

  length = 1.0;
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");
  k = 1.0;
  alpha = 1.0;

  n_Node = n_Cell + 1;
  _n_DOFs = n_Node;

  dx = length / n_Cell;
  dx2 = dx * dx;

  x.resize(_n_DOFs);
  T.resize(_n_DOFs);
  T_old.resize(_n_DOFs);
  T_oldold.resize(_n_DOFs);

  // locations of nodes
  for(unsigned i = 0; i < x.size(); i++)
    x[i] = dx * i;
}

HeatConduction1D::~HeatConduction1D() {}

void
HeatConduction1D::SetupInitialCondition(double * u)
{
  for(unsigned i = 0; i < _n_DOFs; i++)
  { u[i] = 0.0; T[i] = 0.0; }

  T_oldold = T;     T_old = T;
}

void
HeatConduction1D::updateSolution(double * u)
{
  for(unsigned i = 0; i < T.size(); i++)  T[i] = u[i];
}

void
HeatConduction1D::transientResidual(double * res)
{
  // hard-coded BC; T_left = 0; T_right = 0
  res[0] = T[0] - 0;
  res[n_Node-1] = T[_n_DOFs-1] - 0;

  // The remaining nodes
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    // It is typical to use BDF1 for step 1 to startup BDF2
    // however, it is also associated with large error
    // see H. Nishikawa, "On large start-up error of BDF2", Journal of Computational Physics, Vol. 392, 2019, Pages 456-461
    for (unsigned i = 1; i < _n_DOFs-1; i++)
      res[i] = (1.5 * T[i] - 2.0 * T_old[i] + 0.5 * T_oldold[i]) / _dt / alpha;
  }
  else
  {
    for (unsigned i = 1; i < _n_DOFs-1; i++)
      res[i] = (T[i] - T_old[i]) / _dt / alpha;
  }
}

void
HeatConduction1D::RHS(double * rhs)
{
  // no spatial terms for the first and last node (applied Direchlet BC)

  for (unsigned i = 1; i < _n_DOFs-1; i++)
  {
    double qx = k * PI * PI * sin(PI * x[i]);
    // RHS = Laplacian(T) + q/k
    rhs[i] = (T[i-1] - 2*T[i] + T[i+1]) / dx2 + qx / k;
  }
}

void
HeatConduction1D::onTimestepEnd()
{
  // Save old solutions
  T_oldold = T_old; T_old = T;
}

void
HeatConduction1D::writeVTKOutput(FILE * file)
{
  fprintf(file, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_Node, n_Cell);
  fprintf(file, "      <Points>\n");
  fprintf(file, "        <DataArray type = \"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %f 0 0\n", i * dx);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Points>\n");

  fprintf(file, "      <Cells>\n");
  fprintf(file, "        <DataArray type = \"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d %d\n", i, i+1);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type = \"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d\n", 2*(i+1));
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type = \"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d\n", 3);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  fprintf(file, "      <CellData>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "mu");
  double t = _problemSystem->getCurrentTime();
  for (unsigned i = 0; i < n_Cell; i++) // A fake cell data
    fprintf(file, "          %20.6f\n", T[i]*exp(x[i]));
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </CellData>\n");

  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "temperature");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", T[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "T_exact");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", (1.0 - exp(-alpha * PI * PI * t)) * sin(PI * x[i]));
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
HeatConduction1D::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  mnzp->addEntry(0, 0); mnzp->addEntry(0, 1);// (0, 1) not needed because of Direchlet BC, otherwise needed
  mnzp->addEntry(_n_DOFs-1, _n_DOFs-1); mnzp->addEntry(_n_DOFs-1, _n_DOFs-2);
  for (unsigned row = 1; row < _n_DOFs-1; row++)
    for (unsigned col = row-1; col <= row+1; col++)
      mnzp->addEntry(row + _DOF_offset, col + _DOF_offset);
}

void
HeatConduction1D::computeJacobianMatrix(Mat & P_Mat)
{
  PetscInt     row, col;
  PetscScalar  one = 1.0;

  row = 0; col = 0;
  MatSetValues(P_Mat, 1, &row, 1, &col, &one, INSERT_VALUES);
  row = _n_DOFs-1; col = _n_DOFs-1;
  MatSetValues(P_Mat, 1, &row, 1, &col, &one, INSERT_VALUES);

  PetscInt cols[3];
  PetscReal jac[3];
  for (row = 1; row < _n_DOFs-1; row++)
  {
    cols[0] = row - 1; cols[1] = row; cols[2] = row + 1;
    jac[0] = (_time_scheme == CN) ? -0.5 / dx2 : -1.0 / dx2;
    jac[1] = (_time_scheme == CN) ? 1.0 / _dt / alpha + 1.0 / dx2 : 1.0 / _dt / alpha + 2.0 / dx2;
    jac[2] = (_time_scheme == CN) ? -0.5 / dx2 : -1.0 / dx2;

    MatSetValues(P_Mat, 1, &row, 3, cols, jac, INSERT_VALUES);
  }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
}

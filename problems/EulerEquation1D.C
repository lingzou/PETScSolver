#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "EulerEquation1D.h"

EulerEquation1D::EulerEquation1D() :
  PETScProblem()
{
  _order = PetscOptionsGetRequiredInt("-order");
  _gamma = 1.4;

  length = 1.0;
  n_Cell = PetscOptionsGetRequiredInt("-n_cells");
  n_Node = n_Cell + 1;
  n_DOFs = n_Cell * 3;

  dx = length / n_Cell;

  int p_case = PetscOptionsGetRequiredInt("-p_case");
  if (p_case == 1)  // Sod problem
  {
    RHO_L = 1.0;      M_L = 0.0;    E_L = 2.5;
    RHO_R = 0.125;    M_R = 0.0;    E_R = 0.25;
  }
  else if (p_case == 2)  // Lax problem
  {
    RHO_L = 0.445;    M_L = 0.311;    E_L = 8.928;
    RHO_R = 0.500;    M_R = 0.000;    E_R = 1.4275;
  }
  else
    sysError("case = 1 for Sod problem; case = 2 for Lax problem");

  P_L = p_IG(RHO_L, M_L, E_L);
  P_R = p_IG(RHO_R, M_R, E_R);

  rho.resize(n_Cell);       m.resize(n_Cell);       E.resize(n_Cell);     p.resize(n_Cell);
  rho_old.resize(n_Cell);   m_old.resize(n_Cell);   E_old.resize(n_Cell);
  rho_oo.resize(n_Cell);    m_oo.resize(n_Cell);    E_oo.resize(n_Cell);

  flux_rho.resize(n_Cell + 1); flux_m.resize(n_Cell + 1); flux_E.resize(n_Cell + 1);

  if (_order > 1)
  {
    rho_w.resize(n_Cell);       m_w.resize(n_Cell);       E_w.resize(n_Cell);
    rho_e.resize(n_Cell);       m_e.resize(n_Cell);       E_e.resize(n_Cell);
  }
}

EulerEquation1D::~EulerEquation1D()
{
}

void
EulerEquation1D::SetupInitialCondition(double * u)
{
  unsigned int index = 0;

  for(unsigned int i = 0; i < n_Cell; i++)
  {
    rho[i] = (i < n_Cell/2) ? RHO_L : RHO_R;
    m[i]   = (i < n_Cell/2) ? M_L   : M_R;
    E[i]   = (i < n_Cell/2) ? E_L   : E_R;
    p[i]   = (i < n_Cell/2) ? P_L   : P_R;

    u[index++] = rho[i];
    u[index++] = m[i];
    u[index++] = E[i];
  }
  rho_old = rho;  m_old = m; E_old = E;
  rho_oo = rho;   m_oo = m;  E_oo = E;
}

void
EulerEquation1D::updateSolution(double * u, TimeStepIndex index)
{
  unsigned int idx = 0;
  switch (index)
  {
    case NEW:
      for(unsigned int i = 0; i < n_Cell; i++)
      {
        rho[i] = u[idx++];  m[i] = u[idx++];  E[i] = u[idx++];
        p[i] = p_IG(rho[i], m[i], E[i]);
      }
      break;

    case OLD:
      for(unsigned int i = 0; i < n_Cell; i++)
      {
        rho_old[i] = u[idx++];  m_old[i] = u[idx++];  E_old[i] = u[idx++];
      }
      break;

    case OLDOLD:
      for(unsigned int i = 0; i < n_Cell; i++)
      {
        rho_oo[i] = u[idx++];  m_oo[i] = u[idx++];  E_oo[i] = u[idx++];
      }
      break;

    default:
      sysError("Unknown TimeStepIndex");
  }
}

void
EulerEquation1D::transientResidual(double * res)
{
  unsigned int idx = 0;
  if ((_time_scheme == BDF2) && (_step > 1))
  {
    // It is typical to use BDF1 for step 1 to startup BDF2
    // however, it is also associated with large error
    // see H. Nishikawa, "On large start-up error of BDF2", Journal of Computational Physics, Vol. 392, 2019, Pages 456-461
    for(unsigned int i = 0; i < n_Cell; i++)
    {
      res[idx++] = (1.5 * rho[i] - 2.0 * rho_old[i] + 0.5 * rho_oo[i]) / _dt;
      res[idx++] = (1.5 * m[i]   - 2.0 * m_old[i]   + 0.5 * m_oo[i])   / _dt;
      res[idx++] = (1.5 * E[i]   - 2.0 * E_old[i]   + 0.5 * E_oo[i])   / _dt;
    }
  }
  else
  {
    for(unsigned int i = 0; i < n_Cell; i++)
    {
      res[idx++] = (rho[i] - rho_old[i]) / _dt;
      res[idx++] = (m[i]   - m_old[i])   / _dt;
      res[idx++] = (E[i]   - E_old[i])   / _dt;
    }
  }
}

void
EulerEquation1D::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:   updateFluxes();           break;
    case 2:   updateFluxes2ndOrder();   break;
    defaut:   sysError("Spatial order not implemented.");
  }

  unsigned int idx = 0;
  for(unsigned int i = 0; i < n_Cell; i++)
  {
    rhs[idx++] = -(flux_rho[i+1] - flux_rho[i]) / dx;
    rhs[idx++] = -(flux_m[i+1] - flux_m[i]) / dx;
    rhs[idx++] = -(flux_E[i+1] - flux_E[i]) / dx;
  }
}

void
EulerEquation1D::updateFluxes()
{
  for(unsigned int i = 0; i <= n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double rho_left = (i == 0) ? RHO_L : rho[i-1];
    double m_left   = (i == 0) ? M_L   : m[i-1];
    double E_left   = (i == 0) ? E_L   : E[i-1];
    double u_left   = m_left / rho_left;
    double p_left   = p_IG(rho_left, m_left, E_left);
    double c_left   = std::sqrt(_gamma * p_left / rho_left);

    double rho_right= (i == n_Cell) ? RHO_R : rho[i];
    double m_right  = (i == n_Cell) ? M_R   : m[i];
    double E_right  = (i == n_Cell) ? E_R   : E[i];
    double u_right  = m_right / rho_right;
    double p_right  = p_IG(rho_right, m_right, E_right);
    double c_right  = std::sqrt(_gamma * p_right / rho_right);

    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m_left + m_right) + 0.5 * eigen_max * (rho_left - rho_right);
    flux_m[i] = 0.5 * (m_left*m_left/rho_left + p_left + m_right*m_right/rho_right + p_right) + 0.5 * eigen_max * (m_left - m_right);
    flux_E[i] = 0.5 * ((E_left + p_left)*u_left + (E_right + p_right)*u_right) + 0.5 * eigen_max * (E_left - E_right);
  }
}

void
EulerEquation1D::updateFluxes2ndOrder()
{
  // Reference:
  //   A. Kurganov and E. Tadmor, New High-Resolution Central Schemes for Nonlinear Conservation Laws
  //   and Convection–Diffusion Equations, Journal of Computational Physics, 160, 241–282 (2000)
  // Note:
  //   The reconstruction is directly performed on the conservative variables.
  //   Some papers (I cannot find the references now FIXME) claim that it is preferred to use primitive
  //   variables for reconstruction. This has not been tested in this code.
  linearReconstruction(RHO_L, RHO_R, rho, rho_w, rho_e);
  linearReconstruction(M_L, M_R, m, m_w, m_e);
  linearReconstruction(E_L, E_R, E, E_w, E_e);

  for(unsigned int i = 0; i <= n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double rho_left = (i == 0) ? RHO_L : rho_e[i-1];
    double m_left   = (i == 0) ? M_L   : m_e[i-1];
    double E_left   = (i == 0) ? E_L   : E_e[i-1];
    double u_left   = m_left / rho_left;
    double p_left   = p_IG(rho_left, m_left, E_left);
    double c_left   = std::sqrt(_gamma * p_left / rho_left);

    double rho_right= (i == n_Cell) ? RHO_R : rho_w[i];
    double m_right  = (i == n_Cell) ? M_R   : m_w[i];
    double E_right  = (i == n_Cell) ? E_R   : E_w[i];
    double u_right  = m_right / rho_right;
    double p_right  = p_IG(rho_right, m_right, E_right);
    double c_right  = std::sqrt(_gamma * p_right / rho_right);

    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m_left + m_right) + 0.5 * eigen_max * (rho_left - rho_right);
    flux_m[i] = 0.5 * (m_left*m_left/rho_left + p_left + m_right*m_right/rho_right + p_right) + 0.5 * eigen_max * (m_left - m_right);
    flux_E[i] = 0.5 * ((E_left + p_left)*u_left + (E_right + p_right)*u_right) + 0.5 * eigen_max * (E_left - E_right);
  }
}

void
EulerEquation1D::linearReconstruction(double l_ghost, double r_ghost,
  std::vector<double> &u, std::vector<double> &u_w, std::vector<double> &u_e)
{
  for(int i = 0; i < u.size(); i++)
  {
    double u_P = u[i];
    double u_W = (i == 0) ? l_ghost : u[i - 1];
    double u_E = (i == u.size() - 1) ? r_ghost : u[i + 1];

    // I think this is a more decent implementation of TVD limiter
    // Reference:
    // [1] Berger, M., Aftosmis, M.J., Murman, S.M., 2005. Analysis of slope limiters on irreg- ular grids.
    //     In: AIAA Paper 2005-0490, 43rd AIAA Aerospace Sciences Meeting, Jan. 10e13, Reno, NV, 2005.
    double wf = 0.0;
    if ((u_E - u_P)*(u_P - u_W) > 0.0)    // u_P is in between u_W and u_E, i.e., not a local extremum
      if (std::fabs(u_E - u_W) > 1.e-10)  // The difference is large enough, so reconstruction is meaningful
      {
        double f = (u_P - u_W) / (u_E - u_W);
        wf = 2.0 * f * (1.0 - f) / ((1.0 - f) * (1.0 - f) + f * f);   // Eqn. (10) of Ref. [1]
      }
    // Linear reconstructed values
    u_w[i] = u_P - 0.25 * wf * (u_E - u_W);   // Eqn. (5) of Ref. [1]
    u_e[i] = u_P + 0.25 * wf * (u_E - u_W);
  }
}

void
EulerEquation1D::onTimestepEnd()
{
  // March time forward
  _t += _dt;

  // write solution
  writeSolution(_step);
  _step ++;

  // Copy oldold solution into old; new solution into old solution
  if (_time_scheme == BDF2)
  { rho_oo = rho_old; m_oo = m_old; E_oo = E_old; }

  rho_old = rho; m_old = m; E_old = E;
}

void
EulerEquation1D::writeSolution(unsigned int step)
{
  FILE * ptr_File;
  std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".vtk";
  ptr_File = fopen(file_name.c_str(), "w");

  fprintf(ptr_File, "# vtk DataFile Version 4.0\n");
  fprintf(ptr_File, "my data\n");
  fprintf(ptr_File, "ASCII\n");
  fprintf(ptr_File, "DATASET STRUCTURED_GRID\n");
  fprintf(ptr_File, "DIMENSIONS %u 1 1\n", n_Node);
  fprintf(ptr_File, "POINTS %u Float32\n", n_Node);
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f 0 0\n", i * dx);

  // point data
  fprintf(ptr_File, "POINT_DATA %u\n", n_Node);
  // node id
  fprintf(ptr_File, "SCALARS node_id Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE node_id\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%d\n", i);

  // cell data
  fprintf(ptr_File, "CELL_DATA %u\n", n_Cell);
  // rho
  fprintf(ptr_File, "SCALARS rho Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE rho\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", rho[i]);

  // m
  fprintf(ptr_File, "SCALARS m Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE m\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", m[i]);

  // E
  fprintf(ptr_File, "SCALARS E Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE E\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", E[i]);

  // p
  fprintf(ptr_File, "SCALARS p Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", p[i]);

  // velocity u = m/rho
  fprintf(ptr_File, "SCALARS u Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE u\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", m[i]/rho[i]);

  fclose(ptr_File);
}

void
EulerEquation1D::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOFs, n_DOFs, 15, NULL, &P_Mat);

  int n_Var = 3;
  PetscReal v = 1.0;
  PetscInt row, col;
  for (int i = 0; i < n_Cell; i++)                 // loop on cells
    for (int i_var = 0; i_var < n_Var; i_var++)    // loop on the 3 variables on each cells
    {
      row = i * n_Var + i_var;                     // This loops on rows
      for (int j = i-2; j <= i+2; j++)             // loop on its -2 to 2 neighbors; maybe out of bound
        for (int j_var = 0; j_var < n_Var; j_var++)
        {
          // This loops on variables on neighboring cells
          col = j * n_Var + j_var;  // might be smaller than zero
          if ((col >= 0) && (col < n_DOFs))
            MatSetValues(P_Mat, 1, &row, 1, &col, &v, INSERT_VALUES);
        }
    }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
  /*
  MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD);
  */
}

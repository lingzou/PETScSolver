#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "EulerEquation1D.h"

EulerEquation1D::EulerEquation1D() :
  PETScProblem()
{
  _order = 2;
  _gamma = 1.4;

  length = 1.0;
  n_Cell = 100;
  n_Node = n_Cell + 1;
  n_DOFs = n_Cell * 3;

  dx = length / n_Cell;

  rho.resize(n_DOFs);       m.resize(n_DOFs);       E.resize(n_DOFs);     p.resize(n_DOFs);
  rho_old.resize(n_DOFs);   m_old.resize(n_DOFs);   E_old.resize(n_DOFs);
  rho_oo.resize(n_DOFs);    m_oo.resize(n_DOFs);    E_oo.resize(n_DOFs);

  flux_rho.resize(n_Cell + 1); flux_m.resize(n_Cell + 1); flux_E.resize(n_Cell + 1);

  if (_order > 1)
  {
    rho_w.resize(n_DOFs);       m_w.resize(n_DOFs);       E_w.resize(n_DOFs);
    rho_e.resize(n_DOFs);       m_e.resize(n_DOFs);       E_e.resize(n_DOFs);
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
    rho[i] = (i < n_Cell/2) ? 1.0 : 0.125;   // rho
    m[i] = 0.0;  // m; u=0
    p[i] = (i < n_Cell/2) ? 1.0 : 0.1;
    E[i] = p[i] / (_gamma - 1.0);   // E = rho*e + 0.5*u*u = rho*e

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
        p[i] = (_gamma - 1.0) * (E[i] - 0.5 * m[i] * m[i] / rho[i]);
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
      std::cerr << "ERROR.\n";
      exit(1);
      break;
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
    defaut:   std::cerr << "Spatial order not implemented." << std::endl; exit(1);
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
  // Left wall
  flux_rho[0] = 0.0;    flux_m[0] = p[0];   flux_E[0] = 0.0;
  // Right wall
  flux_rho[n_Cell] = 0.0;    flux_m[n_Cell] = p[n_Cell-1];   flux_E[n_Cell] = 0.0;

  // Everything in between
  for(unsigned int i = 1; i < n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double u_left = m[i-1] / rho[i-1];
    double c_left = std::sqrt(_gamma * p[i-1] / rho[i-1]);
    double u_right = m[i] / rho[i];
    double c_right = std::sqrt(_gamma * p[i] / rho[i]);
    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m[i-1] + m[i]) + 0.5 * eigen_max * (rho[i-1] - rho[i]);
    flux_m[i] = 0.5 * (m[i-1]*m[i-1]/rho[i-1] + p[i-1] + m[i]*m[i]/rho[i] + p[i]) + 0.5 * eigen_max * (m[i-1] - m[i]);
    flux_E[i] = 0.5 * ((E[i-1] + p[i-1])*u_left + (E[i] + p[i])*u_right) + 0.5 * eigen_max * (E[i-1] - E[i]);
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
  linearReconstruction(rho[0], rho[n_Cell - 1], rho, rho_w, rho_e);
  linearReconstruction(-m[0], -m[n_Cell - 1], m, m_w, m_e);
  linearReconstruction(E[0], E[n_Cell - 1], E, E_w, E_e);

  // Left wall
  flux_rho[0] = 0.0;    flux_m[0] = p[0];   flux_E[0] = 0.0;
  // Right wall
  flux_rho[n_Cell] = 0.0;    flux_m[n_Cell] = p[n_Cell-1];   flux_E[n_Cell] = 0.0;

  // Everything in between
  for(unsigned int i = 1; i < n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double u_left = m_e[i-1] / rho_e[i-1];
    double p_left = (_gamma - 1.0) * (E_e[i-1] - 0.5 * m_e[i-1] * m_e[i-1] / rho_e[i-1]);
    double c_left = std::sqrt(_gamma * p_left / rho_e[i-1]);
    double u_right = m_w[i] / rho_w[i];
    double p_right = (_gamma - 1.0) * (E_w[i] - 0.5 * m_w[i] * m_w[i] / rho_w[i]);
    double c_right = std::sqrt(_gamma * p_right / rho_w[i]);
    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m_e[i-1] + m_w[i]) + 0.5 * eigen_max * (rho_e[i-1] - rho_w[i]);
    flux_m[i] = 0.5 * (m_e[i-1]*m_e[i-1]/rho_e[i-1] + p_left + m_w[i]*m_w[i]/rho_w[i] + p_right) + 0.5 * eigen_max * (m_e[i-1] - m_w[i]);
    flux_E[i] = 0.5 * ((E_e[i-1] + p_left)*u_left + (E_w[i] + p_right)*u_right) + 0.5 * eigen_max * (E_e[i-1] - E_w[i]);
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

    double r = 0.0;
    if ((u_E - u_P)*(u_P - u_W) < 0.0) // local extremum
      r = 0.0;
    else // u_P is in between u_E and u_W
    {
      if (std::fabs(u_E - u_W) < 1.e-10)  // do not bother to perform reconstruction
        r = 0.0;
      else
        r = (u_P - u_W) / (u_E - u_P);
    }
    // The codes above could be tranformed into a more compact form,
    //   but really not necessary to do so.
    // r = (((u_E - u_P)*(u_P - u_W) > 0.0) && (std::fabs(u_E - u_W) > 1.e-10)) ?
    //        (u_P - u_W) / (u_E - u_P) : 0.0;

    // Van Albada 1
    double phir = (r > 0.0) ? (r + r*r) / (1.0 + r*r) : 0.0;

    // Linear reconstructed values
    u_w[i] = u_P - 0.5 * phir * (u_E - u_P);
    u_e[i] = u_P + 0.5 * phir * (u_E - u_P);
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
  std::string file_name = "output/EulerEquation1D_step_" + std::to_string(step) + ".vtk";
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
  // Temperature point data
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

  fclose(ptr_File);
}

void
EulerEquation1D::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOFs, n_DOFs, 9, NULL, &P_Mat);

  PetscInt row;
  for(unsigned int var = 0; var < 3; var++)
  {
    // First cell
    row = var;
    PetscInt cols_first[6] = {0, 1, 2, 3, 4, 5};
    PetscScalar ones[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    MatSetValues(P_Mat, 1, &row, 6, cols_first, ones, INSERT_VALUES);
    // Last cell
    row = 3 * (n_Cell - 1) + var;
    PetscInt start = 3 * (n_Cell - 2);
    PetscInt cols_last[6] = {start, start + 1, start + 2,
                             start + 3, start + 4, start + 5};
    MatSetValues(P_Mat, 1, &row, 6, cols_last, ones, INSERT_VALUES);
  }

  // Cells in between
  for(unsigned int i = 1; i < n_Cell - 1; i++)
  {
    for(unsigned int var = 0; var < 3; var++)
    {
      row = i * 3 + var;
      PetscInt s = 3 * (i - 1);
      PetscInt cols[9] = {s, s+1, s+2, s+3, s+4, s+5, s+6, s+7, s+8};
      PetscScalar ones[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      MatSetValues(P_Mat, 1, &row, 9, cols, ones, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
  /*
  MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD); */
}

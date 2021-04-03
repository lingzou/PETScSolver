#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "FiveEqnTwoP_StagGrid.h"

// Reference:
//   [1] L. Zou, H. Zhao, and H. Zhang, A revisit to the Hicks' hyperbolic two-pressure
//       two-phase flow model, The 17th International Topical Meeting on Nuclear Reactor
//       Thermal Hydraulics (NURETH-17), Xi’an, China, September 3-8, 2017

/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0    2    5    7                                            5n
   1    3    6    8                                            5n+1
        4         9
*/

FiveEqnTwoP_StagGrid::FiveEqnTwoP_StagGrid() :
  PETScProblem()
{
  ALPHA_INIT =  0.2;
  V_L_INIT   = 10.0;
  V_G_INIT   =  0.0;
  P_INIT     =  1.0e5;
  P_OUTLET   =  1.0e5;
  ALPHA_INLET  = 0.2;
  ALPHA_OUTLET = 0.2;

  g = -9.81;

  C_L = std::sqrt(1.0e7);
  C_G = std::sqrt(1.0e6);

  _order = PetscOptionsGetRequiredInt("-order");
  H_inv = PetscOptionsGetRequiredReal("-H_inv");

  length = 12.0;
  n_Cell = PetscOptionsGetRequiredInt("-n_cells");
  n_Node = n_Cell + 1;
  n_DOFs = n_Cell * 5 + 2;

  dx = length / n_Cell;

  // Primary variables
  alpha.resize(n_Cell);       p_l.resize(n_Cell);       p_g.resize(n_Cell);       v_l.resize(n_Node);       v_g.resize(n_Node);
  alpha_old.resize(n_Cell);   p_l_old.resize(n_Cell);   p_g_old.resize(n_Cell);   v_l_old.resize(n_Node);   v_g_old.resize(n_Node);
  alpha_oo.resize(n_Cell);    p_l_oo.resize(n_Cell);    p_g_oo.resize(n_Cell);    v_l_oo.resize(n_Node);    v_g_oo.resize(n_Node);

  // Helper variables
  rho_l.resize(n_Cell);       rho_g.resize(n_Cell);
  rho_l_old.resize(n_Cell);   rho_g_old.resize(n_Cell);
  rho_l_oo.resize(n_Cell);    rho_g_oo.resize(n_Cell);

  alpha_edge.resize(n_Node);  rho_l_edge.resize(n_Node);  rho_g_edge.resize(n_Node);
  v_l_cell.resize(n_Cell);    v_g_cell.resize(n_Cell);

  v_hat.resize(n_Node);
  p_hat.resize(n_Cell);
  mu.resize(n_Cell);
  p_hat_edge.resize(n_Node);

  // Fluxes
  alpha_flux.resize(n_Node);  rho_l_flux.resize(n_Node);  rho_g_flux.resize(n_Node);
}

FiveEqnTwoP_StagGrid::~FiveEqnTwoP_StagGrid()
{
}

void
FiveEqnTwoP_StagGrid::SetupInitialCondition(double * u)
{
  unsigned int index = 0;
  for(int i = 0; i < n_Cell + 1; i++)
  {
    v_l[i]     = V_L_INIT;
    v_g[i]     = V_G_INIT;
    u[index++] = v_l[i];
    u[index++] = v_g[i];

    if (i < n_Cell)
    {
      alpha[i]   = ALPHA_INIT;
      p_l[i]     = P_INIT;
      p_g[i]     = P_INIT;
      rho_l[i] = rho_l_func(p_l[i]);
      rho_g[i] = rho_g_func(p_g[i]);

      u[index++] = alpha[i];
      u[index++] = p_l[i];
      u[index++] = p_g[i];
    }
  }

  alpha_old = alpha;  p_l_old = p_l;  p_g_old = p_g;  v_l_old = v_l;  v_g_old = v_g;
  alpha_oo  = alpha;  p_l_oo  = p_l;  p_g_oo  = p_g;  v_l_oo  = v_l;  v_g_oo  = v_g;

  rho_l_old = rho_l;  rho_g_old = rho_g;
  rho_l_oo  = rho_l;  rho_g_oo  = rho_g;
}

void
FiveEqnTwoP_StagGrid::updateSolution(double * u, TimeStepIndex index)
{
  unsigned int idx = 0;
  switch (index)
  {
    case NEW:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v_l[i] = u[idx++];  v_g[i] = u[idx++];
        if (i < n_Cell)
        {
          alpha[i] = u[idx++];  p_l[i] = u[idx++];  p_g[i] = u[idx++];
          rho_l[i] = rho_l_func(p_l[i]);
          rho_g[i] = rho_g_func(p_g[i]);
        }
      }
      //update cell values
      for(int i = 0; i < n_Cell; i++)
      {
        mu[i]    = H_inv / (rho_l[i] * C_L + rho_g[i] * C_G);
        p_hat[i] = (p_l[i] * rho_g[i] * C_G + p_g[i] * rho_l[i] * C_L) / (rho_l[i] * C_L + rho_g[i] * C_G);
        v_l_cell[i] = 0.5 * (v_l[i] + v_l[i+1]);
        v_g_cell[i] = 0.5 * (v_g[i] + v_g[i+1]);
      }
      //update edge values
      rho_l_edge[0] = rho_l[0];
      rho_g_edge[0] = rho_g[0];
      alpha_edge[0] = ALPHA_INIT;
      p_hat_edge[0] = p_hat[0];
      for(unsigned int i = 1; i < n_Cell; i++)
      {
        rho_l_edge[i] = 0.5 * (rho_l[i-1] + rho_l[i]);
        rho_g_edge[i] = 0.5 * (rho_g[i-1] + rho_g[i]);
        alpha_edge[i] = 0.5 * (alpha[i-1] + alpha[i]);
        p_hat_edge[i] = 0.5 * (p_hat[i-1] + p_hat[i]);
      }
      rho_l_edge[n_Cell] = rho_l[n_Cell - 1];
      rho_g_edge[n_Cell] = rho_g[n_Cell - 1];
      alpha_edge[n_Cell] = alpha[n_Cell - 1]; // FIXME FIXME
      p_hat_edge[n_Cell] = p_hat[n_Cell - 1];

      for(unsigned int i = 0; i < n_Cell + 1; i++)
        v_hat[i] = alpha_edge[i] * v_g[i] + (1.0 - alpha_edge[i]) * v_l[i];
      break;

    case OLD:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v_l_old[i] = u[idx++];  v_g_old[i] = u[idx++];
        if (i < n_Cell)
        {
          alpha_old[i] = u[idx++];  p_l_old[i] = u[idx++];  p_g_old[i] = u[idx++];
          rho_l_old[i] = rho_l_func(p_l_old[i]);
          rho_g_old[i] = rho_g_func(p_g_old[i]);
        }
      }
      break;

    case OLDOLD:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v_l_oo[i] = u[idx++];  v_g_oo[i] = u[idx++];
        if (i < n_Cell)
        {
          alpha_oo[i] = u[idx++];  p_l_oo[i] = u[idx++];  p_g_oo[i] = u[idx++];
          rho_l_oo[i] = rho_l_func(p_l_oo[i]);
          rho_g_oo[i] = rho_g_func(p_g_oo[i]);
        }
      }
      break;

    default:
      std::cerr << "ERROR.\n";
      exit(1);
      break;
  }
}

void
FiveEqnTwoP_StagGrid::transientResidual(double * res)
{
  unsigned int idx = 0;
  if ((_time_scheme == BDF2) && (_step > 1))
  {
    for(unsigned int i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (1.0 - alpha_edge[i]) * rho_l_edge[i] * (1.5 * v_l[i] - 2.0 * v_l_old[i] + 0.5 * v_l_oo[i]) / _dt;
      res[idx++] = alpha_edge[i] * rho_g_edge[i] * (1.5 * v_g[i] - 2.0 * v_g_old[i] + 0.5 * v_g_oo[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (1.5 * alpha[i] - 2.0 * alpha_old[i] + 0.5 * alpha_oo[i]) / _dt;
        res[idx++] = (1.5 * (1.0 - alpha[i]) * rho_l[i] - 2.0 * (1.0 - alpha_old[i]) * rho_l_old[i] + 0.5 * (1.0 - alpha_oo[i]) * rho_l_oo[i]) / _dt;
        res[idx++] = (1.5 * alpha[i] * rho_g[i] - 2.0 * alpha_old[i] * rho_g_old[i] + 0.5 * alpha_oo[i] * rho_g_oo[i]) / _dt;
      }
    }
  }
  else
  {
    for(unsigned int i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (1.0 - alpha_edge[i]) * rho_l_edge[i] * (v_l[i] - v_l_old[i]) / _dt;
      res[idx++] = alpha_edge[i] * rho_g_edge[i]         * (v_g[i] - v_g_old[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (alpha[i]                    - alpha_old[i])                        / _dt;
        res[idx++] = ((1.0 - alpha[i]) * rho_l[i] - (1.0 - alpha_old[i]) * rho_l_old[i]) / _dt;
        res[idx++] = (alpha[i] * rho_g[i]         - alpha_old[i] * rho_g_old[i])         / _dt;
      }
    }
  }

  // Zero-out transient residual, to apply Dirichlet BC at the inlet later
  res[0] = 0.0;
  res[1] = 0.0;
}

void
FiveEqnTwoP_StagGrid::RHS(double * rhs)
{
  // Boundary values and ghost values
  double p_l_inlet_ghost = p_l[0];
  double p_g_inlet_ghost = p_g[0];

  double rho_l_inlet_bc = rho_l_func(p_l_inlet_ghost);
  double rho_g_inlet_bc = rho_g_func(p_g_inlet_ghost);

  double rho_l_outlet_bc = rho_l_func(P_OUTLET);
  double rho_g_outlet_bc = rho_g_func(P_OUTLET);

  double alpha_outlet_ghost = (v_l[n_Cell] > 0) ? alpha[n_Cell - 1] : ALPHA_OUTLET;

  double v_l_outlet_ghost = v_l[n_Cell];
  double v_g_outlet_ghost = v_g[n_Cell];

  // Upwind donor cell method for void fraction and mass balance equations
  alpha_flux[0] = (v_hat[0] > 0) ? v_hat[0] * ALPHA_INLET : v_hat[0] * alpha[0];
  rho_l_flux[0] = (V_L_INIT > 0) ? V_L_INIT * (1.0 - ALPHA_INLET) * rho_l_inlet_bc : V_L_INIT * (1.0 - alpha[0]) * rho_l[0];
  rho_g_flux[0] = (V_G_INIT > 0) ? V_G_INIT * ALPHA_INLET * rho_g_inlet_bc : V_G_INIT * alpha[0] * rho_g[0];

  alpha_flux[n_Cell] = (v_hat[n_Cell] > 0) ? v_hat[n_Cell] * alpha[n_Cell-1] : v_hat[n_Cell] * ALPHA_OUTLET;
  rho_l_flux[n_Cell] = (v_l[n_Cell] > 0) ? v_l[n_Cell] * (1. - alpha[n_Cell-1]) * rho_l[n_Cell-1] : v_l[n_Cell] * (1.0 - ALPHA_OUTLET) * rho_l_outlet_bc;
  rho_g_flux[n_Cell] = (v_g[n_Cell] > 0) ? v_g[n_Cell] * alpha[n_Cell-1] * rho_g[n_Cell-1] : v_g[n_Cell] * alpha_outlet_ghost * rho_g_outlet_bc;

  for(unsigned int i = 1; i < n_Cell; i++)
  {
    alpha_flux[i] = (v_hat[i] > 0) ? v_hat[i] * alpha[i-1] : v_hat[i] * alpha[i];
    rho_l_flux[i] = (v_l[i] > 0) ? v_l[i] * (1. - alpha[i-1]) * rho_l[i-1] : v_l[i] * (1. - alpha[i]) * rho_l[i];
    rho_g_flux[i] = (v_g[i] > 0) ? v_g[i] * alpha[i-1] * rho_g[i-1] : v_g[i] * alpha[i] * rho_g[i];
  }

  // Momentum equations RHS; Eqn. (22) of Ref. [1]
  rhs[0] = v_l[0] - V_L_INIT;
  rhs[1] = v_g[0] - V_G_INIT;
  for(int i = 1; i < n_Cell + 1; i++) // loop on the remaining edges
  {
    // east and west velocities
    double v_l_east = (i == n_Cell) ? v_l[n_Cell] : v_l[i+1];
    double v_l_west = (i == 1) ? V_L_INIT : v_l[i-1];
    double v_g_east = (i == n_Cell) ? v_g[n_Cell] : v_g[i+1];
    double v_g_west = (i == 1) ? V_G_INIT : v_g[i-1];

    // dv_dx term
    double dv_l_dx = (v_l[i] > 0) ? (v_l[i] - v_l_west)/dx : (v_l_east - v_l[i])/dx;
    double dv_g_dx = (v_g[i] > 0) ? (v_g[i] - v_g_west)/dx : (v_g_east - v_g[i])/dx;

    // d(alpha*p)_dx term
    double d_alphap_l_dx = (i == n_Cell) ? (P_OUTLET * (1. - alpha_outlet_ghost) - p_l[n_Cell-1]*(1. - alpha[n_Cell-1]))/dx*2. : (p_l[i]*(1. - alpha[i]) - p_l[i-1]*(1. - alpha[i-1]))/dx;
    double d_alphap_g_dx = (i == n_Cell) ? (P_OUTLET * alpha_outlet_ghost - p_g[n_Cell-1]*alpha[n_Cell-1])/dx*2. : (p_g[i]*alpha[i] - p_g[i-1]*alpha[i-1])/dx;
    double alpha_grad = (i == n_Cell) ? (alpha_outlet_ghost - alpha[n_Cell-1])/dx*2. : (alpha[i] - alpha[i-1])/dx;

    // assemble RHS for liquid and gas phase momentum equations
    double rhs_v_l = -(1. - alpha_edge[i]) * rho_l_edge[i] * v_l[i] * dv_l_dx;
    double rhs_v_g = -alpha_edge[i] * rho_g_edge[i] * v_g[i] * dv_g_dx;

    rhs_v_l -= d_alphap_l_dx;
    rhs_v_l -= p_hat_edge[i] * alpha_grad;
    rhs_v_g -= d_alphap_g_dx;
    rhs_v_g += p_hat_edge[i] * alpha_grad;

    // gravity terms, alpha * rho * g
    if (i < n_Cell)
    {
      rhs_v_l -= 0.5 * ((1. - alpha[i-1]) * rho_l[i-1] + (1. - alpha[i]) * rho_l[i]) * g;
      rhs_v_g -= 0.5 * (alpha[i-1] * rho_g[i-1] + alpha[i] * rho_g[i]) * g;
    }
    else
    {
      rhs_v_l -= (1. - alpha[n_Cell-1]) * rho_l[n_Cell-1] * g;
      rhs_v_g -= alpha[n_Cell-1] * rho_g[n_Cell-1] * g;
    }

    rhs[5*i] = rhs_v_l;
    rhs[5*i+1] = rhs_v_g;
  }

  // RHS for void fraction and mass conservation equations
  for(int i = 0; i < n_Cell; i++) //loop on cells
  {
    // alpha eqn, Eqn. (23) of Ref. [1]
    double v_hat_grad = (v_hat[i+1] - v_hat[i]) / dx;
    rhs[5*i+2] = -(alpha_flux[i+1] - alpha_flux[i]) / dx + alpha[i] * v_hat_grad + mu[i] * (p_g[i] - p_l[i]);

    // mass conservation equations, Eqn. (2) and (4) of Ref. [1]
    rhs[5*i+3] = -(rho_l_flux[i+1] - rho_l_flux[i]) / dx;
    rhs[5*i+4] = -(rho_g_flux[i+1] - rho_g_flux[i]) / dx;
  }
}

void
FiveEqnTwoP_StagGrid::updateFluxes()
{
}

void
FiveEqnTwoP_StagGrid::updateFluxes2ndOrder()
{
}

void
FiveEqnTwoP_StagGrid::linearReconstruction(double l_ghost, double r_ghost,
  std::vector<double> &u, std::vector<double> &u_w, std::vector<double> &u_e)
{
}

void
FiveEqnTwoP_StagGrid::onTimestepEnd()
{
  // March time forward
  _t += _dt;

  // write solution
  int N_Steps = PetscOptionsGetRequiredInt("-n_steps");
  if ((_step % _output_interval == 0) || (_step == N_Steps))
    writeSolution(_step);
  _step ++;
}

void
FiveEqnTwoP_StagGrid::writeSolution(unsigned int step)
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

  // v_l
  fprintf(ptr_File, "SCALARS v_l Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE v_l\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f\n", v_l[i]);

  // v_g
  fprintf(ptr_File, "SCALARS v_g Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE v_g\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f\n", v_g[i]);

  // v_hat
  fprintf(ptr_File, "SCALARS v_hat Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE v_hat\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f\n", v_hat[i]);

  // cell data
  fprintf(ptr_File, "CELL_DATA %u\n", n_Cell);
  // alpha
  fprintf(ptr_File, "SCALARS alpha Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE alpha\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", alpha[i]);

  // p_l
  fprintf(ptr_File, "SCALARS p_l Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p_l\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", p_l[i]);

  // p_g
  fprintf(ptr_File, "SCALARS p_g Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p_g\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", p_g[i]);

  // rho_l
  fprintf(ptr_File, "SCALARS rho_l Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE rho_l\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", rho_l[i]);

  // rho_g
  fprintf(ptr_File, "SCALARS rho_g Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE rho_g\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", rho_g[i]);

  // mu
  fprintf(ptr_File, "SCALARS mu Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE mu\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", mu[i]);

  // p_hat
  fprintf(ptr_File, "SCALARS p_hat Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p_hat\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", p_hat[i]);

  fclose(ptr_File);
}

void
FiveEqnTwoP_StagGrid::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOFs, n_DOFs, 25, NULL, &P_Mat);

  int n_Var = 5;
  PetscReal v = 1.0;
  for (int i = 0; i < n_Cell + 1; i++)
  {
    for (int var = 0; var < n_Var; var++)
    {
      PetscInt i_dof = i * n_Var + var;
      for (int j_dof = (i - 2) * n_Var; j_dof < (i + 3) * n_Var; j_dof++)
      {
        if ((i_dof >= 0) && (i_dof < n_DOFs) && (j_dof >= 0) && (j_dof < n_DOFs))
          MatSetValues(P_Mat, 1, &i_dof, 1, &j_dof, &v, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
  /*
  MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD);
  */
}

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "SinglePhaseFlow.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and S.J. Kim, Numerical study on the Welander oscillatory natural circulation
//       problem using high-order numerical methods, Progress in Nuclear Energy, 94 (2017) 162-172

/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0(v) 1(p) 3(v) 4(p)                                         3n(p)
        2(T)      5(T)
*/

// This problem has hard-coded boundary conditions: inlet velocity and T, outlet P.

SinglePhaseFlow::SinglePhaseFlow(InputParameterList & pList) :
  PETScProblem(pList)
{
  paramList.readRequiredInputParameter<int>("order");
  paramList.readRequiredInputParameter<int>("n_cells");
  paramList.readRequiredInputParameter<double>("length");

  paramList.readRequiredInputParameter<double>("P_INIT");
  paramList.readRequiredInputParameter<double>("V_INIT");
  paramList.readRequiredInputParameter<double>("T_INIT");

  paramList.readRequiredInputParameter<double>("V_INLET");
  paramList.readRequiredInputParameter<double>("T_INLET");
  paramList.readRequiredInputParameter<double>("P_OUTLET");
  paramList.readRequiredInputParameter<double>("T_OUTLET");

  // initial conditions
  P_INIT     =  paramList.getParameterValue<double>("P_INIT");
  V_INIT     =  paramList.getParameterValue<double>("V_INIT");
  T_INIT     =  paramList.getParameterValue<double>("T_INIT");
  // boundary conditions
  V_INLET    =  paramList.getParameterValue<double>("V_INLET");
  T_INLET    =  paramList.getParameterValue<double>("T_INLET");
  P_OUTLET   =  paramList.getParameterValue<double>("P_OUTLET");
  T_OUTLET   =  paramList.getParameterValue<double>("T_OUTLET");  /* safeguard for reverse flow */

  _order = paramList.getParameterValue<int>("order");
  n_Cell = paramList.getParameterValue<int>("n_cells");
  length = paramList.getParameterValue<double>("length");


  n_Node = n_Cell + 1;
  n_DOFs = n_Cell * 3 + 1;

  dx = length / n_Cell;

  // Primary variables
  p.resize(n_Cell);         v.resize(n_Node);           T.resize(n_Cell);
  p_old.resize(n_Cell);     v_old.resize(n_Node);       T_old.resize(n_Cell);
  p_oo.resize(n_Cell);      v_oo.resize(n_Node);        T_oo.resize(n_Cell);

  // Dependent and Helper variables
  rho.resize(n_Cell);       e.resize(n_Cell);
  rho_old.resize(n_Cell);   e_old.resize(n_Cell);
  rho_oo.resize(n_Cell);    e_oo.resize(n_Cell);
  rho_edge.resize(n_Node);  rho_edge_old.resize(n_Node);  rho_edge_oo.resize(n_Node);
  v_cell.resize(n_Cell);

  // Fluxes
  mass_flux.resize(n_Node);  energy_flux.resize(n_Node);

  // Second-order helper variables
  if (_order == 2)
  {
    p_w.resize(n_Cell);       T_w.resize(n_Cell);
    p_e.resize(n_Cell);       T_e.resize(n_Cell);
  }
}

SinglePhaseFlow::~SinglePhaseFlow()
{
}

void
SinglePhaseFlow::SetupInitialCondition(double * u)
{
  unsigned int index = 0;
  for(int i = 0; i < n_Cell + 1; i++)
  {
    v[i]       = V_INIT;
    u[index++] = v[i];

    if (i < n_Cell)
    {
      p[i]      = P_INIT;
      T[i]      = T_INIT;
      rho[i]    = rho_func(p[i], T[i]);
      e[i]      = e_func(p[i], T[i]);

      u[index++] = p[i];
      u[index++] = T[i];
    }
  }

  rho_old = rho;  v_old = v;  T_old = T;  e_old = e;
  rho_oo  = rho;  v_oo  = v;  T_oo  = T;  e_oo  = e;
}

void
SinglePhaseFlow::updateSolution(double * u, TimeStepIndex index)
{
  unsigned int idx = 0;
  switch (index)
  {
    case NEW:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v[i] = u[idx++];
        if (i < n_Cell)
        {
          p[i] = u[idx++];  T[i] = u[idx++];
          rho[i] = rho_func(p[i], T[i]);
          e[i]   = e_func(p[i], T[i]);
        }
      }
      //update cell values
      for(int i = 0; i < n_Cell; i++)
        v_cell[i] = 0.5 * (v[i] + v[i+1]);
      //update edge values
      rho_edge[0] = rho[0];
      for(unsigned int i = 1; i < n_Cell; i++)    rho_edge[i] = 0.5 * (rho[i-1] + rho[i]);
      rho_edge[n_Cell] = rho[n_Cell - 1];

      break;

    case OLD:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v_old[i] = u[idx++];
        if (i < n_Cell)
        {
          p_old[i] = u[idx++];  T_old[i] = u[idx++];
          rho_old[i] = rho_func(p_old[i], T_old[i]);
          e_old[i]   = e_func(p_old[i], T_old[i]);
        }
      }
      //update edge values
      rho_edge_old[0] = rho_old[0];
      for(unsigned int i = 1; i < n_Cell; i++)    rho_edge_old[i] = 0.5 * (rho_old[i-1] + rho_old[i]);
      rho_edge_old[n_Cell] = rho_old[n_Cell - 1];
      break;

    case OLDOLD:
      for(unsigned int i = 0; i < n_Cell + 1; i++)
      {
        v_oo[i] = u[idx++];
        if (i < n_Cell)
        {
          p_oo[i] = u[idx++];  T_oo[i] = u[idx++];
          rho_oo[i] = rho_func(p_oo[i], T_oo[i]);
          e_oo[i]   = e_func(p_oo[i], T_oo[i]);
        }
      }
      //update edge values
      rho_edge_oo[0] = rho_oo[0];
      for(unsigned int i = 1; i < n_Cell; i++)    rho_edge_oo[i] = 0.5 * (rho_oo[i-1] + rho_oo[i]);
      rho_edge_oo[n_Cell] = rho_oo[n_Cell - 1];
      break;

    default:
      sysError("TimeStepIndex should be one of NEW, OLD, and OLDOLD");
  }
}

void
SinglePhaseFlow::transientResidual(double * res)
{
  unsigned int idx = 0;
  if ((_time_scheme == BDF2) && (_step > 1))
  {
    for(unsigned int i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (1.5 * rho_edge[i] * v[i] - 2.0 * rho_edge_old[i] * v_old[i] + 0.5 * rho_edge_oo[i] * v_oo[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (1.5 * rho[i] - 2.0 * rho_old[i] + 0.5 * rho_oo[i]) / _dt;
        res[idx++] = (1.5 * rho[i] * e[i] - 2.0 * rho_old[i] * e_old[i] + 0.5 * rho_oo[i] * e_oo[i]) / _dt;
      }
    }
  }
  else
  {
    for(unsigned int i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (rho_edge[i] * v[i] - rho_edge_old[i] * v_old[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (rho[i] - rho_old[i]) * dx / _dt;
        res[idx++] = (rho[i] * e[i] - rho_old[i] * e_old[i]) / _dt;
      }
    }
  }

  // Zero-out transient residual, to apply Dirichlet BC at the inlet later
  res[0] = 0.0;
}

void
SinglePhaseFlow::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:     RHS_1st_order(rhs);     break;
    case 2:     RHS_2nd_order(rhs);     break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
SinglePhaseFlow::RHS_1st_order(double * rhs)
{
  // Boundary values and ghost values
  double p_inlet_ghost = p[0]; // projection better?

  // Upwind donor cell method for void fraction and mass balance equations
  mass_flux[0]    = (v[0] > 0) ? v[0] * rho_func(p_inlet_ghost, T_INLET) : v[0] * rho[0];
  energy_flux[0]  = (v[0] > 0) ? v[0] * rho_func(p_inlet_ghost, T_INLET) * e_func(p_inlet_ghost, T_INLET)
                    : v[0] * rho[0] * e[0];

  mass_flux[n_Cell]    = (v[n_Cell] > 0) ? v[n_Cell] * rho[n_Cell-1] : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET);
  energy_flux[n_Cell]  = (v[n_Cell] > 0) ? v[n_Cell] * rho[n_Cell-1] * e[n_Cell-1]
                         : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET) * e_func(P_OUTLET, T_OUTLET);

  for(unsigned int i = 1; i < n_Cell; i++)
  {
    mass_flux[i]    = (v[i] > 0) ? v[i] * rho[i-1]          : v[i] * rho[i];
    energy_flux[i]  = (v[i] > 0) ? v[i] * rho[i-1] * e[i-1] : v[i] * rho[i] * e[i];
  }

  // Momentum equations RHS
  double _f = 0.01; double _dh = 0.01;
  rhs[0] = v[0] - V_INLET;
  for(int i = 1; i < n_Cell + 1; i++) // loop on the remaining edges
  {
    // east and west velocities
    double v_east   = (i == n_Cell) ? v[n_Cell]     : v[i+1];
    double rho_east = (i == n_Cell) ? rho[n_Cell-1] : rho[i];
    double v_west = (i == 1) ? V_INLET : v[i-1];
    double rho_west = rho[i-1]; // always

    // d(rho v v)_dx term
    double drhovv_dx = (v[i] > 0) ? (rho_east * v[i] * v[i] - rho_west * v_west * v_west) / dx
                                    : (rho_east * v_east * v_east - rho_west * v[i] * v[i]) / dx;

    // dp_dx term
    double dp_dx = (i == n_Cell) ? (P_OUTLET - p[n_Cell-1])/dx*2. : (p[i] - p[i-1])/dx;

    // friction term
    double fric = (i == n_Cell) ? 0.5 * _f / _dh * rho[n_Cell-1] * v_cell[n_Cell-1] * std::fabs(v_cell[n_Cell-1])
                                : 0.5 * _f * 0.5 / _dh * rho[i-1] * v_cell[i-1] * std::fabs(v_cell[i-1])
                                  + 0.5 * _f * 0.5 / _dh * rho[i] * v_cell[i] * std::fabs(v_cell[i]);
    //double fric = 0.5 * rho_edge[i] * _f / _dh * v[i] * std::fabs(v[i]);

    // assemble RHS terms
    rhs[3*i] = -drhovv_dx - dp_dx - fric;
  }

  // RHS for mass and energy equations
  double _h = 2000.0; double _aw = 300.0; double Tw = 350.0;
  for(int i = 0; i < n_Cell; i++) //loop on cells
  {
    rhs[3*i+1] = -(mass_flux[i+1] - mass_flux[i]) / dx;
    double p_dv_dx = p[i] * (v[i+1] - v[i]) / dx;
    double src = _h * _aw * (Tw - T[i]);
    rhs[3*i+2] = -(energy_flux[i+1] - energy_flux[i]) / dx - p_dv_dx + src;
  }
}

void
SinglePhaseFlow::RHS_2nd_order(double * rhs)
{
  sysError("Not implemented.");
}

void
SinglePhaseFlow::updateFluxes()
{
}

void
SinglePhaseFlow::updateFluxes2ndOrder()
{
}

void
SinglePhaseFlow::writeVTKOutput(unsigned int step)
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

  // v
  fprintf(ptr_File, "SCALARS v_l Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE v\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f\n", v[i]);

  // cell data
  fprintf(ptr_File, "CELL_DATA %u\n", n_Cell);
  // p
  fprintf(ptr_File, "SCALARS alpha Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", p[i]);

  // T
  fprintf(ptr_File, "SCALARS p_l Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE T\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", T[i]);

  fclose(ptr_File);
}

void
SinglePhaseFlow::writeTextOutput(unsigned int step)
{
  FILE * ptr_File;
  std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".dat";
  ptr_File = fopen(file_name.c_str(), "w");

  // cell data
  fprintf(ptr_File, "Time = %20.6e\n", _t);
  fprintf(ptr_File, "#Cell data\n");
  fprintf(ptr_File, "%20s%20s%20s%20s%20s\n", "x", "p", "T", "rho", "v_cell");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, p[i], T[i], rho[i], v_cell[i]);

  // edge data
  fprintf(ptr_File, "#Edge data\n");
  fprintf(ptr_File, "%20s%20s%20s\n", "x", "v", "rho_edge");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%20.6e%20.6e%20.6e\n", i*dx, v[i], rho_edge[i]);

  fclose(ptr_File);
}

void
SinglePhaseFlow::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOFs, n_DOFs, 25, NULL, &P_Mat);

  int n_Var = 3;
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

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "SinglePhaseFlow.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and S.J. Kim, Numerical study on the Welander oscillatory natural circulation
//       problem using high-order numerical methods, Progress in Nuclear Energy, 94 (2017) 162-172
//
// The momentum equation is simplified to use the primitive form. In the end, we solve:
// 1) Mass Equation:      d(rho)/dt + d(rho*v)/dx = 0
// 2) Momentum Equation:  rho dv/dt + rho*v dv/dx + dp/dx + f/(2dh) * rho*v*|v| = 0
// 3) Energy Equation:    d(rho*e)/dt + d(rho*v*e)/dx + p dv/dx - h_w*a_w*(T_w-T) = 0

/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0(v) 1(p) 3(v) 4(p)                                         3n+1(v)
        2(T)      5(T)
*/

// This problem has hard-coded boundary conditions: inlet velocity and T, outlet P.

SinglePhaseFlow::SinglePhaseFlow(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  // initial conditions
  P_INIT     =  _inputParamList.getValueFromInput<double>("P_INIT");
  V_INIT     =  _inputParamList.getValueFromInput<double>("V_INIT");
  T_INIT     =  _inputParamList.getValueFromInput<double>("T_INIT");
  // boundary conditions
  V_INLET    =  _inputParamList.getValueFromInput<double>("V_INLET");
  T_INLET    =  _inputParamList.getValueFromInput<double>("T_INLET");
  P_OUTLET   =  _inputParamList.getValueFromInput<double>("P_OUTLET");
  T_OUTLET   =  _inputParamList.getValueFromInput<double>("T_OUTLET");  /* safeguard for reverse flow */

  _order = _inputParamList.getValueFromInput<int>("order");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");
  length = _inputParamList.getValueFromInput<double>("length");


  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 3 + 1;

  dx = length / n_Cell;

  // Primary variables
  p.resize(n_Cell);         v.resize(n_Node);           T.resize(n_Cell);
  p_old.resize(n_Cell);     v_old.resize(n_Node);       T_old.resize(n_Cell);
  p_oo.resize(n_Cell);      v_oo.resize(n_Node);        T_oo.resize(n_Cell);

  // Dependent and Helper variables
  rho.resize(n_Cell);       e.resize(n_Cell);
  rho_old.resize(n_Cell);   e_old.resize(n_Cell);
  rho_oo.resize(n_Cell);    e_oo.resize(n_Cell);
  rho_edge.resize(n_Node);

  // Fluxes
  mass_flux.resize(n_Node);  energy_flux.resize(n_Node);

  // Second-order helper variables
  if (_order == 2)
  {
    p_w.resize(n_Cell);       T_w.resize(n_Cell);       rho_w.resize(n_Cell);     e_w.resize(n_Cell);
    p_e.resize(n_Cell);       T_e.resize(n_Cell);       rho_e.resize(n_Cell);     e_e.resize(n_Cell);
  }
}

SinglePhaseFlow::~SinglePhaseFlow() {}

void
SinglePhaseFlow::SetupInitialCondition(double * u)
{
  unsigned index = 0;
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

  rho_edge[0] = rho[0];   rho_edge[n_Cell] = rho[n_Cell - 1];
  for(unsigned i = 1; i < n_Cell; i++)    rho_edge[i] = 0.5 * (rho[i-1] + rho[i]);
}

void
SinglePhaseFlow::updateSolution(double * u)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < n_Cell + 1; i++)
  {
    v[i] = u[idx++];
    if (i < n_Cell)
    {
      p[i] = u[idx++];  T[i] = u[idx++];
      rho[i] = rho_func(p[i], T[i]);
      e[i]   = e_func(p[i], T[i]);
    }
  }
  //update edge values
  rho_edge[0] = rho[0];   rho_edge[n_Cell] = rho[n_Cell - 1];
  for(unsigned i = 1; i < n_Cell; i++)    rho_edge[i] = 0.5 * (rho[i-1] + rho[i]);
}

void
SinglePhaseFlow::transientResidual(double * res)
{
  unsigned idx = 0;
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
    {
      //res[idx++] = rho_edge[i] * (1.5 * v[i] - 2.0 * v_old[i] + 0.5 * v_oo[i]) / _dt;
      res[idx++] = rho_edge[i] * UTILS::BDF2Tran(v[i], v_old[i], v_oo[i], _dt, _dt_old);
      if (i < n_Cell)
      {
        //res[idx++] = (1.5 * rho[i] - 2.0 * rho_old[i] + 0.5 * rho_oo[i]) / _dt;
        //res[idx++] = (1.5 * rho[i] * e[i] - 2.0 * rho_old[i] * e_old[i] + 0.5 * rho_oo[i] * e_oo[i]) / _dt;
        res[idx++] = UTILS::BDF2Tran(rho[i], rho_old[i], rho_oo[i], _dt, _dt_old);
        res[idx++] = UTILS::BDF2Tran(rho[i] * e[i], rho_old[i] * e_old[i], rho_oo[i] * e_oo[i], _dt, _dt_old);
      }
    }
  }
  else
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = rho_edge[i] * (v[i] - v_old[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (rho[i] - rho_old[i]) / _dt;
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
    case 1:     updateFluxes();             break;
    case 2:     updateFluxes2ndOrder();     break;
    default :   sysError("Spatial order not implemented.");
  }

  // Momentum equations RHS
  double _f = 0.01; double _dh = 0.01;
  rhs[0] = v[0] - V_INLET;
  for(int i = 1; i < n_Cell + 1; i++) // loop on the remaining edges
  {
    // east and west velocities
    double v_east = (i == n_Cell) ? v[n_Cell] : v[i+1];
    double v_west = v[i-1];
    double dv_dx = (v[i] > 0) ? (v[i] - v_west) / dx : (v_east - v[i]) / dx;

    // dp_dx term
    double dp_dx = (i == n_Cell) ? (P_OUTLET - p[n_Cell-1])/dx*2. : (p[i] - p[i-1])/dx;

    // friction term
    double fric  = 0.5 * _f / _dh * rho_edge[i] * v[i] * std::fabs(v[i]);

    // assemble RHS terms
    rhs[3*i] = -rho_edge[i] * v[i] * dv_dx - dp_dx - fric;
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
SinglePhaseFlow::updateFluxes()
{
  double rho_inlet = rho_func(p[0], T_INLET); // better to use 1st-order extrapolation for p_inlet?
  double e_inlet = e_func(p[0], T_INLET);

  // Upwind donor cell method for void fraction and mass balance equations
  mass_flux[0]    = (v[0] > 0) ? v[0] * rho_inlet           : v[0] * rho[0];
  energy_flux[0]  = (v[0] > 0) ? v[0] * rho_inlet * e_inlet : v[0] * rho[0] * e[0];

  mass_flux[n_Cell]    = (v[n_Cell] > 0) ? v[n_Cell] * rho[n_Cell-1] : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET);
  energy_flux[n_Cell]  = (v[n_Cell] > 0) ? v[n_Cell] * rho[n_Cell-1] * e[n_Cell-1]
                         : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET) * e_func(P_OUTLET, T_OUTLET);

  for(unsigned i = 1; i < n_Cell; i++)
  {
    mass_flux[i]    = (v[i] > 0) ? v[i] * rho[i-1]          : v[i] * rho[i];
    energy_flux[i]  = (v[i] > 0) ? v[i] * rho[i-1] * e[i-1] : v[i] * rho[i] * e[i];
  }
}

void
SinglePhaseFlow::updateFluxes2ndOrder()
{
  double rho_inlet = rho_func(p[0], T_INLET); // better to use 1st-order extrapolation for p_inlet?
  double e_inlet = e_func(p[0], T_INLET);

  double p_inlet_GHOST = 2 * p[0] - p[1];
  double T_inlet_GHOST = 2 * T_INLET - T[0]; //double T_inlet_GHOST = 2 * T[0] - T[1];
  double p_outlet_GHOST = 2 * P_OUTLET - p[n_Cell-1]; //double p_outlet_GHOST = 2 * p[n_Cell-1] - p[n_Cell - 2];
  double T_outlet_GHOST = (v[n_Cell] > 0) ? 2 * T[n_Cell-1] - T[n_Cell-2] : 2 * T_OUTLET - T[n_Cell-1];
  //double T_outlet_GHOST = 2 * T[n_Cell-1] - T[n_Cell - 2];

  UTILS::linearReconstruction(p_inlet_GHOST, p_outlet_GHOST, p, p_w, p_e);
  UTILS::linearReconstruction(T_inlet_GHOST, T_outlet_GHOST, T, T_w, T_e);

  for (unsigned i = 0; i < n_Cell; i++)
  {
    rho_w[i] = rho_func(p_w[i], T_w[i]);  rho_e[i] = rho_func(p_e[i], T_e[i]);
    e_w[i]   = e_func(p_w[i], T_w[i]);    e_e[i]   = e_func(p_e[i], T_e[i]);
  }

  // Upwind donor cell method for void fraction and mass balance equations
  mass_flux[0]    = (v[0] > 0) ? v[0] * rho_inlet           : v[0] * rho_w[0];
  energy_flux[0]  = (v[0] > 0) ? v[0] * rho_inlet * e_inlet : v[0] * rho_w[0] * e_w[0];

  mass_flux[n_Cell]    = (v[n_Cell] > 0) ? v[n_Cell] * rho_e[n_Cell-1] : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET);
  energy_flux[n_Cell]  = (v[n_Cell] > 0) ? v[n_Cell] * rho_e[n_Cell-1] * e_e[n_Cell-1]
                         : v[n_Cell] * rho_func(P_OUTLET, T_OUTLET) * e_func(P_OUTLET, T_OUTLET);

  for(unsigned i = 1; i < n_Cell; i++)
  {
    mass_flux[i]    = (v[i] > 0) ? v[i] * rho_e[i-1]            : v[i] * rho_w[i];
    energy_flux[i]  = (v[i] > 0) ? v[i] * rho_e[i-1] * e_e[i-1] : v[i] * rho_w[i] * e_w[i];
  }
}

void
SinglePhaseFlow::onTimestepEnd()
{
  // save old solutions
  p_oo = p_old; v_oo = v_old; T_oo = T_old; rho_oo = rho_old; e_oo = e_old;
  p_old = p; v_old = v; T_old = T; rho_old = rho; e_old = e;
}

void
SinglePhaseFlow::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "T");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", T[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </CellData>\n");

  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n", "node_id");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %d\n", i);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", v[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
SinglePhaseFlow::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s%20s\n", "x", "p", "T", "rho", "v_cell");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, p[i], T[i], rho[i], 0.5 * (v[i] + v[i+1]));

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s%20s\n", "x", "v", "rho_edge");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "%20.6e%20.6e%20.6e\n", i*dx, v[i], rho_edge[i]);
}

void
SinglePhaseFlow::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  int n_Var = 3;
  for (int i = 0; i < n_Cell + 1; i++)
  {
    for (int var = 0; var < n_Var; var++)
    {
      int i_dof = i * n_Var + var;
      for (int j_dof = (i - 2) * n_Var; j_dof < (i + 3) * n_Var; j_dof++)
      {
        if ((i_dof >= 0) && (i_dof < _n_DOFs) && (j_dof >= 0) && (j_dof < _n_DOFs))
          mnzp->addEntry(i_dof + _DOF_offset, j_dof + _DOF_offset);
      }
    }
  }
}

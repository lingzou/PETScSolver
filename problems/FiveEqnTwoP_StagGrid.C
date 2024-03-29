#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "FiveEqnTwoP_StagGrid.h"
#include "utils.h"

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

FiveEqnTwoP_StagGrid::FiveEqnTwoP_StagGrid(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
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

  _order = _inputParamList.getValueFromInput<int>("order");
  H_inv  = _inputParamList.getValueFromInput<double>("H_inv");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");

  length = 12.0;
  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 5 + 2;

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

  // Second-order helper variables
  if (_order == 2)
  {
    alpha_w.resize(n_Cell);       p_l_w.resize(n_Cell);       p_g_w.resize(n_Cell);       v_l_w.resize(n_Node);       v_g_w.resize(n_Node);
    alpha_e.resize(n_Cell);       p_l_e.resize(n_Cell);       p_g_e.resize(n_Cell);       v_l_e.resize(n_Node);       v_g_e.resize(n_Node);
  }
}

FiveEqnTwoP_StagGrid::~FiveEqnTwoP_StagGrid()
{
}

void
FiveEqnTwoP_StagGrid::SetupInitialCondition(double * u)
{
  unsigned index = 0;
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
FiveEqnTwoP_StagGrid::updateSolution(double * u)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < n_Cell + 1; i++)
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
  for(unsigned i = 1; i < n_Cell; i++)
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

  for(unsigned i = 0; i < n_Cell + 1; i++)
    v_hat[i] = alpha_edge[i] * v_g[i] + (1.0 - alpha_edge[i]) * v_l[i];
}

void
FiveEqnTwoP_StagGrid::transientResidual(double * res)
{
  unsigned idx = 0;
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (1.0 - alpha_edge[i]) * rho_l_edge[i] * UTILS::BDF2Tran(v_l[i], v_l_old[i], v_l_oo[i], _dt, _dt_old);
      res[idx++] = alpha_edge[i] * rho_g_edge[i] * UTILS::BDF2Tran(v_g[i], v_g_old[i], v_g_oo[i], _dt, _dt_old);
      if (i < n_Cell)
      {
        res[idx++] = UTILS::BDF2Tran(alpha[i], alpha_old[i], alpha_oo[i], _dt, _dt_old);
        res[idx++] = UTILS::BDF2Tran( (1-alpha[i])     * rho_l[i],
                                      (1-alpha_old[i]) * rho_l_old[i],
                                      (1-alpha_oo[i])  * rho_l_oo[i], _dt, _dt_old);
        res[idx++] = UTILS::BDF2Tran(alpha[i]     * rho_g[i],
                                     alpha_old[i] * rho_g_old[i],
                                     alpha_oo[i]  * rho_g_oo[i], _dt, _dt_old);
      }
    }
  }
  else
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
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
  switch (_order)
  {
    case 1:     RHS_1st_order(rhs);     break;
    case 2:     RHS_2nd_order(rhs);     break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
FiveEqnTwoP_StagGrid::RHS_1st_order(double * rhs)
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

  for(unsigned i = 1; i < n_Cell; i++)
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
FiveEqnTwoP_StagGrid::RHS_2nd_order(double * rhs)
{
  // 1. Linear reconstruction with linear-extrapolation boundary ghost cell values
  // 1.1 void fraction
  double alpha_inlet_ghost = 2.0 * ALPHA_INLET - alpha[0];
  double alpha_outlet_ghost = (v_l[n_Cell] > 0) ? 2.0 * alpha[n_Cell - 1] - alpha[n_Cell - 2] : ALPHA_OUTLET;
  UTILS::linearReconstruction(alpha_inlet_ghost, alpha_outlet_ghost, alpha, alpha_w, alpha_e);
  // 1.2 liquid pressures
  double p_l_inlet_ghost = 2.0 * p_l[0] - p_l[1];
  double p_l_out_ghost = 2.0 * P_OUTLET - p_l[n_Cell - 1];
  UTILS::linearReconstruction(p_l_inlet_ghost, p_l_out_ghost, p_l, p_l_w, p_l_e);
  // 1.3 vapor pressures
  double p_g_inlet_ghost = 2.0 * p_g[0] - p_g[1];
  double p_g_out_ghost = 2.0 * P_OUTLET - p_g[n_Cell - 1];
  UTILS::linearReconstruction(p_l_inlet_ghost, p_l_out_ghost, p_l, p_l_w, p_l_e);
  // 1.4 velocities
  //    v_l(g)_inlet_ghost are not used here, as the inlet will be treated differently, see the residual forms
  double v_l_outlet_ghost = 2. * v_l[n_Cell - 1] - v_l[n_Cell - 2];
  double v_g_outlet_ghost = 2. * v_g[n_Cell - 1] - v_g[n_Cell - 2];
  UTILS::linearReconstruction(V_L_INIT, v_l_outlet_ghost, v_l, v_l_w, v_l_e);
  UTILS::linearReconstruction(V_G_INIT, v_g_outlet_ghost, v_g, v_g_w, v_g_e);

  // 2. Void fraction and mass equation fluxes
  // Upwind donor cell method for void fraction and mass balance equations
  double rho_l_inlet_bc = rho_l_func(p_l_inlet_ghost);
  double rho_g_inlet_bc = rho_g_func(p_g_inlet_ghost);
  double rho_l_outlet_bc = rho_l_func(P_OUTLET);
  double rho_g_outlet_bc = rho_g_func(P_OUTLET);
  alpha_flux[0] = (v_hat[0] > 0) ? v_hat[0] * ALPHA_INLET : v_hat[0] * alpha_w[0];
  rho_l_flux[0] = (V_L_INIT > 0) ? V_L_INIT * (1.0 - ALPHA_INLET) * rho_l_inlet_bc : V_L_INIT * (1.0 - alpha_w[0]) * rho_l_func(p_l_w[0]);
  rho_g_flux[0] = (V_G_INIT > 0) ? V_G_INIT * ALPHA_INLET * rho_g_inlet_bc : V_G_INIT * alpha_w[0] * rho_g_func(p_g_w[0]);

  alpha_flux[n_Cell] = (v_hat[n_Cell] > 0) ? v_hat[n_Cell] * alpha_e[n_Cell-1] : v_hat[n_Cell] * ALPHA_OUTLET;
  rho_l_flux[n_Cell] = (v_l[n_Cell] > 0) ? v_l[n_Cell] * (1. - alpha_e[n_Cell-1]) * rho_l_func(p_l_e[n_Cell-1]) : v_l[n_Cell] * (1.0 - ALPHA_OUTLET) * rho_l_outlet_bc;
  rho_g_flux[n_Cell] = (v_g[n_Cell] > 0) ? v_g[n_Cell] * alpha_e[n_Cell-1] * rho_g_func(p_g_e[n_Cell-1]) : v_g[n_Cell] * alpha_outlet_ghost * rho_g_outlet_bc;

  for(unsigned i = 1; i < n_Cell; i++)
  {
    alpha_flux[i] = (v_hat[i] > 0) ? v_hat[i] * alpha_e[i-1] : v_hat[i] * alpha_w[i];
    rho_l_flux[i] = (v_l[i] > 0) ? v_l[i] * (1. - alpha_e[i-1]) * rho_l_func(p_l_e[i-1]) : v_l[i] * (1. - alpha_w[i]) * rho_l_func(p_l_w[i]);
    rho_g_flux[i] = (v_g[i] > 0) ? v_g[i] * alpha_e[i-1] * rho_g_func(p_g_e[i-1]) : v_g[i] * alpha_w[i] * rho_g_func(p_g_w[i]);
  }

  // 3. RHS
  // Momentum equations RHS; Eqn. (22) of Ref. [1]
  rhs[0] = v_l[0] - V_L_INIT;
  rhs[1] = v_g[0] - V_G_INIT;
  for(int i = 1; i < n_Cell + 1; i++) // loop on the remaining edges
  {
    // dv_dx term
    double dv_l_dx = 0.;
    double dv_g_dx = 0.;
    if(v_l[i] > 0)
    {
      double v_l_west = (i == 1) ? 0.5*(V_L_INIT + v_l[1]) : v_l_e[i-1];
      double v_l_east = v_l_e[i];
      dv_l_dx = (v_l_east - v_l_west) / dx;
    }
    else
    {
      double v_l_east = (i == n_Cell-1) ? v_l_outlet_ghost : v_l_w[i+1];
      double v_l_west = v_l_w[i];
      dv_l_dx = (v_l_east - v_l_west) / dx;
    }

    if(v_g[i] > 0)
    {
      double v_g_west = (i == 1) ? 0.5*(V_G_INIT + v_g[1]) : v_g_e[i-1];
      double v_g_east = v_g_e[i];
      dv_g_dx = (v_g_east - v_g_west) / dx;
    }
    else
    {
      double v_g_east = (i == n_Cell-1) ? v_g_outlet_ghost : v_g_w[i+1];
      double v_g_west = v_g_w[i];
      dv_g_dx = (v_g_east - v_g_west) / dx;
    }

    // d(alpha*p)_dx term, the same as 1-st order, central differencing
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
FiveEqnTwoP_StagGrid::onTimestepEnd()
{
  // save old solutions
  alpha_oo  = alpha_old;  p_l_oo  = p_l_old;   p_g_oo  = p_g_old;   v_l_oo  = v_l_old;   v_g_oo  = v_g_old;   rho_l_oo  = rho_l_old;   rho_g_oo  = rho_g_old;
  alpha_old = alpha;      p_l_old = p_l;       p_g_old = p_g;       v_l_old = v_l;       v_g_old = v_g;       rho_l_old = rho_l;       rho_g_old = rho_g;
}

void
FiveEqnTwoP_StagGrid::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "alpha");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", alpha[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p_l");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p_g[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_l");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", rho_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", rho_g[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "mu");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", mu[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p_hat");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p_hat[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </CellData>\n");

  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n", "node_id");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %d\n", i);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v_l");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", v_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v_g");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", v_g[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v_hat");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %20.6f\n", v_hat[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
FiveEqnTwoP_StagGrid::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s%20s%20s%20s%20s\n", "x", "alpha", "p_l", "p_g", "rho_l", "rho_g", "mu", "p_hat");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, alpha[i], p_l[i], p_g[i], rho_l[i], rho_g[i], mu[i], p_hat[i]);

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s%20s%20s\n", "x", "v_l", "v_g", "v_hat");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e\n", i*dx, v_l[i], v_g[i], v_hat[i]);
}

void
FiveEqnTwoP_StagGrid::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  int n_Var = 5;
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

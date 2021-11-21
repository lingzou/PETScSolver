#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "FourEqnOneP.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and H. Zhang,
//        Solving phase appearance/disappearance two-phase flow problems with high resolution staggered grid and
//        fully implicit schemes by the Jacobian-free Newton–Krylov Method
//        Computers & Fluids, Vol. 129, 179-188


/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0    2    4    6                                            4n
   1    3    5    7                                            4n+1
  v_l  alpha
  v_g  p
*/

FourEqnOneP::FourEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _problem_type = _inputParamList.getValueFromInput<int>("problem_type");
  _order = _inputParamList.getValueFromInput<int>("order");
  length = _inputParamList.getValueFromInput<int>("length");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");

  // The sine wave and square wave problems use periodic BC
  _periodic_bc = (_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE);

  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 4 + 2;

  dx = length / n_Cell;

  // Primary variables
  alpha.resize(n_Cell);       p.resize(n_Cell);         v_l.resize(n_Node);       v_g.resize(n_Node);
  alpha_old.resize(n_Cell);   p_old.resize(n_Cell);     v_l_old.resize(n_Node);   v_g_old.resize(n_Node);
  alpha_oo.resize(n_Cell);    p_oo.resize(n_Cell);      v_l_oo.resize(n_Node);    v_g_oo.resize(n_Node);

  // Helper variables
  rho_l.resize(n_Cell);       rho_g.resize(n_Cell);
  rho_l_old.resize(n_Cell);   rho_g_old.resize(n_Cell);
  rho_l_oo.resize(n_Cell);    rho_g_oo.resize(n_Cell);

  alpha_edge.resize(n_Node);  rho_l_edge.resize(n_Node);  rho_g_edge.resize(n_Node);
  v_l_cell.resize(n_Cell);    v_g_cell.resize(n_Cell);

  // Fluxes
  rho_l_flux.resize(n_Node);  rho_g_flux.resize(n_Node);

  // Second-order helper variables
  if (_order == 2)
  {
    alpha_w.resize(n_Cell);       p_w.resize(n_Cell);       v_l_w.resize(n_Node);       v_g_w.resize(n_Node);
    alpha_e.resize(n_Cell);       p_e.resize(n_Cell);       v_l_e.resize(n_Node);       v_g_e.resize(n_Node);
  }
}

FourEqnOneP::~FourEqnOneP()
{
}

void
FourEqnOneP::SetupInitialCondition(double * u)
{
  unsigned index = 0;
  switch (_problem_type)
  {
    case SINE_WAVE:
      for(int i = 0; i < n_Cell + 1; i++)
      {
        v_l[i] = 1.0;
        v_g[i] = 1.0;
        u[index++] = v_l[i];
        u[index++] = v_g[i];

        if (i < n_Cell)
        {
          double xx = (i + 0.5) * dx;
          alpha[i] = 0.5 + 0.2 * std::sin(2 * PI * xx / length);
          p[i] = 1.e5;
          rho_l[i] = rho_l_func(p[i]);
          rho_g[i] = rho_g_func(p[i]);

          u[index++] = alpha[i];
          u[index++] = p[i];
        }
      }
    break;

    case SQUARE_WAVE:
      for(int i = 0; i < n_Cell + 1; i++)
      {
        v_l[i] = 1.0;
        v_g[i] = 1.0;
        u[index++] = v_l[i];
        u[index++] = v_g[i];

        if (i < n_Cell)
        {
          double xx = (i + 0.5) * dx;
          alpha[i] = ((xx < 0.2) || (xx > 0.4)) ? 1 : 0;
          p[i] = 1.e5;
          rho_l[i] = rho_l_func(p[i]);
          rho_g[i] = rho_g_func(p[i]);

          u[index++] = alpha[i];
          u[index++] = p[i];
        }
      }
    break;

    case MANOMETER:
      for(int i = 0; i < n_Cell + 1; i++)
      {
        v_l[i] = -1.0;
        v_g[i] = -1.0;
        u[index++] = v_l[i];
        u[index++] = v_g[i];

        if (i < n_Cell)
        {
          double xx = (i + 0.5) * dx;
          alpha[i] = ((xx < 5) || (xx > 15)) ? 1 : 0;
          p[i] = 1.e5;
          rho_l[i] = rho_l_func(p[i]);
          rho_g[i] = rho_g_func(p[i]);

          u[index++] = alpha[i];
          u[index++] = p[i];
        }
      }
    break;

    case SEDIMENTATION:
    break;

    case WATER_FAUCET:
    break;

    defaut: sysError("Unknown problem_type.");
  }

  alpha_old = alpha;  p_old = p;  v_l_old = v_l;  v_g_old = v_g;
  alpha_oo  = alpha;  p_oo  = p;  v_l_oo  = v_l;  v_g_oo  = v_g;

  rho_l_old = rho_l;  rho_g_old = rho_g;
  rho_l_oo  = rho_l;  rho_g_oo  = rho_g;

  ALPHA_MIN = 1e-6;
}

void
FourEqnOneP::updateSolution(double * u)
{
  unsigned idx = 0;
  for(int i = 0; i < n_Cell + 1; i++)
  {
    v_l[i] = u[idx++];
    v_g[i] = u[idx++];

    if (i < n_Cell)
    {
      alpha[i] = u[idx++];
      p[i] = u[idx++];

      rho_l[i] = rho_l_func(p[i]);
      rho_g[i] = rho_g_func(p[i]);
    }
  }

  //update cell values
  for(int i = 0; i < n_Cell; i++)
  {
    v_l_cell[i] = 0.5 * (v_l[i] + v_l[i+1]);
    v_g_cell[i] = 0.5 * (v_g[i] + v_g[i+1]);
  }

  //update edge values
  rho_l_edge[0] = _periodic_bc ? 0.5 * (rho_l[0] + rho_l[n_Cell-1]) : rho_l[0];
  rho_g_edge[0] = _periodic_bc ? 0.5 * (rho_g[0] + rho_g[n_Cell-1]) : rho_g[0];
  alpha_edge[0] = _periodic_bc ? 0.5 * (alpha[0] + alpha[n_Cell-1]) : alpha[0];
  for(unsigned i = 1; i < n_Cell; i++)
  {
    rho_l_edge[i] = 0.5 * (rho_l[i-1] + rho_l[i]);
    rho_g_edge[i] = 0.5 * (rho_g[i-1] + rho_g[i]);
    alpha_edge[i] = 0.5 * (alpha[i-1] + alpha[i]);
  }
  rho_l_edge[n_Cell] = _periodic_bc ? 0.5 * (rho_l[0] + rho_l[n_Cell-1]) : rho_l[n_Cell - 1];
  rho_g_edge[n_Cell] = _periodic_bc ? 0.5 * (rho_g[0] + rho_g[n_Cell-1]) : rho_g[n_Cell - 1];
  alpha_edge[n_Cell] = _periodic_bc ? 0.5 * (alpha[0] + alpha[n_Cell-1]) : alpha[n_Cell - 1];
}

void
FourEqnOneP::transientResidual(double * res)
{
  unsigned idx = 0;
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (1.5 * v_l[i] - 2.0 * v_l_old[i] + 0.5 * v_l_oo[i]) / _dt;
      res[idx++] = (1.5 * v_g[i] - 2.0 * v_g_old[i] + 0.5 * v_g_oo[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = (1.5 * (1.0 - alpha[i]) * rho_l[i] - 2.0 * (1.0 - alpha_old[i]) * rho_l_old[i] + 0.5 * (1.0 - alpha_oo[i]) * rho_l_oo[i]) / _dt;
        res[idx++] = (1.5 * alpha[i] * rho_g[i] - 2.0 * alpha_old[i] * rho_g_old[i] + 0.5 * alpha_oo[i] * rho_g_oo[i]) / _dt;
      }
    }
  }
  else
  {
    for(unsigned i = 0; i < n_Cell + 1; i++)
    {
      res[idx++] = (v_l[i] - v_l_old[i]) / _dt;
      res[idx++] = (v_g[i] - v_g_old[i]) / _dt;
      if (i < n_Cell)
      {
        res[idx++] = ((1.0 - alpha[i]) * rho_l[i] - (1.0 - alpha_old[i]) * rho_l_old[i]) / _dt;
        res[idx++] = (alpha[i] * rho_g[i]         - alpha_old[i] * rho_g_old[i])         / _dt;
      }
    }
  }
}

void
FourEqnOneP::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:     RHS_1st_order(rhs);     break;
    case 2:     RHS_2nd_order(rhs);     break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
FourEqnOneP::RHS_1st_order(double * rhs)
{
  // Boundary values and ghost values
  double p_inlet_bc = _periodic_bc ? p[n_Cell-1] : p[0];
  double p_outlet_bc = _periodic_bc ? p[0] : p[n_Cell-1];
  double alpha_inlet_bc = _periodic_bc ? alpha[n_Cell-1] : alpha[0];
  double alpha_outlet_bc = _periodic_bc ? alpha[0] : alpha[n_Cell-1];
  double v_l_inlet_bc = _periodic_bc ? v_l[n_Cell-1] : v_l[0];
  double v_l_outlet_bc = _periodic_bc ? v_l[1] : v_l[n_Cell];
  double v_g_inlet_bc = _periodic_bc ? v_g[n_Cell-1] : v_g[0];
  double v_g_outlet_bc = _periodic_bc ? v_g[1] : v_g[n_Cell];

  double rho_l_inlet_bc = rho_l_func(p_inlet_bc);
  double rho_g_inlet_bc = rho_g_func(p_inlet_bc);
  double rho_l_outlet_bc = rho_l_func(p_outlet_bc);
  double rho_g_outlet_bc = rho_g_func(p_outlet_bc);

  // Upwind donor cell method for mass balance equations
  rho_l_flux[0] = (v_l[0] > 0) ? v_l[0] * (1.0 - alpha_inlet_bc) * rho_l_inlet_bc : v_l[0] * (1.0 - alpha[0]) * rho_l[0];
  rho_g_flux[0] = (v_g[0] > 0) ? v_g[0] * alpha_inlet_bc * rho_g_inlet_bc : v_g[0] * alpha[0] * rho_g[0];
  rho_l_flux[n_Cell] = (v_l[n_Cell] > 0) ? v_l[n_Cell] * (1. - alpha[n_Cell - 1]) * rho_l[n_Cell - 1] : v_l[n_Cell] * (1.0 - alpha_outlet_bc) * rho_l_outlet_bc;
  rho_g_flux[n_Cell] = (v_g[n_Cell] > 0) ? v_g[n_Cell] * alpha[n_Cell - 1] * rho_g[n_Cell - 1] : v_g[n_Cell] * alpha_outlet_bc * rho_g_outlet_bc;

  for(unsigned i = 1; i < n_Cell; i++)
  {
    rho_l_flux[i] = (v_l[i] > 0) ? v_l[i] * (1 - alpha[i-1]) * rho_l[i-1] : v_l[i] * (1 - alpha[i]) * rho_l[i];
    rho_g_flux[i] = (v_g[i] > 0) ? v_g[i] * alpha[i-1] * rho_g[i-1] : v_g[i] * alpha[i] * rho_g[i];
  }

  // Momentum equations RHS
  for(int i = 0; i < n_Cell + 1; i++)
  {
    // east and west velocities
    double v_l_east = (i == n_Cell) ? v_l_outlet_bc : v_l[i+1];
    double v_l_west = (i == 0) ? v_l_inlet_bc : v_l[i-1];
    double v_g_east = (i == n_Cell) ? v_g_outlet_bc : v_g[i+1];
    double v_g_west = (i == 0) ? v_g_inlet_bc : v_g[i-1];

    // dv_dx term
    double dv_l_dx = (v_l[i] > 0) ? (v_l[i] - v_l_west)/dx : (v_l_east - v_l[i])/dx;
    double dv_g_dx = (v_g[i] > 0) ? (v_g[i] - v_g_west)/dx : (v_g_east - v_g[i])/dx;

    // dp_dx term
    double dp_dx = 0;
    if (i == 0)               dp_dx = (p[0] - p_inlet_bc) / dx * 2.;
    else if (i == n_Cell)     dp_dx = (p_outlet_bc - p[n_Cell-1])/dx * 2.;
    else                      dp_dx = (p[i] - p[i-1])/dx;

    // assemble RHS for liquid and gas phase momentum equations
    double rhs_v_l = -v_l[i] * dv_l_dx - dp_dx / rho_l_edge[i];
    double rhs_v_g = -v_g[i] * dv_g_dx - dp_dx / rho_g_edge[i];

    // gravity terms, alpha * rho * g
    rhs[4*i] = rhs_v_l;
    rhs[4*i+1] = rhs_v_g;
  }

  // RHS for mass conservation equations
  for(int i = 0; i < n_Cell; i++) //loop on cells
  {
    // mass conservation equations
    rhs[4*i+2] = -(rho_l_flux[i+1] - rho_l_flux[i]) / dx;
    rhs[4*i+3] = -(rho_g_flux[i+1] - rho_g_flux[i]) / dx;
  }
}

void
FourEqnOneP::RHS_2nd_order(double * rhs)
{
  double p_inlet_ghost, p_outlet_ghost, alpha_inlet_ghost, alpha_outlet_ghost, v_l_inlet_ghost, v_l_outlet_ghost, v_g_inlet_ghost, v_g_outlet_ghost;
  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE)) // periodic BCs
  {
    p_inlet_ghost = p[n_Cell-1];
    p_outlet_ghost = p[0];
    alpha_inlet_ghost = alpha[n_Cell-1];
    alpha_outlet_ghost = alpha[0];
    v_l_inlet_ghost = v_l[n_Cell-1];
    v_l_outlet_ghost = v_l[1];
    v_g_inlet_ghost = v_g[n_Cell-1];
    v_g_outlet_ghost = v_g[1];
  }
  else
  {
    p_inlet_ghost = 2 * p[0] - p[1];
    p_outlet_ghost = 2 * p[n_Cell-1] - p[n_Cell-2];
    alpha_inlet_ghost = 2 * alpha[0] - alpha[1];
    alpha_outlet_ghost = 2 * alpha[n_Cell-1] - alpha[n_Cell-2];
    v_l_inlet_ghost = 2 * v_l[0] - v_l[1];
    v_l_outlet_ghost = 2 * v_l[n_Cell] - v_l[n_Cell-1];
    v_g_inlet_ghost = 2 * v_g[0] - v_g[1];
    v_g_outlet_ghost = 2 * v_g[n_Cell] - v_g[n_Cell-1];
  }

  UTILS::linearReconstruction(p_inlet_ghost, p_outlet_ghost, p, p_w, p_e);
  UTILS::linearReconstruction(alpha_inlet_ghost, alpha_outlet_ghost, alpha, alpha_w, alpha_e);
  UTILS::linearReconstruction(v_l_inlet_ghost, v_l_outlet_ghost, v_l, v_l_w, v_l_e);
  UTILS::linearReconstruction(v_g_inlet_ghost, v_g_outlet_ghost, v_g, v_g_w, v_g_e);

  double p_inlet_bc, p_outlet_bc, alpha_inlet_bc, alpha_outlet_bc, rho_l_inlet_bc, rho_l_outlet_bc, rho_g_inlet_bc, rho_g_outlet_bc;
  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
  {
    p_inlet_bc = p[n_Cell-1];
    p_outlet_bc = p[0];
    alpha_inlet_bc = alpha_e[n_Cell-1];
    alpha_outlet_bc = alpha_w[0];
    rho_l_inlet_bc = rho_l_func(p_e[n_Cell-1]);
    rho_l_outlet_bc = rho_l_func(p_w[0]);
    rho_g_inlet_bc = rho_g_func(p_e[n_Cell-1]);
    rho_g_outlet_bc = rho_g_func(p_w[0]);
  }
  else
  {
    p_inlet_bc = 1e5;
    p_outlet_bc = 1e5;
    alpha_inlet_bc = 1;
    alpha_outlet_bc = 1;
    rho_l_inlet_bc = rho_l_func(1e5);
    rho_l_outlet_bc = rho_l_func(1e5);
    rho_g_inlet_bc = rho_g_func(1e5);
    rho_g_outlet_bc = rho_g_func(1e5);
  }

  rho_l_flux[0] = (v_l[0] > 0) ? v_l[0] * (1 - alpha_inlet_bc) * rho_l_inlet_bc : v_l[0] * (1 - alpha_w[0]) * rho_l_func(p_w[0]);
  rho_g_flux[0] = (v_g[0] > 0) ? v_g[0] * alpha_inlet_bc * rho_g_inlet_bc : v_g[0] * alpha_w[0] * rho_g_func(p_w[0]);
  rho_l_flux[n_Cell] = (v_l[n_Cell] > 0) ? v_l[n_Cell] * (1 - alpha_e[n_Cell - 1]) * rho_l_func(p_e[n_Cell -1]) : v_l[n_Cell] * (1.0 - alpha_outlet_ghost) * rho_l_outlet_bc;
  rho_g_flux[n_Cell] = (v_g[n_Cell] > 0) ? v_g[n_Cell] * alpha_e[n_Cell - 1] * rho_g_func(p_e[n_Cell - 1]) : v_g[n_Cell] * alpha_outlet_ghost * rho_g_outlet_bc;

  for(int i = 1; i < n_Cell; i++)
  {
    rho_l_flux[i] = (v_l[i] > 0) ? v_l[i] * (1.0 - alpha_e[i-1]) * rho_l_func(p_e[i-1]) : v_l[i] * (1.0 - alpha_w[i]) * rho_l_func(p_w[i]);
    rho_g_flux[i] = (v_g[i] > 0) ? v_g[i] * alpha_e[i-1] * rho_g_func(p_e[i-1]) : v_g[i] * alpha_w[i] * rho_g_func(p_w[i]);
  }

  // momentum
  for(int i = 0; i < n_Cell + 1; i++)
  {
    double v_l_west = 0, v_l_east = 0;
    if(v_l[i] > 0)
    {
      v_l_west = (i == 0) ? v_l_e[n_Cell-1] : v_l_e[i-1];
      v_l_east = v_l_e[i];
    }
    else
    {
      v_l_east = (i == n_Cell) ? v_l_w[1] : v_l_w[i+1];
      v_l_west = v_l_w[i];
    }
    double dv_l_dx = (v_l_east - v_l_west) / dx;

    double v_g_west = 0, v_g_east = 0;
    if(v_g[i] > 0)
    {
      v_g_west = (i == 0) ? v_g_e[n_Cell-1] : v_g_e[i-1];
      v_g_east = v_g_e[i];
    }
    else
    {
      v_g_east = (i == n_Cell) ? v_g_w[1] : v_g_w[i+1];
      v_g_west = v_g_w[i];
    }
    double dv_g_dx = (v_g_east - v_g_west) / dx;

    // dp_dx term
    double dp_dx = 0;
    if (i == 0)               dp_dx = (p[0] - p_inlet_bc) / dx;
    else if (i == n_Cell)     dp_dx = (p_outlet_bc - p[n_Cell-1]) / dx;
    else                      dp_dx = (p[i] - p[i-1])/dx;

    // gravity term and interfacial friction (drag)
    double gravity_l = 0, gravity_g = 0;
    double fric_l = 0, fric_g = 0;
    if (_problem_type == MANOMETER)
    {
      // gravity
      if (i == 0)
      {
        gravity_l = (1-alpha[0]) * rho_l[0] * gx_vol(0);
        gravity_g = alpha[0] * rho_g[0] * gx_vol(0);
      }
      else if (i == n_Cell)
      {
        gravity_l = (1-alpha[n_Cell-1]) * rho_l[n_Cell-1] * gx_vol(n_Cell-1);
        gravity_g = alpha[n_Cell-1] * rho_g[n_Cell-1] * gx_vol(n_Cell-1);
      }
      else
      {
        gravity_l = 0.5 * ((1-alpha[i-1]) * rho_l[i-1] * gx_vol(i-1) + (1-alpha[i]) * rho_l[i] * gx_vol(i));
        gravity_g = 0.5 * (alpha[i-1] * rho_g[i-1] * gx_vol(i-1) + alpha[i] * rho_g[i] * gx_vol(i));
      }

      double alpha_edge_l = std::max(1 - alpha_edge[i], ALPHA_MIN);
      gravity_l = gravity_l / alpha_edge_l / rho_l_edge[i];
      double alpha_edge_g = std::max(alpha_edge[i], ALPHA_MIN);
      gravity_g = gravity_g / alpha_edge_g / rho_g_edge[i];

      // drag
      double rp = 5e-4;
      double cd = 0.44;
      double a_int = 3 * std::max(alpha_edge[i] * (1 - alpha_edge[i]), ALPHA_MIN * (1 - ALPHA_MIN)) / rp;
      double rho_mix = alpha_edge[i] * rho_g_edge[i] + (1 - alpha_edge[i]) * rho_l_edge[i];
      double v_diff_sqr = (v_g[i] - v_l[i]) * std::fabs(v_g[i] - v_l[i]);

      fric_l =  0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_l / rho_l_edge[i];
      fric_g = -0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_g / rho_g_edge[i];
      /*
      double rp = 5.e-4;
      double cd = 0.44; //2.2; //4.4;
      double a_int = 3 * alpha_edge[i] * (1 - alpha_edge[i]) / rp;
      double a_int_min = 3 * ALPHA_MIN * (1 - ALPHA_MIN) / rp;
      double rho_mix = alpha_edge[i] * rho_g_edge[i] + (1 - alpha_edge[i]) * rho_l_edge[i];
      double v_diff_sqr = (v_g[i] - v_l[i]) * std::fabs(v_g[i] - v_l[i]);

      if((1 - alpha_edge[i]) > ALPHA_MIN)
        fric_l = 0.125 * cd * a_int * rho_mix * v_diff_sqr / (1 - alpha_edge[i]) / rho_l_edge[i];
      else
        fric_l = 0.125 * cd * a_int_min * rho_mix * v_diff_sqr / ALPHA_MIN / rho_l_edge[i];

      if(alpha_edge[i] > ALPHA_MIN)
        fric_g = -0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge[i] / rho_g_edge[i];
      else
        fric_g = -0.125 * cd * a_int_min * rho_mix * v_diff_sqr / ALPHA_MIN / rho_g_edge[i];*/
    }

    // residuals
    rhs[4*i]   = -v_l[i] * dv_l_dx - dp_dx / rho_l_edge[i] + gravity_l + fric_l;
    rhs[4*i+1] = -v_g[i] * dv_g_dx - dp_dx / rho_g_edge[i] + gravity_g + fric_g;
  }

  // RHS for mass conservation equations
  for(int i = 0; i < n_Cell; i++) //loop on cells
  {
    // mass conservation equations
    rhs[4*i+2] = -(rho_l_flux[i+1] - rho_l_flux[i]) / dx;
    rhs[4*i+3] = -(rho_g_flux[i+1] - rho_g_flux[i]) / dx;
  }
}

double
FourEqnOneP::gx_vol(unsigned i)
{
  if (_problem_type == MANOMETER)     return (i < n_Cell/2) ? 9.81 : -9.81;
  else                                return 0;
}

void
FourEqnOneP::onTimestepEnd()
{
  // save old solutions
  alpha_oo  = alpha_old;  p_oo  = p_old;     v_l_oo  = v_l_old;   v_g_oo  = v_g_old;   rho_l_oo  = rho_l_old;   rho_g_oo  = rho_g_old;
  alpha_old = alpha;      p_old = p;         v_l_old = v_l;       v_g_old = v_g;       rho_l_old = rho_l;       rho_g_old = rho_g;

  if (_problem_type == MANOMETER)
  {
    time.push_back(_problemSystem->getCurrentTime());
    v_l_bottom.push_back(v_l[n_Cell/2]);
    p_bottom.push_back(0.5 * (p[n_Cell/2-1]+p[n_Cell/2]));
  }
}

void
FourEqnOneP::onLastTimestepEnd()
{
  if (_problem_type == MANOMETER)
  {
    FILE * file;
    std::string file_name = "output/manometer_results.csv";
    file = fopen(file_name.c_str(), "w");

    fprintf(file, "%20s%20s%20s\n", "time,", "v_l_bottom,", "p_bottom");
    for (unsigned i = 0; i < time.size(); i++)
      fprintf(file, "%20.6e,%20.6e,%20.6e\n", time[i], v_l_bottom[i], p_bottom[i]);

    fclose(file);
  }
}

void
FourEqnOneP::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_l");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", rho_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", rho_g[i]);
  fprintf(file, "        </DataArray>\n");
  if (_problem_type == SINE_WAVE) // analytical solution
  {
    fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "alpha_ana");
    double t = _problemSystem->getCurrentTime();
    for (unsigned i = 0; i < n_Cell; i++)
    {
      double xx = (i + 0.5) * dx;
      double alpha_ana = 0.5 + 0.2 * std::sin(2 * PI * (xx - t) / length);
      fprintf(file, "          %20.6f\n", alpha_ana);
    }
    fprintf(file, "        </DataArray>\n");
  }
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
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
FourEqnOneP::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s%20s\n", "x", "alpha", "p", "rho_l", "rho_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, alpha[i], p[i], rho_l[i], rho_g[i]);

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s%20s\n", "x", "v_l", "v_g");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "%20.6e%20.6e%20.6e\n", i*dx, v_l[i], v_g[i]);
}

void
FourEqnOneP::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  int n_Var = 4;
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

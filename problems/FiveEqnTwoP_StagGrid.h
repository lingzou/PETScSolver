#ifndef FIVE_EQN_TWOP_STAGGRID_H
#define FIVE_EQN_TWOP_STAGGRID_H

#include <vector>
#include "PETScProblem.h"

class FiveEqnTwoP_StagGrid : public PETScProblem
{
public:
  FiveEqnTwoP_StagGrid(InputParameterList & pList);
  ~FiveEqnTwoP_StagGrid();

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void onTimestepEnd();
  virtual void writeVTKOutput(unsigned int step) final;
  virtual void writeTextOutput(unsigned int step) final;

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) final;

  void RHS_1st_order(double * rhs);
  void RHS_2nd_order(double * rhs);
protected:
  void updateFluxes();
  void updateFluxes2ndOrder();
  void linearReconstruction(double, double, std::vector<double> &, std::vector<double> &, std::vector<double> &);

  double rho_l_func(double p) { return 1000.0 + 1.0e-7 * (p - 1.0e5); }
  double rho_g_func(double p) { return    0.5 + 1.0e-6 * (p - 1.0e5); }

protected:
  double ALPHA_INIT, V_L_INIT, V_G_INIT, P_INIT, C_L, C_G;
  double P_OUTLET, ALPHA_INLET, ALPHA_OUTLET;

  unsigned int _order;
  double H_inv;
  double g;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // State variables
  // 1) Primary variables
  std::vector<double> alpha,     p_l,     p_g,     v_l,     v_g;
  std::vector<double> alpha_old, p_l_old, p_g_old, v_l_old, v_g_old;
  std::vector<double> alpha_oo,  p_l_oo,  p_g_oo,  v_l_oo,  v_g_oo;
  // 2) Dependent variables (rho appears in d/dt terms, needs appropriately initialized)
  std::vector<double> rho_l, rho_g, rho_l_old, rho_g_old, rho_l_oo, rho_g_oo;

  // Helper variables
  std::vector<double> alpha_edge, rho_l_edge, rho_g_edge;
  std::vector<double> v_l_cell, v_g_cell;
  std::vector<double> v_hat, p_hat, mu, p_hat_edge;

  // Fluxes
  std::vector<double> alpha_flux, rho_l_flux, rho_g_flux;

  // Second-order helper variables
  std::vector<double> alpha_w, alpha_e, p_l_w, p_l_e, p_g_w, p_g_e;
  std::vector<double> v_l_w, v_l_e, v_g_w, v_g_e;
};

#endif /*FIVE_EQN_TWOP_STAGGRID_H*/

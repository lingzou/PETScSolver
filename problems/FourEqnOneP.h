#ifndef FOUR_EQN_ONEP_H
#define FOUR_EQN_ONEP_H

#include <vector>
#include "PETScProblem.h"

class FourEqnOneP : public PETScProblem
{
public:
  FourEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  ~FourEqnOneP();

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE * file) override final;
  virtual void writeTextOutput(FILE * file) override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;

  void RHS_1st_order(double * rhs);
  void RHS_2nd_order(double * rhs);

  virtual void onLastTimestepEnd() override final;

  enum ProblemType
  {
    SINE_WAVE       = 1,      // A gentel advection problem (sine wave)
    SQUARE_WAVE     = 2,      // Advection problem with extremely small alpha_l and alpha_g values, 1.e-6
    MANOMETER       = 3,      // Manometer problem
    SEDIMENTATION   = 4,      // Sedimentation problem
    WATER_FAUCET    = 5
  };

protected:
  double rho_l_func(double p) { return 1000.0 + 1.0e-7 * (p - 1.0e5); }
  double rho_g_func(double p) { return    0.5 + 1.0e-5 * (p - 1.0e5); }
  double gx_vol(unsigned);

protected:
  int _problem_type;
  bool _periodic_bc;
  double P_INLET, P_OUTLET, ALPHA_INLET, ALPHA_OUTLET, V_L_INLET, V_G_INLET;

  unsigned int _order;
  double g;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;
  double ALPHA_MIN;

  // State variables
  // 1) Primary variables
  std::vector<double> alpha,     p,     v_l,     v_g;
  std::vector<double> alpha_old, p_old, v_l_old, v_g_old;
  std::vector<double> alpha_oo,  p_oo,  v_l_oo,  v_g_oo;
  // 2) Dependent variables (rho appears in d/dt terms, needs appropriately initialized)
  std::vector<double> rho_l, rho_g, rho_l_old, rho_g_old, rho_l_oo, rho_g_oo;

  // Helper variables
  std::vector<double> alpha_edge, rho_l_edge, rho_g_edge;
  std::vector<double> v_l_cell, v_g_cell;

  // Fluxes
  std::vector<double> rho_l_flux, rho_g_flux;

  // Second-order helper variables
  std::vector<double> alpha_w, alpha_e, p_w, p_e;
  std::vector<double> v_l_w, v_l_e, v_g_w, v_g_e;

  // Post-processing stuff, e.g., in Manometer problem, v_l_bottom vs. time
  std::vector<double> time;
  std::vector<double> v_l_bottom, p_bottom;
};

#endif /*FOUR_EQN_ONEP_H*/

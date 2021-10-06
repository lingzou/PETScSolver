#ifndef SINGLE_PHASE_FLOW_H
#define SINGLE_PHASE_FLOW_H

#include <vector>
#include "PETScProblem.h"

class SinglePhaseFlow : public PETScProblem
{
public:
  SinglePhaseFlow(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  ~SinglePhaseFlow();

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void onTimestepEnd();
  virtual void writeVTKOutput(unsigned int step) final;
  virtual void writeTextOutput(unsigned int step) final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) final;

protected:
  void updateFluxes();
  void updateFluxes2ndOrder();

  double rho_func(double p, double T) { return 1.e3 + 4.e-7 * (p - 1.e5) - 0.46 * (T - 300.0); }
  double e_func(double /*p*/, double T) { return 112.55e3 + 4.e3 * (T - 300.0); }

protected:
  double P_INIT, V_INIT, T_INIT;
  double P_OUTLET, V_INLET, T_INLET, T_OUTLET;

  unsigned int _order;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // State variables
  // 1) Primary variables
  std::vector<double> p,        v,        T;
  std::vector<double> p_old,    v_old,    T_old;
  std::vector<double> p_oo,     v_oo,     T_oo;
  // 2) Dependent variables
  std::vector<double> rho, rho_old, rho_oo, e, e_old, e_oo;

  // Helper variables
  std::vector<double> rho_edge;

  // Fluxes
  std::vector<double> mass_flux, energy_flux;

  // Second-order helper variables
  std::vector<double> p_w, p_e, T_w, T_e, rho_w, rho_e, e_w, e_e;
};

#endif /*SINGLE_PHASE_FLOW_H*/

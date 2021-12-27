#pragma once

#include <vector>
#include "PETScProblem.h"

class EulerEquation1D : public PETScProblem
{
public:
  EulerEquation1D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~EulerEquation1D();

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE *) override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;
protected:
  void updateFluxes();
  void updateFluxes2ndOrder();

  double p_IG(double rho, double m, double E) { return (_gamma - 1.0) * (E - 0.5 * m * m / rho); }

protected:
  unsigned int _order;
  double _gamma;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // Initial conditions
  double RHO_L, M_L, E_L, P_L;
  double RHO_R, M_R, E_R, P_R;

  std::vector<double> rho, m, E;              // rho; m=rho*u; E=rho*e+0.5*u*u
  std::vector<double> rho_old, m_old, E_old;  // old solutions
  std::vector<double> rho_oo, m_oo, E_oo;     // old-old solutions
  std::vector<double> p;

  std::vector<double> rho_w, rho_e, m_w, m_e, E_w, E_e;
  std::vector<double> flux_rho, flux_m, flux_E; // Fluxes
};

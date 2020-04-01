#ifndef EULER_EQUATION_1D_H
#define EULER_EQUATION_1D_H

#include <vector>
#include "PETScProblem.h"

class EulerEquation1D : public PETScProblem
{
public:
  EulerEquation1D();
  ~EulerEquation1D();

  virtual void onTimestepEnd() final;

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u, TimeStepIndex index) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void writeSolution(unsigned int step) final;

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) final;
protected:
  void updateFluxes();
  void updateFluxes2ndOrder();
  void linearReconstruction(double, double, std::vector<double> &, std::vector<double> &, std::vector<double> &);

protected:
  unsigned int _order;
  double _gamma;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  std::vector<double> rho, m, E;              // rho; m=rho*u; E=rho*e+0.5*u*u
  std::vector<double> rho_old, m_old, E_old;  // old solutions
  std::vector<double> rho_oo, m_oo, E_oo;     // old-old solutions
  std::vector<double> p;

  std::vector<double> rho_w, rho_e, m_w, m_e, E_w, E_e;
  std::vector<double> flux_rho, flux_m, flux_E; // Fluxes
};

#endif /*EULER_EQUATION_1D_H*/
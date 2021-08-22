#ifndef HEAT_CONDUCTION_1D_H
#define HEAT_CONDUCTION_1D_H

#include <vector>
#include "PETScProblem.h"

static const double PI = asin(1) * 2;

class HeatConduction1D : public PETScProblem
{
public:
  HeatConduction1D(InputParameterList & pList);
  virtual ~HeatConduction1D();

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u, TimeStepIndex index) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void writeVTKOutput(unsigned int step) final;

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) final;
  virtual void computeJacobianMatrix(Mat & P_Mat) final;

protected:
  double length;
  double dx, dx2;
  unsigned int n_Cell;
  unsigned int n_Node;

  std::vector<double> x;          // x position
  std::vector<double> T;          // Temperature
  std::vector<double> T_old;      // Old Temperature
  std::vector<double> T_oldold;   // OldOld Temperature

  double k;
  double alpha;
};

#endif /*HEAT_CONDUCTION_1D_H*/

#ifndef HEAT_CONDUCTION_1D_H
#define HEAT_CONDUCTION_1D_H

#include <vector>
#include "PETScProblem.h"

class HeatConduction1D : public PETScProblem
{
public:
  HeatConduction1D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~HeatConduction1D();

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void writeVTKOutput(FILE *) override final;
  virtual void onTimestepEnd() override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;
  virtual void computeJacobianMatrix(Mat & P_Mat) override final;

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

#ifndef HEAT_CONDUCTION_1D_H
#define HEAT_CONDUCTION_1D_H

#include <petsc.h>
#include <vector>

enum TimeStepIndex
{
  NEW = 0,
  OLD = 1,
  OLDOLD = 2
};

static const double PI = asin(1) * 2;

class HeatConduction1D
{
public:
  HeatConduction1D();
  ~HeatConduction1D();

  unsigned int getNDOF() { return n_DOFs; }

  void SetupInitialCondition(double * u);
  void updateSolution(double *u, TimeStepIndex index);

  void transientResidual(double * res);
  void RHS(double * rhs);
  void writeSolution();

  void FillJacobianMatrix(Mat P_mat);

protected:
  double _dt;

  double length;
  double dx, dx2;
  unsigned int n_Cell;
  unsigned int n_Node;
  unsigned int n_DOFs;

  std::vector<double> x;          // x position
  std::vector<double> T;          // Temperature
  std::vector<double> T_old;      // Old Temperature
  std::vector<double> T_oldold;   // OldOld Temperature

  double k;
  double alpha;
};

#endif /*HEAT_CONDUCTION_1D_H*/

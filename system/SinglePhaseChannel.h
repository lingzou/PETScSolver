#ifndef SINGLE_PHASE_CHANNEL_H
#define SINGLE_PHASE_CHANNEL_H

#include <vector>
#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class SinglePhaseChannel : public PETScProblem
{
public:
  SinglePhaseChannel(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  ~SinglePhaseChannel();

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void onTimestepEnd();
  virtual void writeVTKOutput(unsigned int step) final;
  virtual void writeTextOutput(unsigned int step) final;

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) final;

protected:
  double P_INIT, V_INIT, T_INIT;
  double P_OUTLET, V_INLET, T_INLET, T_OUTLET;

  unsigned int _order;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // base data structures
  std::vector<SPCell*> _cells;
  std::vector<EdgeBase*> _edges;
};

#endif /*SINGLE_PHASE_CHANNEL_H*/

#ifndef PBC_H
#define PBC_H

#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class pBC : public PETScProblem
{
public:
  pBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~pBC();

  virtual void setDOFoffset(unsigned int offset) override;
  virtual void setupConnections() final;
  virtual void setupExtendedConnections() override { _edge->setExtendedNghbrs(); _edge->printConnection(); }

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u) final;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void onTimestepEnd() override;
  virtual void writeVTKOutput(unsigned int) override {}
  virtual void writeTextOutput(unsigned int) override {}

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) final;
  virtual void updateEdgeCellHelperVar() override;

protected:
  double P_OUTLET, T_OUTLET, V_INIT;

  unsigned int _order;
  double _dx;

  EdgeBase* _edge;
};

#endif /*PBC_H*/

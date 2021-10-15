#ifndef PBC_H
#define PBC_H

#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class pBC : public PETScProblem
{
public:
  pBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~pBC();

  virtual void setDOFoffset(unsigned int offset) override final;
  virtual void setupConnections() override final;
  virtual void setupExtendedConnections() override final { _edge->setExtendedNghbrs(); _edge->printConnection(); }

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE *) override final { /* do nothing */ }
  virtual void writeTextOutput(FILE *) override final { /* do nothing */ }

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;
  virtual void updateEdgeCellHelperVar() override final;

protected:
  unsigned int _order;
  double _dx;

  EdgeBase* _edge;
};

#endif /*PBC_H*/

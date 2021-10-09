#ifndef VBC_H
#define VBC_H

#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class vBC : public PETScProblem
{
public:
  vBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~vBC();

  virtual void setupConnections() override final;
  virtual void setupExtendedConnections() override final { _edge->setExtendedNghbrs(); _edge->printConnection(); }

  virtual void setDOFoffset(unsigned offset) override final;
  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(unsigned int) override final {}
  virtual void writeTextOutput(unsigned int) override final {}

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;
  virtual void updateEdgeCellHelperVar() override final;

protected:
  double V_INLET, T_INLET;

  unsigned _order;
  double _dx;

  EdgeBase* _edge;
  SPCell* _EAST_cell;
};

#endif /*VBC_H*/

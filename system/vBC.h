#ifndef VBC_H
#define VBC_H

#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class vBC : public PETScProblem
{
public:
  vBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~vBC();

  virtual void setupConnections() final;
  virtual void setupExtendedConnections() override { _edge->setExtendedNghbrs(); _edge->printConnection(); }

  virtual void setDOFoffset(unsigned int offset) override;
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
  double V_INLET, T_INLET;

  unsigned _order;
  double _dx;

  EdgeBase* _edge;
  SPCell* _EAST_cell;
};

#endif /*VBC_H*/

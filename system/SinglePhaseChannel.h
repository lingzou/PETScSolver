#ifndef SINGLE_PHASE_CHANNEL_H
#define SINGLE_PHASE_CHANNEL_H

#include <vector>
#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class SinglePhaseChannel : public PETScProblem
{
public:
  SinglePhaseChannel(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~SinglePhaseChannel();

  virtual void SetupInitialCondition(double * u) final;
  virtual void updateSolution(double *u) final;
  virtual void setDOFoffset(unsigned int offset) override;

  virtual void transientResidual(double * res) final;
  virtual void RHS(double * rhs) final;
  virtual void onTimestepEnd() override;
  virtual void writeVTKOutput(unsigned int step) override;
  virtual void writeTextOutput(unsigned int step) override;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) final;

  // connection
  virtual double getDx() { return dx; }
  virtual SPCell* getFirstCell() { return _cells.front(); }
  virtual SPCell* getLastCell()  { return _cells.back();  }
  virtual void acceptConnections(EdgeBase*, std::string);
  virtual void setupExtendedConnections() override;
  virtual void linearReconstruction() override;
  virtual void updateEdgeCellHelperVar() override;
protected:
  double P_INIT, V_INIT, T_INIT;

  unsigned int _order;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // base data structures
  std::vector<SPCell*> _cells;
  std::vector<EdgeBase*> _edges;
  EdgeBase *edge_begin, *edge_end;
};

#endif /*SINGLE_PHASE_CHANNEL_H*/

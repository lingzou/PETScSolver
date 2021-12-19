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

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;
  virtual void setDOFoffset(unsigned offset) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE *) override final;
  virtual void writeTextOutput(FILE *) override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;

  // connection
  virtual double getDx() const { return dx; }
  virtual SPCell* getFirstCell() const { return _cells.front(); }
  virtual SPCell* getLastCell()  const { return _cells.back();  }
  virtual void acceptConnections(EdgeBase*, std::string);
  virtual void setupExtendedConnections() override final;
  virtual void linearReconstruction() override final;
  virtual void updateEdgeCellHelperVar() override final;
protected:
  double P_INIT, V_INIT, T_INIT;

  unsigned int _order;

  double length, dx;
  double f, dh, gx;
  double h, aw, Tw;
  unsigned int n_Cell, n_Node;

  // base data structures
  std::vector<SPCell*> _cells;
  std::vector<EdgeBase*> _edges;
  EdgeBase *edge_begin, *edge_end;
};

#endif /*SINGLE_PHASE_CHANNEL_H*/

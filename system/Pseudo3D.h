#pragma once

#include <vector>
#include "PETScProblem.h"
#include "SinglePhaseFields.h"

class Channel;
class Pseudo3D;

class CrossFlow
{
public:
  CrossFlow(Channel* west_chan, Channel* east_chan, InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~CrossFlow();
  virtual void setDOFoffset(unsigned int offset);
  virtual unsigned int getDOFoffset() const final { return _DOF_offset; }
  virtual unsigned int getNDOF() const final { return _n_DOFs; }

  virtual void onTimestepEnd();
  virtual void SetupInitialCondition(double * u);
  virtual void updateSolution(double * u);

  virtual void transientResidual(double * res);
  virtual void RHS(double * rhs);

  virtual void writeVTKOutput(FILE *) {}
  virtual void writeTextOutput(FILE *) {}

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp);
  //virtual void computeJacobianMatrix(Mat & P_Mat);

  //virtual void setupExtendedConnections();
  //virtual void linearReconstruction();
  virtual void updateEdgeCellHelperVar();

  // temp functions
  virtual void setParent(Pseudo3D* parent) { _parent = parent; }
  virtual void printSolution() { for(auto& edge : _edges) std::cout << "v_cross = " << edge->v() << std::endl; }

protected:
  InputParameterList & _globalParamList;
  InputParameterList & _inputParamList;
  ProblemSystem * _problemSystem;
  Pseudo3D * _parent;

  const TimeScheme& _time_scheme;
  const double& _dt;
  const double& _dt_old;

  unsigned _order;
  double length;
  double dx;
  unsigned n_Cell, n_Node;

  unsigned int _DOF_offset;
  unsigned int _n_DOFs;

  // base data structures
  std::vector<crossEdge*> _edges;
  Channel *WEST_CHAN, *EAST_CHAN;
};

class Channel
{
public:
  Channel(unsigned i, unsigned j, InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~Channel();
  virtual void setDOFoffset(unsigned int offset);
  virtual unsigned int getDOFoffset() const final { return _DOF_offset; }
  virtual unsigned int getNDOF() const final { return _n_DOFs; }

  virtual void onTimestepEnd();
  virtual void SetupInitialCondition(double * u);
  virtual void updateSolution(double * u);

  virtual void transientResidual(double * res);
  virtual void RHS(double * rhs);

  virtual void writeVTKOutput(FILE *);
  virtual void writeTextOutput(FILE *) {}

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp);
  //virtual void computeJacobianMatrix(Mat & P_Mat);

  virtual void setupExtendedConnections();
  virtual void linearReconstruction();
  virtual void updateEdgeCellHelperVar(double dh);

  // temp functions
  virtual void setParent(Pseudo3D* parent) { _parent = parent; }
  virtual EdgeBase * getFirstEdge() { return _edges.front(); }
  virtual SPCell * getFirstCell() { return _cells.front(); }
  virtual SPCell * getCell(unsigned i) { return _cells[i]; }

protected:
  InputParameterList & _globalParamList;
  InputParameterList & _inputParamList;
  ProblemSystem * _problemSystem;
  Pseudo3D * _parent;

  double P_INIT, V_INIT, T_INIT;

  const TimeScheme& _time_scheme;
  const double& _dt;
  const double& _dt_old;

  unsigned _order;
  double length;
  double dx;
  unsigned n_Cell, n_Node;

  unsigned int _DOF_offset;
  unsigned int _n_DOFs;

  // base data structures
  std::vector<SPCell*> _cells;
  std::vector<EdgeBase*> _edges;

  unsigned _row, _col;
};

class Pseudo3D : public PETScProblem
{
public:
  Pseudo3D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~Pseudo3D();

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
  virtual void setupExtendedConnections() override final;
  virtual void linearReconstruction() override final;
  virtual void updateEdgeCellHelperVar() override final;

  // temp functions
  double getPInlet() { return _p_inlet; }
protected:
  double P_INIT, V_INIT, T_INIT;

  double _p_inlet, _p_inlet_o, _p_inlet_oo;

  unsigned int _order;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // base data structures
  std::vector<Channel*> _channels;
  std::vector<CrossFlow*> _crossflows;
  std::vector<double> _Areas, _dh;
  double _PI;
};

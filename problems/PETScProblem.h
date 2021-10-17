#ifndef PETSC_PROBLEM_H
#define PETSC_PROBLEM_H

#include <petsc.h>
#include "ParameterList.h"
#include "ProblemSystem.h"

class PETScProblem
{
public:
  PETScProblem(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  virtual ~PETScProblem();

  virtual void setDOFoffset(unsigned int offset) { _DOF_offset = offset; }
  virtual unsigned int getDOFoffset() const final { return _DOF_offset; }
  virtual unsigned int getNDOF() const final { return _n_DOFs; }

  virtual void onTimestepEnd() = 0;
  virtual void SetupInitialCondition(double * u) = 0;
  virtual void updateSolution(double *u) = 0;

  virtual void transientResidual(double * res) = 0;
  virtual void RHS(double * rhs) = 0;

  virtual void writeVTKOutput(FILE *) = 0;
  virtual void writeTextOutput(FILE *);

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) = 0;
  virtual void computeJacobianMatrix(Mat & P_Mat);

  virtual void setupConnections() {}
  virtual void setupExtendedConnections() {}
  virtual void linearReconstruction() {}
  virtual void updateEdgeCellHelperVar() {}

protected:
  InputParameterList & _globalParamList;
  InputParameterList & _inputParamList;
  ProblemSystem * _problemSystem;

  TimeScheme _time_scheme;
  double _dt;

  unsigned _n_DOFs, _DOF_offset;
};
#endif /*PETSC_PROBLEM_H*/

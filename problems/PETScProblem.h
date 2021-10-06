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

  virtual void setDOFoffset(long int offset) final { _DOF_offset = offset; }
  virtual unsigned int getDOFoffset() const final { return _DOF_offset; }
  virtual unsigned int getNDOF() const final { return _n_DOFs; }

  virtual void onTimestepEnd() = 0;
  virtual void SetupInitialCondition(double * u) = 0;
  virtual void updateSolution(double *u) = 0;

  virtual void transientResidual(double * res) = 0;
  virtual void RHS(double * rhs) = 0;

  virtual void writeVTKOutput(unsigned int step) = 0;
  virtual void writeTextOutput(unsigned int step);

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) = 0;
  virtual void computeJacobianMatrix(Mat & P_Mat);

protected:
  InputParameterList & _globalParamList;
  InputParameterList & _inputParamList;
  ProblemSystem * _problemSystem;

  std::string _input_file_name;
  TimeScheme _time_scheme;
  double _dt;

  unsigned int _DOF_offset;
  unsigned int _n_DOFs;
};
#endif /*PETSC_PROBLEM_H*/

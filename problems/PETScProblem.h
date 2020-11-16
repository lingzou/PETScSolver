#ifndef PETSC_PROBLEM_H
#define PETSC_PROBLEM_H

#include <petsc.h>
#include "PETScProblemInterface.h"

enum TimeStepIndex
{
  NEW = 0,
  OLD = 1,
  OLDOLD = 2
};

enum TimeScheme
{
  BDF1    = 0,
  BDF2    = 1,
  CN      = 2,
  INVALID = 99
};

TimeScheme StringToEnum(std::string str);

std::string trim_file_name(std::string file_name_with_path);

class PETScProblem
{
public:
  PETScProblem();
  virtual ~PETScProblem();

  virtual unsigned int getNDOF() { return n_DOFs; }
  virtual TimeScheme getTimeScheme() { return _time_scheme; }
  virtual double getCurrentTime() { return _t; }
  virtual double getDt() { return _dt; }

  virtual void onTimestepEnd() = 0;

  virtual void SetupInitialCondition(double * u) = 0;
  virtual void updateSolution(double *u, TimeStepIndex index) = 0;

  virtual void transientResidual(double * res) = 0;
  virtual void RHS(double * rhs) = 0;
  virtual void writeSolution(unsigned int step) = 0;

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) = 0;
  virtual void computeJacobianMatrix(Mat & P_Mat);

protected:
  std::string _input_file_name;
  TimeScheme _time_scheme;
  double _t;
  double _dt;
  unsigned int _step;

  unsigned int n_DOFs;
};
#endif /*PETSC_PROBLEM_H*/

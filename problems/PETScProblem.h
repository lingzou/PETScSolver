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

class PETScProblem
{
public:
  PETScProblem(InputParameterList & pList);
  virtual ~PETScProblem();

  virtual unsigned int getNDOF() { return n_DOFs; }
  virtual TimeScheme getTimeScheme() { return _time_scheme; }
  virtual double getCurrentTime() { return _t; }
  virtual double getDt() { return _dt; }

  virtual void onTimestepEnd();

  virtual void SetupInitialCondition(double * u) = 0;
  virtual void updateSolution(double *u, TimeStepIndex index) = 0;

  virtual void transientResidual(double * res) = 0;
  virtual void RHS(double * rhs) = 0;

  void writeOutput(unsigned int step);
  virtual void writeVTKOutput(unsigned int step) = 0;
  virtual void writeTextOutput(unsigned int step);

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat) = 0;
  virtual void computeJacobianMatrix(Mat & P_Mat);

protected:
  InputParameterList & paramList;
  std::string _input_file_name;
  TimeScheme _time_scheme;
  double _t;
  double _dt;
  int _n_steps;
  unsigned int _step;
  int _output_interval;
  bool _text_output;

  unsigned int n_DOFs;
};
#endif /*PETSC_PROBLEM_H*/

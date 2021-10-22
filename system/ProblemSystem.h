#ifndef PROBLEM_SYSTEM_H
#define PROBLEM_SYSTEM_H

#include <petsc.h>
#include "ParameterList.h"
#include "PETScProblemInterface.h"
#include "SinglePhaseFluid.h"

class PETScProblem;
class MatrixNonZeroPattern;

class ProblemSystem
{
public:
  ProblemSystem(InputParser* input_parser);
  virtual ~ProblemSystem();

  virtual unsigned int getNDOF() { return _n_DOFs; }
  virtual TimeScheme getTimeScheme() { return _time_scheme; }
  virtual double getCurrentTime() { return _t; }
  virtual unsigned int getCurrentTimeStep() { return _step; }
  virtual double getDt() { return _dt; }

  virtual void onTimestepEnd();

  virtual void SetupInitialCondition(double * u);
  virtual void updateSolution(double *u);

  virtual void transientResidual(double * res);
  virtual void RHS(double * rhs);

  virtual void writeOutput(unsigned int step);

  virtual void FillJacobianMatrixNonZeroPattern(Mat & P_Mat);
  virtual void computeJacobianMatrix(Mat & P_Mat);

  virtual PETScProblem* getProblem(std::string prob_name);
  virtual SinglePhaseFluid* getDefaultFluid();

protected:
  std::map <std::string, PETScProblem*> problem_system;
  std::map <std::string, SinglePhaseFluid*> fluid_system;

  InputParameterList & _globalParamList;
  std::map<std::string, InputParameterList *>& _problemParamList_map;
  std::map<std::string, InputParameterList *>& _fluidParamList_map;

  std::string _input_file_name;
  TimeScheme _time_scheme;
  double _t;
  double _dt;
  int _n_steps;
  unsigned int _step;
  int _output_interval;
  bool _text_output;

  unsigned int _n_DOFs;
};
#endif /*PROBLEM_SYSTEM_H*/

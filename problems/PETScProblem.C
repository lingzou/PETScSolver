#include <iostream>
#include "PETScProblem.h"
#include "PETScProblemInterface.h"

PETScProblem::PETScProblem(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  _globalParamList(globalParamList),
  _inputParamList(inputParamList),
  _problemSystem(problemSystem),
  //_input_file_name(_globalParamList.getParameterValue<std::string>("input_file_name")),
  _time_scheme(_globalParamList.getParameterValue<TimeScheme>("ts")),
  _dt(_globalParamList.getParameterValue<double>("dt")),
  _n_DOFs(0)
{
}

PETScProblem::~PETScProblem()
{}

void
PETScProblem::computeJacobianMatrix(Mat & P_Mat)
{
  // By default this is not required
  sysError("Exact Jacobian has not been implemented.");
}

void
PETScProblem::writeTextOutput(FILE * /**/)
{
  sysError("writeTextOutput() function has not been implemented.");
}

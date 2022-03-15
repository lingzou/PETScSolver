#include <iostream>
#include "PETScProblem.h"
#include "PETScProblemInterface.h"

PETScProblem::PETScProblem(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  _globalParamList(globalParamList),
  _inputParamList(inputParamList),
  _problemSystem(problemSystem),
  _time_scheme(_problemSystem->getTimeScheme()),
  _dt(_problemSystem->getDt()),
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

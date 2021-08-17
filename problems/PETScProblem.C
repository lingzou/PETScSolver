#include <iostream>
#include "PETScProblem.h"
#include "PETScProblemInterface.h"

PETScProblem::PETScProblem(ParameterList & pList) :
  paramList(pList),
  _input_file_name(paramList.getParameterValue<std::string>("input_file_name")),
  _time_scheme(paramList.getParameterValue<TimeScheme>("ts")),
  _t(0.0),
  _dt(paramList.getParameterValue<double>("dt")),
  _step(1),
  _n_steps(paramList.getParameterValue<int>("n_steps")),
  _output_interval(paramList.getParameterValue<int>("output_interval")),
  _text_output(paramList.getParameterValue<bool>("text_output")),
  n_DOFs(1)
{
}

PETScProblem::~PETScProblem()
{}

void
PETScProblem::onTimestepEnd()
{
  // March time forward
  _t += _dt;

  // write solution
  int N_Steps = PetscOptionsGetRequiredInt("-n_steps");
  if ((_step % _output_interval == 0) || (_step == N_Steps))
    writeOutput(_step);

  _step ++;
}

void
PETScProblem::computeJacobianMatrix(Mat & P_Mat)
{
  // By default this is not required
  sysError("Exact Jacobian has not been implemented.");
}

void
PETScProblem::writeOutput(unsigned int step)
{
  writeVTKOutput(step);
  if (_text_output) writeTextOutput(step);
}

void
PETScProblem::writeTextOutput(unsigned int step)
{
  sysError("writeTextOutput() function has not been implemented.");
}

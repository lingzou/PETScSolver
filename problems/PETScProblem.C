#include <iostream>
#include "PETScProblem.h"
#include "PETScProblemInterface.h"

TimeScheme
StringToEnum(std::string str)
{
  if      (str.compare("BDF1") == 0)  return BDF1;
  else if (str.compare("BDF2") == 0)  return BDF2;
  else if (str.compare("CN")   == 0)  return CN;
  else    {sysError("ERROR: UNKNOWN TimeScheme: " + str); return INVALID;}
}

std::string
trim_file_name(std::string full_file_name)
{
  unsigned int pos_of_point = full_file_name.find_last_of(".");
  unsigned int pos_of_slash = full_file_name.find_last_of("\\/");
  unsigned int length = pos_of_point - pos_of_slash - 1;
  if(length < 1)
  {
    sysError("ERROR: The input file name, '" + full_file_name + "', cannot be properly trimmed.");
    return std::string("");
  }
  else
    return full_file_name.substr(pos_of_slash + 1, length);
}

PETScProblem::PETScProblem() :
  _time_scheme(INVALID), _t(0.0), _dt(0.0), _step(1), n_DOFs(1)
{
  std::string full_name = PetscOptionsGetRequiredString("-input_file_name");
  _input_file_name = trim_file_name(full_name);

  _dt = PetscOptionsGetRequiredReal("-dt");
  std::string ts_str = PetscOptionsGetRequiredString("-ts");
  _time_scheme = StringToEnum(ts_str);

  _output_interval = PetscOptionsGetOptionalInt("-output_interval", 1);
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
    writeVTKOutput(_step);
  _step ++;
}

void
PETScProblem::computeJacobianMatrix(Mat & P_Mat)
{
  // By default this is not required
  sysError("Exact Jacobian has not been implemented.");
}

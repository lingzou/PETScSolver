#include <iostream>
#include "PETScProblem.h"
#include "PETScProblemInterface.h"

TimeScheme
StringToEnum(std::string str)
{
  if      (str.compare("BDF1") == 0)  return BDF1;
  else if (str.compare("BDF2") == 0)  return BDF2;
  else if (str.compare("CN")   == 0)  return CN;
  else { std::cerr << "ERROR: UNKNOWN TimeScheme: " << str << std::endl; exit(1); return INVALID; }
}

PETScProblem::PETScProblem() :
  _time_scheme(INVALID), _t(0.0), _dt(0.0), _step(1), n_DOFs(1)
{
  _dt = PetscOptionsGetRequiredReal("-dt");
  std::string ts_str = PetscOptionsGetRequiredString("-ts");
  _time_scheme = StringToEnum(ts_str);
}

PETScProblem::~PETScProblem()
{}

void
PETScProblem::computeJacobianMatrix(Mat & P_Mat)
{
  // By default this is not required
  std::cerr << "Exact Jacobian has not been implemented." << std::endl;
  exit(1);
}

#include <iostream>
#include "PETScProblem.h"

TimeScheme
StringToEnum(const char* str)
{
  if      (strcmp("BDF1", str) == 0)  return BDF1;
  else if (strcmp("BDF2", str) == 0)  return BDF2;
  else if (strcmp("CN", str) == 0)    return CN;
  else { std::cerr << "ERROR: UNKNOWN TimeScheme: " << str << std::endl; exit(1); return INVALID; }
}

PETScProblem::PETScProblem(TimeScheme ts, double t_start, double dt) :
  _time_scheme(ts), _t(t_start), _dt(dt), _step(1)
{}

PETScProblem::~PETScProblem()
{}

void
PETScProblem::computeJacobianMatrix(Mat & P_Mat)
{
  // By default this is not required
  std::cerr << "Exact Jacobian has not been implemented." << std::endl;
}

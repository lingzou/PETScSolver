#include <iostream>
#include "PETScProblem.h"

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

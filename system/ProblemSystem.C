#include <iostream>
#include "ProblemSystem.h"
#include "PETScProblem.h"

#include "HeatConduction1D.h"
#include "EulerEquation1D.h"
#include "FiveEqnTwoP_StagGrid.h"
#include "SinglePhaseFlow.h"
#include "SinglePhaseChannel.h"

ProblemSystem::ProblemSystem(InputParameterList & globalParameterList, std::map<std::string, InputParameterList *>& problemParamList_map) :
  _globalParamList(globalParameterList),
  _problemParamList_map(problemParamList_map),
  _input_file_name(_globalParamList.getParameterValue<std::string>("input_file_name")),
  _time_scheme(_globalParamList.getParameterValue<TimeScheme>("ts")),
  _t(0.0),
  _dt(_globalParamList.getParameterValue<double>("dt")),
  _step(1),
  _n_steps(_globalParamList.getParameterValue<int>("n_steps")),
  _output_interval(_globalParamList.getParameterValue<int>("output_interval")),
  _text_output(_globalParamList.getParameterValue<bool>("text_output")),
  _n_DOFs(0)
{
  for (auto& it : _problemParamList_map)
  {
    std::string problem_name = it.second->getParameterValue<std::string>("type");

    if (problem_name.compare("HeatConduction1D") == 0)
      problem_system[it.first] = new HeatConduction1D(_globalParamList, *(it.second), this);
    else if (problem_name.compare("EulerEquation1D") == 0)
      problem_system[it.first] = new EulerEquation1D(_globalParamList, *(it.second), this);
    else if (problem_name.compare("FiveEqnTwoP_StagGrid") == 0)
      problem_system[it.first] = new FiveEqnTwoP_StagGrid(_globalParamList, *(it.second), this);
    else if (problem_name.compare("SinglePhaseFlow") == 0)
      problem_system[it.first] = new SinglePhaseFlow(_globalParamList, *(it.second), this);
    else if (problem_name.compare("SinglePhaseChannel") == 0)
      problem_system[it.first] = new SinglePhaseChannel(_globalParamList, *(it.second), this);
    else
      sysError("ERROR: UNKNOWN problem: " + problem_name);
  }

  // Survey each subProblems to compute total n_DOFs
  for (auto& it : problem_system)
  {
    it.second->setDOFoffset(_n_DOFs);
    _n_DOFs += it.second->getNDOF();
  }
}

ProblemSystem::~ProblemSystem()
{
  for (auto& it : problem_system)  delete it.second;
}

void
ProblemSystem::onTimestepEnd()
{
  // March time forward
  _t += _dt;

  // write solution
  if ((_step % _output_interval == 0) || (_step == _n_steps))     writeOutput(_step);

  // Let subProblems act
  for (auto& it : problem_system)   it.second->onTimestepEnd();

  // March step forward
  _step ++;
}

void
ProblemSystem::SetupInitialCondition(double * u)
{
  for (auto& it : problem_system)
  {
    unsigned int offset = it.second->getDOFoffset();
    double * u_local = &u[offset];
    it.second->SetupInitialCondition(u_local);
  }
}

void
ProblemSystem::updateSolution(double *u)
{
  for (auto& it : problem_system)
  {
    unsigned int offset = it.second->getDOFoffset();
    double * u_local = &u[offset];
    it.second->updateSolution(u_local);
  }
}

void
ProblemSystem::transientResidual(double * res)
{
  for (auto& it : problem_system)
  {
    unsigned int offset = it.second->getDOFoffset();
    double * res_local = &res[offset];
    it.second->transientResidual(res_local);
  }
}

void
ProblemSystem::RHS(double * rhs)
{
  for (auto& it : problem_system)
  {
    unsigned int offset = it.second->getDOFoffset();
    double * rhs_local = &rhs[offset];
    it.second->RHS(rhs_local);
  }
}

void
ProblemSystem::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  // FIXME
  for (auto& it : problem_system)
    it.second->FillJacobianMatrixNonZeroPattern(P_Mat);
}

void
ProblemSystem::computeJacobianMatrix(Mat & P_Mat)
{
  // FIXME
  for (auto& it : problem_system)
    it.second->computeJacobianMatrix(P_Mat);
}

void
ProblemSystem::writeOutput(unsigned int step)
{
  for (auto& it : problem_system)
  {
    it.second->writeVTKOutput(step);
    if (_text_output) it.second->writeTextOutput(step);
  }
}

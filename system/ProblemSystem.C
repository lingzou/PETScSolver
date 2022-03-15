#include <iostream>
#include "ProblemSystem.h"
#include "PETScProblem.h"

#include "HeatConduction1D.h"
#include "EulerEquation1D.h"
#include "FiveEqnTwoP_StagGrid.h"
#include "FourEqnOneP.h"
#include "SixEqnOneP.h"
#include "SinglePhaseFlow.h"
#include "SinglePhaseChannel.h"
#include "vBC.h"
#include "pBC.h"
#include "Pseudo3D.h"

ProblemSystem::ProblemSystem(InputParser* input_parser) :
  _globalParamList(input_parser->getGlobalParamList()),
  _problemParamList_map(input_parser->getProblemSystemParamList()),
  _fluidParamList_map(input_parser->getFluidParamList()),
  _input_file_name(_globalParamList.getParameterValue<std::string>("input_file_name")),
  _time_scheme(_globalParamList.getParameterValue<TimeScheme>("ts")),
  _t(0.0),
  _dt(_globalParamList.getParameterValue<double>("dt")),
  _dt_max(_dt),
  _step(1),
  _n_steps(_globalParamList.getParameterValue<int>("n_steps")),
  _output_interval(_globalParamList.getParameterValue<int>("output_interval")),
  _text_output(_globalParamList.getParameterValue<bool>("text_output")),
  _n_DOFs(0)
{
  for (auto& it : _fluidParamList_map)
  {
    std::string fluid_name = it.second->getParameterValue<std::string>("type");
    if (fluid_name == "linearFluid")
      fluid_system[it.first] = new linearFluid(*(it.second));
    else
      sysError("ERROR: UNKNOWN fluid: " + fluid_name);
  }

  for (auto& it : _problemParamList_map)
  {
    std::string problem_name = it.second->getParameterValue<std::string>("type");

    if (problem_name == "HeatConduction1D")
      problem_system[it.first] = new HeatConduction1D(_globalParamList, *(it.second), this);
    else if (problem_name == "EulerEquation1D")
      problem_system[it.first] = new EulerEquation1D(_globalParamList, *(it.second), this);
    else if (problem_name == "FiveEqnTwoP_StagGrid")
      problem_system[it.first] = new FiveEqnTwoP_StagGrid(_globalParamList, *(it.second), this);
    else if (problem_name == "FourEqnOneP")
      problem_system[it.first] = new FourEqnOneP(_globalParamList, *(it.second), this);
    else if (problem_name == "SixEqnOneP")
      problem_system[it.first] = new SixEqnOneP(_globalParamList, *(it.second), this);
    else if (problem_name == "SinglePhaseFlow")
      problem_system[it.first] = new SinglePhaseFlow(_globalParamList, *(it.second), this);
    else if (problem_name == "SinglePhaseChannel")
      problem_system[it.first] = new SinglePhaseChannel(_globalParamList, *(it.second), this);
    else if (problem_name == "vBC")
      problem_system[it.first] = new vBC(_globalParamList, *(it.second), this);
    else if (problem_name == "pBC")
      problem_system[it.first] = new pBC(_globalParamList, *(it.second), this);
    else if (problem_name == "Pseudo3D")
      problem_system[it.first] = new Pseudo3D(_globalParamList, *(it.second), this);
    else
      sysError("ERROR: UNKNOWN problem: " + problem_name);
  }

  // Survey each subProblems to compute total n_DOFs
  for (auto& it : problem_system)
  {
    it.second->setDOFoffset(_n_DOFs);
    _n_DOFs += it.second->getNDOF();
  }

  for (auto& it : problem_system) it.second->setupConnections();
  for (auto& it : problem_system) it.second->setupExtendedConnections();
}

ProblemSystem::~ProblemSystem()
{
  for (auto& it : problem_system)   delete it.second;
  for (auto& it : fluid_system)     delete it.second;
}

void
ProblemSystem::adjustTimeStepSize(double ratio)
{
  _dt = std::min(_dt * ratio, _dt_max);
  for (auto& it : problem_system)   it.second->updateTimeStepSize(_dt);
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

  if (_step == _n_steps)
    for (auto& it : problem_system)   it.second->onLastTimestepEnd();

  // March step forward
  _step ++;
}

void
ProblemSystem::SetupInitialCondition(double * u)
{
  for (auto& it : problem_system)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->SetupInitialCondition(u + offset);
  }
}

void
ProblemSystem::updateSolution(double *u)
{
  for (auto& it : problem_system)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->updateSolution(u + offset);
  }

  for (auto& it : problem_system)   it.second->linearReconstruction();
  for (auto& it : problem_system)   it.second->updateEdgeCellHelperVar();
}

void
ProblemSystem::transientResidual(double * res)
{
  for (auto& it : problem_system)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->transientResidual(res + offset);
  }
}

void
ProblemSystem::RHS(double * rhs)
{
  for (auto& it : problem_system)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->RHS(rhs + offset);
  }
}

void
ProblemSystem::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatrixNonZeroPattern * mnzp = new MatrixNonZeroPattern(_n_DOFs);

  for (auto& it : problem_system)     it.second->FillJacobianMatrixNonZeroPattern(mnzp);
  std::vector<std::set<unsigned int>> &nzp = mnzp->getNonZeroPattern();

  PetscInt *nnz;
  PetscMalloc(_n_DOFs * sizeof(PetscInt), &nnz);
  for(unsigned i = 0; i < nzp.size(); i++)    nnz[i] = nzp[i].size();

  MatCreateSeqAIJ(PETSC_COMM_SELF, _n_DOFs, _n_DOFs, 0, nnz, &P_Mat);

  PetscReal zero = 0.0;
  for(unsigned int i = 0; i < nzp.size(); i++)  // loop on rows
  {
    PetscInt row = i;
    for (auto j : nzp[i]) // for each row, loop on its columns
    {
      PetscInt col = j;
      MatSetValues(P_Mat, 1, &row, 1, &col, &zero, INSERT_VALUES);
    }
  }
  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);

  // MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  // MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD);

  delete mnzp;
  PetscFree(nnz);
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
  // vtk (vtu) output
  FILE * vtk_File;
  std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".vtu";
  vtk_File = fopen(file_name.c_str(), "w");

  fprintf(vtk_File, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(vtk_File, "  <UnstructuredGrid>\n");
  for (auto& it : problem_system)     it.second->writeVTKOutput(vtk_File);
  fprintf(vtk_File, "  </UnstructuredGrid>\n");
  fprintf(vtk_File, "</VTKFile>\n");
  fclose(vtk_File);

  // txt output
  if (_text_output)
  {
    FILE * txt_File;
    std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".dat";
    txt_File = fopen(file_name.c_str(), "w");
    for (auto& it : problem_system)   it.second->writeTextOutput(txt_File);
    fclose(txt_File);
  }
}

PETScProblem*
ProblemSystem::getProblem(std::string name)
{
  if (problem_system.find(name) == problem_system.end())    sysError("Problem " + name + " does not exist.");
  return problem_system[name];
}

SinglePhaseFluid*
ProblemSystem::getDefaultFluid()
{
  if (fluid_system.empty())     sysError("There is no fluid properties in the system.");
  if (fluid_system.size() > 1)  sysError("There is more than one fluid properties in the system.");
  return fluid_system.begin()->second;
}

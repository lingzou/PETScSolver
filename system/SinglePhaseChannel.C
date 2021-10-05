#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "SinglePhaseChannel.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and S.J. Kim, Numerical study on the Welander oscillatory natural circulation
//       problem using high-order numerical methods, Progress in Nuclear Energy, 94 (2017) 162-172
//
// The momentum equation is simplified to use the primitive form. In the end, we solve:
// 1) Mass Equation:      d(rho)/dt + d(rho*v)/dx = 0
// 2) Momentum Equation:  rho dv/dt + rho*v dv/dx + dp/dx + f/(2dh) * rho*v*|v| = 0
// 3) Energy Equation:    d(rho*e)/dt + d(rho*v*e)/dx + p dv/dx - h_w*a_w*(T_w-T) = 0

/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0(v) 1(p) 3(v) 4(p)                                         3n+1(v)
        2(T)      5(T)
*/

// This problem has hard-coded boundary conditions: inlet velocity and T, outlet P.

SinglePhaseChannel::SinglePhaseChannel(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _inputParamList.readRequiredInputParameter<int>("order");
  _inputParamList.readRequiredInputParameter<int>("n_cells");
  _inputParamList.readRequiredInputParameter<double>("length");

  _inputParamList.readRequiredInputParameter<double>("P_INIT");
  _inputParamList.readRequiredInputParameter<double>("V_INIT");
  _inputParamList.readRequiredInputParameter<double>("T_INIT");

  _inputParamList.readRequiredInputParameter<double>("V_INLET");
  _inputParamList.readRequiredInputParameter<double>("T_INLET");
  _inputParamList.readRequiredInputParameter<double>("P_OUTLET");
  _inputParamList.readRequiredInputParameter<double>("T_OUTLET");

  // initial conditions
  P_INIT     =  _inputParamList.getParameterValue<double>("P_INIT");
  V_INIT     =  _inputParamList.getParameterValue<double>("V_INIT");
  T_INIT     =  _inputParamList.getParameterValue<double>("T_INIT");
  // boundary conditions
  V_INLET    =  _inputParamList.getParameterValue<double>("V_INLET");
  T_INLET    =  _inputParamList.getParameterValue<double>("T_INLET");
  P_OUTLET   =  _inputParamList.getParameterValue<double>("P_OUTLET");
  T_OUTLET   =  _inputParamList.getParameterValue<double>("T_OUTLET");  /* safeguard for reverse flow */

  _order = _inputParamList.getParameterValue<int>("order");
  n_Cell = _inputParamList.getParameterValue<int>("n_cells");
  length = _inputParamList.getParameterValue<double>("length");


  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 3 + 1;

  dx = length / n_Cell;

  // Create cells
  for (unsigned int i = 0; i < n_Cell; i++)
    _cells.push_back(new SPCell("CELL_"+std::to_string(i)));

  // Create edges
  _edges.push_back(new vBndryEdge("Inlet", V_INLET, T_INLET));
  for (unsigned int i = 1; i < n_Cell; i++)
    _edges.push_back(new IntEdge("EDGE_"+std::to_string(i)));
  _edges.push_back(new pBndryEdge("Outlet", P_OUTLET, T_OUTLET));

  // Setup neighboring edges for cells
  for (unsigned int i = 0; i < n_Cell; i++)
    _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);

  // Setup neighboring cells/edges for edges (has to be done after cells setup edges)
  _edges[0]->setNghbrCells(NULL, _cells[0]);
  _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  for (unsigned int i = 1; i < n_Cell; i++)
    _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);

  // debug:
  // for(auto& itr : _cells)   itr->printConnection();
  // for(auto& itr : _edges)   itr->printConnection();
}

SinglePhaseChannel::~SinglePhaseChannel()
{
  for (unsigned int i = 0; i < n_Cell + 1; i++)   delete _edges[i];
  for (unsigned int i = 0; i < n_Cell; i++)       delete _cells[i];
}

void
SinglePhaseChannel::SetupInitialCondition(double * u)
{
  for(int i = 0; i < n_Cell; i++)
  {
    _cells[i]->initialize(P_INIT, T_INIT);
    u[3*i+1] = P_INIT;
    u[3*i+2] = T_INIT;
  }

  for(int i = 0; i < n_Cell + 1; i++)
  {
    _edges[i]->initialize(V_INIT);
    u[3*i] = V_INIT;
  }
}

void
SinglePhaseChannel::updateSolution(double * u)
{
  for(int i = 0; i < n_Cell; i++)
    _cells[i]->updateSolution(u[3*i+1], u[3*i+2]);

  for(int i = 0; i < n_Cell + 1; i++)
    _edges[i]->updateSolution(u[3*i]);
}

void
SinglePhaseChannel::transientResidual(double * res)
{
  unsigned int idx = 0;
  unsigned int time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(int i = 0; i < n_Cell; i++)
    {
      res[3*i+1] = _cells[i]->massTranResBDF2(_dt);
      res[3*i+2] = _cells[i]->energyTranResBDF2(_dt);
    }

    for(int i = 0; i < n_Cell + 1; i++)
      res[3*i] = _edges[i]->computeTranResBDF2(_dt);
  }
  else
  {
    for(int i = 0; i < n_Cell; i++)
    {
      res[3*i+1] = _cells[i]->massTranRes(_dt);
      res[3*i+2] = _cells[i]->energyTranRes(_dt);
    }

    for(int i = 0; i < n_Cell + 1; i++)
      res[3*i] = _edges[i]->computeTranRes(_dt);
  }
}

void
SinglePhaseChannel::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:
      for(auto& itr : _edges)   itr->computeFluxes();
      break;
    case 2:
    {
      for(int i = 0; i < n_Cell; i++)
      {
        double p_W = (i == 0)        ? 2 * _cells[0]->p() - _cells[1]->p()                  : _cells[i-1]->p();
        double p_E = (i == n_Cell-1) ? 2 * _cells[n_Cell-1]->p() - _cells[n_Cell - 2]->p()  : _cells[i+1]->p();
        double T_W = (i == 0)        ? 2 * _cells[0]->T() - _cells[1]->T()                  : _cells[i-1]->T();
        double T_E = (i == n_Cell-1) ? 2 * _cells[n_Cell-1]->T() - _cells[n_Cell - 2]->T()  : _cells[i+1]->T();
        _cells[i]->linearReconstruction(p_W, p_E, T_W, T_E);
      }
      for(auto& itr : _edges)   itr->computeFluxes2nd();
    }
      break;
    default :   sysError("Spatial order not implemented.");
  }

  for (unsigned int i = 0; i < n_Cell; i++)
  {
    rhs[3*i+1] = _cells[i]->computeMassRHS(dx);
    rhs[3*i+2] = _cells[i]->computeEnergyRHS(dx);
  }
  for (unsigned int i = 0; i < n_Cell + 1; i++)
    rhs[3*i] = _edges[i]->computeRHS(dx);
}

void
SinglePhaseChannel::onTimestepEnd()
{
  // save old solutions
  for(auto& itr : _cells)   itr->saveOldSlns();
  for(auto& itr : _edges)   itr->saveOldSlns();
}

void
SinglePhaseChannel::writeVTKOutput(unsigned int step)
{
  FILE * ptr_File;
  std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".vtk";
  ptr_File = fopen(file_name.c_str(), "w");

  fprintf(ptr_File, "# vtk DataFile Version 4.0\n");
  fprintf(ptr_File, "my data\n");
  fprintf(ptr_File, "ASCII\n");
  fprintf(ptr_File, "DATASET STRUCTURED_GRID\n");
  fprintf(ptr_File, "DIMENSIONS %u 1 1\n", n_Node);
  fprintf(ptr_File, "POINTS %u Float32\n", n_Node);
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f 0 0\n", i * dx);

  // point data
  fprintf(ptr_File, "POINT_DATA %u\n", n_Node);
  // node id
  fprintf(ptr_File, "SCALARS node_id Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE node_id\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%d\n", i);

  // v
  fprintf(ptr_File, "SCALARS v Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE v\n");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%f\n", _edges[i]->v());

  // cell data
  fprintf(ptr_File, "CELL_DATA %u\n", n_Cell);
  // p
  fprintf(ptr_File, "SCALARS p Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE p\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", _cells[i]->p());

  // T
  fprintf(ptr_File, "SCALARS T Float32 1\n");
  fprintf(ptr_File, "LOOKUP_TABLE T\n");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%f\n", _cells[i]->T());

  fclose(ptr_File);
}

void
SinglePhaseChannel::writeTextOutput(unsigned int step)
{
  FILE * ptr_File;
  std::string file_name = "output/" + _input_file_name + "_step_" + std::to_string(step) + ".dat";
  ptr_File = fopen(file_name.c_str(), "w");

  // cell data
  fprintf(ptr_File, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(ptr_File, "#Cell data\n");
  fprintf(ptr_File, "%20s%20s%20s%20s%20s\n", "x", "p", "T", "rho", "v_cell");
  for (unsigned int i = 0; i < n_Cell; i++)
    fprintf(ptr_File, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, _cells[i]->p(), _cells[i]->T(), _cells[i]->rho(), 0.5 * (_edges[i]->v() + _edges[i+1]->v()));

  // edge data
  fprintf(ptr_File, "#Edge data\n");
  fprintf(ptr_File, "%20s%20s\n", "x", "v");
  for (unsigned int i = 0; i < n_Node; i++)
    fprintf(ptr_File, "%20.6e%20.6e\n", i*dx, _edges[i]->v());

  fclose(ptr_File);
}

void
SinglePhaseChannel::FillJacobianMatrixNonZeroPattern(Mat & P_Mat)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, _n_DOFs, _n_DOFs, 15, NULL, &P_Mat);

  int n_Var = 3;
  PetscReal v = 1.0;
  for (int i = 0; i < n_Cell + 1; i++)
  {
    for (int var = 0; var < n_Var; var++)
    {
      PetscInt i_dof = i * n_Var + var;
      for (int j_dof = (i - 2) * n_Var; j_dof < (i + 3) * n_Var; j_dof++)
      {
        if ((i_dof >= 0) && (i_dof < _n_DOFs) && (j_dof >= 0) && (j_dof < _n_DOFs))
          MatSetValues(P_Mat, 1, &i_dof, 1, &j_dof, &v, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);
  /*
  MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD);
  */
}

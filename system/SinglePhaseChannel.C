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
   *---------|---------|---------|---------|---------|---------*
  (BC)  0(p) 2(v) 3(p) ...                                    (BC)
        1(T)      4(T)
*/

// This problem has hard-coded boundary conditions: inlet velocity and T, outlet P.

SinglePhaseChannel::SinglePhaseChannel(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _inputParamList.addOptionalParamFromInput<double>("gx", 0.0);
  _inputParamList.addOptionalParamFromInput<double>("h", 0.0);
  _inputParamList.addOptionalParamFromInput<double>("aw", 0.0);
  _inputParamList.addOptionalParamFromInput<double>("Tw", 0.0);

  // initial conditions
  P_INIT =  _inputParamList.getValueFromInput<double>("P_INIT");
  V_INIT =  _inputParamList.getValueFromInput<double>("V_INIT");
  T_INIT =  _inputParamList.getValueFromInput<double>("T_INIT");

  _order = _inputParamList.getValueFromInput<int>("order");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");
  length = _inputParamList.getValueFromInput<double>("length");
  f      = _inputParamList.getValueFromInput<double>("f");
  dh     = _inputParamList.getValueFromInput<double>("dh");

  gx     = _inputParamList.getParameterValue<double>("gx");
  h      = _inputParamList.getParameterValue<double>("h");
  aw     = _inputParamList.getParameterValue<double>("aw");
  Tw     = _inputParamList.getParameterValue<double>("Tw");

  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 3-1;

  dx = length / n_Cell;

  SinglePhaseFluid* fluid = _problemSystem->getDefaultFluid();
  // Create cells/edges
  for (unsigned i = 0; i < n_Cell; i++)       _cells.push_back(new SPCell("CELL_"+std::to_string(i), fluid));
  for (unsigned i = 0; i < n_Cell - 1; i++)   _edges.push_back(new IntEdge("EDGE_"+std::to_string(i), fluid));

  // Setup neighboring edges for cells
  _cells[0]->setNghbrEdges(NULL, _edges[0]);
  for (unsigned i = 1; i < n_Cell - 1; i++)   _cells[i]->setNghbrEdges(_edges[i-1], _edges[i]);
  _cells[n_Cell-1]->setNghbrEdges(_edges[n_Cell-2], NULL);

  // Setup neighboring cells for edges
  for (unsigned i = 0; i < n_Cell-1; i++)     _edges[i]->setNghbrCells(_cells[i], _cells[i+1]);
}

void
SinglePhaseChannel::setDOFoffset(unsigned offset)
{
  _DOF_offset = offset;
  for(int i = 0; i < n_Cell; i++)     _cells[i]->setDOF(_DOF_offset + 3*i, _DOF_offset + 3*i + 1);
  for(int i = 0; i < n_Cell-1; i++)   _edges[i]->setDOF(_DOF_offset + 3*i + 2);
}

SinglePhaseChannel::~SinglePhaseChannel()
{
  for(auto& itr : _edges)   delete itr;
  for(auto& itr : _cells)   delete itr;
}

void
SinglePhaseChannel::acceptConnections(EdgeBase* edge, std::string type)
{
  if (type == "begin")      { edge_begin = edge; _cells[0]->setNghbrEdges(edge_begin, _edges[0]); }
  else if (type == "end")   { edge_end   = edge; _cells[n_Cell-1]->setNghbrEdges(_edges[n_Cell-2], edge_end); }
  else                      sysError("Incorrect keyword.");
}

void
SinglePhaseChannel::setupExtendedConnections()
{
  for(auto& itr : _edges)   itr->setExtendedNghbrs();
  for(auto& itr : _cells)   itr->setExtendedNghbrs();

  // debug
  // for(auto& itr : _cells)   itr->printConnection();
  // for(auto& itr : _edges)   itr->printConnection();
}

void
SinglePhaseChannel::SetupInitialCondition(double * u)
{
  for(int i = 0; i < n_Cell; i++)
  {
    _cells[i]->initialize(P_INIT, T_INIT);
    u[3*i] = P_INIT;
    u[3*i+1] = T_INIT;
  }

  for(int i = 0; i < n_Cell-1; i++)
  {
    _edges[i]->initialize(V_INIT);
    u[3*i+2] = V_INIT;
  }
}

void
SinglePhaseChannel::updateSolution(double * u)
{
  for(int i = 0; i < n_Cell; i++)     _cells[i]->updateSolution(u[3*i], u[3*i+1]);
  for(int i = 0; i < n_Cell-1; i++)   _edges[i]->updateSolution(u[3*i+2]);
}

void
SinglePhaseChannel::transientResidual(double * res)
{
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(int i = 0; i < n_Cell; i++)
    {
      res[3*i] = _cells[i]->massTranResBDF2(_dt, _dt_old);
      res[3*i+1] = _cells[i]->energyTranResBDF2(_dt, _dt_old);
    }
    for(int i = 0; i < n_Cell-1; i++)
      res[3*i+2] = _edges[i]->computeTranResBDF2(_dt, _dt_old);
  }
  else
  {
    for(int i = 0; i < n_Cell; i++)
    {
      res[3*i] = _cells[i]->massTranRes(_dt);
      res[3*i+1] = _cells[i]->energyTranRes(_dt);
    }
    for(int i = 0; i < n_Cell-1; i++)
      res[3*i+2] = _edges[i]->computeTranRes(_dt);
  }
}

void
SinglePhaseChannel::linearReconstruction()
{
  if (_order == 2)
    for(int i = 0; i < n_Cell; i++)
    {
      double p_W = (i == 0)        ? 2 * _cells[0]->p() - _cells[1]->p()                  : _cells[i-1]->p();
      double p_E = (i == n_Cell-1) ? 2 * _cells[n_Cell-1]->p() - _cells[n_Cell - 2]->p()  : _cells[i+1]->p();
      double T_W = (i == 0)        ? 2 * _cells[0]->T() - _cells[1]->T()                  : _cells[i-1]->T();
      double T_E = (i == n_Cell-1) ? 2 * _cells[n_Cell-1]->T() - _cells[n_Cell - 2]->T()  : _cells[i+1]->T();
      _cells[i]->linearReconstruction(p_W, p_E, T_W, T_E);
    }
}
void
SinglePhaseChannel::updateEdgeCellHelperVar()
{
  for (auto & cell : _cells)    cell->computeDP(f, dh, gx);

  switch (_order)
  {
    case 1:     for(auto& itr : _edges)   itr->computeFluxes();       break;
    case 2:     for(auto& itr : _edges)   itr->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
SinglePhaseChannel::RHS(double * rhs)
{
  for (unsigned i = 0; i < n_Cell; i++)
  {
    rhs[3*i] = _cells[i]->computeMassRHS(dx);
    rhs[3*i+1] = _cells[i]->computeEnergyRHS(dx, h, aw, Tw);
  }
  for (unsigned i = 0; i < n_Cell-1; i++)
    rhs[3*i+2] = _edges[i]->computeRHS(dx);
}

void
SinglePhaseChannel::onTimestepEnd()
{
  // save old solutions
  for(auto& itr : _cells)   itr->saveOldSlns();
  for(auto& itr : _edges)   itr->saveOldSlns();
}

void
SinglePhaseChannel::writeVTKOutput(FILE * file)
{
  fprintf(file, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_Node, n_Cell);
  fprintf(file, "      <Points>\n");
  fprintf(file, "        <DataArray type = \"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %f 0 0\n", i * dx);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Points>\n");

  fprintf(file, "      <Cells>\n");
  fprintf(file, "        <DataArray type = \"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d %d\n", i, i+1);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type = \"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d\n", 2*(i+1));
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type = \"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %d\n", 3);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  fprintf(file, "      <CellData>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->p());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "T");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->T());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </CellData>\n");

  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n", "node_id");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %d\n", i);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v");
  fprintf(file, "          %f\n", edge_begin->v());
  for (unsigned i = 0; i < n_Cell-1; i++)
    fprintf(file, "          %f\n", _edges[i]->v());
  fprintf(file, "          %f\n", edge_end->v());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
SinglePhaseChannel::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s\n", "x", "p", "T", "rho");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, _cells[i]->p(), _cells[i]->T(), _cells[i]->rho());

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s\n", "x", "v");
  fprintf(file, "%20.6e%20.6e\n", 0.0, edge_begin->v());
  for (unsigned i = 0; i < n_Cell - 1; i++)
    fprintf(file, "%20.6e%20.6e\n", (i+1)*dx, _edges[i]->v());
  fprintf(file, "%20.6e%20.6e\n", length, edge_end->v());
}

void
SinglePhaseChannel::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  for(auto& cell : _cells)   { mnzp->addRow(cell->pDOF(), cell->getConnectedDOFs()); mnzp->addRow(cell->TDOF(), cell->getConnectedDOFs());}
  for(auto& edge : _edges)   { mnzp->addRow(edge->vDOF(), edge->getConnectedDOFs()); }
}

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>

#include "Pseudo3D.h"
#include "utils.h"

CrossFlow::CrossFlow(Channel* west_chan, Channel* east_chan, InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  WEST_CHAN(west_chan),
  EAST_CHAN(east_chan),
  _globalParamList(globalParamList),
  _inputParamList(inputParamList),
  _problemSystem(problemSystem)
{
  _order = _inputParamList.getParameterValue<int>("order");
  n_Cell = _inputParamList.getParameterValue<int>("n_cells");
  length = _inputParamList.getParameterValue<double>("length");

  dx = length / n_Cell;
  n_Node = n_Cell + 1;
  _dt = _globalParamList.getParameterValue<double>("dt");
  _time_scheme = _globalParamList.getParameterValue<TimeScheme>("ts");

  _n_DOFs = n_Cell;

  SinglePhaseFluid* fluid = _problemSystem->getDefaultFluid();

  for (unsigned i = 0; i < n_Cell; i++)
    _edges.push_back(new crossEdge("CROSS_EDGE_"+std::to_string(i), fluid));

  for (unsigned i = 0; i < n_Cell; i++)
  {
    _edges[i]->setNghbrCells(WEST_CHAN->getCell(i), EAST_CHAN->getCell(i));
    WEST_CHAN->getCell(i)->addCrossEdge(_edges[i], -1);
    EAST_CHAN->getCell(i)->addCrossEdge(_edges[i], 1);
  }

  _edges[0]->setNghbrEdges(NULL, _edges[1]);
  for(unsigned i = 1; i < n_Cell - 1; i++)  _edges[i]->setNghbrEdges(_edges[i-1], _edges[i+1]);
  _edges[n_Cell - 1]->setNghbrEdges(_edges[n_Cell-2], NULL);
}

CrossFlow::~CrossFlow()
{
  for(auto& edge : _edges)   delete edge;
}

void
CrossFlow::setDOFoffset(unsigned offset)
{
  _DOF_offset = offset;
  for(unsigned i = 0; i < n_Cell; i++)     _edges[i]->setDOF(_DOF_offset + i);
}

void
CrossFlow::SetupInitialCondition(double * u)
{
  for(unsigned i = 0; i < n_Cell; i++)
  {
    _edges[i]->initialize(0);
    u[i] = 0;
  }
}

void
CrossFlow::updateSolution(double * u)
{
  for(unsigned i = 0; i < n_Cell; i++)   _edges[i]->updateSolution(u[i]);
}

void
CrossFlow::transientResidual(double * res)
{
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(int i = 0; i < n_Cell; i++)
      res[i] = _edges[i]->computeTranResBDF2(_dt);
  }
  else
  {
    for(int i = 0; i < n_Cell; i++)
      res[i] = _edges[i]->computeTranRes(_dt);
  }
}

void
CrossFlow::updateEdgeCellHelperVar()
{
  for(auto& itr : _edges)   itr->computeFluxes();
}

void
CrossFlow::RHS(double * rhs)
{
  for (unsigned int i = 0; i < n_Cell; i++)
    rhs[i] = _edges[i]->computeRHS(0.738 * 0.0254);
}

void
CrossFlow::onTimestepEnd()
{
  // save old solutions
  for(auto& edge : _edges)   edge->saveOldSlns();
}

void
CrossFlow::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  for(auto& edge : _edges)   { mnzp->addRow(edge->vDOF(), edge->getConnectedDOFs()); }
}

Channel::Channel(unsigned row, unsigned col, InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  _globalParamList(globalParamList),
  _inputParamList(inputParamList),
  _problemSystem(problemSystem),
  _row(row),
  _col(col)
{
  // std::cout << "Creating Channel" << std::endl;
  P_INIT =  _inputParamList.getParameterValue<double>("P_INIT");
  V_INIT =  _inputParamList.getParameterValue<double>("V_INIT");
  T_INIT =  _inputParamList.getParameterValue<double>("T_INIT");
  _order = _inputParamList.getParameterValue<int>("order");
  n_Cell = _inputParamList.getParameterValue<int>("n_cells");
  length = _inputParamList.getParameterValue<double>("length");

  dx = length / n_Cell;
  n_Node = n_Cell + 1;
  _dt = _globalParamList.getParameterValue<double>("dt");
  _time_scheme = _globalParamList.getParameterValue<TimeScheme>("ts");

  _n_DOFs = n_Cell * 3 + 1;

  SinglePhaseFluid* fluid = _problemSystem->getDefaultFluid();
  // Create cells
  for (unsigned i = 0; i < n_Cell; i++)
    _cells.push_back(new SPCell("CELL_"+std::to_string(i), fluid));

  // Create edges
  _edges.push_back(new pPseudo3DBndryEdge("Inlet", fluid, P_INIT+1e3, T_INIT));
  for (unsigned int i = 1; i < n_Cell; i++)
    _edges.push_back(new IntEdge("EDGE_"+std::to_string(i), fluid));
  _edges.push_back(new pBndryEdge("Outlet", fluid, P_INIT, T_INIT));

  // Setup neighboring edges for cells
  for (unsigned i = 0; i < n_Cell; i++)
    _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);
  // Setup neighboring cells for edges
  _edges[0]->setNghbrCells(NULL, _cells[0]);
  _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  for (unsigned i = 1; i < n_Cell; i++)
    _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
}

void
Channel::setDOFoffset(unsigned offset)
{
  _DOF_offset = offset;
  // std::cout << "Channel::setDOFoffset. _DOF_offset = " << _DOF_offset << std::endl;
  for(int i = 0; i < n_Node; i++)     _edges[i]->setDOF(_DOF_offset + 3*i);
  for(int i = 0; i < n_Cell; i++)     _cells[i]->setDOF(_DOF_offset + 3*i + 1, _DOF_offset + 3*i + 2);
}

Channel::~Channel()
{
  for(auto& edge : _edges)   delete edge;
  for(auto& cell : _cells)   delete cell;
}

void
Channel::setupExtendedConnections()
{
  for(auto& edge : _edges)   edge->setExtendedNghbrs();
  for(auto& cell : _cells)   cell->setExtendedNghbrs();
}

void
Channel::SetupInitialCondition(double * u)
{
  for(unsigned i = 0; i < n_Cell + 1; i++)
  {
    _edges[i]->initialize(V_INIT);
    u[3*i] = V_INIT;
  }

  for(unsigned i = 0; i < n_Cell; i++)
  {
    _cells[i]->initialize(P_INIT, T_INIT);
    u[3*i+1] = P_INIT;
    u[3*i+2] = T_INIT;
  }
}

void
Channel::updateSolution(double * u)
{
  for(unsigned i = 0; i < n_Cell + 1; i++)   _edges[i]->updateSolution(u[3*i]);
  for(unsigned i = 0; i < n_Cell; i++)       _cells[i]->updateSolution(u[3*i+1], u[3*i+2]);
}

void
Channel::transientResidual(double * res)
{
  unsigned time_step = _problemSystem->getCurrentTimeStep();
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
Channel::linearReconstruction()
{
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
Channel::updateEdgeCellHelperVar(double dh)
{
  for (auto & cell : _cells)    cell->computeDP(0.01, dh, -9.8);

  switch (_order)
  {
    case 1:     for(auto& itr : _edges)   itr->computeFluxes();       break;
    case 2:     for(auto& itr : _edges)   itr->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
Channel::RHS(double * rhs)
{
  for (unsigned i = 0; i < n_Cell; i++)
  {
    rhs[3*i+1] = _cells[i]->computeMassRHS(dx);
    rhs[3*i+2] = _cells[i]->computeEnergyRHS(dx, 2000, 300, 350);
  }
  for (unsigned int i = 0; i < n_Cell + 1; i++)
    rhs[3*i] = _edges[i]->computeRHS(dx);
}

void
Channel::writeVTKOutput(FILE * file)
{
  fprintf(file, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_Node, n_Cell);
  fprintf(file, "      <Points>\n");
  fprintf(file, "        <DataArray type = \"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %f %f %f\n", _row * 0.01, _col * 0.01, i * dx);
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
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %f\n", _edges[i]->v());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
Channel::onTimestepEnd()
{
  // save old solutions
  for(auto& edge : _edges)   edge->saveOldSlns();
  for(auto& cell : _cells)   cell->saveOldSlns();
}

void
Channel::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  for(auto& cell : _cells)   { mnzp->addRow(cell->pDOF(), cell->getConnectedDOFs()); mnzp->addRow(cell->TDOF(), cell->getConnectedDOFs());}
  for(auto& edge : _edges)   { mnzp->addRow(edge->vDOF(), edge->getConnectedDOFs()); }
}


/*************************************
*************************************/
Pseudo3D::Pseudo3D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _inputParamList.addRequiredParamFromInput<int>("order");
  _inputParamList.addRequiredParamFromInput<int>("n_cells");
  _inputParamList.addRequiredParamFromInput<double>("length");

  _inputParamList.addRequiredParamFromInput<double>("P_INIT");
  _inputParamList.addRequiredParamFromInput<double>("V_INIT");
  _inputParamList.addRequiredParamFromInput<double>("T_INIT");

  // initial conditions
  P_INIT     =  _inputParamList.getParameterValue<double>("P_INIT");
  V_INIT     =  _inputParamList.getParameterValue<double>("V_INIT");
  T_INIT     =  _inputParamList.getParameterValue<double>("T_INIT");

  _order = _inputParamList.getParameterValue<int>("order");
  n_Cell = _inputParamList.getParameterValue<int>("n_cells");
  length = _inputParamList.getParameterValue<double>("length");

  n_Node = n_Cell + 1;

  dx = length / n_Cell;

  for (unsigned i = 0; i < 4; i++)
    for (unsigned j = 0; j < 4; j++)
      _channels.push_back(new Channel(i, j, _globalParamList, _inputParamList, _problemSystem));

  double unit_area = (0.25 * 0.0254 * 0.25 * 0.0254) * 0.25;
  _Areas.resize(16, 0);
  _Areas[4*0 + 0] = unit_area;
  _Areas[4*0 + 1] = unit_area * 2;
  _Areas[4*0 + 2] = unit_area * 2;
  _Areas[4*0 + 3] = unit_area;
  _Areas[4*1 + 0] = unit_area * 2;
  _Areas[4*1 + 1] = unit_area * 4;
  _Areas[4*1 + 2] = unit_area * 4;
  _Areas[4*1 + 3] = unit_area * 2;
  _Areas[4*2 + 0] = unit_area * 2;
  _Areas[4*2 + 1] = unit_area * 4;
  _Areas[4*2 + 2] = unit_area * 4;
  _Areas[4*2 + 3] = unit_area * 2;
  _Areas[4*3 + 0] = unit_area;
  _Areas[4*3 + 1] = unit_area * 2;
  _Areas[4*3 + 2] = unit_area * 2;
  _Areas[4*3 + 3] = unit_area;

  _dh.resize(16, 0.02);

  double A_total = 0;
  for(auto & area : _Areas)  A_total += area;
  //for(auto & dh : _dh)       std::cout << "dh = " << dh << std::endl;
  std::cout << "A_total = " << A_total << std::endl;

  // create cross flows
  for (unsigned i = 0; i < 4; i++)
    for (unsigned j = 0; j < 4; j++)
    {
      unsigned this_idx = 4*i + j;
      if (j < 3) // has right
      {
        _crossflows.push_back(new CrossFlow(_channels[this_idx], _channels[this_idx+1], _globalParamList, _inputParamList, _problemSystem));
        //std::cout << "Cross flow connecting: " << this_idx << " to " << this_idx+1 << std::endl;
      }

      if (i < 3)
      {
        _crossflows.push_back(new CrossFlow(_channels[this_idx], _channels[this_idx+4], _globalParamList, _inputParamList, _problemSystem));
        //std::cout << "Cross flow connecting: " << this_idx << " to " << this_idx+4 << std::endl;
      }
    }

  _n_DOFs = 1; // inlet p
  for(auto& ch : _channels)     _n_DOFs += ch->getNDOF();
  for(auto& cf : _crossflows)   _n_DOFs += cf->getNDOF();

  for(auto& ch : _channels)     ch->setParent(this);
  for(auto& cf : _crossflows)   cf->setParent(this);

  std::cout << "Total nDOF = " << _n_DOFs << std::endl;
}

void
Pseudo3D::setDOFoffset(unsigned offset)
{
  _DOF_offset = offset;

  unsigned local_offset = offset + 1; // skip the inlet P
  for(auto& ch : _channels)
  {
    //std::cout << "Pseudo3D::setDOFoffset. 1) local_offset = " << local_offset << std::endl;
    ch->setDOFoffset(local_offset);
    local_offset += ch->getNDOF();
    //std::cout << "Pseudo3D::setDOFoffset. 2) local_offset = " << local_offset << std::endl;
  }
  for(auto& cf : _crossflows)
  {
    //std::cout << "CrossFlow::setDOFoffset. 1) local_offset = " << local_offset << std::endl;
    cf->setDOFoffset(local_offset);
    local_offset += cf->getNDOF();
    //std::cout << "CrossFlow::setDOFoffset. 2) local_offset = " << local_offset << std::endl;
  }
}

Pseudo3D::~Pseudo3D()
{
  for(auto& ch : _channels)     delete ch;
  for(auto& cf : _crossflows)   delete cf;
}

void
Pseudo3D::setupExtendedConnections()
{
  for(auto& ch : _channels)   ch->setupExtendedConnections();
  // CrossFlow no need setupExtendedConnections
}

void
Pseudo3D::SetupInitialCondition(double * u)
{
  _p_inlet = P_INIT; u[0] = P_INIT;
  unsigned offset = 1;
  for(auto& ch : _channels)     { ch->SetupInitialCondition(u + offset); offset += ch->getNDOF(); }
  for(auto& cf : _crossflows)   { cf->SetupInitialCondition(u + offset); offset += cf->getNDOF(); }
}

void
Pseudo3D::updateSolution(double * u)
{
  _p_inlet = u[0];

  unsigned offset = 1;
  for(auto& ch : _channels)
  {
    pPseudo3DBndryEdge * inlet_edge = dynamic_cast<pPseudo3DBndryEdge*>(ch->getFirstEdge());
    inlet_edge->updatePBC(_p_inlet);
    ch->updateSolution(u + offset); offset += ch->getNDOF();
  }

  for(auto& cf : _crossflows)
  {
    cf->updateSolution(u + offset); offset += cf->getNDOF();
  }
}

void
Pseudo3D::transientResidual(double * res)
{
  res[0] = 0; // no transient
  unsigned offset = 1;
  for(auto& ch : _channels)     { ch->transientResidual(res + offset); offset += ch->getNDOF(); }
  for(auto& cf : _crossflows)   { cf->transientResidual(res + offset); offset += cf->getNDOF(); }
}

void
Pseudo3D::linearReconstruction()
{
  if (_order == 2)
    for(auto& ch : _channels)   ch->linearReconstruction();
}
void
Pseudo3D::updateEdgeCellHelperVar()
{
  for(unsigned i = 0; i < 16; i++)    _channels[i]->updateEdgeCellHelperVar(_dh[i]);

  for(auto& cf : _crossflows) cf->updateEdgeCellHelperVar();
}

void
Pseudo3D::RHS(double * rhs)
{
  rhs[0] = 2.5; //_p_inlet;
  unsigned offset = 1;
  for (unsigned i = 0; i < 16; i++)
  {
    _channels[i]->RHS(rhs + offset); offset += _channels[i]->getNDOF();

    EdgeBase * inlet_edge = _channels[i]->getFirstEdge();
    rhs[0] -= inlet_edge->mass_flux() * _Areas[i];
  }

  for(auto& cf : _crossflows) { cf->RHS(rhs + offset); offset += cf->getNDOF();}
}

void
Pseudo3D::onTimestepEnd()
{
  _p_inlet_oo = _p_inlet_o; _p_inlet_o = _p_inlet;
  // save old solutions
  double total_mfr = 0;
  for (unsigned i = 0; i < 16; i++)
  {
    _channels[i]->onTimestepEnd();

    EdgeBase * inlet_edge = _channels[i]->getFirstEdge();
    total_mfr += inlet_edge->mass_flux() * _Areas[i];
  }
  std::cout << "total_mfr = " << total_mfr << std::endl;

  for(auto& cf : _crossflows) cf->onTimestepEnd();
}

void
Pseudo3D::writeVTKOutput(FILE * file)
{
  for(auto& ch : _channels)   ch->writeVTKOutput(file);
}

void
Pseudo3D::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());

  for(auto& ch : _channels)   ch->writeTextOutput(file);
}

void
Pseudo3D::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  mnzp->addEntry(0, 0);
  for(auto& ch : _channels)
  {
    mnzp->addEntry(0, ch->getFirstEdge()->vDOF());
    mnzp->addEntry(0, ch->getFirstCell()->pDOF());
    mnzp->addEntry(0, ch->getFirstCell()->TDOF());
    mnzp->addEntry(ch->getFirstEdge()->vDOF(), 0);
    mnzp->addEntry(ch->getFirstCell()->pDOF(), 0);
    mnzp->addEntry(ch->getFirstCell()->TDOF(), 0);

    ch->FillJacobianMatrixNonZeroPattern(mnzp);
  }

  for(auto& cf : _crossflows)   cf->FillJacobianMatrixNonZeroPattern(mnzp);
}

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "vBC.h"
#include "utils.h"
#include "SinglePhaseChannel.h"

vBC::vBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  // boundary conditions
  V_BC = _inputParamList.getValueFromInput<double>("V_BC");
  T_BC = _inputParamList.getValueFromInput<double>("T_BC");
  _order = _inputParamList.getValueFromInput<int>("order");

  _n_DOFs = 1;

  // Create edge
  SinglePhaseFluid* fluid = _problemSystem->getDefaultFluid();
  _edge = new vBndryEdge("vBC", fluid, V_BC, T_BC);
}

vBC::~vBC() { delete _edge; }

void
vBC::setDOFoffset(unsigned int offset) { _DOF_offset = offset; _edge->setDOF(_DOF_offset); }

void
vBC::setupConnections()
{
  std::string connection = _inputParamList.getValueFromInput<std::string>("Connection");
  std::string prob_name = connection.substr(0, connection.find(':'));
  std::string type = connection.substr(connection.find(':')+1, connection.size());

  if (type == "begin")
  {
    SinglePhaseChannel* spc = dynamic_cast<SinglePhaseChannel*>(_problemSystem->getProblem(prob_name));
    _edge->setNghbrCells(NULL, spc->getFirstCell());
    _dx = spc->getDx();
    spc->acceptConnections(_edge, "begin");
  }
  else sysError("vBC can only connect the 'begin' of flow channel.");
}

void
vBC::SetupInitialCondition(double * u)
{
  _edge->initialize(V_BC);
  u[0] = V_BC;
}

void
vBC::updateSolution(double * u)
{
  _edge->updateSolution(u[0]);
}

void
vBC::transientResidual(double * res)
{
  unsigned int time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))      res[0] = _edge->computeTranResBDF2(_dt, _dt_old);
  else                                                res[0] = _edge->computeTranRes(_dt);
}

void
vBC::updateEdgeCellHelperVar()
{
  switch (_order)
  {
    case 1:     _edge->computeFluxes();       break;
    case 2:     _edge->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
vBC::RHS(double * rhs)
{
  rhs[0] = _edge->computeRHS(_dx);
}

void
vBC::onTimestepEnd()
{
  // save old solutions
  _edge->saveOldSlns();
}

void
vBC::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  mnzp->addRow(_edge->vDOF(), _edge->getConnectedDOFs());
}

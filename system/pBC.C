#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "pBC.h"
#include "utils.h"
#include "SinglePhaseChannel.h"

pBC::pBC(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  // boundary conditions
  double p_BC   =  _inputParamList.getValueFromInput<double>("P_BC");
  double T_BC   =  _inputParamList.getValueFromInput<double>("T_BC");
  _order = _inputParamList.getValueFromInput<int>("order");

  _n_DOFs = 1;

  // Create edge
  SinglePhaseFluid* fluid = _problemSystem->getDefaultFluid();
  _edge = new pBndryEdge("Outlet", fluid, p_BC, T_BC);
}

pBC::~pBC() { delete _edge; }

void
pBC::setDOFoffset(unsigned int offset) { _DOF_offset = offset; _edge->setDOF(_DOF_offset); }

void
pBC::setupConnections()
{
  std::string connection = _inputParamList.getValueFromInput<std::string>("Connection");
  std::string prob_name = connection.substr(0, connection.find(':'));
  std::string type = connection.substr(connection.find(':')+1, connection.size());
  SinglePhaseChannel* spc = dynamic_cast<SinglePhaseChannel*>(_problemSystem->getProblem(prob_name));
  _dx = spc->getDx();

  if (type == "begin")      _edge->setNghbrCells(NULL, spc->getFirstCell());
  else if (type == "end")   _edge->setNghbrCells(spc->getLastCell(), NULL);
  else                      sysError("pBC can only connect the 'end' of flow channel.");

  spc->acceptConnections(_edge, type);
}

void
pBC::SetupInitialCondition(double * u)
{
  double v_initial = _inputParamList.getValueFromInput<double>("V_INIT");
  _edge->initialize(v_initial);
  u[0] = v_initial;
}

void
pBC::updateSolution(double * u)
{
  _edge->updateSolution(u[0]);
}

void
pBC::transientResidual(double * res)
{
  unsigned int time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))      res[0] = _edge->computeTranResBDF2(_dt, _dt_old);
  else                                                res[0] = _edge->computeTranRes(_dt);
}

void
pBC::updateEdgeCellHelperVar()
{
  switch (_order)
  {
    case 1:     _edge->computeFluxes();       break;
    case 2:     _edge->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
pBC::RHS(double * rhs)
{
  rhs[0] = _edge->computeRHS(_dx);
}

void
pBC::onTimestepEnd()
{
  // save old solutions
  _edge->saveOldSlns();
}

void
pBC::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  mnzp->addRow(_edge->vDOF(), _edge->getConnectedDOFs());
}

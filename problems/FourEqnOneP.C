#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "FourEqnOneP.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and H. Zhang,
//        Solving phase appearance/disappearance two-phase flow problems with high resolution staggered grid and
//        fully implicit schemes by the Jacobian-free Newton–Krylov Method
//        Computers & Fluids, Vol. 129, 179-188


/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0    2    4    6                                            4n
   1    3    5    7                                            4n+1
  v_l  alpha
  v_g  p
*/


/**** Data Structure for 4E1P problems ****/
double rho_l_func(double p) { return 1000.0 + 1.0e-7 * (p - 1.0e5); }
double rho_g_func(double p) { return    0.5 + 1.0e-5 * (p - 1.0e5); }

EdgeBase4E1P*
Cell4E1P::getOtherSideEdge(EdgeBase4E1P* edge)
{
  return (edge == WEST_EDGE) ? EAST_EDGE : WEST_EDGE;
}

void
Cell4E1P::setExtendedNghbrs()
{
  if (WEST_EDGE != NULL) WEST_CELL = WEST_EDGE->getOtherSideCell(this);
  if (EAST_EDGE != NULL) EAST_CELL = EAST_EDGE->getOtherSideCell(this);
}

void
Cell4E1P::initialize(double alpha, double p)
{
  _p = p;           _p_o = p;           _p_oo = p;
  _alpha = alpha;   _alpha_o = alpha;   _alpha_oo = alpha;
  _rho_l = rho_l_func(p);
  _rho_l_o = _rho_l; _rho_l_oo = _rho_l;
  _rho_g = rho_g_func(p);
  _rho_g_o = _rho_g; _rho_g_oo = _rho_g;
}

void
Cell4E1P::linearReconstruction(double alpha_W, double alpha_E, double p_W, double p_E)
{
  UTILS::linearReconstruction(alpha_W, alpha_E, _alpha, _alpha_w, _alpha_e);
  UTILS::linearReconstruction(p_W, p_E, _p, _p_w, _p_e);

  _rho_l_w = rho_l_func(_p_w);
  _rho_l_e = rho_l_func(_p_e);
  _rho_g_w = rho_g_func(_p_w);
  _rho_g_e = rho_g_func(_p_e);
}

void
Cell4E1P::updateSolution(double alpha, double p)
{
  _alpha = alpha; _p = p;
  _rho_l = rho_l_func(p);
  _rho_g = rho_g_func(p);
}

double
Cell4E1P::liquidMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux_l() - WEST_EDGE->mass_flux_l()) / dx;
}

double
Cell4E1P::gasMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux_g() - WEST_EDGE->mass_flux_g()) / dx;
}

void
Cell4E1P::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout          << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << std::endl;
}

std::vector<unsigned>
Cell4E1P::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.push_back(_aDOF);  dofs.push_back(_pDOF);
  if (WEST_CELL != NULL) { dofs.push_back(WEST_CELL->aDOF());   dofs.push_back(WEST_CELL->pDOF());}
  if (WEST_EDGE != NULL) { dofs.push_back(WEST_EDGE->vlDOF());  dofs.push_back(WEST_EDGE->vgDOF());}
  if (EAST_CELL != NULL) { dofs.push_back(EAST_CELL->aDOF());   dofs.push_back(EAST_CELL->pDOF());}
  if (EAST_EDGE != NULL) { dofs.push_back(EAST_EDGE->vlDOF());  dofs.push_back(EAST_EDGE->vgDOF());}
  return dofs;
}

void
Cell4E1P::saveOldSlns()
{
  _p_oo = _p_o; _alpha_oo = _alpha_o; _rho_l_oo = _rho_l_o; _rho_g_oo = _rho_g_o;
  _p_o  = _p;   _alpha_o  = _alpha;   _rho_l_o  = _rho_l;   _rho_g_o  = _rho_g;
}

void
EdgeBase4E1P::setNghbrCells(Cell4E1P * west, Cell4E1P * east)
{
  WEST_CELL = west;
  EAST_CELL = east;
}

Cell4E1P*
EdgeBase4E1P::getOtherSideCell(Cell4E1P* cell)
{
  return (cell == WEST_CELL) ? EAST_CELL : WEST_CELL;
}

void
EdgeBase4E1P::setExtendedNghbrs()
{
  if (WEST_CELL != NULL) WEST_EDGE = WEST_CELL->getOtherSideEdge(this);
  if (EAST_CELL != NULL) EAST_EDGE = EAST_CELL->getOtherSideEdge(this);
}

void
EdgeBase4E1P::linearReconstruction(double vl_W, double vl_E, double vg_W, double vg_E)
{
  UTILS::linearReconstruction(vl_W, vl_E, _vl, _vl_w, _vl_e);
  UTILS::linearReconstruction(vg_W, vg_E, _vg, _vg_w, _vg_e);
}

void
EdgeBase4E1P::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout          << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << std::endl;
}

std::vector<unsigned>
EdgeBase4E1P::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.push_back(_vlDOF); dofs.push_back(_vgDOF);
  if (WEST_CELL != NULL) { dofs.push_back(WEST_CELL->aDOF());   dofs.push_back(WEST_CELL->pDOF());}
  if (WEST_EDGE != NULL) { dofs.push_back(WEST_EDGE->vlDOF());  dofs.push_back(WEST_EDGE->vgDOF());}
  if (EAST_CELL != NULL) { dofs.push_back(EAST_CELL->aDOF());   dofs.push_back(EAST_CELL->pDOF());}
  if (EAST_EDGE != NULL) { dofs.push_back(EAST_EDGE->vlDOF());  dofs.push_back(EAST_EDGE->vgDOF());}
  return dofs;
}

void
vBndryEdge4E1P::computeFluxes()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p() : WEST_CELL->p();
  double rho_l_ghost = rho_l_func(p_ghost);
  double rho_g_ghost = rho_g_func(p_ghost);
  _mass_flux_l = (1 - _alpha_bc) * rho_l_ghost * _vl;
  _mass_flux_g = _alpha_bc * rho_g_ghost * _vg;
}

void
vBndryEdge4E1P::computeFluxes2nd()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p_w() : WEST_CELL->p_e();
  double rho_l_ghost = rho_l_func(p_ghost);
  double rho_g_ghost = rho_g_func(p_ghost);
  _mass_flux_l = (1 - _alpha_bc) * rho_l_ghost * _vl;
  _mass_flux_g = _alpha_bc * rho_g_ghost * _vg;
}

void
pBndryEdge4E1P::computeFluxes()
{
  double alpha_west = (WEST_CELL == NULL) ? _alpha_bc : WEST_CELL->alpha();
  double alpha_east = (EAST_CELL == NULL) ? _alpha_bc : EAST_CELL->alpha();
  double rho_l_west = (WEST_CELL == NULL) ? rho_l_func(_p_bc) : WEST_CELL->rho_l();
  double rho_l_east = (EAST_CELL == NULL) ? rho_l_func(_p_bc) : EAST_CELL->rho_l();
  double rho_g_west = (WEST_CELL == NULL) ? rho_g_func(_p_bc) : WEST_CELL->rho_g();
  double rho_g_east = (EAST_CELL == NULL) ? rho_g_func(_p_bc) : EAST_CELL->rho_g();

  _mass_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west : _vl * (1 - alpha_east) * rho_l_east;
  _mass_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west : _vg * alpha_east * rho_g_east;
}

void
pBndryEdge4E1P::computeFluxes2nd()
{
  double alpha_west = (WEST_CELL == NULL) ? _alpha_bc : WEST_CELL->alpha_e();
  double alpha_east = (EAST_CELL == NULL) ? _alpha_bc : EAST_CELL->alpha_w();
  double rho_l_west = (WEST_CELL == NULL) ? rho_l_func(_p_bc) : WEST_CELL->rho_l_e();
  double rho_l_east = (EAST_CELL == NULL) ? rho_l_func(_p_bc) : EAST_CELL->rho_l_w();
  double rho_g_west = (WEST_CELL == NULL) ? rho_g_func(_p_bc) : WEST_CELL->rho_g_e();
  double rho_g_east = (EAST_CELL == NULL) ? rho_g_func(_p_bc) : EAST_CELL->rho_g_w();

  _mass_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west : _vl * (1 - alpha_east) * rho_l_east;
  _mass_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west : _vg * alpha_east * rho_g_east;
}

double
pBndryEdge4E1P::liquidVelTranRes(double dt)
{
  return (_vl - _vl_o) / dt;
}

double
pBndryEdge4E1P::gasVelTranRes(double dt)
{
  return (_vg - _vg_o) / dt;
}

double
pBndryEdge4E1P::liquidVelTranResBDF2(double dt)
{
  return (1.5 * _vl - 2 * _vl_o + 0.5 * _vl_oo) / dt;
}

double
pBndryEdge4E1P::gasVelTranResBDF2(double dt)
{
  return (1.5 * _vg - 2 * _vg_o + 0.5 * _vg_oo) / dt;
}

void
pWESTBndryEdge4E1P::computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g)
{
  double rho_l_edge = EAST_CELL->rho_l();
  double rho_g_edge = EAST_CELL->rho_g();
  double alpha_edge = EAST_CELL->alpha();
  double alpha_edge_l = std::max(1 - alpha_edge, 1e-6);
  double alpha_edge_g = std::max(alpha_edge, 1e-6);

  double dv_l_dx = (_vl < 0) ? (EAST_EDGE->v_l() - _vl) / dx : 0;
  double dp_dx = 2 * (EAST_CELL->p() - _p_bc) / dx;
  double gravity_l = (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l() * EAST_CELL->g();
  gravity_l = gravity_l / rho_l_edge / alpha_edge_l;

  double dv_g_dx = (_vg < 0) ? (EAST_EDGE->v_g() - _vg) / dx : 0;
  double gravity_g = EAST_CELL->alpha() * EAST_CELL->rho_g() * EAST_CELL->g();
  gravity_g = gravity_g / rho_g_edge / alpha_edge_g;

  // friction
  double rp = 5e-4;
  double cd = 0.44;
  double a_int = 3 * std::max(alpha_edge * (1 - alpha_edge), 1e-6 * (1 - 1e-6)) / rp;
  double rho_mix = alpha_edge * rho_g_edge + (1 - alpha_edge) * rho_l_edge;
  double v_diff_sqr = (_vg - _vl) * std::fabs(_vg - _vl);

  double fric_l =  0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_l / rho_l_edge;
  double fric_g = -0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_g / rho_g_edge;

  rhs_l = -_vl * dv_l_dx - dp_dx / rho_l_edge + gravity_l + fric_l;
  rhs_g = -_vg * dv_g_dx - dp_dx / rho_g_edge + gravity_g + fric_g;
}

void
pEASTBndryEdge4E1P::computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g)
{
  double rho_l_edge = WEST_CELL->rho_l();
  double rho_g_edge = WEST_CELL->rho_g();
  double alpha_edge = WEST_CELL->alpha();
  double alpha_edge_l = std::max(1 - alpha_edge, 1e-6);
  double alpha_edge_g = std::max(alpha_edge, 1e-6);

  double dv_l_dx = (_vl > 0) ? (_vl - WEST_EDGE->v_l()) / dx : 0;
  double dp_dx = 2 * (_p_bc - WEST_CELL->p()) / dx;
  double gravity_l = (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l() * WEST_CELL->g();
  gravity_l = gravity_l / rho_l_edge / alpha_edge_l;

  double dv_g_dx = (_vg > 0) ? (_vg - WEST_EDGE->v_g()) / dx : 0;
  double gravity_g = WEST_CELL->alpha() * WEST_CELL->rho_g() * WEST_CELL->g();
  gravity_g = gravity_g / rho_g_edge / alpha_edge_g;

  // friction
  double rp = 5e-4;
  double cd = 0.44;
  double a_int = 3 * std::max(alpha_edge * (1 - alpha_edge), 1e-6 * (1 - 1e-6)) / rp;
  double rho_mix = alpha_edge * rho_g_edge + (1 - alpha_edge) * rho_l_edge;
  double v_diff_sqr = (_vg - _vl) * std::fabs(_vg - _vl);

  double fric_l =  0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_l / rho_l_edge;
  double fric_g = -0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_g / rho_g_edge;

  rhs_l = -_vl * dv_l_dx - dp_dx / rho_l_edge + gravity_l + fric_l;
  rhs_g = -_vg * dv_g_dx - dp_dx / rho_g_edge + gravity_g + fric_g;
}

void
IntEdge4E1P::computeFluxes()
{
  _mass_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l()
                           : _vl * (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l();

  _mass_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha() * WEST_CELL->rho_g()
                           : _vg * EAST_CELL->alpha() * EAST_CELL->rho_g();
}

void
IntEdge4E1P::computeFluxes2nd()
{
  _mass_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha_e()) * WEST_CELL->rho_l_e()
                           : _vl * (1 - EAST_CELL->alpha_w()) * EAST_CELL->rho_l_w();

  _mass_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha_e() * WEST_CELL->rho_g_e()
                           : _vg * EAST_CELL->alpha_w() * EAST_CELL->rho_g_w();
}

double
IntEdge4E1P::liquidVelTranRes(double dt)
{
  return (_vl - _vl_o) / dt;
}

double
IntEdge4E1P::gasVelTranRes(double dt)
{
  return (_vg - _vg_o) / dt;
}

double
IntEdge4E1P::liquidVelTranResBDF2(double dt)
{
  return (1.5 * _vl - 2 * _vl_o + 0.5 * _vl_oo) / dt;
}

double
IntEdge4E1P::gasVelTranResBDF2(double dt)
{
  return (1.5 * _vg - 2 * _vg_o + 0.5 * _vg_oo) / dt;
}

void
IntEdge4E1P::computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g)
{
  double rho_l_edge = 0.5 * (WEST_CELL->rho_l() + EAST_CELL->rho_l());
  double rho_g_edge = 0.5 * (WEST_CELL->rho_g() + EAST_CELL->rho_g());
  double alpha_edge = 0.5 * (WEST_CELL->alpha() + EAST_CELL->alpha());
  double alpha_edge_l = std::max(1 - alpha_edge, 1e-6);
  double alpha_edge_g = std::max(alpha_edge, 1e-6);

  double dv_l_dx = 0, dv_g_dx = 0;
  if (order == 1)
  {
    dv_l_dx = (_vl > 0) ? (_vl - WEST_EDGE->v_l()) / dx : (EAST_EDGE->v_l() - _vl) / dx;
    dv_g_dx = (_vg > 0) ? (_vg - WEST_EDGE->v_g()) / dx : (EAST_EDGE->v_g() - _vg) / dx;
  }
  else if (order == 2)
  {
    dv_l_dx = (_vl > 0) ? (_vl_e - WEST_EDGE->vl_e()) / dx : (EAST_EDGE->vl_w() - _vl_w) / dx;
    dv_g_dx = (_vg > 0) ? (_vg_e - WEST_EDGE->vg_e()) / dx : (EAST_EDGE->vg_w() - _vg_w) / dx;
  }
  else
    sysError("order not supported in IntEdge4E1P::computeLiquidRHS");

  double dp_dx = (EAST_CELL->p() - WEST_CELL->p()) / dx;

  double gravity_l = (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l() * EAST_CELL->g()
                   + (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l() * WEST_CELL->g();
  gravity_l *= 0.5;
  gravity_l = gravity_l / rho_l_edge / alpha_edge_l;

  double gravity_g = EAST_CELL->alpha() * EAST_CELL->rho_g() * EAST_CELL->g()
                   + WEST_CELL->alpha() * WEST_CELL->rho_g() * WEST_CELL->g();
  gravity_g *= 0.5;
  gravity_g = gravity_g / rho_g_edge / alpha_edge_g;

  // friction
  double rp = 5e-4;
  double cd = 0.44;
  double a_int = 3 * std::max(alpha_edge * (1 - alpha_edge), 1e-6 * (1 - 1e-6)) / rp;
  double rho_mix = alpha_edge * rho_g_edge + (1 - alpha_edge) * rho_l_edge;
  double v_diff_sqr = (_vg - _vl) * std::fabs(_vg - _vl);

  double fric_l =  0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_l / rho_l_edge;
  double fric_g = -0.125 * cd * a_int * rho_mix * v_diff_sqr / alpha_edge_g / rho_g_edge;

  rhs_l = -_vl * dv_l_dx - dp_dx / rho_l_edge + gravity_l + fric_l;
  rhs_g = -_vg * dv_g_dx - dp_dx / rho_g_edge + gravity_g + fric_g;
}
/**** End of Data Structure for 4E1P problems ****/

FourEqnOneP::FourEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _problem_type = _inputParamList.getValueFromInput<int>("problem_type");
  _order = _inputParamList.getValueFromInput<int>("order");
  length = _inputParamList.getValueFromInput<int>("length");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");

  // The sine wave and square wave problems use periodic BC
  _periodic_bc = (_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE);

  n_Node = n_Cell + 1;
  dx = length / n_Cell;

  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
  {
    _n_DOFs = n_Cell * 4;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell4E1P("CELL_"+std::to_string(i)));
    for (unsigned i = 0; i < n_Cell; i++)   _edges.push_back(new IntEdge4E1P("EDGE_"+std::to_string(i)));

    for (unsigned i = 0; i < n_Cell-1; i++) _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);
    _cells[n_Cell-1]->setNghbrEdges(_edges[n_Cell-1], _edges[0]);

    _edges[0]->setNghbrCells(_cells[n_Cell-1], _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++) _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
  }
  else if (_problem_type == MANOMETER)
  {
    _n_DOFs = n_Cell * 4 + 2;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell4E1P("CELL_"+std::to_string(i)));

    _edges.push_back(new pWESTBndryEdge4E1P("INLET", 1e5, 1));
    for (unsigned i = 0; i < n_Cell-1; i++)   _edges.push_back(new IntEdge4E1P("EDGE_"+std::to_string(i)));
    _edges.push_back(new pEASTBndryEdge4E1P("OUTLET", 1e5, 1));

    for (unsigned i = 0; i < n_Cell; i++)   _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);

    _edges[0]->setNghbrCells(NULL, _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++)   _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
    _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  }

  ALPHA_MIN = 1e-6;
}

FourEqnOneP::~FourEqnOneP()
{
  for(auto& itr : _edges)   delete itr;
  for(auto& itr : _cells)   delete itr;
}

void
FourEqnOneP::setDOFoffset(unsigned offset)
{
  _DOF_offset = offset;
  for(unsigned i = 0; i < _edges.size(); i++)   _edges[i]->setDOF(_DOF_offset + 4*i,      _DOF_offset + 4*i + 1);
  for(unsigned i = 0; i < _cells.size(); i++)   _cells[i]->setDOF(_DOF_offset + 4*i + 2,  _DOF_offset + 4*i + 3);
}

void
FourEqnOneP::setupExtendedConnections()
{
  for(auto& itr : _edges)   itr->setExtendedNghbrs();
  for(auto& itr : _cells)   itr->setExtendedNghbrs();
  /*
  for(auto& itr : _cells)   itr->printConnection();
  for(auto& itr : _edges)   itr->printConnection();
  */
}

void
FourEqnOneP::SetupInitialCondition(double * u)
{
  unsigned index = 0;
  switch (_problem_type)
  {
    case SINE_WAVE:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(1, 1);
        u[4*i] = 1;
        u[4*i+1] = 1;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double val = 0.5 + 0.2 * std::sin(2 * PI * xx / length);
        _cells[i]->initialize(val, 1e5);
        u[4*i+2] = val;
        u[4*i+3] = 1e5;

        _cells[i]->set_g(0);
      }
    break;

    case SQUARE_WAVE:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(1, 1);
        u[4*i] = 1;
        u[4*i+1] = 1;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double val = ((xx < 0.2) || (xx > 0.4)) ? 1 : 0;
        _cells[i]->initialize(val, 1e5);
        u[4*i+2] = val;
        u[4*i+3] = 1e5;

        _cells[i]->set_g(0);
      }
    break;

    case MANOMETER:
    /*
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(-1, -1);
        u[4*i] = -1;
        u[4*i+1] = -1;
      }

      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double val = ((xx < 5) || (xx > 15)) ? 1 : 0;
        _cells[i]->initialize(val, 1e5);

        u[4*i+2] = val;
        u[4*i+3] = 1e5;

        if (i < n_Cell/2)   _cells[i]->set_g(9.81);
        else                _cells[i]->set_g(-9.81);
      }*/
    break;

    case SEDIMENTATION:
    break;

    case WATER_FAUCET:
    break;

    defaut: sysError("Unknown problem_type.");
  }
}

void
FourEqnOneP::updateSolution(double * u)
{
  for(unsigned i = 0; i < _edges.size(); i++)   _edges[i]->updateSolution(u[4*i], u[4*i+1]);    // vl, vg
  for(unsigned i = 0; i < _cells.size(); i++)   _cells[i]->updateSolution(u[4*i+2], u[4*i+3]);  // alpha, p
}

void
FourEqnOneP::transientResidual(double * res)
{
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(unsigned i = 0; i < _edges.size(); i++)
    {
      res[4*i] = _edges[i]->liquidVelTranResBDF2(_dt);
      res[4*i+1] = _edges[i]->gasVelTranResBDF2(_dt);
    }
    for(unsigned i = 0; i < _cells.size(); i++)
    {
      res[4*i+2] = _cells[i]->liquidMassTranResBDF2(_dt);
      res[4*i+3] = _cells[i]->gasMassTranResBDF2(_dt);
    }
  }
  else
  {
    for(unsigned i = 0; i < _edges.size(); i++)
    {
      res[4*i] = _edges[i]->liquidVelTranRes(_dt);
      res[4*i+1] = _edges[i]->gasVelTranRes(_dt);
    }
    for(unsigned i = 0; i < _cells.size(); i++)
    {
      res[4*i+2] = _cells[i]->liquidMassTranRes(_dt);
      res[4*i+3] = _cells[i]->gasMassTranRes(_dt);
    }
  }
}

void
FourEqnOneP::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:     RHS_1st_order(rhs);     break;
    case 2:     RHS_2nd_order(rhs);     break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
FourEqnOneP::linearReconstruction()
{
  if (_order == 2)
  {
    if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {
        double alpha_W  = (i == 0) ? _cells.back()->alpha() : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? _cells.back()->p()     : _cells[i-1]->p();
        double alpha_E  = (i == n_Cell-1) ? _cells.front()->alpha() : _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? _cells.front()->p()     : _cells[i+1]->p();
        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E);
      }
      for (unsigned i = 0; i < _edges.size(); i++)
      {
        double vl_W     = (i == 0) ? _edges.back()->v_l()  : _edges[i-1]->v_l();
        double vg_W     = (i == 0) ? _edges.back()->v_g()  : _edges[i-1]->v_g();
        double vl_E     = (i == n_Cell-1) ? _edges.front()->v_l() : _edges[i+1]->v_l();
        double vg_E     = (i == n_Cell-1) ? _edges.front()->v_g() : _edges[i+1]->v_g();
        _edges[i]->linearReconstruction(vl_W, vl_E, vg_W, vg_E);
      }
    }
    /*
    else if (_problem_type == MANOMETER)
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {

        double alpha_W  = (i == 0) ? 2*_cells[0]->alpha()-_cells[1]->alpha()  : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? 2*_cells[0]->p()-_cells[1]->p()          : _cells[i-1]->p();
        double alpha_E  = (i == n_Cell-1) ? 2*_cells[n_Cell-1]->alpha()-_cells[n_Cell-2]->alpha() : _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? 2*_cells[n_Cell-1]->p()-_cells[n_Cell-2]->p()         : _cells[i+1]->p();
        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E);
      }
      for (unsigned i = 0; i < _edges.size(); i++)
      {

        double vl_W     = (i == 0) ? 2*_edges[0]->v_l()-_edges[1]->v_l()  : _edges[i-1]->v_l();
        double vg_W     = (i == 0) ? 2*_edges[0]->v_g()-_edges[1]->v_g()  : _edges[i-1]->v_g();
        double vl_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_l()-_edges[n_Cell-1]->v_l() : _edges[i+1]->v_l();
        double vg_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_g()-_edges[n_Cell-1]->v_g() : _edges[i+1]->v_g();
        _edges[i]->linearReconstruction(vl_W, vl_E, vg_W, vg_E);
      }
    }*/
    else
      sysError("Not implemented");
  }
}

void
FourEqnOneP::updateEdgeCellHelperVar()
{
  switch (_order)
  {
    case 1:     for(auto& edge : _edges)   edge->computeFluxes();       break;
    case 2:     for(auto& edge : _edges)   edge->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
FourEqnOneP::RHS_1st_order(double * rhs)
{
  for(unsigned i = 0; i < _edges.size(); i++)
    _edges[i]->computeRHS(1, dx, rhs[4*i], rhs[4*i+1]);

  for(unsigned i = 0; i < _cells.size(); i++)
  {
    rhs[4*i+2] = _cells[i]->liquidMassRHS(dx);
    rhs[4*i+3] = _cells[i]->gasMassRHS(dx);
  }
}

void
FourEqnOneP::RHS_2nd_order(double * rhs)
{
  for(unsigned i = 0; i < _edges.size(); i++)
    _edges[i]->computeRHS(2, dx, rhs[4*i], rhs[4*i+1]);

  for(unsigned i = 0; i < _cells.size(); i++)
  {
    rhs[4*i+2] = _cells[i]->liquidMassRHS(dx);
    rhs[4*i+3] = _cells[i]->gasMassRHS(dx);
  }
}

double
FourEqnOneP::gx_vol(unsigned i)
{
  if (_problem_type == MANOMETER)     return (i < n_Cell/2) ? 9.81 : -9.81;
  else                                return 0;
}

void
FourEqnOneP::onTimestepEnd()
{
  // save old solutions
  for(auto& itr : _cells)   itr->saveOldSlns();
  for(auto& itr : _edges)   itr->saveOldSlns();
  /*
  if (_problem_type == MANOMETER)
  {
    time.push_back(_problemSystem->getCurrentTime());
    v_l_bottom.push_back(v_l[n_Cell/2]);
    p_bottom.push_back(0.5 * (p[n_Cell/2-1]+p[n_Cell/2]));
  }*/
}

void
FourEqnOneP::onLastTimestepEnd()
{
  /*
  if (_problem_type == MANOMETER)
  {
    FILE * file;
    std::string file_name = "output/manometer_results.csv";
    file = fopen(file_name.c_str(), "w");

    fprintf(file, "%20s%20s%20s\n", "time,", "v_l_bottom,", "p_bottom");
    for (unsigned i = 0; i < time.size(); i++)
      fprintf(file, "%20.6e,%20.6e,%20.6e\n", time[i], v_l_bottom[i], p_bottom[i]);

    fclose(file);
  }*/
}

void
FourEqnOneP::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "alpha");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->alpha());//alpha[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->p());//p[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_l");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->rho_l());//rho_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->rho_g());//rho_g[i]);
  fprintf(file, "        </DataArray>\n");
  if (_problem_type == SINE_WAVE) // analytical solution
  {
    fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "alpha_ana");
    double t = _problemSystem->getCurrentTime();
    for (unsigned i = 0; i < n_Cell; i++)
    {
      double xx = (i + 0.5) * dx;
      double alpha_ana = 0.5 + 0.2 * std::sin(2 * PI * (xx - t) / length);
      fprintf(file, "          %20.6f\n", alpha_ana);
    }
    fprintf(file, "        </DataArray>\n");
  }
  fprintf(file, "      </CellData>\n");

  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n", "node_id");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "          %d\n", i);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v_l");
  for(unsigned i = 0; i < _edges.size(); i++)
    fprintf(file, "          %20.6f\n", _edges[i]->v_l());// v_l[i]);
  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
    fprintf(file, "          %20.6f\n", _edges[0]->v_l());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "v_g");
  for(unsigned i = 0; i < _edges.size(); i++)
    fprintf(file, "          %20.6f\n", _edges[i]->v_g());//v_g[i]);
  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
    fprintf(file, "          %20.6f\n", _edges[0]->v_g());
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");

  fprintf(file, "    </Piece>\n");
}

void
FourEqnOneP::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s%20s\n", "x", "alpha", "p", "rho_l", "rho_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, alpha[i], p[i], rho_l[i], rho_g[i]);

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s%20s\n", "x", "v_l", "v_g");
  for (unsigned i = 0; i < n_Node; i++)
    fprintf(file, "%20.6e%20.6e%20.6e\n", i*dx, v_l[i], v_g[i]);
}

void
FourEqnOneP::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  for(auto& cell : _cells)   { mnzp->addRow(cell->aDOF(), cell->getConnectedDOFs());  mnzp->addRow(cell->pDOF(), cell->getConnectedDOFs());}
  for(auto& edge : _edges)   { mnzp->addRow(edge->vlDOF(), edge->getConnectedDOFs()); mnzp->addRow(edge->vgDOF(), edge->getConnectedDOFs());}
}

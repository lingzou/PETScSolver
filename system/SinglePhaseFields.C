#include <cmath>
#include "SinglePhaseFields.h"

EdgeBase*
SPCell::getOtherSideEdge(EdgeBase* edge)
{
  return (edge == WEST_EDGE) ? EAST_EDGE : WEST_EDGE;
}

void
SPCell::initialize(double p, double T)
{
  _p = p;   _p_o = p;   _p_oo = p;
  _T = T;   _T_o = T;   _T_oo = T;
  _rho = rho_func(p, T);
  _rho_o = _rho; _rho_oo = _rho;
  _e = e_func(p, T);
  _e_o = _e; _e_oo = _e;
}

void
SPCell::updateSolution(double p, double T)
{
  _p = p; _T = T;
  _rho = rho_func(p, T);
  _e = e_func(p, T);
}

double
SPCell::computeMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux() - WEST_EDGE->mass_flux()) / dx;
}

double
SPCell::computeEnergyRHS(double dx)
{
  double p_dv_dx = _p * (EAST_EDGE->v() - WEST_EDGE->v()) / dx;

  double _h = 2000.0; double _aw = 300.0; double Tw = 350.0;
  double src = _h * _aw * (Tw - _T);

  return -(EAST_EDGE->energy_flux() - WEST_EDGE->energy_flux()) / dx - p_dv_dx + src;
}

void
SPCell::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  " << WEST_EDGE->name() << ", " << EAST_EDGE->name() << std::endl;
}

void
SPCell::saveOldSlns()
{
  _p_oo = _p_o; _T_oo = _T_o; _rho_oo = _rho_o; _e_oo = _e_o;
  _p_o  = _p;   _T_o  = _T;   _rho_o  = _rho;   _e_o = _e;
}

void
EdgeBase::setNghbrCells(SPCell * west, SPCell * east)
{
  WEST_CELL = west;
  EAST_CELL = east;
  if (WEST_CELL != NULL) WEST_EDGE = WEST_CELL->getOtherSideEdge(this);
  if (EAST_CELL != NULL) EAST_EDGE = EAST_CELL->getOtherSideEdge(this);
}

void
EdgeBase::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout          << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << std::endl;
}

void
vBndryEdge::computeFluxes()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p() : WEST_CELL->p();
  double rho_ghost = rho_func(p_ghost, _T_bc);
  double e_ghost = e_func(p_ghost, _T_bc);
  _mass_flux = rho_ghost * _v;
  _energy_flux = rho_ghost * e_ghost * _v;
}

void
pBndryEdge::computeFluxes()
{
  double rho_west = (WEST_CELL == NULL) ? rho_func(_p_bc, _T_bc) : WEST_CELL->rho();
  double rho_east = (EAST_CELL == NULL) ? rho_func(_p_bc, _T_bc) : EAST_CELL->rho();
  double e_west = (WEST_CELL == NULL) ? e_func(_p_bc, _T_bc) : WEST_CELL->e();
  double e_east = (EAST_CELL == NULL) ? e_func(_p_bc, _T_bc) : EAST_CELL->e();

  _mass_flux = (_v > 0) ? _v * rho_west : _v * rho_east;
  _energy_flux = (_v > 0) ? _v * rho_west * e_west : _v * rho_east * e_east;
}

double
pBndryEdge::computeTranRes(double dt)
{
  double rho_edge = (WEST_CELL == NULL) ? EAST_CELL->rho() : WEST_CELL->rho();
  return rho_edge * (_v - _v_o) / dt;
}

double
pBndryEdge::computeRHS(double dx)
{
  double dp_dx = 0;
  double dv_dx = 0;
  if (WEST_CELL == NULL)
  {
    dp_dx = (EAST_CELL->p() - _p_bc) / dx * 2;
    dv_dx = (_v < 0) ? (EAST_EDGE->v() - _v) / dx : 0;
  }
  else
  {
    dp_dx = (_p_bc - WEST_CELL->p()) / dx * 2;
    dv_dx = (_v > 0) ? (_v - WEST_EDGE->v()) / dx : 0;
  }

  double _f = 0.01; double _dh = 0.01;
  double rho_edge = (WEST_CELL == NULL) ? EAST_CELL->rho() : WEST_CELL->rho();
  double fric = 0.5 * _f / _dh * rho_edge * _v * std::fabs(_v);

  return -rho_edge * _v * dv_dx - dp_dx - fric;
}

void
IntEdge::computeFluxes()
{
  _mass_flux = (_v > 0) ? _v * WEST_CELL->rho() : _v * EAST_CELL->rho();
  _energy_flux = (_v > 0) ? _v * WEST_CELL->rho() * WEST_CELL->e() : _v * EAST_CELL->rho() * EAST_CELL->e();
}

double
IntEdge::computeTranRes(double dt)
{
  double rho_edge = 0.5 * (WEST_CELL->rho() + EAST_CELL->rho());
  return rho_edge * (_v - _v_o) / dt;
}

double
IntEdge::computeRHS(double dx)
{
  double dp_dx = (EAST_CELL->p() - WEST_CELL->p()) / dx;

  double _f = 0.01; double _dh = 0.01;
  double rho_edge = 0.5 * (WEST_CELL->rho() + EAST_CELL->rho());
  double fric = 0.5 * _f / _dh * rho_edge * _v * std::fabs(_v);

  double dv_dx = (_v > 0) ? (_v - WEST_EDGE->v()) / dx : (EAST_EDGE->v() - _v) / dx;

  return -rho_edge * _v * dv_dx - dp_dx - fric;
}

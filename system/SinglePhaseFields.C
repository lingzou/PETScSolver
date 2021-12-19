#include <cmath>
#include "utils.h"
#include "SinglePhaseFields.h"

EdgeBase*
SPCell::getOtherSideEdge(EdgeBase* edge)
{
  return (edge == WEST_EDGE) ? EAST_EDGE : WEST_EDGE;
}

void
SPCell::setExtendedNghbrs()
{
  if (WEST_EDGE != NULL) WEST_CELL = WEST_EDGE->getOtherSideCell(this);
  if (EAST_EDGE != NULL) EAST_CELL = EAST_EDGE->getOtherSideCell(this);
}

void
SPCell::initialize(double p, double T)
{
  _p = p;   _p_o = p;   _p_oo = p;
  _T = T;   _T_o = T;   _T_oo = T;
  _rho = _fluid->rho(p, T);
  _rho_o = _rho; _rho_oo = _rho;
  _e = _fluid->e(p, T);
  _e_o = _e; _e_oo = _e;
}

void
SPCell::linearReconstruction(double p_W, double p_E, double T_W, double T_E)
{
  UTILS::linearReconstruction(p_W, p_E, _p, _p_w, _p_e);
  UTILS::linearReconstruction(T_W, T_E, _T, _T_w, _T_e);

  _rho_w = _fluid->rho(_p_w, _T_w);
  _rho_e = _fluid->rho(_p_e, _T_e);
  _e_w = _fluid->e(_p_w, _T_w);
  _e_e = _fluid->e(_p_e, _T_e);
}

void
SPCell::updateSolution(double p, double T)
{
  _p = p; _T = T;
  _rho = _fluid->rho(p, T);
  _e = _fluid->e(p, T);
}

void
SPCell::SPCell::computeDP(double f, double dh, double gx)
{
  double v_cell = 0.5 * (EAST_EDGE->v() + WEST_EDGE->v());

  _dpdx_fric = 0.5 * f / dh * _rho * v_cell * std::fabs(v_cell);
  _dpdx_gravity = _rho * gx;
}

double
SPCell::computeMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux() - WEST_EDGE->mass_flux()) / dx;
}

double
SPCell::computeEnergyRHS(double dx, double h, double aw, double Tw)
{
  double p_dv_dx = _p * (EAST_EDGE->v() - WEST_EDGE->v()) / dx;

  double src = h * aw * (Tw - _T);

  return -(EAST_EDGE->energy_flux() - WEST_EDGE->energy_flux()) / dx - p_dv_dx + src;
}

void
SPCell::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout          << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << std::endl;
}

std::vector<unsigned>
SPCell::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.push_back(_pDOF);  dofs.push_back(_TDOF);
  if (WEST_CELL != NULL) { dofs.push_back(WEST_CELL->pDOF()); dofs.push_back(WEST_CELL->TDOF());}
  if (WEST_EDGE != NULL) { dofs.push_back(WEST_EDGE->vDOF()); }
  if (EAST_CELL != NULL) { dofs.push_back(EAST_CELL->pDOF()); dofs.push_back(EAST_CELL->TDOF());}
  if (EAST_EDGE != NULL) { dofs.push_back(EAST_EDGE->vDOF()); }
  return dofs;
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
}

SPCell*
EdgeBase::getOtherSideCell(SPCell* cell)
{
  return (cell == WEST_CELL) ? EAST_CELL : WEST_CELL;
}

void
EdgeBase::setExtendedNghbrs()
{
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

std::vector<unsigned>
EdgeBase::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.push_back(_vDOF);
  if (WEST_CELL != NULL) { dofs.push_back(WEST_CELL->pDOF()); dofs.push_back(WEST_CELL->TDOF());}
  if (WEST_EDGE != NULL) { dofs.push_back(WEST_EDGE->vDOF()); }
  if (EAST_CELL != NULL) { dofs.push_back(EAST_CELL->pDOF()); dofs.push_back(EAST_CELL->TDOF());}
  if (EAST_EDGE != NULL) { dofs.push_back(EAST_EDGE->vDOF()); }
  return dofs;
}

void
vBndryEdge::computeFluxes()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p() : WEST_CELL->p();
  double rho_ghost = _fluid->rho(p_ghost, _T_bc);
  double e_ghost = _fluid->e(p_ghost, _T_bc);
  _mass_flux = rho_ghost * _v;
  _energy_flux = rho_ghost * e_ghost * _v;
}

void
vBndryEdge::computeFluxes2nd()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p_w() : WEST_CELL->p_e();
  double rho_ghost = _fluid->rho(p_ghost, _T_bc);
  double e_ghost = _fluid->e(p_ghost, _T_bc);
  _mass_flux = rho_ghost * _v;
  _energy_flux = rho_ghost * e_ghost * _v;
}

void
pBndryEdge::computeFluxes()
{
  double rho_west = (WEST_CELL == NULL) ? _fluid->rho(_p_bc, _T_bc) : WEST_CELL->rho();
  double rho_east = (EAST_CELL == NULL) ? _fluid->rho(_p_bc, _T_bc) : EAST_CELL->rho();
  double e_west = (WEST_CELL == NULL) ? _fluid->e(_p_bc, _T_bc) : WEST_CELL->e();
  double e_east = (EAST_CELL == NULL) ? _fluid->e(_p_bc, _T_bc) : EAST_CELL->e();

  _mass_flux = (_v > 0) ? _v * rho_west : _v * rho_east;
  _energy_flux = (_v > 0) ? _v * rho_west * e_west : _v * rho_east * e_east;
}

void
pBndryEdge::computeFluxes2nd()
{
  double rho_west = (WEST_CELL == NULL) ? _fluid->rho(_p_bc, _T_bc) : WEST_CELL->rho_e();
  double rho_east = (EAST_CELL == NULL) ? _fluid->rho(_p_bc, _T_bc) : EAST_CELL->rho_w();
  double e_west = (WEST_CELL == NULL) ? _fluid->e(_p_bc, _T_bc) : WEST_CELL->e_e();
  double e_east = (EAST_CELL == NULL) ? _fluid->e(_p_bc, _T_bc) : EAST_CELL->e_w();

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
pBndryEdge::computeTranResBDF2(double dt)
{
  double rho_edge = (WEST_CELL == NULL) ? EAST_CELL->rho() : WEST_CELL->rho();
  return rho_edge * (1.5 * _v - 2.0 * _v_o + 0.5 * _v_oo) / dt;
}

double
pBndryEdge::computeRHS(double dx)
{
  double dp_dx = 0, dv_dx = 0, rho_edge = 0, fric = 0, gravity = 0;
  if (WEST_CELL == NULL)
  {
    dp_dx = (EAST_CELL->p() - _p_bc) / dx * 2;
    dv_dx = (_v < 0) ? (EAST_EDGE->v() - _v) / dx : 0;
    rho_edge = EAST_CELL->rho();
    fric = EAST_CELL->dpdx_fric();
    gravity = EAST_CELL->dpdx_gravity();
  }
  else
  {
    dp_dx = (_p_bc - WEST_CELL->p()) / dx * 2;
    dv_dx = (_v > 0) ? (_v - WEST_EDGE->v()) / dx : 0;
    rho_edge = WEST_CELL->rho();
    fric = WEST_CELL->dpdx_fric();
    gravity = WEST_CELL->dpdx_gravity();
  }

  return -rho_edge * _v * dv_dx - dp_dx - fric + gravity;
}

void
IntEdge::computeFluxes()
{
  _mass_flux = (_v > 0) ? _v * WEST_CELL->rho() : _v * EAST_CELL->rho();
  _energy_flux = (_v > 0) ? _v * WEST_CELL->rho() * WEST_CELL->e() : _v * EAST_CELL->rho() * EAST_CELL->e();
}

void
IntEdge::computeFluxes2nd()
{
  _mass_flux = (_v > 0) ? _v * WEST_CELL->rho_e() : _v * EAST_CELL->rho_w();
  _energy_flux = (_v > 0) ? _v * WEST_CELL->rho_e() * WEST_CELL->e_e() : _v * EAST_CELL->rho_w() * EAST_CELL->e_w();
}

double
IntEdge::computeTranRes(double dt)
{
  double rho_edge = 0.5 * (WEST_CELL->rho() + EAST_CELL->rho());
  return rho_edge * (_v - _v_o) / dt;
}

double
IntEdge::computeTranResBDF2(double dt)
{
  double rho_edge = 0.5 * (WEST_CELL->rho() + EAST_CELL->rho());
  return rho_edge * (1.5 * _v - 2.0 * _v_o + 0.5 * _v_oo) / dt;
}

double
IntEdge::computeRHS(double dx)
{
  double dp_dx = (EAST_CELL->p() - WEST_CELL->p()) / dx;

  double rho_edge = 0.5 * (WEST_CELL->rho() + EAST_CELL->rho());
  double fric = 0.5 * (WEST_CELL->dpdx_fric() + EAST_CELL->dpdx_fric());
  double gravity = 0.5 * (WEST_CELL->dpdx_gravity() + EAST_CELL->dpdx_gravity());

  double dv_dx = (_v > 0) ? (_v - WEST_EDGE->v()) / dx : (EAST_EDGE->v() - _v) / dx;

  return -rho_edge * _v * dv_dx - dp_dx - fric + gravity;
}

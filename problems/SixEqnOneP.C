#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "SixEqnOneP.h"
#include "utils.h"

// Reference:
//   [1] L. Zou, H. Zhao, and H. Zhang,
//        Solving phase appearance/disappearance two-phase flow problems with high resolution staggered grid and
//        fully implicit schemes by the Jacobian-free Newton–Krylov Method
//        Computers & Fluids, Vol. 129, 179-188
//
//   [2] L. Zou, H. Zhao, and H. Zhang,
//        Implicitly solving phase appearance and disappearance problems using two-fluid six-equation model,
//        Progress in Nuclear Energy 88 (2016) 198-210
//
//   The implementation in this code follows referene [2].

/* Staggered-grid mesh arrangement

     cell 0       1         2                            n-1
   |---------|---------|---------|---------|---------|---------|
   0    2    6    8                                            6n
   1    3    7    9                                            6n+1
        4         10
        5         11
  v_l  alpha
  v_g  p
       T_l
       T_g
*/


/**** Data Structure for 6E1P problems ****/
EdgeBase6E1P*
Cell6E1P::getOtherSideEdge(EdgeBase6E1P* edge)
{
  return (edge == WEST_EDGE) ? EAST_EDGE : WEST_EDGE;
}

void
Cell6E1P::setExtendedNghbrs()
{
  if (WEST_EDGE != NULL) WEST_CELL = WEST_EDGE->getOtherSideCell(this);
  if (EAST_EDGE != NULL) EAST_CELL = EAST_EDGE->getOtherSideCell(this);
}

void
Cell6E1P::initialize(double alpha, double p, double T_l, double T_g)
{
  _alpha_oo = _alpha_o    = _alpha  = alpha;
  _p_oo     = _p_o        = _p      = p;
  _T_l_oo   = _T_l_o      = _T_l    = T_l;
  _T_g_oo   = _T_g_o      = _T_g    = T_g;
  _rho_l_oo = _rho_l_o    = _rho_l  = _water->rho(p, T_l);
  _rho_g_oo = _rho_g_o    = _rho_g  = _steam->rho(p, T_g);
  _e_l_oo   = _e_l_o      = _e_l    = _water->e(p, T_l);
  _e_g_oo   = _e_g_o      = _e_g    = _steam->e(p, T_g);
}

void
Cell6E1P::linearReconstruction(double alpha_W, double alpha_E, double p_W, double p_E, double T_l_W, double T_l_E, double T_g_W, double T_g_E)
{
  UTILS::linearReconstruction(alpha_W, alpha_E, _alpha, _alpha_w, _alpha_e);
  UTILS::linearReconstruction(p_W, p_E, _p, _p_w, _p_e);
  UTILS::linearReconstruction(T_l_W, T_l_E, _T_l, _T_l_w, _T_l_e);
  UTILS::linearReconstruction(T_g_W, T_g_E, _T_g, _T_g_w, _T_g_e);

  _rho_l_w = _water->rho(_p_w, _T_l_w);
  _rho_l_e = _water->rho(_p_e, _T_l_e);
  _rho_g_w = _steam->rho(_p_w, _T_g_w);
  _rho_g_e = _steam->rho(_p_e, _T_g_e);
  _e_l_w = _water->e(_p_w, _T_l_w);
  _e_l_e = _water->e(_p_e, _T_l_e);
  _e_g_w = _steam->e(_p_w, _T_g_w);
  _e_g_e = _steam->e(_p_e, _T_g_e);
}

void
Cell6E1P::updateSolution(double alpha, double p, double T_l, double T_g)
{
  _alpha = alpha; _p = p; _T_l = T_l; _T_g = T_g;
  _rho_l = _water->rho(_p, _T_l);
  _rho_g = _steam->rho(_p, _T_g);
  _e_l = _water->e(_p, _T_l);
  _e_g = _steam->e(_p, _T_g);
}

double
Cell6E1P::liquidEnergyTranRes(double dt)
{
  double dt_term = ((1-_alpha) * _rho_l * _e_l - (1-_alpha_o) * _rho_l_o * _e_l_o) / dt;
  double dalpha_dt = (_alpha - _alpha_o) / dt;
  return dt_term - _p * dalpha_dt;
}

double
Cell6E1P::gasEnergyTranRes(double dt)
{
  double dt_term = (_alpha * _rho_g * _e_g - _alpha_o * _rho_g_o * _e_g_o) / dt;
  double dalpha_dt = (_alpha - _alpha_o) / dt;
  return dt_term + _p * dalpha_dt;
}

double
Cell6E1P::liquidEnergyTranResBDF2(double dt)
{
  double dt_term = (1.5 * (1-_alpha) * _rho_l * _e_l - 2.0 * (1-_alpha_o) * _rho_l_o * _e_l_o + 0.5 * (1-_alpha_oo) * _rho_l_oo * _e_l_oo) / dt;
  double dalpha_dt = (1.5 * _alpha - 2.0 * _alpha_o + 0.5 * _alpha_oo) / dt;

  return dt_term - _p * dalpha_dt;
}

double
Cell6E1P::gasEnergyTranResBDF2(double dt)
{
  double dt_term = (1.5 * _alpha * _rho_g * _e_g - 2.0 * _alpha_o * _rho_g_o * _e_g_o + 0.5 * _alpha_oo * _rho_g_oo * _e_g_oo) / dt;
  double dalpha_dt = (1.5 * _alpha - 2.0 * _alpha_o + 0.5 * _alpha_oo) / dt;
  return dt_term + _p * dalpha_dt;
}

double
Cell6E1P::liquidMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux_l() - WEST_EDGE->mass_flux_l()) / dx;
}

double
Cell6E1P::gasMassRHS(double dx)
{
  return -(EAST_EDGE->mass_flux_g() - WEST_EDGE->mass_flux_g()) / dx;
}

double
Cell6E1P::liquidEnergyRHS(double dx)
{
  double flux_term = -(EAST_EDGE->energy_flux_l() - WEST_EDGE->energy_flux_l()) / dx;
  double d_alphal_ul_dx = (EAST_EDGE->alphal_ul() - WEST_EDGE->alphal_ul()) / dx;

  return flux_term - _p * d_alphal_ul_dx; // + 10000 * (_T_g - _T_l);
}

double
Cell6E1P::gasEnergyRHS(double dx)
{
  double flux_term = -(EAST_EDGE->energy_flux_g() - WEST_EDGE->energy_flux_g()) / dx;
  double d_alphag_ug_dx = (EAST_EDGE->alphag_ug() - WEST_EDGE->alphag_ug()) / dx;

  return flux_term - _p * d_alphag_ug_dx; // + 10000 * (_T_l - _T_g);
}

void
Cell6E1P::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout          << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << std::endl;
}

std::vector<unsigned>
Cell6E1P::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.insert(dofs.end(), {_aDOF, _pDOF, _TlDOF, _TgDOF});
  if (WEST_CELL != NULL)
  {
    dofs.insert(dofs.end(), {WEST_CELL->aDOF(), WEST_CELL->pDOF(), WEST_CELL->TlDOF(), WEST_CELL->TgDOF()});

    Cell6E1P * WW_CELL = WEST_CELL->west_cell();
    if (WW_CELL != NULL)
      dofs.insert(dofs.end(), {WW_CELL->aDOF(), WW_CELL->pDOF(), WW_CELL->TlDOF(), WW_CELL->TgDOF()});
  }
  if (WEST_EDGE != NULL) dofs.insert(dofs.end(), {WEST_EDGE->vlDOF(), WEST_EDGE->vgDOF()});
  if (EAST_CELL != NULL)
  {
    dofs.insert(dofs.end(), {EAST_CELL->aDOF(), EAST_CELL->pDOF(), EAST_CELL->TlDOF(), EAST_CELL->TgDOF()});

    Cell6E1P * EE_CELL = EAST_CELL->east_cell();
    if (EE_CELL != NULL)
      dofs.insert(dofs.end(), {EE_CELL->aDOF(), EE_CELL->pDOF(), EE_CELL->TlDOF(), EE_CELL->TgDOF()});
  }
  if (EAST_EDGE != NULL) dofs.insert(dofs.end(), {EAST_EDGE->vlDOF(), EAST_EDGE->vgDOF()});
  return dofs;
}

void
Cell6E1P::saveOldSlns()
{
  _p_oo = _p_o; _alpha_oo = _alpha_o; _T_l_oo = _T_l_o;   _T_g_oo = _T_g_o; _rho_l_oo = _rho_l_o; _rho_g_oo = _rho_g_o; _e_l_oo = _e_l_o; _e_g_oo = _e_g_o;
  _p_o  = _p;   _alpha_o  = _alpha;   _T_l_o  = _T_l;     _T_g_o  = _T_g;   _rho_l_o  = _rho_l;   _rho_g_o  = _rho_g;   _e_l_o  = _e_l;   _e_g_o  = _e_g;
}

void
EdgeBase6E1P::setNghbrCells(Cell6E1P * west, Cell6E1P * east)
{
  WEST_CELL = west;
  EAST_CELL = east;
}

Cell6E1P*
EdgeBase6E1P::getOtherSideCell(Cell6E1P* cell)
{
  return (cell == WEST_CELL) ? EAST_CELL : WEST_CELL;
}

void
EdgeBase6E1P::setExtendedNghbrs()
{
  if (WEST_CELL != NULL) WEST_EDGE = WEST_CELL->getOtherSideEdge(this);
  if (EAST_CELL != NULL) EAST_EDGE = EAST_CELL->getOtherSideEdge(this);
}

void
EdgeBase6E1P::linearReconstruction(double vl_W, double vl_E, double vg_W, double vg_E)
{
  UTILS::linearReconstruction(vl_W, vl_E, _vl, _vl_w, _vl_e);
  UTILS::linearReconstruction(vg_W, vg_E, _vg, _vg_w, _vg_e);
}

void
EdgeBase6E1P::printConnection()
{
  std::cout << _name << std::endl;
  std::cout << "  "  << ((WEST_EDGE == NULL) ? "NULL" : WEST_EDGE->name()) << " -> ";
  std::cout          << ((WEST_CELL == NULL) ? "NULL" : WEST_CELL->name()) << " -> ";
  std::cout << _name << " -> ";
  std::cout          << ((EAST_CELL == NULL) ? "NULL" : EAST_CELL->name()) << " -> ";
  std::cout          << ((EAST_EDGE == NULL) ? "NULL" : EAST_EDGE->name()) << std::endl;
}

std::vector<unsigned>
EdgeBase6E1P::getConnectedDOFs()
{
  std::vector<unsigned> dofs;
  dofs.push_back(_vlDOF); dofs.push_back(_vgDOF);
  if (WEST_CELL != NULL)
  {
    dofs.push_back(WEST_CELL->aDOF());    dofs.push_back(WEST_CELL->pDOF());
    dofs.push_back(WEST_CELL->TlDOF());   dofs.push_back(WEST_CELL->TgDOF());
  }
  if (WEST_EDGE != NULL)
  {
    dofs.push_back(WEST_EDGE->vlDOF());  dofs.push_back(WEST_EDGE->vgDOF());
    EdgeBase6E1P * WW_EDGE = WEST_EDGE->west_edge();
    if (WW_EDGE != NULL ) { dofs.push_back(WW_EDGE->vlDOF());  dofs.push_back(WW_EDGE->vgDOF()); }
  }
  if (EAST_CELL != NULL)
  {
    dofs.push_back(EAST_CELL->aDOF());    dofs.push_back(EAST_CELL->pDOF());
    dofs.push_back(EAST_CELL->TlDOF());   dofs.push_back(EAST_CELL->TgDOF());
  }
  if (EAST_EDGE != NULL)
  {
    dofs.push_back(EAST_EDGE->vlDOF());  dofs.push_back(EAST_EDGE->vgDOF());
    EdgeBase6E1P * EE_EDGE = EAST_EDGE->east_edge();
    if (EE_EDGE != NULL ) { dofs.push_back(EE_EDGE->vlDOF());  dofs.push_back(EE_EDGE->vgDOF()); }
  }
  return dofs;
}

double
EdgeBase6E1P::interfacial_drag(double alpha_edge, double rho_l_edge, double rho_g_edge)
{
  double rp = 5e-4;
  double cd = 0.44;
  double ALPHA_MIN = 1e-6;
  double a_int = 3 * std::max(alpha_edge * (1 - alpha_edge), ALPHA_MIN * (1 - ALPHA_MIN)) / rp;
  double rho_mix = alpha_edge * rho_g_edge + (1 - alpha_edge) * rho_l_edge;
  double v_diff_sqr = (_vg - _vl) * std::fabs(_vg - _vl);

  return 0.125 * cd * a_int * rho_mix * v_diff_sqr;
}

void
vBndryEdge6E1P::computeFluxes()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p() : WEST_CELL->p();
  double rho_l_ghost = _water->rho(p_ghost, _T_l_bc);
  double rho_g_ghost = _steam->rho(p_ghost, _T_g_bc);
  double e_l_ghost = _water->e(p_ghost, _T_l_bc);
  double e_g_ghost = _steam->e(p_ghost, _T_g_bc);

  _mass_flux_l = (1 - _alpha_bc) * rho_l_ghost * _vl;
  _mass_flux_g = _alpha_bc * rho_g_ghost * _vg;
  _energy_flux_l = (1 - _alpha_bc) * rho_l_ghost * e_l_ghost * _vl;
  _energy_flux_g = _alpha_bc * rho_g_ghost * e_g_ghost * _vg;

  _alphal_ul = (1 - _alpha_bc) * _vl;
  _alphag_ug = _alpha_bc * _vg;
}

void
vBndryEdge6E1P::computeFluxes2nd()
{
  double p_ghost = (WEST_CELL == NULL) ? EAST_CELL->p_w() : WEST_CELL->p_e();
  double rho_l_ghost = _water->rho(p_ghost, _T_l_bc);
  double rho_g_ghost = _steam->rho(p_ghost, _T_g_bc);
  double e_l_ghost = _water->e(p_ghost, _T_l_bc);
  double e_g_ghost = _steam->e(p_ghost, _T_g_bc);

  _mass_flux_l = (1 - _alpha_bc) * rho_l_ghost * _vl;
  _mass_flux_g = _alpha_bc * rho_g_ghost * _vg;
  _energy_flux_l = (1 - _alpha_bc) * rho_l_ghost * e_l_ghost * _vl;
  _energy_flux_g = _alpha_bc * rho_g_ghost * e_g_ghost * _vg;

  _alphal_ul = (1 - _alpha_bc) * _vl;
  _alphag_ug = _alpha_bc * _vg;
}

void
pBndryEdge6E1P::computeFluxes()
{
  double alpha_west = 0, alpha_east = 0;
  double rho_l_west = 0, rho_l_east = 0;
  double rho_g_west = 0, rho_g_east = 0;
  double e_l_west = 0, e_l_east = 0;
  double e_g_west = 0, e_g_east = 0;
  if (WEST_CELL == NULL)
  {
    alpha_west = _alpha_bc;
    rho_l_west = _water->rho(_p_bc, _T_l_bc);
    rho_g_west = _steam->rho(_p_bc, _T_g_bc);
    e_l_west   = _water->e(_p_bc, _T_l_bc);
    e_g_west   = _steam->e(_p_bc, _T_g_bc);

    alpha_east = EAST_CELL->alpha();
    rho_l_east = EAST_CELL->rho_l();
    rho_g_east = EAST_CELL->rho_g();
    e_l_east   = EAST_CELL->e_l();
    e_g_east   = EAST_CELL->e_g();
  }
  else
  {
    alpha_west = WEST_CELL->alpha();
    rho_l_west = WEST_CELL->rho_l();
    rho_g_west = WEST_CELL->rho_g();
    e_l_west   = WEST_CELL->e_l();
    e_g_west   = WEST_CELL->e_g();

    alpha_east = _alpha_bc;
    rho_l_east = _water->rho(_p_bc, _T_l_bc);
    rho_g_east = _steam->rho(_p_bc, _T_g_bc);
    e_l_east   = _water->e(_p_bc, _T_l_bc);
    e_g_east   = _steam->e(_p_bc, _T_g_bc);
  }

  _mass_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west : _vl * (1 - alpha_east) * rho_l_east;
  _mass_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west : _vg * alpha_east * rho_g_east;
  _energy_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west * e_l_west
                             : _vl * (1 - alpha_east) * rho_l_east * e_l_east;
  _energy_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west * e_g_west
                             : _vg * alpha_east * rho_g_east * e_g_east;

  _alphal_ul = (_vl > 0) ? _vl * (1 - alpha_west) : _vl * (1 - alpha_east);
  _alphag_ug = (_vg > 0) ? _vg * alpha_west : _vg * alpha_east;
}

void
pBndryEdge6E1P::computeFluxes2nd()
{
  double alpha_west = 0, alpha_east = 0;
  double rho_l_west = 0, rho_l_east = 0;
  double rho_g_west = 0, rho_g_east = 0;
  double e_l_west = 0, e_l_east = 0;
  double e_g_west = 0, e_g_east = 0;
  if (WEST_CELL == NULL)
  {
    alpha_west = _alpha_bc;
    rho_l_west = _water->rho(_p_bc, _T_l_bc);
    rho_g_west = _steam->rho(_p_bc, _T_g_bc);
    e_l_west   = _water->e(_p_bc, _T_l_bc);
    e_g_west   = _steam->e(_p_bc, _T_g_bc);

    alpha_east = EAST_CELL->alpha_w();
    rho_l_east = EAST_CELL->rho_l_w();
    rho_g_east = EAST_CELL->rho_g_w();
    e_l_east   = EAST_CELL->e_l_w();
    e_g_east   = EAST_CELL->e_g_w();
  }
  else
  {
    alpha_west = WEST_CELL->alpha_e();
    rho_l_west = WEST_CELL->rho_l_e();
    rho_g_west = WEST_CELL->rho_g_e();
    e_l_west   = WEST_CELL->e_l_e();
    e_g_west   = WEST_CELL->e_g_e();

    alpha_east = _alpha_bc;
    rho_l_east = _water->rho(_p_bc, _T_l_bc);
    rho_g_east = _steam->rho(_p_bc, _T_g_bc);
    e_l_east   = _water->e(_p_bc, _T_l_bc);
    e_g_east   = _steam->e(_p_bc, _T_g_bc);
  }

  _mass_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west : _vl * (1 - alpha_east) * rho_l_east;
  _mass_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west : _vg * alpha_east * rho_g_east;
  _energy_flux_l = (_vl > 0) ? _vl * (1 - alpha_west) * rho_l_west * e_l_west
                             : _vl * (1 - alpha_east) * rho_l_east * e_l_east;
  _energy_flux_g = (_vg > 0) ? _vg * alpha_west * rho_g_west * e_g_west
                             : _vg * alpha_east * rho_g_east * e_g_east;

  _alphal_ul = (_vl > 0) ? _vl * (1 - alpha_west) : _vl * (1 - alpha_east);
  _alphag_ug = (_vg > 0) ? _vg * alpha_west : _vg * alpha_east;
}

double
pBndryEdge6E1P::liquidVelTranRes(double dt)
{
  double rho_l_edge = (WEST_CELL == NULL) ? EAST_CELL->rho_l() : WEST_CELL->rho_l();
  double alpha_edge = (WEST_CELL == NULL) ? EAST_CELL->alpha() : WEST_CELL->alpha();

  return (1 - alpha_edge) * rho_l_edge * (_vl - _vl_o) / dt;
}

double
pBndryEdge6E1P::gasVelTranRes(double dt)
{
  double rho_g_edge = (WEST_CELL == NULL) ? EAST_CELL->rho_g() : WEST_CELL->rho_g();
  double alpha_edge = (WEST_CELL == NULL) ? EAST_CELL->alpha() : WEST_CELL->alpha();

  return alpha_edge * rho_g_edge * (_vg - _vg_o) / dt;
}

double
pBndryEdge6E1P::liquidVelTranResBDF2(double dt)
{
  double rho_l_edge = (WEST_CELL == NULL) ? EAST_CELL->rho_l() : WEST_CELL->rho_l();
  double alpha_edge = (WEST_CELL == NULL) ? EAST_CELL->alpha() : WEST_CELL->alpha();

  return (1 - alpha_edge) * rho_l_edge * (1.5 * _vl - 2 * _vl_o + 0.5 * _vl_oo) / dt;
}

double
pBndryEdge6E1P::gasVelTranResBDF2(double dt)
{
  double rho_g_edge = (WEST_CELL == NULL) ? EAST_CELL->rho_g() : WEST_CELL->rho_g();
  double alpha_edge = (WEST_CELL == NULL) ? EAST_CELL->alpha() : WEST_CELL->alpha();

  return alpha_edge * rho_g_edge * (1.5 * _vg - 2 * _vg_o + 0.5 * _vg_oo) / dt;
}

void
pWESTBndryEdge6E1P::computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g)
{
  // FIXME
  /*
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

  double int_drag = compute_int_drag ? interfacial_drag(alpha_edge, rho_l_edge, rho_g_edge) : 0;
  rhs_l = -_vl * dv_l_dx - dp_dx / rho_l_edge + gravity_l + int_drag / alpha_edge_l / rho_l_edge;
  rhs_g = -_vg * dv_g_dx - dp_dx / rho_g_edge + gravity_g - int_drag / alpha_edge_g / rho_g_edge;*/
}

void
pEASTBndryEdge6E1P::computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g)
{
  double rho_l_edge = WEST_CELL->rho_l();
  double rho_g_edge = WEST_CELL->rho_g();
  double alpha_edge = WEST_CELL->alpha();
  double alpha_edge_l = std::max(1 - alpha_edge, 1e-6);
  double alpha_edge_g = std::max(alpha_edge, 1e-6);
  /*
  double dv_l_dx = 0;
  if (order == 1)
    dv_l_dx = (_vl > 0) ? (_vl - WEST_EDGE->v_l()) / dx : 0;
  else
    dv_l_dx = (_vl > 0) ? (_vl_e - WEST_EDGE->vl_e()) / dx : (_vl - WEST_EDGE->v_l()) / dx;
  */
  double dv_l_dx = (_vl > 0) ? (_vl - WEST_EDGE->v_l()) / dx : 0;
  double dp_dx = 2 * (_p_bc - WEST_CELL->p()) / dx;
  double gravity_l = (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l() * WEST_CELL->g();

  double dv_g_dx = (_vg > 0) ? (_vg - WEST_EDGE->v_g()) / dx : 0;
  /*
  double dv_g_dx = 0;
  if (order == 1)
    dv_g_dx = (_vg > 0) ? (_vg - WEST_EDGE->v_g()) / dx : 0;
  else
    dv_g_dx = (_vg > 0) ? (_vg_e - WEST_EDGE->vg_e()) / dx : (_vg - WEST_EDGE->v_g()) / dx;*/
  double gravity_g = WEST_CELL->alpha() * WEST_CELL->rho_g() * WEST_CELL->g();

  double int_drag = compute_int_drag ? interfacial_drag(alpha_edge, rho_l_edge, rho_g_edge) : 0;
  rhs_l = -(1 - alpha_edge) * rho_l_edge * _vl * dv_l_dx - (1 - alpha_edge) * dp_dx + gravity_l + int_drag;
  rhs_g = -alpha_edge       * rho_g_edge * _vg * dv_g_dx -       alpha_edge * dp_dx + gravity_g - int_drag;
}

void
IntEdge6E1P::computeFluxes()
{
  _mass_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l()
                           : _vl * (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l();

  _mass_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha() * WEST_CELL->rho_g()
                           : _vg * EAST_CELL->alpha() * EAST_CELL->rho_g();

  _energy_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l() * WEST_CELL->e_l()
                             : _vl * (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l() * EAST_CELL->e_l();

  _energy_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha() * WEST_CELL->rho_g() * WEST_CELL->e_g()
                             : _vg * EAST_CELL->alpha() * EAST_CELL->rho_g() * EAST_CELL->e_g();

  _alphal_ul = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha()) : _vl * (1 - EAST_CELL->alpha());
  _alphag_ug = (_vg > 0) ? _vg * WEST_CELL->alpha() : _vg * EAST_CELL->alpha();
}

void
IntEdge6E1P::computeFluxes2nd()
{
  _mass_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha_e()) * WEST_CELL->rho_l_e()
                           : _vl * (1 - EAST_CELL->alpha_w()) * EAST_CELL->rho_l_w();

  _mass_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha_e() * WEST_CELL->rho_g_e()
                           : _vg * EAST_CELL->alpha_w() * EAST_CELL->rho_g_w();

  _energy_flux_l = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha_e()) * WEST_CELL->rho_l_e() * WEST_CELL->e_l_e()
                             : _vl * (1 - EAST_CELL->alpha_w()) * EAST_CELL->rho_l_w() * EAST_CELL->e_l_w();

  _energy_flux_g = (_vg > 0) ? _vg * WEST_CELL->alpha_e() * WEST_CELL->rho_g_e() * WEST_CELL->e_g_e()
                             : _vg * EAST_CELL->alpha_w() * EAST_CELL->rho_g_w() * EAST_CELL->e_g_w();

  _alphal_ul = (_vl > 0) ? _vl * (1 - WEST_CELL->alpha_e()) : _vl * (1 - EAST_CELL->alpha_w());
  _alphag_ug = (_vg > 0) ? _vg * WEST_CELL->alpha_e() : _vg * EAST_CELL->alpha_w();
}

double
IntEdge6E1P::liquidVelTranRes(double dt)
{
  double rho_l_edge = 0.5 * (WEST_CELL->rho_l() + EAST_CELL->rho_l());
  double alpha_edge = 0.5 * (WEST_CELL->alpha() + EAST_CELL->alpha());

  return (1 - alpha_edge) * rho_l_edge * (_vl - _vl_o) / dt;
}

double
IntEdge6E1P::gasVelTranRes(double dt)
{
  double rho_g_edge = 0.5 * (WEST_CELL->rho_g() + EAST_CELL->rho_g());
  double alpha_edge = 0.5 * (WEST_CELL->alpha() + EAST_CELL->alpha());

  return alpha_edge * rho_g_edge * (_vg - _vg_o) / dt;
}

double
IntEdge6E1P::liquidVelTranResBDF2(double dt)
{
  double rho_l_edge = 0.5 * (WEST_CELL->rho_l() + EAST_CELL->rho_l());
  double alpha_edge = 0.5 * (WEST_CELL->alpha() + EAST_CELL->alpha());

  return (1 - alpha_edge) * rho_l_edge * (1.5 * _vl - 2 * _vl_o + 0.5 * _vl_oo) / dt;
}

double
IntEdge6E1P::gasVelTranResBDF2(double dt)
{
  double rho_g_edge = 0.5 * (WEST_CELL->rho_g() + EAST_CELL->rho_g());
  double alpha_edge = 0.5 * (WEST_CELL->alpha() + EAST_CELL->alpha());

  return alpha_edge * rho_g_edge * (1.5 * _vg - 2 * _vg_o + 0.5 * _vg_oo) / dt;
}

void
IntEdge6E1P::computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g)
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
    sysError("order not supported in IntEdge6E1P::computeLiquidRHS");

  double dp_dx = (EAST_CELL->p() - WEST_CELL->p()) / dx;

  double gravity_l = (1 - EAST_CELL->alpha()) * EAST_CELL->rho_l() * EAST_CELL->g()
                   + (1 - WEST_CELL->alpha()) * WEST_CELL->rho_l() * WEST_CELL->g();
  gravity_l *= 0.5;

  double gravity_g = EAST_CELL->alpha() * EAST_CELL->rho_g() * EAST_CELL->g()
                   + WEST_CELL->alpha() * WEST_CELL->rho_g() * WEST_CELL->g();
  gravity_g *= 0.5;

  double int_drag = compute_int_drag ? interfacial_drag(alpha_edge, rho_l_edge, rho_g_edge) : 0;
  rhs_l = -(1 - alpha_edge) * rho_l_edge * _vl * dv_l_dx - (1 - alpha_edge) * dp_dx + gravity_l + int_drag;
  rhs_g = -alpha_edge       * rho_g_edge * _vg * dv_g_dx -       alpha_edge * dp_dx + gravity_g - int_drag;
}
/**** End of Data Structure for 6E1P problems ****/

SixEqnOneP::SixEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _problem_type = _inputParamList.getValueFromInput<int>("problem_type");
  _order = _inputParamList.getValueFromInput<int>("order");
  length = _inputParamList.getValueFromInput<int>("length");
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");

  n_Node = n_Cell + 1;
  dx = length / n_Cell;

  _water = new Water_IF97(inputParamList); //inputParamList a dummy
  _steam = new Steam_IF97(inputParamList); //inputParamList a dummy

/*
  if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
  {
    _n_DOFs = n_Cell * 6;
    _compute_int_drag = false;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell6E1P("CELL_"+std::to_string(i)));
    for (unsigned i = 0; i < n_Cell; i++)   _edges.push_back(new IntEdge6E1P("EDGE_"+std::to_string(i)));

    for (unsigned i = 0; i < n_Cell-1; i++) _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);
    _cells[n_Cell-1]->setNghbrEdges(_edges[n_Cell-1], _edges[0]);

    _edges[0]->setNghbrCells(_cells[n_Cell-1], _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++) _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
  }*/
  if (_problem_type == MANOMETER)
  {
    _n_DOFs = n_Cell * 6 + 2;
    _compute_int_drag = true;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell6E1P("CELL_"+std::to_string(i)));

    _edges.push_back(new pWESTBndryEdge6E1P("INLET", 1e5, 1, 300, 300));
    for (unsigned i = 0; i < n_Cell-1; i++)   _edges.push_back(new IntEdge6E1P("EDGE_"+std::to_string(i)));
    _edges.push_back(new pEASTBndryEdge6E1P("OUTLET", 1e5, 1, 300, 300));

    for (unsigned i = 0; i < n_Cell; i++)   _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);

    _edges[0]->setNghbrCells(NULL, _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++)   _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
    _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  }
  else if (_problem_type == SEDIMENTATION)
  {
    _n_DOFs = n_Cell * 6 + 2;
    _compute_int_drag = true;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell6E1P("CELL_"+std::to_string(i)));

    _edges.push_back(new vBndryEdge6E1P("TOP", 0, 0, 1, 300, 300));
    for (unsigned i = 0; i < n_Cell-1; i++)   _edges.push_back(new IntEdge6E1P("EDGE_"+std::to_string(i)));
    _edges.push_back(new vBndryEdge6E1P("OUTLET", 0, 0, 1, 300, 300));

    for (unsigned i = 0; i < n_Cell; i++)   _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);

    _edges[0]->setNghbrCells(NULL, _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++)   _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
    _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  }
  else if (_problem_type == WATER_FAUCET)
  {
    _n_DOFs = n_Cell * 6 + 2;
    _compute_int_drag = false;

    for (unsigned i = 0; i < n_Cell; i++)   _cells.push_back(new Cell6E1P("CELL_"+std::to_string(i)));

    _edges.push_back(new vBndryEdge6E1P("TOP", 10, 0, 0.2, 450, 450));
    for (unsigned i = 0; i < n_Cell-1; i++)   _edges.push_back(new IntEdge6E1P("EDGE_"+std::to_string(i)));
    _edges.push_back(new pEASTBndryEdge6E1P("OUTLET", 1e5, 0.2, 450, 450));

    for (unsigned i = 0; i < n_Cell; i++)   _cells[i]->setNghbrEdges(_edges[i], _edges[i+1]);

    _edges[0]->setNghbrCells(NULL, _cells[0]);
    for (unsigned i = 1; i < n_Cell; i++)   _edges[i]->setNghbrCells(_cells[i-1], _cells[i]);
    _edges[n_Cell]->setNghbrCells(_cells[n_Cell-1], NULL);
  }
  else
    sysError("Problem type not supported.");

  for(auto& edge : _edges)   edge->setEOS(_water, _steam);
  for(auto& cell : _cells)   cell->setEOS(_water, _steam);
  std::cout << "SixEqnOneP::SixEqnOneP end" << std::endl;
}

SixEqnOneP::~SixEqnOneP()
{
  for(auto& itr : _edges)   delete itr;
  for(auto& itr : _cells)   delete itr;

  delete _water;
  delete _steam;
}

void
SixEqnOneP::setDOFoffset(unsigned offset)
{
  std::cout << "SixEqnOneP::setDOFoffset begin" << std::endl;
  _DOF_offset = offset;
  for(unsigned i = 0; i < _edges.size(); i++)   _edges[i]->setDOF(_DOF_offset + 6*i,      _DOF_offset + 6*i + 1);
  for(unsigned i = 0; i < _cells.size(); i++)   _cells[i]->setDOF(_DOF_offset + 6*i + 2,  _DOF_offset + 6*i + 3,  _DOF_offset + 6*i + 4,  _DOF_offset + 6*i + 5);
  std::cout << "SixEqnOneP::setDOFoffset end" << std::endl;
}

void
SixEqnOneP::setupExtendedConnections()
{
  std::cout << "SixEqnOneP::setupExtendedConnections begin" << std::endl;
  for(auto& itr : _edges)   itr->setExtendedNghbrs();
  for(auto& itr : _cells)   itr->setExtendedNghbrs();
  std::cout << "SixEqnOneP::setupExtendedConnections end" << std::endl;
  /*
  for(auto& itr : _cells)   itr->printConnection();
  for(auto& itr : _edges)   itr->printConnection();
  */
}

void
SixEqnOneP::SetupInitialCondition(double * u)
{
  std::cout << "SixEqnOneP::SetupInitialCondition begin" << std::endl;
  unsigned index = 0;
  switch (_problem_type)
  {
    case SINE_WAVE:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(1, 1);
        u[6*i] = 1;
        u[6*i+1] = 1;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double val = 0.5 + 0.2 * std::sin(2 * PI * xx / length);
        _cells[i]->initialize(val, 1e5, 300, 300);
        u[6*i+2] = val;
        u[6*i+3] = 1e5;
        u[6*i+4] = 300;
        u[6*i+5] = 300;

        _cells[i]->set_g(0);
      }
    break;

    case SQUARE_WAVE:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(1, 1);
        u[6*i] = 1;
        u[6*i+1] = 1;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double val = ((xx < 0.2) || (xx > 0.4)) ? 1 : 0;
        _cells[i]->initialize(val, 1e5, 300, 300);
        u[6*i+2] = val;
        u[6*i+3] = 1e5;
        u[6*i+4] = 300;
        u[6*i+5] = 300;

        _cells[i]->set_g(0);
      }
    break;

    case MANOMETER:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(-1, -1);
        u[6*i] = -1;
        u[6*i+1] = -1;
      }

      for(unsigned i = 0; i < _cells.size(); i++)
      {
        double xx = (i + 0.5) * dx;
        double alpha_init = 0, p_init = 0, depth = 0;
        if (xx < 5) // left gas column
        {
          _cells[i]->set_g(9.81);
          alpha_init = 1;
          p_init = 1e5 + 0.5 * 9.81 * xx;
        }
        else if (xx < 10) // left water column
        {
          _cells[i]->set_g(9.81);
          depth = xx - 5;
          alpha_init = 0;
          p_init = 1e5 + 0.5 * 9.81 * 5 + 1000 * 9.81 * depth;
        }
        else if (xx < 15) // right water column
        {
          _cells[i]->set_g(-9.81);
          depth = 15 - xx;
          alpha_init = 0;
          p_init = 1e5 + 0.5 * 9.81 * 5 + 1000 * 9.81 * depth;
        }
        else // right gas column
        {
          _cells[i]->set_g(-9.81);
          depth = 20 - xx;
          alpha_init = 1;
          p_init = 1e5 + 0.5 * 9.81 * xx;
        }

        _cells[i]->initialize(alpha_init, p_init, 300, 300);

        u[6*i+2] = alpha_init;
        u[6*i+3] = p_init;
        u[6*i+4] = 300;
        u[6*i+5] = 300;
      }
    break;

    case SEDIMENTATION:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(0, 0);
        u[6*i] = 0;
        u[6*i+1] = 0;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        _cells[i]->initialize(0.5, 1e5, 300, 300);
        u[6*i+2] = 0.5;
        u[6*i+3] = 1e5;
        u[6*i+4] = 300;
        u[6*i+5] = 300;

        _cells[i]->set_g(-9.81);
      }
    break;

    case WATER_FAUCET:
      for(unsigned i = 0; i < _edges.size(); i++)
      {
        _edges[i]->initialize(10, 0);
        u[6*i] = 10;
        u[6*i+1] = 0;
      }
      for(unsigned i = 0; i < _cells.size(); i++)
      {
        _cells[i]->initialize(0.2, 1e5, 450, 450);
        u[6*i+2] = 0.2;
        u[6*i+3] = 1e5;
        u[6*i+4] = 450;
        u[6*i+5] = 450;

        _cells[i]->set_g(9.81);
      }
    break;

    defaut: sysError("Unknown problem_type.");
  }
  std::cout << "SixEqnOneP::SetupInitialCondition end" << std::endl;
}

void
SixEqnOneP::updateSolution(double * u)
{
  for(unsigned i = 0; i < _edges.size(); i++)   _edges[i]->updateSolution(u[6*i], u[6*i+1]);    // vl, vg
  for(unsigned i = 0; i < _cells.size(); i++)   _cells[i]->updateSolution(u[6*i+2], u[6*i+3], u[6*i+4], u[6*i+5]);  // alpha, p, Tl, Tg
}

void
SixEqnOneP::transientResidual(double * res)
{
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    for(unsigned i = 0; i < _edges.size(); i++)
    {
      res[6*i] = _edges[i]->liquidVelTranResBDF2(_dt);
      res[6*i+1] = _edges[i]->gasVelTranResBDF2(_dt);
    }
    for(unsigned i = 0; i < _cells.size(); i++)
    {
      res[6*i+2] = _cells[i]->liquidMassTranResBDF2(_dt);
      res[6*i+3] = _cells[i]->gasMassTranResBDF2(_dt);
      res[6*i+4] = _cells[i]->liquidEnergyTranResBDF2(_dt);
      res[6*i+5] = _cells[i]->gasEnergyTranResBDF2(_dt);
    }
  }
  else
  {
    for(unsigned i = 0; i < _edges.size(); i++)
    {
      res[6*i] = _edges[i]->liquidVelTranRes(_dt);
      res[6*i+1] = _edges[i]->gasVelTranRes(_dt);
    }
    for(unsigned i = 0; i < _cells.size(); i++)
    {
      res[6*i+2] = _cells[i]->liquidMassTranRes(_dt);
      res[6*i+3] = _cells[i]->gasMassTranRes(_dt);
      res[6*i+4] = _cells[i]->liquidEnergyTranRes(_dt);
      res[6*i+5] = _cells[i]->gasEnergyTranRes(_dt);
    }
  }
}

void
SixEqnOneP::RHS(double * rhs)
{
  for(unsigned i = 0; i < _edges.size(); i++)
    _edges[i]->computeRHS(_compute_int_drag, _order, dx, rhs[6*i], rhs[6*i+1]);

  for(unsigned i = 0; i < _cells.size(); i++)
  {
    rhs[6*i+2] = _cells[i]->liquidMassRHS(dx);
    rhs[6*i+3] = _cells[i]->gasMassRHS(dx);
    rhs[6*i+4] = _cells[i]->liquidEnergyRHS(dx);
    rhs[6*i+5] = _cells[i]->gasEnergyRHS(dx);
  }
}

void
SixEqnOneP::linearReconstruction()
{
  if (_order == 2)
  {
    if ((_problem_type == SINE_WAVE) || (_problem_type == SQUARE_WAVE))
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {
        double alpha_W  = (i == 0) ? _cells.back()->alpha() : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? _cells.back()->p()     : _cells[i-1]->p();
        double T_l_W    = (i == 0) ? _cells.back()->T_l()   : _cells[i-1]->T_l();
        double T_g_W    = (i == 0) ? _cells.back()->T_g()   : _cells[i-1]->T_g();
        double alpha_E  = (i == n_Cell-1) ? _cells.front()->alpha() : _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? _cells.front()->p()     : _cells[i+1]->p();
        double T_l_E    = (i == n_Cell-1) ? _cells.front()->T_l()   : _cells[i+1]->T_l();
        double T_g_E    = (i == n_Cell-1) ? _cells.front()->T_g()   : _cells[i+1]->T_g();
        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E, T_l_W, T_l_E, T_g_W, T_g_E);
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
    else if (_problem_type == MANOMETER)
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {
        double alpha_W  = (i == 0) ? 1            : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? 1e5          : _cells[i-1]->p();
        double T_l_W    = (i == 0) ? 300          : _cells[i-1]->T_l();
        double T_g_W    = (i == 0) ? 300          : _cells[i-1]->T_g();
        double alpha_E  = (i == n_Cell-1) ? 1     : _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? 1e5   : _cells[i+1]->p();
        double T_l_E    = (i == n_Cell-1) ? 300   : _cells[i+1]->T_l();
        double T_g_E    = (i == n_Cell-1) ? 300   : _cells[i+1]->T_g();

        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E, T_l_W, T_l_E, T_g_W, T_g_E);
      }
      for (unsigned i = 0; i < _edges.size(); i++)
      {
        double vl_W     = (i == 0) ? 2*_edges[0]->v_l()-_edges[1]->v_l()  : _edges[i-1]->v_l();
        double vg_W     = (i == 0) ? 2*_edges[0]->v_g()-_edges[1]->v_g()  : _edges[i-1]->v_g();
        double vl_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_l()-_edges[n_Cell-1]->v_l() : _edges[i+1]->v_l();
        double vg_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_g()-_edges[n_Cell-1]->v_g() : _edges[i+1]->v_g();
        _edges[i]->linearReconstruction(vl_W, vl_E, vg_W, vg_E);
      }
    }
    else if (_problem_type == SEDIMENTATION)
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {
        double alpha_W  = (i == 0) ? _cells[0]->alpha()      : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? _cells[0]->p()          : _cells[i-1]->p();
        double T_l_W    = (i == 0) ? _cells[0]->T_l()        : _cells[i-1]->T_l();
        double T_g_W    = (i == 0) ? _cells[0]->T_g()        : _cells[i-1]->T_g();
        double alpha_E  = (i == n_Cell-1) ? _cells[n_Cell-1]->alpha()     : _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? _cells[n_Cell-1]->p()         : _cells[i+1]->p();
        double T_l_E    = (i == n_Cell-1) ? _cells[n_Cell-1]->T_l()       : _cells[i+1]->T_l();
        double T_g_E    = (i == n_Cell-1) ? _cells[n_Cell-1]->T_g()       : _cells[i+1]->T_g();

        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E, T_l_W, T_l_E, T_g_W, T_g_E);
      }
      for (unsigned i = 0; i < _edges.size(); i++)
      {
        double vl_W     = (i == 0) ? 2*_edges[0]->v_l()-_edges[1]->v_l()  : _edges[i-1]->v_l();
        double vg_W     = (i == 0) ? 2*_edges[0]->v_g()-_edges[1]->v_g()  : _edges[i-1]->v_g();
        double vl_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_l()-_edges[n_Cell-1]->v_l() : _edges[i+1]->v_l();
        double vg_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_g()-_edges[n_Cell-1]->v_g() : _edges[i+1]->v_g();
        _edges[i]->linearReconstruction(vl_W, vl_E, vg_W, vg_E);
      }
    }
    else if (_problem_type == WATER_FAUCET)
    {
      for (unsigned i = 0; i < _cells.size(); i++)
      {
        double alpha_W  = (i == 0) ? 0.2                : _cells[i-1]->alpha();
        double p_W      = (i == 0) ? _cells[0]->p()     : _cells[i-1]->p();
        double T_l_W    = (i == 0) ? 450   : _cells[i-1]->T_l();
        double T_g_W    = (i == 0) ? 450   : _cells[i-1]->T_g();

        double alpha_E  = 0;
        if (i == n_Cell-1)    alpha_E = (_edges[n_Cell]->v_l() > 0) ? _cells[n_Cell-1]->alpha() : 0.2;
        else                  alpha_E = _cells[i+1]->alpha();
        double p_E      = (i == n_Cell-1) ? 1e5   : _cells[i+1]->p();
        double T_l_E    = (i == n_Cell-1) ? 450   : _cells[i+1]->T_l();
        double T_g_E    = (i == n_Cell-1) ? 450   : _cells[i+1]->T_g();

        _cells[i]->linearReconstruction(alpha_W, alpha_E, p_W, p_E, T_l_W, T_l_E, T_g_W, T_g_E);
      }
      for (unsigned i = 0; i < _edges.size(); i++)
      {
        double vl_W     = (i == 0) ? 2*_edges[0]->v_l()-_edges[1]->v_l()  : _edges[i-1]->v_l();
        double vg_W     = (i == 0) ? 2*_edges[0]->v_g()-_edges[1]->v_g()  : _edges[i-1]->v_g();
        double vl_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_l()-_edges[n_Cell-1]->v_l() : _edges[i+1]->v_l();
        double vg_E     = (i == n_Cell) ? 2*_edges[n_Cell]->v_g()-_edges[n_Cell-1]->v_g() : _edges[i+1]->v_g();
        _edges[i]->linearReconstruction(vl_W, vl_E, vg_W, vg_E);
      }
    }
    else
      sysError("Not implemented");
  }
}

void
SixEqnOneP::updateEdgeCellHelperVar()
{
  switch (_order)
  {
    case 1:     for(auto& edge : _edges)    edge->computeFluxes();      break;
    case 2:     for(auto& edge : _edges)   edge->computeFluxes2nd();    break;
    default :   sysError("Spatial order not implemented.");
  }
}

void
SixEqnOneP::onTimestepEnd()
{
  // save old solutions
  for(auto& itr : _cells)   itr->saveOldSlns();
  for(auto& itr : _edges)   itr->saveOldSlns();

  if (_problem_type == MANOMETER)
  {
    time.push_back(_problemSystem->getCurrentTime());
    v_l_bottom.push_back(_edges[n_Cell/2]->v_l());
    p_bottom.push_back(0.5 * (_cells[n_Cell/2-1]->p()+_cells[n_Cell/2]->p()));

    double h_l, h_r;
    for(unsigned i = 0; i < n_Cell - 1; i++)
    {
      if ((_cells[i]->alpha() > 0.5) && (_cells[i+1]->alpha() <= 0.5))
        h_l = (i + 0.5) * dx + dx * (0.5 - _cells[i]->alpha()) / (_cells[i+1]->alpha() - _cells[i]->alpha());

      if ((_cells[i]->alpha() <= 0.5) && (_cells[i+1]->alpha() > 0.5))
        h_r = (i + 0.5) * dx + dx * (0.5 - _cells[i]->alpha()) / (_cells[i+1]->alpha() - _cells[i]->alpha());
    }
    h_left.push_back(10 - h_l);
    h_right.push_back(h_r - 10);
  }
}

void
SixEqnOneP::onLastTimestepEnd()
{
  if (_problem_type == MANOMETER)
  {
    FILE * file;
    std::string file_name = "output/manometer_results.csv";
    file = fopen(file_name.c_str(), "w");

    // Reference analytical solutions see reference:
    // Section 9.2 of:
    // H. St\"{a}dtke, Gasdynamic Aspects of Two-Phase Flow, Wiley-VCH, 2006.
    double L = 10; // length of total water column
    double T = sqrt(2 * PI * PI * L / 9.81);
    double h_max = sqrt(L / 2 / 9.81);

    fprintf(file, "%20s,%20s,%20s,%20s,%20s,%20s,%20s\n", "time", "v_l_bottom", "p_bottom", "h_left", "h_right", "v_ana", "h_left_ana");
    for (unsigned i = 0; i < time.size(); i++)
    {
      double tt = time[i];
      double v_ana = -cos(2 * PI * tt / T);
      double h_left_ana = 5 + h_max * sin(2 * PI * tt / T);
      fprintf(file, "%20.6e,%20.6e,%20.6e,%20.6e,%20.6e,%20.6e,%20.6e\n", tt, v_l_bottom[i], p_bottom[i], h_left[i], h_right[i], v_ana, h_left_ana);
    }

    fclose(file);
  }
}

void
SixEqnOneP::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "T_l");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->T_l());//rho_l[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "T_g");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", _cells[i]->T_g());//rho_g[i]);
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
SixEqnOneP::writeTextOutput(FILE * file)
{
  // cell data
  fprintf(file, "Time = %20.6e\n", _problemSystem->getCurrentTime());
  fprintf(file, "#Cell data\n");
  fprintf(file, "%20s%20s%20s%20s%20s\n", "x", "alpha", "p", "rho_l", "rho_g");
  for (unsigned i = 0; i < _cells.size(); i++)
    fprintf(file, "%20.6e%20.6e%20.6e%20.6e%20.6e\n", (i+0.5)*dx, _cells[i]->alpha(), _cells[i]->p(), _cells[i]->rho_l(), _cells[i]->rho_g());

  // edge data
  fprintf(file, "#Edge data\n");
  fprintf(file, "%20s%20s%20s\n", "x", "v_l", "v_g");
  for (unsigned i = 0; i < _edges.size(); i++)
    fprintf(file, "%20.6e%20.6e%20.6e\n", i*dx, _edges[i]->v_l(), _edges[i]->v_g());
}

void
SixEqnOneP::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  std::cout << "SixEqnOneP::FillJacobianMatrixNonZeroPattern begin" << std::endl;
  for(auto& cell : _cells)
  {
    mnzp->addRow(cell->aDOF(), cell->getConnectedDOFs());
    mnzp->addRow(cell->pDOF(), cell->getConnectedDOFs());
    mnzp->addRow(cell->TlDOF(), cell->getConnectedDOFs());
    mnzp->addRow(cell->TgDOF(), cell->getConnectedDOFs());
  }
  for(auto& edge : _edges)
  {
    mnzp->addRow(edge->vlDOF(), edge->getConnectedDOFs());
    mnzp->addRow(edge->vgDOF(), edge->getConnectedDOFs());
  }
  std::cout << "SixEqnOneP::FillJacobianMatrixNonZeroPattern end" << std::endl;
}

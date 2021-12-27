#ifndef SINGLE_PHASE_FIELDS_H
#define SINGLE_PHASE_FIELDS_H

#include <iostream>
#include "SinglePhaseFluid.h"

class EdgeBase;
class crossEdge;

class SPCell
{
public:
  SPCell(std::string name, SinglePhaseFluid* fluid) : _name(name), _fluid(fluid), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~SPCell() {}
  virtual std::string name() { return _name; }

  // Data access
  virtual double p()    const { return _p;    }
  virtual double T()    const { return _T;    }
  virtual double rho()  const { return _rho;  }
  virtual double e()    const { return _e;    }
  virtual double v_cell()         const { return _v_cell; }
  virtual double dpdx_fric()      const { return _dpdx_fric; }
  virtual double dpdx_gravity()   const { return _dpdx_gravity; }
  //   2nd-order related
  virtual double p_w()    const { return _p_w;    }
  virtual double p_e()    const { return _p_e;    }
  virtual double T_w()    const { return _T_w;    }
  virtual double T_e()    const { return _T_e;    }
  virtual double rho_w()  const { return _rho_w;  }
  virtual double rho_e()  const { return _rho_e;  }
  virtual double e_w()    const { return _e_w;    }
  virtual double e_e()    const { return _e_e;    }
  virtual void linearReconstruction(double, double, double, double);

  // Connection related functions
  virtual EdgeBase * getOtherSideEdge(EdgeBase* edge);
  virtual void setNghbrEdges(EdgeBase * west, EdgeBase * east) { WEST_EDGE = west; EAST_EDGE = east; }
  virtual void setExtendedNghbrs();

  // Residual related functions
  virtual void setDOF(unsigned pDOF, unsigned TDOF) { _pDOF = pDOF; _TDOF = TDOF; }
  virtual unsigned pDOF() const { return _pDOF; }
  virtual unsigned TDOF() const { return _TDOF; }
  virtual std::vector<unsigned> getConnectedDOFs();
  virtual void initialize(double p, double T);
  virtual void updateSolution(double p, double T);
  virtual double massTranRes(double dt) { return (_rho - _rho_o) / dt; }
  virtual double energyTranRes(double dt) { return (_rho * _e - _rho_o * _e_o) / dt; }
  virtual double massTranResBDF2(double dt) { return (1.5 * _rho - 2.0 * _rho_o + 0.5 * _rho_oo) / dt; }
  virtual double energyTranResBDF2(double dt) { return (1.5 * _rho * _e - 2.0 * _rho_o * _e_o + 0.5 * _rho_oo * _e_oo) / dt; }
  virtual double computeMassRHS(double dx);
  virtual double computeEnergyRHS(double dx, double h, double aw, double Tw);
  virtual void computeDP(double f, double dh, double gx);
  virtual void saveOldSlns();

  // Debug functions
  virtual void printConnection();

  // temp functions
  virtual void addCrossEdge(crossEdge* edge, double norm) { _cross_edges.push_back(edge); _cross_edges_in_norm.push_back(norm); }

protected:
  std::string _name;
  SinglePhaseFluid* _fluid;

  unsigned _pDOF, _TDOF;

  double _p,    _T,     _rho,     _e;
  double _p_o,  _T_o,   _rho_o,   _e_o;
  double _p_oo, _T_oo,  _rho_oo,  _e_oo;

  // second-order related
  double _p_w, _p_e, _T_w, _T_e;
  double _rho_w, _rho_e, _e_w, _e_e;

  // helper variables
  double _v_cell, _dpdx_fric, _dpdx_gravity;

  EdgeBase *WEST_EDGE, *EAST_EDGE;
  SPCell *WEST_CELL, *EAST_CELL;

  // cross flow
  std::vector<crossEdge*> _cross_edges;
  std::vector<double> _cross_edges_in_norm;
};

class EdgeBase
{
public:
  EdgeBase(std::string name, SinglePhaseFluid* fluid) : _name(name), _fluid(fluid), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~EdgeBase() {}
  virtual std::string name() final { return _name; }

  // Data access
  virtual double v()           const final { return _v; }
  virtual double mass_flux()   const final { return _mass_flux; }
  virtual double energy_flux() const final { return _energy_flux; }

  // Connection related functions
  virtual void setNghbrCells(SPCell * west, SPCell * east) final;
  virtual void setExtendedNghbrs();
  virtual SPCell * getOtherSideCell(SPCell* cell) final;

  // Residual related functions
  virtual void setDOF(unsigned vDOF) final { _vDOF = vDOF; }
  virtual unsigned vDOF() const final { return _vDOF; }
  virtual std::vector<unsigned> getConnectedDOFs() final;
  virtual void initialize(double v)     final { _v = v; _v_o = v; _v_oo = v; }
  virtual void updateSolution(double v) final { _v = v; }
  virtual void saveOldSlns()            final { _v_oo = _v_o; _v_o = _v; }
  virtual double computeTranRes(double dt) = 0;
  virtual double computeTranResBDF2(double dt) = 0;
  virtual double computeRHS(double dx) = 0;
  virtual void computeFluxes() = 0;
  virtual void computeFluxes2nd() = 0;

  // Debug functions
  virtual void printConnection() final;

protected:
  std::string _name;
  SinglePhaseFluid* _fluid;

  unsigned _vDOF;

  double _rho_edge;
  double _v, _v_o, _v_oo;
  double _mass_flux, _energy_flux;

  SPCell *WEST_CELL, *EAST_CELL;
  EdgeBase *WEST_EDGE, *EAST_EDGE;
};

class vBndryEdge : public EdgeBase
{
public:
  vBndryEdge(std::string name, SinglePhaseFluid* fluid, double v_bc, double T_bc) : EdgeBase(name, fluid), _v_bc(v_bc), _T_bc(T_bc) {}
  virtual ~vBndryEdge() {}

  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double computeTranRes(double /*dt*/) override final { return 0; }
  virtual double computeTranResBDF2(double /*dt*/) override final { return 0; }
  virtual double computeRHS(double dx) override final { return _v - _v_bc; }

protected:
  double _v_bc, _T_bc;
};

class pBndryEdge : public EdgeBase
{
public:
  pBndryEdge(std::string name, SinglePhaseFluid* fluid, double p_bc, double T_bc) : EdgeBase(name, fluid), _p_bc(p_bc), _T_bc(T_bc) {}
  virtual ~pBndryEdge() {}

  virtual void computeFluxes() override;
  virtual void computeFluxes2nd() override;
  virtual double computeTranRes(double dt) override;
  virtual double computeTranResBDF2(double dt) override;
  virtual double computeRHS(double dx) override;

protected:
  double _p_bc, _T_bc;
};

class pPseudo3DBndryEdge : public pBndryEdge
{
public:
  pPseudo3DBndryEdge(std::string name, SinglePhaseFluid* fluid, double p_bc, double T_bc) : pBndryEdge(name, fluid, p_bc, T_bc) {}
  virtual ~pPseudo3DBndryEdge() {}

  virtual void updatePBC(double p) { _p_bc = p; }
};

class IntEdge : public EdgeBase
{
public:
  IntEdge(std::string name, SinglePhaseFluid* fluid) : EdgeBase(name, fluid) {}
  ~IntEdge() {}
  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double computeTranRes(double dt) override final;
  virtual double computeTranResBDF2(double dt) override final;
  virtual double computeRHS(double dx) override;
};

class crossEdge : public IntEdge
{
public:
  crossEdge(std::string name, SinglePhaseFluid* fluid) : IntEdge(name, fluid) {}
  ~crossEdge() {}

  virtual void setNghbrEdges(EdgeBase * west, EdgeBase * east) { WEST_EDGE = west; EAST_EDGE = east; }

  virtual void setExtendedNghbrs() override final { /* do not find extended neighbors*/ }

  virtual double computeRHS(double dx) override;
};

#endif /*SINGLE_PHASE_FIELDS_H*/

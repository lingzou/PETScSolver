#ifndef SINGLE_PHASE_FIELDS_H
#define SINGLE_PHASE_FIELDS_H

#include <iostream>

class EdgeBase;

class SPCell
{
public:
  SPCell(std::string name) : _name(name) {}
  virtual ~SPCell() {}
  virtual std::string name() { return _name; }

  // Data access
  virtual double p()    const { return _p;    }
  virtual double T()    const { return _T;    }
  virtual double rho()  const { return _rho;  }
  virtual double e()    const { return _e;    }

  // Connection related functions
  virtual EdgeBase * getOtherSideEdge(EdgeBase* edge);
  virtual void setNghbrEdges(EdgeBase * west, EdgeBase * east) { WEST_EDGE = west; EAST_EDGE = east; }

  // Residual related functions
  virtual void initialize(double p, double T);
  virtual void updateSolution(double p, double T);
  virtual double massTranRes(double dt) { return (_rho - _rho_o) / dt; }
  virtual double energyTranRes(double dt) { return (_rho * _e - _rho_o * _e_o) / dt; }
  virtual double computeMassRHS(double dx);
  virtual double computeEnergyRHS(double dx);
  virtual void saveOldSlns();

  // Debug functions
  virtual void printConnection();

protected:
  double rho_func(double p, double T) { return 1.e3 + 4.e-7 * (p - 1.e5) - 0.46 * (T - 300.0); }
  double e_func(double /*p*/, double T) { return 112.55e3 + 4.e3 * (T - 300.0); }

protected:
  std::string _name;

  double _p,    _T,     _rho,     _e;
  double _p_o,  _T_o,   _rho_o,   _e_o;
  double _p_oo, _T_oo,  _rho_oo,  _e_oo;

  EdgeBase *WEST_EDGE, *EAST_EDGE;
};

class EdgeBase
{
public:
  EdgeBase(std::string name) : _name(name), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~EdgeBase() {}
  virtual std::string name() { return _name; }

  // Data access
  virtual double v()           const { return _v; }
  virtual double mass_flux()   const { return _mass_flux; }
  virtual double energy_flux() const { return _energy_flux; }

  // Connection related functions
  virtual void setNghbrCells(SPCell * west, SPCell * east);

  // Residual related functions
  virtual void initialize(double v)     { _v = v; _v_o = v; _v_oo = v; }
  virtual void updateSolution(double v) { _v = v; }
  virtual void saveOldSlns()            { _v_oo = _v_o; _v_o = _v; }
  virtual double computeTranRes(double dt) = 0;
  virtual double computeRHS(double dx) = 0;
  virtual void computeFluxes() = 0;

  // Debug functions
  virtual void printConnection();

protected:
  double rho_func(double p, double T) { return 1.e3 + 4.e-7 * (p - 1.e5) - 0.46 * (T - 300.0); }
  double e_func(double /*p*/, double T) { return 112.55e3 + 4.e3 * (T - 300.0); }

protected:
  std::string _name;

  double _rho_edge;
  double _v, _v_o, _v_oo;
  double _mass_flux, _energy_flux;

  SPCell *WEST_CELL, *EAST_CELL;
  EdgeBase *WEST_EDGE, *EAST_EDGE;
};

class vBndryEdge : public EdgeBase
{
public:
  vBndryEdge(std::string name, double v_bc, double T_bc) : EdgeBase(name), _v_bc(v_bc), _T_bc(T_bc) {}
  virtual ~vBndryEdge() {}

  virtual void computeFluxes() final;
  virtual double computeTranRes(double /*dt*/) final { return 0; }
  virtual double computeRHS(double dx) final { return _v - _v_bc; }

protected:
  double _v_bc, _T_bc;
};

class pBndryEdge : public EdgeBase
{
public:
  pBndryEdge(std::string name, double p_bc, double T_bc) : EdgeBase(name), _p_bc(p_bc), _T_bc(T_bc) {}
  virtual ~pBndryEdge() {}

  virtual void computeFluxes() final;
  virtual double computeTranRes(double dt) final;
  virtual double computeRHS(double dx) final;

protected:
  double _p_bc, _T_bc;
};

class IntEdge : public EdgeBase
{
public:
  IntEdge(std::string name) : EdgeBase(name) {}
  ~IntEdge() {}
  virtual void computeFluxes() final;
  virtual double computeTranRes(double dt) final;
  virtual double computeRHS(double dx) final;
};

#endif /*SINGLE_PHASE_FIELDS_H*/

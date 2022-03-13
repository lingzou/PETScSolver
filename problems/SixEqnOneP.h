#pragma once

#include <vector>
#include "PETScProblem.h"
#include "WaterIF97.h"
#include "SteamIF97.h"

/**** Data Structure for 6E1P problems ****/
class EdgeBase6E1P;
class Cell6E1P;
/*
class SimpleWater
{
public:
  SimpleWater() {}
  virtual ~SimpleWater() {}
  double rho(double p, double T) { return 1000.0 + 1.0e-7 * (p - 1.0e5) - 1.0e-7 * (T - 450); }
  double e(double p, double T)   { return 4200 * T; }
};

class SimpleSteam
{
public:
  SimpleSteam() {}
  virtual ~SimpleSteam() {}
  double rho(double p, double T) { return 0.5 + 1.0e-5 * (p - 1.0e5) - 1.0e-7 * (T - 450); }
  double e(double p, double T)   { return 1000 * T; }
};*/

class Cell6E1P
{
public:
  Cell6E1P(std::string name) : _name(name), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~Cell6E1P() {}
  virtual std::string name() { return _name; }

  // Data access
  virtual double p()      const { return _p;      }
  virtual double T_l()    const { return _T_l;    }
  virtual double T_g()    const { return _T_g;    }
  virtual double alpha()  const { return _alpha;  }
  virtual double rho_l()  const { return _rho_l;    }
  virtual double rho_g()  const { return _rho_g;    }
  virtual double e_l()    const { return _e_l;    }
  virtual double e_g()    const { return _e_g;    }
  virtual double g()      const { return _g;        }
  void set_g(double g) { _g = g; }
  //   2nd-order related
  virtual double p_w()      const { return _p_w;      }
  virtual double p_e()      const { return _p_e;      }
  virtual double T_l_w()    const { return _T_l_w;    }
  virtual double T_l_e()    const { return _T_l_e;    }
  virtual double T_g_w()    const { return _T_g_w;    }
  virtual double T_g_e()    const { return _T_g_e;    }
  virtual double alpha_w()  const { return _alpha_w;  }
  virtual double alpha_e()  const { return _alpha_e;  }
  virtual double rho_l_w()  const { return _rho_l_w;  }
  virtual double rho_l_e()  const { return _rho_l_e;  }
  virtual double rho_g_w()  const { return _rho_g_w;  }
  virtual double rho_g_e()  const { return _rho_g_e;  }
  virtual double e_l_w()    const { return _e_l_w;  }
  virtual double e_l_e()    const { return _e_l_e;  }
  virtual double e_g_w()    const { return _e_g_w;  }
  virtual double e_g_e()    const { return _e_g_e;  }
  virtual void linearReconstruction(double, double, double, double, double, double, double, double);

  // Connection related functions
  virtual EdgeBase6E1P * getOtherSideEdge(EdgeBase6E1P* edge);
  virtual void setNghbrEdges(EdgeBase6E1P * west, EdgeBase6E1P * east) { WEST_EDGE = west; EAST_EDGE = east; }
  virtual void setExtendedNghbrs();

  virtual Cell6E1P * west_cell() final { return WEST_CELL; }
  virtual Cell6E1P * east_cell() final { return EAST_CELL; }

  // Residual related functions
  virtual void setDOF(unsigned aDOF, unsigned pDOF, unsigned TlDOF, unsigned TgDOF) { _aDOF = aDOF; _pDOF = pDOF; _TlDOF = TlDOF; _TgDOF = TgDOF; }
  virtual unsigned aDOF()   const { return _aDOF; }
  virtual unsigned pDOF()   const { return _pDOF; }
  virtual unsigned TlDOF()  const { return _TlDOF; }
  virtual unsigned TgDOF()  const { return _TgDOF; }
  virtual std::vector<unsigned> getConnectedDOFs();
  virtual void initialize(double alpha, double p, double Tl, double Tg);
  virtual void updateSolution(double alpha, double p, double Tl, double Tg);
  virtual double liquidMassTranRes(double dt) { return ((1-_alpha) * _rho_l - (1-_alpha_o) * _rho_l_o) / dt; }
  virtual double gasMassTranRes(double dt) { return (_alpha * _rho_g - _alpha_o * _rho_g_o) / dt; }
  virtual double liquidEnergyTranRes(double dt);
  virtual double gasEnergyTranRes(double dt);
  virtual double liquidMassTranResBDF2(double dt) { return (1.5 * (1-_alpha) * _rho_l - 2.0 * (1-_alpha_o) * _rho_l_o + 0.5 * (1-_alpha_oo) * _rho_l_oo) / dt; }
  virtual double gasMassTranResBDF2(double dt) { return (1.5 * _alpha * _rho_g - 2.0 * _alpha_o * _rho_g_o + 0.5 * _alpha_oo * _rho_g_oo) / dt; }
  virtual double liquidEnergyTranResBDF2(double dt);
  virtual double gasEnergyTranResBDF2(double dt);
  virtual double liquidMassRHS(double dx);
  virtual double gasMassRHS(double dx);
  virtual double liquidEnergyRHS(double dx);
  virtual double gasEnergyRHS(double dx);
  virtual void saveOldSlns();
  virtual void setEOS(Water_IF97* water, Steam_IF97 *steam) { _water = water; _steam = steam; }

  // Debug functions
  virtual void printConnection();

protected:
  std::string _name;
  Water_IF97 * _water;
  Steam_IF97 * _steam;

  double _g;

  unsigned _aDOF, _pDOF, _TlDOF, _TgDOF;

  double _alpha,    _p,     _T_l,     _T_g,     _rho_l,     _rho_g,     _e_l,     _e_g;
  double _alpha_o,  _p_o,   _T_l_o,   _T_g_o,   _rho_l_o,   _rho_g_o,   _e_l_o,   _e_g_o;
  double _alpha_oo, _p_oo,  _T_l_oo,  _T_g_oo,  _rho_l_oo,  _rho_g_oo,  _e_l_oo,  _e_g_oo;

  // second-order related
  double _alpha_w, _alpha_e, _p_w, _p_e, _T_l_w, _T_l_e, _T_g_w, _T_g_e;
  double _rho_l_w, _rho_l_e, _rho_g_w, _rho_g_e;
  double _e_l_w, _e_l_e, _e_g_w, _e_g_e;

  EdgeBase6E1P *WEST_EDGE, *EAST_EDGE;
  Cell6E1P *WEST_CELL, *EAST_CELL;
};

class EdgeBase6E1P
{
public:
  EdgeBase6E1P(std::string name) : _name(name), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~EdgeBase6E1P() {}
  virtual std::string name() final { return _name; }

  // Data access
  virtual double v_l()          const final { return _vl; }
  virtual double v_g()          const final { return _vg; }
  virtual double vl_w()         const final { return _vl_w;  }
  virtual double vl_e()         const final { return _vl_e;  }
  virtual double vg_w()         const final { return _vg_w;  }
  virtual double vg_e()         const final { return _vg_e;  }
  virtual double mass_flux_l()  const final { return _mass_flux_l; }
  virtual double mass_flux_g()  const final { return _mass_flux_g; }
  virtual double energy_flux_l()  const final { return _energy_flux_l; }
  virtual double energy_flux_g()  const final { return _energy_flux_g; }
  virtual double alphal_ul()    const final { return _alphal_ul; }
  virtual double alphag_ug()    const final { return _alphag_ug; }

  virtual void linearReconstruction(double, double, double, double);

  // Connection related functions
  virtual void setNghbrCells(Cell6E1P * west, Cell6E1P * east) final;
  virtual void setExtendedNghbrs() final ;
  virtual Cell6E1P * getOtherSideCell(Cell6E1P* cell) final;

  virtual EdgeBase6E1P * west_edge() final { return WEST_EDGE; }
  virtual EdgeBase6E1P * east_edge() final { return EAST_EDGE; }

  // Residual related functions
  virtual void setDOF(unsigned vlDOF, unsigned vgDOF) final { _vlDOF = vlDOF; _vgDOF = vgDOF; }
  virtual unsigned vlDOF() const final { return _vlDOF; }
  virtual unsigned vgDOF() const final { return _vgDOF; }
  virtual std::vector<unsigned> getConnectedDOFs() final;
  virtual void initialize(double vl, double vg) final
  { _vl = vl; _vl_o = vl; _vl_oo = vl;
    _vg = vg; _vg_o = vg; _vg_oo = vg; }
  virtual void updateSolution(double vl, double vg) final { _vl = vl; _vg = vg; }
  virtual void saveOldSlns() final
  { _vl_oo = _vl_o; _vl_o = _vl;   _vg_oo = _vg_o; _vg_o = _vg;}
  virtual double liquidVelTranRes(double dt) = 0;
  virtual double gasVelTranRes(double dt) = 0;
  virtual double liquidVelTranResBDF2(double dt) = 0;
  virtual double gasVelTranResBDF2(double dt) = 0;
  virtual void computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g) = 0;
  virtual void computeFluxes() = 0;
  virtual void computeFluxes2nd() = 0;
  virtual double interfacial_drag(double alpha_edge, double rho_l_edge, double rho_g_edge) final;

  // Debug functions
  virtual void printConnection() final;

  virtual void setEOS(Water_IF97* water, Steam_IF97 *steam) { _water = water; _steam = steam; }

protected:
  std::string _name;
  Water_IF97 * _water;
  Steam_IF97 * _steam;

  unsigned _vlDOF, _vgDOF;

  double _rho_l_edge, _rho_g_edge;
  double _vl, _vl_o, _vl_oo;
  double _vg, _vg_o, _vg_oo;
  double _mass_flux_l, _mass_flux_g, _energy_flux_l, _energy_flux_g;
  double _alphal_ul, _alphag_ug;

  // second-order related
  double _vl_w, _vl_e, _vg_w, _vg_e;

  Cell6E1P *WEST_CELL, *EAST_CELL;
  EdgeBase6E1P *WEST_EDGE, *EAST_EDGE;
};

class vBndryEdge6E1P : public EdgeBase6E1P
{
public:
  vBndryEdge6E1P(std::string name, double vl_bc, double vg_bc, double alpha_bc, double T_l_bc, double T_g_bc) :
    EdgeBase6E1P(name), _vl_bc(vl_bc), _vg_bc(vg_bc), _alpha_bc(alpha_bc), _T_l_bc(T_l_bc), _T_g_bc(T_g_bc) {}
  virtual ~vBndryEdge6E1P() {}

  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double /*dt*/) override final { return 0; }
  virtual double liquidVelTranResBDF2(double /*dt*/) override final { return 0; }
  virtual double gasVelTranRes(double /*dt*/) override final { return 0; }
  virtual double gasVelTranResBDF2(double /*dt*/) override final { return 0; }
  virtual void computeRHS(bool /*compute_int_drag*/, unsigned /*order*/, double /*dx*/, double & rhs_l, double & rhs_g) override final
  { rhs_l = _vl_bc - _vl; rhs_g = _vg_bc - _vg; } // rhs will flip sign

protected:
  double _vl_bc, _vg_bc, _alpha_bc, _T_l_bc, _T_g_bc;
};

class pBndryEdge6E1P : public EdgeBase6E1P
{
public:
  pBndryEdge6E1P(std::string name, double p_bc, double alpha_bc, double T_l_bc, double T_g_bc) :
    EdgeBase6E1P(name), _p_bc(p_bc), _alpha_bc(alpha_bc), _T_l_bc(T_l_bc), _T_g_bc(T_g_bc) {}
  virtual ~pBndryEdge6E1P() {}

  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double dt) override final;
  virtual double liquidVelTranResBDF2(double dt) override final;
  virtual double gasVelTranRes(double dt) override final;
  virtual double gasVelTranResBDF2(double dt) override final;

protected:
  double _p_bc, _alpha_bc, _T_l_bc, _T_g_bc;
};

class pWESTBndryEdge6E1P : public pBndryEdge6E1P
{
public:
  pWESTBndryEdge6E1P(std::string name, double p_bc, double alpha_bc, double T_l_bc, double T_g_bc) :
    pBndryEdge6E1P(name, p_bc, alpha_bc, T_l_bc, T_g_bc) {}
  virtual ~pWESTBndryEdge6E1P() {}

  virtual void computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

class pEASTBndryEdge6E1P : public pBndryEdge6E1P
{
public:
  pEASTBndryEdge6E1P(std::string name, double p_bc, double alpha_bc, double T_l_bc, double T_g_bc) :
    pBndryEdge6E1P(name, p_bc, alpha_bc, T_l_bc, T_g_bc) {}
  virtual ~pEASTBndryEdge6E1P() {}

  virtual void computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

class IntEdge6E1P : public EdgeBase6E1P
{
public:
  IntEdge6E1P(std::string name) : EdgeBase6E1P(name) {}
  ~IntEdge6E1P() {}
  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double dt) override final;
  virtual double liquidVelTranResBDF2(double dt) override final;
  virtual double gasVelTranRes(double dt) override final;
  virtual double gasVelTranResBDF2(double dt) override final;
  virtual void computeRHS(bool compute_int_drag, unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

/**** End of Data Structure for 6E1P problems ****/

class SixEqnOneP : public PETScProblem
{
public:
  SixEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  ~SixEqnOneP();

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE * file) override final;
  virtual void writeTextOutput(FILE * file) override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;

  virtual void onLastTimestepEnd() override final;

  virtual void setupExtendedConnections() override final;
  virtual void updateEdgeCellHelperVar() override final;
  virtual void setDOFoffset(unsigned) override final;
  virtual void linearReconstruction() override final;

  enum ProblemType
  {
    SINE_WAVE       = 1,      // A gentel advection problem (sine wave)
    SQUARE_WAVE     = 2,      // Advection problem with extremely small alpha_l and alpha_g values, 1.e-6
    MANOMETER       = 3,      // Manometer problem
    SEDIMENTATION   = 4,      // Sedimentation problem
    WATER_FAUCET    = 5
  };

protected:
  double gx_vol(unsigned);

protected:
  int _problem_type;
  bool _compute_int_drag;

  unsigned int _order;
  double length;
  double dx;
  unsigned int n_Cell, n_Node;

  // Post-processing stuff, e.g., in Manometer problem, v_l_bottom vs. time
  std::vector<double> time;
  std::vector<double> v_l_bottom, p_bottom, h_left, h_right;

  //
  std::vector<EdgeBase6E1P*> _edges;
  std::vector<Cell6E1P*> _cells;

  Water_IF97 * _water;
  Steam_IF97 * _steam;
};

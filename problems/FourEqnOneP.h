#ifndef FOUR_EQN_ONEP_H
#define FOUR_EQN_ONEP_H

#include <vector>
#include "PETScProblem.h"

//double rho_l_func(double p) { return 1000.0 + 1.0e-7 * (p - 1.0e5); }
//double rho_g_func(double p) { return    0.5 + 1.0e-5 * (p - 1.0e5); }

/**** Data Structure for 4E1P problems ****/
class EdgeBase4E1P;

class Cell4E1P
{
public:
  Cell4E1P(std::string name) : _name(name), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~Cell4E1P() {}
  virtual std::string name() { return _name; }

  // Data access
  virtual double p()      const { return _p;      }
  virtual double alpha()  const { return _alpha;  }
  virtual double rho_l()  const { return _rho_l;    }
  virtual double rho_g()  const { return _rho_g;    }
  virtual double g()      const { return _g;        }
  void set_g(double g) { _g = g; }
  //virtual double dpdx_fric()      const { return _dpdx_fric; }
  //virtual double dpdx_gravity()   const { return _dpdx_gravity; }
  //   2nd-order related
  virtual double p_w()      const { return _p_w;      }
  virtual double p_e()      const { return _p_e;      }
  virtual double alpha_w()  const { return _alpha_w;  }
  virtual double alpha_e()  const { return _alpha_e;  }
  virtual double rho_l_w()  const { return _rho_l_w;  }
  virtual double rho_l_e()  const { return _rho_l_e;  }
  virtual double rho_g_w()  const { return _rho_g_w;  }
  virtual double rho_g_e()  const { return _rho_g_e;  }
  virtual void linearReconstruction(double, double, double, double);

  // Connection related functions
  virtual EdgeBase4E1P * getOtherSideEdge(EdgeBase4E1P* edge);
  virtual void setNghbrEdges(EdgeBase4E1P * west, EdgeBase4E1P * east) { WEST_EDGE = west; EAST_EDGE = east; }
  virtual void setExtendedNghbrs();

  // Residual related functions
  virtual void setDOF(unsigned aDOF, unsigned pDOF) { _aDOF = aDOF; _pDOF = pDOF; }
  virtual unsigned aDOF() const { return _aDOF; }
  virtual unsigned pDOF() const { return _pDOF; }
  virtual std::vector<unsigned> getConnectedDOFs();
  virtual void initialize(double alpha, double p);
  virtual void updateSolution(double alpha, double p);
  virtual double liquidMassTranRes(double dt) { return ((1-_alpha) * _rho_l - (1-_alpha_o) * _rho_l_o) / dt; }
  virtual double gasMassTranRes(double dt) { return (_alpha * _rho_g - _alpha_o * _rho_g_o) / dt; }
  virtual double liquidMassTranResBDF2(double dt) { return (1.5 * (1-_alpha) * _rho_l - 2.0 * (1-_alpha_o) * _rho_l_o + 0.5 * (1-_alpha_oo) * _rho_l_oo) / dt; }
  virtual double gasMassTranResBDF2(double dt) { return (1.5 * _alpha * _rho_g - 2.0 * _alpha_o * _rho_g_o + 0.5 * _alpha_oo * _rho_g_oo) / dt; }
  virtual double liquidMassRHS(double dx);
  virtual double gasMassRHS(double dx);
  //virtual void computeDP();
  virtual void saveOldSlns();

  // Debug functions
  virtual void printConnection();

protected:
  std::string _name;
  //SinglePhaseFluid* _fluid;

  double _g;

  unsigned _aDOF, _pDOF;

  double _alpha,    _p,     _rho_l,     _rho_g;
  double _alpha_o,  _p_o,   _rho_l_o,   _rho_g_o;
  double _alpha_oo, _p_oo,  _rho_l_oo,  _rho_g_oo;

  // second-order related
  double _alpha_w, _alpha_e, _p_w, _p_e;
  double _rho_l_w, _rho_l_e, _rho_g_w, _rho_g_e;

  // helper variables
  //double _dpdx_fric, _dpdx_gravity;

  EdgeBase4E1P *WEST_EDGE, *EAST_EDGE;
  Cell4E1P *WEST_CELL, *EAST_CELL;
};

class EdgeBase4E1P
{
public:
  EdgeBase4E1P(std::string name) : _name(name), WEST_CELL(NULL), EAST_CELL(NULL), WEST_EDGE(NULL), EAST_EDGE(NULL) {}
  virtual ~EdgeBase4E1P() {}
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

  virtual void linearReconstruction(double, double, double, double);

  // Connection related functions
  virtual void setNghbrCells(Cell4E1P * west, Cell4E1P * east) final;
  virtual void setExtendedNghbrs() final ;
  virtual Cell4E1P * getOtherSideCell(Cell4E1P* cell) final;

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
  //virtual double computeLiquidRHS(unsigned order, double dx) = 0;
  //virtual double computeGasRHS(unsigned order, double dx) = 0;
  virtual void computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g) = 0;
  virtual void computeFluxes() = 0;
  virtual void computeFluxes2nd() = 0;

  // Debug functions
  virtual void printConnection() final;

protected:
  std::string _name;
  //SinglePhaseFluid* _fluid;

  unsigned _vlDOF, _vgDOF;

  double _rho_l_edge, _rho_g_edge;
  double _vl, _vl_o, _vl_oo;
  double _vg, _vg_o, _vg_oo;
  double _mass_flux_l, _mass_flux_g;

  // second-order related
  double _vl_w, _vl_e, _vg_w, _vg_e;

  Cell4E1P *WEST_CELL, *EAST_CELL;
  EdgeBase4E1P *WEST_EDGE, *EAST_EDGE;
};

class vBndryEdge4E1P : public EdgeBase4E1P
{
public:
  vBndryEdge4E1P(std::string name, double vl_bc, double vg_bc, double alpha_bc) : EdgeBase4E1P(name), _vl_bc(vl_bc), _vg_bc(vg_bc), _alpha_bc(alpha_bc) {}
  virtual ~vBndryEdge4E1P() {}

  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double /*dt*/) override final { return 0; }
  virtual double liquidVelTranResBDF2(double /*dt*/) override final { return 0; }
  virtual double gasVelTranRes(double /*dt*/) override final { return 0; }
  virtual double gasVelTranResBDF2(double /*dt*/) override final { return 0; }
  //virtual double computeLiquidRHS(unsigned /*order*/, double /*dx*/) override final { return _vl - _vl_bc; }
  //virtual double computeGasRHS(unsigned /*order*/, double /*dx*/) override final { return _vg - _vg_bc; }
  virtual void computeRHS(unsigned /*order*/, double /*dx*/, double & rhs_l, double & rhs_g) override final
  { rhs_l = _vl - _vl_bc; rhs_g = _vg - _vg_bc; }

protected:
  double _vl_bc, _vg_bc, _alpha_bc;
};

class pBndryEdge4E1P : public EdgeBase4E1P
{
public:
  pBndryEdge4E1P(std::string name, double p_bc, double alpha_bc) : EdgeBase4E1P(name), _p_bc(p_bc), _alpha_bc(alpha_bc) {}
  virtual ~pBndryEdge4E1P() {}

  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double dt) override final;
  virtual double liquidVelTranResBDF2(double dt) override final;
  virtual double gasVelTranRes(double dt) override final;
  virtual double gasVelTranResBDF2(double dt) override final;
  //virtual double computeLiquidRHS(unsigned order, double dx) = 0;
  //virtual double computeGasRHS(unsigned order, double dx) = 0;

protected:
  double _p_bc, _alpha_bc;
};

class pWESTBndryEdge4E1P : public pBndryEdge4E1P
{
public:
  pWESTBndryEdge4E1P(std::string name, double p_bc, double alpha_bc) : pBndryEdge4E1P(name, p_bc, alpha_bc) {}
  virtual ~pWESTBndryEdge4E1P() {}

  //virtual double computeLiquidRHS(unsigned order, double dx) override final;
  //virtual double computeGasRHS(unsigned order, double dx) override final;
  virtual void computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

class pEASTBndryEdge4E1P : public pBndryEdge4E1P
{
public:
  pEASTBndryEdge4E1P(std::string name, double p_bc, double alpha_bc) : pBndryEdge4E1P(name, p_bc, alpha_bc) {}
  virtual ~pEASTBndryEdge4E1P() {}

  //virtual double computeLiquidRHS(unsigned order, double dx) override final;
  //virtual double computeGasRHS(unsigned order, double dx) override final;
  virtual void computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

class IntEdge4E1P : public EdgeBase4E1P
{
public:
  IntEdge4E1P(std::string name) : EdgeBase4E1P(name) {}
  ~IntEdge4E1P() {}
  virtual void computeFluxes() override final;
  virtual void computeFluxes2nd() override final;
  virtual double liquidVelTranRes(double dt) override final;
  virtual double liquidVelTranResBDF2(double dt) override final;
  virtual double gasVelTranRes(double dt) override final;
  virtual double gasVelTranResBDF2(double dt) override final;
  //virtual double computeLiquidRHS(unsigned order, double dx) override final;
  //virtual double computeGasRHS(unsigned order, double dx) override final;
  virtual void computeRHS(unsigned order, double dx, double & rhs_l, double & rhs_g) override final;
};

/**** End of Data Structure for 4E1P problems ****/

class FourEqnOneP : public PETScProblem
{
public:
  FourEqnOneP(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem);
  ~FourEqnOneP();

  virtual void SetupInitialCondition(double * u) override final;
  virtual void updateSolution(double *u) override final;

  virtual void transientResidual(double * res) override final;
  virtual void RHS(double * rhs) override final;
  virtual void onTimestepEnd() override final;
  virtual void writeVTKOutput(FILE * file) override final;
  virtual void writeTextOutput(FILE * file) override final;

  virtual void FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp) override final;

  void RHS_1st_order(double * rhs);
  void RHS_2nd_order(double * rhs);

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
  //double rho_l_func(double p) { return 1000.0 + 1.0e-7 * (p - 1.0e5); }
  //double rho_g_func(double p) { return    0.5 + 1.0e-5 * (p - 1.0e5); }
  double gx_vol(unsigned);

protected:
  int _problem_type;
  bool _periodic_bc;
  double P_INLET, P_OUTLET, ALPHA_INLET, ALPHA_OUTLET, V_L_INLET, V_G_INLET;

  unsigned int _order;
  double g;

  double length;
  double dx;
  unsigned int n_Cell, n_Node;
  double ALPHA_MIN;

  // State variables
  // 1) Primary variables
  std::vector<double> alpha,     p,     v_l,     v_g;
  std::vector<double> alpha_old, p_old, v_l_old, v_g_old;
  std::vector<double> alpha_oo,  p_oo,  v_l_oo,  v_g_oo;
  // 2) Dependent variables (rho appears in d/dt terms, needs appropriately initialized)
  std::vector<double> rho_l, rho_g, rho_l_old, rho_g_old, rho_l_oo, rho_g_oo;

  // Helper variables
  std::vector<double> alpha_edge, rho_l_edge, rho_g_edge;
  std::vector<double> v_l_cell, v_g_cell;

  // Fluxes
  std::vector<double> rho_l_flux, rho_g_flux;

  // Second-order helper variables
  std::vector<double> alpha_w, alpha_e, p_w, p_e;
  std::vector<double> v_l_w, v_l_e, v_g_w, v_g_e;

  // Post-processing stuff, e.g., in Manometer problem, v_l_bottom vs. time
  std::vector<double> time;
  std::vector<double> v_l_bottom, p_bottom;

  //
  std::vector<EdgeBase4E1P*> _edges;
  std::vector<Cell4E1P*> _cells;
};

#endif /*FOUR_EQN_ONEP_H*/

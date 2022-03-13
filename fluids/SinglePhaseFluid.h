#pragma once

#include "ParameterList.h"

class SinglePhaseFluid
{
public:
  SinglePhaseFluid(InputParameterList & inputParamList);
  virtual ~SinglePhaseFluid() {}

  virtual double rho(double p, double T) const = 0;
  virtual double e(double p, double T) const = 0;

protected:
  InputParameterList & _inputParamList;
};

class linearFluid : public SinglePhaseFluid
{
public:
  linearFluid(InputParameterList & inputParamList);
  virtual ~linearFluid() {}

  virtual double rho(double p, double T)   const override final { return _rho0 + (p - _p0) * _drho_dp + (T - _T0) * _drho_dT; }
  virtual double e(double /*p*/, double T) const override final { return _e0 + _cv * T; }

protected:
  double _rho0, _p0, _T0, _e0;
  double _drho_dp, _drho_dT, _cv;
};

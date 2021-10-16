#include "SinglePhaseFluid.h"

SinglePhaseFluid::SinglePhaseFluid(InputParameterList & inputParamList) :
  _inputParamList(inputParamList)
{}

linearFluid::linearFluid(InputParameterList & inputParamList) :
  SinglePhaseFluid(inputParamList)
{
  _inputParamList.readRequiredInputParameter<double>("rho0");
  _inputParamList.readRequiredInputParameter<double>("p0");
  _inputParamList.readRequiredInputParameter<double>("T0");
  _inputParamList.readRequiredInputParameter<double>("e0");
  _inputParamList.readRequiredInputParameter<double>("drho_dp");
  _inputParamList.readRequiredInputParameter<double>("drho_dT");
  _inputParamList.readRequiredInputParameter<double>("cv");

  _rho0 = _inputParamList.getParameterValue<double>("rho0");
  _p0 = _inputParamList.getParameterValue<double>("p0");
  _T0 = _inputParamList.getParameterValue<double>("T0");
  _e0 = _inputParamList.getParameterValue<double>("e0");
  _drho_dp = _inputParamList.getParameterValue<double>("drho_dp");
  _drho_dT = _inputParamList.getParameterValue<double>("drho_dT");
  _cv = _inputParamList.getParameterValue<double>("cv");
}

#pragma once

#include "SinglePhaseFluid.h"

// Based on data generated from IF97 library.
// https://github.com/lingzou/IAPWS-1995-IF97
// Region 3 metastable data were genereated from extrapolation of stable region data,
// so they are not accurate and may even cause convergence issue.
// This needs future improvement.

class Water_IF97 : public SinglePhaseFluid
{
public:
  Water_IF97(InputParameterList & inputParamList);
  ~Water_IF97() {}

  virtual double rho(double p, double T) const override final;
  virtual double e(double p, double T) const override final;
  double k(double p, double T) const;
  double mu(double p, double T) const;
  double cp(double p, double T) const;

  void computeProperties(double p, double T, double & rho, double & e, double & k, double & mu, double &cp, double &beta);

protected:
  int find_T_lower_bound(double T) const;
  int find_p_lower_bound(double p) const;

protected:
  // 1k, 3k, 5k, 7k;              # 4     data pts
  // 10k -> 100k every 5k         # 19    data pts
  // 200k -> 22000k every 100k    # 219   data pts
  float p_array[242];

  // 274 -> 647 every 1Kelvin
  float T_array[374];
  float rho_matrix[242][374];   // matrix[p_idx][T_idx]
  float e_matrix[242][374];
  float k_matrix[242][374];
  float mu_matrix[242][374];
  float cp_matrix[242][374];
};

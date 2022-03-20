#include <stdio.h>
#include <math.h>
#include <petscsnes.h>

#include "WaterIF97.h"
#include "utils.h"

Water_IF97::Water_IF97(InputParameterList & inputParamList) :
  SinglePhaseFluid(inputParamList)
{
  FILE * pFile_rho  = fopen("fluids/rho_l.txt", "r");
  FILE * pFile_e    = fopen("fluids/e_l.txt", "r");
  FILE * pFile_k    = fopen("fluids/k_l.txt", "r");
  FILE * pFile_mu   = fopen("fluids/mu_l.txt", "r");
  FILE * pFile_cp   = fopen("fluids/cp_l.txt", "r");

  if(pFile_rho == NULL)       sysError("Cannot open fluids/rho_l.txt.");
  if(pFile_e == NULL)         sysError("Cannot open fluids/e_l.txt.");
  if(pFile_k == NULL)         sysError("Cannot open fluids/k_l.txt.");
  if(pFile_mu == NULL)        sysError("Cannot open fluids/mu_l.txt.");
  if(pFile_cp == NULL)        sysError("Cannot open fluids/cp_l.txt.");

  // Read in the first line, temperatures (just once)
  float dummy;
  for(int i = 0; i < 374; i++)
  {
    fscanf(pFile_rho, "%f", &T_array[i]);
    fscanf(pFile_e,   "%f", &dummy);
    fscanf(pFile_k,   "%f", &dummy);
    fscanf(pFile_mu,  "%f", &dummy);
    fscanf(pFile_cp,  "%f", &dummy);
  }

  // Read in data line by line. 1st data is pressure, the remaining numbers are rho (e, etc.)
  for(int i = 0; i < 242; i++) // loop on p
  {
    fscanf(pFile_rho, "%f", &p_array[i]);
    fscanf(pFile_e,   "%f", &dummy);
    fscanf(pFile_k,   "%f", &dummy);
    fscanf(pFile_mu,  "%f", &dummy);
    fscanf(pFile_cp,  "%f", &dummy);

    for(int j = 0; j < 374; j++) // loop on T
    {
      fscanf(pFile_rho, "%f", &rho_matrix[i][j]);
      fscanf(pFile_e,   "%f", &e_matrix[i][j]);
      fscanf(pFile_k,   "%f", &k_matrix[i][j]);
      fscanf(pFile_mu,  "%f", &mu_matrix[i][j]);
      fscanf(pFile_cp,  "%f", &cp_matrix[i][j]);
    }
  }

  fclose (pFile_rho);
  fclose (pFile_e);
  fclose (pFile_k);
  fclose (pFile_mu);
  fclose (pFile_cp);
}

double
Water_IF97::rho(double p, double T) const
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  return UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], rho_matrix[i][j], rho_matrix[i][j+1], rho_matrix[i+1][j], rho_matrix[i+1][j+1]);
}

double
Water_IF97::e(double p, double T) const
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  return UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], e_matrix[i][j], e_matrix[i][j+1], e_matrix[i+1][j], e_matrix[i+1][j+1]);
}

double
Water_IF97::k(double p, double T) const
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  return UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], k_matrix[i][j], k_matrix[i][j+1], k_matrix[i+1][j], k_matrix[i+1][j+1]);
}

double
Water_IF97::mu(double p, double T) const
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  return UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], mu_matrix[i][j], mu_matrix[i][j+1], mu_matrix[i+1][j], mu_matrix[i+1][j+1]);
}

double
Water_IF97::cp(double p, double T) const
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  return UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], cp_matrix[i][j], cp_matrix[i][j+1], cp_matrix[i+1][j], cp_matrix[i+1][j+1]);
}

void
Water_IF97::computeProperties(double p, double T, double & rho, double & e, double & k, double & mu, double & cp, double & beta)
{
  int i = find_p_lower_bound(p);
  int j = find_T_lower_bound(T);

  rho = UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1], rho_matrix[i][j], rho_matrix[i][j+1], rho_matrix[i+1][j], rho_matrix[i+1][j+1]);
  e   = UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1],   e_matrix[i][j],   e_matrix[i][j+1],   e_matrix[i+1][j],   e_matrix[i+1][j+1]);
  k   = UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1],   k_matrix[i][j],   k_matrix[i][j+1],   k_matrix[i+1][j],   k_matrix[i+1][j+1]);
  mu  = UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1],  mu_matrix[i][j],  mu_matrix[i][j+1],  mu_matrix[i+1][j],  mu_matrix[i+1][j+1]);
  cp  = UTILS::bilinear(p, T, p_array[i], p_array[i+1], T_array[j], T_array[j+1],  cp_matrix[i][j],  cp_matrix[i][j+1],  cp_matrix[i+1][j],  cp_matrix[i+1][j+1]);

  double dp = p_array[i+1] - p_array[i];
  double dT = T_array[i+1] - T_array[i];
  double dp_1 = p_array[i+1] - p;
  double dp_0 = p - p_array[i];
  double d_rho_0 = rho_matrix[i][j+1] - rho_matrix[i][j];     // 01 - 00
  double d_rho_1 = rho_matrix[i+1][j+1] - rho_matrix[i+1][j]; // 11 - 10

  double d_rho_dT = (dp_1 * d_rho_0 + dp_0 * d_rho_1) / (dp * dT);
  beta = - d_rho_dT / rho;
}

int
Water_IF97::find_T_lower_bound(double T) const
{
  if (T < T_array[0])           SETERRQ(PETSC_COMM_WORLD, 103, "103"); //throw std::runtime_error("103"); //sysException("T out of lower bound." + std::to_string(T), 103);
  else if (T < T_array[373])    return floor(T - T_array[0]);
  else                          SETERRQ(PETSC_COMM_WORLD, 104, "104"); //throw std::runtime_error("104"); //ysException("T out of upper bound." + std::to_string(T), 104);

  return -1;
}

int
Water_IF97::find_p_lower_bound(double p) const
{
  // 1k, 3k, 5k, 7k;  (4)
  // 10k -> 95k every 5k (18)
  // 100k -> 22000k every 100k (220)
  if (p < p_array[0])
    //sysError("p out of lower bound.");
    SETERRQ(PETSC_COMM_WORLD, 203, "203");
  else if (p < 3e3)
    return 0;
  else if (p < 5e3)
    return 1;
  else if (p < 7e3)
    return 2;
  else if (p < 10e3)
    return 3;
  else if (p < 100e3)
    return floor((p - 10.e3) / 5.e3) + 4;
  else if (p < 22e6)
    return floor((p - 100.e3) / 100.e3) + 22;
  else
    SETERRQ(PETSC_COMM_WORLD, 204, "204");
    //sysError("p out of upper bound.");

  return -1;
}

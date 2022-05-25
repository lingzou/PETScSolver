#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "EulerEquation1D.h"
#include "utils.h"

EulerEquation1D::EulerEquation1D(InputParameterList & globalParamList, InputParameterList & inputParamList, ProblemSystem * problemSystem) :
  PETScProblem(globalParamList, inputParamList, problemSystem)
{
  _order = _inputParamList.getValueFromInput<int>("order");
  _gamma = 1.4;

  length = 1.0;
  n_Cell = _inputParamList.getValueFromInput<int>("n_cells");
  n_Node = n_Cell + 1;
  _n_DOFs = n_Cell * 3;

  dx = length / n_Cell;

  int p_case = _inputParamList.getValueFromInput<int>("p_case");
  if (p_case == 1)  // Sod problem
  {
    RHO_L = 1.0;      M_L = 0.0;    E_L = 2.5;
    RHO_R = 0.125;    M_R = 0.0;    E_R = 0.25;
  }
  else if (p_case == 2)  // Lax problem
  {
    RHO_L = 0.445;    M_L = 0.311;    E_L = 8.928;
    RHO_R = 0.500;    M_R = 0.000;    E_R = 1.4275;
  }
  else
    sysError("case = 1 for Sod problem; case = 2 for Lax problem");

  P_L = p_IG(RHO_L, M_L, E_L);
  P_R = p_IG(RHO_R, M_R, E_R);

  rho.resize(n_Cell);       m.resize(n_Cell);       E.resize(n_Cell);     p.resize(n_Cell);
  rho_old.resize(n_Cell);   m_old.resize(n_Cell);   E_old.resize(n_Cell);
  rho_oo.resize(n_Cell);    m_oo.resize(n_Cell);    E_oo.resize(n_Cell);

  flux_rho.resize(n_Cell + 1); flux_m.resize(n_Cell + 1); flux_E.resize(n_Cell + 1);

  if (_order > 1)
  {
    rho_w.resize(n_Cell);       m_w.resize(n_Cell);       E_w.resize(n_Cell);
    rho_e.resize(n_Cell);       m_e.resize(n_Cell);       E_e.resize(n_Cell);
  }
}

EulerEquation1D::~EulerEquation1D() {}

void
EulerEquation1D::SetupInitialCondition(double * u)
{
  unsigned index = 0;

  for(unsigned i = 0; i < n_Cell; i++)
  {
    rho[i] = (i < n_Cell/2) ? RHO_L : RHO_R;
    m[i]   = (i < n_Cell/2) ? M_L   : M_R;
    E[i]   = (i < n_Cell/2) ? E_L   : E_R;
    p[i]   = (i < n_Cell/2) ? P_L   : P_R;

    u[index++] = rho[i];
    u[index++] = m[i];
    u[index++] = E[i];
  }
  rho_old = rho;  m_old = m; E_old = E;
  rho_oo = rho;   m_oo = m;  E_oo = E;
}

void
EulerEquation1D::updateSolution(double * u)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < n_Cell; i++)
  {
    rho[i] = u[idx++];  m[i] = u[idx++];  E[i] = u[idx++];
    p[i] = p_IG(rho[i], m[i], E[i]);
  }
}

void
EulerEquation1D::transientResidual(double * res)
{
  unsigned idx = 0;
  unsigned time_step = _problemSystem->getCurrentTimeStep();
  if ((_time_scheme == BDF2) && (time_step > 1))
  {
    // It is typical to use BDF1 for step 1 to startup BDF2
    // however, it is also associated with large error
    // see H. Nishikawa, "On large start-up error of BDF2", Journal of Computational Physics, Vol. 392, 2019, Pages 456-461
    for(unsigned i = 0; i < n_Cell; i++)
    {
      res[idx++] = UTILS::BDF2Tran(rho[i], rho_old[i], rho_oo[i], _dt, _dt_old);
      res[idx++] = UTILS::BDF2Tran(m[i],   m_old[i],   m_oo[i],   _dt, _dt_old);
      res[idx++] = UTILS::BDF2Tran(E[i],   E_old[i],   E_oo[i],   _dt, _dt_old);
    }
  }
  else
  {
    for(unsigned i = 0; i < n_Cell; i++)
    {
      res[idx++] = (rho[i] - rho_old[i]) / _dt;
      res[idx++] = (m[i]   - m_old[i])   / _dt;
      res[idx++] = (E[i]   - E_old[i])   / _dt;
    }
  }
}

void
EulerEquation1D::RHS(double * rhs)
{
  switch (_order)
  {
    case 1:   updateFluxes();           break;
    case 2:   updateFluxes2ndOrder();   break;
    defaut:   sysError("Spatial order not implemented.");
  }

  unsigned idx = 0;
  for(unsigned i = 0; i < n_Cell; i++)
  {
    rhs[idx++] = -(flux_rho[i+1] - flux_rho[i]) / dx;
    rhs[idx++] = -(flux_m[i+1] - flux_m[i]) / dx;
    rhs[idx++] = -(flux_E[i+1] - flux_E[i]) / dx;
  }
}

void
EulerEquation1D::updateFluxes()
{
  for(unsigned i = 0; i <= n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double rho_left = (i == 0) ? RHO_L : rho[i-1];
    double m_left   = (i == 0) ? M_L   : m[i-1];
    double E_left   = (i == 0) ? E_L   : E[i-1];
    double u_left   = m_left / rho_left;
    double p_left   = p_IG(rho_left, m_left, E_left);
    double c_left   = std::sqrt(_gamma * p_left / rho_left);

    double rho_right= (i == n_Cell) ? RHO_R : rho[i];
    double m_right  = (i == n_Cell) ? M_R   : m[i];
    double E_right  = (i == n_Cell) ? E_R   : E[i];
    double u_right  = m_right / rho_right;
    double p_right  = p_IG(rho_right, m_right, E_right);
    double c_right  = std::sqrt(_gamma * p_right / rho_right);

    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m_left + m_right) + 0.5 * eigen_max * (rho_left - rho_right);
    flux_m[i] = 0.5 * (m_left*m_left/rho_left + p_left + m_right*m_right/rho_right + p_right) + 0.5 * eigen_max * (m_left - m_right);
    flux_E[i] = 0.5 * ((E_left + p_left)*u_left + (E_right + p_right)*u_right) + 0.5 * eigen_max * (E_left - E_right);
  }
}

void
EulerEquation1D::updateFluxes2ndOrder()
{
  // Reference:
  //   A. Kurganov and E. Tadmor, New High-Resolution Central Schemes for Nonlinear Conservation Laws
  //   and Convection–Diffusion Equations, Journal of Computational Physics, 160, 241–282 (2000)
  // Note:
  //   The reconstruction is directly performed on the conservative variables.
  //   Some papers (I cannot find the references now FIXME) claim that it is preferred to use primitive
  //   variables for reconstruction. This has not been tested in this code.
  UTILS::linearReconstruction(RHO_L, RHO_R, rho, rho_w, rho_e);
  UTILS::linearReconstruction(M_L, M_R, m, m_w, m_e); // should it be linearReconstruction(-M_L, -M_R, m, m_w, m_e) ?
  UTILS::linearReconstruction(E_L, E_R, E, E_w, E_e);

  for(unsigned i = 0; i <= n_Cell; i++)
  {
    // Find the local max abs(eigenvalue), i.e., Spectral radius
    double rho_left = (i == 0) ? RHO_L : rho_e[i-1];
    double m_left   = (i == 0) ? M_L   : m_e[i-1];
    double E_left   = (i == 0) ? E_L   : E_e[i-1];
    double u_left   = m_left / rho_left;
    double p_left   = p_IG(rho_left, m_left, E_left);
    double c_left   = std::sqrt(_gamma * p_left / rho_left);

    double rho_right= (i == n_Cell) ? RHO_R : rho_w[i];
    double m_right  = (i == n_Cell) ? M_R   : m_w[i];
    double E_right  = (i == n_Cell) ? E_R   : E_w[i];
    double u_right  = m_right / rho_right;
    double p_right  = p_IG(rho_right, m_right, E_right);
    double c_right  = std::sqrt(_gamma * p_right / rho_right);

    double eigen_max = std::max({std::fabs(u_left),  std::fabs(u_left + c_left),   std::fabs(u_left - c_left),
                                 std::fabs(u_right), std::fabs(u_right + c_right), std::fabs(u_right - c_right)});

    flux_rho[i] = 0.5 * (m_left + m_right) + 0.5 * eigen_max * (rho_left - rho_right);
    flux_m[i] = 0.5 * (m_left*m_left/rho_left + p_left + m_right*m_right/rho_right + p_right) + 0.5 * eigen_max * (m_left - m_right);
    flux_E[i] = 0.5 * ((E_left + p_left)*u_left + (E_right + p_right)*u_right) + 0.5 * eigen_max * (E_left - E_right);
  }
}

void
EulerEquation1D::onTimestepEnd()
{
  // save old solutions
  rho_oo  = rho_old;  m_oo  = m_old;  E_oo  = E_old;
  rho_old = rho;      m_old = m;      E_old = E;
}

void
EulerEquation1D::writeVTKOutput(FILE * file)
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
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "rho");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", rho[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "m");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", m[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "E");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", E[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "p");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", p[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", "u");
  for (unsigned i = 0; i < n_Cell; i++)
    fprintf(file, "          %20.6f\n", m[i]/rho[i]);
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </CellData>\n");

  fprintf(file, "    </Piece>\n");
}

void
EulerEquation1D::FillJacobianMatrixNonZeroPattern(MatrixNonZeroPattern * mnzp)
{
  int n_Var = 3;
  for (int i = 0; i < n_Cell; i++)                 // loop on cells
    for (int i_var = 0; i_var < n_Var; i_var++)    // loop on the 3 variables on each cells
    {
      int row = i * n_Var + i_var;                     // This loops on rows
      for (int j = i-2; j <= i+2; j++)             // loop on its -2 to 2 neighbors; maybe out of bound
        for (int j_var = 0; j_var < n_Var; j_var++)
        {
          // This loops on variables on neighboring cells
          int col = j * n_Var + j_var;  // might be smaller than zero
          if ((col >= 0) && (col < _n_DOFs))
            mnzp->addEntry(row + _DOF_offset, col + _DOF_offset);
        }
    }
}

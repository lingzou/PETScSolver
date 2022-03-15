#pragma once

#include <set>
#include "InputParser.h"
#include "ProblemSystem.h"

PetscErrorCode SNESFormFunction(SNES, Vec, Vec, void*);
PetscErrorCode SNESMonitor(SNES, PetscInt, PetscReal, void*);
PetscErrorCode KSPMonitor(KSP, PetscInt, PetscReal, void*);
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);

class ProblemSystem;

class MatrixNonZeroPattern
{
public:
  MatrixNonZeroPattern(unsigned int size) { _non_zero_entries.resize(size); }
  virtual ~MatrixNonZeroPattern() {}
  virtual void addEntry(unsigned row, unsigned col) { _non_zero_entries[row].insert(col); }
  virtual void addRow(unsigned row, std::vector<unsigned> cols) { for (auto& col : cols) addEntry(row, col); }
  std::vector<std::set<unsigned>> &getNonZeroPattern() { return _non_zero_entries; }
protected:
  std::vector<std::set<unsigned>> _non_zero_entries;
};

struct ApplicationCtx
{
  InputParser *   _input_parser;
  ProblemSystem * myProblemSystem;
  PetscInt        N_DOFs;     // Number of degrees of freedom

  SNES            snes;       // SNES
  KSP             ksp;        // KSP

  PetscBool       hasFDColoring;
  MatFDColoring   fdcoloring;   // MatFDColoring
  Mat             J_Mat;        // Jacobian matrix
  Mat             P_Mat;        // Preconditioning matrix
  Mat             J_MatrixFree; // Jacobian-free

  Vec             u;          // unknown vector
  Vec             u_backup;   // unknown vector backup (in case retry a failed time step)

  Vec             r;               // total residual = res_transient - res_RHS
  Vec             res_transient;   /* residual from transient term */
  Vec             res_RHS;         /* residual from RHS */
  Vec             res_RHS_old;     /* residual from RHS (old time step)*/

  void initializePETScApp(InputParser*);
  void setupPETScWorkSpace();
  void setupInitialConditions();
  void FreePETScWorkSpace();
  void setupMatrices();
};

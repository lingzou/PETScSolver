# p_case = 1 for Sod problem; = 2 for Lax problem
[Global]
  ts        = BDF2
  dt        = 1.0e-3
  n_steps   = 40
  output_interval   = 5
  petsc_options = '-snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 1
    length    = 1
    order     = 2
    n_cells   = 96
  [../]
[]

#-snes_converged_reason
#-ksp_converged_reason
#-snes_view
#-pc_factor_shift_type NONZERO
#-pc_factor_shift_amount 1.0e-8
#-snes_mffd_type ds
#-snes_fd_type ds

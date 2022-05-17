[Global]
  ts        = BDF1
  dt        = 1.0e-2
  n_steps   = 50
  output_interval   = 10
  # use tighter convergence for regression test
  petsc_options = '-snes_rtol 1e-10 -snes_atol 1e-8 -snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = SixEqnOneP
    problem_type = 5
    length    = 12
    order     = 1
    n_cells   = 100
  [../]
[]

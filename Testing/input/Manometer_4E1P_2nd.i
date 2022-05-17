[Global]
  ts        = BDF2
  dt        = 5.0e-2
  n_steps   = 600
  output_interval   = 100
  # use tighter convergence for regression test
  petsc_options = '-snes_rtol 1e-10 -snes_atol 1e-8 -snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 3
    length    = 20
    order     = 2
    n_cells   = 40
  [../]
[]

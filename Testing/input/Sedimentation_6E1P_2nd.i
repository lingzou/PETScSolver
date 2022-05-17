[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 1000
  output_interval   = 20
  # use tighter convergence for regression test
  petsc_options = '-snes_linesearch_type basic -snes_rtol 1e-10 -snes_atol 1e-8 -snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = SixEqnOneP
    problem_type = 4
    length    = 2
    order     = 2
    n_cells   = 50
  [../]
[]

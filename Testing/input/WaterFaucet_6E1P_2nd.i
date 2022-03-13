[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 50
  output_interval   = 10
  linear_rtol = 1.e-4
  # default line search has difficulty for dt = 0.01
  petsc_options = '-snes_linesearch_type basic'
[]
[System]
  [./problem]
    type      = SixEqnOneP
    problem_type = 5
    length    = 12
    order     = 2
    n_cells   = 100
  [../]
[]

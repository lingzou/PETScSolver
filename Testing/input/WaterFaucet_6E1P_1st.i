[Global]
  ts        = BDF1
  dt        = 1.0e-2
  n_steps   = 50
  output_interval   = 10
  # default line search has difficulty for dt = 0.01
  petsc_options = '-snes_linesearch_type basic'
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

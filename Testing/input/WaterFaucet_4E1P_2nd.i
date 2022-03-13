[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 50
  output_interval   = 10
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 5
    length    = 12
    order     = 2
    n_cells   = 100
  [../]
[]

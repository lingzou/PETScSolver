[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 1000
  output_interval   = 20
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

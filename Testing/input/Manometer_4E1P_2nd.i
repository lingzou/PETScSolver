[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 60
  output_interval   = 10
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 3
    length    = 20
    order     = 2
    n_cells   = 200
  [../]
[]

[Global]
  ts        = BDF2
  dt        = 5e-3
  n_steps   = 2000
  output_interval   = 100
  linear_rtol = 1e-4
  linear_max_its = 100
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 4
    length    = 2
    order     = 2
    n_cells   = 50
  [../]
[]

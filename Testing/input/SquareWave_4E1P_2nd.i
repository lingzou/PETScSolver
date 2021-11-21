[Global]
  ts        = BDF2
  dt        = 1.0e-3
  n_steps   = 40
  output_interval   = 10
  #petsc_options = '-snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 2
    length    = 1
    order     = 2
    n_cells   = 100
  [../]
[]

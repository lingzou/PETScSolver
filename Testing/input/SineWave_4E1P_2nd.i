[Global]
  ts        = BDF2
  dt        = 1.0e-3
  n_steps   = 40
  output_interval   = 5
  petsc_options = '-snes_stol 1e-10'
[]
[System]
  [./problem]
    type      = FourEqnOneP
    problem_type = 1
    length    = 1
    order     = 2
    n_cells   = 96
  [../]
[]

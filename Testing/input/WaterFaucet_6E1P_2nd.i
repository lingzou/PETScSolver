[Global]
  ts        = BDF2
  dt        = 1.0e-2
  n_steps   = 50
  output_interval   = 10
  linear_rtol = 1.e-4
  # default line search has difficulty for dt = 0.01
  # Using '-snes_linesearch_type basic' would help it converge
  # but here I'd like to test the time step size cut feature
  petsc_options = '-snes_max_it 15'
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

[Global]
  ts        = BDF2
  dt        = 1.0e-2
  #dt_min = 1
  n_steps   = 600
  output_interval   = 100

  linear_rtol = 1.e-4
  linear_max_its = 60

  petsc_options = '-mat_mffd_type ds -mat_fd_type ds'
[]
[System]
  [./problem]
    type      = SixEqnOneP
    problem_type = 3
    length    = 20
    order     = 2
    n_cells   = 40
  [../]
[]

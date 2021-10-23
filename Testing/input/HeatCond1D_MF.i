[Global]
  ts       = CN
  dt       = 0.02
  n_steps  = 40
  solver_option = 1

  fake_int = 2 # This generates a warning as an unused input
[]
[System]
  [./problem]
    type    = HeatConduction1D
    n_cells = 100
  [../]
[]

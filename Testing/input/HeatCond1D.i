[Global]
  ts       = CN
  dt       = 0.02
  n_steps  = 40
[]
[System]
  [./problem]
    type    = HeatConduction1D
    n_cells = 100
  [../]
[]

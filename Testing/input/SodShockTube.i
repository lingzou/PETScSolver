# p_case = 1 for Sod problem; = 2 for Lax problem
[Global]
  ts        = CN
  dt        = 1.0e-3
  n_steps   = 40
[]
[System]
  [./problem]
    type      = EulerEquation1D
    p_case    = 1
    order     = 2
    n_cells   = 100
  [../]
[]

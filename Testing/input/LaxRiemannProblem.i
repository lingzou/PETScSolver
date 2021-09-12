# p_case = 1 for Sod problem; = 2 for Lax problem
[Global]
  ts        = CN
  dt        = 1.0e-3
  n_steps   = 160
[]
[System]
  [./problem]
    type      = EulerEquation1D
    p_case    = 2
    order     = 2
    n_cells   = 200
  [../]
[]

[Global]
  ts       = BDF1
  dt       = 1
  n_steps  = 20
  text_output = true
  output_interval = 5
[]
[System]
  [./problem]
    type     = SinglePhaseFlow
    order    = 1
    n_cells  = 10
    length   = 1
    P_INIT   = 1e5
    V_INIT   = 0.0
    T_INIT   = 300
    V_INLET  = 1
    T_INLET  = 300
    P_OUTLET = 1e5
    T_OUTLET = 300
  [../]
[]

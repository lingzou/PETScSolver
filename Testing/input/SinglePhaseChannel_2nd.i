[Global]
  ts       = BDF2
  dt       = 1
  n_steps  = 20
  text_output = true
  output_interval = 5
[]
[System]
  [./vBC]
    type     = vBC
    Connection = channel-1:begin
    order    = 2
    V_INLET  = 1
    T_INLET  = 300
  [../]
  [./channel-1]
    type     = SinglePhaseChannel
    order    = 2
    n_cells  = 10
    length   = 1
    P_INIT   = 1e5
    V_INIT   = 0.0
    T_INIT   = 300
  [../]
  [./pBC]
    type     = pBC
    Connection = channel-1:end
    order    = 2
    P_OUTLET = 1e5
    T_OUTLET = 300
    V_INIT   = 0
  [../]
[]

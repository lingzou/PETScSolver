[Global]
  ts       = BDF1
  dt       = 1
  n_steps  = 20
  text_output = true
  output_interval = 5
[]

[Fluids]
  [./simpleWater]
    type    = linearFluid
    rho0    = 1000
    p0      = 1e5
    T0      = 300
    e0      = 112.55e3
    drho_dp = 4e-7
    drho_dT = -0.46
    cv      = 4e3
  [../]
[]

[System]
  [./vBC]
    type     = vBC
    Connection = channel-1:begin
    order    = 1
    V_BC     = 1
    T_BC     = 300
  [../]
  [./channel-1]
    type     = SinglePhaseChannel
    order    = 1
    n_cells  = 10
    length   = 1
    P_INIT   = 1e5
    V_INIT   = 0.0
    T_INIT   = 300
    f        = 0.01
    dh       = 0.01
    h        = 2000
    aw       = 300
    Tw       = 350
  [../]
  [./pBC]
    type     = pBC
    Connection = channel-1:end
    order    = 1
    P_BC     = 1e5
    T_BC     = 300
    V_INIT   = 0
  [../]
[]

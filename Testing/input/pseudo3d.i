[Global]
  ts       = BDF2
  dt       = 1
  n_steps  = 20
  text_output = true
  output_interval = 1
  linear_rtol = 1e-4
  linear_max_its = 100
  petsc_options = '-pc_factor_mat_ordering_type rcm'
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
  [./channels]
    type     = Pseudo3D
    order    = 2
    n_cells  = 10
    length   = 2
    P_INIT   = 1e5
    V_INIT   = 1.0
    T_INIT   = 300
  [../]
[]

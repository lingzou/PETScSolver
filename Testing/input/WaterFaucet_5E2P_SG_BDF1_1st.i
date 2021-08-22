# p_case = 1 for Sod problem; = 2 for Lax problem
problem   = FiveEqnTwoP_StagGrid
H_inv     = 1000
order     = 1
ts        = BDF1
dt        = 1.0e-3
n_steps   = 40
n_cells   = 96

output_interval   = 20

#-snes_converged_reason
#-ksp_converged_reason
#-snes_view
#-pc_factor_shift_type NONZERO
#-pc_factor_shift_amount 1.0e-8
#-snes_mffd_type ds
#-snes_fd_type ds

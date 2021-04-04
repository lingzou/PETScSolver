# p_case = 1 for Sod problem; = 2 for Lax problem
-problem  FiveEqnTwoP_StagGrid
-H_inv    2000        # 1/H
-order    1           # First-order spatial scheme
-ts       BDF1        # Backward Euler (aka BDF1)
-dt       5.0e-4      # Time step size
-n_steps  100         # Total number of time steps
-n_cells  768         # Number of finite volumes (cells)

-output_interval 100  # Output interval, write output file every 100 time steps
-text_output true

#-snes_converged_reason
#-ksp_converged_reason
#-snes_view
#-pc_factor_shift_type NONZERO
#-pc_factor_shift_amount 1.0e-8
#-snes_mffd_type ds
#-snes_fd_type ds

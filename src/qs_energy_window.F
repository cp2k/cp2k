!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Does all kind of post scf calculations for GPW/GAPW
!> \par History
!>      Started as a copy from the relevant part of qs_scf
!>      Start to adapt for k-points [07.2015, JGH]
!> \author Joost VandeVondele (10.2003)
! **************************************************************************************************
MODULE qs_energy_window
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_copy_into_existing, dbcsr_create, dbcsr_desymmetrize, &
        dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_type
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE dbcsr_iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_integrate_function
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_window'

   PUBLIC :: energy_windows

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE energy_windows(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'energy_windows', routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: debug = .FALSE.
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      CHARACTER(len=40)                                  :: ext, title
      INTEGER                                            :: handle, i, lanzcos_max_iter, last, nao, &
                                                            nelectron_total, newton_schulz_order, &
                                                            next, nwindows, print_unit, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: stride(:)
      LOGICAL                                            :: mpi_io, print_cube, restrict_range
      REAL(KIND=dp) :: bin_width, density_ewindow_total, density_total, energy_range, fermi_level, &
         filter_eps, frob_norm, lanzcos_threshold, lower_bound, occupation, total_rho, upper_bound
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, P_eigenvalues, &
                                                            window_eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct, window_fm_struct
      TYPE(cp_fm_type), POINTER :: eigenvectors, eigenvectors_nonorth, matrix_ks_fm, &
         P_eigenvectors, P_window_fm, rho_ao_ortho_fm, S_minus_half_fm, tmp_fm, &
         window_eigenvectors, window_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, rho_ao
      TYPE(dbcsr_type)                                   :: matrix_ks_nosym, S_half, S_minus_half, &
                                                            tmp
      TYPE(dbcsr_type), POINTER                          :: rho_ao_ortho, window
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), POINTER                           :: rho_g, rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input, ls_scf_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, matrix_ks=matrix_ks, pw_env=pw_env, rho=rho, &
                      input=input, nelectron_total=nelectron_total, subsys=subsys, ks_env=ks_env, matrix_s=matrix_s)
      CALL qs_subsys_get(subsys, particles=particles)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao)
      IF (SIZE(rho_ao) > 1) CALL cp_warn(__LOCATION__, &
                                         "The printing of energy windows is currently only implemented for clsoe shell systems")
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      !reading the input
      dft_section => section_vals_get_subs_vals(input, "DFT")
      ls_scf_section => section_vals_get_subs_vals(input, "DFT%LS_SCF")
      CALL section_vals_val_get(dft_section, "PRINT%ENERGY_WINDOWS%N_WINDOWS", i_val=nwindows)
      CALL section_vals_val_get(dft_section, "PRINT%ENERGY_WINDOWS%PRINT_CUBES", l_val=print_cube)
      CALL section_vals_val_get(dft_section, "PRINT%ENERGY_WINDOWS%RESTRICT_RANGE", l_val=restrict_range)
      CALL section_vals_val_get(dft_section, "PRINT%ENERGY_WINDOWS%RANGE", r_val=energy_range)
      NULLIFY (stride)
      ALLOCATE (stride(3))
      stride = section_get_ivals(dft_section, "PRINT%ENERGY_WINDOWS%STRIDE")
      CALL section_vals_val_get(dft_section, "PRINT%ENERGY_WINDOWS%EPS_FILTER", r_val=filter_eps)
      CALL section_vals_val_get(ls_scf_section, "EPS_LANCZOS", r_val=lanzcos_threshold)
      CALL section_vals_val_get(ls_scf_section, "MAX_ITER_LANCZOS", i_val=lanzcos_max_iter)
      CALL section_vals_val_get(ls_scf_section, "SIGN_SQRT_ORDER", i_val=newton_schulz_order)

      !Initialize data
      ALLOCATE (window, rho_ao_ortho)
      CALL dbcsr_get_info(matrix=matrix_ks(1)%matrix, nfullrows_total=nao)
      ALLOCATE (eigenvalues(nao))
      CALL dbcsr_create(tmp, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_create(S_minus_half, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_create(S_half, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_create(window, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_create(rho_ao_ortho, template=matrix_s(1)%matrix, matrix_type="N")
      CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, context=blacs_env, nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(P_window_fm, ao_ao_fmstruct)
      CALL cp_fm_create(matrix_ks_fm, ao_ao_fmstruct)
      CALL cp_fm_create(rho_ao_ortho_fm, ao_ao_fmstruct)
      CALL cp_fm_create(S_minus_half_fm, ao_ao_fmstruct)
      CALL cp_fm_create(eigenvectors, ao_ao_fmstruct)
      CALL cp_fm_create(eigenvectors_nonorth, ao_ao_fmstruct)
      NULLIFY (rho_r, rho_g)
      ALLOCATE (rho_r, rho_g)
      NULLIFY (rho_r%pw, rho_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      !calculate S_minus_half
      CALL matrix_sqrt_Newton_Schulz(S_half, S_minus_half, matrix_s(1)%matrix, filter_eps, &
                                     newton_schulz_order, lanzcos_threshold, lanzcos_max_iter)

      !get the full ks matrix
      CALL dbcsr_desymmetrize(matrix_ks(1)%matrix, matrix_ks_nosym)

      !switching to orthonormal basis
      CALL dbcsr_multiply("N", "N", one, S_minus_half, matrix_ks_nosym, zero, tmp, filter_eps=filter_eps)
      CALL dbcsr_multiply("N", "N", one, tmp, S_minus_half, zero, matrix_ks_nosym, filter_eps=filter_eps)
      CALL copy_dbcsr_to_fm(matrix_ks_nosym, matrix_ks_fm)
      CALL dbcsr_multiply("N", "N", one, S_half, rho_ao(1)%matrix, zero, tmp, filter_eps=filter_eps)
      CALL dbcsr_multiply("N", "N", one, tmp, S_half, zero, rho_ao_ortho, filter_eps=filter_eps)
      CALL copy_dbcsr_to_fm(rho_ao_ortho, rho_ao_ortho_fm)

      !diagonalize the full ks matrix
      CALL choose_eigv_solver(matrix_ks_fm, eigenvectors, eigenvalues)
      fermi_level = eigenvalues((nelectron_total + MOD(nelectron_total, 2))/2)
      IF (restrict_range) THEN
         lower_bound = MAX(fermi_level - energy_range, eigenvalues(1))
         upper_bound = MIN(fermi_level + energy_range, eigenvalues(SIZE(eigenvalues)))
      ELSE
         lower_bound = eigenvalues(1)
         upper_bound = eigenvalues(SIZE(eigenvalues))
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) " Creating energy windows. Fermi level: ", fermi_level
         WRITE (unit_nr, *) " Printing Energy Levels from ", lower_bound, " to ", upper_bound
      END IF
      !Rotate the eigenvectors back out of the orthonormal basis
      CALL copy_dbcsr_to_fm(S_minus_half, S_minus_half_fm)
      !calculate the density caused by the mos in the energy window
      CALL cp_gemm("N", "N", nao, nao, nao, one, S_minus_half_fm, eigenvectors, zero, eigenvectors_nonorth)

      IF (debug) THEN
         !check difference to actual density
         CALL cp_fm_struct_create(fmstruct=window_fm_struct, context=blacs_env, nrow_global=nao, &
                                  ncol_global=nelectron_total/2)
         CALL cp_fm_create(window_fm, window_fm_struct)
         CALL cp_fm_to_fm(eigenvectors_nonorth, window_fm, nelectron_total/2, 1, 1)
         CALL cp_gemm("N", "T", nao, nao, nelectron_total/2, 2*one, window_fm, window_fm, zero, P_window_fm)
         !ensure the correct sparsity
         CALL copy_fm_to_dbcsr(P_window_fm, tmp)
         CALL dbcsr_copy(window, matrix_ks(1)%matrix)
         CALL dbcsr_copy_into_existing(window, tmp)
         CALL calculate_rho_elec(matrix_p=window, &
                                 rho=rho_r, &
                                 rho_gspace=rho_g, &
                                 total_rho=total_rho, &
                                 ks_env=ks_env)
         density_total = pw_integrate_function(rho_r%pw)
         IF (unit_nr > 0) WRITE (unit_nr, *) " Ground-state density: ", density_total
         frob_norm = dbcsr_frobenius_norm(window)
         IF (unit_nr > 0) WRITE (unit_nr, *) " Frob norm of calculated ground-state density matrix: ", frob_norm
         CALL dbcsr_add(window, rho_ao(1)%matrix, one, -one)
         frob_norm = dbcsr_frobenius_norm(rho_ao(1)%matrix)
         IF (unit_nr > 0) WRITE (unit_nr, *) " Frob norm of current density matrix: ", frob_norm
         frob_norm = dbcsr_frobenius_norm(window)
         IF (unit_nr > 0) WRITE (unit_nr, *) " Difference between calculated ground-state density and current density: ", frob_norm
         CALL cp_fm_struct_release(window_fm_struct)
         CALL cp_fm_create(tmp_fm, ao_ao_fmstruct)
         CALL cp_fm_to_fm(rho_ao_ortho_fm, tmp_fm)
         CALL cp_fm_create(P_eigenvectors, ao_ao_fmstruct)
         ALLOCATE (P_eigenvalues(nao))
         CALL choose_eigv_solver(tmp_fm, P_eigenvectors, P_eigenvalues)
         CALL cp_fm_create(window_eigenvectors, ao_ao_fmstruct)
         ALLOCATE (window_eigenvalues(nao))
         CALL cp_fm_to_fm(eigenvectors, window_fm, nelectron_total/2, 1, 1)
         CALL cp_gemm("N", "T", nao, nao, nelectron_total/2, 2*one, window_fm, window_fm, zero, P_window_fm)
         CALL choose_eigv_solver(P_window_fm, window_eigenvectors, window_eigenvalues)
         DO i = 1, nao
            IF (unit_nr > 0) THEN
              WRITE (unit_nr, *) i, "H:", eigenvalues(i), "P:", P_eigenvalues(nao - i + 1), "Pnew:", window_eigenvalues(nao - i + 1)
            ENDIF
         END DO
         DEALLOCATE (P_eigenvalues)
         CALL cp_fm_release(tmp_fm)
         CALL cp_fm_release(P_eigenvectors)
         DEALLOCATE (window_eigenvalues)
         CALL cp_fm_release(window_eigenvectors)
         CALL cp_fm_release(window_fm)
      END IF

      !create energy windows
      bin_width = (upper_bound - lower_bound)/nwindows
      next = 0

      DO i = 1, nwindows
         DO WHILE (eigenvalues(next + 1) < lower_bound)
            next = next + 1
         END DO
         last = next
         DO WHILE (eigenvalues(next + 1) < lower_bound + i*bin_width)
            next = next + 1
            IF (next == SIZE(eigenvalues)) EXIT
         END DO
         !calculate the occupation
         !not sure how bad this is now load balanced due to using the same blacs_env
         CALL cp_fm_struct_create(fmstruct=window_fm_struct, context=blacs_env, nrow_global=nao, ncol_global=next - last)
         CALL cp_fm_create(window_fm, window_fm_struct)
         !copy the mos in the energy window into a separate matrix
         CALL cp_fm_to_fm(eigenvectors, window_fm, next - last, last + 1, 1)
         CALL cp_gemm("N", "T", nao, nao, next - last, one, window_fm, window_fm, zero, P_window_fm)
         CALL cp_fm_trace(P_window_fm, rho_ao_ortho_fm, occupation)
         IF (print_cube) THEN
            CALL cp_fm_to_fm(eigenvectors_nonorth, window_fm, next - last, last + 1, 1)
            !print the energy window to a cube file
            !calculate the density caused by the mos in the energy window
            CALL cp_gemm("N", "T", nao, nao, next - last, one, window_fm, window_fm, zero, P_window_fm)
            CALL copy_fm_to_dbcsr(P_window_fm, tmp)
            !ensure the correct sparsity
            CALL dbcsr_copy(window, matrix_ks(1)%matrix)
            CALL dbcsr_copy_into_existing(window, tmp)
            CALL calculate_rho_elec(matrix_p=window, &
                                    rho=rho_r, &
                                    rho_gspace=rho_g, &
                                    total_rho=total_rho, &
                                    ks_env=ks_env)
            WRITE (ext, "(A14,I5.5,A)") "-ENERGY-WINDOW", i, TRIM(cp_iter_string(logger%iter_info))//".cube"
            mpi_io = .TRUE.
            print_unit = cp_print_key_unit_nr(logger, dft_section, "PRINT%ENERGY_WINDOWS", &
                                              extension=ext, file_status="REPLACE", file_action="WRITE", &
                                              log_filename=.FALSE., mpi_io=mpi_io)
            WRITE (title, "(A14,I5)") "ENERGY WINDOW ", i
            CALL cp_pw_to_cube(rho_r%pw, print_unit, title, particles=particles, stride=stride, mpi_io=mpi_io)
            CALL cp_print_key_finished_output(print_unit, logger, dft_section, &
                                              "PRINT%ENERGY_WINDOWS", mpi_io=mpi_io)
            density_ewindow_total = pw_integrate_function(rho_r%pw)
            IF (unit_nr > 0) WRITE (unit_nr, "(A,F16.10,A,I5,A,F20.14,A,F20.14)") " Energy Level: ", &
               lower_bound + (i - 0.5_dp)*bin_width, " Number of states: ", next - last, " Occupation: ", &
               occupation, " Grid Density ", density_ewindow_total
         ELSE
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, "(A,F16.10,A,I5,A,F20.14)") " Energy Level: ", lower_bound + (i - 0.5_dp)*bin_width, &
                  " Number of states: ", next - last, " Occupation: ", occupation
            END IF
         END IF
         CALL cp_fm_release(window_fm)
         CALL cp_fm_struct_release(window_fm_struct)
      END DO

      !clean up
      CALL dbcsr_release(matrix_ks_nosym)
      CALL dbcsr_release(tmp)
      CALL dbcsr_release(window)
      CALL dbcsr_release(S_minus_half)
      CALL dbcsr_release(S_half)
      CALL dbcsr_release(rho_ao_ortho)
      DEALLOCATE (window, rho_ao_ortho)
      CALL cp_fm_struct_release(ao_ao_fmstruct)
      CALL cp_fm_release(matrix_ks_fm)
      CALL cp_fm_release(rho_ao_ortho_fm)
      CALL cp_fm_release(eigenvectors)
      CALL cp_fm_release(P_window_fm)
      CALL cp_fm_release(eigenvectors_nonorth)
      CALL cp_fm_release(S_minus_half_fm)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
      DEALLOCATE (rho_r, rho_g)
      DEALLOCATE (eigenvalues)
      DEALLOCATE (STRIDE)

      CALL timestop(handle)

   END SUBROUTINE energy_windows

!**************************************************************************************************

END MODULE qs_energy_window

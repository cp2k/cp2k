!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the propagation via RT-BSE method.
!> \note  The control is handed directly from cp2k_runs
!> \author Maximilian Graml (03.26)
! **************************************************************************************************

MODULE rt_bse_linearized
   USE bibliography, ONLY: Marek2025, &
                           cite_reference
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE force_env_types, ONLY: force_env_type
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_to_fm, &
                          cp_fm_to_fm_submat_general, &
                          cp_fm_create, &
                          cp_fm_set_all, &
                          cp_fm_release, &
                           cp_fm_get_info
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_create, &
                           cp_cfm_get_info, &
                           cp_cfm_release, &
                           cp_cfm_set_all
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm
   USE kinds, ONLY: dp
   USE cp_dbcsr_api, ONLY: dbcsr_p_type, &
                           dbcsr_type, &
                           dbcsr_create, &
                           dbcsr_release, &
                           dbcsr_copy, &
                           dbcsr_add, &
                           dbcsr_set, &
                           dbcsr_clear, &
                           dbcsr_iterator_type, &
                           dbcsr_iterator_start, &
                           dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block, &
                           dbcsr_reserve_blocks, &
                           dbcsr_get_num_blocks, &
                           dbcsr_get_block_p, &
                           dbcsr_get_info
   USE dbt_api, ONLY: dbt_clear, &
                      dbt_contract, &
                      dbt_copy_matrix_to_tensor, &
                      dbt_copy_tensor_to_matrix, &
                      dbt_type
   USE libint_2c_3c, ONLY: libint_potential_type
   USE qs_tensors, ONLY: build_2c_integrals, &
                         build_2c_neighbor_lists
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, &
                                     release_neighbor_list_sets
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr, &
                                  cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg, ONLY: cp_fm_scale, &
                                 cp_fm_invert, &
                                 cp_fm_transpose, &
                                 cp_fm_column_scale, &
                                 cp_fm_scale_and_add
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_scale_and_add, &
                                  cp_cfm_scale, &
                                  cp_cfm_transpose, &
                                  cp_cfm_norm, &
                                  cp_cfm_trace, &
                                  cp_cfm_column_scale
   USE cp_cfm_diag, ONLY: cp_cfm_geeig
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE qs_moments, ONLY: build_local_moment_matrix
   USE moments_utils, ONLY: get_reference_point
   USE force_env_methods, ONLY: force_env_calc_energy_force
   USE efield_utils, ONLY: make_field
   USE message_passing, ONLY: mp_para_env_type
   USE input_constants, ONLY: rtp_bse_ham_g0w0, &
                              use_mom_ref_zero, &
                              do_bch, &
                              do_exact, &
                              use_rt_restart
    USE rt_bse, ONLY:             initialize_rtbse_env, &
            get_idempotence_deviation, &
            antiherm_metric, &
            init_hartree, &
            get_hartree, &
            get_sigma, &
            rho_metric, &
            cp_cfm_gexp, &
            propagate_density
   USE rt_bse_types, ONLY: rtbse_env_type, &
                           create_rtbse_env, &
                           release_rtbse_env, &
                           multiply_fm_cfm, &
                           multiply_cfm_fm
   USE rt_bse_io, ONLY: output_moments, &
                        read_field, &
                        output_field, &
                        output_mos_contravariant, &
                        read_restart, &
                        output_restart, &
                        print_timestep_info, &
                        print_etrs_info_header, &
                        print_etrs_info, &
                        print_rtbse_header_info
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_add_iter_level, &
                                 cp_rm_iter_level, &
                                 cp_iterate
   USE rt_propagation_output, ONLY: print_ft
   USE rt_propagation_utils, ONLY: read_moments
   USE mathconstants, ONLY: twopi
   USE physcon,                         ONLY: evolt, &
   seconds

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse"

   PUBLIC :: run_propagation_linearized_bse

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time propagation of the linearized BSE
!> \param force_env Force environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_linearized_bse(force_env)
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_linearized_bse'
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER                                            :: i, j, k, handle
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric, enum_re, enum_im, &
                                                            idempotence_dev, a_metric_1, a_metric_2, &
                                                            full_max, oo_max, vv_max, ov_max, vo_max, &
                                                            full_max_old, oo_max_old, vv_max_old, ov_max_old, vo_max_old, &
                                                            ov_growth, vo_growth
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! To Do: Bibliography information

      logger => cp_get_default_logger()

      ! Run the initial SCF calculation / read SCF restart information
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)

      ! Allocate all persistant storage and read input that does not need further processing
      CALL create_rtbse_env(rtbse_env, force_env, linearized=.TRUE.)

      CALL initialize_maximum_timestep(rtbse_env)

      CALL print_linrtbse_header_info(rtbse_env)

      ! Build the truncated MO coefficient slabs C_active(:, first_active_mo..last_active_mo)
      ! used by all AO<->MO transforms in the linearized path.
      CALL populate_C_active(rtbse_env)

      ! Initiate iteration level "MD" in order to copy the structure of other RTP codes
      CALL cp_add_iter_level(logger%iter_info, "MD")
      ! Initialize non-trivial values
      !  - calculates the moment operators
      CALL initialize_moments(rtbse_env)
      !  - populates the initial density matrix
      !     - reads the restart density if requested
      !  - reads the field and moment trace from previous runs
      !  - populates overlap and inverse overlap matrices
      !  - prints some info about loaded files into the output
      CALL initialize_rtbse_env(rtbse_env)

      !  - populates the initial density matrix
      !     - reads the restart density if requested
      CALL initialize_density_matrix(rtbse_env)
      !  - calculates/populates the G0W0/KS Hamiltonian, respectively
      CALL initialize_singleparticle_hamiltonian(rtbse_env)
      ! Transform initial density matrix to AO basis for use in Hartree and self-energy calculations
      DO i = 1, rtbse_env%n_spin
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_orig(i), rtbse_env%rho_ao_lrrtbse(i), i)
      END DO
      !  - calculates the Hartree reference potential
      CALL initialize_hartree_potential(rtbse_env)
      !  - calculates the COHSEX reference self-energy
      CALL initialize_cohsex_selfenergy(rtbse_env)

      ! Debug-only diagnostic: estimate the spectral radius of the linearized
      ! Liouvillian L by power iteration on a random Hermitian seed restricted
      ! to the OV+VO subspace. Prints rho(L) and the corresponding RK4-stable
      ! dt = 2*sqrt(2)/rho(L) next to the gap-based estimate. Compile-time flag.
      IF (.TRUE.) THEN
         CALL estimate_omega_max_krylov(rtbse_env, n_iter=15, rel_tol=1.0E-2_dp)
      END IF

      ! Setup the time based on the starting step
      ! Assumes identical dt between two runs
      rtbse_env%sim_time = REAL(rtbse_env%sim_start, dp)*rtbse_env%sim_dt
      ! Output 0 time moments and field
      IF (.NOT. rtbse_env%restart_extracted) THEN
         CALL output_field(rtbse_env)
         CALL output_moments(rtbse_env, rtbse_env%rho)
      END IF

      ! Do not apply the delta kick if we are doing a restart calculation
        IF (rtbse_env%dft_control%rtp_control%apply_delta_pulse .AND. (.NOT. rtbse_env%restart_extracted)) THEN
           CALL apply_delta_pulse_MO(rtbse_env)
        END IF


      ! ********************** Start the time loop **********************
      ! NOTE : Time-loop starts at index sim_start = 0, unless restarted or configured otherwise
      DO i = rtbse_env%sim_start, rtbse_env%sim_nsteps - 1

         ! Update the simulation time
         rtbse_env%sim_time = REAL(i, dp)*rtbse_env%sim_dt
         rtbse_env%sim_step = i

         CALL solve_rk4_timestep(rtbse_env, rtbse_env%rho, rtbse_env%rho_new, converged, k, metric)
         DO j = 1, rtbse_env%n_spin
            !MG for debugging
            CALL get_rho_block_maxnorms(rtbse_env, rtbse_env%rho(j), j, &
                                        full_max_old, oo_max_old, vv_max_old, ov_max_old, vo_max_old)
            CALL get_rho_block_maxnorms(rtbse_env, rtbse_env%rho_new(j), j, &
                                        full_max, oo_max, vv_max, ov_max, vo_max)
            ov_growth = ov_max/MAX(ov_max_old, TINY(1.0_dp))
            vo_growth = vo_max/MAX(vo_max_old, TINY(1.0_dp))
            ! IF (rtbse_env%unit_nr > 0) THEN
            !    WRITE (rtbse_env%unit_nr, '(A,I3, A, 5ES12.3)') 'MG: Block max norms for spin ', j, ': ', &
            !       full_max, oo_max, vv_max, ov_max, vo_max 
            !    WRITE (rtbse_env%unit_nr, '(A,I3,A,2ES12.3,A,2ES12.3)') 'MG: OV/VO growth for spin ', j, &
            !       ' (new/old): ', ov_growth, vo_growth, ' | old OV/VO: ', ov_max_old, vo_max_old
            ! END IF
         END DO
         CALL get_electron_number_MO(rtbse_env, rtbse_env%rho_new, enum_re, enum_im)
         IF (.FALSE.) THEN
            ! Not all of these are used, but they are all good metrics to check the convergence in problematic cases
            ! TODO : Allow for conditional warning
            CALL get_idempotence_deviation(rtbse_env, rtbse_env%rho_new, idempotence_dev)
            DO j = 1, rtbse_env%n_spin
               CALL cp_cfm_to_fm(rtbse_env%sigma_SEX(j), rtbse_env%real_workspace(1), rtbse_env%real_workspace(2))
               CALL antiherm_metric(real_fm=rtbse_env%real_workspace(1), imag_fm=rtbse_env%real_workspace(2), &
                                    workspace=rtbse_env%rho_workspace, metric=a_metric_1)
               CALL antiherm_metric(real_fm=rtbse_env%hartree_curr(j), &
                                    workspace=rtbse_env%rho_workspace, metric=a_metric_2)
            END DO
         END IF
         CALL print_timestep_info(rtbse_env, i, enum_re)
         !IF (.NOT. converged) CPABORT("ETRS did not converge")
         CALL cp_iterate(logger%iter_info, iter_nr=i, last=(i == rtbse_env%sim_nsteps-1))


         ! Update rho
         DO j = 1, rtbse_env%n_spin
            CALL cp_cfm_to_cfm(rtbse_env%rho_new(j), rtbse_env%rho(j))
         END DO
         ! Print the updated field
         CALL output_field(rtbse_env)
         ! If needed, print out the density matrix in MO basis
         CALL output_mos_contravariant(rtbse_env, rtbse_env%rho, rtbse_env%rho_section)
         ! Also handles outputting to memory
         CALL output_moments(rtbse_env, rtbse_env%rho)
         ! Output restart files, so that the restart starts at the following time index
         CALL output_restart(rtbse_env, rtbse_env%rho, i + 1)
      END DO
      ! ********************** End the time loop **********************

      CALL cp_rm_iter_level(logger%iter_info, "MD")

      ! Carry out the FT
      CALL print_ft(rtbse_env%rtp_section, &
                    rtbse_env%moments_trace, &
                    rtbse_env%time_trace, &
                    rtbse_env%field_trace, &
                    rtbse_env%dft_control%rtp_control, &
                    info_opt=rtbse_env%unit_nr)

      ! Deallocate everything
      CALL release_rtbse_env(rtbse_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_linearized_bse

! **************************************************************************************************
!> \brief Writes the header and basic info (incl maximum timestep inside stability region) 
!>        to the standard output
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
      SUBROUTINE initialize_maximum_timestep(rtbse_env)
         TYPE(rtbse_env_type), POINTER                      :: rtbse_env
         CHARACTER(len=*), PARAMETER                        :: routineN = "initialize_maximum_timestep"
         CHARACTER(len=256)                                 :: hint_msg
         REAL(kind=dp)                                      :: omega_max, sim_dt_as, total_time
         INTEGER                                            :: handle, i_first, i_last, n_steps_old, n_steps_new

         CALL timeset(routineN, handle)

         i_first = rtbse_env%first_active_mo
         i_last = rtbse_env%last_active_mo

         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
            omega_max = MAXVAL(rtbse_env%bs_env%eigenval_G0W0(i_first:i_last, :, :)) - &
                        MINVAL(rtbse_env%bs_env%eigenval_G0W0(i_first:i_last, :, :))
         ELSE
            omega_max = MAXVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:i_last, :)) - &
                        MINVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:i_last, :))
         END IF

         rtbse_env%omega_max = omega_max

         IF (omega_max > 0.0_dp) THEN
            rtbse_env%maximum_timestep = 2.0_dp*SQRT(2.0_dp)/omega_max
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "Error in estimating maximum timestep: largest KS/GW gap is "// &
                          "non-positive. Check the active MO window (cutoffs) and the "// &
                          "eigenvalues.")
         END IF

         IF (rtbse_env%sim_dt <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "TIMESTEP must be positive for linearized RT-BSE. Use RTBSE%ENFORCE_MAX_DT "// &
                          "with a positive TIMESTEP to automatically rewrite TIMESTEP and STEPS.")
         END IF

         n_steps_old = MAX(0, rtbse_env%sim_nsteps)

         IF (rtbse_env%dft_control%rtp_control%enforce_max_dt) THEN
            total_time = REAL(n_steps_old, dp)*rtbse_env%sim_dt
            ! Grace factor of 4
            n_steps_new = MAX(1, CEILING(total_time/(0.25_dp*rtbse_env%maximum_timestep)))
            rtbse_env%sim_dt = total_time/REAL(n_steps_new, dp)
            rtbse_env%sim_nsteps = n_steps_new
            sim_dt_as = rtbse_env%sim_dt*seconds*1e18_dp
            WRITE (hint_msg, '(A,F16.4,A,I0,A)') &
               'ENFORCE_MAX_DT enabled. Resetting TIMESTEP to ', sim_dt_as, &
               ' as and STEPS to ', n_steps_new, '.'
            CALL cp_hint(__LOCATION__, TRIM(hint_msg))
         END IF

         CALL timestop(handle)
      END SUBROUTINE initialize_maximum_timestep

   ! **************************************************************************************************
   !> \brief Writes the header and basic info (incl maximum timestep inside stability region) 
   !>        to the standard output
   !> \param rtbse_env Entry point - rtbse environment
   ! **************************************************************************************************
   SUBROUTINE print_linrtbse_header_info(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp)                                      :: total_time, fft_resolution, nyquist_frequency
      INTEGER                                            :: n_steps

      logger => cp_get_default_logger()
      n_steps = MAX(0, rtbse_env%sim_nsteps)
      total_time = REAL(n_steps, dp)*rtbse_env%sim_dt
      fft_resolution = 0.0_dp
      nyquist_frequency = 0.0_dp
      IF (total_time > 0.0_dp) fft_resolution = twopi/total_time
      IF (rtbse_env%sim_dt > 0.0_dp) nyquist_frequency = twopi/(2.0_dp*rtbse_env%sim_dt)

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, *) ''
         WRITE (rtbse_env%unit_nr, '(A)') ' /-----------------------------------------------'// &
            '------------------------------\'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |              Linearized Real Time Bethe-Salpeter Propagation'// &
            '                |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' \-----------------------------------------------'// &
            '------------------------------/'
         WRITE (rtbse_env%unit_nr, *) ''


         WRITE (rtbse_env%unit_nr, '(A18,L62)') ' Apply delta pulse', &
            rtbse_env%dft_control%rtp_control%apply_delta_pulse

         WRITE (rtbse_env%unit_nr, '(A)') ''

         WRITE (rtbse_env%unit_nr, '(A,T65,F16.4)') ' Estimated maximum timestep within stability region [as]:', &
            rtbse_env%maximum_timestep*seconds*1e18_dp
         WRITE (rtbse_env%unit_nr, '(A,T65,F16.4)') ' Applied timestep [as]:', &
            rtbse_env%sim_dt*seconds*1e18_dp
         WRITE (rtbse_env%unit_nr, '(A,T71,I10)') ' Number of propagation steps:', n_steps
         WRITE (rtbse_env%unit_nr, '(A,T65,F16.4)') ' Total propagation time [as]:', &
            total_time*seconds*1e18_dp
         WRITE (rtbse_env%unit_nr, '(A,T65,F16.6)') ' Estimated FFT frequency resolution without interpolation [eV]:', &
            fft_resolution*evolt
         WRITE (rtbse_env%unit_nr, '(A,T65,F16.6)') ' Nyquist frequency [eV]:', &
            nyquist_frequency*evolt
         WRITE (rtbse_env%unit_nr, '(A,T65,F16.6)') ' Estimated maximum oscillation frequency (gap-based) [eV]:', &
            rtbse_env%omega_max*evolt

         ! Active MO window (energy-cutoff truncation) for linearized RT-BSE
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
            WRITE (rtbse_env%unit_nr, '(A,T75,A6)') ' Active-window single-particle spectrum source:', '  G0W0'
         ELSE
            WRITE (rtbse_env%unit_nr, '(A,T75,A6)') ' Active-window single-particle spectrum source:', '    KS'
         END IF
         IF (rtbse_env%rtbse_energy_cutoff_occ > 0.0_dp) THEN
            WRITE (rtbse_env%unit_nr, '(A,T71,F10.3)') ' Active-window occupied energy cutoff [eV]:', &
               rtbse_env%rtbse_energy_cutoff_occ*evolt
         ELSE
            WRITE (rtbse_env%unit_nr, '(A,T71,A10)') ' Active-window occupied energy cutoff [eV]:', '  disabled'
         END IF
         IF (rtbse_env%rtbse_energy_cutoff_empty > 0.0_dp) THEN
            WRITE (rtbse_env%unit_nr, '(A,T71,F10.3)') ' Active-window virtual energy cutoff [eV]:', &
               rtbse_env%rtbse_energy_cutoff_empty*evolt
         ELSE
            WRITE (rtbse_env%unit_nr, '(A,T71,A10)') ' Active-window virtual energy cutoff [eV]:', '  disabled'
         END IF
         WRITE (rtbse_env%unit_nr, '(A,T71,I10)') ' First active occupied MO index:', rtbse_env%first_active_mo
         WRITE (rtbse_env%unit_nr, '(A,T71,I10)') ' Last active virtual MO index:', rtbse_env%last_active_mo
         WRITE (rtbse_env%unit_nr, '(A,T71,I10)') ' Number of active MOs:', rtbse_env%mo_active
      END IF

   END SUBROUTINE print_linrtbse_header_info

! **************************************************************************************************
!> \brief Populates rtbse_env%C_active(i_spin) (n_ao x mo_active) by extracting columns
!>        first_active_mo..last_active_mo from bs_env%fm_mo_coeff_Gamma(i_spin).
!> \param rtbse_env RT-BSE environment
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE populate_C_active(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "populate_C_active"
      INTEGER                                      :: i, handle

      CALL timeset(routineN, handle)

      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_set_all(rtbse_env%C_active(i), 0.0_dp)
         CALL cp_fm_to_fm_submat_general( &
            rtbse_env%bs_env%fm_mo_coeff_Gamma(i), rtbse_env%C_active(i), &
            rtbse_env%n_ao, rtbse_env%mo_active, &
            1, rtbse_env%first_active_mo, &
            1, 1, &
            rtbse_env%bs_env%fm_mo_coeff_Gamma(i)%matrix_struct%context)
      END DO

      CALL timestop(handle)
   END SUBROUTINE populate_C_active

! **************************************************************************************************
!> \brief Calculates the single particle Hamiltonian
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - refactor in prep. of linearized propagation
! **************************************************************************************************
   SUBROUTINE initialize_moments(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_moments"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: moments_dbcsr_p
      INTEGER                                      :: i,k, handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: matrix_s
      REAL(kind=dp), DIMENSION(3)                  :: rpoint
      TYPE(cp_fm_type)                             :: tmp_ao

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env, matrix_s=matrix_s)

      ! AO-sized scratch buffer for moment matrices before transform to MO basis
      CALL cp_fm_create(tmp_ao, bs_env%fm_s_Gamma%matrix_struct)

      ! ****** START MOMENTS OPERATOR CALCULATION
      ! Construct moments from dbcsr
      NULLIFY (moments_dbcsr_p)
      ALLOCATE (moments_dbcsr_p(3))
      DO k = 1, 3
         ! Make sure the pointer is empty
         NULLIFY (moments_dbcsr_p(k)%matrix)
         ! Allocate a new matrix that the pointer points to
         ALLOCATE (moments_dbcsr_p(k)%matrix)
         ! Create the matrix storage - matrix copies the structure of overlap matrix
         CALL dbcsr_copy(moments_dbcsr_p(k)%matrix, matrix_s(1)%matrix)
      END DO
      ! Run the moment calculation
      ! check for presence to prevent memory errors
      rpoint(:) = 0.0_dp
      CALL get_reference_point(rpoint, qs_env=rtbse_env%qs_env, &
                               reference=rtbse_env%moment_ref_type, ref_point=rtbse_env%user_moment_ref_point)
      CALL build_local_moment_matrix(rtbse_env%qs_env, moments_dbcsr_p, 1, rpoint)
      ! Copy to AO scratch then transform to MO-active
      DO k = 1, 3
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, tmp_ao)
         CALL transform_ao_to_mo_covariant_fm(rtbse_env, tmp_ao, rtbse_env%moments(k), 1)
      END DO
      ! MG to do: remove moments_field; only needed for TDDFT comparison
      ! Now, repeat without reference point to get the moments for field
      CALL get_reference_point(rpoint, qs_env=rtbse_env%qs_env, &
                               reference=use_mom_ref_zero)
      CALL build_local_moment_matrix(rtbse_env%qs_env, moments_dbcsr_p, 1, rpoint)
      DO k = 1, 3
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, tmp_ao)
         CALL transform_ao_to_mo_covariant_fm(rtbse_env, tmp_ao, rtbse_env%moments_field(k), 1)
      END DO

      ! Now can deallocate dbcsr matrices
      DO k = 1, 3
         CALL dbcsr_release(moments_dbcsr_p(k)%matrix)
         DEALLOCATE (moments_dbcsr_p(k)%matrix)
      END DO
      DEALLOCATE (moments_dbcsr_p)
      CALL cp_fm_release(tmp_ao)
      ! ****** END MOMENTS OPERATOR CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_moments

   
! **************************************************************************************************
!> \brief Calculates the initial density matrix, based on the SCF density
!> \param rtbse_env RT-BSE environment
!> \author Maximilian Graml (03/26) 
! **************************************************************************************************
   SUBROUTINE initialize_density_matrix(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_density_matrix"
      INTEGER                                      :: i, handle, ii, jj, nrow_global, ncol_global, &
                                                      nrow_local, ncol_local, &
                                                      i_row_global, j_col_global
      INTEGER, DIMENSION(:), POINTER               :: row_indices, col_indices


      CALL timeset(routineN, handle)

      ! Get distribution of MO-active workspace
      CALL cp_fm_get_info(rtbse_env%real_workspace_mo(1), &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      ! Iterate over both spins
      DO i = 1, rtbse_env%n_spin
         !Ensure that workspace is set to 0
         CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
               IF (i_row_global == j_col_global .AND. &
                   (i_row_global + rtbse_env%first_active_mo - 1) <= rtbse_env%n_occ(i)) THEN
                  rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = 1.0_dp
               END IF
            END DO
         END DO
         ! Sets imaginary part to zero
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), mtarget=rtbse_env%rho(i))
         ! Save the reference value for the case of delta kick
         CALL cp_cfm_to_cfm(rtbse_env%rho(i), rtbse_env%rho_orig(i))
      END DO
      
      !MG to do: restart reading of density matrix

      CALL timestop(handle)
   END SUBROUTINE initialize_density_matrix

! **************************************************************************************************
!> \brief Calculates the single particle Hamiltonian in the MO basis
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - refactor in prep. of linearized propagation
! **************************************************************************************************
   SUBROUTINE initialize_singleparticle_hamiltonian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_singleparticle_hamiltonian"
      !MG to remove - H2 energies from N^4 GW
      ! REAL(KIND=dp), DIMENSION(12), PARAMETER      :: last_column = [ &
      !    -15.7854_dp,   9.0657_dp,  21.3465_dp,  33.0050_dp, &
      !     37.8467_dp,  37.8467_dp,  57.2756_dp,  56.9618_dp, &
      !     56.9619_dp,  83.8232_dp,  99.2828_dp, 116.7532_dp  &
      ! ]
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      INTEGER                                      :: i, ii,jj, handle, nrow_global, ncol_global, &
                                                      nrow_local, ncol_local, &
                                                      i_row_global, j_col_global, abs_mo_idx
      INTEGER, DIMENSION(:), POINTER               :: row_indices, col_indices

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)

      ! Get distribution of MO-active workspace
      CALL cp_fm_get_info(rtbse_env%real_workspace_mo(1), &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      !Ensure that workspace is set to 0
      CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
      ! ****** START SINGLE PARTICLE HAMILTONIAN CALCULATION
      DO i = 1, rtbse_env%n_spin
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
                  IF (i_row_global == j_col_global) THEN
                     abs_mo_idx = i_row_global + rtbse_env%first_active_mo - 1
                     IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
                        ! G0W0 Hamiltonian
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = bs_env%eigenval_G0W0(abs_mo_idx, 1, i)
                     ELSE
                        ! KS Hamiltonian
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = bs_env%eigenval_scf_Gamma(abs_mo_idx, i)
                     END IF
                     ! For testing purposes, use the hacked eigenvalues to check the correctness of the implementation
                     !MG to remove
                     !rtbse_env%real_workspace(1)%local_data(ii, jj) = last_column(i_row_global)/evolt
                  END IF
            END DO
         END DO
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), mtarget=rtbse_env%ham_reference_singleparticle(i))
      END DO
      ! ****** END SINGLE PARTICLE HAMILTONIAN CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_singleparticle_hamiltonian

! **************************************************************************************************
!> \brief Calculates the Hartree potential
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_hartree_potential(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_hartree_potential"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      INTEGER                                      :: i, handle

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! Calculate Coulomb RI elements, necessary for Hartree calculation
      CALL init_hartree(rtbse_env, rtbse_env%v_dbcsr)
      ! Calculate the original Hartree potential
      ! Uses rho_orig - same as rho for initial run but different for continued run
      DO i = 1, rtbse_env%n_spin
         CALL get_hartree(rtbse_env, rtbse_env%rho_ao_lrrtbse(i), rtbse_env%hartree_curr_ao(i))
         ! Scaling by spin degeneracy
         CALL cp_fm_scale(rtbse_env%spin_degeneracy, rtbse_env%hartree_curr_ao(i))
         ! Transform to MO basis (AO scratch -> MO-active result)
         CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%hartree_curr_ao(i), rtbse_env%hartree_curr(i), i)
         ! Apply occupation factor f_n-f_m
         CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%hartree_curr(i), i)
         ! Subtract the reference from the reference Hamiltonian
         ! following H_eff = ... + V_Hartree(delta_rho) - V_Hartree(rho_0), 
         CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(i), mtarget=rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(0.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                   CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
      END DO
      ! ****** END HARTREE POTENTIAL REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_hartree_potential

! **************************************************************************************************
!> \brief Calculates the single particle Hamiltonian
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_cohsex_selfenergy(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_cohsex_selfenergy"
      INTEGER                                      :: i, handle

      CALL timeset(routineN, handle)

      ! ****** START COHSEX REFERENCE CALCULATION
      ! w_dbcsr and screened_dbt are needed for get_sigma routines
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
         ! In a non-HF calculation, copy the actual correlation part of the interaction
         CALL copy_fm_to_dbcsr(rtbse_env%bs_env%fm_W_MIC_freq_zero, rtbse_env%w_dbcsr)
      ELSE
         ! In HF, correlation is set to zero
         CALL dbcsr_set(rtbse_env%w_dbcsr, 0.0_dp)
      END IF
      ! Add the Hartree to the screened_dbt tensor - now W = V + W^c
      CALL dbcsr_add(rtbse_env%w_dbcsr, rtbse_env%v_dbcsr, 1.0_dp, 1.0_dp)
      CALL dbt_copy_matrix_to_tensor(rtbse_env%w_dbcsr, rtbse_env%screened_dbt)
      ! Calculate the COHSEX starting energies
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
         ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
         DO i = 1, rtbse_env%n_spin
            ! TODO : Allow no COH calculation for static screening
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH_ao(i), -0.5_dp, rtbse_env%S_inv_fm)
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%sigma_COH_ao(i), rtbse_env%sigma_COH(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%sigma_COH(i), i)
            ! Copy and subtract from the complex reference hamiltonian
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(i), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
            ! Calculate exchange part - TODO : should this be applied for different spins? - TEST with O2 HF propagation?
            ! So far only closed shell tested
            ! Uses rho_orig - same as rho for initial run but different for continued run
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(i), -1.0_dp, rtbse_env%rho_ao_lrrtbse(i))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(i), rtbse_env%sigma_SEX(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(i), i)
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      ELSE
         ! KS Hamiltonian - use time-dependent Fock exchange
         DO i = 1, rtbse_env%n_spin
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(i), -1.0_dp, rtbse_env%rho_ao_lrrtbse(i))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(i), rtbse_env%sigma_SEX(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(i), i)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      END IF
      ! ****** END COHSEX REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_cohsex_selfenergy


! **************************************************************************************************
!> \brief Propagates the density to the next timestep using Runge-Kutta 4th order (RK4) method
!> \param rho_start Initial density matrix
!> \param rho_end Final density matrix 
!> \param converged Whether the resulting rho_end is self-consistent
!> \param k How many SC iterations were done
!> \param metric The difference metric from the last self-consistent iteration (for printing/evaluation)
! **************************************************************************************************
   SUBROUTINE solve_rk4_timestep(rtbse_env, rho_start, rho_end, converged, k, metric)
      TYPE(rtbse_env_type), POINTER                     :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: rho_start, &
                                                           rho_mid, &
                                                           rho_end
      LOGICAL                                           :: converged
      INTEGER                                           :: k
      REAL(kind=dp)                                     :: metric
      CHARACTER(len=*), PARAMETER                       :: routineN = "solve_rk4_timestep"
      INTEGER                                           :: i, handle
      REAL(kind=dp)                                     :: norm_k1, norm_rho, omega_runtime
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER        :: a_ptr, b_ptr, c_ptr
      TYPE(mp_para_env_type), POINTER                   :: para_env

      CALL timeset(routineN, handle)
      ! RK4 follows the typical scheme
      ! d/dt ρ = f(t, ρ)
      ! f(t,ρ) = -i( Δε Δρ + (f_n-f_m) V_Hartree(Δρ) + ΔΣ(Δρ) ) 
      ! i.e. RK4 reads
      ! k_1 = f(t, ρ_start)
      ! k_2 = f(t + dt/2, ρ_start + dt/2 * k_1)
      ! k_3 = f(t + dt/2, ρ_start + dt/2 * k_2)
      ! k_4 = f(t + dt, ρ_start + dt * k_3)
      ! ρ_end = ρ_start + dt/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
      ! Note that the effective Hamiltonian needs to be updated for each evaluation of f, 
      ! as it depends on the density matrix at the respective time

      DO i = 1, rtbse_env%n_spin
         ! k_1
         CALL update_effective_ham_MO(rtbse_env, rho_start(i), rtbse_env%rk4_coefficients(1), i)
         ! Runtime estimate of the highest oscillating frequency contained in the
         ! linearised RT-BSE RHS f(t,rho) = -i [H_eff, rho]. The static
         ! occupied-diagonal block of rho does not oscillate, so we measure the
         ! deviation from the reference density:
         !   ||f||_F = ||d(drho)/dt||_F  <=  ||L||_2 ||drho||_F,    drho = rho - rho_orig
         ! the ratio ||k_1||_F / ||drho||_F is a sharp instantaneous lower bound
         ! on the spectral radius of the linear propagator L acting on drho and
         ! can be compared against the gap-based omega_max estimate to refine
         ! the stability region. Printed only at the first integrated step to
         ! keep the runtime overhead negligible.
         IF (rtbse_env%sim_step == rtbse_env%sim_start) THEN
            ! NOTE: para_env%sum is collective. ALL ranks must enter this branch,
            ! not just the I/O rank, otherwise the run deadlocks inside MPI.
            para_env => rho_start(i)%matrix_struct%para_env
            a_ptr => rtbse_env%rk4_coefficients(1)%local_data
            b_ptr => rho_start(i)%local_data
            c_ptr => rtbse_env%rho_orig(i)%local_data

            norm_k1 = SUM(REAL(a_ptr*CONJG(a_ptr), kind=dp))
            norm_rho = SUM(REAL((b_ptr - c_ptr)*CONJG(b_ptr - c_ptr), kind=dp))

            CALL para_env%sum(norm_k1)
            CALL para_env%sum(norm_rho)

            norm_k1 = SQRT(norm_k1)
            norm_rho = SQRT(norm_rho)
            IF (norm_rho > 0.0_dp) THEN
               omega_runtime = norm_k1/norm_rho
            ELSE
               omega_runtime = 0.0_dp
            END IF

            IF (rtbse_env%unit_nr > 0) THEN
               WRITE (rtbse_env%unit_nr, '(A,I3,A,F16.6)') &
                  ' Runtime instantaneous max oscillation frequency (spin', i, &
                  ') ||drho/dt||_F/||delta_rho||_F [eV]:', omega_runtime*evolt
            END IF
         END IF
         ! MG TODO: dt/2 needs to implemented in field
         !          Problem: sim_step is needed as integer for envelope
         ! k_2 - needs to be evaluated at t + dt/2 and rho + dt/2 * k_1
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(0.5_dp*rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(1))
         CALL update_effective_ham_MO(rtbse_env,  rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(2), i)
         ! k_3 - needs to be evaluated at t + dt/2 and rho + dt/2 * k_2
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(0.5_dp*rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(2))
         CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(3), i)
         ! k_4 - needs to be evaluated at t + dt and rho + dt * k_3
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(3))
         CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(4), i)
         ! Now calculate the final rho
         CALL cp_cfm_to_cfm(rho_start(i), rho_end(i))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rho_end(i), &
                                   CMPLX(rtbse_env%sim_dt/6.0_dp, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rho_end(i), &
                                   CMPLX(rtbse_env%sim_dt/3.0_dp, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(2))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rho_end(i), &
                                   CMPLX(rtbse_env%sim_dt/3.0_dp, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(3))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rho_end(i), &
                                   CMPLX(rtbse_env%sim_dt/6.0_dp, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(4))
      END DO

      ! Update bookkeeping to the next timestep similar to logic of etrs_scf_loop
      rtbse_env%sim_step = rtbse_env%sim_step + 1
      rtbse_env%sim_time = rtbse_env%sim_time + rtbse_env%sim_dt
      
      CALL timestop(handle)
   END SUBROUTINE solve_rk4_timestep

   ! **************************************************************************************************
!> \brief Updates the effective Hamiltonian in the MO basis, given a density matrix rho
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
!> \param ham_effective Effective Hamiltonian in the MO basis that is updated in this routine
! **************************************************************************************************
   SUBROUTINE update_effective_ham_MO(rtbse_env, rho, ham_effective, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)           :: rho, ham_effective
      INTEGER                                           :: ispin
      CHARACTER(len=*), PARAMETER                        :: routineN = "update_effective_ham_MO"
      LOGICAL, PARAMETER                                 :: debug_disable_sex = .FALSE.
      LOGICAL, PARAMETER                                 :: debug_disable_hartree = .FALSE.
      LOGICAL, PARAMETER                                 :: debug_disable_coh = .FALSE.
      INTEGER                                            :: k, j, nspin, handle

      CALL timeset(routineN, handle)
      ! Shorthand
      !nspin = rtbse_env%n_spin

      ! Update density matrix in AO from latest rho as it is needed for Hartree and SEX contributions
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rho, rtbse_env%rho_ao_lrrtbse(ispin), ispin)

      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0 - reference COHSEX - reference Hartree
      !DO j = 1, nspin
         ! Sets the imaginary part to zero
         CALL cp_cfm_to_cfm(rtbse_env%ham_reference(ispin), ham_effective)
         ! Multiply energy differences by delta rho and add to effective Hamiltonian
         ! MG TODO: Code that in local_data chunks later on
         CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                         CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference_singleparticle(ispin), &
                         rho, CMPLX(1.0, 0.0, kind=dp), ham_effective )
         CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                         CMPLX(-1.0, 0.0, kind=dp), rho, rtbse_env%ham_reference_singleparticle(ispin), &
                         CMPLX(1.0, 0.0, kind=dp), ham_effective )
      !END DO
      ! Determine the field at current time
      IF (rtbse_env%dft_control%apply_efield_field) THEN
         CPABORT("E(t) NYI")
         ! MG needs trafo to MO
         !CALL make_field(rtbse_env%dft_control, rtbse_env%field, rtbse_env%sim_step, rtbse_env%sim_time)
      ELSE
         ! No field
         rtbse_env%field(:) = 0.0_dp
      END IF
      !DO j = 1, nspin
         DO k = 1, 3
         ! MG needs trafo to MO
            ! Minus sign due to charge of electrons
            ! CALL cp_fm_to_cfm(msourcer=rtbse_env%moments_field(k), mtarget=rtbse_env%ham_workspace(1))
            ! CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
            !                           CMPLX(rtbse_env%field(k), 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END DO
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0 .AND. .NOT. debug_disable_coh) THEN
            ! Add the COH part - so far static but can be dynamic in principle through the W updates
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH_ao(ispin), -0.5_dp, rtbse_env%S_inv_fm)
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%sigma_COH_ao(ispin), rtbse_env%sigma_COH(ispin), ispin)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%sigma_COH(ispin), ispin)
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(ispin), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                      CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END IF
         ! Calculate the (S)EX part - based on provided rho
         ! iGW = - rho W
         IF (.NOT. debug_disable_sex) THEN
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, rtbse_env%rho_ao_lrrtbse(ispin))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), rtbse_env%sigma_SEX(ispin), ispin)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(ispin), ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                      CMPLX(1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(ispin))
         END IF
         ! Calculate Hartree potential
         ! Hartree potential is scaled by number of electrons in each MO - spin degeneracy
         IF (.NOT. debug_disable_hartree) THEN
            CALL get_hartree(rtbse_env, rtbse_env%rho_ao_lrrtbse(ispin), &
                             rtbse_env%hartree_curr_ao(ispin))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%hartree_curr_ao(ispin), rtbse_env%hartree_curr(ispin), ispin)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%hartree_curr(ispin), ispin)
            CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(ispin), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                      CMPLX(rtbse_env%spin_degeneracy, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END IF
         ! Keep only the anti-Hermitian part so the RHS remains structure-preserving.
         ! CALL cp_cfm_transpose(ham_effective, 'C', rtbse_env%ham_workspace(1))
         ! CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, kind=dp), ham_effective, &
         !                           CMPLX(-0.5_dp, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1))
         ! Return the actual RHS f(t,rho) = -i * A(rho) for RK4
         CALL cp_cfm_scale(CMPLX(0.0_dp, -1.0_dp, kind=dp), ham_effective)
      !END DO

      CALL timestop(handle)
   END SUBROUTINE update_effective_ham_MO

! **************************************************************************************************
!> \brief Transforms from AO to MO basis
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_covariant_fm(rtbse_env, fm_ao, fm_mo, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                         :: fm_ao, fm_mo
      INTEGER, INTENT(IN)                                :: i_spin
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_ao_to_mo_covariant_fm"
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! M_uv * C_active_vn  -> (n_ao x mo_active)
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, fm_ao, rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      ! (C_active_um)^T * (n_ao x mo_active) -> (mo_active x mo_active)
      CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%ao_mo_workspace(1), &
                         0.0_dp, fm_mo)

      CALL timestop(handle)
   END SUBROUTINE transform_ao_to_mo_covariant_fm

! **************************************************************************************************
!> \brief Transforms from AO to MO basis for complex matrices
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_covariant_cfm(rtbse_env, fm_ao, fm_mo, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                         :: fm_ao, fm_mo
      INTEGER, INTENT(IN)                                :: i_spin
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_ao_to_mo_covariant_cfm"
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Decompose into real/imag AO-sized parts
      CALL cp_cfm_to_fm(msource=fm_ao, mtargetr=rtbse_env%real_workspace(1), &
                        mtargeti=rtbse_env%real_workspace(2))
      ! REAL : C_active^T * (real_workspace(1) * C_active) -> mo_active x mo_active
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%real_workspace(1), rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%ao_mo_workspace(1), &
                         0.0_dp, rtbse_env%real_workspace_mo(1))
      ! IMAG
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%real_workspace(2), rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%ao_mo_workspace(1), &
                         0.0_dp, rtbse_env%real_workspace_mo(2))
      ! Reassemble into MO-sized cfm
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), &
                        msourcei=rtbse_env%real_workspace_mo(2), &
                        mtarget=fm_mo)

      CALL timestop(handle)
   END SUBROUTINE transform_ao_to_mo_covariant_cfm

! **************************************************************************************************
!> \brief Transforms from MO to AO basis for contravariant complex matrices (e.g. density matrix) 
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_mo_to_ao_contravariant_cfm(rtbse_env, fm_mo, fm_ao, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                         :: fm_mo, fm_ao
      INTEGER, INTENT(IN)                                :: i_spin
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_mo_to_ao_contravariant_cfm"
      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: tmp_ao_mo

      CALL timeset(routineN, handle)

      ! Rectangular n_ao x mo_active complex intermediate
      CALL cp_cfm_create(tmp_ao_mo, rtbse_env%fm_struct_ao_mo_active)

      ! C_active * fm_mo  -> (n_ao x mo_active)
      CALL multiply_fm_cfm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%mo_active, &
                           1.0_dp, rtbse_env%C_active(i_spin), fm_mo, &
                           0.0_dp, tmp_ao_mo)
      ! (n_ao x mo_active) * C_active^T -> (n_ao x n_ao)
      CALL multiply_cfm_fm("N", "T", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%mo_active, &
                           1.0_dp, tmp_ao_mo, rtbse_env%C_active(i_spin), &
                           0.0_dp, fm_ao)

      CALL cp_cfm_release(tmp_ao_mo)

      CALL timestop(handle)
   END SUBROUTINE transform_mo_to_ao_contravariant_cfm

! **************************************************************************************************
!> \brief Applies the occupation factor difference f_n - f_m 
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_mo_occupation_factor_diff_cfm(rtbse_env, cfm, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                         :: cfm
      INTEGER                                           :: i_spin
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_mo_occupation_factor_diff_cfm"
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_cfm_to_fm(msource=cfm, mtargetr=rtbse_env%real_workspace_mo(1), &
                        mtargeti=rtbse_env%real_workspace_mo(2))
      ! REAL
      CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%real_workspace_mo(1), i_spin)
      ! IMAGINARY
      CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%real_workspace_mo(2), i_spin)
      ! Copy back to cfm
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), &
                        msourcei=rtbse_env%real_workspace_mo(2), &
                        mtarget=cfm)

      CALL timestop(handle)
   END SUBROUTINE transform_mo_occupation_factor_diff_cfm
   

   SUBROUTINE transform_mo_occupation_factor_diff_fm(rtbse_env, fm, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                         :: fm
      INTEGER                                                :: i_spin
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_mo_occupation_factor_diff_fm"
      INTEGER                                            :: handle, i_local, j_local, i_global, j_global, n_occ, &
                                                            i_global_mo, j_global_mo, shift
      INTEGER                                            :: nrow_local, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices, col_indices
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data
      REAL(kind=dp)                                      :: occ_factor

      CALL timeset(routineN, handle)

      n_occ = rtbse_env%n_occ(i_spin)
      ! Shift mapping local active-window index to absolute MO index
      shift = rtbse_env%first_active_mo - 1

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      local_data => fm%local_data

      DO i_local = 1, nrow_local
         i_global = row_indices(i_local)
         i_global_mo = i_global + shift
         DO j_local = 1, ncol_local
            j_global = col_indices(j_local)
            j_global_mo = j_global + shift

            IF (i_global_mo <= n_occ .AND. j_global_mo > n_occ) THEN
               occ_factor = -1.0_dp
            ELSE IF (i_global_mo > n_occ .AND. j_global_mo <= n_occ) THEN
               occ_factor = 1.0_dp
            ELSE
               occ_factor = 0.0_dp
            END IF
            
            local_data(i_local, j_local) = occ_factor*local_data(i_local, j_local)
         END DO
      END DO


      CALL timestop(handle)
   END SUBROUTINE transform_mo_occupation_factor_diff_fm

! **************************************************************************************************
!> \brief Custom reimplementation of the delta pulse routines
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml - trafo to MO and linearized version following 10.1021/acs.jctc.2c00644 (03.26)
! **************************************************************************************************
   SUBROUTINE apply_delta_pulse_MO(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "apply_delta_pulse_MO"
      REAL(kind=dp)                                :: intensity, metric
      REAL(kind=dp), DIMENSION(3)                  :: kvec
      INTEGER                                      :: i, k, handle

      CALL timeset(routineN, handle)

      ! Report application
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A28)') ' RTBSE| Applying delta pulse'
      ! Extra minus for the propagation of density
      intensity = -rtbse_env%dft_control%rtp_control%delta_pulse_scale
      metric = 0.0_dp
      kvec(:) = rtbse_env%dft_control%rtp_control%delta_pulse_direction(:)
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A38,E14.4E3,E14.4E3,E14.4E3)') &
         " RTBSE| Delta pulse elements (a.u.) : ", intensity*kvec(:)
      ! So far no spin dependence, but can be added by different structure of delta pulse
      CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
      DO k = 1, 3
         CALL cp_fm_scale_and_add(1.0_dp, rtbse_env%real_workspace_mo(1), &
                                  kvec(k), rtbse_env%moments_field(k))
      END DO
      ! enforce hermiticity of the effective Hamiltonian
      CALL cp_fm_transpose(rtbse_env%real_workspace_mo(1), rtbse_env%real_workspace_mo(2))
      CALL cp_fm_scale_and_add(0.5_dp, rtbse_env%real_workspace_mo(1), &
                               0.5_dp, rtbse_env%real_workspace_mo(2))
      ! Prepare the exponential/exponent for propagation
      ! MG: For now only use BCH, where Overlap matrix can easily be removed
      ! multiply by intensity
      ! MG to do: linearize delta pulse
      CALL cp_fm_scale(intensity, rtbse_env%real_workspace_mo(1))
      DO i = 1, rtbse_env%n_spin
         ! Sets real part to zero
         CALL cp_fm_to_cfm(msourcei=rtbse_env%real_workspace_mo(1), mtarget=rtbse_env%ham_workspace(i))
      END DO
      ! Propagate the density by the effect of the delta pulse
      CALL propagate_density(rtbse_env, rtbse_env%ham_workspace, rtbse_env%rho, rtbse_env%rho_new)
      metric = rho_metric(rtbse_env%rho_new, rtbse_env%rho, rtbse_env%n_spin)
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, ('(A42,E38.8E3)')) " RTBSE| Metric difference after delta kick", metric
      ! Copy the new density to the old density
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rtbse_env%rho_new(i), rtbse_env%rho(i))
      END DO

      CALL timestop(handle)
   END SUBROUTINE apply_delta_pulse_MO

! **************************************************************************************************
!> \brief Determines maxnorm of rho and its OO/VV/OV/VO blocks for one spin channel
!> \param rtbse_env Entry point - rtbse environment
!> \param rho Density matrix in MO basis for one spin channel
!> \param i_spin Spin index
!> \param full_max Maxnorm over all matrix elements
!> \param oo_max Maxnorm in occupied-occupied block
!> \param vv_max Maxnorm in virtual-virtual block
!> \param ov_max Maxnorm in occupied-virtual block
!> \param vo_max Maxnorm in virtual-occupied block
! **************************************************************************************************
   SUBROUTINE get_rho_block_maxnorms(rtbse_env, rho, i_spin, full_max, oo_max, vv_max, ov_max, vo_max)
      TYPE(rtbse_env_type), INTENT(IN)                    :: rtbse_env
      TYPE(cp_cfm_type), INTENT(IN)                       :: rho
      INTEGER, INTENT(IN)                                 :: i_spin
      REAL(kind=dp), INTENT(OUT)                          :: full_max, oo_max, vv_max, ov_max, vo_max

      INTEGER                                             :: i_local, j_local, i_global, j_global, n_occ, &
                                                             i_global_mo, j_global_mo, shift
      INTEGER                                             :: nrow_local, ncol_local
      INTEGER, DIMENSION(:), POINTER                      :: row_indices, col_indices
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER          :: local_data
      REAL(kind=dp)                                       :: aij
      TYPE(mp_para_env_type), POINTER                     :: para_env

      full_max = 0.0_dp
      oo_max = 0.0_dp
      vv_max = 0.0_dp
      ov_max = 0.0_dp
      vo_max = 0.0_dp

      n_occ = rtbse_env%n_occ(i_spin)
      shift = rtbse_env%first_active_mo - 1

      CALL cp_cfm_get_info(matrix=rho, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      local_data => rho%local_data

      DO i_local = 1, nrow_local
         i_global = row_indices(i_local)
         i_global_mo = i_global + shift
         DO j_local = 1, ncol_local
            j_global = col_indices(j_local)
            j_global_mo = j_global + shift

            aij = ABS(local_data(i_local, j_local))
            full_max = MAX(full_max, aij)

            IF (i_global_mo <= n_occ .AND. j_global_mo <= n_occ) THEN
               oo_max = MAX(oo_max, aij)
            ELSE IF (i_global_mo > n_occ .AND. j_global_mo > n_occ) THEN
               vv_max = MAX(vv_max, aij)
            ELSE IF (i_global_mo <= n_occ .AND. j_global_mo > n_occ) THEN
               ov_max = MAX(ov_max, aij)
            ELSE
               vo_max = MAX(vo_max, aij)
            END IF
         END DO
      END DO

      CALL get_qs_env(rtbse_env%qs_env, para_env=para_env)
      CALL para_env%max(full_max)
      CALL para_env%max(oo_max)
      CALL para_env%max(vv_max)
      CALL para_env%max(ov_max)
      CALL para_env%max(vo_max)
   END SUBROUTINE get_rho_block_maxnorms

! **************************************************************************************************
!> \brief Outputs the number of electrons in the system from the density matrix
!> \param rtbse_env Entry point - rtbse environment
!> \param rho Density matrix in MO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_electron_number_MO(rtbse_env, rho, electron_n_re, electron_n_im)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      REAL(kind=dp), INTENT(OUT)                         :: electron_n_re, electron_n_im

      CHARACTER(len=*), PARAMETER                  :: routineN = "get_electron_number_MO"


      INTEGER                                            :: handle, j, i_local, j_local, i_global, j_global
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: nrow_local, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices, col_indices

      electron_n_re = 0.0_dp
      electron_n_im = 0.0_dp
      DO j = 1, rtbse_env%n_spin
         CALL cp_cfm_get_info(matrix=rho(j), &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)
         local_data => rho(j)%local_data
         ! Extract and sum diagonal elements
         DO i_local = 1, nrow_local
            i_global = row_indices(i_local)
            ! Search column indices for the diagonal position
            DO j_local = 1, ncol_local
               j_global = col_indices(j_local)
               IF (j_global == i_global) THEN
                  ! Found diagonal element
                  electron_n_re = electron_n_re + REAL(local_data(i_local, j_local), kind=dp)
                  electron_n_im = electron_n_im + AIMAG(local_data(i_local, j_local))
                  EXIT
               END IF
            END DO
         END DO
      END DO
      ! Scale by spin degeneracy
      electron_n_re = electron_n_re*rtbse_env%spin_degeneracy
      electron_n_im = electron_n_im*rtbse_env%spin_degeneracy
   END SUBROUTINE get_electron_number_MO

! **************************************************************************************************
!> \brief Debug routine: estimate the spectral radius rho(L) of the linearized
!>        Liouvillian L by power iteration on a random Hermitian seed delta-rho.
!>        Reuses update_effective_ham_MO as the matvec; subtracts the constant
!>        offset f(rho=0) to expose the purely linear action L*delta_rho.
!>        For the BSE Liouvillian the eigenvalues sit on (or near) the imaginary
!>        axis, so RK4 stability is governed by rho(L) directly: dt < 2*sqrt(2)/rho(L).
!>        Note: the operator 2-norm sigma_max(L) cannot be obtained from
!>        L-only matvecs without a true adjoint (which is not available here).
!>        The empirical safety factor needed in practice (~4) reflects
!>        non-normal transient amplification, not a missing factor between
!>        rho(L) and sigma_max(L); see option B (embedded RK45) for a proper
!>        per-step error driver.
!>        Restores rho_ao_lrrtbse from rho_orig before returning.
!> \param rtbse_env Entry point - rtbse environment
!> \param n_iter Maximum number of power iterations
!> \param rel_tol Relative tolerance on rho(L) for early break
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE estimate_omega_max_krylov(rtbse_env, n_iter, rel_tol)
      TYPE(rtbse_env_type), POINTER                     :: rtbse_env
      INTEGER, INTENT(IN)                               :: n_iter
      REAL(kind=dp), INTENT(IN)                         :: rel_tol
      CHARACTER(len=*), PARAMETER                       :: routineN = "estimate_omega_max_krylov"
      INTEGER                                           :: handle, it, ispin, &
                                                           i_local, j_local, &
                                                           nrow_local, ncol_local, seed_size
      INTEGER, DIMENSION(:), POINTER                    :: row_indices, col_indices
      INTEGER, ALLOCATABLE, DIMENSION(:)                :: rng_seed
      REAL(kind=dp)                                     :: omega, omega_prev, dt_max, &
                                                           omega_gap, dt_max_gap, norm_seed, rr, ri
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER        :: seed_data
      TYPE(cp_cfm_type)                                 :: seed, v, k0, zero_rho, work_h

      CALL timeset(routineN, handle)

      ispin = 1
      omega_gap = rtbse_env%omega_max
      dt_max_gap = 0.0_dp
      IF (omega_gap > 0.0_dp) dt_max_gap = 2.0_dp*SQRT(2.0_dp)/omega_gap

      ! Allocate scratch on MO-active struct
      CALL cp_cfm_create(seed,     rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(v,        rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(k0,       rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(zero_rho, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(work_h,   rtbse_env%fm_struct_mo_active)

      ! 1) Constant offset k0 = update_effective_ham_MO(rho = 0). update_effective_ham_MO
      !    is affine in rho; subtracting k0 yields the pure linear action L*seed.
      CALL cp_cfm_set_all(zero_rho, CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL update_effective_ham_MO(rtbse_env, zero_rho, k0, ispin)

      ! 2) Random Hermitian seed (no subspace projection; let the iteration
      !    converge to the dominant eigendirection of L).
      CALL cp_cfm_get_info(matrix=seed, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      CALL RANDOM_SEED(SIZE=seed_size)
      ALLOCATE (rng_seed(seed_size))
      rng_seed = 20260506
      CALL RANDOM_SEED(PUT=rng_seed)
      DEALLOCATE (rng_seed)
      seed_data => seed%local_data
      DO j_local = 1, ncol_local
         DO i_local = 1, nrow_local
            CALL RANDOM_NUMBER(rr)
            CALL RANDOM_NUMBER(ri)
            seed_data(i_local, j_local) = CMPLX(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
         END DO
      END DO
      ! Initial Hermitization (BSE delta-rho is Hermitian)
      CALL cp_cfm_transpose(seed, 'C', work_h)
      CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, kind=dp), seed, &
                                CMPLX(0.5_dp, 0.0_dp, kind=dp), work_h)
      norm_seed = cp_cfm_norm(seed, 'F')
      IF (norm_seed <= 0.0_dp) THEN
         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A)') &
               ' KRYLOV-OMEGA: seed norm is zero on this distribution; skipping.'
         END IF
         GOTO 999
      END IF
      CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), seed)

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' ============= KRYLOV rho(L) DEBUG ESTIMATE ============='
         WRITE (rtbse_env%unit_nr, '(A,F16.6,A,F16.4,A)') &
            ' Gap-based  omega_max = ', omega_gap*evolt, ' eV   ->   dt_max = ', &
            dt_max_gap*seconds*1e18_dp, ' as'
         WRITE (rtbse_env%unit_nr, '(A,A4,3X,A16,3X,A16,3X,A16)') &
            ' ', 'iter', 'rho(L) [eV]', 'dt_max [as]', 'rel.change'
      END IF

      ! 3) Power iteration: rho(L) = lim_k ||L^k seed||^{1/k}.
      omega = 0.0_dp
      omega_prev = 0.0_dp
      dt_max = 0.0_dp
      DO it = 1, n_iter
         ! v = L * seed
         CALL update_effective_ham_MO(rtbse_env, seed, v, ispin)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                   CMPLX(-1.0_dp, 0.0_dp, kind=dp), k0)

         omega = cp_cfm_norm(v, 'F')
         IF (omega > 0.0_dp) dt_max = 2.0_dp*SQRT(2.0_dp)/omega

         IF (rtbse_env%unit_nr > 0) THEN
            IF (it == 1) THEN
               WRITE (rtbse_env%unit_nr, '(A,I4,3X,F16.6,3X,F16.4,3X,A16)') &
                  ' ', it, omega*evolt, dt_max*seconds*1e18_dp, '       --'
            ELSE
               WRITE (rtbse_env%unit_nr, '(A,I4,3X,F16.6,3X,F16.4,3X,ES16.3)') &
                  ' ', it, omega*evolt, dt_max*seconds*1e18_dp, &
                  ABS(omega - omega_prev)/MAX(omega, TINY(1.0_dp))
            END IF
         END IF

         IF (it > 1 .AND. ABS(omega - omega_prev) <= rel_tol*MAX(omega, TINY(1.0_dp))) EXIT
         omega_prev = omega

         ! Next seed: re-Hermitize (round-off cleanup) and renormalize
         CALL cp_cfm_to_cfm(v, seed)
         CALL cp_cfm_transpose(seed, 'C', work_h)
         CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, kind=dp), seed, &
                                   CMPLX(0.5_dp, 0.0_dp, kind=dp), work_h)
         norm_seed = cp_cfm_norm(seed, 'F')
         IF (norm_seed <= 0.0_dp) EXIT
         CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), seed)
      END DO

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A,F16.6,A,F16.4,A)') &
            ' Krylov     rho(L)    = ', omega*evolt, ' eV   ->   dt_max = ', &
            dt_max*seconds*1e18_dp, ' as'
         IF (omega_gap > 0.0_dp) THEN
            WRITE (rtbse_env%unit_nr, '(A,F16.6)') &
               ' Ratio rho(L) / omega_gap = ', omega/omega_gap
         END IF
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' ========================================================'
      END IF

999   CONTINUE
      ! Restore rho_ao_lrrtbse from rho_orig (update_effective_ham_MO has
      ! overwritten it with the seed AO image during the iteration)
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_orig(ispin), &
                                                rtbse_env%rho_ao_lrrtbse(ispin), ispin)

      CALL cp_cfm_release(seed)
      CALL cp_cfm_release(v)
      CALL cp_cfm_release(k0)
      CALL cp_cfm_release(zero_rho)
      CALL cp_cfm_release(work_h)

      CALL timestop(handle)
   END SUBROUTINE estimate_omega_max_krylov


END MODULE
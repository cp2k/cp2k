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
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_norm,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: &
        cp_cfm_create, cp_cfm_get_info, cp_cfm_get_submatrix, cp_cfm_release, cp_cfm_set_all, &
        cp_cfm_set_element, cp_cfm_set_submatrix, cp_cfm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
        cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE dbt_api,                         ONLY: dbt_copy_matrix_to_tensor
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_type
   USE input_constants,                 ONLY: rtp_bse_ham_g0w0,&
                                              use_mom_ref_zero
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt,&
                                              seconds
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power
   USE rt_bse,                          ONLY: antiherm_metric,&
                                              get_hartree,&
                                              get_idempotence_deviation,&
                                              get_sigma,&
                                              init_hartree,&
                                              initialize_rtbse_env,&
                                              propagate_density,&
                                              rho_metric
   USE rt_bse_io,                       ONLY: output_field,&
                                              output_moments,&
                                              output_mos_contravariant,&
                                              output_restart,&
                                              print_timestep_info
   USE rt_bse_ri_rs,                    ONLY: compute_hartree_ri_rs,&
                                              compute_hartree_ri_rs_complex,&
                                              compute_sigma_ri_rs,&
                                              compute_sigma_ri_rs_complex,&
                                              rt_bse_ri_rs_ensure_V_grid,&
                                              rt_bse_ri_rs_ensure_W0_grid
   USE rt_bse_types,                    ONLY: create_rtbse_env,&
                                              multiply_cfm_fm,&
                                              multiply_fm_cfm,&
                                              release_rtbse_env,&
                                              rtbse_env_type
   USE rt_propagation_output,           ONLY: print_ft
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_linearized'

   ! Set .TRUE. to print the Krylov rho(L) diagnostic next to the analytical
   ! gap-based omega_max estimate. Diagnostic only
   LOGICAL, PARAMETER :: debug_krylov_estimate = .FALSE.

   ! Validation-only: H2 G0W0 QPs from a bse_full_diag.F TDA+ABBA reference run (2026-05-26,
   ! single SCF, single GW for both TDA and ABBA eigenvalues). Used by both the propagator's
   ! reference H and the diagnostic's eps source to align QP inputs across the two pipelines
   ! for a clean cross-check. REMOVE this PARAMETER plus both hack call sites after validation.
   REAL(KIND=dp), DIMENSION(12), PARAMETER, PRIVATE      :: last_column = [ &
                                                            -15.7837_dp, 9.0662_dp, 21.4040_dp, 33.0524_dp, &
                                                            37.7741_dp, 37.7848_dp, 59.7905_dp, 55.5754_dp, &
                                                            55.5908_dp, 82.1245_dp, 112.9021_dp, 114.1125_dp &
                                                            ]

   PUBLIC :: run_propagation_linearized_bse

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time propagation of the linearized BSE
!> \param force_env Force environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_linearized_bse(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'run_propagation_linearized_bse'

      INTEGER                                            :: handle, i, j, k
      LOGICAL                                            :: converged
      REAL(kind=dp) :: a_metric_1, a_metric_2, enum_im, enum_re, full_max, full_max_old, &
         idempotence_dev, metric, oo_max, oo_max_old, ov_growth, ov_max, ov_max_old, t_phys, &
         timestep_walltime, timestep_walltime_start, vo_growth, vo_max, vo_max_old, vv_max, &
         vv_max_old
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_lab
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

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
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_set_all(rtbse_env%ham_reference(i), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      END DO
      CALL initialize_hartree_potential(rtbse_env)
      !  - calculates the COHSEX reference self-energy
      CALL initialize_cohsex_selfenergy(rtbse_env)

      ! Krylov rho(L) diagnostics: print the kernel-inclusive spectral radius
      ! of the linearized propagator alongside the bare-gap omega_max so the
      ! analytical bound can be checked. Power-iteration (Krylov) variant is a
      ! slow lower bound; the Lanczos variant gives the largest-magnitude Ritz
      ! value of the tridiagonal projection and converges superlinearly under
      ! TDA (operator is self-adjoint on Hermitian OV+VO drho). Both are
      ! diagnostic only - they do NOT overwrite omega_max / maximum_timestep.
      ! refine_omega_shift_via_krylov remains DEAD CODE (degenerates for
      ! H2-like single-peak systems).
      ! [INTERNAL] dev_notes/rt_bse_pitfalls_physics.md P4 (Krylov lower-bound issue)
      IF (debug_krylov_estimate .AND. rtbse_env%tda_active .AND. &
          rtbse_env%tda_shift_to_first_peak) THEN
         CALL estimate_omega_max_krylov(rtbse_env, n_iter=30, rel_tol=1.0e-5_dp)
         ! Hermitian Lanczos variant: in TDA the kernel response is self-adjoint on
         ! Hermitian OV+VO drho, so plain Lanczos applies (Arnoldi reduces to Lanczos)
         ! and Ritz values converge superlinearly toward |lambda_max(L)| - typically
         ! an order of magnitude faster than power iteration on clustered top spectra.
         CALL estimate_omega_max_lanczos(rtbse_env, n_iter=30, rel_tol=1.0e-8_dp)
      END IF

      ! Liouvillian eigenvalue diagnostic (one-shot at init, TDA or ABBA via dispatcher).
      ! Detached from the propagator; safe to call after the reference-init routines.
      ! [INTERNAL] dev_notes/rt_bse_observables.md O-D10; dev_notes/rt_bse_equations.md §8.5
      IF (rtbse_env%diagnose_liouvillian_eig) THEN
         CALL diagnose_liouvillian_eigenvalues(rtbse_env)
      END IF

      ! Setup the time based on the starting step
      ! Assumes identical dt between two runs
      rtbse_env%sim_time = REAL(rtbse_env%sim_start, dp)*rtbse_env%sim_dt
      NULLIFY (rho_lab)
      ! Output 0 time moments and field
      IF (.NOT. rtbse_env%restart_extracted) THEN
         CALL output_field(rtbse_env)
         CALL build_rho_lab(rtbse_env, rtbse_env%rho, rtbse_env%sim_time, rho_lab)
         CALL output_moments(rtbse_env, rho_lab)
      END IF

      ! Do not apply the delta kick if we are doing a restart calculation
      IF (rtbse_env%dft_control%rtp_control%apply_delta_pulse .AND. (.NOT. rtbse_env%restart_extracted)) THEN
         CALL apply_delta_pulse_MO(rtbse_env)
      END IF

      ! ********************** Start the time loop **********************
      ! NOTE : Time-loop starts at index sim_start = 0, unless restarted or configured otherwise
      DO i = rtbse_env%sim_start, rtbse_env%sim_nsteps - 1
         timestep_walltime_start = m_walltime()

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
         timestep_walltime = m_walltime() - timestep_walltime_start
         CALL print_timestep_info(rtbse_env, i, enum_re, step_walltime=timestep_walltime)
         !IF (.NOT. converged) CPABORT("ETRS did not converge")
         CALL cp_iterate(logger%iter_info, iter_nr=i, last=(i == rtbse_env%sim_nsteps - 1))

         ! Update rho
         DO j = 1, rtbse_env%n_spin
            CALL cp_cfm_to_cfm(rtbse_env%rho_new(j), rtbse_env%rho(j))
         END DO
         ! Print the updated field
         CALL output_field(rtbse_env)
         ! rho is the rotating-frame density at physical time t_phys = (i+1)*dt.
         ! Build a lab-frame copy once and feed it to all observable/restart sinks.
         t_phys = REAL(i + 1, dp)*rtbse_env%sim_dt
         CALL build_rho_lab(rtbse_env, rtbse_env%rho, t_phys, rho_lab)
         ! If needed, print out the density matrix in MO basis
         CALL output_mos_contravariant(rtbse_env, rho_lab, rtbse_env%rho_section)
         ! Also handles outputting to memory
         CALL output_moments(rtbse_env, rho_lab)
         ! Output restart files, so that the restart starts at the following time index
         CALL output_restart(rtbse_env, rho_lab, i + 1)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_maximum_timestep'

      CHARACTER(len=256)                                 :: hint_msg
      INTEGER                                            :: handle, i_first, i_last, ispin, &
                                                            n_steps_new, n_steps_old
      REAL(kind=dp) :: eps_max_ai, eps_min_ai, eps_occ_max, eps_occ_min, eps_virt_max, &
         eps_virt_min, ev_tmp, grace_factor, omega_max, sim_dt_as, total_time

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

      ! First-peak shift (TDA only): set Omega_0 = eps_min_ai so the lowest
      ! active OV mode oscillates at zero frequency in the rotating frame
      ! (RK4-exact for peak 1). omega_max is the full active OV width
      ! Delta = eps_max_ai - eps_min_ai, where eps_ai = eps_a - eps_i runs
      ! over the active OV pairs only (i in active occupied, a in active
      ! virtual).
      rtbse_env%omega_shift = 0.0_dp
      IF (rtbse_env%tda_active .AND. rtbse_env%tda_shift_to_first_peak) THEN
         eps_occ_min = HUGE(0.0_dp)
         eps_occ_max = -HUGE(0.0_dp)
         eps_virt_min = HUGE(0.0_dp)
         eps_virt_max = -HUGE(0.0_dp)
         DO ispin = 1, rtbse_env%n_spin
            ! Active occupied window: first_active_mo .. n_occ(ispin)
            IF (rtbse_env%n_occ(ispin) >= i_first) THEN
               IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_G0W0(i_first:rtbse_env%n_occ(ispin), :, ispin))
                  eps_occ_min = MIN(eps_occ_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_G0W0(i_first:rtbse_env%n_occ(ispin), :, ispin))
                  eps_occ_max = MAX(eps_occ_max, ev_tmp)
               ELSE
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:rtbse_env%n_occ(ispin), ispin))
                  eps_occ_min = MIN(eps_occ_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:rtbse_env%n_occ(ispin), ispin))
                  eps_occ_max = MAX(eps_occ_max, ev_tmp)
               END IF
            END IF
            ! Active virtual window: n_occ(ispin)+1 .. last_active_mo
            IF (rtbse_env%n_occ(ispin) < i_last) THEN
               IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_G0W0(rtbse_env%n_occ(ispin) + 1:i_last, :, ispin))
                  eps_virt_min = MIN(eps_virt_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_G0W0(rtbse_env%n_occ(ispin) + 1:i_last, :, ispin))
                  eps_virt_max = MAX(eps_virt_max, ev_tmp)
               ELSE
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%n_occ(ispin) + 1:i_last, ispin))
                  eps_virt_min = MIN(eps_virt_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%n_occ(ispin) + 1:i_last, ispin))
                  eps_virt_max = MAX(eps_virt_max, ev_tmp)
               END IF
            END IF
         END DO

         IF (eps_occ_max > -HUGE(0.0_dp) .AND. eps_virt_min < HUGE(0.0_dp)) THEN
            eps_min_ai = eps_virt_min - eps_occ_max
            eps_max_ai = eps_virt_max - eps_occ_min
            rtbse_env%omega_shift = eps_min_ai
            omega_max = eps_max_ai - eps_min_ai
            IF (rtbse_env%unit_nr > 0) THEN
               WRITE (rtbse_env%unit_nr, '(A)') &
                  " RTBSE| ---------- First-peak shift diagnostics (TDA, active OV pairs) ----------"
               WRITE (rtbse_env%unit_nr, '(A,F14.6,A,F14.6)') &
                  " RTBSE|   eps_occ  [eV] min / max =", eps_occ_min*evolt, &
                  " /", eps_occ_max*evolt
               WRITE (rtbse_env%unit_nr, '(A,F14.6,A,F14.6)') &
                  " RTBSE|   eps_virt [eV] min / max =", eps_virt_min*evolt, &
                  " /", eps_virt_max*evolt
               WRITE (rtbse_env%unit_nr, '(A,F14.6,A,F14.6)') &
                  " RTBSE|   eps_ai   [eV] min / max =", eps_min_ai*evolt, &
                  " /", eps_max_ai*evolt
               WRITE (rtbse_env%unit_nr, '(A,F14.6)') &
                  " RTBSE|   omega_shift [eV]        =", rtbse_env%omega_shift*evolt
               WRITE (rtbse_env%unit_nr, '(A,F14.6)') &
                  " RTBSE|   omega_max   [eV] (full) =", omega_max*evolt
               WRITE (rtbse_env%unit_nr, '(A)') &
                  " RTBSE| ------------------------------------------------------------------------"
            END IF
         ELSE
            ! Active window has no genuine OV pair - fall back to no shift
            rtbse_env%omega_shift = 0.0_dp
            rtbse_env%tda_shift_to_first_peak = .FALSE.
         END IF
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

      IF (rtbse_env%enforce_max_dt) THEN
         total_time = REAL(n_steps_old, dp)*rtbse_env%sim_dt
         ! Full code needs a grace factor of 4
         IF (rtbse_env%tda_active) THEN
            grace_factor = 1.0_dp
         ELSE
            grace_factor = 4.0_dp
         END IF
         n_steps_new = MAX(1, CEILING(total_time/(rtbse_env%maximum_timestep/grace_factor)))
         rtbse_env%sim_dt = total_time/REAL(n_steps_new, dp)
         rtbse_env%sim_nsteps = n_steps_new
         sim_dt_as = rtbse_env%sim_dt*seconds*1e18_dp
         WRITE (hint_msg, '(A,F16.4,A,I0,A)') &
            'ENFORCE_MAX_DT enabled. Resetting TIMESTEP to ', sim_dt_as, &
            ' as and STEPS to ', n_steps_new, '.'
         CALL cp_hint(__LOCATION__, TRIM(hint_msg))
      END IF

      IF (rtbse_env%sim_nsteps /= n_steps_old) THEN
         CALL reallocate_ft_traces(rtbse_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE initialize_maximum_timestep

! **************************************************************************************************
!> \brief Reallocate FT trace buffers after ENFORCE_MAX_DT rewrites the step count.
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE reallocate_ft_traces(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_ft_traces'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(rtbse_env%moments_trace)) DEALLOCATE (rtbse_env%moments_trace)
      IF (ASSOCIATED(rtbse_env%field_trace)) DEALLOCATE (rtbse_env%field_trace)
      IF (ASSOCIATED(rtbse_env%time_trace)) DEALLOCATE (rtbse_env%time_trace)

      ALLOCATE (rtbse_env%moments_trace(rtbse_env%n_spin, 3, rtbse_env%sim_nsteps + 1), &
                source=CMPLX(0.0_dp, 0.0_dp, kind=dp))
      ALLOCATE (rtbse_env%field_trace(3, rtbse_env%sim_nsteps + 1), &
                source=CMPLX(0.0_dp, 0.0_dp, kind=dp))
      ALLOCATE (rtbse_env%time_trace(rtbse_env%sim_nsteps + 1), source=0.0_dp)

      CALL timestop(handle)
   END SUBROUTINE reallocate_ft_traces

   ! **************************************************************************************************
   !> \brief Writes the header and basic info (incl maximum timestep inside stability region)
   !>        to the standard output
   !> \param rtbse_env Entry point - rtbse environment
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param rtbse_env ...
! **************************************************************************************************
   SUBROUTINE print_linrtbse_header_info(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env

      INTEGER                                            :: n_steps
      REAL(kind=dp)                                      :: fft_resolution, nyquist_frequency, &
                                                            total_time
      TYPE(cp_logger_type), POINTER                      :: logger

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
         WRITE (rtbse_env%unit_nr, '(A18,L62)') ' Use Tamm-Dancoff approximation', &
            rtbse_env%tda_active
         IF (rtbse_env%tda_active) THEN
            WRITE (rtbse_env%unit_nr, '(A,T71,L10)') ' TDA first-peak shift active', &
               rtbse_env%tda_shift_to_first_peak
            IF (rtbse_env%tda_shift_to_first_peak) THEN
               WRITE (rtbse_env%unit_nr, '(A,T65,F16.6)') ' TDA first-peak shift Omega_0 [eV]:', &
                  rtbse_env%omega_shift*evolt
            END IF
         END IF

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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'populate_C_active'

      INTEGER                                            :: handle, i

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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_moments'

      INTEGER                                            :: handle, k
      REAL(kind=dp), DIMENSION(3)                        :: rpoint
      TYPE(cp_fm_type)                                   :: tmp_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments_dbcsr_p
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_density_matrix'

      INTEGER                                            :: handle, i, i_row_global, ii, &
                                                            j_col_global, jj, ncol_global, &
                                                            ncol_local, nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_singleparticle_hamiltonian'

      INTEGER                                            :: abs_mo_idx, handle, i, i_row_global, ii, &
                                                            j_col_global, jj, ncol_global, &
                                                            ncol_local, nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

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
                  ! First-peak shift (TDA only): rotate the single-particle Hamiltonian
                  ! into a frame where the lowest active OV mode oscillates at zero, with
                  ! Omega_0 = eps_min_ai. Adds +Omega_0/2 on active occupied diagonals and
                  ! -Omega_0/2 on active virtual diagonals so that [h', rho]_OV picks up an
                  ! overall (- eps_ai + Omega_0) and OO/VV blocks remain commutator-free.
                  IF (rtbse_env%tda_active .AND. rtbse_env%tda_shift_to_first_peak) THEN
                     IF (abs_mo_idx <= rtbse_env%n_occ(i)) THEN
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                           rtbse_env%real_workspace_mo(1)%local_data(ii, jj) + 0.5_dp*rtbse_env%omega_shift
                     ELSE
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                           rtbse_env%real_workspace_mo(1)%local_data(ii, jj) - 0.5_dp*rtbse_env%omega_shift
                     END IF
                  END IF
                  ! Debug-only QP override (paired with the matching hack in
                  ! apply_liouvillian_to_drho). Uncomment BOTH sites together to inject the
                  ! module-level last_column QPs (bse_full_diag.F values) for LRBSE cross-check
                  ! parity. Validated 2026-05-26: TDA + ABBA both within 55 meV of LRBSE.
                  !rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = last_column(i_row_global)/evolt
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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_hartree_potential'

      INTEGER                                            :: handle, i
      LOGICAL                                            :: use_hartree_reference, use_rirs_hartree, &
                                                            use_rirs_sex
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)
      use_hartree_reference = (.NOT. rtbse_env%tda_active) .AND. &
                              (.NOT. rtbse_env%debug_disable_hartree)
      use_rirs_hartree = rtbse_env%rirs_hartree
      use_rirs_sex = rtbse_env%rirs_sex

      ! Make sure the RI-RS V_grid kernel is available if we need it for Hartree.
      IF (use_rirs_hartree) CALL rt_bse_ri_rs_ensure_V_grid(bs_env, rtbse_env%qs_env)

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! v_dbcsr is needed by either AO-RI Hartree (here) or by AO-RI SEX (W = V + W^c assembly).
      IF ((.NOT. use_rirs_hartree) .OR. (.NOT. use_rirs_sex)) THEN
         CALL init_hartree(rtbse_env, rtbse_env%v_dbcsr)
      END IF
      ! Always zero ham_reference here (this routine is the first to touch it).
      ! The Hartree reference subtraction is then conditionally added below.
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_set_all(rtbse_env%ham_reference(i), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      END DO
      ! Calculate the original Hartree potential
      ! Uses rho_orig - same as rho for initial run but different for continued run
      ! In TDA the propagator evaluates separate OV/VO-projected kernel passes.
      ! rho_orig is OO-diagonal in the MO basis, so its OV/VO projections vanish
      ! and the corresponding Hartree reference is identically zero.
      IF (use_hartree_reference) THEN
         DO i = 1, rtbse_env%n_spin
            IF (use_rirs_hartree) THEN
               ! AO-RI get_hartree uses only Re(rho); mirror that on the RI-RS path.
               CALL cp_cfm_to_fm(msource=rtbse_env%rho_ao_lrrtbse(i), &
                                 mtargetr=rtbse_env%real_workspace(1))
               CALL compute_hartree_ri_rs(bs_env, rtbse_env%real_workspace(1), &
                                          rtbse_env%hartree_curr_ao(i))
            ELSE
               CALL get_hartree(rtbse_env, rtbse_env%rho_ao_lrrtbse(i), rtbse_env%hartree_curr_ao(i))
            END IF
            ! Scaling by spin degeneracy
            CALL cp_fm_scale(rtbse_env%spin_degeneracy, rtbse_env%hartree_curr_ao(i))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%hartree_curr_ao(i), rtbse_env%hartree_curr(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%hartree_curr(i), i)
            ! Subtract the reference from the reference Hamiltonian
            ! following H_eff = ... + V_Hartree(rho) - V_Hartree(rho_0),
            CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(i), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END DO
      END IF
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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_cohsex_selfenergy'

      INTEGER                                            :: handle, i
      LOGICAL                                            :: use_coh_reference, use_rirs_sex, &
                                                            use_sex_reference

      CALL timeset(routineN, handle)
      use_coh_reference = .NOT. rtbse_env%debug_disable_coh
      use_sex_reference = (.NOT. rtbse_env%tda_active) .AND. &
                          (.NOT. rtbse_env%debug_disable_sex)
      use_rirs_sex = rtbse_env%rirs_sex

      ! Make sure the RI-RS W0_grid kernel is available if we need it for COH/SEX.
      IF (use_rirs_sex) CALL rt_bse_ri_rs_ensure_W0_grid(rtbse_env%bs_env, rtbse_env%qs_env)

      ! ****** START COHSEX REFERENCE CALCULATION
      ! w_dbcsr and screened_dbt are needed for get_sigma routines (AO-RI path only).
      ! For RI-RS the W = V + W^c kernel is on the real-space grid in mat_W0_grid_rtbse.
      IF (.NOT. use_rirs_sex) THEN
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
      END IF
      ! Calculate the COHSEX starting energies
      DO i = 1, rtbse_env%n_spin
         ! COH is only meaningful with screening (g0w0 path); KS reference has W^c = 0.
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0 .AND. use_coh_reference) THEN
            ! COH is currently a static no-op: this reference subtraction is cancelled
            ! exactly by the per-step addition in update_effective_ham_MO (both build
            ! Sigma_COH from S_inv_fm, no rho dependence). DEBUG_DISABLE_COH is inert.
            ! [INTERNAL] dev_notes/rt_bse_open_questions.md Q6
            CPWARN("linRTBSE COH is currently a static no-op")
            ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
            ! TODO : Allow no COH calculation for static screening
            IF (use_rirs_sex) THEN
               CALL compute_sigma_ri_rs(rtbse_env%bs_env, rtbse_env%sigma_COH_ao(i), -0.5_dp, &
                                        rtbse_env%S_inv_fm)
            ELSE
               CALL get_sigma(rtbse_env, rtbse_env%sigma_COH_ao(i), -0.5_dp, rtbse_env%S_inv_fm)
            END IF
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%sigma_COH_ao(i), rtbse_env%sigma_COH(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%sigma_COH(i), i)
            ! Copy and subtract from the complex reference hamiltonian
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(i), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END IF
         ! Calculate exchange part - TODO : should this be applied for different spins? - TEST with O2 HF propagation?
         ! So far only closed shell tested
         ! Uses rho_orig - same as rho for initial run but different for continued run
         ! For KS reference this is the time-dependent Fock exchange (w_dbcsr = v only).
         ! In TDA the propagator evaluates separate OV/VO-projected kernel passes.
         ! rho_orig is OO-diagonal in the MO basis, so its OV/VO projections
         ! vanish and the SEX reference must remain zero in TDA.
         IF (use_sex_reference) THEN
            IF (use_rirs_sex) THEN
               CALL compute_sigma_ri_rs_complex(rtbse_env%bs_env, rtbse_env%sigma_SEX_ao(i), -1.0_dp, &
                                                rtbse_env%rho_ao_lrrtbse(i))
            ELSE
               CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(i), -1.0_dp, rtbse_env%rho_ao_lrrtbse(i))
            END IF
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(i), rtbse_env%sigma_SEX(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(i), i)
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END IF
      END DO
      ! ****** END COHSEX REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_cohsex_selfenergy

! **************************************************************************************************
!> \brief Propagates the density to the next timestep using Runge-Kutta 4th order (RK4) method
!> \param rtbse_env ...
!> \param rho_start Initial density matrix
!> \param rho_end Final density matrix
!> \param converged Whether the resulting rho_end is self-consistent
!> \param k How many SC iterations were done
!> \param metric The difference metric from the last self-consistent iteration (for printing/evaluation)
! **************************************************************************************************
   SUBROUTINE solve_rk4_timestep(rtbse_env, rho_start, rho_end, converged, k, metric)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_start, rho_end
      LOGICAL                                            :: converged
      INTEGER                                            :: k
      REAL(kind=dp)                                      :: metric

      CHARACTER(len=*), PARAMETER :: routineN = 'solve_rk4_timestep'

      INTEGER                                            :: handle, i

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
         IF (rtbse_env%tda_active) &
            CALL project_drho_to_ov(rtbse_env, rtbse_env%rk4_coefficients(1), i)
         ! MG TODO: dt/2 needs to implemented in field
         !          Problem: sim_step is needed as integer for envelope
         ! k_2 - needs to be evaluated at t + dt/2 and rho + dt/2 * k_1
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(0.5_dp*rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(1))
         CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(2), i)
         IF (rtbse_env%tda_active) &
            CALL project_drho_to_ov(rtbse_env, rtbse_env%rk4_coefficients(2), i)
         ! k_3 - needs to be evaluated at t + dt/2 and rho + dt/2 * k_2
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(0.5_dp*rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(2))
         CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(3), i)
         IF (rtbse_env%tda_active) &
            CALL project_drho_to_ov(rtbse_env, rtbse_env%rk4_coefficients(3), i)
         ! k_4 - needs to be evaluated at t + dt and rho + dt * k_3
         CALL cp_cfm_to_cfm(rho_start(i), rtbse_env%rho_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%rho_workspace(1), &
                                   CMPLX(rtbse_env%sim_dt, 0.0, kind=dp), rtbse_env%rk4_coefficients(3))
         CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho_workspace(1), rtbse_env%rk4_coefficients(4), i)
         IF (rtbse_env%tda_active) &
            CALL project_drho_to_ov(rtbse_env, rtbse_env%rk4_coefficients(4), i)
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
!> \param rho Real and imaginary parts ( + spin) of the density at current time
!> \param ham_effective Effective Hamiltonian in the MO basis that is updated in this routine
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE update_effective_ham_MO(rtbse_env, rho, ham_effective, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho, ham_effective
      INTEGER                                            :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'update_effective_ham_MO'

      INTEGER                                            :: handle, k
      LOGICAL                                            :: use_coh, use_hartree, use_rirs_hartree, &
                                                            use_rirs_sex, use_sex

      CALL timeset(routineN, handle)
      use_coh = .NOT. rtbse_env%debug_disable_coh
      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex
      use_rirs_hartree = rtbse_env%rirs_hartree
      use_rirs_sex = rtbse_env%rirs_sex

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
                       rho, CMPLX(1.0, 0.0, kind=dp), ham_effective)
      CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                       CMPLX(-1.0, 0.0, kind=dp), rho, rtbse_env%ham_reference_singleparticle(ispin), &
                       CMPLX(1.0, 0.0, kind=dp), ham_effective)
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
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0 .AND. use_coh) THEN
         ! Add the COH part - so far static but can be dynamic in principle through the W updates
         IF (use_rirs_sex) THEN
            CALL compute_sigma_ri_rs(rtbse_env%bs_env, rtbse_env%sigma_COH_ao(ispin), -0.5_dp, &
                                     rtbse_env%S_inv_fm)
         ELSE
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH_ao(ispin), -0.5_dp, rtbse_env%S_inv_fm)
         END IF
         ! Transform to MO basis (AO scratch -> MO-active result)
         CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%sigma_COH_ao(ispin), rtbse_env%sigma_COH(ispin), ispin)
         ! Apply occupation factor f_n-f_m
         CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%sigma_COH(ispin), ispin)
         CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(ispin), mtarget=rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                   CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
      END IF
      IF (.NOT. rtbse_env%tda_active) THEN
         ! ----- Standard (non-TDA) kernel: single AO density built once at routine entry -----
         ! Calculate the (S)EX part - based on provided rho
         ! iGW = - rho W
         IF (use_sex) THEN
            IF (use_rirs_sex) THEN
               CALL compute_sigma_ri_rs_complex(rtbse_env%bs_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, &
                                                rtbse_env%rho_ao_lrrtbse(ispin))
            ELSE
               CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, rtbse_env%rho_ao_lrrtbse(ispin))
            END IF
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                  rtbse_env%sigma_SEX(ispin), ispin)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(ispin), ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                      CMPLX(1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(ispin))
         END IF
         ! Calculate Hartree potential
         ! Hartree potential is scaled by number of electrons in each MO - spin degeneracy
         IF (use_hartree) THEN
            IF (use_rirs_hartree) THEN
               ! AO-RI get_hartree uses only Re(rho); mirror that on the RI-RS path.
               CALL cp_cfm_to_fm(msource=rtbse_env%rho_ao_lrrtbse(ispin), &
                                 mtargetr=rtbse_env%real_workspace(1))
               CALL compute_hartree_ri_rs(rtbse_env%bs_env, rtbse_env%real_workspace(1), &
                                          rtbse_env%hartree_curr_ao(ispin))
            ELSE
               CALL get_hartree(rtbse_env, rtbse_env%rho_ao_lrrtbse(ispin), &
                                rtbse_env%hartree_curr_ao(ispin))
            END IF
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, rtbse_env%hartree_curr_ao(ispin), &
                                                 rtbse_env%hartree_curr(ispin), ispin)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%hartree_curr(ispin), ispin)
            CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(ispin), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), ham_effective, &
                                      CMPLX(rtbse_env%spin_degeneracy, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END IF
      ELSE
         ! ----- Tamm-Dancoff approximation: drop B-coupling between OV and VO
         ! blocks. Single OV pass + Hermitian conjugate for VO.
         !   Delta rho_VO = (Delta rho_OV)^C (Delta rho is Hermitian).
         !   The (mu nu|lambda sigma) and W kernels are symmetric in their AO-pair
         !   indices and C_active is real, so the AO<->MO transforms commute with
         !   conjugate transpose: K_MO[Delta rho_VO] = (K_MO[Delta rho_OV])^C for
         !   both K = V_H and K = Sigma_SEX. We therefore evaluate each kernel only
         !   once on Delta rho_OV, mask the MO result to OV, and add the VO-block
         !   contribution as the conjugate transpose of the OV result.
         !   rho_orig is OO-diagonal in the MO basis, so masking rho gives the same
         !   OV input as masking rho - rho_orig (TDA reference remains 0 by
         !   construction).
         !
         ! Kernel inputs are non-Hermitian (Delta rho_OV alone has no VO entries),
         ! so the kernels MUST be complex-linear: get_sigma_complex (Re/Im handled
         ! separately) and get_hartree_complex (Re/Im fed to real-input get_hartree
         ! separately, then reassembled). Real-input get_hartree on a non-Hermitian
         ! input would silently drop Im and break Hermitian conjugacy of the
         ! OV+VO contributions to ham_effective, giving the propagator real
         ! eigenvalues -> exponential blow-up.
         !
         ! Hartree runs before SEX: sigma_SEX_ao(ispin) is reused as the AO
         ! scratch holding the complex Hartree result, then overwritten by
         ! get_sigma for the SEX contribution.
         !
         ! Signs: (f_n - f_m) is -1 on OV and +1 on VO; applied explicitly.

         ! Build Delta rho_OV (in MO) and project to AO
         CALL cp_cfm_to_cfm(rho, rtbse_env%rho_delta_mo(ispin))
         CALL mask_mo_block_cfm(rtbse_env, rtbse_env%rho_delta_mo(ispin), ispin, keep_OV=.TRUE.)
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_delta_mo(ispin), &
                                                   rtbse_env%rho_ao_lrrtbse(ispin), ispin)

         ! Hartree must run before SEX: sigma_SEX_ao(ispin) is reused as the AO
         ! scratch for the complex Hartree result, then overwritten by get_sigma.
         IF (use_hartree) THEN
            IF (use_rirs_hartree) THEN
               CALL compute_hartree_ri_rs_complex(rtbse_env%bs_env, rtbse_env%rho_ao_lrrtbse(ispin), &
                                                  rtbse_env%sigma_SEX_ao(ispin))
            ELSE
               CALL get_hartree_complex(rtbse_env, rtbse_env%rho_ao_lrrtbse(ispin), &
                                        rtbse_env%sigma_SEX_ao(ispin), ispin)
            END IF
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                  rtbse_env%ham_workspace(1), ispin)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%ham_workspace(1), ispin, keep_OV=.TRUE.)
            ! OV contribution: sign -1 * spin_degeneracy
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(-rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), &
                                      rtbse_env%ham_workspace(1))
            ! VO contribution: (OV result)^C, sign +1 * spin_degeneracy.
            ! Use rho_delta_mo(ispin) as scratch (idle since the AO projection above).
            CALL cp_cfm_transpose(rtbse_env%ham_workspace(1), 'C', rtbse_env%rho_delta_mo(ispin))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), &
                                      rtbse_env%rho_delta_mo(ispin))
         END IF
         IF (use_sex) THEN
            IF (use_rirs_sex) THEN
               CALL compute_sigma_ri_rs_complex(rtbse_env%bs_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                -1.0_dp, rtbse_env%rho_ao_lrrtbse(ispin))
            ELSE
               CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, rtbse_env%rho_ao_lrrtbse(ispin))
            END IF
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                  rtbse_env%sigma_SEX(ispin), ispin)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%sigma_SEX(ispin), ispin, keep_OV=.TRUE.)
            ! OV contribution: sign -1
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%sigma_SEX(ispin))
            ! VO contribution: (OV result)^C, sign +1
            CALL cp_cfm_transpose(rtbse_env%sigma_SEX(ispin), 'C', rtbse_env%rho_delta_mo(ispin))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_delta_mo(ispin))
         END IF

         ! Restore rho_ao_lrrtbse to the AO image of the current full rho, so any
         ! post-routine consumer (output_mos, COH refresh) sees the same invariant
         ! as on the non-TDA path.
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rho, &
                                                   rtbse_env%rho_ao_lrrtbse(ispin), ispin)
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
!> \brief Apply the linearized RT-BSE Liouvillian (TDA) to an input Delta rho_MO and return
!>        L * Delta rho_MO in MO basis, OV+VO blocks populated. Detached from the propagator
!>        (does NOT touch rho / rho_orig / ham_effective / ham_reference); only env scratches
!>        drho_probe, L_drho, rho_ao_lrrtbse, sigma_SEX_ao, ham_workspace(1), sigma_SEX,
!>        real_workspace_mo(1) are used.
!>
!>        Matrix elements correspond to the Casida-A matrix:
!>            L_{ia,jb} = (eps_a - eps_i) delta_{ij} delta_{ab} + (ia|jb) - W_{ij,ab}
!>        with no (f_n - f_m) factor (the propagator path applies -1 on OV; we want raw +K).
!> [INTERNAL] dev_notes/rt_bse_equations.md §8 (derivation of L_TDA = A_Casida)
!>
!> \param rtbse_env RT-BSE environment (TDA, n_spin=1).
!> \param drho_in Input Delta rho_MO (mo_active x mo_active complex).
!> \param L_drho_out Output L * Delta rho_MO; OV+VO blocks populated; OO/VV zero.
!> \param ispin Spin index.
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE apply_liouvillian_to_drho(rtbse_env, drho_in, L_drho_out, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: drho_in
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: L_drho_out
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_liouvillian_to_drho'

      INTEGER                                            :: abs_mo_idx, handle, i_row_global, ii, &
                                                            j_col_global, jj, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: use_hartree, use_sex

      CALL timeset(routineN, handle)

      ! Mirror update_effective_ham_MO:1048-1050. No use_coh: COH has no LR action
      ! (Q6) and is not applied in this routine; only Hartree + SX are gated.
      ! [INTERNAL] dev_notes/rt_bse_open_questions.md Q5
      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex

      ! 1. Stage drho_in into drho_probe. TDA: mask to OV (propagator carries only OV
      !    and the VO contribution is recovered later via the (.)^C shortcut). ABBA: keep
      !    the full OV+VO content - drho_in carries both blocks independently.
      ! [INTERNAL] dev_notes/rt_bse_equations.md §8.5.5
      CALL cp_cfm_to_cfm(drho_in, rtbse_env%drho_probe)
      IF (rtbse_env%tda_active) THEN
         CALL mask_mo_block_cfm(rtbse_env, rtbse_env%drho_probe, ispin, keep_OV=.TRUE.)
      END IF

      ! 2. Project Delta rho_OV (MO -> AO, contravariant).
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%drho_probe, &
                                                rtbse_env%rho_ao_lrrtbse(ispin), ispin)

      ! 3. Initialize the output accumulator.
      CALL cp_cfm_set_all(L_drho_out, CMPLX(0.0_dp, 0.0_dp, kind=dp))

      ! 4. eps_{ai} diagonal contribution. MUST run BEFORE the kernel steps (5/6) below:
      !    the kernel path scratches drho_probe via cp_cfm_transpose for the VO Hermitian
      !    conjugate, and drho_in IS rtbse_env%drho_probe (aliased at the caller) - reading
      !    drho_in for the [drho, H_eps] commutator after the kernel ran corrupts the eps
      !    part. Build H_eps as a real diagonal fm in real_workspace_mo(1), convert to cfm
      !    in ham_workspace(1), accumulate [drho_in, H_eps] = drho * H - H * drho.
      !    On OV: ([drho, H])_{ia} = (eps_a - eps_i) * drho_{ia} = +eps_{ai} * drho_{ia}.
      !    On VO: -eps_{ai} * drho_{ai} (sign flips); irrelevant - driver only reads OV.
      !    Bare GW eigenvalues (lab frame) so eigenvalues compare 1:1 to bse_full.F.
      CALL cp_fm_get_info(rtbse_env%real_workspace_mo(1), &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
      DO ii = 1, nrow_local
         i_row_global = row_indices(ii)
         DO jj = 1, ncol_local
            j_col_global = col_indices(jj)
            IF (i_row_global == j_col_global) THEN
               abs_mo_idx = i_row_global + rtbse_env%first_active_mo - 1
               rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                  rtbse_env%bs_env%eigenval_G0W0(abs_mo_idx, 1, ispin)
               ! Debug-only QP override (paired with the matching hack in
               ! initialize_singleparticle_hamiltonian). Uncomment BOTH sites together to inject
               ! the module-level last_column QPs for LRBSE cross-check parity.
               !rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = last_column(i_row_global)/evolt
            END IF
         END DO
      END DO
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), &
                        mtarget=rtbse_env%ham_workspace(1))
      ! drho * H_eps
      CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                       CMPLX(1.0_dp, 0.0_dp, kind=dp), drho_in, rtbse_env%ham_workspace(1), &
                       CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out)
      ! -H_eps * drho
      CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                       CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1), drho_in, &
                       CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out)

      ! 5. Hartree contribution: get_hartree_complex (Re/Im split) on AO Delta rho, then
      !    AO->MO covariant, mask OV+VO, accumulate +spin_degeneracy on OV and on VO.
      !    No (f_n - f_m) factor: we want raw +K (Casida convention), not the propagator's -K_OV.
      !    VO contribution comes from the Hermitian conjugate of the OV result (P2 in
      !    rt_bse_pitfalls_physics.md): for real C_active + AO-pair-symmetric kernel,
      !    K_MO[Delta rho_VO] = (K_MO[Delta rho_OV])^C.
      IF (use_hartree) THEN
         CALL get_hartree_complex(rtbse_env, rtbse_env%rho_ao_lrrtbse(ispin), &
                                  rtbse_env%sigma_SEX_ao(ispin), ispin)
         CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                               rtbse_env%ham_workspace(1), ispin)
         CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%ham_workspace(1), L_drho_out, &
                                 CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), ispin)
      END IF

      ! 6. Screened-exchange contribution: get_sigma is already complex-aware. The
      !    -1.0_dp factor passed to get_sigma builds the -W contribution; we then add
      !    +1.0 here (no occupation-factor flip), giving raw K^SX = -W as required by
      !    K = (ia|jb) - W_{ij,ab} (Casida convention with +K).
      IF (use_sex) THEN
         CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, rtbse_env%rho_ao_lrrtbse(ispin))
         CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                               rtbse_env%sigma_SEX(ispin), ispin)
         CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%sigma_SEX(ispin), L_drho_out, &
                                 CMPLX(1.0_dp, 0.0_dp, kind=dp), ispin)
      END IF

      CALL timestop(handle)
   END SUBROUTINE apply_liouvillian_to_drho

! **************************************************************************************************
!> \brief Adds the MO-domain kernel contribution K_MO * Delta rho to L_drho_out, branching on
!>        rtbse_env%tda_active. TDA: masks K_MO to OV in place, adds scale*K_MO on OV, then
!>        builds the VO contribution as (K_MO_OV)^C and adds it (valid only under the TDA
!>        assumption Delta rho_VO = (Delta rho_OV)^C). ABBA: adds the full mo_active x mo_active
!>        K_MO directly - OV+VO blocks are computed naturally from the full OV+VO input.
!>        K_MO is INTENT(INOUT); in the TDA branch it is masked in place (treat as scratch
!>        after this call). drho_probe is also scratched in the TDA branch (the caller's eps
!>        step must already have consumed drho_in).
!> [INTERNAL] dev_notes/rt_bse_equations.md §8.5.5
!> \param rtbse_env RT-BSE environment.
!> \param K_MO Kernel contribution in MO basis (mo_active x mo_active). Scratched in TDA branch.
!> \param L_drho_out Accumulator (mo_active x mo_active).
!> \param scale Complex scale factor (spin_degeneracy for Hartree, 1.0 for SX).
!> \param ispin Spin index.
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE add_K_MO_to_L_drho(rtbse_env, K_MO, L_drho_out, scale, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: K_MO, L_drho_out
      COMPLEX(kind=dp), INTENT(IN)                       :: scale
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'add_K_MO_to_L_drho'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (rtbse_env%tda_active) THEN
         CALL mask_mo_block_cfm(rtbse_env, K_MO, ispin, keep_OV=.TRUE.)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out, scale, K_MO)
         CALL cp_cfm_transpose(K_MO, 'C', rtbse_env%drho_probe)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out, scale, rtbse_env%drho_probe)
      ELSE
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out, scale, K_MO)
      END IF

      CALL timestop(handle)
   END SUBROUTINE add_K_MO_to_L_drho

! **************************************************************************************************
!> \brief Public entry for the Liouvillian eigenvalue diagnostic (observable O-D10).
!>        Dispatches to the TDA branch (Casida-A via cp_cfm_heevd) or the ABBA branch
!>        (Furche reduction via cp_cfm_power) based on rtbse_env%tda_active.
!>        Called once at job init from run_propagation_linearized_bse, gated on
!>        rtbse_env%diagnose_liouvillian_eig (n_spin = 1 enforced at env creation).
!> [INTERNAL] dev_notes/rt_bse_equations.md §8.5; dev_notes/rt_bse_observables.md O-D10
!> \param rtbse_env RT-BSE environment with diagnostic scratch already allocated.
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE diagnose_liouvillian_eigenvalues(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagnose_liouvillian_eigenvalues'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (rtbse_env%tda_active) THEN
         CALL diagnose_TDA_liouvillian(rtbse_env)
      ELSE
         CALL diagnose_ABBA_liouvillian(rtbse_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE diagnose_liouvillian_eigenvalues

! **************************************************************************************************
!> \brief TDA branch of the Liouvillian eigenvalue diagnostic (O-D10). Probes
!>        apply_liouvillian_to_drho with canonical OV basis vectors, assembles the Casida-A
!>        matrix as L_pairs (N_OV x N_OV), and diagonalizes via cp_cfm_heevd. Eigenvalues are
!>        written to stdout (line prefix RTBSE|) and to the LIOUVILLIAN_EIG print key file.
!>        Iterates over rtbse_env%n_spin channels (n_spin = 1 enforced at env creation).
!>        Detached from RK4 state. Called from the dispatcher diagnose_liouvillian_eigenvalues
!>        when tda_active=.TRUE.
!> [INTERNAL] dev_notes/rt_bse_equations.md §8.2 (Casida-A identification)
!> \param rtbse_env RT-BSE environment with TDA diagnostic scratch already allocated.
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE diagnose_TDA_liouvillian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagnose_TDA_liouvillian'

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: ov_block
      INTEGER                                            :: b, eig_unit, handle, ispin, j, k_col, n, &
                                                            n_act_occ, n_act_virt, n_ov
      REAL(kind=dp)                                      :: residual_max
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      ! Banner (per-routine; spin sections follow below).
      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE| ----- TDA Liouvillian diagnostic -----'
      END IF

      ! File output: open once (REWIND), append per-spin sections, close after the loop.
      ! ignore_should_output=.TRUE. so this fires at init regardless of MD-iteration cadence.
      eig_unit = cp_print_key_unit_nr(logger, rtbse_env%eig_section, &
                                      extension=".dat", &
                                      file_form="FORMATTED", &
                                      file_position="REWIND", &
                                      ignore_should_output=.TRUE.)
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '# TDA Liouvillian eigenvalues (one section per spin)'
      END IF

      ! Per-spin loop. n_spin = 1 enforced at env creation today; structure is
      ! multi-spin-ready for when that restriction lifts.
      DO ispin = 1, rtbse_env%n_spin
         n_act_occ = rtbse_env%n_occ(ispin) - rtbse_env%first_active_mo + 1
         n_act_virt = rtbse_env%last_active_mo - rtbse_env%n_occ(ispin)
         n_ov = n_act_occ*n_act_virt
         ALLOCATE (ov_block(n_act_occ, n_act_virt))

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,I0)') ' RTBSE| Spin channel ', ispin
            WRITE (rtbse_env%unit_nr, '(A,I0,A,I0,A,I0)') &
               ' RTBSE| n_act_occ = ', n_act_occ, &
               ', n_act_virt = ', n_act_virt, &
               ', N_OV = ', n_ov
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
         END IF

         ! Column-major OV index: k = (b_local - 1) * n_act_occ + j_local
         !   j_local = j - first_active_mo + 1  (1 .. n_act_occ)
         !   b_local = b - n_occ(ispin)         (1 .. n_act_virt)
         ! Probe with each canonical unit vector e_{(j,b)} and read out the OV block of L_drho
         ! as the k_col-th column of L_pairs. Block-level transfers (cp_cfm_{get,set}_submatrix)
         ! keep this at O(N_OV) collective ops total instead of O(N_OV^2) per-element ones.
         ! The Fortran column-major layout of ov_block (n_act_occ, n_act_virt) matches the
         ! column-major OV index by construction, so RESHAPE without padding gives the right
         ! (N_OV, 1) column for L_pairs.
         DO b = rtbse_env%n_occ(ispin) + 1, rtbse_env%last_active_mo
            DO j = rtbse_env%first_active_mo, rtbse_env%n_occ(ispin)
               k_col = (b - rtbse_env%n_occ(ispin) - 1)*n_act_occ + (j - rtbse_env%first_active_mo + 1)

               CALL cp_cfm_set_all(rtbse_env%drho_probe, CMPLX(0.0_dp, 0.0_dp, kind=dp))
               CALL cp_cfm_set_element(rtbse_env%drho_probe, &
                                       j - rtbse_env%first_active_mo + 1, &
                                       b - rtbse_env%first_active_mo + 1, &
                                       CMPLX(1.0_dp, 0.0_dp, kind=dp))

               CALL apply_liouvillian_to_drho(rtbse_env, rtbse_env%drho_probe, rtbse_env%L_drho, ispin)

               ! Extract OV block of L_drho: rows 1..n_act_occ (i_local), cols n_act_occ+1..mo_active (a_local).
               CALL cp_cfm_get_submatrix(rtbse_env%L_drho, ov_block, &
                                         start_row=1, start_col=n_act_occ + 1, &
                                         n_rows=n_act_occ, n_cols=n_act_virt)

               ! Pack OV block into column k_col of L_pairs (column-major reshape).
               CALL cp_cfm_set_submatrix(rtbse_env%L_pairs, &
                                         RESHAPE(ov_block, [n_ov, 1]), &
                                         start_row=1, start_col=k_col, &
                                         n_rows=n_ov, n_cols=1)
            END DO
         END DO

         ! Hermitian residual: real-orbital BSE => L is real symmetric, so ||L - L^H||_max
         ! should be at the floating-point floor (~1e-12). Compute the difference matrix
         ! D = L^H - L into eigvecs_pairs (idle here; overwritten by cp_cfm_heevd below),
         ! then take its max-element norm via the BLACS-native pzlange path - O(N_OV^2)
         ! local work + one MPI reduce, instead of 2*N_OV^2 collective get_element calls.
         CALL cp_cfm_transpose(rtbse_env%L_pairs, 'C', rtbse_env%eigvecs_pairs)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%eigvecs_pairs, &
                                   CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%L_pairs)
         residual_max = cp_cfm_norm(rtbse_env%eigvecs_pairs, 'M')
         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,ES16.6)') &
               ' RTBSE| Hermitian residual ||L - L^H||_max = ', residual_max
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
         END IF
         IF (residual_max > 1.0E-6_dp) THEN
            CPABORT("Liouvillian Hermitian residual > 1e-6 - check kernel signs / symmetry.")
         END IF

         ! Diagonalize. cp_cfm_heevd returns real eigenvalues + complex eigenvectors.
         CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, &
                           rtbse_env%eigenvalues_liouvillian)

         ! Stdout table for this spin. Atomic-units Omega goes to the .dat file below;
         ! stdout shows only eV in F12.4 right-aligned to col 80 (L-7 convention).
         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,T26,A,T59,A)') &
               ' RTBSE|', "Excitation index n", "Excitation energy (eV)"
            DO n = 1, n_ov
               WRITE (rtbse_env%unit_nr, '(A,T40,I4,T69,F12.4)') &
                  ' RTBSE|', n, rtbse_env%eigenvalues_liouvillian(n)*evolt
            END DO
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
         END IF

         ! File output for this spin (a.u. + eV).
         IF (eig_unit > 0) THEN
            WRITE (eig_unit, '(/,A,I0)') '# Spin channel ', ispin
            WRITE (eig_unit, '(A)') '#   n          Omega [a.u.]               Omega [eV]'
            DO n = 1, n_ov
               WRITE (eig_unit, '(I5,4X,ES24.14E3,4X,ES24.14E3)') n, &
                  rtbse_env%eigenvalues_liouvillian(n), &
                  rtbse_env%eigenvalues_liouvillian(n)*evolt
            END DO
         END IF

         DEALLOCATE (ov_block)
      END DO

      CALL cp_print_key_finished_output(eig_unit, logger, rtbse_env%eig_section)

      CALL timestop(handle)
   END SUBROUTINE diagnose_TDA_liouvillian

! **************************************************************************************************
!> \brief ABBA branch of the Liouvillian eigenvalue diagnostic. Probes apply_liouvillian_to_drho
!>        with canonical OV unit vectors, reads OV -> A and VO -> -B^* (one probe pass), recovers
!>        B by sign-flip + complex conjugation, forms (A+/-B), enforces (A-B) > 0 via an
!>        independent cp_cfm_heevd pre-check, builds (A-B)^{1/2} via cp_cfm_power, forms
!>        C = (A-B)^{1/2}(A+B)(A-B)^{1/2}, diagonalizes via cp_cfm_heevd, takes Omega_n = sqrt(C).
!>        Per-spin loop; n_spin=1 enforced at env creation. Output mirrors the TDA branch:
!>        stdout table prefixed RTBSE| plus a .dat section per spin in the LIOUVILLIAN_EIG file.
!>        Residual + lambda_min(A-B) diagnostics go to stdout only (matches TDA pattern).
!> [INTERNAL] dev_notes/rt_bse_equations.md §8.5
!> \param rtbse_env RT-BSE environment with ABBA diagnostic scratch allocated.
!> \author Maximilian Graml (05/26)
! **************************************************************************************************
   SUBROUTINE diagnose_ABBA_liouvillian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagnose_ABBA_liouvillian'

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: ov_block, vo_block
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: B_local
      INTEGER                                            :: b, eig_unit, handle, ispin, j, k_col, n, &
                                                            n_act_occ, n_act_virt, n_ov
      REAL(kind=dp)                                      :: lambda_min_AmB, residual_A, residual_B
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE| ----- ABBA Liouvillian diagnostic -----'
      END IF

      eig_unit = cp_print_key_unit_nr(logger, rtbse_env%eig_section, &
                                      extension=".dat", &
                                      file_form="FORMATTED", &
                                      file_position="REWIND", &
                                      ignore_should_output=.TRUE.)
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '# ABBA Liouvillian eigenvalues (one section per spin)'
      END IF

      DO ispin = 1, rtbse_env%n_spin
         n_act_occ = rtbse_env%n_occ(ispin) - rtbse_env%first_active_mo + 1
         n_act_virt = rtbse_env%last_active_mo - rtbse_env%n_occ(ispin)
         n_ov = n_act_occ*n_act_virt
         ALLOCATE (ov_block(n_act_occ, n_act_virt))
         ALLOCATE (vo_block(n_act_virt, n_act_occ))

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,I0)') ' RTBSE| Spin channel ', ispin
            WRITE (rtbse_env%unit_nr, '(A,I0,A,I0,A,I0)') &
               ' RTBSE| n_act_occ = ', n_act_occ, &
               ', n_act_virt = ', n_act_virt, &
               ', N_OV = ', n_ov
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
         END IF

         ! Probe with e_{ia}^OV; OV block of output -> column k_col of A; VO block -> column k_col
         ! of -B^* (recovery one line below). Column-major OV index matches the TDA probe.
         DO b = rtbse_env%n_occ(ispin) + 1, rtbse_env%last_active_mo
            DO j = rtbse_env%first_active_mo, rtbse_env%n_occ(ispin)
               k_col = (b - rtbse_env%n_occ(ispin) - 1)*n_act_occ + (j - rtbse_env%first_active_mo + 1)

               CALL cp_cfm_set_all(rtbse_env%drho_probe, CMPLX(0.0_dp, 0.0_dp, kind=dp))
               CALL cp_cfm_set_element(rtbse_env%drho_probe, &
                                       j - rtbse_env%first_active_mo + 1, &
                                       b - rtbse_env%first_active_mo + 1, &
                                       CMPLX(1.0_dp, 0.0_dp, kind=dp))

               CALL apply_liouvillian_to_drho(rtbse_env, rtbse_env%drho_probe, rtbse_env%L_drho, ispin)

               CALL cp_cfm_get_submatrix(rtbse_env%L_drho, ov_block, &
                                         start_row=1, start_col=n_act_occ + 1, &
                                         n_rows=n_act_occ, n_cols=n_act_virt)
               CALL cp_cfm_set_submatrix(rtbse_env%A_mat, &
                                         RESHAPE(ov_block, [n_ov, 1]), &
                                         start_row=1, start_col=k_col, &
                                         n_rows=n_ov, n_cols=1)

               ! TRANSPOSE puts vo_block in (n_act_occ, n_act_virt) layout, matching the OV pack.
               CALL cp_cfm_get_submatrix(rtbse_env%L_drho, vo_block, &
                                         start_row=n_act_occ + 1, start_col=1, &
                                         n_rows=n_act_virt, n_cols=n_act_occ)
               CALL cp_cfm_set_submatrix(rtbse_env%B_mat, &
                                         RESHAPE(TRANSPOSE(vo_block), [n_ov, 1]), &
                                         start_row=1, start_col=k_col, &
                                         n_rows=n_ov, n_cols=1)
            END DO
         END DO

         ! Recover B from -B^* via rank-local pass on the cfm's MPI-local data (documented
         ! exception to the fm/cfm-routines-only rule); B_recovered = -CONJG(stored).
         B_local => rtbse_env%B_mat%local_data
         B_local = -CONJG(B_local)

         ! Block-symmetry residuals: A Hermitian, B real-symmetric (B = B^T under real orbitals).
         CALL cp_cfm_transpose(rtbse_env%A_mat, 'C', rtbse_env%eigvecs_pairs)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%eigvecs_pairs, &
                                   CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%A_mat)
         residual_A = cp_cfm_norm(rtbse_env%eigvecs_pairs, 'M')

         CALL cp_cfm_transpose(rtbse_env%B_mat, 'T', rtbse_env%eigvecs_pairs)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%eigvecs_pairs, &
                                   CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%B_mat)
         residual_B = cp_cfm_norm(rtbse_env%eigvecs_pairs, 'M')

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,ES16.6)') &
               ' RTBSE| Hermitian residual ||A - A^H||_max = ', residual_A
            WRITE (rtbse_env%unit_nr, '(A,ES16.6)') &
               ' RTBSE| Symmetry residual  ||B - B^T||_max = ', residual_B
         END IF
         IF (residual_A > 1.0E-6_dp) THEN
            CPABORT("A is not Hermitian within 1e-6 - check kernel signs / symmetry.")
         END IF
         IF (residual_B > 1.0E-6_dp) THEN
            CPABORT("B is not symmetric within 1e-6 - check kernel signs / symmetry.")
         END IF

         ! A +/- B in scratches (both Hermitian).
         CALL cp_cfm_to_cfm(rtbse_env%A_mat, rtbse_env%AmB_scratch)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%AmB_scratch, &
                                   CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%B_mat)
         CALL cp_cfm_to_cfm(rtbse_env%A_mat, rtbse_env%ApB_scratch)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%ApB_scratch, &
                                   CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%B_mat)

         ! (A-B) positivity gate: heevd on a copy (heevd destroys input); abort if lambda_min < 0.
         CALL cp_cfm_to_cfm(rtbse_env%AmB_scratch, rtbse_env%L_pairs)
         CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, rtbse_env%eigenvalues_liouvillian)
         lambda_min_AmB = rtbse_env%eigenvalues_liouvillian(1)
         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,ES16.6,A,F12.6,A)') &
               ' RTBSE| lambda_min(A - B)                       = ', &
               lambda_min_AmB, '  a.u.  (', lambda_min_AmB*evolt, ' eV)'
         END IF
         IF (lambda_min_AmB < 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "(A - B) not positive definite. Likely causes: "// &
                          "(1) B-readout sign error (see dev_notes/rt_bse_equations.md §8.5.4), "// &
                          "(2) wrong sign on a Coulomb/SX contribution to B, "// &
                          "(3) genuine triplet/charge-transfer instability.")
         END IF

         ! In-place AmB_scratch -> (A-B)^{1/2} via cp_cfm_power. threshold=0 substitution
         ! codepath is unreachable here (lambda_min_AmB > 0 already enforced above).
         CALL cp_cfm_power(rtbse_env%AmB_scratch, threshold=0.0_dp, exponent=0.5_dp)

         ! C = (A-B)^{1/2}(A+B)(A-B)^{1/2}. B_mat free after step "A+/-B" -> reuse as T scratch.
         CALL cp_cfm_gemm('N', 'N', n_ov, n_ov, n_ov, &
                          CMPLX(1.0_dp, 0.0_dp, kind=dp), &
                          rtbse_env%AmB_scratch, rtbse_env%ApB_scratch, &
                          CMPLX(0.0_dp, 0.0_dp, kind=dp), &
                          rtbse_env%B_mat)
         CALL cp_cfm_gemm('N', 'N', n_ov, n_ov, n_ov, &
                          CMPLX(1.0_dp, 0.0_dp, kind=dp), &
                          rtbse_env%B_mat, rtbse_env%AmB_scratch, &
                          CMPLX(0.0_dp, 0.0_dp, kind=dp), &
                          rtbse_env%L_pairs)

         ! Diagonalize C -> Omega_n^2; take +sqrt. Safety clamp on tiny-negative noise.
         CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, rtbse_env%eigenvalues_liouvillian)
         DO n = 1, n_ov
            IF (rtbse_env%eigenvalues_liouvillian(n) < 0.0_dp) THEN
               rtbse_env%eigenvalues_liouvillian(n) = 0.0_dp
            ELSE
               rtbse_env%eigenvalues_liouvillian(n) = SQRT(rtbse_env%eigenvalues_liouvillian(n))
            END IF
         END DO

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
            WRITE (rtbse_env%unit_nr, '(A,T26,A,T59,A)') &
               ' RTBSE|', "Excitation index n", "Excitation energy (eV)"
            DO n = 1, n_ov
               WRITE (rtbse_env%unit_nr, '(A,T40,I4,T69,F12.4)') &
                  ' RTBSE|', n, rtbse_env%eigenvalues_liouvillian(n)*evolt
            END DO
            WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
         END IF
         IF (eig_unit > 0) THEN
            WRITE (eig_unit, '(/,A,I0)') '# Spin channel ', ispin
            WRITE (eig_unit, '(A)') '#   n          Omega [a.u.]               Omega [eV]'
            DO n = 1, n_ov
               WRITE (eig_unit, '(I5,4X,ES24.14E3,4X,ES24.14E3)') n, &
                  rtbse_env%eigenvalues_liouvillian(n), &
                  rtbse_env%eigenvalues_liouvillian(n)*evolt
            END DO
         END IF

         DEALLOCATE (ov_block, vo_block)
      END DO

      CALL cp_print_key_finished_output(eig_unit, logger, rtbse_env%eig_section)

      CALL timestop(handle)
   END SUBROUTINE diagnose_ABBA_liouvillian

! **************************************************************************************************
!> \brief Transforms from AO to MO basis
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param fm_ao ...
!> \param fm_mo ...
!> \param i_spin ...
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_covariant_fm(rtbse_env, fm_ao, fm_mo, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: fm_ao, fm_mo
      INTEGER, INTENT(IN)                                :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_ao_to_mo_covariant_fm'

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
!> \param fm_ao ...
!> \param fm_mo ...
!> \param i_spin ...! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_covariant_cfm(rtbse_env, fm_ao, fm_mo, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: fm_ao, fm_mo
      INTEGER, INTENT(IN)                                :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_ao_to_mo_covariant_cfm'

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
!> \param fm_mo ...
!> \param fm_ao ...
!> \param i_spin ...
! **************************************************************************************************
   SUBROUTINE transform_mo_to_ao_contravariant_cfm(rtbse_env, fm_mo, fm_ao, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: fm_mo, fm_ao
      INTEGER, INTENT(IN)                                :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_mo_to_ao_contravariant_cfm'

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
!> \param cfm ...
!> \param i_spin ...
! **************************************************************************************************
   SUBROUTINE transform_mo_occupation_factor_diff_cfm(rtbse_env, cfm, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: cfm
      INTEGER                                            :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_mo_occupation_factor_diff_cfm'

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

! **************************************************************************************************
!> \brief ...
!> \param rtbse_env ...
!> \param fm ...
!> \param i_spin ...
! **************************************************************************************************
   SUBROUTINE transform_mo_occupation_factor_diff_fm(rtbse_env, fm, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: fm
      INTEGER                                            :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_mo_occupation_factor_diff_fm'

      INTEGER                                            :: handle, i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(kind=dp)                                      :: occ_factor
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data

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
!> \brief Mask an MO-active cfm: keep either OV or VO block, zero everything else.
!> \param rtbse_env RT-BSE environment
!> \param cfm MO-active cfm to mask in place
!> \param i_spin Spin index
!> \param keep_OV .TRUE. keeps the (occ row, virt col) block; .FALSE. keeps (virt row, occ col)
! **************************************************************************************************
   SUBROUTINE mask_mo_block_cfm(rtbse_env, cfm, i_spin, keep_OV)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: cfm
      INTEGER, INTENT(IN)                                :: i_spin
      LOGICAL, INTENT(IN)                                :: keep_OV

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mask_mo_block_cfm'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: handle, i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: keep

      CALL timeset(routineN, handle)

      n_occ = rtbse_env%n_occ(i_spin)
      shift = rtbse_env%first_active_mo - 1

      CALL cp_cfm_get_info(matrix=cfm, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      local_data => cfm%local_data

      DO i_local = 1, nrow_local
         i_global = row_indices(i_local)
         i_global_mo = i_global + shift
         DO j_local = 1, ncol_local
            j_global = col_indices(j_local)
            j_global_mo = j_global + shift
            IF (keep_OV) THEN
               keep = (i_global_mo <= n_occ .AND. j_global_mo > n_occ)
            ELSE
               keep = (i_global_mo > n_occ .AND. j_global_mo <= n_occ)
            END IF
            IF (.NOT. keep) local_data(i_local, j_local) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE mask_mo_block_cfm

! **************************************************************************************************
!> \brief Mask an MO-active fm: keep either OV or VO block, zero everything else.
!> \param rtbse_env RT-BSE environment
!> \param fm MO-active fm to mask in place
!> \param i_spin Spin index
!> \param keep_OV .TRUE. keeps the (occ row, virt col) block; .FALSE. keeps (virt row, occ col)
! **************************************************************************************************
   SUBROUTINE mask_mo_block_fm(rtbse_env, fm, i_spin, keep_OV)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: fm
      INTEGER, INTENT(IN)                                :: i_spin
      LOGICAL, INTENT(IN)                                :: keep_OV

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mask_mo_block_fm'

      INTEGER                                            :: handle, i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: keep
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data

      CALL timeset(routineN, handle)

      n_occ = rtbse_env%n_occ(i_spin)
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
            IF (keep_OV) THEN
               keep = (i_global_mo <= n_occ .AND. j_global_mo > n_occ)
            ELSE
               keep = (i_global_mo > n_occ .AND. j_global_mo <= n_occ)
            END IF
            IF (.NOT. keep) local_data(i_local, j_local) = 0.0_dp
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE mask_mo_block_fm

! **************************************************************************************************
!> \brief Multiply a MO-active cfm rho by the TDA symmetric-shift rotation:
!>        rho_OV *= exp(i * phase),  rho_VO *= exp(-i * phase).
!>        OO/VV blocks are left untouched. For direction='to_lab' pass phase = -Omega_0*t;
!>        for direction='to_rotating' pass phase = +Omega_0*t. No-op when omega_shift = 0.
!>        Preserves Hermiticity since the two phases are complex conjugates of each other.
!> \param rtbse_env RT-BSE environment
!> \param rho MO-active cfm rotated in place
!> \param i_spin Spin index
!> \param phase Real phase argument (radians); typically +/- Omega_0 * t
! **************************************************************************************************
   SUBROUTINE rotate_rho_phase(rtbse_env, rho, i_spin, phase)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho
      INTEGER, INTENT(IN)                                :: i_spin
      REAL(kind=dp), INTENT(IN)                          :: phase

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rotate_rho_phase'

      COMPLEX(kind=dp)                                   :: phase_ov, phase_vo
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: handle, i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      IF (rtbse_env%omega_shift == 0.0_dp) RETURN

      CALL timeset(routineN, handle)

      n_occ = rtbse_env%n_occ(i_spin)
      shift = rtbse_env%first_active_mo - 1
      phase_ov = CMPLX(COS(phase), SIN(phase), kind=dp)
      phase_vo = CONJG(phase_ov)

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
            IF (i_global_mo <= n_occ .AND. j_global_mo > n_occ) THEN
               ! OV block
               local_data(i_local, j_local) = phase_ov*local_data(i_local, j_local)
            ELSE IF (i_global_mo > n_occ .AND. j_global_mo <= n_occ) THEN
               ! VO block
               local_data(i_local, j_local) = phase_vo*local_data(i_local, j_local)
            END IF
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE rotate_rho_phase

! **************************************************************************************************
!> \brief Build a lab-frame copy of the (possibly rotating-frame) density rho for I/O.
!>        On the TDA + symmetric-shift path returns rho_lab(t) by multiplying OV/VO by
!>        exp(-/+ i Omega_0 t). Otherwise returns a plain copy. Writes into rho_new_last
!>        (mo_struct-sized, idle outside ETRS) and returns a pointer to it; falls back to
!>        the input rho when no scratch is available.
!> \param rtbse_env RT-BSE environment
!> \param rho_in Rotating-frame density (per spin)
!> \param t_phys Physical time associated with rho_in
!> \param rho_lab On exit, points to a per-spin cfm array holding rho in the lab frame.
! **************************************************************************************************
   SUBROUTINE build_rho_lab(rtbse_env, rho_in, t_phys, rho_lab)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_in
      REAL(kind=dp), INTENT(IN)                          :: t_phys
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_lab

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_rho_lab'

      INTEGER                                            :: handle, i

      IF (rtbse_env%omega_shift == 0.0_dp .OR. .NOT. ASSOCIATED(rtbse_env%rho_new_last)) THEN
         rho_lab => rho_in
         RETURN
      END IF

      CALL timeset(routineN, handle)

      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rho_in(i), rtbse_env%rho_new_last(i))
         ! Lab-frame rho_OV(t) = exp(+i*Omega_0*t) * rho_tilde_OV(t)
         ! (sign derived from [h_shifted, rho_tilde]_OV = (-eps_ai + Omega_0) rho_tilde_OV).
         CALL rotate_rho_phase(rtbse_env, rtbse_env%rho_new_last(i), i, rtbse_env%omega_shift*t_phys)
      END DO
      rho_lab => rtbse_env%rho_new_last

      CALL timestop(handle)
   END SUBROUTINE build_rho_lab

! **************************************************************************************************
!> \brief Complex-linear Hartree contraction. Calls the real-input get_hartree on
!>        Re(rho_AO) and on Im(rho_AO) separately and assembles
!>           v_AO = V_H[Re(rho_AO)] + i * V_H[Im(rho_AO)] .
!>        Required by the TDA propagator where the per-pass input Delta rho_OV (or
!>        Delta rho_VO) is non-Hermitian, so the imaginary part must be carried.
!> \param rtbse_env RT-BSE environment
!> \param rho_cfm AO complex input density
!> \param v_cfm AO complex Hartree output (overwritten)
!> \param ispin Spin index (selects scratch slots in rtbse_env)
! **************************************************************************************************
   SUBROUTINE get_hartree_complex(rtbse_env, rho_cfm, v_cfm, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_cfm
      TYPE(cp_cfm_type)                                  :: v_cfm
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'get_hartree_complex'

      INTEGER                                            :: handle

      MARK_USED(ispin)

      CALL timeset(routineN, handle)

      ! Mirrors get_sigma_complex: split rho_cfm into real and imaginary fm parts,
      ! call the real-input Hartree contraction on each, then assemble
      ! v_cfm = V_H[Re(rho_cfm)] + i*V_H[Im(rho_cfm)].
      ! Scratch usage:
      !   real_workspace(1) - holds Re(rho) then V_H[Re]
      !   real_workspace(2) - holds Im(rho) then V_H[Im]
      !   sigma_complex_workspace(1) - cfm wrapper feeding the real-input slot of get_hartree

      ! Real part
      CALL cp_cfm_to_fm(msource=rho_cfm, mtargetr=rtbse_env%real_workspace(1))
      CALL cp_cfm_set_all(rtbse_env%sigma_complex_workspace(1), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), &
                        mtarget=rtbse_env%sigma_complex_workspace(1))
      CALL get_hartree(rtbse_env, rtbse_env%sigma_complex_workspace(1), &
                       rtbse_env%real_workspace(1))

      ! Imaginary part: extract Im(rho_cfm) into real_workspace(2)
      CALL cp_cfm_to_fm(msource=rho_cfm, mtargeti=rtbse_env%real_workspace(2))
      CALL cp_cfm_set_all(rtbse_env%sigma_complex_workspace(1), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(2), &
                        mtarget=rtbse_env%sigma_complex_workspace(1))
      CALL get_hartree(rtbse_env, rtbse_env%sigma_complex_workspace(1), &
                       rtbse_env%real_workspace(2))

      ! Assemble v_cfm = real_workspace(1) + i * real_workspace(2)
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), &
                        msourcei=rtbse_env%real_workspace(2), &
                        mtarget=v_cfm)

      CALL timestop(handle)
   END SUBROUTINE get_hartree_complex

! **************************************************************************************************
!> \brief Custom reimplementation of the delta pulse routines
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml - trafo to MO and linearized version following 10.1021/acs.jctc.2c00644 (03.26)
! **************************************************************************************************
   SUBROUTINE apply_delta_pulse_MO(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_MO'

      INTEGER                                            :: handle, i, k
      REAL(kind=dp)                                      :: intensity, metric
      REAL(kind=dp), DIMENSION(3)                        :: kvec

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
      TYPE(rtbse_env_type), INTENT(IN)                   :: rtbse_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho
      INTEGER, INTENT(IN)                                :: i_spin
      REAL(kind=dp), INTENT(OUT)                         :: full_max, oo_max, vv_max, ov_max, vo_max

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(kind=dp)                                      :: aij
      TYPE(mp_para_env_type), POINTER                    :: para_env

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
!> \brief Zero the OO and VV blocks of an MO-basis derivative cfm in place. Used to enforce
!>        strict linear response on the RK4 derivatives so that the linearized propagator
!>        only carries the OV/VO branches and the OO/VV orbital-energy spreads do not enter
!>        the RK4 stability bound.
!> \param rtbse_env Entry point - rtbse environment
!> \param cfm Derivative-like cfm in MO basis for one spin channel (modified in place)
!> \param i_spin Spin index
! **************************************************************************************************
   SUBROUTINE project_drho_to_ov(rtbse_env, cfm, i_spin)
      TYPE(rtbse_env_type), INTENT(IN)                   :: rtbse_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm
      INTEGER, INTENT(IN)                                :: i_spin

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      n_occ = rtbse_env%n_occ(i_spin)
      shift = rtbse_env%first_active_mo - 1

      CALL cp_cfm_get_info(matrix=cfm, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      local_data => cfm%local_data

      DO j_local = 1, ncol_local
         j_global = col_indices(j_local)
         j_global_mo = j_global + shift
         DO i_local = 1, nrow_local
            i_global = row_indices(i_local)
            i_global_mo = i_global + shift
            IF ((i_global_mo <= n_occ .AND. j_global_mo <= n_occ) .OR. &
                (i_global_mo > n_occ .AND. j_global_mo > n_occ)) THEN
               local_data(i_local, j_local) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
            END IF
         END DO
      END DO
   END SUBROUTINE project_drho_to_ov

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

      CHARACTER(len=*), PARAMETER :: routineN = 'get_electron_number_MO'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: i_global, i_local, j, j_global, j_local, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

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
!> \brief Add +/- delta_omega/2 to the diagonal of ham_reference_singleparticle:
!>        +delta/2 on active occupied diagonals, -delta/2 on active virtual diagonals.
!>        Used to move between lab frame (omega_shift = 0) and rotating frame
!>        (omega_shift = Omega_0). Per spin. No-op when delta_omega = 0.
!> \param rtbse_env RT-BSE environment
!> \param delta_omega Frequency-shift increment to add (subtract: pass negative)
! **************************************************************************************************
   SUBROUTINE shift_href_diagonal(rtbse_env, delta_omega)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      REAL(kind=dp), INTENT(IN)                          :: delta_omega

      CHARACTER(len=*), PARAMETER :: routineN = 'shift_href_diagonal'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: ld
      INTEGER                                            :: abs_mo_idx, handle, i, i_row_global, ii, &
                                                            j_col_global, jj, n_occ_i, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(kind=dp)                                      :: half

      IF (delta_omega == 0.0_dp) RETURN

      CALL timeset(routineN, handle)
      half = 0.5_dp*delta_omega

      DO i = 1, rtbse_env%n_spin
         n_occ_i = rtbse_env%n_occ(i)
         CALL cp_cfm_get_info(matrix=rtbse_env%ham_reference_singleparticle(i), &
                              nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices)
         ld => rtbse_env%ham_reference_singleparticle(i)%local_data
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
               IF (i_row_global /= j_col_global) CYCLE
               abs_mo_idx = i_row_global + rtbse_env%first_active_mo - 1
               IF (abs_mo_idx <= n_occ_i) THEN
                  ld(ii, jj) = ld(ii, jj) + CMPLX(half, 0.0_dp, kind=dp)
               ELSE
                  ld(ii, jj) = ld(ii, jj) - CMPLX(half, 0.0_dp, kind=dp)
               END IF
            END DO
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE shift_href_diagonal

! **************************************************************************************************
!> \brief Fill seed with a deterministic-seeded random complex matrix, Hermitize it,
!>        and return its Frobenius norm. Used to start Krylov power iterations.
!> \param seed cfm to fill (must be allocated on the MO-active struct)
!> \param work_h Hermitization scratch cfm (same struct as seed)
!> \param rng_seed_val Integer seed for RANDOM_SEED (deterministic across runs)
!> \param nrow_local Local row count for seed
!> \param ncol_local Local column count for seed
!> \param norm_seed Output: Frobenius norm of the Hermitized seed
! **************************************************************************************************
   SUBROUTINE init_random_hermitian_seed(seed, work_h, rng_seed_val, &
                                         nrow_local, ncol_local, norm_seed)
      TYPE(cp_cfm_type)                                  :: seed, work_h
      INTEGER, INTENT(IN)                                :: rng_seed_val, nrow_local, ncol_local
      REAL(kind=dp), INTENT(OUT)                         :: norm_seed

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: seed_data
      INTEGER                                            :: i_local, j_local, seed_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rng_seed
      REAL(kind=dp)                                      :: ri, rr

      CALL RANDOM_SEED(SIZE=seed_size)
      ALLOCATE (rng_seed(seed_size))
      rng_seed = rng_seed_val
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
      CALL cp_cfm_transpose(seed, 'C', work_h)
      CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, kind=dp), seed, &
                                CMPLX(0.5_dp, 0.0_dp, kind=dp), work_h)
      norm_seed = cp_cfm_norm(seed, 'F')
   END SUBROUTINE init_random_hermitian_seed

! **************************************************************************************************
!> \brief Two-sweep Krylov estimator of the TDA excitation-spectrum extremes in lab frame.
!>
!>        The linearized RT-BSE Liouvillian L (returned by update_effective_ham_MO,
!>        with the affine constant k0 = L*0 subtracted) has eigenvalues +/- i*Omega_n
!>        for the BSE-A excitation energies Omega_n. Therefore L^2 has real negative
!>        eigenvalues -Omega_n^2.
!>
!>        Sweep 1: power iteration on -L^2 (matvec = -L*(L*v)) converges to the
!>          largest eigenvalue Omega_max^2.
!>        Sweep 2: spectral folding on M = Omega_max^2 * I + L^2. Eigenvalues are
!>          Omega_max^2 - Omega_n^2 (>=0); the largest gives Omega_max^2 - Omega_min^2.
!>          Hence Omega_min = sqrt(Omega_max^2 - dominant fold eigenvalue).
!>
!>        Caller MUST have placed ham_reference_singleparticle in the desired (typically
!>        lab) frame; this routine does not shift. Each iteration costs two
!>        update_effective_ham_MO matvecs. rho_ao_lrrtbse is restored on exit.
!> \param rtbse_env RT-BSE environment
!> \param n_iter Maximum number of power iterations per sweep
!> \param rel_tol Relative tolerance on the iterated eigenvalue for early break
!> \param omega_min Output: smallest TDA excitation energy
!> \param omega_max Output: largest TDA excitation energy
! **************************************************************************************************
   SUBROUTINE compute_spectral_extremes_krylov(rtbse_env, n_iter, rel_tol, omega_min, omega_max)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER, INTENT(IN)                                :: n_iter
      REAL(kind=dp), INTENT(IN)                          :: rel_tol
      REAL(kind=dp), INTENT(OUT)                         :: omega_min, omega_max

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_spectral_extremes_krylov'

      INTEGER                                            :: handle, ispin, it, ncol_local, nrow_local
      LOGICAL                                            :: seed_ok
      REAL(kind=dp)                                      :: lam, lam_prev, norm_seed, omega_max_sq
      TYPE(cp_cfm_type)                                  :: k0, l_v, seed, v, work_h, zero_rho

      CALL timeset(routineN, handle)
      ispin = 1
      omega_min = 0.0_dp
      omega_max = 0.0_dp
      omega_max_sq = 0.0_dp

      CALL cp_cfm_create(seed, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(v, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(l_v, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(k0, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(zero_rho, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(work_h, rtbse_env%fm_struct_mo_active)

      ! Affine offset k0 = f(rho = 0); L*v = f(v) - k0.
      CALL cp_cfm_set_all(zero_rho, CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL update_effective_ham_MO(rtbse_env, zero_rho, k0, ispin)

      CALL cp_cfm_get_info(matrix=seed, &
                           nrow_local=nrow_local, ncol_local=ncol_local)

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' RTBSE| ===== Krylov spectral extremes (TDA, lab frame) ====='
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' RTBSE| Sweep 1: power iteration on -L^2  ->  Omega_max'
         WRITE (rtbse_env%unit_nr, '(A,T13,A,T35,A,T57,A)') &
            ' RTBSE|', 'iter', 'Omega_max [eV]', 'rel.change'
      END IF

      ! ----- Sweep 1: power iteration on -L^2 -> Omega_max^2 -----
      CALL init_random_hermitian_seed(seed, work_h, 20260507, nrow_local, ncol_local, &
                                      norm_seed)
      CALL project_drho_to_ov(rtbse_env, seed, ispin)
      norm_seed = cp_cfm_norm(seed, 'F')
      seed_ok = (norm_seed > 0.0_dp)
      IF (seed_ok) THEN
         CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), seed)
         lam = 0.0_dp
         lam_prev = 0.0_dp
         DO it = 1, n_iter
            ! v = L * seed
            CALL update_effective_ham_MO(rtbse_env, seed, v, ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), k0)
            CALL project_drho_to_ov(rtbse_env, v, ispin)
            ! l_v = L * v = L^2 * seed
            CALL update_effective_ham_MO(rtbse_env, v, l_v, ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), l_v, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), k0)
            CALL project_drho_to_ov(rtbse_env, l_v, ispin)
            ! Operate with -L^2 to expose positive eigenvalues +Omega_n^2
            CALL cp_cfm_scale(CMPLX(-1.0_dp, 0.0_dp, kind=dp), l_v)

            lam = cp_cfm_norm(l_v, 'F')
            IF (rtbse_env%unit_nr > 0) THEN
               IF (it == 1) THEN
                  WRITE (rtbse_env%unit_nr, '(A,T13,I4,T35,F10.4,T57,A)') &
                     ' RTBSE|', it, SQRT(MAX(lam, 0.0_dp))*evolt, '--'
               ELSE
                  WRITE (rtbse_env%unit_nr, '(A,T13,I4,T35,F10.4,T55,ES12.3)') &
                     ' RTBSE|', it, SQRT(MAX(lam, 0.0_dp))*evolt, &
                     ABS(lam - lam_prev)/MAX(lam, TINY(1.0_dp))
               END IF
            END IF
            IF (it > 1 .AND. ABS(lam - lam_prev) <= &
                rel_tol*MAX(lam, TINY(1.0_dp))) EXIT
            lam_prev = lam
            ! Re-Hermitize and renormalize
            CALL cp_cfm_to_cfm(l_v, seed)
            CALL cp_cfm_transpose(seed, 'C', work_h)
            CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, kind=dp), seed, &
                                      CMPLX(0.5_dp, 0.0_dp, kind=dp), work_h)
            norm_seed = cp_cfm_norm(seed, 'F')
            IF (norm_seed <= 0.0_dp) EXIT
            CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), seed)
         END DO
         omega_max_sq = MAX(lam, 0.0_dp)
         omega_max = SQRT(omega_max_sq)
      END IF

      ! Sweep 2 (folded power iteration on Omega_max^2 I + L^2 -> Omega_min) intentionally
      ! skipped: it relies on catastrophic cancellation Omega_max^2 - lam, which gave a
      ! noisy/under-estimated Omega_min and pushed RK4 past its stability disk. The
      ! conservative choice Omega_min = 0 gives the symmetric worst-case shift around zero.
      omega_min = 0.0_dp

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A,T63,F18.4)') &
            ' RTBSE| Omega_min [eV] (forced to 0)', omega_min*evolt
         WRITE (rtbse_env%unit_nr, '(A,T63,F18.4)') &
            ' RTBSE| Omega_max [eV]', omega_max*evolt
      END IF

      ! Restore rho_ao_lrrtbse from rho_orig (update_effective_ham_MO overwrote it).
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_orig(ispin), &
                                                rtbse_env%rho_ao_lrrtbse(ispin), ispin)

      CALL cp_cfm_release(seed)
      CALL cp_cfm_release(v)
      CALL cp_cfm_release(l_v)
      CALL cp_cfm_release(k0)
      CALL cp_cfm_release(zero_rho)
      CALL cp_cfm_release(work_h)

      CALL timestop(handle)
   END SUBROUTINE compute_spectral_extremes_krylov

! **************************************************************************************************
!> \brief Refine the TDA symmetric frequency shift using the two-step Krylov
!>        spectral-extremes estimator. Undoes the bare-gap-based shift on
!>        ham_reference_singleparticle, measures Omega_min/Omega_max of the full
!>        (kernel-inclusive) TDA spectrum, sets Omega_0 = (Omega_min+Omega_max)/2
!>        and omega_max = (Omega_max-Omega_min)/2, reapplies the optimal shift,
!>        and updates maximum_timestep. If ENFORCE_MAX_DT is on, rewrites
!>        sim_dt and sim_nsteps accordingly.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE refine_omega_shift_via_krylov(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'refine_omega_shift_via_krylov'

      CHARACTER(len=256)                                 :: hint_msg
      INTEGER                                            :: handle, n_steps_new, n_steps_old
      REAL(kind=dp)                                      :: omega_max, omega_max_new, omega_min, &
                                                            omega_shift_new, omega_shift_old, &
                                                            sim_dt_as, total_time

      CALL timeset(routineN, handle)

      omega_shift_old = rtbse_env%omega_shift

      ! 1) Undo the bare-gap shift -> lab frame
      CALL shift_href_diagonal(rtbse_env, -omega_shift_old)

      ! 2) Measure Omega_min, Omega_max with the kernel included
      CALL compute_spectral_extremes_krylov(rtbse_env, n_iter=30, rel_tol=1.0E-5_dp, &
                                            omega_min=omega_min, omega_max=omega_max)

      ! 3) Optimal symmetric shift
      omega_shift_new = 0.5_dp*(omega_min + omega_max)
      omega_max_new = 0.5_dp*(omega_max - omega_min)

      ! 4) Apply the new shift to ham_reference_singleparticle
      CALL shift_href_diagonal(rtbse_env, omega_shift_new)

      ! 5) Update env
      rtbse_env%omega_shift = omega_shift_new
      IF (omega_max_new > 0.0_dp) THEN
         rtbse_env%omega_max = omega_max_new
         rtbse_env%maximum_timestep = 2.0_dp*SQRT(2.0_dp)/omega_max_new
      END IF

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A,F16.6,A,F16.6,A)') &
            ' Krylov-refined Omega_0   = ', omega_shift_new*evolt, &
            ' eV   (gap-based was ', omega_shift_old*evolt, ' eV)'
         WRITE (rtbse_env%unit_nr, '(A,F16.6,A,F16.4,A)') &
            ' Krylov-refined omega_max = ', omega_max_new*evolt, &
            ' eV   ->   dt_max = ', rtbse_env%maximum_timestep*seconds*1e18_dp, ' as'
      END IF

      ! 6) If requested, rewrite sim_dt/sim_nsteps to fit the (now tighter) bound
      IF (rtbse_env%enforce_max_dt .AND. rtbse_env%sim_dt > 0.0_dp .AND. omega_max_new > 0.0_dp) THEN
         n_steps_old = MAX(0, rtbse_env%sim_nsteps)
         total_time = REAL(n_steps_old, dp)*rtbse_env%sim_dt
         n_steps_new = MAX(1, CEILING(total_time/rtbse_env%maximum_timestep))
         rtbse_env%sim_dt = total_time/REAL(n_steps_new, dp)
         IF (rtbse_env%sim_nsteps /= n_steps_new) THEN
            rtbse_env%sim_nsteps = n_steps_new
            CALL reallocate_ft_traces(rtbse_env)
         END IF
         sim_dt_as = rtbse_env%sim_dt*seconds*1e18_dp
         WRITE (hint_msg, '(A,F16.4,A,I0,A)') &
            'Krylov-refined ENFORCE_MAX_DT: TIMESTEP -> ', sim_dt_as, &
            ' as,  STEPS -> ', n_steps_new, '.'
         CALL cp_hint(__LOCATION__, TRIM(hint_msg))
      END IF

      CALL timestop(handle)
   END SUBROUTINE refine_omega_shift_via_krylov

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
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER, INTENT(IN)                                :: n_iter
      REAL(kind=dp), INTENT(IN)                          :: rel_tol

      CHARACTER(len=*), PARAMETER :: routineN = 'estimate_omega_max_krylov'

      INTEGER                                            :: handle, ispin, it, ncol_local, nrow_local
      LOGICAL                                            :: seed_ok
      REAL(kind=dp)                                      :: dt_max, dt_max_gap, norm_seed, omega, &
                                                            omega_gap, omega_prev
      TYPE(cp_cfm_type)                                  :: k0, seed, v, work_h, zero_rho

      CALL timeset(routineN, handle)

      ispin = 1
      omega_gap = rtbse_env%omega_max
      dt_max_gap = 0.0_dp
      IF (omega_gap > 0.0_dp) dt_max_gap = 2.0_dp*SQRT(2.0_dp)/omega_gap

      ! Allocate scratch on MO-active struct
      CALL cp_cfm_create(seed, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(v, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(k0, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(zero_rho, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(work_h, rtbse_env%fm_struct_mo_active)

      ! 1) Constant offset k0 = update_effective_ham_MO(rho = 0). update_effective_ham_MO
      !    is affine in rho; subtracting k0 yields the pure linear action L*seed.
      CALL cp_cfm_set_all(zero_rho, CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL update_effective_ham_MO(rtbse_env, zero_rho, k0, ispin)

      ! 2) Random Hermitian seed projected onto the OV/VO subspace so the iteration
      !    measures rho(P L P) - the spectral radius of the operator that the
      !    OV-projected RK4 actually integrates.
      CALL cp_cfm_get_info(matrix=seed, &
                           nrow_local=nrow_local, ncol_local=ncol_local)
      CALL init_random_hermitian_seed(seed, work_h, 20260506, nrow_local, ncol_local, &
                                      norm_seed)
      CALL project_drho_to_ov(rtbse_env, seed, ispin)
      norm_seed = cp_cfm_norm(seed, 'F')
      seed_ok = (norm_seed > 0.0_dp)
      IF (.NOT. seed_ok .AND. rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' RTBSE| Krylov rho(L): seed norm is zero on this distribution; skipping.'
      END IF

      IF (seed_ok) THEN
         CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), seed)

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A)') &
               ' RTBSE| ===== Krylov rho(L) debug estimate ====='
            WRITE (rtbse_env%unit_nr, '(A,T35,F10.4,T55,F12.4)') &
               ' RTBSE| Gap omega_max [eV] / dt [as]', &
               omega_gap*evolt, dt_max_gap*seconds*1e18_dp
            WRITE (rtbse_env%unit_nr, '(A,T13,A,T28,A,T45,A,T62,A)') &
               ' RTBSE|', 'iter', 'rho(L) [eV]', 'dt_max [as]', 'rel.change'
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
            IF (rtbse_env%tda_active) &
               CALL project_drho_to_ov(rtbse_env, v, ispin)

            omega = cp_cfm_norm(v, 'F')
            IF (omega > 0.0_dp) dt_max = 2.0_dp*SQRT(2.0_dp)/omega

            IF (rtbse_env%unit_nr > 0) THEN
               IF (it == 1) THEN
                  WRITE (rtbse_env%unit_nr, '(A,T13,I4,T28,F10.4,T45,F12.4,T62,A)') &
                     ' RTBSE|', it, omega*evolt, dt_max*seconds*1e18_dp, '--'
               ELSE
                  WRITE (rtbse_env%unit_nr, &
                         '(A,T13,I4,T28,F10.4,T45,F12.4,T60,ES12.3)') &
                     ' RTBSE|', it, omega*evolt, dt_max*seconds*1e18_dp, &
                     ABS(omega - omega_prev)/MAX(omega, TINY(1.0_dp))
               END IF
            END IF

            IF (it > 1 .AND. ABS(omega - omega_prev) <= &
                rel_tol*MAX(omega, TINY(1.0_dp))) EXIT
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
            WRITE (rtbse_env%unit_nr, '(A,T35,F10.4,T55,F12.4)') &
               ' RTBSE| Krylov rho(L) [eV] / dt [as]', &
               omega*evolt, dt_max*seconds*1e18_dp
            IF (omega_gap > 0.0_dp) THEN
               WRITE (rtbse_env%unit_nr, '(A,T61,F18.4)') &
                  ' RTBSE| rho(L) / omega_gap', omega/omega_gap
            END IF
         END IF
      END IF

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

! **************************************************************************************************
!> \brief Real Frobenius inner product <A,B>_R = Re Tr(A^H B), distributed.
!> \param a left operand
!> \param b right operand
!> \param res returned real inner product (summed across para_env)
! **************************************************************************************************
   SUBROUTINE cfm_real_inner_product(a, b, res)
      TYPE(cp_cfm_type), INTENT(IN)                      :: a, b
      REAL(kind=dp), INTENT(OUT)                         :: res

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: ad, bd
      INTEGER                                            :: i, j, ncol_local, nrow_local
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL cp_cfm_get_info(matrix=a, nrow_local=nrow_local, ncol_local=ncol_local)
      ad => a%local_data
      bd => b%local_data
      res = 0.0_dp
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            res = res + REAL(CONJG(ad(i, j))*bd(i, j), dp)
         END DO
      END DO
      para_env => a%matrix_struct%para_env
      CALL para_env%sum(res)
   END SUBROUTINE cfm_real_inner_product

! **************************************************************************************************
!> \brief Hermitian Lanczos diagnostic for the kernel-inclusive linearized propagator.
!>        TDA-only: in TDA the kernel response is self-adjoint on the Hermitian OV+VO
!>        drho subspace, so plain Lanczos applies. Ritz values converge superlinearly
!>        to lambda_max(L); a faster alternative to the power-iteration variant on
!>        clustered top spectra. Diagnostic only - does NOT overwrite omega_max.
!> \param rtbse_env environment
!> \param n_iter maximum Lanczos iterations
!> \param rel_tol relative tolerance on the largest Ritz value for early exit
! **************************************************************************************************
   SUBROUTINE estimate_omega_max_lanczos(rtbse_env, n_iter, rel_tol)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER, INTENT(IN)                                :: n_iter
      REAL(kind=dp), INTENT(IN)                          :: rel_tol

      CHARACTER(len=*), PARAMETER :: routineN = 'estimate_omega_max_lanczos'

      INTEGER                                            :: handle, info, ispin, it, j, ncol_local, &
                                                            nrow_local
      LOGICAL                                            :: seed_ok
      REAL(kind=dp)                                      :: alpha_k, beta_k, c, dt_max, dt_max_gap, &
                                                            norm_seed, omega, omega_gap, omega_prev
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: alpha, beta, d_loc, e_loc, work
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: z_dummy
      TYPE(cp_cfm_type)                                  :: k0, v, work_h, zero_rho
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: q

      CALL timeset(routineN, handle)

      ispin = 1
      omega_gap = rtbse_env%omega_max
      dt_max_gap = 0.0_dp
      IF (omega_gap > 0.0_dp) dt_max_gap = 2.0_dp*SQRT(2.0_dp)/omega_gap

      ! Lanczos basis q(1..n_iter+1) + scratch
      ALLOCATE (q(n_iter + 1))
      DO it = 1, n_iter + 1
         CALL cp_cfm_create(q(it), rtbse_env%fm_struct_mo_active)
      END DO
      CALL cp_cfm_create(v, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(k0, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(zero_rho, rtbse_env%fm_struct_mo_active)
      CALL cp_cfm_create(work_h, rtbse_env%fm_struct_mo_active)

      ALLOCATE (alpha(n_iter), beta(n_iter))
      alpha = 0.0_dp
      beta = 0.0_dp

      ! Constant offset: linearise the affine action of update_effective_ham_MO.
      CALL cp_cfm_set_all(zero_rho, CMPLX(0.0_dp, 0.0_dp, kind=dp))
      CALL update_effective_ham_MO(rtbse_env, zero_rho, k0, ispin)

      ! Random Hermitian seed, OV-projected so the iteration probes rho(P L P).
      CALL cp_cfm_get_info(matrix=q(1), nrow_local=nrow_local, ncol_local=ncol_local)
      CALL init_random_hermitian_seed(q(1), work_h, 20260512, nrow_local, ncol_local, &
                                      norm_seed)
      CALL project_drho_to_ov(rtbse_env, q(1), ispin)
      norm_seed = cp_cfm_norm(q(1), 'F')
      seed_ok = (norm_seed > 0.0_dp)
      IF (.NOT. seed_ok .AND. rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') &
            ' RTBSE| Lanczos rho(L): seed norm is zero on this distribution; skipping.'
      END IF

      omega = 0.0_dp
      omega_prev = 0.0_dp
      dt_max = 0.0_dp

      IF (seed_ok) THEN
         CALL cp_cfm_scale(CMPLX(1.0_dp/norm_seed, 0.0_dp, kind=dp), q(1))

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A)') &
               ' RTBSE| ===== Lanczos rho(L) debug estimate ====='
            WRITE (rtbse_env%unit_nr, '(A,T35,F10.4,T55,F12.4)') &
               ' RTBSE| Gap omega_max [eV] / dt [as]', &
               omega_gap*evolt, dt_max_gap*seconds*1e18_dp
            WRITE (rtbse_env%unit_nr, '(A,T13,A,T28,A,T45,A,T62,A)') &
               ' RTBSE|', 'iter', 'rho(L) [eV]', 'dt_max [as]', 'rel.change'
         END IF

         DO it = 1, n_iter
            ! v = L * q_it = update_effective_ham_MO(q_it) - k0, OV-projected.
            CALL update_effective_ham_MO(rtbse_env, q(it), v, ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), k0)
            CALL project_drho_to_ov(rtbse_env, v, ispin)

            ! alpha_k = Re <q_it, v>
            CALL cfm_real_inner_product(q(it), v, alpha_k)
            alpha(it) = alpha_k

            ! v <- v - alpha_k q_it - beta_{k-1} q_{k-1}
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                      CMPLX(-alpha_k, 0.0_dp, kind=dp), q(it))
            IF (it > 1) THEN
               CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                         CMPLX(-beta(it - 1), 0.0_dp, kind=dp), q(it - 1))
            END IF

            ! Full reorthogonalisation against q(1..it) for numerical safety.
            DO j = 1, it
               CALL cfm_real_inner_product(q(j), v, c)
               CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), v, &
                                         CMPLX(-c, 0.0_dp, kind=dp), q(j))
            END DO

            beta_k = cp_cfm_norm(v, 'F')
            beta(it) = beta_k

            ! Largest |Ritz| of T_it = tridiag(alpha(1:it), beta(1:it-1)).
            IF (it == 1) THEN
               omega = ABS(alpha(1))
            ELSE
               ALLOCATE (d_loc(it), e_loc(it - 1), work(MAX(1, 2*it - 2)), z_dummy(1, 1))
               d_loc(1:it) = alpha(1:it)
               e_loc(1:it - 1) = beta(1:it - 1)
               info = 0
               CALL dstev('N', it, d_loc, e_loc, z_dummy, 1, work, info)
               IF (info == 0) THEN
                  omega = MAX(ABS(d_loc(1)), ABS(d_loc(it)))
               ELSE
                  omega = 0.0_dp
               END IF
               DEALLOCATE (d_loc, e_loc, work, z_dummy)
            END IF

            IF (omega > 0.0_dp) dt_max = 2.0_dp*SQRT(2.0_dp)/omega

            IF (rtbse_env%unit_nr > 0) THEN
               IF (it == 1) THEN
                  WRITE (rtbse_env%unit_nr, '(A,T13,I4,T28,F10.4,T45,F12.4,T62,A)') &
                     ' RTBSE|', it, omega*evolt, dt_max*seconds*1e18_dp, '--'
               ELSE
                  WRITE (rtbse_env%unit_nr, &
                         '(A,T13,I4,T28,F10.4,T45,F12.4,T60,ES12.3)') &
                     ' RTBSE|', it, omega*evolt, dt_max*seconds*1e18_dp, &
                     ABS(omega - omega_prev)/MAX(omega, TINY(1.0_dp))
               END IF
            END IF

            IF (it > 1 .AND. ABS(omega - omega_prev) <= &
                rel_tol*MAX(omega, TINY(1.0_dp))) EXIT
            omega_prev = omega

            IF (beta_k <= 0.0_dp) EXIT
            IF (it < n_iter) THEN
               CALL cp_cfm_to_cfm(v, q(it + 1))
               CALL cp_cfm_scale(CMPLX(1.0_dp/beta_k, 0.0_dp, kind=dp), q(it + 1))
            END IF
         END DO

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,T35,F10.4,T55,F12.4)') &
               ' RTBSE| Lanczos rho(L) [eV] / dt [as]', &
               omega*evolt, dt_max*seconds*1e18_dp
            IF (omega_gap > 0.0_dp) THEN
               WRITE (rtbse_env%unit_nr, '(A,T61,F18.4)') &
                  ' RTBSE| rho(L) / omega_gap', omega/omega_gap
            END IF
         END IF
      END IF

      ! Restore rho_ao_lrrtbse from rho_orig (update_effective_ham_MO overwrites it).
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_orig(ispin), &
                                                rtbse_env%rho_ao_lrrtbse(ispin), ispin)

      DO it = 1, n_iter + 1
         CALL cp_cfm_release(q(it))
      END DO
      DEALLOCATE (q)
      CALL cp_cfm_release(v)
      CALL cp_cfm_release(k0)
      CALL cp_cfm_release(zero_rho)
      CALL cp_cfm_release(work_h)
      DEALLOCATE (alpha, beta)

      CALL timestop(handle)
   END SUBROUTINE estimate_omega_max_lanczos

END MODULE rt_bse_linearized

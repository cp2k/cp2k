!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the propagation of the linearized RT-BSE equations of motion.
!>        Propagates the first-order density matrix response Δρ within the active MO window,
!>        in the Tamm-Dancoff approximation or with the full (A, B) coupling, instead of the
!>        lesser Green's function propagated by rt_bse. Also provides the Liouvillian eigenvalue
!>        diagnostic, which builds the Liouvillian by probing the kernel with canonical basis
!>        vectors and diagonalizes it.
!> \note  The control is handed directly from cp2k_runs
!>        The initialization and delta-kick routines are adapted from the full RT-BSE
!>        propagator in rt_bse.F.
!> \author Maximilian Graml (03.26)
!> \author Stepan Marek (09.24) - original RT-BSE routines adapted here
! **************************************************************************************************

MODULE rt_bse_linearized
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_norm,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: &
        cp_cfm_get_info, cp_cfm_get_submatrix, cp_cfm_set_all, cp_cfm_set_element, &
        cp_cfm_set_submatrix, cp_cfm_to_cfm, cp_cfm_to_fm, cp_cfm_type, cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
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
   USE input_constants,                 ONLY: evGW0,&
                                              rtp_bse_ham_gw,&
                                              use_mom_ref_zero,&
                                              use_rt_restart
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt,&
                                              seconds
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power
   USE rt_bse,                          ONLY: get_hartree,&
                                              get_sigma,&
                                              init_hartree,&
                                              initialize_rtbse_env,&
                                              propagate_density,&
                                              rho_metric
   USE rt_bse_io,                       ONLY: &
        check_restart_eps_consistency, output_field, output_moments, output_mos_contravariant, &
        output_restart_linearized, print_timestep_info, read_restart_C, read_restart_density, &
        read_restart_info, read_restart_trace
   USE rt_bse_ri_rs,                    ONLY: compute_hartree_ri_rs,&
                                              compute_hartree_ri_rs_complex,&
                                              compute_hartree_ri_rs_from_diag,&
                                              rt_bse_ri_rs_ensure_V_grid,&
                                              rt_bse_ri_rs_ensure_W0_grid
   USE rt_bse_types,                    ONLY: create_rtbse_env,&
                                              release_rtbse_env,&
                                              rtbse_env_type
   USE rt_propagation_output,           ONLY: print_ft
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_linearized'

   ! build_shared_sex_and_hartree input-convention selector. The mask choice also fixes the input
   ! Hermiticity, which gates the Hartree imaginary channel (Im computed iff non-Hermitian = OV only).
   INTEGER, PARAMETER, PRIVATE :: kernel_input_ov = 1, kernel_input_ovvo = 2, kernel_input_full = 3

   PUBLIC :: run_propagation_linearized_bse

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time propagation of the linearized BSE
!> \param force_env Force environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_linearized_bse(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'run_propagation_linearized_bse'

      INTEGER                                            :: handle, i, j
      REAL(kind=dp)                                      :: t_phys, t_start, timestep_walltime, &
                                                            timestep_walltime_start
      REAL(kind=dp), DIMENSION(2)                        :: enum_im, enum_re
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_lab
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      ! Per-spin (alpha/beta) electron numbers; only 1:n_spin entries are used.

      CALL timeset(routineN, handle)

      CALL cp_warn(__LOCATION__, &
                   "Linearized RT-BSE is under active development. Make sure you understand "// &
                   "the method and validate results before using it for production calculations.")

      ! To Do: Bibliography information

      logger => cp_get_default_logger()

      ! Run the initial SCF calculation / read SCF restart information
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)

      ! Allocate all persistant storage and read input that does not need further processing
      CALL create_rtbse_env(rtbse_env, force_env, linearized=.TRUE.)

      ! Restart phase 1a: read sim_start + the original run's dt (from the trace header) BEFORE
      ! ENFORCE_MAX_DT, so the continuation inherits that dt instead of a window-dependent one.
      IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
         CALL read_restart_info(rtbse_env)
      END IF

      CALL initialize_maximum_timestep(rtbse_env)

      ! Restart phase 1b: load the trace prefix now that ENFORCE_MAX_DT has sized the trace arrays.
      IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
         IF (rtbse_env%sim_start >= rtbse_env%sim_nsteps) THEN
            CPABORT("RT_RESTART: restart step >= STEPS - increase MOTION%MD%STEPS")
         END IF
         CALL read_restart_trace(rtbse_env)
      END IF

      CALL print_linrtbse_header_info(rtbse_env)

      ! Build the truncated MO coefficient slabs C_active(:, first_active_mo..last_active_mo)
      ! used by all AO<->MO transforms in the linearized path.
      CALL populate_C_active(rtbse_env)

      ! Initiate iteration level "MD" in order to copy the structure of other RTP codes
      CALL cp_add_iter_level(logger%iter_info, "MD")
      ! Initialize non-trivial values
      !  - calculates the moment operators
      CALL initialize_moments(rtbse_env)
      !  - populates overlap and inverse overlap matrices
      CALL initialize_rtbse_env(rtbse_env)

      !  - populates the fresh SCF density matrix rho^0 (and the rho_orig reference for delta rho)
      CALL initialize_density_matrix(rtbse_env)

      ! Restart phase 2: overwrite rho from the lab-frame restart files, bridge into this run's MO
      ! gauge, then enter this run's rotating frame (rotate_rho_phase is a no-op when omega_shift=0)
      IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
         CALL read_restart_density(rtbse_env)
         IF (rtbse_env%restart_extracted) THEN
            CALL apply_restart_basis_bridge(rtbse_env)
            t_start = REAL(rtbse_env%sim_start, dp)*rtbse_env%sim_dt
            DO i = 1, rtbse_env%n_spin
               CALL rotate_rho_phase(rtbse_env, rtbse_env%rho(i), i, -rtbse_env%omega_shift*t_start)
            END DO
         END IF
      END IF
      !  - calculates/populates the G0W0/KS Hamiltonian, respectively
      CALL initialize_singleparticle_hamiltonian(rtbse_env)
      ! Restart Hamiltonian-consistency heads-up: eps_active exists only now, so compare here
      CALL check_restart_eps_consistency(rtbse_env)
      ! Transform initial density matrix to AO basis for use in Hartree and self-energy calculations
      DO i = 1, rtbse_env%n_spin
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_orig(i), rtbse_env%rho_ao_scratch(i), i)
      END DO
      !  - calculates the Hartree reference potential
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_set_all(rtbse_env%ham_reference(i), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      END DO
      CALL initialize_hartree_potential(rtbse_env)
      !  - calculates the SEX reference self-energy
      CALL initialize_sex_selfenergy(rtbse_env)

      ! Liouvillian eigenvalue diagnostic (one-shot at init, TDA or ABBA via dispatcher).
      ! Detached from the propagator; safe to call after the reference-init routines.
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

         CALL solve_rk4_timestep(rtbse_env, rtbse_env%rho, rtbse_env%rho_new)
         CALL get_electron_number_MO(rtbse_env, rtbse_env%rho_new, &
                                     enum_re(1:rtbse_env%n_spin), enum_im(1:rtbse_env%n_spin))
         timestep_walltime = m_walltime() - timestep_walltime_start
         CALL print_timestep_info(rtbse_env, i, enum_re(1:rtbse_env%n_spin), step_walltime=timestep_walltime)
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
         ! Output restart files, so that the restart resumes at the step recorded in .info
         CALL output_restart_linearized(rtbse_env, rho_lab)
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
!> \brief Computes the analytic RK4 stability bound t* = 2√2 / Ω_max (a.u.) from the largest active
!>        GW/KS gap, and (TDA + first-peak) the rotating-frame shift Ω_0 = ε^ai_min.
!>        Writes the timestep diagnostics to stdout; optionally rewrites TIMESTEP/STEPS under
!>        ENFORCE_MAX_DT; on restart it inherits the original dt from the trace and only rescales STEPS.
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

      IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
         omega_max = MAXVAL(rtbse_env%bs_env%eigenval_GW(i_first:i_last, :, :)) - &
                     MINVAL(rtbse_env%bs_env%eigenval_GW(i_first:i_last, :, :))
      ELSE
         omega_max = MAXVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:i_last, :)) - &
                     MINVAL(rtbse_env%bs_env%eigenval_scf_Gamma(i_first:i_last, :))
      END IF

      ! First-peak shift (TDA only): set Ω_0 = eps_min_ai so the lowest
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
               IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_GW(i_first:rtbse_env%n_occ(ispin), :, ispin))
                  eps_occ_min = MIN(eps_occ_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_GW(i_first:rtbse_env%n_occ(ispin), :, ispin))
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
               IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
                  ev_tmp = MINVAL(rtbse_env%bs_env%eigenval_GW(rtbse_env%n_occ(ispin) + 1:i_last, :, ispin))
                  eps_virt_min = MIN(eps_virt_min, ev_tmp)
                  ev_tmp = MAXVAL(rtbse_env%bs_env%eigenval_GW(rtbse_env%n_occ(ispin) + 1:i_last, :, ispin))
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
         ! t* = 2√2 / Ω_max  (a.u.; imaginary-axis RK4 bound |R(iy)| ≤ 1 at y = 2√2)
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
         IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart .AND. &
             rtbse_env%sim_dt_restart > 0.0_dp) THEN
            ! Continuation: dt is frozen in the trace, so inherit it and only rescale the step count
            ! to the requested window. Recomputing dt from the (longer) window would desync the trace
            ! time-grid and trip the continuation guard in read_restart_trace.
            rtbse_env%sim_dt = rtbse_env%sim_dt_restart
            n_steps_new = MAX(1, NINT(total_time/rtbse_env%sim_dt))
            rtbse_env%sim_nsteps = n_steps_new
            sim_dt_as = rtbse_env%sim_dt*seconds*1e18_dp
            WRITE (hint_msg, '(A,F16.4,A,I0,A)') &
               'ENFORCE_MAX_DT on restart: inheriting original TIMESTEP ', sim_dt_as, &
               ' as and setting STEPS to ', n_steps_new, '.'
            CALL cp_hint(__LOCATION__, TRIM(hint_msg))
            ! The inherited dt was stable in the original run; warn only if this run's stability
            ! window shrank below it (e.g. the recomputed GW eigenvalues shifted the gap).
            IF (rtbse_env%sim_dt > rtbse_env%maximum_timestep/grace_factor) THEN
               CALL cp_warn(__LOCATION__, &
                            "ENFORCE_MAX_DT restart: inherited dt exceeds this run's stability "// &
                            "limit - the recomputed Hamiltonian may make the propagation unstable.")
            END IF
         ELSE
            n_steps_new = MAX(1, CEILING(total_time/(rtbse_env%maximum_timestep/grace_factor)))
            rtbse_env%sim_dt = total_time/REAL(n_steps_new, dp)
            rtbse_env%sim_nsteps = n_steps_new
            sim_dt_as = rtbse_env%sim_dt*seconds*1e18_dp
            WRITE (hint_msg, '(A,F16.4,A,I0,A)') &
               'ENFORCE_MAX_DT enabled. Resetting TIMESTEP to ', sim_dt_as, &
               ' as and STEPS to ', n_steps_new, '.'
            CALL cp_hint(__LOCATION__, TRIM(hint_msg))
         END IF
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
!> \brief Prints the linRTBSE run header to stdout: active-MO window (first/last/count) and the
!>        occupied/virtual energy cutoffs.
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_linrtbse_header_info(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env

      INTEGER                                            :: ispin, n_steps
      REAL(kind=dp)                                      :: e_first, e_last, fft_resolution, &
                                                            nyquist_frequency, total_time
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

         ! Which single-particle eigenvalues the propagation uses. bs_env%eigenval_GW carries the
         ! result of the highest GW flavour requested.
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
            IF (rtbse_env%bs_env%gw_flavour == evGW0) THEN
               WRITE (rtbse_env%unit_nr, '(A,T75,A6)') &
                  ' GW flavor for computing GW eigenvalues used in RT-BSE:', ' evGW0'
            ELSE
               WRITE (rtbse_env%unit_nr, '(A,T75,A6)') &
                  ' GW flavor for computing GW eigenvalues used in RT-BSE:', '  G0W0'
            END IF
         ELSE
            WRITE (rtbse_env%unit_nr, '(A,T75,A6)') &
               ' Single-particle eigenvalues used in RT-BSE:', '    KS'
         END IF

         ! Active MO window (energy-cutoff truncation) for linearized RT-BSE
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
         IF (rtbse_env%active_mo_truncation) THEN
            ! The window is cut on the DFT axis but propagated on the QP axis, so the QP edges may
            ! exceed the nominal cutoff. Print both so the window can be checked against the input.
            DO ispin = 1, rtbse_env%n_spin
               e_first = (rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%first_active_mo, ispin) - &
                          rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%n_occ(ispin), ispin))*evolt
               e_last = (rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%last_active_mo, ispin) - &
                         rtbse_env%bs_env%eigenval_scf_Gamma(rtbse_env%n_occ(ispin) + 1, ispin))*evolt
               WRITE (rtbse_env%unit_nr, '(A,I1,A,T71,F10.3)') ' Spin ', ispin, &
                  ' first active MO, E - E_HOMO (KS) [eV]:', e_first
               WRITE (rtbse_env%unit_nr, '(A,I1,A,T71,F10.3)') ' Spin ', ispin, &
                  ' last active MO, E - E_LUMO (KS) [eV]:', e_last
               IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
                  e_first = (rtbse_env%bs_env%eigenval_GW(rtbse_env%first_active_mo, 1, ispin) - &
                             rtbse_env%bs_env%eigenval_GW(rtbse_env%n_occ(ispin), 1, ispin))*evolt
                  e_last = (rtbse_env%bs_env%eigenval_GW(rtbse_env%last_active_mo, 1, ispin) - &
                            rtbse_env%bs_env%eigenval_GW(rtbse_env%n_occ(ispin) + 1, 1, ispin))*evolt
                  WRITE (rtbse_env%unit_nr, '(A,I1,A,T71,F10.3)') ' Spin ', ispin, &
                     ' first active MO, E - E_HOMO (QP) [eV]:', e_first
                  WRITE (rtbse_env%unit_nr, '(A,I1,A,T71,F10.3)') ' Spin ', ispin, &
                     ' last active MO, E - E_LUMO (QP) [eV]:', e_last
               END IF
            END DO
         END IF
      END IF

   END SUBROUTINE print_linrtbse_header_info

! **************************************************************************************************
!> \brief Populates rtbse_env%C_active(i_spin) (n_ao x mo_active) by extracting columns
!>        first_active_mo..last_active_mo from bs_env%fm_mo_coeff_Gamma(i_spin).
!> \param rtbse_env RT-BSE environment
!> \author Maximilian Graml (05.26)
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
!> \brief Builds the dipole moment operators r_mn in the active-MO basis (per axis, per spin) from
!>        the AO moment matrices: moments(k,σ) at the reference point, moments_field(k,σ) at origin.
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03.26) - refactor in prep. of linearized propagation
! **************************************************************************************************
   SUBROUTINE initialize_moments(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_moments'

      INTEGER                                            :: handle, i_spin, k
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
         DO i_spin = 1, rtbse_env%n_spin
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, tmp_ao, rtbse_env%moments(k, i_spin), i_spin)
         END DO
      END DO
      ! TODO: remove moments_field (only needed for the TDDFT comparison)
      ! Now, repeat without reference point to get the moments for field
      CALL get_reference_point(rpoint, qs_env=rtbse_env%qs_env, &
                               reference=use_mom_ref_zero)
      CALL build_local_moment_matrix(rtbse_env%qs_env, moments_dbcsr_p, 1, rpoint)
      DO k = 1, 3
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, tmp_ao)
         DO i_spin = 1, rtbse_env%n_spin
            CALL transform_ao_to_mo_covariant_fm(rtbse_env, tmp_ao, rtbse_env%moments_field(k, i_spin), i_spin)
         END DO
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
!> \brief Initial MO density ρ^0_mn = f_m δ_mn (f_m = 1 on active occupied, 0 on virtual), copied to
!>        rho_orig as the reference for the δ-kick. Imaginary part zero.
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03.26) - adapted to the linearized active-MO path
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
      ! rho_orig stays the SCF reference for delta rho; a restart overwrites only rho, downstream in
      ! the driver (read_restart_density + apply_restart_basis_bridge + rotate_rho_phase).

      CALL timestop(handle)
   END SUBROUTINE initialize_density_matrix

! **************************************************************************************************
!> \brief Single-particle reference Hamiltonian in the active-MO basis: H^0_mn = ε^GW_m δ_mn (or KS
!>        ε^scf), diagonal. TDA first-peak adds +Ω_0/2 on occupied, -Ω_0/2 on virtual diagonals.
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03.26) - refactor in prep. of linearized propagation
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
      rtbse_env%eps_active(:, :) = 0.0_dp
      ! ****** START SINGLE PARTICLE HAMILTONIAN CALCULATION
      DO i = 1, rtbse_env%n_spin
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
               IF (i_row_global == j_col_global) THEN
                  abs_mo_idx = i_row_global + rtbse_env%first_active_mo - 1
                  IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
                     ! GW Hamiltonian
                     rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = bs_env%eigenval_GW(abs_mo_idx, 1, i)
                  ELSE
                     ! KS Hamiltonian
                     rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = bs_env%eigenval_scf_Gamma(abs_mo_idx, i)
                  END IF
                  ! First-peak shift (TDA only): rotate the single-particle Hamiltonian
                  ! into a frame where the lowest active OV mode oscillates at zero, with
                  ! Ω_0 = eps_min_ai. Adds +Ω_0/2 on active occupied diagonals and
                  ! -Ω_0/2 on active virtual diagonals so that [h', rho]_OV picks up an
                  ! overall (- eps_ai + Ω_0) and OO/VV blocks remain commutator-free.
                  IF (rtbse_env%tda_active .AND. rtbse_env%tda_shift_to_first_peak) THEN
                     IF (abs_mo_idx <= rtbse_env%n_occ(i)) THEN
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                           rtbse_env%real_workspace_mo(1)%local_data(ii, jj) + 0.5_dp*rtbse_env%omega_shift
                     ELSE
                        rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                           rtbse_env%real_workspace_mo(1)%local_data(ii, jj) - 0.5_dp*rtbse_env%omega_shift
                     END IF
                  END IF
                  ! Mirror the finalized diagonal into the replicated active-energy array
                  rtbse_env%eps_active(i_row_global, i) = rtbse_env%real_workspace_mo(1)%local_data(ii, jj)
               END IF
            END DO
         END DO
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), mtarget=rtbse_env%ham_reference_singleparticle(i))
      END DO
      ! Each diagonal element was set on its single owner rank; sum to replicate eps_active.
      CALL rtbse_env%real_workspace_mo(1)%matrix_struct%para_env%sum(rtbse_env%eps_active)
      ! ****** END SINGLE PARTICLE HAMILTONIAN CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_singleparticle_hamiltonian

! **************************************************************************************************
!> \brief Reference Hartree subtraction: builds V^H[ρ^0] (AO-RI or RI-RS) and subtracts it into
!>        ham_reference, realizing H_eff = ... + V_H[ρ] - V_H[ρ^0]. Only for non-TDA n_spin=1 (in
!>        TDA the OV/VO projection of the OO-diagonal ρ^0 vanishes, so the reference is zero).
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03.26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_hartree_potential(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_hartree_potential'

      INTEGER                                            :: handle, i, n_grid
      LOGICAL                                            :: use_hartree_reference, use_rirs_kernel
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)
      use_hartree_reference = (.NOT. rtbse_env%tda_active) .AND. (rtbse_env%n_spin == 1) .AND. &
                              (.NOT. rtbse_env%debug_disable_hartree)
      use_rirs_kernel = rtbse_env%rirs_kernel

      ! Make sure the RI-RS V_grid kernel is available if we need it for Hartree.
      IF (use_rirs_kernel) CALL rt_bse_ri_rs_ensure_V_grid(bs_env, rtbse_env%qs_env)

      ! Spin-summed grid-density accumulators for the RI-RS Hartree reuse (diag(φρφ^T) harvested in SEX;
      ! mat_phi_mu_l is grid x AO, so its row count is n_grid). Allocated ONLY for RI-RS + Hartree, so the
      ! get_sigma harvest calls (passed unconditionally below) see an absent optional and self-disable on
      ! AO-RI / Hartree-off (F2008 unallocated-allocatable -> absent optional). Sole allocator, runs once.
      IF (use_rirs_kernel .AND. (.NOT. rtbse_env%debug_disable_hartree)) THEN
         CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, nfullrows_total=n_grid)
         ALLOCATE (rtbse_env%hartree_diag_re(n_grid), rtbse_env%hartree_diag_im(n_grid))
      END IF

      ! The RI-RS Hartree normally reuses the grid density harvested by the SEX kernel. With SEX
      ! disabled but Hartree on, that diagonal is never produced, so the Hartree must rebuild the full
      ! real-space density grid itself every RK4 stage - much slower. Warn once (this is a debug-only
      ! configuration); the rebuild fallback lives in the use_sex branches of the Hartree kernels.
      IF (use_rirs_kernel .AND. (.NOT. rtbse_env%debug_disable_hartree) .AND. rtbse_env%debug_disable_sex) THEN
         CALL cp_warn(__LOCATION__, &
                      "RI-RS Hartree rebuilds the full density grid every RK4 stage because SEX is "// &
                      "disabled (DEBUG_DISABLE_SEX) and no SEX-harvested diagonal is available to reuse. "// &
                      "This slows down the Hartree computation considerably.")
      END IF

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! v_dbcsr is needed by either AO-RI Hartree (here) or by AO-RI SEX (W = V + W^c assembly).
      IF (.NOT. use_rirs_kernel) THEN
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
            IF (use_rirs_kernel) THEN
               ! V^H_λσ = sum_l φ_λ(r_l) v_l φ_σ(r_l), v_l = sum_l' V_ll' n_l'  (RI-RS)
               ! AO-RI get_hartree uses only Re(rho); mirror that on the RI-RS path.
               CALL cp_cfm_to_fm(msource=rtbse_env%rho_ao_scratch(i), &
                                 mtargetr=rtbse_env%real_workspace(1))
               CALL compute_hartree_ri_rs(bs_env, rtbse_env%real_workspace(1), &
                                          rtbse_env%hartree_curr_ao(i))
            ELSE
               ! V^H_λσ = sum_PQ (λσ|P) V_PQ [sum_µν (µν|Q) ρ^0_µν]  (AO-RI; reference density)
               CALL get_hartree(rtbse_env, rtbse_env%rho_ao_scratch(i), rtbse_env%hartree_curr_ao(i))
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
!> \brief Reference SEX self-energy subtraction into ham_reference (non-TDA, n_spin=1). Assembles
!>        W = V + W^c, then Σ^SX = -W ρ^0, (f_n - f_m)-weighted, subtracted.
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03.26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_sex_selfenergy(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_sex_selfenergy'

      INTEGER                                            :: handle, i
      LOGICAL                                            :: use_rirs_kernel, use_sex_reference

      CALL timeset(routineN, handle)
      use_sex_reference = (.NOT. rtbse_env%tda_active) .AND. (rtbse_env%n_spin == 1) .AND. &
                          (.NOT. rtbse_env%debug_disable_sex)
      use_rirs_kernel = rtbse_env%rirs_kernel

      ! Make sure the RI-RS W0_grid kernel is available if we need it for SEX.
      IF (use_rirs_kernel) CALL rt_bse_ri_rs_ensure_W0_grid(rtbse_env%bs_env, rtbse_env%qs_env)

      ! ****** START SEX REFERENCE CALCULATION
      ! w_dbcsr and screened_dbt are needed for get_sigma routines (AO-RI path only).
      ! For RI-RS the W = V + W^c kernel is on the real-space grid in mat_W0_grid_rtbse.
      IF (.NOT. use_rirs_kernel) THEN
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_gw) THEN
            ! W(w=0) is built by the GW step only under its RTBSE rtp_method gate; reaching this
            ! consumer without it means gate and consumer disagree. The HF branch below needs no W.
            IF (.NOT. ASSOCIATED(rtbse_env%bs_env%fm_W_MIC_freq_zero%matrix_struct)) THEN
               CALL cp_abort(__LOCATION__, &
                             "RT-BSE AO-RI kernel needs the screened interaction W(w=0), which the "// &
                             "GW step did not build. Select the RT-BSE propagator with '&RTBSE' or "// &
                             "'&RTBSE RTBSE', not '&RTBSE TDDFT'.")
            END IF
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
      ! Calculate the SEX starting energies
      DO i = 1, rtbse_env%n_spin
         ! Calculate the exchange (SEX) part for this spin channel
         ! Uses rho_orig - same as rho for initial run but different for continued run
         ! For KS reference this is the time-dependent Fock exchange (w_dbcsr = v only).
         ! In TDA the propagator evaluates separate OV/VO-projected kernel passes.
         ! rho_orig is OO-diagonal in the MO basis, so its OV/VO projections
         ! vanish and the SEX reference must remain zero in TDA.
         IF (use_sex_reference) THEN
            ! Σ^SX_λσ = -sum_νQ [sum_µ (λµ|Q) ρ^0_µν][sum_P (νσ|P) W_PQ]  (reference)
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(i), -1.0_dp, rtbse_env%rho_ao_scratch(i))
            ! Transform to MO basis (AO scratch -> MO-active result)
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(i), rtbse_env%sigma_SEX(i), i)
            ! Apply occupation factor f_n-f_m
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(i), i)
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END IF
      END DO
      ! ****** END SEX REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_sex_selfenergy

! **************************************************************************************************
!> \brief Propagates the density one timestep by RK4 for ∂_t ρ = f(t,ρ) =
!>        -i( Δε Δρ + (f_n - f_m) V_Hartree(Δρ) + ΔΣ(Δρ) ); see body for the 4-stage scheme.
!>        Spin loop is inner to each stage (cross-spin Hartree coupling).
!> \param rtbse_env Entry point - rtbse environment
!> \param rho_start Initial density matrix
!> \param rho_end Final density matrix
! **************************************************************************************************
   SUBROUTINE solve_rk4_timestep(rtbse_env, rho_start, rho_end)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_start, rho_end

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

      ! Spin loop is INNER to each RK4 stage (inside do_rk4_stage): cross-spin-coupled kernels (Hartree in
      ! open shell) need every spin's stage density before any spin advances. rk4_coefficients(i) holds
      ! spin i's CURRENT-stage k (indexed by spin, not stage - see its allocation in create_rtbse_env),
      ! reused across stages; rho_workspace(i) holds spin i's stage density. Bit-identical for n_spin=1.
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rho_start(i), rho_end(i))
      END DO
      ! Each stage: do_rk4_stage evaluates k = f(stage density) for all spins, accumulates
      ! rho_end += result_weight*dt*k (Butcher b = 1/6, 1/3, 1/3, 1/6) and forms the next stage density
      ! rho_workspace = rho_start + advance_weight*dt*k (node c = 1/2, 1/2, 1; omitted on k_4, which only
      ! accumulates). The dt factor is applied inside do_rk4_stage, so the calls show the bare weights.
      ! k_1 = f(t, rho_start)
      CALL do_rk4_stage(rtbse_env, rho_start, rho_start, rho_end, &
                        result_weight=1.0_dp/6.0_dp, advance_weight=0.5_dp)
      ! k_2 = f(t + dt/2, rho_start + dt/2 * k_1)
      CALL do_rk4_stage(rtbse_env, rtbse_env%rho_workspace, rho_start, rho_end, &
                        result_weight=1.0_dp/3.0_dp, advance_weight=0.5_dp)
      ! k_3 = f(t + dt/2, rho_start + dt/2 * k_2)
      CALL do_rk4_stage(rtbse_env, rtbse_env%rho_workspace, rho_start, rho_end, &
                        result_weight=1.0_dp/3.0_dp, advance_weight=1.0_dp)
      ! k_4 = f(t + dt, rho_start + dt * k_3)
      CALL do_rk4_stage(rtbse_env, rtbse_env%rho_workspace, rho_start, rho_end, &
                        result_weight=1.0_dp/6.0_dp)

      ! Update bookkeeping to the next timestep similar to logic of etrs_scf_loop
      rtbse_env%sim_step = rtbse_env%sim_step + 1
      rtbse_env%sim_time = rtbse_env%sim_time + rtbse_env%sim_dt

      CALL timestop(handle)
   END SUBROUTINE solve_rk4_timestep

! **************************************************************************************************
!> \brief Takes one RK4 stage. Evaluates k = f(t, rho_eval) for every spin into
!>        rtbse_env%rk4_coefficients, accumulates it into the running result (rho_end += result_weight*dt*k)
!>        and - unless this is the last stage - forms the next stage density
!>        (rho_workspace = rho_base + advance_weight*dt*k). Builds the stage's shared kernels first: the
!>        cross-spin Hartree is built once per stage and consumed by every spin in update_effective_ham_MO.
!>        All shell/kernel combinations go through build_shared_sex_and_hartree (mask_mode selects the
!>        input convention); update_effective_ham_MO is then a pure consumer.
!>        No timeset/timestop: the callees are individually timed and this runs 4x per RK4 timestep.
!> \param rtbse_env RT-BSE environment
!> \param rho_eval Per-spin MO density f is evaluated at (rho_start for k_1, rho_workspace otherwise)
!> \param rho_base Per-spin MO density the next stage advances from (the step's rho_start)
!> \param rho_end Per-spin RK4 result accumulator (= rho_start + dt/6*(k1+2k2+2k3+k4) after all 4 stages)
!> \param result_weight RK4 Butcher weight b (dt factor applied inside) for accumulating k into rho_end
!> \param advance_weight RK4 node c (dt factor applied inside) for the next stage density; ABSENT on the last stage
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE do_rk4_stage(rtbse_env, rho_eval, rho_base, rho_end, result_weight, advance_weight)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_eval, rho_base, rho_end
      REAL(kind=dp), INTENT(IN)                          :: result_weight
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: advance_weight

      INTEGER                                            :: i, mask_mode

      ! Input convention for this stage's shared kernel build:
      IF (rtbse_env%tda_active) THEN
         mask_mode = kernel_input_ov          ! TDA (any shell): OV only, non-Hermitian
      ELSE IF (rtbse_env%n_spin > 1) THEN
         mask_mode = kernel_input_ovvo        ! open-shell ABBA: OV+VO, Hermitian
      ELSE
         mask_mode = kernel_input_full        ! closed-shell ABBA: full ρ, Hermitian
      END IF
      ! Build this stage's shared SEX + bare Hartree; consumers are per-spin in update_effective_ham_MO.
      CALL build_shared_sex_and_hartree(rtbse_env, rho_eval, mask_mode)
      DO i = 1, rtbse_env%n_spin
         CALL update_effective_ham_MO(rtbse_env, rho_eval(i), rtbse_env%rk4_coefficients(i), i)
         IF (rtbse_env%tda_active .OR. rtbse_env%n_spin > 1) THEN
            CALL project_drho_to_ov(rtbse_env, rtbse_env%rk4_coefficients(i), i)
         END IF
      END DO
      ! Fold each spin's k into the RK4 result and (unless last stage) form the next stage density.
      ! The dt factor lives here so the call sites carry the bare RK4 weights.
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rho_end(i), &
                                   CMPLX(result_weight*rtbse_env%sim_dt, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(i))
         IF (PRESENT(advance_weight)) THEN
            CALL cp_cfm_to_cfm(rho_base(i), rtbse_env%rho_workspace(i))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_workspace(i), &
                                      CMPLX(advance_weight*rtbse_env%sim_dt, 0.0_dp, kind=dp), rtbse_env%rk4_coefficients(i))
         END IF
      END DO
   END SUBROUTINE do_rk4_stage

! **************************************************************************************************
!> \brief Builds the spin-summed complex Hartree potential in AO, once per RK4 stage (AO-RI or RI-RS).
!>        rho_total = spin_degeneracy * sum_sigma (OV-masked rho^sigma -> AO); hartree_total_ao =
!>        V_H[rho_total]. update_effective_ham_MO back-transforms it with each spin's C, so this
!>        single call feeds every spin block - the cross-spin Hartree coupling of open shell.
!>        Side effect: leaves rho_ao_scratch(sigma) = OV-masked AO density (recomputed per spin in
!>        update_effective_ham_MO; built here only to form the sum).
!> \param rtbse_env RT-BSE environment
!> \param rho_stage Per-spin MO density at the current RK4 stage
!> \param keep_ovvo .FALSE. = OV source mask (TDA); .TRUE. = OV+VO (open-shell ABBA).
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE build_shared_hartree_ao(rtbse_env, rho_stage, keep_ovvo)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_stage
      LOGICAL, INTENT(IN)                                :: keep_ovvo

      CHARACTER(len=*), PARAMETER :: routineN = 'build_shared_hartree_ao'

      INTEGER                                            :: handle, isp

      CALL timeset(routineN, handle)
      ! Per spin: mask the stage density and project to AO.
      DO isp = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rho_stage(isp), rtbse_env%rho_delta_mo(isp))
         CALL mask_mo_block_cfm(rtbse_env, rtbse_env%rho_delta_mo(isp), isp, &
                                keep_OV=.TRUE., keep_ovvo=keep_ovvo)
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_delta_mo(isp), &
                                                   rtbse_env%rho_ao_scratch(isp), isp)
      END DO
      ! Sum: rho_total = spin_degeneracy * sum_isp rho_ao_scratch(isp).
      CALL cp_cfm_set_all(rtbse_env%rho_total_ao_scratch, CMPLX(0.0_dp, 0.0_dp, kind=dp))
      DO isp = 1, rtbse_env%n_spin
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_total_ao_scratch, &
                                   CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), &
                                   rtbse_env%rho_ao_scratch(isp))
      END DO
      ! One complex Hartree contraction on the summed density. RI backend = kernel choice: both
      ! get_hartree_complex and compute_hartree_ri_rs_complex take AO x AO in/out and are spin-blind,
      ! so the spin-summed cross-spin density feeds whichever kernel is active.
      IF (rtbse_env%rirs_kernel) THEN
         ! V^H_λσ = sum_l φ_λ(r_l) v_l φ_σ(r_l), v_l = sum_l' V_ll' n_l'[ρ^total]  (RI-RS)
         CALL compute_hartree_ri_rs_complex(rtbse_env%bs_env, rtbse_env%rho_total_ao_scratch, &
                                            rtbse_env%hartree_total_ao)
      ELSE
         ! V^H_λσ = sum_PQ (λσ|P) V_PQ [sum_µν (µν|Q) ρ^total_µν]  (AO-RI)
         CALL get_hartree_complex(rtbse_env, rtbse_env%rho_total_ao_scratch, &
                                  rtbse_env%hartree_total_ao, 1)
      END IF
      CALL timestop(handle)
   END SUBROUTINE build_shared_hartree_ao

! **************************************************************************************************
!> \brief Unified cross-spin kernel builder, once per RK4 stage for every shell. Computes the
!>        per-spin SEX self-energy (stashed in sigma_SEX_ao(σ)) and the single shared bare Hartree
!>        (hartree_total_ao), so update_effective_ham_MO only consumes them.
!>        mask_mode selects the source-density convention (OV / OV+VO / full-ρ); it also determines
!>        the input Hermiticity, which gates the Hartree imaginary channel: Im computed only when
!>        mask_mode = kernel_input_ov (non-Hermitian TDA input); for Hermitian input Im ≡ 0 analytically
!>        and is skipped. Hartree emitted BARE (no spin_degeneracy); consumer scales g on the MO output.
!> \param rtbse_env RT-BSE environment
!> \param rho_stage Per-spin MO density at the current RK4 stage
!> \param mask_mode Input-convention selector: kernel_input_ov / kernel_input_ovvo / kernel_input_full
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE build_shared_sex_and_hartree(rtbse_env, rho_stage, mask_mode)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_stage
      INTEGER, INTENT(IN)                                :: mask_mode

      CHARACTER(len=*), PARAMETER :: routineN = 'build_shared_sex_and_hartree'

      INTEGER                                            :: handle, isp
      LOGICAL                                            :: harvest_im, use_hartree, &
                                                            use_rirs_kernel, use_sex

      CALL timeset(routineN, handle)
      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex
      use_rirs_kernel = rtbse_env%rirs_kernel
      ! Im channel iff the masked input is non-Hermitian (OV-only, TDA). For Hermitian input (OV+VO
      ! or full ρ) Im(ρ) is antisymmetric and Coulomb factors are symmetric, so V_H[Im] ≡ 0.
      harvest_im = (mask_mode == kernel_input_ov)

      ! Pre-zero the RI-RS spin-summed grid-diagonal accumulators (the SEX harvest target).
      IF (use_rirs_kernel .AND. use_sex .AND. use_hartree) THEN
         rtbse_env%hartree_diag_re(:) = 0.0_dp
         IF (harvest_im) rtbse_env%hartree_diag_im(:) = 0.0_dp
      END IF
      ! Per spin: stage ρ into AO (masked per mask_mode); SEX (stash sigma_SEX_ao(σ)); RI-RS harvests
      ! the diagonal.
      DO isp = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rho_stage(isp), rtbse_env%rho_delta_mo(isp))
         SELECT CASE (mask_mode)
         CASE (kernel_input_ov)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%rho_delta_mo(isp), isp, &
                                   keep_OV=.TRUE., keep_ovvo=.FALSE.)
         CASE (kernel_input_ovvo)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%rho_delta_mo(isp), isp, &
                                   keep_OV=.TRUE., keep_ovvo=.TRUE.)
         CASE (kernel_input_full)
            ! no mask: full ρ (closed-shell ABBA; reference subtracted later via ham_reference)
         CASE DEFAULT
            CPABORT("Unknown mask_mode in build_shared_sex_and_hartree")
         END SELECT
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%rho_delta_mo(isp), &
                                                   rtbse_env%rho_ao_scratch(isp), isp)
         IF (use_sex) THEN
            ! Σ^SX_λσ = -sum_νQ [sum_µ (λµ|Q) Δρ_µν][sum_P (νσ|P) W_PQ]. Im accumulator passed only
            ! when harvesting; absent (unallocated optional) on AO-RI / Hartree-off / Hermitian.
            IF (harvest_im) THEN
               CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(isp), -1.0_dp, rtbse_env%rho_ao_scratch(isp), &
                              grid_diag_re_accum=rtbse_env%hartree_diag_re, &
                              grid_diag_im_accum=rtbse_env%hartree_diag_im)
            ELSE
               CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(isp), -1.0_dp, rtbse_env%rho_ao_scratch(isp), &
                              grid_diag_re_accum=rtbse_env%hartree_diag_re)
            END IF
         END IF
      END DO
      ! Shared cross-spin Hartree, emitted BARE (no spin_degeneracy — consumer scales g on output).
      ! RI-RS reuses the spin-summed grid diagonal; AO-RI contracts rho_total = sum_σ ρ_σ.
      ! Hermitian input: real Hartree path (Im ≡ 0); non-Hermitian: complex.
      IF (use_hartree) THEN
         IF (use_rirs_kernel .AND. use_sex) THEN
            IF (harvest_im) THEN
               CALL compute_hartree_ri_rs_from_diag(rtbse_env%bs_env, rtbse_env%hartree_diag_re, &
                                                    rtbse_env%hartree_total_ao, n_im=rtbse_env%hartree_diag_im)
            ELSE
               CALL compute_hartree_ri_rs_from_diag(rtbse_env%bs_env, rtbse_env%hartree_diag_re, &
                                                    rtbse_env%hartree_total_ao)
            END IF
         ELSE
            CALL cp_cfm_set_all(rtbse_env%rho_total_ao_scratch, CMPLX(0.0_dp, 0.0_dp, kind=dp))
            DO isp = 1, rtbse_env%n_spin
               CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_total_ao_scratch, &
                                         CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_ao_scratch(isp))
            END DO
            IF (use_rirs_kernel) THEN
               IF (harvest_im) THEN
                  CALL compute_hartree_ri_rs_complex(rtbse_env%bs_env, rtbse_env%rho_total_ao_scratch, &
                                                     rtbse_env%hartree_total_ao)
               ELSE
                  ! Hermitian input: real RI-RS on Re(ρ_total) -> cfm with Im ≡ 0.
                  CALL cp_cfm_to_fm(msource=rtbse_env%rho_total_ao_scratch, &
                                    mtargetr=rtbse_env%real_workspace(1))
                  CALL compute_hartree_ri_rs(rtbse_env%bs_env, rtbse_env%real_workspace(1), &
                                             rtbse_env%hartree_curr_ao(1))
                  CALL cp_cfm_set_all(rtbse_env%hartree_total_ao, CMPLX(0.0_dp, 0.0_dp, kind=dp))
                  CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr_ao(1), &
                                    mtarget=rtbse_env%hartree_total_ao)
               END IF
            ELSE
               IF (harvest_im) THEN
                  CALL get_hartree_complex(rtbse_env, rtbse_env%rho_total_ao_scratch, &
                                           rtbse_env%hartree_total_ao, 1)
               ELSE
                  ! Hermitian input: real AO-RI on Re(ρ_total) -> cfm with Im ≡ 0.
                  CALL get_hartree(rtbse_env, rtbse_env%rho_total_ao_scratch, &
                                   rtbse_env%hartree_curr_ao(1))
                  CALL cp_cfm_set_all(rtbse_env%hartree_total_ao, CMPLX(0.0_dp, 0.0_dp, kind=dp))
                  CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr_ao(1), &
                                    mtarget=rtbse_env%hartree_total_ao)
               END IF
            END IF
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE build_shared_sex_and_hartree

   ! **************************************************************************************************
!> \brief Assembles the linearized RT-BSE right-hand side in the MO basis and returns it scaled by -i:
!>        ham_effective <- -i( [H^0, ρ] + (f_n - f_m)(ΔΣ^SX[Δρ] + ΔV^H[Δρ]) ), i.e. the
!>        f(t,ρ) of ∂_t Δρ_mn = -i( (ε_m - ε_n)Δρ_mn + (f_n - f_m)(V^H_mn + Σ^SX_mn) ).
!>        ham_reference already carries KS+G0W0 minus the reference SEX/Hartree, so the kernels
!>        enter as differences vs the reference. Forks: tda_active (drop B-coupling: one OV kernel
!>        pass + VO conjugate) vs full ABBA (OV+VO); n_spin and the KERNEL_RI (rirs_kernel) flag
!>        select the AO-RI or RI-RS backend per term.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param rho Real and imaginary parts ( + spin) of the density at current time
!> \param ham_effective Effective Hamiltonian in the MO basis that is updated in this routine
!> \param ispin Spin channel σ being assembled
! **************************************************************************************************
   SUBROUTINE update_effective_ham_MO(rtbse_env, rho, ham_effective, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho, ham_effective
      INTEGER                                            :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'update_effective_ham_MO'

      INTEGER                                            :: handle, i_global, i_loc, j_global, &
                                                            j_loc, ncl, nrl
      INTEGER, DIMENSION(:), POINTER                     :: c_idx, r_idx
      LOGICAL                                            :: use_hartree, use_sex

      CALL timeset(routineN, handle)
      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex

      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0 - reference SEX - reference Hartree
      ! Sets the imaginary part to zero
      CALL cp_cfm_to_cfm(rtbse_env%ham_reference(ispin), ham_effective)
      ! [H^0, ρ]_mn = (ε_m - ε_n) ρ_mn exactly (H^0 diagonal in active-MO basis). Element-wise
      ! local pass replacing the two gemms; ham_effective and rho share fm_struct_mo_active so
      ! their local layouts coincide. Covers all (m,n) (OO/VV retained for closed-shell ABBA).
      CALL cp_cfm_get_info(matrix=ham_effective, nrow_local=nrl, ncol_local=ncl, &
                           row_indices=r_idx, col_indices=c_idx)
      DO j_loc = 1, ncl
         j_global = c_idx(j_loc)
         DO i_loc = 1, nrl
            i_global = r_idx(i_loc)
            ham_effective%local_data(i_loc, j_loc) = ham_effective%local_data(i_loc, j_loc) &
                                            + CMPLX(rtbse_env%eps_active(i_global, ispin) - rtbse_env%eps_active(j_global, ispin), &
                                                             0.0_dp, kind=dp)*rho%local_data(i_loc, j_loc)
         END DO
      END DO
      ! Determine the field at current time
      IF (rtbse_env%dft_control%apply_efield_field) THEN
         CALL cp_abort(__LOCATION__, &
                       "Continuous/pulsed E(t) field coupling is not implemented for linearized "// &
                       "RT-BSE. Only the delta-kick (impulsive) absorption spectrum is supported; "// &
                       "use APPLY_DELTA_PULSE.")
      ELSE
         ! No field
         rtbse_env%field(:) = 0.0_dp
      END IF
      IF (.NOT. rtbse_env%tda_active) THEN
         ! ===== ABBA: consume the prebuilt per-spin SEX (sigma_SEX_ao(σ)) and shared bare Hartree
         ! (hartree_total_ao). (f_n-f_m) zeros OO/VV and sets OV/VO signs.  Closed shell uses full-ρ
         ! input (no mask, reference subtracted via ham_reference); open shell uses OV+VO mask. =====
         IF (use_sex) THEN
            ! Σ^SX_λσ = -sum_νQ [sum_µ (λµ|Q) Δρ_µν][sum_P (νσ|P) W_PQ]
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                  rtbse_env%sigma_SEX(ispin), ispin)
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%sigma_SEX(ispin), ispin)
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%sigma_SEX(ispin))
         END IF
         IF (use_hartree) THEN
            ! Builder emits bare V_H (no spin_degeneracy); fold g here (post-occ-factor).
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%hartree_total_ao, &
                                                  rtbse_env%ham_workspace(1), ispin)
            CALL transform_mo_occupation_factor_diff_cfm(rtbse_env, rtbse_env%ham_workspace(1), ispin)
            CALL cp_cfm_scale(CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1))
         END IF
      ELSE
         ! ----- TDA: drop B-coupling. Consume the OV-input kernels prebuilt by
         ! build_shared_sex_and_hartree. K_MO[Δρ_VO] = (K_MO[Δρ_OV])^C (real C + symmetric AO kernels),
         ! so evaluate on OV, mask MO to OV, add VO as conjugate transpose.
         ! Signs: (f_n - f_m) = -1 on OV, +1 on VO; applied explicitly. -----

         ! SEX: AO->MO, mask OV, stash (assembled after Hartree).
         IF (use_sex) THEN
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                                  rtbse_env%sigma_SEX(ispin), ispin)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%sigma_SEX(ispin), ispin, keep_OV=.TRUE.)
         END IF

         ! Hartree: hartree_total_ao is bare (no spin_degeneracy); fold g here (post-mask).
         ! VO = (OV)^C; rho_delta_mo(ispin) is idle scratch.
         IF (use_hartree) THEN
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%hartree_total_ao, &
                                                  rtbse_env%ham_workspace(1), ispin)
            CALL mask_mo_block_cfm(rtbse_env, rtbse_env%ham_workspace(1), ispin, keep_OV=.TRUE.)
            CALL cp_cfm_scale(CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1))
            CALL cp_cfm_transpose(rtbse_env%ham_workspace(1), 'C', rtbse_env%rho_delta_mo(ispin))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_delta_mo(ispin))
         END IF

         ! SEX assembled (stashed sigma_SEX, OV-masked). OV sign -1; VO = (OV)^C sign +1.
         IF (use_sex) THEN
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%sigma_SEX(ispin))
            CALL cp_cfm_transpose(rtbse_env%sigma_SEX(ispin), 'C', rtbse_env%rho_delta_mo(ispin))
            CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), ham_effective, &
                                      CMPLX(1.0_dp, 0.0_dp, kind=dp), rtbse_env%rho_delta_mo(ispin))
         END IF

         ! Restore rho_ao_scratch to the AO image of the current full rho, so any
         ! post-routine consumer (output_mos) sees the same invariant.
         CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rho, rtbse_env%rho_ao_scratch(ispin), ispin)
      END IF
      ! Return the actual RHS f(t,rho) = -i * A(rho) for RK4
      CALL cp_cfm_scale(CMPLX(0.0_dp, -1.0_dp, kind=dp), ham_effective)

      CALL timestop(handle)
   END SUBROUTINE update_effective_ham_MO

! **************************************************************************************************
!> \brief Single-spin Liouvillian matvec: apply L^sigma_sigma to one spin's Delta rho_MO and return
!>        L * Delta rho_MO in MO basis, OV+VO blocks populated. Drives the n_spin=1 TDA diagnostic
!>        and the (n_spin=1) ABBA diagnostic; the open-shell TDA path uses the array routine
!>        apply_liouvillian_to_drho instead. Detached from the propagator (does NOT touch
!>        rho / rho_orig / ham_effective / ham_reference); only env scratches drho_probe(ispin),
!>        rho_ao_scratch, sigma_SEX_ao, hartree_total_ao, ham_workspace(1), sigma_SEX,
!>        rho_delta_mo(ispin), real_workspace_mo(1) are used.
!>
!>        Matrix elements correspond to the Casida-A matrix:
!>            L_{ia,jb} = (eps_a - eps_i) delta_{ij} delta_{ab} + (ia|jb) - W_{ij,ab}
!>        with no (f_n - f_m) factor (the propagator path applies -1 on OV; we want raw +K).
!>        Honors rtbse_env%rirs_kernel to dispatch Hartree to the RI-RS grid kernel
!>        (compute_hartree_ri_rs_complex); the SX RI-RS dispatch also reads rirs_kernel inside get_sigma.
!>        AO-RI path: get_hartree_complex / get_sigma.
!>
!> \param rtbse_env RT-BSE environment (TDA, n_spin=1).
!> \param drho_in Input Delta rho_MO (mo_active x mo_active complex).
!> \param L_drho_out Output L * Delta rho_MO; OV+VO blocks populated; OO/VV zero.
!> \param ispin Spin index.
!> \author Maximilian Graml (05.26)
! **************************************************************************************************
   SUBROUTINE apply_liouvillian_to_drho_spin(rtbse_env, drho_in, L_drho_out, ispin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: drho_in
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: L_drho_out
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_liouvillian_to_drho_spin'

      INTEGER                                            :: abs_mo_idx, handle, i_row_global, ii, &
                                                            j_col_global, jj, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: use_hartree, use_sex

      CALL timeset(routineN, handle)

      ! Mirror update_effective_ham_MO kernel gating: only Hartree + SX are gated
      ! (no COH kernel in linRTBSE).
      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex
      ! RI-RS dispatch reads rtbse_env%rirs_kernel directly below.
      ! V_grid / W0_grid are populated by initialize_hartree_potential /
      ! initialize_sex_selfenergy, which run before this diagnostic.

      ! 1. Stage drho_in into drho_probe. TDA: mask to OV (propagator carries only OV
      !    and the VO contribution is recovered later via the (.)^C shortcut). ABBA: keep
      !    the full OV+VO content - drho_in carries both blocks independently.
      CALL cp_cfm_to_cfm(drho_in, rtbse_env%drho_probe(ispin))
      IF (rtbse_env%tda_active) THEN
         CALL mask_mo_block_cfm(rtbse_env, rtbse_env%drho_probe(ispin), ispin, keep_OV=.TRUE.)
      END IF

      ! 2. Project Delta rho_OV (MO -> AO, contravariant).
      CALL transform_mo_to_ao_contravariant_cfm(rtbse_env, rtbse_env%drho_probe(ispin), &
                                                rtbse_env%rho_ao_scratch(ispin), ispin)

      ! 3. Initialize the output accumulator.
      CALL cp_cfm_set_all(L_drho_out, CMPLX(0.0_dp, 0.0_dp, kind=dp))

      ! 4. eps_{ai} diagonal contribution, kept before the kernel steps (5/6). drho_in aliases
      !    rtbse_env%drho_probe(ispin) at the caller; the kernels' VO Hermitian-conjugate scratch
      !    is rho_delta_mo(ispin) (NOT drho_probe), so they no longer corrupt drho_in - but eps
      !    first is the clean ordering. Build H_eps as a real diagonal fm in real_workspace_mo(1),
      !    convert to cfm in ham_workspace(1), accumulate [drho_in, H_eps] = drho * H - H * drho.
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
                  rtbse_env%bs_env%eigenval_GW(abs_mo_idx, 1, ispin)
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
         IF (rtbse_env%rirs_kernel) THEN
            ! V^H_λσ = sum_l φ_λ(r_l) v_l φ_σ(r_l), v_l = sum_l' V_ll' n_l'  (RI-RS)
            CALL compute_hartree_ri_rs_complex(rtbse_env%bs_env, rtbse_env%rho_ao_scratch(ispin), &
                                               rtbse_env%hartree_total_ao)
         ELSE
            ! V^H_λσ = sum_PQ (λσ|P) V_PQ [sum_µν (µν|Q) Δρ_µν]  (AO-RI)
            CALL get_hartree_complex(rtbse_env, rtbse_env%rho_ao_scratch(ispin), &
                                     rtbse_env%hartree_total_ao, ispin)
         END IF
         CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%hartree_total_ao, &
                                               rtbse_env%ham_workspace(1), ispin)
         CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%ham_workspace(1), L_drho_out, &
                                 CMPLX(rtbse_env%spin_degeneracy, 0.0_dp, kind=dp), ispin)
      END IF

      ! 6. Screened-exchange contribution: get_sigma is already complex-aware. The
      !    -1.0_dp factor passed to get_sigma builds the -W contribution; we then add
      !    +1.0 here (no occupation-factor flip), giving raw K^SX = -W as required by
      !    K = (ia|jb) - W_{ij,ab} (Casida convention with +K).
      IF (use_sex) THEN
         ! Σ^SX_λσ = -sum_νQ [sum_µ (λµ|Q) Δρ_µν][sum_P (νσ|P) W_PQ]
         CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), -1.0_dp, rtbse_env%rho_ao_scratch(ispin))
         CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(ispin), &
                                               rtbse_env%sigma_SEX(ispin), ispin)
         CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%sigma_SEX(ispin), L_drho_out, &
                                 CMPLX(1.0_dp, 0.0_dp, kind=dp), ispin)
      END IF

      CALL timestop(handle)
   END SUBROUTINE apply_liouvillian_to_drho_spin

! **************************************************************************************************
!> \brief Adds the MO-domain kernel contribution K_MO * Delta rho to L_drho_out, branching on
!>        rtbse_env%tda_active. TDA: masks K_MO to OV in place, adds scale*K_MO on OV, then
!>        builds the VO contribution as (K_MO_OV)^C and adds it (valid only under the TDA
!>        assumption Delta rho_VO = (Delta rho_OV)^C). ABBA: adds the full mo_active x mo_active
!>        K_MO directly - OV+VO blocks are computed naturally from the full OV+VO input.
!>        K_MO is INTENT(INOUT); in the TDA branch it is masked in place (treat as scratch
!>        after this call). rho_delta_mo(ispin) is used as the VO-transpose scratch in the TDA
!>        branch (NOT drho_probe, which the open-shell array driver keeps as the live probe).
!> \param rtbse_env RT-BSE environment.
!> \param K_MO Kernel contribution in MO basis (mo_active x mo_active). Scratched in TDA branch.
!> \param L_drho_out Accumulator (mo_active x mo_active).
!> \param scale Complex scale factor (spin_degeneracy for Hartree, 1.0 for SX).
!> \param ispin Spin index.
!> \author Maximilian Graml (05.26)
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
         CALL cp_cfm_transpose(K_MO, 'C', rtbse_env%rho_delta_mo(ispin))
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out, scale, rtbse_env%rho_delta_mo(ispin))
      ELSE
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out, scale, K_MO)
      END IF

      CALL timestop(handle)
   END SUBROUTINE add_K_MO_to_L_drho

! **************************************************************************************************
!> \brief Open-shell TDA Liouvillian matvec: apply the joint spin-block L_TDA to a per-spin probe
!>        and return L * Delta rho on every spin block. A probe on spin sigma feeds the diagonal
!>        block A^{sigma,sigma} (eps + Coulomb + SX) AND the off-diagonal Coulomb block
!>        A^{sigma',sigma} for every other output spin sigma' (cross-spin Hartree - the term that
!>        produces the singlet/triplet split). AO-RI only: the Phase-D guard forbids RIRS
!>        for n_spin>1, so there is no RIRS branch here; the n_spin=1 diagnostic routes
!>        through apply_liouvillian_to_drho_spin instead (which keeps the RIRS path).
!> \param rtbse_env RT-BSE environment (TDA).
!> \param drho_in Per-spin probe Delta rho_MO (mo_active x mo_active each); zero on non-probed spins.
!> \param L_drho_out Per-spin output L * Delta rho; OV+VO blocks populated, OO/VV zero.
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE apply_liouvillian_to_drho(rtbse_env, drho_in, L_drho_out)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: drho_in, L_drho_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_liouvillian_to_drho'

      INTEGER                                            :: abs_mo_idx, handle, i_row_global, ii, &
                                                            isp, isp_out, j_col_global, jj, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: use_hartree, use_sex

      CALL timeset(routineN, handle)

      use_hartree = .NOT. rtbse_env%debug_disable_hartree
      use_sex = .NOT. rtbse_env%debug_disable_sex

      ! Zero every output spin block before accumulating.
      DO isp = 1, rtbse_env%n_spin
         CALL cp_cfm_set_all(L_drho_out(isp), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      END DO

      ! eps^sigma commutator -> diagonal block only (eps is spin-diagonal). MUST precede the
      ! kernel steps (add_K_MO_to_L_drho scratches MO buffers). Lab-frame bare GW eigenvalues
      ! so they compare 1:1 to bse_full_diag.F.
      CALL cp_fm_get_info(rtbse_env%real_workspace_mo(1), &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      DO isp = 1, rtbse_env%n_spin
         CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
               IF (i_row_global == j_col_global) THEN
                  abs_mo_idx = i_row_global + rtbse_env%first_active_mo - 1
                  rtbse_env%real_workspace_mo(1)%local_data(ii, jj) = &
                     rtbse_env%bs_env%eigenval_GW(abs_mo_idx, 1, isp)
               END IF
            END DO
         END DO
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), &
                           mtarget=rtbse_env%ham_workspace(1))
         ! [drho, H_eps] = drho * H_eps - H_eps * drho
         CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                          CMPLX(1.0_dp, 0.0_dp, kind=dp), drho_in(isp), rtbse_env%ham_workspace(1), &
                          CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out(isp))
         CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                          CMPLX(-1.0_dp, 0.0_dp, kind=dp), rtbse_env%ham_workspace(1), drho_in(isp), &
                          CMPLX(1.0_dp, 0.0_dp, kind=dp), L_drho_out(isp))
      END DO

      ! Project per-spin OV-masked rho_ao(sigma) (consumed by SX) and build the single cross-spin
      ! Hartree V_H[spin_degeneracy * sum_sigma rho_ao(sigma)] (Coulomb is spin-blind).
      ! build_shared_hartree_ao does both; skip only when neither kernel is active.
      IF (use_hartree .OR. use_sex) THEN
         CALL build_shared_hartree_ao(rtbse_env, drho_in, keep_ovvo=.FALSE.)
      END IF

      ! Hartree: the one shared V_H read back with each output spin's C fills the diagonal
      ! A^{sigma,sigma} AND the off-diagonal A^{sigma',sigma}. Coeff 1.0 (spin_degeneracy lives in
      ! the summed density); no (f_n - f_m) factor (Casida +K convention).
      IF (use_hartree) THEN
         DO isp_out = 1, rtbse_env%n_spin
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%hartree_total_ao, &
                                                  rtbse_env%ham_workspace(1), isp_out)
            CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%ham_workspace(1), L_drho_out(isp_out), &
                                    CMPLX(1.0_dp, 0.0_dp, kind=dp), isp_out)
         END DO
      END IF

      ! Screened exchange: spin-diagonal (W^sigma acts only within spin sigma). get_sigma builds
      ! -W; add raw +1 (no occupation-factor flip), giving K^SX = -W.
      IF (use_sex) THEN
         DO isp = 1, rtbse_env%n_spin
            ! Σ^SX_λσ = -sum_νQ [sum_µ (λµ|Q) Δρ_µν][sum_P (νσ|P) W_PQ]  (per spin)
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX_ao(isp), -1.0_dp, rtbse_env%rho_ao_scratch(isp))
            CALL transform_ao_to_mo_covariant_cfm(rtbse_env, rtbse_env%sigma_SEX_ao(isp), &
                                                  rtbse_env%sigma_SEX(isp), isp)
            CALL add_K_MO_to_L_drho(rtbse_env, rtbse_env%sigma_SEX(isp), L_drho_out(isp), &
                                    CMPLX(1.0_dp, 0.0_dp, kind=dp), isp)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE apply_liouvillian_to_drho

! **************************************************************************************************
!> \brief Public entry for the Liouvillian eigenvalue diagnostic.
!>        Dispatches to the TDA branch (Casida-A via cp_cfm_heevd) or the ABBA branch
!>        (Furche reduction via cp_cfm_power) based on rtbse_env%tda_active.
!>        Called once at job init from run_propagation_linearized_bse, gated on
!>        rtbse_env%diagnose_liouvillian_eig (n_spin = 1 enforced at env creation).
!> \param rtbse_env RT-BSE environment with diagnostic scratch already allocated.
!> \author Maximilian Graml (05.26)
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
!> \brief TDA branch of the Liouvillian eigenvalue diagnostic. Probes the Liouvillian with
!>        canonical OV unit vectors and assembles the joint spin-block Casida-A matrix as L_pairs
!>        (N_OV_joint x N_OV_joint, spin blocks stacked), then diagonalizes via cp_cfm_heevd.
!>        A probe on spin sigma fills its column block and the response on every output spin lands
!>        in that spin's row block (the off-diagonal blocks carry the cross-spin Coulomb that
!>        splits singlet/triplet). The matvec is dispatched on n_spin: n_spin=1 uses the
!>        single-spin apply_liouvillian_to_drho_spin (keeps the RIRS path, bit-identical to the
!>        closed-shell baseline); n_spin=2 uses the AO-RI array apply_liouvillian_to_drho.
!>        Eigenvalues go to stdout (RTBSE|) and to the LIOUVILLIAN_EIG .dat file. Detached from
!>        RK4 state. Called from the dispatcher when tda_active=.TRUE.
!> \param rtbse_env RT-BSE environment with TDA diagnostic scratch already allocated.
!> \author Maximilian Graml (05.26)
! **************************************************************************************************
   SUBROUTINE diagnose_TDA_liouvillian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagnose_TDA_liouvillian'

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: ov_block
      INTEGER                                            :: b, eig_unit, handle, j, k_col, k_local, &
                                                            n, n_ov_joint, sigma_out, sigma_probe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_act_occ, n_act_virt, n_ov, off
      REAL(kind=dp)                                      :: residual_max
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      ! Per-spin OV counts + offsets into the stacked joint Liouvillian. off(1)=0,
      ! off(2)=n_ov(1); n_ov_joint = sum_sigma n_ov(sigma). n_spin=1 -> single block.
      ALLOCATE (n_act_occ(rtbse_env%n_spin), n_act_virt(rtbse_env%n_spin), &
                n_ov(rtbse_env%n_spin), off(rtbse_env%n_spin))
      n_ov_joint = 0
      DO sigma_probe = 1, rtbse_env%n_spin
         n_act_occ(sigma_probe) = rtbse_env%n_occ(sigma_probe) - rtbse_env%first_active_mo + 1
         n_act_virt(sigma_probe) = rtbse_env%last_active_mo - rtbse_env%n_occ(sigma_probe)
         n_ov(sigma_probe) = n_act_occ(sigma_probe)*n_act_virt(sigma_probe)
         off(sigma_probe) = n_ov_joint
         n_ov_joint = n_ov_joint + n_ov(sigma_probe)
      END DO

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE| ----- TDA Liouvillian diagnostic -----'
         WRITE (rtbse_env%unit_nr, '(A,I0,A,I0)') &
            ' RTBSE| n_spin = ', rtbse_env%n_spin, ', joint N_OV = ', n_ov_joint
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
      END IF

      ! Single joint block to the .dat file (REWIND). ignore_should_output=.TRUE. so this fires
      ! at init regardless of MD-iteration cadence.
      eig_unit = cp_print_key_unit_nr(logger, rtbse_env%eig_section, &
                                      extension=".dat", &
                                      file_form="FORMATTED", &
                                      file_position="REWIND", &
                                      ignore_should_output=.TRUE.)
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '# Joint spin-block TDA Liouvillian eigenvalues'
         IF (rtbse_env%n_spin == 1) THEN
            WRITE (eig_unit, '(A,I0,A,I0)') '# n_spin = ', rtbse_env%n_spin, ', N_OV = ', n_ov(1)
         ELSE
            WRITE (eig_unit, '(A,I0,A,I0,A,I0)') '# n_spin = ', rtbse_env%n_spin, &
               ', N_OV(1) = ', n_ov(1), ', N_OV(2) = ', n_ov(2)
         END IF
      END IF

      ! Assemble the joint Casida-A matrix. A probe on spin sigma_probe with canonical OV unit
      ! vector e_{(j,b)} fills column off(sigma_probe)+k_local; the response on each output spin
      ! sigma_out lands in its row block [off(sigma_out)+1 .. +n_ov(sigma_out)]. Block-level
      ! transfers keep this at O(N_OV_joint) collective ops. Column-major OV index
      ! k_local = (b_local-1)*n_act_occ + j_local matches the Fortran layout of ov_block, so
      ! RESHAPE without padding gives the right (n_ov, 1) column.
      DO sigma_probe = 1, rtbse_env%n_spin
         DO b = rtbse_env%n_occ(sigma_probe) + 1, rtbse_env%last_active_mo
            DO j = rtbse_env%first_active_mo, rtbse_env%n_occ(sigma_probe)
               k_local = (b - rtbse_env%n_occ(sigma_probe) - 1)*n_act_occ(sigma_probe) + &
                         (j - rtbse_env%first_active_mo + 1)
               k_col = off(sigma_probe) + k_local

               DO sigma_out = 1, rtbse_env%n_spin
                  CALL cp_cfm_set_all(rtbse_env%drho_probe(sigma_out), CMPLX(0.0_dp, 0.0_dp, kind=dp))
               END DO
               CALL cp_cfm_set_element(rtbse_env%drho_probe(sigma_probe), &
                                       j - rtbse_env%first_active_mo + 1, &
                                       b - rtbse_env%first_active_mo + 1, &
                                       CMPLX(1.0_dp, 0.0_dp, kind=dp))

               ! n_spin=1 keeps the single-spin matvec (RIRS-capable, bit-identical baseline);
               ! n_spin=2 is AO-RI-only (Phase-D guard) -> cross-spin array matvec.
               IF (rtbse_env%n_spin == 1) THEN
                  CALL apply_liouvillian_to_drho_spin(rtbse_env, rtbse_env%drho_probe(1), &
                                                      rtbse_env%L_drho(1), 1)
               ELSE
                  CALL apply_liouvillian_to_drho(rtbse_env, rtbse_env%drho_probe, rtbse_env%L_drho)
               END IF

               ! Stack each output spin's OV response into the joint column.
               DO sigma_out = 1, rtbse_env%n_spin
                  ALLOCATE (ov_block(n_act_occ(sigma_out), n_act_virt(sigma_out)))
                  CALL cp_cfm_get_submatrix(rtbse_env%L_drho(sigma_out), ov_block, &
                                            start_row=1, start_col=n_act_occ(sigma_out) + 1, &
                                            n_rows=n_act_occ(sigma_out), n_cols=n_act_virt(sigma_out))
                  CALL cp_cfm_set_submatrix(rtbse_env%L_pairs, &
                                            RESHAPE(ov_block, [n_ov(sigma_out), 1]), &
                                            start_row=off(sigma_out) + 1, start_col=k_col, &
                                            n_rows=n_ov(sigma_out), n_cols=1)
                  DEALLOCATE (ov_block)
               END DO
            END DO
         END DO
      END DO

      ! Hermitian residual on the joint matrix: real-orbital BSE => L real symmetric, so
      ! ||L - L^H||_max should be at the FP floor. D = L^H - L into eigvecs_pairs (idle here;
      ! heevd overwrites it), then its max-element norm via the BLACS-native pzlange path.
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

      ! Diagonalize the joint matrix. cp_cfm_heevd returns ascending real eigenvalues.
      CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, &
                        rtbse_env%eigenvalues_liouvillian)

      ! Stdout table (eV, F12.4 right-aligned to col 80, L-7) + .dat (a.u. + eV).
      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A,T26,A,T59,A)') &
            ' RTBSE|', "Excitation index n", "Excitation energy (eV)"
         DO n = 1, n_ov_joint
            WRITE (rtbse_env%unit_nr, '(A,T40,I4,T69,F12.4)') &
               ' RTBSE|', n, rtbse_env%eigenvalues_liouvillian(n)*evolt
         END DO
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
      END IF
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '#   n          Omega [a.u.]               Omega [eV]'
         DO n = 1, n_ov_joint
            WRITE (eig_unit, '(I5,4X,ES24.14E3,4X,ES24.14E3)') n, &
               rtbse_env%eigenvalues_liouvillian(n), &
               rtbse_env%eigenvalues_liouvillian(n)*evolt
         END DO
      END IF

      CALL cp_print_key_finished_output(eig_unit, logger, rtbse_env%eig_section)

      DEALLOCATE (n_act_occ, n_act_virt, n_ov, off)

      CALL timestop(handle)
   END SUBROUTINE diagnose_TDA_liouvillian

! **************************************************************************************************
!> \brief ABBA branch of the Liouvillian eigenvalue diagnostic. Assembles the joint spin-block
!>        A and B by probing apply_liouvillian_to_drho with canonical OV unit vectors: a probe on
!>        spin sigma_probe fills joint column off(sigma_probe)+k_local; each output spin's OV
!>        response -> A, VO response -> -B^* (recovered by sign-flip + conjugation). One joint
!>        Furche reduction follows: (A-B)>0 gate (independent cp_cfm_heevd), (A-B)^{1/2} via
!>        cp_cfm_power, C = (A-B)^{1/2}(A+B)(A-B)^{1/2}, cp_cfm_heevd, Ω_n = √(C). The
!>        matvec is dispatched on n_spin: n_spin=1 -> apply_liouvillian_to_drho_spin (RIRS-capable,
!>        bit-identical to the closed-shell baseline); n_spin=2 -> the AO-RI cross-spin array
!>        apply_liouvillian_to_drho. For n_spin=1 the routine reduces to the single-block path.
!>        Output: a single joint spectrum (stdout RTBSE| + LIOUVILLIAN_EIG .dat). All eigenvalues
!>        are retained, including optically dark triplet modes (the kernel-correctness gate).
!> \param rtbse_env RT-BSE environment with ABBA diagnostic scratch allocated.
!> \author Maximilian Graml (05.26)
! **************************************************************************************************
   SUBROUTINE diagnose_ABBA_liouvillian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagnose_ABBA_liouvillian'

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: ov_block, vo_block
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: B_local
      INTEGER                                            :: b, eig_unit, handle, j, k_col, k_local, &
                                                            n, n_ov_joint, sigma_out, sigma_probe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_act_occ, n_act_virt, n_ov, off
      REAL(kind=dp)                                      :: lambda_min_AmB, residual_A, residual_B
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      ! Per-spin OV counts + offsets into the stacked joint A/B (same layout as the TDA
      ! diagnostic). off(1)=0, off(2)=n_ov(1); n_ov_joint = sum_sigma n_ov(sigma). For
      ! n_spin=1 this is the single block, bit-identical to the closed-shell ABBA path.
      ALLOCATE (n_act_occ(rtbse_env%n_spin), n_act_virt(rtbse_env%n_spin), &
                n_ov(rtbse_env%n_spin), off(rtbse_env%n_spin))
      n_ov_joint = 0
      DO sigma_probe = 1, rtbse_env%n_spin
         n_act_occ(sigma_probe) = rtbse_env%n_occ(sigma_probe) - rtbse_env%first_active_mo + 1
         n_act_virt(sigma_probe) = rtbse_env%last_active_mo - rtbse_env%n_occ(sigma_probe)
         n_ov(sigma_probe) = n_act_occ(sigma_probe)*n_act_virt(sigma_probe)
         off(sigma_probe) = n_ov_joint
         n_ov_joint = n_ov_joint + n_ov(sigma_probe)
      END DO

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE| ----- ABBA Liouvillian diagnostic -----'
         WRITE (rtbse_env%unit_nr, '(A,I0,A,I0)') &
            ' RTBSE| n_spin = ', rtbse_env%n_spin, ', joint N_OV = ', n_ov_joint
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
      END IF

      eig_unit = cp_print_key_unit_nr(logger, rtbse_env%eig_section, &
                                      extension=".dat", &
                                      file_form="FORMATTED", &
                                      file_position="REWIND", &
                                      ignore_should_output=.TRUE.)
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '# Joint spin-block ABBA Liouvillian eigenvalues'
         IF (rtbse_env%n_spin == 1) THEN
            WRITE (eig_unit, '(A,I0,A,I0)') '# n_spin = ', rtbse_env%n_spin, ', N_OV = ', n_ov(1)
         ELSE
            WRITE (eig_unit, '(A,I0,A,I0,A,I0)') '# n_spin = ', rtbse_env%n_spin, &
               ', N_OV(1) = ', n_ov(1), ', N_OV(2) = ', n_ov(2)
         END IF
      END IF

      ! Assemble joint A and B. Probe on spin sigma_probe, OV pair (j,b) -> joint column
      ! k_col = off(sigma_probe)+k_local (column-major k_local, as TDA). Each output spin's
      ! OV response -> A rows [off(sigma_out)+1 ..]; VO response -> B rows (TRANSPOSE to OV
      ! layout), recovered as -B^* below. The n_spin=2 array matvec produces the full OV+VO
      ! readout (add_K_MO_to_L_drho .NOT.tda_active branch); one cross-spin Hartree fills
      ! both A^{s,s'} and B^{s,s'} Coulomb; SX stays spin-diagonal.
      DO sigma_probe = 1, rtbse_env%n_spin
         DO b = rtbse_env%n_occ(sigma_probe) + 1, rtbse_env%last_active_mo
            DO j = rtbse_env%first_active_mo, rtbse_env%n_occ(sigma_probe)
               k_local = (b - rtbse_env%n_occ(sigma_probe) - 1)*n_act_occ(sigma_probe) + &
                         (j - rtbse_env%first_active_mo + 1)
               k_col = off(sigma_probe) + k_local

               DO sigma_out = 1, rtbse_env%n_spin
                  CALL cp_cfm_set_all(rtbse_env%drho_probe(sigma_out), CMPLX(0.0_dp, 0.0_dp, kind=dp))
               END DO
               CALL cp_cfm_set_element(rtbse_env%drho_probe(sigma_probe), &
                                       j - rtbse_env%first_active_mo + 1, &
                                       b - rtbse_env%first_active_mo + 1, &
                                       CMPLX(1.0_dp, 0.0_dp, kind=dp))

               IF (rtbse_env%n_spin == 1) THEN
                  CALL apply_liouvillian_to_drho_spin(rtbse_env, rtbse_env%drho_probe(1), &
                                                      rtbse_env%L_drho(1), 1)
               ELSE
                  CALL apply_liouvillian_to_drho(rtbse_env, rtbse_env%drho_probe, rtbse_env%L_drho)
               END IF

               DO sigma_out = 1, rtbse_env%n_spin
                  ALLOCATE (ov_block(n_act_occ(sigma_out), n_act_virt(sigma_out)))
                  ALLOCATE (vo_block(n_act_virt(sigma_out), n_act_occ(sigma_out)))
                  CALL cp_cfm_get_submatrix(rtbse_env%L_drho(sigma_out), ov_block, &
                                            start_row=1, start_col=n_act_occ(sigma_out) + 1, &
                                            n_rows=n_act_occ(sigma_out), n_cols=n_act_virt(sigma_out))
                  CALL cp_cfm_set_submatrix(rtbse_env%A_mat, &
                                            RESHAPE(ov_block, [n_ov(sigma_out), 1]), &
                                            start_row=off(sigma_out) + 1, start_col=k_col, &
                                            n_rows=n_ov(sigma_out), n_cols=1)
                  ! TRANSPOSE puts vo_block in (n_act_occ, n_act_virt) layout, matching the OV pack.
                  CALL cp_cfm_get_submatrix(rtbse_env%L_drho(sigma_out), vo_block, &
                                            start_row=n_act_occ(sigma_out) + 1, start_col=1, &
                                            n_rows=n_act_virt(sigma_out), n_cols=n_act_occ(sigma_out))
                  CALL cp_cfm_set_submatrix(rtbse_env%B_mat, &
                                            RESHAPE(TRANSPOSE(vo_block), [n_ov(sigma_out), 1]), &
                                            start_row=off(sigma_out) + 1, start_col=k_col, &
                                            n_rows=n_ov(sigma_out), n_cols=1)
                  DEALLOCATE (ov_block, vo_block)
               END DO
            END DO
         END DO
      END DO

      ! Recover B from -B^* via rank-local pass on the cfm's MPI-local data (documented
      ! exception to the fm/cfm-routines-only rule); B_recovered = -CONJG(stored).
      B_local => rtbse_env%B_mat%local_data
      B_local = -CONJG(B_local)

      ! Block-symmetry residuals on the JOINT matrices: A Hermitian, B real-symmetric.
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

      ! (A-B) positivity gate on the joint matrix: heevd on a copy; abort if lambda_min < 0.
      CALL cp_cfm_to_cfm(rtbse_env%AmB_scratch, rtbse_env%L_pairs)
      CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, rtbse_env%eigenvalues_liouvillian)
      lambda_min_AmB = rtbse_env%eigenvalues_liouvillian(1)
      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A,ES16.6,A,F12.6,A)') &
            ' RTBSE| lambda_min(A - B)                       = ', &
            lambda_min_AmB, '  a.u.  (', lambda_min_AmB*evolt, ' eV)'
      END IF
      ! Hard abort by design: a non-positive (A-B) breaks the Furche reduction.
      IF (lambda_min_AmB < 0.0_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "(A - B) not positive definite - this may hint at a triplet or "// &
                       "charge-transfer instability of the reference state.")
      END IF

      ! In-place AmB_scratch -> (A-B)^{1/2} via cp_cfm_power. threshold=0 substitution
      ! codepath is unreachable here (lambda_min_AmB > 0 already enforced above).
      CALL cp_cfm_power(rtbse_env%AmB_scratch, threshold=0.0_dp, exponent=0.5_dp)

      ! C = (A-B)^{1/2}(A+B)(A-B)^{1/2}. B_mat free after step "A+/-B" -> reuse as T scratch.
      CALL cp_cfm_gemm('N', 'N', n_ov_joint, n_ov_joint, n_ov_joint, &
                       CMPLX(1.0_dp, 0.0_dp, kind=dp), &
                       rtbse_env%AmB_scratch, rtbse_env%ApB_scratch, &
                       CMPLX(0.0_dp, 0.0_dp, kind=dp), &
                       rtbse_env%B_mat)
      CALL cp_cfm_gemm('N', 'N', n_ov_joint, n_ov_joint, n_ov_joint, &
                       CMPLX(1.0_dp, 0.0_dp, kind=dp), &
                       rtbse_env%B_mat, rtbse_env%AmB_scratch, &
                       CMPLX(0.0_dp, 0.0_dp, kind=dp), &
                       rtbse_env%L_pairs)

      ! Diagonalize C -> Ω_n^2; take +√. Safety clamp on tiny-negative noise.
      CALL cp_cfm_heevd(rtbse_env%L_pairs, rtbse_env%eigvecs_pairs, rtbse_env%eigenvalues_liouvillian)
      DO n = 1, n_ov_joint
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
         DO n = 1, n_ov_joint
            WRITE (rtbse_env%unit_nr, '(A,T40,I4,T69,F12.4)') &
               ' RTBSE|', n, rtbse_env%eigenvalues_liouvillian(n)*evolt
         END DO
         WRITE (rtbse_env%unit_nr, '(A)') ' RTBSE|'
      END IF
      IF (eig_unit > 0) THEN
         WRITE (eig_unit, '(A)') '#   n          Omega [a.u.]               Omega [eV]'
         DO n = 1, n_ov_joint
            WRITE (eig_unit, '(I5,4X,ES24.14E3,4X,ES24.14E3)') n, &
               rtbse_env%eigenvalues_liouvillian(n), &
               rtbse_env%eigenvalues_liouvillian(n)*evolt
         END DO
      END IF

      CALL cp_print_key_finished_output(eig_unit, logger, rtbse_env%eig_section)

      DEALLOCATE (n_act_occ, n_act_virt, n_ov, off)

      CALL timestop(handle)
   END SUBROUTINE diagnose_ABBA_liouvillian

! **************************************************************************************************
!> \brief Covariant AO->MO transform of an operator (real): M^MO_mn = sum_µν C_µm M^AO_µν C_νn.
!>        For operator-like kernels (Σ^SX, V^H); the density uses the contravariant routine.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param fm_ao operator in the AO basis (n_ao x n_ao), input
!> \param fm_mo operator in the active-MO basis (mo_active x mo_active), output
!> \param i_spin spin channel σ; selects C_active(σ)
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_covariant_fm(rtbse_env, fm_ao, fm_mo, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: fm_ao, fm_mo
      INTEGER, INTENT(IN)                                :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_ao_to_mo_covariant_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! step 1: T_µn = sum_ν M^AO_µν C_νn  (n_ao x mo_active)
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, fm_ao, rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      ! step 2: M^MO_mn = sum_µ C_µm T_µn  (mo_active x mo_active)
      CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%ao_mo_workspace(1), &
                         0.0_dp, fm_mo)

      CALL timestop(handle)
   END SUBROUTINE transform_ao_to_mo_covariant_fm

! **************************************************************************************************
!> \brief Covariant AO->MO transform (complex): M^MO_mn = sum_µν C_µm M^AO_µν C_νn, applied to the
!>        real and imaginary AO parts separately (2x the real cost).
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param fm_ao operator in the AO basis (n_ao x n_ao) cfm, input
!> \param fm_mo operator in the active-MO basis (mo_active x mo_active) cfm, output
!> \param i_spin spin channel σ; selects C_active(σ)
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
      ! Re(M^MO)_mn = sum_µν C_µm Re(M^AO)_µν C_νn  (two gemms via ao_mo_workspace)
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%real_workspace(1), rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%ao_mo_workspace(1), &
                         0.0_dp, rtbse_env%real_workspace_mo(1))
      ! Im(M^MO)_mn = sum_µν C_µm Im(M^AO)_µν C_νn
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
!> \brief Contravariant MO->AO transform of the density (complex): Δρ^AO_µν = sum_mn C_µm Δρ^MO_mn C_νn.
!>        Density-like (expands MO indices), unlike the covariant operator transform.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param fm_mo density in the active-MO basis (mo_active x mo_active) cfm, input
!> \param fm_ao density in the AO basis (n_ao x n_ao) cfm, output
!> \param i_spin spin channel σ; selects C_active(σ)
! **************************************************************************************************
   SUBROUTINE transform_mo_to_ao_contravariant_cfm(rtbse_env, fm_mo, fm_ao, i_spin)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: fm_mo, fm_ao
      INTEGER, INTENT(IN)                                :: i_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_mo_to_ao_contravariant_cfm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Re/Im split of Δρ^MO (mo_active x mo_active) into persistent MO-sized scratch
      CALL cp_cfm_to_fm(msource=fm_mo, mtargetr=rtbse_env%real_workspace_mo(1), &
                        mtargeti=rtbse_env%real_workspace_mo(2))
      ! Re(Δρ^AO)_µν = sum_mn C_µm Re(Δρ^MO)_mn C_νn  (C·ρ via ao_mo_workspace, then ·C^T)
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%mo_active, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%real_workspace_mo(1), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      CALL parallel_gemm("N", "T", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%mo_active, &
                         1.0_dp, rtbse_env%ao_mo_workspace(1), rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%real_workspace(1))
      ! Im(Δρ^AO)_µν = sum_mn C_µm Im(Δρ^MO)_mn C_νn
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%mo_active, &
                         1.0_dp, rtbse_env%C_active(i_spin), rtbse_env%real_workspace_mo(2), &
                         0.0_dp, rtbse_env%ao_mo_workspace(1))
      CALL parallel_gemm("N", "T", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%mo_active, &
                         1.0_dp, rtbse_env%ao_mo_workspace(1), rtbse_env%C_active(i_spin), &
                         0.0_dp, rtbse_env%real_workspace(2))
      ! Reassemble into AO-sized cfm
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), &
                        msourcei=rtbse_env%real_workspace(2), mtarget=fm_ao)

      CALL timestop(handle)
   END SUBROUTINE transform_mo_to_ao_contravariant_cfm

! **************************************************************************************************
!> \brief Scales each MO-active element (m,n) by the occupation prefactor (f_n - f_m), f in {0,1}
!>        (complex): OV -> -1, VO -> +1, OO/VV -> 0. The only trace of "ρ^0 diagonal in MO".
!>        Applied to whatever kernel the caller passes (Σ^SX, V^H in MO) - bound by the caller.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param cfm MO-active kernel matrix (mo_active x mo_active) cfm, scaled in place
!> \param i_spin spin channel σ; OV/VO boundary set by n_occ(σ)
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
      ! (f_n - f_m) applied to the real part
      CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%real_workspace_mo(1), i_spin)
      ! (f_n - f_m) applied to the imaginary part
      CALL transform_mo_occupation_factor_diff_fm(rtbse_env, rtbse_env%real_workspace_mo(2), i_spin)
      ! Copy back to cfm
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), &
                        msourcei=rtbse_env%real_workspace_mo(2), &
                        mtarget=cfm)

      CALL timestop(handle)
   END SUBROUTINE transform_mo_occupation_factor_diff_cfm

! **************************************************************************************************
!> \brief Scales each MO-active element (m,n) by the occupation prefactor (f_n - f_m), f in {0,1}
!>        (real): OV -> -1, VO -> +1, OO/VV -> 0. Real-input worker for the cfm variant; the only
!>        trace of "ρ^0 diagonal in MO". Bound by the caller to the kernel being scaled.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param fm MO-active kernel matrix (mo_active x mo_active) fm, scaled in place
!> \param i_spin spin channel σ; OV/VO boundary set by n_occ(σ)
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
!> \param keep_ovvo if present and .TRUE., keep BOTH off-diagonal blocks (OV and VO) and zero
!>        OO/VV; overrides keep_OV. Absent/false reproduces the keep_OV behaviour.
! **************************************************************************************************
   SUBROUTINE mask_mo_block_cfm(rtbse_env, cfm, i_spin, keep_OV, keep_ovvo)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: cfm
      INTEGER, INTENT(IN)                                :: i_spin
      LOGICAL, INTENT(IN)                                :: keep_OV
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_ovvo

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mask_mo_block_cfm'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: handle, i_global, i_global_mo, i_local, &
                                                            j_global, j_global_mo, j_local, n_occ, &
                                                            ncol_local, nrow_local, shift
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: keep, l_keep_ovvo

      CALL timeset(routineN, handle)

      l_keep_ovvo = .FALSE.
      IF (PRESENT(keep_ovvo)) l_keep_ovvo = keep_ovvo

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
            IF (l_keep_ovvo) THEN
               ! keep both off-diagonal blocks (OV and VO); drop OO/VV
               keep = ((i_global_mo <= n_occ) .NEQV. (j_global_mo <= n_occ))
            ELSE IF (keep_OV) THEN
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
!>        OO/VV blocks are left untouched. For direction='to_lab' pass phase = -Ω_0*t;
!>        for direction='to_rotating' pass phase = +Ω_0*t. No-op when omega_shift = 0.
!>        Preserves Hermiticity since the two phases are complex conjugates of each other.
!> \param rtbse_env RT-BSE environment
!> \param rho MO-active cfm rotated in place
!> \param i_spin Spin index
!> \param phase Real phase argument (radians); typically +/- Ω_0 * t
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
               ! ρ_OV *= e^{+iφ}
               local_data(i_local, j_local) = phase_ov*local_data(i_local, j_local)
            ELSE IF (i_global_mo > n_occ .AND. j_global_mo <= n_occ) THEN
               ! ρ_VO *= e^{-iφ}
               local_data(i_local, j_local) = phase_vo*local_data(i_local, j_local)
            END IF
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE rotate_rho_phase

! **************************************************************************************************
!> \brief Build a lab-frame copy of the (possibly rotating-frame) density rho for I/O.
!>        On the TDA + symmetric-shift path returns rho_lab(t) by multiplying OV/VO by
!>        exp(-/+ i Ω_0 t). Otherwise returns a plain copy. Writes into rho_new_last
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
         ! Lab-frame rho_OV(t) = exp(+i*Ω_0*t) * rho_tilde_OV(t)
         ! (sign derived from [h_shifted, rho_tilde]_OV = (-eps_ai + Ω_0) rho_tilde_OV).
         CALL rotate_rho_phase(rtbse_env, rtbse_env%rho_new_last(i), i, rtbse_env%omega_shift*t_phys)
      END DO
      rho_lab => rtbse_env%rho_new_last

      CALL timestop(handle)
   END SUBROUTINE build_rho_lab

! **************************************************************************************************
!> \brief Bridge the restart density from the previous run's active-MO gauge into this run's:
!>        overlap-metric basis change  U_mn = sum_µν C2_µm S_µν C1_νn   (mo_active × mo_active),
!>        then  ρ_mn ← sum_pq U_mp ρ_pq U_nq   (ρ ← U ρ U^T, U real orthogonal up to FP).
!>        Exact under per-MO sign flips and degenerate-subspace rotations of the SCF solution.
!>        Diagnostics per spin: max|U−1| (total gauge correction), sign-flip count, max off-diag
!>        (degenerate rotation), max|U^T U−1| (representability loss; warn ≥1e-10, abort ≥1e-3),
!>        max|U_OV| (occ/virt mixing; warn ≥1e-6). No-op (U=1) when the two gauges agree.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE apply_restart_basis_bridge(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_restart_basis_bridge'
      COMPLEX(kind=dp), PARAMETER :: c_one = CMPLX(1.0_dp, 0.0_dp, kind=dp), &
         c_zero = CMPLX(0.0_dp, 0.0_dp, kind=dp)

      INTEGER                                            :: handle, i, i_glob, i_mo, ii, j_glob, &
                                                            j_mo, jj, n_flip, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: found
      REAL(kind=dp)                                      :: dev_ident, dev_offdiag, dev_ov, &
                                                            dev_unitary
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: u_diag
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: u_data
      TYPE(cp_fm_type)                                   :: SC_old
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: C_old

      CALL timeset(routineN, handle)

      NULLIFY (C_old)
      ALLOCATE (C_old(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_create(C_old(i), rtbse_env%fm_struct_ao_mo_active)
      END DO
      CALL read_restart_C(rtbse_env, C_old, found)
      IF (.NOT. found) THEN
         DO i = 1, rtbse_env%n_spin
            CALL cp_fm_release(C_old(i))
         END DO
         DEALLOCATE (C_old)
         CALL timestop(handle)
         RETURN
      END IF

      CALL cp_fm_create(SC_old, rtbse_env%fm_struct_ao_mo_active)
      ALLOCATE (u_diag(rtbse_env%mo_active))

      DO i = 1, rtbse_env%n_spin
         ! S C1 : [S C1]_µn = sum_ν S_µν C1_νn
         CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%mo_active, rtbse_env%n_ao, &
                            1.0_dp, rtbse_env%S_fm, C_old(i), 0.0_dp, SC_old)
         ! U = C2^T (S C1) : U_mn = sum_µ C2_µm [S C1]_µn   -> real_workspace_mo(1)
         CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%n_ao, &
                            1.0_dp, rtbse_env%C_active(i), SC_old, 0.0_dp, rtbse_env%real_workspace_mo(1))

         ! Diagnostics on U (local blocks + global MAX reduction; diagonal is gathered globally)
         CALL cp_fm_get_diag(rtbse_env%real_workspace_mo(1), u_diag)
         n_flip = COUNT(u_diag < 0.0_dp)
         CALL cp_fm_get_info(rtbse_env%real_workspace_mo(1), nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=u_data)
         dev_ident = 0.0_dp; dev_offdiag = 0.0_dp; dev_ov = 0.0_dp
         DO ii = 1, nrow_local
            i_glob = row_indices(ii)
            i_mo = i_glob + rtbse_env%first_active_mo - 1
            DO jj = 1, ncol_local
               j_glob = col_indices(jj)
               j_mo = j_glob + rtbse_env%first_active_mo - 1
               IF (i_glob == j_glob) THEN
                  dev_ident = MAX(dev_ident, ABS(u_data(ii, jj) - 1.0_dp))
               ELSE
                  dev_ident = MAX(dev_ident, ABS(u_data(ii, jj)))
                  dev_offdiag = MAX(dev_offdiag, ABS(u_data(ii, jj)))
               END IF
               IF ((i_mo <= rtbse_env%n_occ(i)) .NEQV. (j_mo <= rtbse_env%n_occ(i))) THEN
                  dev_ov = MAX(dev_ov, ABS(u_data(ii, jj)))
               END IF
            END DO
         END DO
         ! U^T U − 1 : [U^T U]_mn = sum_p U_pm U_pn   -> real_workspace_mo(2)
         CALL parallel_gemm("T", "N", rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                            1.0_dp, rtbse_env%real_workspace_mo(1), rtbse_env%real_workspace_mo(1), &
                            0.0_dp, rtbse_env%real_workspace_mo(2))
         CALL cp_fm_get_info(rtbse_env%real_workspace_mo(2), nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=u_data)
         dev_unitary = 0.0_dp
         DO ii = 1, nrow_local
            i_glob = row_indices(ii)
            DO jj = 1, ncol_local
               j_glob = col_indices(jj)
               dev_unitary = MAX(dev_unitary, ABS(u_data(ii, jj) - MERGE(1.0_dp, 0.0_dp, i_glob == j_glob)))
            END DO
         END DO
         CALL rtbse_env%real_workspace_mo(1)%matrix_struct%para_env%max(dev_ident)
         CALL rtbse_env%real_workspace_mo(1)%matrix_struct%para_env%max(dev_offdiag)
         CALL rtbse_env%real_workspace_mo(1)%matrix_struct%para_env%max(dev_ov)
         CALL rtbse_env%real_workspace_mo(1)%matrix_struct%para_env%max(dev_unitary)

         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(A,I3,A)') " RTBSE| Restart basis bridge U = C2^T S C1 (spin ", i, "):"
            WRITE (rtbse_env%unit_nr, '(A,ES12.3,A)') " RTBSE|   max |U - 1|         ", dev_ident, &
               "  (total gauge correction)"
            WRITE (rtbse_env%unit_nr, '(A,I12)') " RTBSE|   sign flips (U_ii<0)", n_flip
            WRITE (rtbse_env%unit_nr, '(A,ES12.3,A)') " RTBSE|   max offdiag |U_ij| ", dev_offdiag, &
               "  (degenerate-subspace rotation)"
            WRITE (rtbse_env%unit_nr, '(A,ES12.3,A)') " RTBSE|   max |U^T U - 1|    ", dev_unitary, &
               "  (representability loss)"
            WRITE (rtbse_env%unit_nr, '(A,ES12.3,A)') " RTBSE|   max |U_OV|         ", dev_ov, &
               "  (occ/virt structure change)"
         END IF
         IF (dev_unitary >= 1.0e-3_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "Restart basis bridge: active spaces of the two runs differ severely (|U^T U - 1| >= 1e-3)")
         END IF
         IF (dev_unitary >= 1.0e-10_dp .AND. dev_unitary < 1.0e-3_dp) THEN
            CALL cp_warn(__LOCATION__, &
                         "Restart basis bridge: representability loss above 1e-10 - active windows differ slightly.")
         END IF
         ! 1e-6 floor: benign SCF reconvergence gives ~1e-9 occ/virt gauge noise (the bridge maps it
         ! correctly either way); only a genuine occupation-structure change reaches this threshold.
         IF (dev_ov >= 1.0e-6_dp) THEN
            CALL cp_warn(__LOCATION__, &
                         "Restart basis bridge: occupied/virtual mixing above 1e-6 - occupation structure changed.")
         END IF

         ! ρ ← U ρ U^T : lift U to complex, [Uρ]_mn = sum_p U_mp ρ_pn, then ρ_mn = sum_q [Uρ]_mq U_nq
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace_mo(1), mtarget=rtbse_env%rho_workspace(1))
         CALL cp_cfm_gemm('N', 'N', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                          c_one, rtbse_env%rho_workspace(1), rtbse_env%rho(i), c_zero, rtbse_env%rho_workspace(2))
         CALL cp_cfm_gemm('N', 'C', rtbse_env%mo_active, rtbse_env%mo_active, rtbse_env%mo_active, &
                          c_one, rtbse_env%rho_workspace(2), rtbse_env%rho_workspace(1), c_zero, rtbse_env%rho(i))
      END DO

      DEALLOCATE (u_diag)
      CALL cp_fm_release(SC_old)
      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_release(C_old(i))
      END DO
      DEALLOCATE (C_old)
      CALL timestop(handle)
   END SUBROUTINE apply_restart_basis_bridge

! **************************************************************************************************
!> \brief Complex-linear Hartree contraction. Calls the real-input get_hartree on
!>        Re(rho_AO) and on Im(rho_AO) separately and assembles
!>           v_AO = V_H[Re(rho_AO)] + i * V_H[Im(rho_AO)] .
!>        Required by the TDA propagator where the per-pass input Delta rho_OV (or
!>        Delta rho_VO) is non-Hermitian, so the imaginary part must be carried.
!>        The real kernel get_hartree realises V^H_λσ = sum_PQ (λσ|P) V_PQ [sum_µν (µν|Q) Δρ_µν].
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

      ! V^H[Re(Δρ^AO)] -> real_workspace(1)
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
!> \brief δ-kick (Marek2025) seeding the linearized EOM: builds the MO-active dipole operator
!>        A = intensity * sum_k kvec_k r_k (MO basis), Hermitizes it, and propagates ρ by exp(-iA),
!>        so Δρ^+_nm = i(f_n - f_m) A_nm excites only the OV/VO blocks.
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
      ! Per-spin kick: each spin uses its own MO-active dipole operator (C_active(i_spin) basis)
      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_set_all(rtbse_env%real_workspace_mo(1), 0.0_dp)
         DO k = 1, 3
            CALL cp_fm_scale_and_add(1.0_dp, rtbse_env%real_workspace_mo(1), &
                                     kvec(k), rtbse_env%moments_field(k, i))
         END DO
         ! enforce hermiticity of the effective Hamiltonian
         CALL cp_fm_transpose(rtbse_env%real_workspace_mo(1), rtbse_env%real_workspace_mo(2))
         CALL cp_fm_scale_and_add(0.5_dp, rtbse_env%real_workspace_mo(1), &
                                  0.5_dp, rtbse_env%real_workspace_mo(2))
         ! multiply by intensity, set as the imaginary exponent for this spin
         CALL cp_fm_scale(intensity, rtbse_env%real_workspace_mo(1))
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

      ! keep OV/VO, zero OO/VV: project Δρ onto the δ-kick sectors
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
!> \brief Per-spin electron numbers from the MO density: N_e^σ = spin_degeneracy * Re Tr[ρ^σ],
!>        returned as one entry per spin channel (alpha/beta). The node-local diagonal partial
!>        sums are reduced over the BLACS grid before scaling; the imaginary trace is a
!>        non-Hermiticity diagnostic.
!> \param rtbse_env Entry point - rtbse environment
!> \param rho Density matrix in MO basis (per spin)
!> \param electron_n_re Real electron number per spin channel (size n_spin)
!> \param electron_n_im Imaginary electron number per spin channel (numerical non-hermiticity)
! **************************************************************************************************
   SUBROUTINE get_electron_number_MO(rtbse_env, rho, electron_n_re, electron_n_im)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      REAL(kind=dp), DIMENSION(:), INTENT(OUT)           :: electron_n_re, electron_n_im

      CHARACTER(len=*), PARAMETER :: routineN = 'get_electron_number_MO'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER                                            :: handle, i_global, i_local, j, j_global, &
                                                            j_local, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)
      electron_n_re(:) = 0.0_dp
      electron_n_im(:) = 0.0_dp
      DO j = 1, rtbse_env%n_spin
         CALL cp_cfm_get_info(matrix=rho(j), &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)
         local_data => rho(j)%local_data
         ! accumulate Tr[ρ^σ] = sum_m ρ^σ_mm (real and imaginary parts separately)
         DO i_local = 1, nrow_local
            i_global = row_indices(i_local)
            ! Search column indices for the diagonal position
            DO j_local = 1, ncol_local
               j_global = col_indices(j_local)
               IF (j_global == i_global) THEN
                  ! Found diagonal element
                  electron_n_re(j) = electron_n_re(j) + REAL(local_data(i_local, j_local), kind=dp)
                  electron_n_im(j) = electron_n_im(j) + AIMAG(local_data(i_local, j_local))
                  EXIT
               END IF
            END DO
         END DO
         ! reduce the per-rank partial traces over the process grid (the MO diagonal is distributed)
         CALL rho(j)%matrix_struct%para_env%sum(electron_n_re(j))
         CALL rho(j)%matrix_struct%para_env%sum(electron_n_im(j))
         ! N_e^σ = spin_degeneracy * Tr[ρ^σ] (g=2 closed shell; 1 per channel open shell)
         electron_n_re(j) = electron_n_re(j)*rtbse_env%spin_degeneracy
         electron_n_im(j) = electron_n_im(j)*rtbse_env%spin_degeneracy
      END DO

      CALL timestop(handle)
   END SUBROUTINE get_electron_number_MO

END MODULE rt_bse_linearized

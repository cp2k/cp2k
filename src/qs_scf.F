!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the Quickstep SCF run.
!> \par History
!>      - Joost VandeVondele (02.2002)
!>           added code for: incremental (pab and gvg) update
!>                            initialisation (init_cube, l_info)
!>      - Joost VandeVondele (02.2002)
!>           called the poisson code of the classical part
!>           this takes into account the spherical cutoff and allows for
!>           isolated systems
!>      - Joost VandeVondele (02.2002)
!>           added multiple grid feature
!>           changed to spherical cutoff consistently (?)
!>           therefore removed the gradient correct functionals
!>      - updated with the new QS data structures (10.04.02,MK)
!>      - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!>      - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!>      - set_mo_occupation for smearing of the MO occupation numbers
!>        (17.04.02,MK)
!>      - MO level shifting added (22.04.02,MK)
!>      - Usage of TYPE mo_set_p_type
!>      - Joost VandeVondele (05.2002)
!>            added cholesky based diagonalisation
!>      - 05.2002 added pao method [fawzi]
!>      - parallel FFT (JGH 22.05.2002)
!>      - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!>      - started to include more LSD (01.2003,Joost VandeVondele)
!>      - 02.2003 scf_env [fawzi]
!>      - got rid of nrebuild (01.2004, Joost VandeVondele)
!>      - 10.2004 removed pao [fawzi]
!>      - 03.2006 large cleaning action [Joost VandeVondele]
!>      - High-spin ROKS added (05.04.06,MK)
!>      - Mandes (10.2013)
!>        intermediate energy communication with external communicator added
!>      - kpoints (08.2014, JGH)
!>      - unified k-point and gamma-point code (2014.11) [Ole Schuett]
!>      - added extra SCF loop for CDFT constraints (12.2015) [Nico Holmberg]
!> \author Matthias Krack (30.04.2001)
! **************************************************************************************************
MODULE qs_scf
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, &
        dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release_p, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_norm
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        cdft2ot, cholesky_dbcsr, diag_update_method_adiis, history_guess, ot2cdft, &
        ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
        outer_scf_becke_constraint, outer_scf_hirshfeld_constraint, outer_scf_optimizer_broyden, &
        outer_scf_optimizer_newton_ls, restart_guess, sic_eo, smear_fermi_dirac, smear_gaussian, &
        smear_mp, smear_mv, tblite_scc_mixer_tblite
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_io,                       ONLY: write_kpoints_restart
   USE kpoint_methods,                  ONLY: kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos,&
                                              kpoint_set_mo_occupation
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: invert_matrix
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE preconditioner,                  ONLY: make_preconditioner_complex_full_all,&
                                              make_preconditioner_complex_full_kinetic,&
                                              make_preconditioner_complex_full_s_inverse,&
                                              make_preconditioner_complex_full_single_inverse,&
                                              prepare_preconditioner,&
                                              restart_preconditioner
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
   USE qs_cdft_scf_utils,               ONLY: build_diagonal_jacobian,&
                                              create_tmp_logger,&
                                              initialize_inverse_jacobian,&
                                              prepare_jacobian_stencil,&
                                              print_inverse_jacobian,&
                                              restart_inverse_jacobian
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_charge_mixing,                ONLY: charge_mixing_scc_error
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_density_mixing_types,         ONLY: gspace_mixing_nr
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_clear_kp,&
                                              qs_diis_b_create,&
                                              qs_diis_b_create_kp
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_methods,                   ONLY: evaluate_core_matrix_traces,&
                                              qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release
   USE qs_mo_io,                        ONLY: write_mo_set_to_restart
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_restrict,&
                                              mo_set_type,&
                                              reassign_allocated_mos
   USE qs_ot,                           ONLY: qs_ot_get_orbitals_ref_complex,&
                                              qs_ot_get_p,&
                                              qs_ot_get_p_complex,&
                                              qs_ot_new_preconditioner
   USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                              ot_scf_init,&
                                              ot_scf_read_input
   USE qs_ot_types,                     ONLY: &
        qs_ot_allocate, qs_ot_channel_index, qs_ot_check_channel_context, qs_ot_init, &
        qs_ot_kpoint_preconditioner_solver_supported, qs_ot_kpoint_preconditioner_supported, &
        qs_ot_number_of_channels, qs_ot_set_context, qs_ot_type
   USE qs_outer_scf,                    ONLY: outer_loop_gradient,&
                                              outer_loop_optimize,&
                                              outer_loop_purge_history,&
                                              outer_loop_switch,&
                                              outer_loop_update_qs_env
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals
   USE qs_rho_atom_types,               ONLY: zero_rho_atom_integrals
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp,&
                                              qs_scf_get_ot_kpoint_operator
   USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
   USE qs_scf_loop_utils,               ONLY: &
        qs_scf_candidate_density_delta, qs_scf_check_inner_exit, qs_scf_check_outer_exit, &
        qs_scf_commit_density_candidate, qs_scf_density_mixing, qs_scf_inner_finalize, &
        qs_scf_kp_search_endpoint, qs_scf_new_mos, qs_scf_new_mos_kp, qs_scf_rho_update, &
        qs_scf_set_loop_flags
   USE qs_scf_oda,                      ONLY: qs_scf_oda_apply
   USE qs_scf_output,                   ONLY: qs_scf_cdft_info,&
                                              qs_scf_cdft_initial_info,&
                                              qs_scf_gce_info,&
                                              qs_scf_loop_info,&
                                              qs_scf_loop_print,&
                                              qs_scf_outer_loop_info,&
                                              qs_scf_write_mos
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_subspace,                 ONLY: qs_scf_subspace_build,&
                                              qs_scf_subspace_push,&
                                              qs_scf_subspace_restart
   USE qs_scf_subspace_types,           ONLY: qs_scf_subspace_buffer_clear,&
                                              qs_scf_subspace_buffer_create,&
                                              qs_scf_subspace_buffer_release
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_scf_env_type, &
        smeagol_method_nr, special_diag_method_nr
   USE qs_vxc_atom,                     ONLY: gapw_cdft_one_center
   USE qs_wf_history_methods,           ONLY: wfi_purge_history,&
                                              wfi_update
   USE scf_control_types,               ONLY: scf_control_type
   USE smeagol_interface,               ONLY: run_smeagol_bulktrans,&
                                              run_smeagol_emtrans
   USE tblite_interface,                ONLY: tb_get_energy,&
                                              tb_native_scc_mixer_active,&
                                              tb_scf_mixer_error,&
                                              tb_update_charges
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'
   LOGICAL, PRIVATE                     :: reuse_precond = .FALSE.
   LOGICAL, PRIVATE                     :: used_history = .FALSE.

   PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf, cdft_scf, init_scf_loop

CONTAINS

! **************************************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param has_converged ...
!> \param total_scf_steps ...
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
!> \note
! **************************************************************************************************
   SUBROUTINE scf(qs_env, has_converged, total_scf_steps)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_converged
      INTEGER, INTENT(OUT), OPTIONAL                     :: total_scf_steps

      INTEGER                                            :: ihistory, max_scf_tmp, tsteps
      LOGICAL                                            :: converged, outer_scf_loop, should_stop
      LOGICAL, SAVE                                      :: first_step_flag = .TRUE.
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, variable_history
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))
      IF (PRESENT(has_converged)) THEN
         has_converged = .FALSE.
      END IF
      IF (PRESENT(total_scf_steps)) THEN
         total_scf_steps = 0
      END IF
      CALL get_qs_env(qs_env, scf_env=scf_env, input=input, &
                      dft_control=dft_control, scf_control=scf_control)
      IF (scf_control%max_scf > 0) THEN

         dft_section => section_vals_get_subs_vals(input, "DFT")
         scf_section => section_vals_get_subs_vals(dft_section, "SCF")

         IF (.NOT. ASSOCIATED(scf_env)) THEN
            CALL qs_scf_env_initialize(qs_env, scf_env)
            ! Moved here from qs_scf_env_initialize to be able to have more scf_env
            CALL set_qs_env(qs_env, scf_env=scf_env)
         ELSE
            CALL qs_scf_env_initialize(qs_env, scf_env)
         END IF

         IF ((scf_control%density_guess == history_guess) .AND. (first_step_flag)) THEN
            max_scf_tmp = scf_control%max_scf
            scf_control%max_scf = 1
            outer_scf_loop = scf_control%outer_scf%have_scf
            scf_control%outer_scf%have_scf = .FALSE.
         END IF

         IF (.NOT. dft_control%qs_control%cdft) THEN
            CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ELSE
            ! Third SCF loop needed for CDFT with OT to properly restart OT inner loop
            CALL cdft_scf(qs_env=qs_env, should_stop=should_stop)
         END IF

         ! If SCF has not converged, then we should not start MP2
         IF (ASSOCIATED(qs_env%mp2_env)) qs_env%mp2_env%hf_fail = .NOT. converged

         ! Add the converged outer_scf SCF gradient(s)/variable(s) to history
         IF (scf_control%outer_scf%have_scf) THEN
            ihistory = scf_env%outer_scf%iter_count
            CALL get_qs_env(qs_env, gradient_history=gradient_history, &
                            variable_history=variable_history)
            ! We only store the latest two values
            gradient_history(:, 1) = gradient_history(:, 2)
            gradient_history(:, 2) = scf_env%outer_scf%gradient(:, ihistory)
            variable_history(:, 1) = variable_history(:, 2)
            variable_history(:, 2) = scf_env%outer_scf%variables(:, ihistory)
            ! Reset flag
            IF (used_history) used_history = .FALSE.
            ! Update a counter and check if the Jacobian should be deallocated
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               scf_control%outer_scf%cdft_opt_control%ijacobian(2) = scf_control%outer_scf%cdft_opt_control%ijacobian(2) + 1
               IF (scf_control%outer_scf%cdft_opt_control%ijacobian(2) >= &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) .AND. &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) > 0) THEN
                  scf_env%outer_scf%deallocate_jacobian = .TRUE.
               END IF
            END IF
         END IF
         !   *** add the converged wavefunction to the wavefunction history
         IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
             ((scf_control%density_guess /= history_guess) .OR. &
              (.NOT. first_step_flag))) THEN
            IF (.NOT. dft_control%qs_control%cdft) THEN
               CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
            ELSE
               IF (dft_control%qs_control%cdft_control%should_purge) THEN
                  CALL wfi_purge_history(qs_env)
                  CALL outer_loop_purge_history(qs_env)
                  dft_control%qs_control%cdft_control%should_purge = .FALSE.
               ELSE
                  CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
               END IF
            END IF
         ELSE IF ((scf_control%density_guess == history_guess) .AND. &
                  (first_step_flag)) THEN
            scf_control%max_scf = max_scf_tmp
            scf_control%outer_scf%have_scf = outer_scf_loop
            first_step_flag = .FALSE.
         END IF

         ! *** compute properties that depend on the converged wavefunction
         IF (.NOT. (should_stop)) CALL qs_scf_compute_properties(qs_env)

         ! *** SMEAGOL interface ***
         IF (.NOT. (should_stop)) THEN
            ! compute properties that depend on the converged wavefunction ..
            CALL run_smeagol_emtrans(qs_env, last=.TRUE., iter=0)
            ! .. or save matrices related to bulk leads
            CALL run_smeagol_bulktrans(qs_env)
         END IF

         ! *** cleanup
         CALL scf_env_cleanup(scf_env)
         IF (dft_control%qs_control%cdft) THEN
            CALL cdft_control_cleanup(dft_control%qs_control%cdft_control)
         END IF

         IF (PRESENT(has_converged)) THEN
            has_converged = converged
         END IF
         IF (PRESENT(total_scf_steps)) THEN
            total_scf_steps = tsteps
         END IF

      END IF

   END SUBROUTINE scf

! **************************************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param should_stop ...
!> \param total_scf_steps ...
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!>      10.2019 switch_surf_dip [SGh]
!> \author Matthias Krack
!> \note
! **************************************************************************************************
   SUBROUTINE scf_env_do_scf(scf_env, scf_control, qs_env, converged, should_stop, total_scf_steps)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: converged, should_stop
      INTEGER, INTENT(OUT)                               :: total_scf_steps

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scf_env_do_scf'
      INTEGER, PARAMETER :: max_ot_kp_subspace_refreshes = 4
      REAL(KIND=dp), PARAMETER :: adiis_stagnation_ratio = 1.0E-2_dp

      CHARACTER(LEN=default_string_length)               :: description, name
      INTEGER :: accepted_ot_kp_searches, ext_master_id, handle, handle2, i_tmp, ic, irep, ispin, &
         iter_count, n_rep_hf, oda_base_slot, ot_kp_subspace_refresh_count, &
         ot_kp_subspace_refresh_iter_count, output_unit, scf_energy_message_tag, total_steps
      LOGICAL :: ace_active, added_mos_auto_grow, adiis_pushed, adiis_restarted, adiis_stagnated, &
         adiis_validation, density_full_step, diis_step, do_kpoints, energy_only, exit_inner_loop, &
         exit_outer_loop, inner_loop_converged, internal_tblite_density_full_step, &
         internal_tblite_mixer, just_energy, oda_applied, oda_eligible, oda_globalization_active, &
         oda_interior, oda_predictor_valid, oda_restart_trigger, oda_state_evaluated, oda_trigger, &
         ot_kp_subspace_refresh, ot_kp_subspace_refresh_pending, outer_loop_converged, &
         tblite_native_mixer
      REAL(KIND=dp)                                      :: adiis_raw_delta, adiis_step_delta, &
                                                            oda_predictor_lambda, t1, t2
      REAL(KIND=dp), DIMENSION(3)                        :: res_val_3
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_last_converged
      TYPE(mp_comm_type)                                 :: external_comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type)                               :: oda_base_energy, oda_evaluated_energy
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: ace_section, dft_section, hfx_sections, &
                                                            input, scf_section

! Weak orbital metrics can need several bounded REF rebuilds to close the canonical density tail.

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, rho, energy, &
               logger, qs_charges, ks_env, mos, atomic_kind_set, qs_kind_set, &
               particle_set, ace_section, dft_section, hfx_sections, input, &
               scf_section, para_env, results, kpoints, pw_env, matrix_ks, &
               matrix_ks_kp, rho_ao_kp, mos_last_converged)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      energy=energy, &
                      particle_set=particle_set, &
                      qs_charges=qs_charges, &
                      ks_env=ks_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      mos=mos, &
                      matrix_ks_kp=matrix_ks_kp, &
                      input=input, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      results=results, &
                      pw_env=pw_env, &
                      para_env=para_env)
      tblite_native_mixer = dft_control%qs_control%xtb_control%do_tblite .AND. &
                            scf_env%method /= ot_method_nr .AND. &
                            tb_native_scc_mixer_active(dft_control)
      internal_tblite_mixer = (dft_control%qs_control%dftb .AND. &
                               dft_control%qs_control%dftb_control%tblite_scc_mixer == tblite_scc_mixer_tblite) .OR. &
                              (dft_control%qs_control%xtb .AND. &
                               .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                               dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite)
      internal_tblite_density_full_step = dft_control%qs_control%xtb .AND. &
                                          .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                                          dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      ! ODA uses the STANDARD diagonalization endpoint for Gamma-point and k-point runs.
      ! Smearing remains excluded by ADIIS; other ADIIS-compatible variants are admitted
      ! only when their raw endpoint is a descent direction and its evaluated trial lowers the energy.
      oda_eligible = scf_control%diagonalization%update_method == diag_update_method_adiis .AND. &
                     scf_control%diagonalization%do_oda .AND. &
                     scf_env%method == general_diag_method_nr .AND. &
                     .NOT. (dft_control%smear .OR. scf_control%smear%do_smear)
      IF (oda_eligible) THEN
         ace_active = .FALSE.
         hfx_sections => section_vals_get_subs_vals(dft_section, "XC%HF")
         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
         DO irep = 1, n_rep_hf
            ace_section => section_vals_get_subs_vals(hfx_sections, "ACE", i_rep_section=irep)
            IF (ASSOCIATED(ace_section)) THEN
               CALL section_vals_val_get(ace_section, "ACTIVE", l_val=ace_active)
               IF (ace_active) EXIT
            END IF
         END DO
         IF (ace_active) THEN
            CALL cp_warn(__LOCATION__, &
                         "ODA is disabled because HF%ACE is active; ODA trial densities have no matching "// &
                         "MOs for the ACE projector.")
            oda_eligible = .FALSE.
         END IF
      END IF

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      IF (scf_control%gce%do_gce .AND. output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,78('-'))")
         WRITE (UNIT=output_unit, FMT="(T31,A)") &
            "GRAND-CANONICAL SCF"
         WRITE (UNIT=output_unit, FMT="(T20,A,F12.6,A)") &
            "Target work function (TWF):", &
            evolt*scf_control%gce%target_workfunction, " eV"
         WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
      END IF

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
         "SCF WAVEFUNCTION OPTIMIZATION"

      ! when switch_surf_dip is switched on, indicate storing mos from the last converged step
      IF (dft_control%switch_surf_dip) THEN
         CALL get_qs_env(qs_env, mos_last_converged=mos_last_converged)
         DO ispin = 1, dft_control%nspins
            CALL reassign_allocated_mos(mos(ispin), mos_last_converged(ispin))
         END DO
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "COPIED mos_last_converged ---> mos"
      END IF

      IF ((output_unit > 0) .AND. (.NOT. scf_control%use_ot)) THEN
         WRITE (UNIT=output_unit, &
                FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)") &
            "Step", "Update method", "Time", "Convergence", "Total energy", "Change", &
            REPEAT("-", 78)
      END IF
      CALL cp_add_iter_level(logger%iter_info, "QS_SCF")

      ! check for external communicator and if the intermediate energy should be sent
      res_val_3(:) = -1.0_dp
      description = "[EXT_SCF_ENER_COMM]"
      IF (test_for_result(results, description=description)) THEN
         CALL get_results(results, description=description, &
                          values=res_val_3, n_entries=i_tmp)
         CPASSERT(i_tmp == 3)
         IF (ALL(res_val_3(:) <= 0.0)) THEN
            CALL cp_abort(__LOCATION__, &
                          " Trying to access result ("//TRIM(description)// &
                          ") which is not correctly stored.")
         END IF
         CALL external_comm%set_handle(NINT(res_val_3(1)))
      END IF
      ext_master_id = NINT(res_val_3(2))
      scf_energy_message_tag = NINT(res_val_3(3))

      ! *** outer loop of the scf, can treat other variables,
      ! *** such as lagrangian multipliers
      scf_env%outer_scf%iter_count = 0
      accepted_ot_kp_searches = 0
      iter_count = 0
      ot_kp_subspace_refresh_count = 0
      ot_kp_subspace_refresh_iter_count = 0
      ot_kp_subspace_refresh_pending = .FALSE.
      total_steps = 0
      energy%tot_old = 0.0_dp

      scf_outer_loop: DO

         CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
                            scf_section=scf_section)

         oda_globalization_active = oda_eligible
         oda_predictor_lambda = 1.0_dp
         oda_predictor_valid = .FALSE.

         CALL qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)
         IF (ot_kp_subspace_refresh_pending) THEN
            scf_env%iter_count = ot_kp_subspace_refresh_iter_count
            ot_kp_subspace_refresh_pending = .FALSE.
         END IF

         ! decide whether to switch off dipole correction for convergence purposes
         dft_control%surf_dip_correct_switch = dft_control%correct_surf_dip
         IF ((dft_control%correct_surf_dip) .AND. (scf_control%outer_scf%have_scf) .AND. &
             (scf_env%outer_scf%iter_count > FLOOR(scf_control%outer_scf%max_scf/2.0_dp))) THEN
            IF (dft_control%switch_surf_dip) THEN
               dft_control%surf_dip_correct_switch = .FALSE.
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "SURFACE DIPOLE CORRECTION switched off"
            END IF
         END IF

         scf_loop: DO

            CALL timeset(routineN//"_inner_loop", handle2)

            IF (.NOT. just_energy) scf_env%iter_count = scf_env%iter_count + 1
            iter_count = iter_count + 1
            CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter_count)

            IF (output_unit > 0) CALL m_flush(output_unit)

            total_steps = total_steps + 1
            just_energy = energy_only

            CALL qs_ks_update_qs_env(qs_env, just_energy=just_energy, &
                                     calculate_forces=.FALSE.)

            scf_env%raw_map_delta = 0.0_dp
            scf_env%oda_lambda = 0.0_dp
            scf_env%oda_energy = 0.0_dp
            scf_env%oda_gradient0 = 0.0_dp
            scf_env%oda_gradient1 = 0.0_dp
            scf_env%oda_evaluations = 0
            scf_env%oda_status = ""
            scf_env%raw_map_delta_valid = .FALSE.
            oda_applied = .FALSE.
            oda_interior = .FALSE.
            oda_state_evaluated = .FALSE.
            ! print 'heavy weight' or relatively expensive quantities
            CALL qs_scf_loop_print(qs_env, scf_env, para_env)

            added_mos_auto_grow = .FALSE.
            ot_kp_subspace_refresh = .FALSE.
            adiis_validation = scf_control%diagonalization%update_method == diag_update_method_adiis .AND. &
                               scf_env%adiis_check_next
            adiis_step_delta = scf_env%step_norm
            scf_env%adiis_validated = .FALSE.
            oda_base_energy = energy

            IF (scf_control%diagonalization%update_method == diag_update_method_adiis) THEN
               ! ADIIS is a Fock-space SCF method and never uses density mixing.  A small
               ! accelerated step is checked against the unmodified F[P] map before
               ! convergence is accepted. The check detects an artificial subspace zero;
               ! it is not a second, equally strict convergence criterion.
               CPASSERT(ASSOCIATED(scf_env%scf_subspace_buffer))
               scf_env%scf_subspace_buffer%last_restart = .FALSE.
               scf_env%scf_subspace_buffer%use_combined_fock = .FALSE.
               scf_env%scf_subspace_buffer%last_old_fock_weight = 0.0_dp
               oda_trigger = .FALSE.
               oda_restart_trigger = .FALSE.
               IF (.NOT. adiis_validation .AND. scf_env%iter_count > 1) THEN
                  ! Only accepted densities whose raw F[P] has just been evaluated enter history.
                  ! The initial guess and raw-validation trial endpoint are deliberately excluded.
                  IF (scf_env%scf_subspace_buffer%diis_weight < 1.0_dp .OR. &
                      scf_env%iter_delta >= scf_control%eps_diis) THEN
                     CALL qs_scf_subspace_push(scf_env%scf_subspace_buffer, matrix_ks_kp, rho_ao_kp, &
                                               energy%total, adiis_pushed)
                     IF (.NOT. adiis_pushed) CPABORT("Failed to append the accepted ADIIS SCF state")
                     CALL qs_scf_subspace_build(scf_env%scf_subspace_buffer)
                  END IF
                  oda_restart_trigger = scf_env%scf_subspace_buffer%last_restart
                  ! ODA globalizes the initial raw Roothaan phase. Once ADIIS has formed
                  ! a genuine history combination, remain in ADIIS unless its history is
                  ! explicitly restarted; a later raw-vertex solution alone is not a new
                  ! globalization event.
                  IF (oda_restart_trigger) THEN
                     oda_globalization_active = .TRUE.
                  END IF
                  IF (scf_env%scf_subspace_buffer%use_combined_fock) THEN
                     oda_globalization_active = .FALSE.
                  END IF
                  IF (oda_restart_trigger .OR. scf_env%scf_subspace_buffer%use_combined_fock) THEN
                     oda_predictor_lambda = 1.0_dp
                     oda_predictor_valid = .FALSE.
                  END IF
                  oda_trigger = oda_eligible .AND. &
                                scf_env%iter_count < scf_control%max_scf .AND. &
                                oda_globalization_active .AND. &
                                .NOT. scf_env%scf_subspace_buffer%use_combined_fock
               END IF
               IF (do_kpoints) THEN
                  CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, &
                                         added_mos_auto_grow=added_mos_auto_grow)
               ELSE
                  CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only)
               END IF

               IF (oda_trigger .OR. adiis_validation) THEN
                  CALL qs_scf_candidate_density_delta(scf_env, rho, para_env, adiis_raw_delta)
                  scf_env%raw_map_delta = adiis_raw_delta
                  scf_env%raw_map_delta_valid = .TRUE.
               END IF

               IF (adiis_validation) THEN
                  oda_predictor_lambda = 1.0_dp
                  oda_predictor_valid = .FALSE.
                  IF (adiis_raw_delta < scf_control%eps_scf) THEN
                     scf_env%iter_delta = adiis_raw_delta
                     scf_env%iter_param = 0.0_dp
                     scf_env%iter_method = "ADIIS/Chk."
                     scf_env%adiis_validated = .TRUE.
                     scf_env%adiis_check_next = .FALSE.
                  ELSE
                     adiis_stagnated = adiis_step_delta <= adiis_stagnation_ratio*adiis_raw_delta
                     IF (adiis_stagnated) THEN
                        ! A nearly vanishing ADIIS step can be a fixed point of the current
                        ! subspace without being a fixed point of the raw SCF map.  Discard
                        ! the stale subspace, retain the current paired P,F[P] state, and
                        ! accept the already computed raw candidate as the restart step.
                        ! This avoids duplicate history entries and a second diagonalization.
                        CALL qs_scf_subspace_restart(scf_env%scf_subspace_buffer, matrix_ks_kp, rho_ao_kp, &
                                                     energy%total, adiis_restarted)
                        IF (.NOT. adiis_restarted) THEN
                           CPABORT("Failed to restart the ADIIS history from the current SCF state")
                        END IF
                        oda_restart_trigger = .TRUE.
                        oda_trigger = oda_eligible .AND. scf_env%iter_count < scf_control%max_scf
                        IF (oda_trigger) THEN
                           oda_globalization_active = .TRUE.
                        ELSE
                           scf_env%iter_delta = adiis_raw_delta
                           scf_env%iter_param = 0.0_dp
                           scf_env%iter_method = "ADIIS/Rst."
                           scf_env%adiis_check_next = .FALSE.
                        END IF
                     ELSE
                        ! The raw map is not converged, but the accelerated step is not an
                        ! artificial subspace zero. Preserve the useful history, discard the
                        ! raw trial, and repeat the normal ADIIS step in this iteration.
                        CALL qs_scf_subspace_push(scf_env%scf_subspace_buffer, matrix_ks_kp, rho_ao_kp, &
                                                  energy%total, adiis_pushed)
                        IF (.NOT. adiis_pushed) CPABORT("Failed to append the accepted ADIIS SCF state")
                        CALL qs_scf_subspace_build(scf_env%scf_subspace_buffer)
                        IF (scf_env%scf_subspace_buffer%last_restart) THEN
                           oda_globalization_active = .TRUE.
                        END IF
                        IF (scf_env%scf_subspace_buffer%use_combined_fock) THEN
                           oda_globalization_active = .FALSE.
                        END IF
                        IF (scf_env%scf_subspace_buffer%last_restart .OR. &
                            scf_env%scf_subspace_buffer%use_combined_fock) THEN
                           oda_predictor_lambda = 1.0_dp
                           oda_predictor_valid = .FALSE.
                        END IF
                        IF (do_kpoints) THEN
                           CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, &
                                                  added_mos_auto_grow=added_mos_auto_grow)
                        ELSE
                           CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only)
                        END IF
                        CALL qs_scf_candidate_density_delta(scf_env, rho, para_env, scf_env%iter_delta)
                        scf_env%adiis_validated = scf_env%iter_delta < scf_control%eps_scf
                        IF (scf_env%adiis_validated .AND. .NOT. diis_step) THEN
                           scf_env%iter_param = 0.0_dp
                           scf_env%iter_method = "ADIIS/Chk."
                        END IF
                        ! A successful sanity check must not recursively demand another raw
                        ! endpoint check. A later, independently small step can trigger one.
                        scf_env%adiis_check_next = .FALSE.
                     END IF
                  END IF
               END IF

               IF (oda_trigger) THEN
                  CPASSERT(ALLOCATED(scf_env%scf_subspace_buffer%generation))
                  CPASSERT(SIZE(scf_env%scf_subspace_buffer%generation) > 0)
                  CPASSERT(ASSOCIATED(scf_env%scf_subspace_buffer%fock))
                  CPASSERT(scf_env%scf_subspace_buffer%nstored > 0)
                  oda_base_slot = MAXLOC(scf_env%scf_subspace_buffer%generation, DIM=1)
                  CALL qs_scf_oda_apply(qs_env, scf_env, rho, ks_env, &
                                        scf_env%scf_subspace_buffer%fock(oda_base_slot, :, :), &
                                        matrix_ks_kp, rho_ao_kp, oda_base_energy, &
                                        oda_predictor_lambda, oda_predictor_valid, oda_applied, &
                                        oda_state_evaluated, oda_evaluated_energy)
                  IF (oda_applied) THEN
                     oda_interior = scf_env%oda_lambda < 1.0_dp - 64.0_dp*EPSILON(1.0_dp)
                  END IF
                  scf_env%iter_delta = scf_env%oda_lambda*adiis_raw_delta
                  scf_env%iter_param = scf_env%oda_lambda
                  IF (.NOT. oda_applied) THEN
                     scf_env%iter_method = "ODA/Fallback"
                  ELSE IF (oda_restart_trigger) THEN
                     scf_env%iter_method = "ODA/Rst."
                  ELSE
                     scf_env%iter_method = "ODA/Diag."
                  END IF
                  scf_env%adiis_check_next = scf_env%iter_delta < scf_control%eps_scf
               ELSE IF (.NOT. adiis_validation) THEN
                  CALL qs_scf_candidate_density_delta(scf_env, rho, para_env, scf_env%iter_delta)
                  scf_env%adiis_check_next = .NOT. diis_step .AND. &
                                             scf_env%iter_delta < scf_control%eps_scf
               END IF
            ELSE
               ! Keep the existing density-mixing SCF path unchanged.
               IF (do_kpoints) THEN
                  ! kpoints
                  IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                     scf_control%smear%do_smear = .FALSE.
                     CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, dft_control%probe, &
                                            energy_only=energy_only)
                  ELSE
                     CALL qs_scf_new_mos_kp( &
                        qs_env, scf_env, scf_control, diis_step, &
                        ot_kp_subspace_refresh=ot_kp_subspace_refresh, &
                        allow_ot_kp_subspace_refresh= &
                        ot_kp_subspace_refresh_count == 0, &
                        allow_ot_kp_exit_refresh= &
                        ot_kp_subspace_refresh_count > 0 .AND. &
                        ot_kp_subspace_refresh_count < max_ot_kp_subspace_refreshes, &
                        accepted_ot_kp_searches=accepted_ot_kp_searches, &
                        added_mos_auto_grow=added_mos_auto_grow, energy_only=energy_only)
                  END IF
               ELSE
                  ! Gamma points only
                  IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                     scf_control%smear%do_smear = .FALSE.
                     CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only, &
                                         dft_control%probe)
                  ELSE
                     CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only)
                  END IF
               END IF
            END IF

            IF (added_mos_auto_grow) THEN
               CALL qs_scf_grow_added_mos_auto_kp(qs_env, scf_env, scf_control, output_unit)
               total_steps = MAX(0, total_steps - 1)
               iter_count = MAX(0, iter_count - 1)
               IF (.NOT. just_energy) scf_env%iter_count = MAX(0, scf_env%iter_count - 1)
               CALL timestop(handle2)
               CYCLE scf_outer_loop
            END IF

            IF (ot_kp_subspace_refresh) THEN
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T2,A)") &
                     "K-point OT: rebuilding physical virtual subspace by full KS diagonalization."
               END IF
               IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
                  DO ispin = 1, SIZE(scf_env%qs_ot_env)
                     CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
                  END DO
                  DEALLOCATE (scf_env%qs_ot_env)
                  NULLIFY (scf_env%qs_ot_env)
               END IF
               IF (ASSOCIATED(kpoints%scf_diis_buffer)) THEN
                  CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
               END IF
               reuse_precond = .FALSE.
               accepted_ot_kp_searches = 0
               ot_kp_subspace_refresh_count = ot_kp_subspace_refresh_count + 1
               total_steps = MAX(0, total_steps - 1)
               iter_count = MAX(0, iter_count - 1)
               IF (.NOT. just_energy) scf_env%iter_count = MAX(0, scf_env%iter_count - 1)
               ot_kp_subspace_refresh_iter_count = scf_env%iter_count
               ot_kp_subspace_refresh_pending = .TRUE.
               CALL timestop(handle2)
               CYCLE scf_outer_loop
            END IF

            IF (do_kpoints .AND. scf_env%method == ot_method_nr .AND. &
                qs_scf_kp_search_endpoint(scf_env%iter_method)) THEN
               accepted_ot_kp_searches = accepted_ot_kp_searches + 1
            END IF

            ! Print requested MO information (can be computationally expensive with OT)
            ! The accepted interior ODA density is relaxed and has no matching MO set;
            ! suppress iteration-local MO properties until the next idempotent endpoint.
            IF (.NOT. oda_interior) THEN
               CALL qs_scf_write_mos(qs_env, scf_env, final_mos=.FALSE.)
            END IF

            IF (dft_control%qs_control%xtb_control%do_tblite) THEN
               IF (scf_env%method == ot_method_nr) THEN
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .TRUE.)
                  CALL evaluate_core_matrix_traces(qs_env)
               ELSE
                  CPASSERT(scf_env%mixing_method > 0)
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .FALSE.)
                  CALL evaluate_core_matrix_traces(qs_env, rho_ao_ext=scf_env%p_mix_new)
               END IF
               CALL tb_get_energy(qs_env, qs_env%tb_tblite, energy)
            END IF

            IF (scf_control%diagonalization%update_method == diag_update_method_adiis) THEN
               ! ADIIS never mixes densities. ODA trials already update rho; only an
               ! unevaluated raw endpoint still needs to be committed from p_mix_new.
               IF (.NOT. oda_state_evaluated) CALL qs_scf_commit_density_candidate(scf_env, rho)
               ! An accepted interior ODA density has no matching MO set.  Admit the
               ! state to CDIIS only after a subsequent full diagonalization step.
               scf_env%scf_subspace_buffer%diis_state_valid = .NOT. oda_interior
               density_full_step = .TRUE.
            ELSE
               density_full_step = diis_step .OR. tblite_native_mixer .OR. internal_tblite_density_full_step
               CALL qs_scf_density_mixing(scf_env, rho, para_env, density_full_step)
            END IF
            IF (dft_control%qs_control%xtb_control%do_tblite .AND. &
                .NOT. (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        tb_scf_mixer_error(dft_control, qs_env%tb_tblite, &
                                                           scf_control%eps_scf))
            END IF
            IF (dft_control%qs_control%dftb .OR. &
                (dft_control%qs_control%xtb .AND. .NOT. dft_control%qs_control%xtb_control%do_tblite)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        charge_mixing_scc_error(scf_env%mixing_store, scf_control%eps_scf))
            END IF
            IF (tblite_native_mixer) THEN
               scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               scf_env%iter_method = "TBLite/Diag"
            ELSE IF (internal_tblite_mixer) THEN
               scf_env%iter_method = "TBLite/Diag"
               IF (dft_control%qs_control%dftb) THEN
                  scf_env%iter_param = dft_control%qs_control%dftb_control%tblite_mixer_damping
               ELSE
                  scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               END IF
            END IF
            scf_env%step_norm = scf_env%iter_delta

            t2 = m_walltime()

            CALL qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy, &
                                  scf_control%diagonalization%adiis_verbose)

            IF (scf_control%gce%do_gce) THEN
               CALL qs_scf_gce_info(output_unit, qs_env, just_energy)
            END IF

            IF (.NOT. just_energy) energy%tot_old = energy%total

            ! check for external communicator and if the intermediate energy should be sent
            IF (scf_energy_message_tag > 0) THEN
               CALL external_comm%send(energy%total, ext_master_id, scf_energy_message_tag)
            END IF

            CALL qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, just_energy, &
                                         exit_inner_loop, inner_loop_converged, output_unit)

            ! In case we decide to exit we perform few more check to see if this one
            ! is really the last SCF step
            IF (exit_inner_loop) THEN

               CALL qs_scf_inner_finalize(scf_env, qs_env, density_full_step, output_unit)

               CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                            outer_loop_converged, exit_outer_loop)

               ! Let's tag the last SCF cycle so we can print informations only of the last step
               IF (exit_outer_loop) CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter_count)

            END IF

            IF (do_kpoints) THEN
               CALL write_kpoints_restart(rho_ao_kp, kpoints, scf_env, dft_section, particle_set, qs_kind_set)
            ELSE
               ! A damped ODA density is generally non-idempotent and cannot be represented by
               ! the raw-endpoint MOs. Defer the WFN restart until a later idempotent step.
               IF (.NOT. dft_control%mtlr_dft_with_perturbation .AND. &
                   .NOT. oda_interior) THEN
                  ! Write wavefunction restart file
                  IF (scf_env%method == ot_method_nr) THEN
                     ! With OT: provide the Kohn-Sham matrix for the calculation of the MO eigenvalues
                     CALL get_ks_env(ks_env=ks_env, matrix_ks=matrix_ks)
                     CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set, &
                                                  matrix_ks=matrix_ks)
                  ELSE
                     CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set)
                  END IF
               END IF
            END IF

            ! Exit if we have finished with the SCF inner loop
            IF (exit_inner_loop) THEN
               CALL timestop(handle2)
               EXIT scf_loop
            END IF

            IF (oda_state_evaluated) THEN
               ! ODA has already updated the accepted AO density, real-space density, Fock
               ! matrix, and all energy components.  Preserve the ordinary SCF row above,
               ! then carry that evaluated state directly into the next iteration.
               oda_evaluated_energy%tot_old = energy%tot_old
               energy = oda_evaluated_energy
            END IF

            IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       scf_section, "PRINT%ITERATION_INFO/TIME_CUMUL"), cp_p_file)) THEN
               t1 = m_walltime()
            END IF

            ! mixing methods have the new density matrix in p_mix_new
            IF (scf_env%mixing_method > 0) THEN
               DO ic = 1, SIZE(rho_ao_kp, 2)
                  DO ispin = 1, dft_control%nspins
                     CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                     CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
                  END DO
               END DO
            END IF

            IF (.NOT. oda_state_evaluated) THEN
               CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, &
                                      mix_rho=scf_env%mixing_method >= gspace_mixing_nr)
            END IF

            CALL timestop(handle2)

         END DO scf_loop

         IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

         ! In case we use the OUTER SCF loop let's print some info..
         CALL qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)

         ! Save MOs to converged MOs if outer_loop_converged and surf_dip_correct_switch is true
         IF (exit_outer_loop) THEN
            IF ((dft_control%switch_surf_dip) .AND. (outer_loop_converged) .AND. &
                (dft_control%surf_dip_correct_switch)) THEN
               DO ispin = 1, dft_control%nspins
                  CALL reassign_allocated_mos(mos_last_converged(ispin), mos(ispin))
               END DO
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "COPIED mos ---> mos_last_converged"
            END IF
         END IF

         IF (exit_outer_loop) EXIT scf_outer_loop

         !
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)

      END DO scf_outer_loop

      converged = inner_loop_converged .AND. outer_loop_converged
      total_scf_steps = total_steps

      IF (dft_control%qs_control%cdft) THEN
         dft_control%qs_control%cdft_control%total_steps = &
            dft_control%qs_control%cdft_control%total_steps + total_steps
      END IF

      IF (.NOT. converged) THEN
         IF (scf_control%ignore_convergence_failure .OR. should_stop) THEN
            CALL cp_warn(__LOCATION__, "SCF run NOT converged")
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "SCF run NOT converged. To continue the calculation "// &
                          "regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.")
         END IF
      END IF

      ! Skip Harris functional calculation if ground-state is NOT converged
      IF (qs_env%energy_correction) THEN
         CALL get_qs_env(qs_env, ec_env=ec_env)
         ec_env%do_skip = .FALSE.
         IF (ec_env%skip_ec .AND. .NOT. converged) ec_env%do_skip = .TRUE.
      END IF

      ! if needed copy mo_coeff dbcsr->fm for later use in post_scf!fm->dbcsr
      DO ispin = 1, SIZE(mos) !fm -> dbcsr
         IF (mos(ispin)%use_mo_coeff_b) THEN !fm->dbcsr
            IF (.NOT. ASSOCIATED(mos(ispin)%mo_coeff_b)) THEN
               !fm->dbcsr
               CPABORT("mo_coeff_b is not allocated")
            END IF !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, & !fm->dbcsr
                                  mos(ispin)%mo_coeff) !fm -> dbcsr
         END IF !fm->dbcsr
      END DO !fm -> dbcsr

      CALL cp_rm_iter_level(logger%iter_info, level_name="QS_SCF")
      CALL timestop(handle)

   END SUBROUTINE scf_env_do_scf

! **************************************************************************************************
!> \brief grow an automatic k-point smearing virtual space and rebuild dimensioned state
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_grow_added_mos_auto_kp(qs_env, scf_env, scf_control, output_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_grow_added_mos_auto_kp'

      INTEGER :: base_nmo, common_base_nmo, current_added, current_target_nmo, grow_by, handle, &
         ic, ik, ispin, nmo_mat, nspins, old_nmo, target_nmo
      INTEGER, DIMENSION(2)                              :: base_nmo_spin, new_added, new_nmo
      LOGICAL                                            :: diag_step, has_unit_metric, need_resize, &
                                                            shared_spin_auto
      REAL(KIND=dp)                                      :: energy_step, flexible_electron_count, &
                                                            maxocc, n_el_f
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, &
                                                            old_eigenvalues, old_occupation_numbers
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: old_mos
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, blacs_env, dft_control, eigenvalues, kpoints, matrix_ks, matrix_s, &
               mo_coeff, mo_derivs, mos, occupation_numbers, old_eigenvalues, &
               old_occupation_numbers, para_env)

      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      kpoints=kpoints, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      mo_derivs=mo_derivs, mos=mos, para_env=para_env, &
                      has_unit_metric=has_unit_metric)

      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(scf_control))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))

      nspins = dft_control%nspins
      CPASSERT(nspins >= 1 .AND. nspins <= SIZE(new_nmo))
      new_added = scf_control%added_mos
      new_nmo(:) = 0
      need_resize = .FALSE.
      shared_spin_auto = nspins == 2 .AND. ALL(scf_control%added_mos_auto(1:2))

      IF (shared_spin_auto) THEN
         DO ispin = 1, nspins
            current_added = MAX(0, scf_control%added_mos(ispin))
            base_nmo_spin(ispin) = MAX(0, mos(ispin)%nmo - current_added)
         END DO
         common_base_nmo = MAXVAL(base_nmo_spin(1:nspins))
         current_target_nmo = MAXVAL(mos(1:nspins)%nmo)
         current_added = MAX(0, current_target_nmo - common_base_nmo)
         grow_by = MAX(4, current_added)
         target_nmo = MIN(MINVAL(mos(1:nspins)%nao), current_target_nmo + grow_by)
         DO ispin = 1, nspins
            new_nmo(ispin) = target_nmo
            new_added(ispin) = target_nmo - base_nmo_spin(ispin)
            need_resize = need_resize .OR. new_nmo(ispin) > mos(ispin)%nmo
         END DO
      ELSE
         DO ispin = 1, nspins
            old_nmo = mos(ispin)%nmo
            new_nmo(ispin) = old_nmo
            IF (.NOT. scf_control%added_mos_auto(ispin)) CYCLE
            current_added = MAX(0, scf_control%added_mos(ispin))
            base_nmo = MAX(0, old_nmo - current_added)
            grow_by = MAX(4, current_added)
            new_added(ispin) = MIN(mos(ispin)%nao - base_nmo, current_added + grow_by)
            new_nmo(ispin) = base_nmo + new_added(ispin)
            need_resize = need_resize .OR. new_nmo(ispin) > old_nmo
         END DO

         IF (nspins == 2) THEN
            target_nmo = MAXVAL(new_nmo(1:nspins))
            DO ispin = 1, nspins
               IF (target_nmo > mos(ispin)%nao) THEN
                  CALL cp_abort(__LOCATION__, &
                                "K-point ADDED_MOS AUTO exhausted the AO basis while matching spin bands.")
               END IF
               IF (target_nmo > mos(ispin)%nmo .AND. .NOT. scf_control%added_mos_auto(ispin)) THEN
                  CALL cp_abort(__LOCATION__, &
                                "K-point ADDED_MOS AUTO needs to grow a spin channel with explicit ADDED_MOS.")
               END IF
               base_nmo = MAX(0, mos(ispin)%nmo - MAX(0, scf_control%added_mos(ispin)))
               new_nmo(ispin) = target_nmo
               new_added(ispin) = target_nmo - base_nmo
               need_resize = need_resize .OR. new_nmo(ispin) > mos(ispin)%nmo
            END DO
         END IF
      END IF

      IF (.NOT. need_resize) THEN
         CALL cp_abort(__LOCATION__, &
                       "K-point ADDED_MOS AUTO cannot grow further although the highest band is still occupied.")
      END IF

      scf_control%added_mos(1:nspins) = new_added(1:nspins)
      IF (output_unit > 0) THEN
         IF (nspins == 2) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,2I5)") &
               "K-point ADDED_MOS AUTO: growing virtual-space buffer to:", &
               scf_control%added_mos(1:nspins)
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "K-point ADDED_MOS AUTO: growing virtual-space buffer to: ", &
               scf_control%added_mos(1)
         END IF
      END IF

      ALLOCATE (old_mos(nspins))
      DO ispin = 1, nspins
         CALL duplicate_mo_set(old_mos(ispin), mos(ispin))
         CALL get_mo_set(old_mos(ispin), maxocc=maxocc, n_el_f=n_el_f, &
                         flexible_electron_count=flexible_electron_count)
         CALL deallocate_mo_set(mos(ispin))
         CALL allocate_mo_set(mo_set=mos(ispin), nao=old_mos(ispin)%nao, nmo=new_nmo(ispin), &
                              nelectron=old_mos(ispin)%nelectron, n_el_f=n_el_f, maxocc=maxocc, &
                              flexible_electron_count=flexible_electron_count)
         mos(ispin)%use_mo_coeff_b = old_mos(ispin)%use_mo_coeff_b
      END DO

      CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=mos, blacs_env=blacs_env, para_env=para_env)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)

      DO ispin = 1, nspins
         CALL init_mo_set(mos(ispin), fm_pool=ao_mo_fm_pools(ispin)%pool, &
                          name="qs_env%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                         occupation_numbers=occupation_numbers)
         CALL get_mo_set(old_mos(ispin), eigenvalues=old_eigenvalues, &
                         occupation_numbers=old_occupation_numbers)

         old_nmo = old_mos(ispin)%nmo
         CALL cp_fm_init_random(mo_coeff, mos(ispin)%nmo)
         IF (has_unit_metric) THEN
            CALL make_basis_simple(mo_coeff, mos(ispin)%nmo)
         ELSE
            CALL make_basis_sm(mo_coeff, mos(ispin)%nmo, matrix_s(1, 1)%matrix)
         END IF

         eigenvalues(1:old_nmo) = old_eigenvalues(1:old_nmo)
         occupation_numbers(:) = 0.0_dp
         occupation_numbers(1:old_nmo) = old_occupation_numbers(1:old_nmo)
         IF (mos(ispin)%nmo > old_nmo) THEN
            energy_step = MAX(scf_control%smear%electronic_temperature, 1.0E-3_dp)
            DO ic = old_nmo + 1, mos(ispin)%nmo
               eigenvalues(ic) = eigenvalues(old_nmo) + energy_step*REAL(ic - old_nmo, KIND=dp)
            END DO
         END IF
         mos(ispin)%homo = old_mos(ispin)%homo
         mos(ispin)%lfomo = old_mos(ispin)%lfomo
         mos(ispin)%kTS = old_mos(ispin)%kTS
         mos(ispin)%mu = old_mos(ispin)%mu
         mos(ispin)%uniform_occupation = old_mos(ispin)%uniform_occupation
      END DO

      IF (dft_control%restricted) CALL mo_set_restrict(mos)

      DO ispin = 1, nspins
         IF (.NOT. mos(ispin)%use_mo_coeff_b) CYCLE
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL dbcsr_init_p(mos(ispin)%mo_coeff_b)
         CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_coeff_b, template=matrix_s(1, 1)%matrix, &
                                                n=mos(ispin)%nmo, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(mo_coeff, mos(ispin)%mo_coeff_b)
      END DO

      IF (ASSOCIATED(mo_derivs)) THEN
         DO ispin = 1, SIZE(mo_derivs)
            IF (ASSOCIATED(mo_derivs(ispin)%matrix)) CALL dbcsr_release_p(mo_derivs(ispin)%matrix)
         END DO
         DEALLOCATE (mo_derivs)
         NULLIFY (mo_derivs)
         CALL set_qs_env(qs_env, mo_derivs=mo_derivs)
      END IF
      IF (qs_env%requires_mo_derivs) THEN
         nmo_mat = MERGE(1, nspins, dft_control%restricted)
         ALLOCATE (mo_derivs(nmo_mat))
         DO ispin = 1, nmo_mat
            NULLIFY (mo_derivs(ispin)%matrix)
            CALL dbcsr_init_p(mo_derivs(ispin)%matrix)
            CALL dbcsr_create(mo_derivs(ispin)%matrix, template=mos(ispin)%mo_coeff_b, &
                              name="mo_derivs", matrix_type=dbcsr_type_no_symmetry)
         END DO
         CALL set_qs_env(qs_env, mo_derivs=mo_derivs)
      END IF

      DO ik = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ik)%kpoint_env
         IF (ASSOCIATED(kp%mos)) THEN
            DO ispin = 1, SIZE(kp%mos, 2)
               DO ic = 1, SIZE(kp%mos, 1)
                  CALL deallocate_mo_set(kp%mos(ic, ispin))
               END DO
            END DO
            DEALLOCATE (kp%mos)
            NULLIFY (kp%mos)
         END IF
         CALL cp_fm_release(kp%pmat)
         CALL cp_fm_release(kp%wmat)
      END DO
      CALL mpools_release(kpoints%mpools)
      CALL kpoint_initialize_mos(kpoints, mos)
      CALL kpoint_initialize_mo_set(kpoints)

      diag_step = .FALSE.
      CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .FALSE., diag_step)

      IF (ASSOCIATED(kpoints%scf_diis_buffer)) CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
      IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
         NULLIFY (scf_env%qs_ot_env)
      END IF
      reuse_precond = .FALSE.

      DO ispin = 1, nspins
         CALL deallocate_mo_set(old_mos(ispin))
      END DO
      DEALLOCATE (old_mos)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_grow_added_mos_auto_kp

! **************************************************************************************************
!> \brief inits those objects needed if you want to restart the scf with, say
!>        only a new initial guess, or different density functional or ...
!>        this will happen just before the scf loop starts
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE init_scf_loop(scf_env, qs_env, scf_section)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_scf_loop'

      INTEGER                                            :: handle, ikind, ispin, nkpoint, nmo, &
                                                            nspin_ot
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: do_adiis, do_kpoints, do_rotation, &
                                                            has_unit_metric, kp_diis_step, &
                                                            kpoint_mos_initialized
      REAL(KIND=dp)                                      :: u_ramping
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h_kp, matrix_ks_kp, matrix_s_kp, &
                                                            matrix_t_kp
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (scf_control, matrix_h_kp, matrix_ks_kp, matrix_s, matrix_s_kp, matrix_t_kp, matrix_ks, &
               dft_control, mos, mo_coeff, kpoints, qs_kind_set, wkp)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env=qs_env, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      mos=mos, &
                      qs_kind_set=qs_kind_set)

      nkpoint = 1
      kp_range = 0
      nspin_ot = dft_control%nspins
      IF (dft_control%restricted) nspin_ot = 1
      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, nkp=nkpoint, kp_range=kp_range, wkp=wkp)
      END IF

      ! if using mo_coeff_b then copy to fm
      DO ispin = 1, SIZE(mos) !fm->dbcsr
         IF (mos(1)%use_mo_coeff_b) THEN !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff) !fm->dbcsr
         END IF !fm->dbcsr
      END DO !fm->dbcsr

      ! this just guarantees that all mo_occupations match the eigenvalues, if smear
      DO ispin = 1, dft_control%nspins
         ! do not reset mo_occupations if the maximum overlap method is in use
         IF (.NOT. scf_control%diagonalization%mom) THEN
            !if the hair probes section is present, this sends hairy_probes to set_mo_occupation subroutine
            !and switches off the standard smearing
            IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         probe=dft_control%probe)
               END IF
            ELSE
               IF (.NOT. scf_control%gce%do_gce) THEN
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         emit_warnings=.NOT. do_kpoints)
               ELSE
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         gce=scf_control%gce, &
                                         emit_warnings=.NOT. do_kpoints)
               END IF
            END IF
         END IF
      END DO

      do_adiis = scf_control%diagonalization%update_method == diag_update_method_adiis
      IF (do_adiis) THEN
         IF (scf_env%method /= general_diag_method_nr) THEN
            CPABORT("ADIIS currently requires SCF%DIAGONALIZATION ALGORITHM STANDARD")
         END IF
         IF (scf_env%cholesky_method == cholesky_dbcsr) THEN
            CPABORT("ADIIS is not yet compatible with CHOLESKY INVERSE_DBCSR")
         END IF
         IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb .OR. &
             dft_control%qs_control%semi_empirical) THEN
            CPABORT("ADIIS currently requires a Quickstep DFT KS matrix")
         END IF
         IF (dft_control%do_admm_mo) THEN
            CPABORT("ADIIS is not yet compatible with ADMM-MO")
         END IF
         IF (dft_control%sic_method_id == sic_eo) THEN
            CPABORT("ADIIS is not yet compatible with explicit-orbital SIC")
         END IF
         IF (dft_control%apply_period_efield) THEN
            CPABORT("ADIIS is not yet compatible with PERIODIC_EFIELD")
         END IF
         IF (dft_control%dft_plus_u) THEN
            CPASSERT(ASSOCIATED(qs_kind_set))
            DO ikind = 1, SIZE(qs_kind_set)
               CALL get_qs_kind(qs_kind_set(ikind), u_ramping=u_ramping)
               IF (u_ramping > 0.0_dp) THEN
                  CPABORT("ADIIS is not yet compatible with DFT+U U_RAMPING")
               END IF
            END DO
         END IF
         IF (dft_control%smear .OR. scf_control%smear%do_smear) THEN
            CPABORT("ADIIS finite-temperature smearing support is not implemented yet")
         END IF
         IF (dft_control%roks .OR. scf_control%diagonalization%mom .OR. &
             dft_control%hairy_probes .OR. scf_control%gce%do_gce) THEN
            CPABORT("ADIIS currently supports ordinary RKS/UKS occupations only")
         END IF
         IF (scf_control%do_diag_sub) THEN
            CPABORT("ADIIS is not yet compatible with DIAG_SUB_SCF")
         END IF
         IF (scf_control%diagonalization%max_history < 1) THEN
            CPABORT("ADIIS requires SCF%ADIIS%MAX_HISTORY >= 1")
         END IF
         IF (.NOT. ASSOCIATED(scf_env%scf_subspace_buffer)) THEN
            ALLOCATE (scf_env%scf_subspace_buffer)
            CALL qs_scf_subspace_buffer_create(scf_env%scf_subspace_buffer, &
                                               scf_control%diagonalization%max_history)
         ELSE IF (scf_env%scf_subspace_buffer%nbuffer /= scf_control%diagonalization%max_history) THEN
            CALL qs_scf_subspace_buffer_release(scf_env%scf_subspace_buffer)
            CALL qs_scf_subspace_buffer_create(scf_env%scf_subspace_buffer, &
                                               scf_control%diagonalization%max_history)
         END IF
         CALL qs_scf_subspace_buffer_clear(scf_env%scf_subspace_buffer)
      END IF

      SELECT CASE (scf_env%method)
      CASE DEFAULT

         CPABORT("Unknown SCF method <"//TRIM(cp_to_string(scf_env%method))//"> found. Check the code!")

      CASE (filter_matrix_diag_method_nr)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

      CASE (general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr, smeagol_method_nr)
         IF (.NOT. scf_env%skip_diis) THEN
            IF (do_kpoints) THEN
               IF (.NOT. ASSOCIATED(kpoints%scf_diis_buffer)) THEN
                  ALLOCATE (kpoints%scf_diis_buffer)
                  CALL qs_diis_b_create_kp(kpoints%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
            ELSE
               IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
                  ALLOCATE (scf_env%scf_diis_buffer)
                  CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
            END IF
         END IF

      CASE (ot_diag_method_nr)
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

         ! disable DFTB and SE for now
         IF (dft_control%qs_control%dftb .OR. &
             dft_control%qs_control%xtb .OR. &
             dft_control%qs_control%semi_empirical) THEN
            CPABORT("DFTB and SE not available with OT/DIAG")
         END IF

         ! if an old preconditioner is still around (i.e. outer SCF is active),
         ! remove it if this could be worthwhile
         CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     dft_control%nspins)

         CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     scf_control%diagonalization%ot_settings%precond_solver_type, &
                                     scf_control%diagonalization%ot_settings%energy_gap, dft_control%nspins, &
                                     chebyshev_degree=scf_control%diagonalization%ot_settings%chebyshev_degree)

      CASE (block_davidson_diag_method_nr)
         ! Preconditioner initialized within the loop, when required
      CASE (ot_method_nr)
         IF (do_kpoints) THEN
            IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN
               CALL get_qs_env(qs_env, matrix_h_kp=matrix_h_kp, matrix_s_kp=matrix_s_kp, &
                               kinetic_kp=matrix_t_kp)
               CPASSERT(ASSOCIATED(matrix_h_kp))
               CPASSERT(ASSOCIATED(matrix_s_kp))
               CALL allocate_qs_ot_envs(scf_env=scf_env, &
                                        scf_control=scf_control, &
                                        dft_control=dft_control, &
                                        scf_section=scf_section, &
                                        do_kpoints=do_kpoints, &
                                        nkpoint=nkpoint, &
                                        nspin_ot=nspin_ot, &
                                        kp_range=kp_range, &
                                        wkp=wkp)

               kpoint_mos_initialized = qs_ot_kpoint_mos_initialized(kpoints)
               IF (.NOT. kpoint_mos_initialized) THEN
                  kp_diis_step = .FALSE.
                  CALL do_general_diag_kp(matrix_h_kp, matrix_s_kp, kpoints, scf_env, &
                                          scf_control, .FALSE., kp_diis_step)
                  IF (SIZE(matrix_h_kp, 1) < nspin_ot) THEN
                     CPASSERT(SIZE(matrix_h_kp, 1) == 1)
                     CALL qs_ot_copy_kpoint_spin_mos(kpoints, nspin_ot)
                  END IF
               END IF
               CALL kpoint_set_mo_occupation(kpoints, scf_control%smear)
               IF (.NOT. (dft_control%smear .OR. scf_control%smear%do_smear)) THEN
                  CALL qs_ot_set_fixed_kpoint_occupations(kpoints, nspin_ot)
               END IF

               ! Build the initial KS matrix from the requested density guess.
               CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                                        calculate_forces=.FALSE.)
               CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks_kp)
               CPASSERT(ASSOCIATED(matrix_ks_kp))

               ! Variable-occupation OT needs its fixed REF complement from the
               ! physical Hamiltonian of the requested density. A restart needs the
               ! same synchronization for both fixed and variable occupations.
               IF (scf_env%qs_ot_env(1)%settings%do_ener .OR. &
                   scf_control%density_guess == restart_guess) THEN
                  kp_diis_step = .FALSE.
                  CALL do_general_diag_kp(matrix_ks_kp, matrix_s_kp, kpoints, scf_env, &
                                          scf_control, .FALSE., kp_diis_step)
                  IF (SIZE(matrix_ks_kp, 1) < nspin_ot) THEN
                     CPASSERT(SIZE(matrix_ks_kp, 1) == 1)
                     CALL qs_ot_copy_kpoint_spin_mos(kpoints, nspin_ot)
                  END IF
                  CALL kpoint_set_mo_occupation(kpoints, scf_control%smear)
                  IF (.NOT. (dft_control%smear .OR. scf_control%smear%do_smear)) THEN
                     CALL qs_ot_set_fixed_kpoint_occupations(kpoints, nspin_ot)
                  END IF
               END IF

               ! Redistribute the physical operators to their owning k-point groups.  OT must not
               ! call the collective real-space transform independently with group-local k values.
               kp_diis_step = .FALSE.
               CALL do_general_diag_kp(matrix_ks_kp, matrix_s_kp, kpoints, scf_env, &
                                       scf_control, .FALSE., kp_diis_step, matrix_t=matrix_t_kp, &
                                       store_ot_matrices=.TRUE., transfer_only=.TRUE.)

               CALL allocate_qs_ot_kpoint_state(qs_env=qs_env, &
                                                scf_env=scf_env, &
                                                dft_control=dft_control, &
                                                kpoints=kpoints, &
                                                matrix_ks_kp=matrix_ks_kp, &
                                                matrix_s_kp=matrix_s_kp, &
                                                matrix_t_kp=matrix_t_kp, &
                                                nspin_ot=nspin_ot)
            END IF

            CALL timestop(handle)
            RETURN
         END IF

         CALL get_qs_env(qs_env, &
                         has_unit_metric=has_unit_metric, &
                         matrix_s=matrix_s, &
                         matrix_ks=matrix_ks)

         ! reortho the wavefunctions if we are having an outer scf and
         ! this is not the first iteration
         ! this is useful to avoid the build-up of numerical noise
         ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
         IF (scf_control%do_outer_scf_reortho) THEN
            IF (scf_control%outer_scf%have_scf .AND. .NOT. dft_control%restricted) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  DO ispin = 1, dft_control%nspins
                     CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
                     IF (has_unit_metric) THEN
                        CALL make_basis_simple(mo_coeff, nmo)
                     ELSE
                        CALL make_basis_sm(mo_coeff, nmo, matrix_s(1)%matrix)
                     END IF
                  END DO
               END IF
            END IF
         ELSE
            ! dont need any dirty trick for the numerically stable irac algorithm.
         END IF

         IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN

            CALL allocate_qs_ot_envs(scf_env=scf_env, &
                                     scf_control=scf_control, &
                                     dft_control=dft_control, &
                                     scf_section=scf_section, &
                                     do_kpoints=do_kpoints, &
                                     nkpoint=nkpoint, &
                                     nspin_ot=nspin_ot, &
                                     kp_range=kp_range, &
                                     wkp=wkp)

            ! might need the KS matrix to init properly
            CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                                     calculate_forces=.FALSE.)

            ! if an old preconditioner is still around (i.e. outer SCF is active),
            ! remove it if this could be worthwhile
            IF (.NOT. reuse_precond) THEN
               CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           dft_control%nspins)
            END IF

            !
            ! preconditioning still needs to be done correctly with has_unit_metric
            ! notice that a big part of the preconditioning (S^-1) is fine anyhow
            !
            IF (has_unit_metric) THEN
               NULLIFY (orthogonality_metric)
            ELSE
               orthogonality_metric => matrix_s(1)%matrix
            END IF

            IF (.NOT. reuse_precond) THEN
               CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                           scf_env%qs_ot_env(1)%settings%energy_gap, dft_control%nspins, &
                                           has_unit_metric=has_unit_metric, &
                                           chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type, &
                                           chebyshev_degree=scf_env%qs_ot_env(1)%settings%chebyshev_degree)
            END IF
            IF (reuse_precond) reuse_precond = .FALSE.

            CALL ot_scf_init(mo_array=mos, matrix_s=orthogonality_metric, &
                             broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma, &
                             qs_ot_env=scf_env%qs_ot_env, matrix_ks=matrix_ks(1)%matrix)

            SELECT CASE (scf_env%qs_ot_env(1)%settings%preconditioner_type)
            CASE (ot_precond_none)
            CASE (ot_precond_full_all, ot_precond_full_single_inverse)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(ispin)%preconditioner)
               END DO
            CASE (ot_precond_s_inverse, ot_precond_full_single)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            CASE DEFAULT
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            END SELECT
         END IF

         ! if we have non-uniform occupations we should be using rotation
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         DO ispin = 1, SIZE(mos)
            IF (.NOT. mos(ispin)%uniform_occupation) THEN
               CPASSERT(do_rotation)
            END IF
         END DO
      END SELECT

      ! another safety check
      IF (dft_control%low_spin_roks) THEN
         CPASSERT(scf_env%method == ot_method_nr)
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         CPASSERT(do_rotation)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_scf_loop

! **************************************************************************************************
!> \brief allocate OT environments and label their spin/k-point channel identity
!> \param scf_env ...
!> \param scf_control ...
!> \param dft_control ...
!> \param scf_section ...
!> \param do_kpoints ...
!> \param nkpoint ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_envs(scf_env, scf_control, dft_control, scf_section, &
                                  do_kpoints, nkpoint, nspin_ot, kp_range, wkp)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: scf_section
      LOGICAL, INTENT(IN)                                :: do_kpoints
      INTEGER, INTENT(IN)                                :: nkpoint, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp

      CHARACTER(len=*), PARAMETER :: kpoint_precond_error = &
         "Complex K-point OT supports PRECONDITIONER NONE, FULL_S_INVERSE, "// &
         "FULL_KINETIC, FULL_SINGLE_INVERSE, or FULL_ALL.", kpoint_solver_error = &
         "Complex K-point FULL_ALL supports PRECOND_SOLVER DEFAULT; inverse "// &
         "preconditioners support DEFAULT or INVERSE_CHOLESKY."

      INTEGER                                            :: ikpoint, ispin, local_kpoint, &
                                                            number_of_OT_envs, ot_channel
      LOGICAL                                            :: do_rotation, do_smear, is_full_all
      REAL(KIND=dp)                                      :: kpoint_weight

      CPASSERT(.NOT. ASSOCIATED(scf_env%qs_ot_env))

      ! Restricted calculations require just one set of OT orbitals per k-point.
      number_of_OT_envs = qs_ot_number_of_channels(dft_control%nspins, &
                                                   nkpoint=nkpoint, &
                                                   restricted=dft_control%restricted)

      ALLOCATE (scf_env%qs_ot_env(number_of_OT_envs))

      ! XXX Joost XXX should disentangle reading input from this part
      IF (scf_env%outer_scf%iter_count > 0) THEN
         IF (scf_env%iter_delta < scf_control%eps_diis) THEN
            scf_env%qs_ot_env(1)%settings%ot_state = 1
         END IF
      END IF

      CALL ot_scf_read_input(scf_env%qs_ot_env, scf_section, do_kpoints)

      IF (do_kpoints) THEN
         do_smear = dft_control%smear .OR. scf_control%smear%do_smear
         IF (do_smear) THEN
            IF (.NOT. scf_env%qs_ot_env(1)%settings%do_ener) THEN
               CPABORT("K-point OT smearing requires OT%ENERGIES and OT%ROTATION.")
            END IF
            IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) THEN
               CPABORT("K-point OT smearing requires OT%ROTATION.")
            END IF
            SELECT CASE (scf_control%smear%method)
            CASE (smear_fermi_dirac, smear_gaussian, smear_mp, smear_mv)
               CONTINUE
            CASE DEFAULT
               CPABORT("K-point Mermin OT does not support the selected smearing method.")
            END SELECT
            IF (scf_env%qs_ot_env(1)%settings%ot_algorithm /= "TOD" .AND. &
                scf_env%qs_ot_env(1)%settings%ot_algorithm /= "REF") THEN
               CPABORT("K-point Mermin OT currently supports OT%ALGORITHM STRICT or IRAC.")
            END IF
            IF (scf_env%qs_ot_env(1)%settings%OT_METHOD /= "CG" .AND. &
                scf_env%qs_ot_env(1)%settings%OT_METHOD /= "SD" .AND. &
                scf_env%qs_ot_env(1)%settings%OT_METHOD /= "DIIS" .AND. &
                scf_env%qs_ot_env(1)%settings%OT_METHOD /= "BROY" .AND. &
                scf_env%qs_ot_env(1)%settings%OT_METHOD /= "LBFG") THEN
               CPABORT("K-point Mermin OT supports OT%MINIMIZER CG, SD, DIIS, BROYDEN, or LBFGS.")
            END IF
         ELSE IF (scf_env%qs_ot_env(1)%settings%do_ener) THEN
            CPABORT("OT%ENERGIES requires smearing in the complex K-point path.")
         END IF
         IF (scf_env%qs_ot_env(1)%settings%ot_algorithm /= "TOD" .AND. &
             scf_env%qs_ot_env(1)%settings%ot_algorithm /= "REF") THEN
            CPABORT("K-point OT currently supports OT%ALGORITHM STRICT or IRAC.")
         END IF
         IF (scf_env%qs_ot_env(1)%settings%OT_METHOD /= "CG" .AND. &
             scf_env%qs_ot_env(1)%settings%OT_METHOD /= "SD" .AND. &
             scf_env%qs_ot_env(1)%settings%OT_METHOD /= "DIIS" .AND. &
             scf_env%qs_ot_env(1)%settings%OT_METHOD /= "BROY" .AND. &
             scf_env%qs_ot_env(1)%settings%OT_METHOD /= "LBFG") THEN
            CPABORT("K-point OT currently supports OT%MINIMIZER CG, SD, DIIS, BROYDEN, or LBFGS.")
         END IF
         IF (.NOT. qs_ot_kpoint_preconditioner_supported( &
             scf_env%qs_ot_env(1)%settings%preconditioner_type, .FALSE.)) THEN
            CPABORT(kpoint_precond_error)
         END IF
         IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_kinetic .AND. &
             (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
              dft_control%qs_control%xtb)) THEN
            CPABORT("FULL_KINETIC is unavailable for semi-empirical methods")
         END IF
         IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
             scf_env%qs_ot_env(1)%settings%preconditioner_type, &
             scf_env%qs_ot_env(1)%settings%precond_solver_type, .FALSE.)) THEN
            CPABORT(kpoint_solver_error)
         END IF
         IF (scf_env%qs_ot_env(1)%settings%occupation_preconditioner .AND. .NOT. do_smear) THEN
            CPABORT("K-point OCCUPATION_PRECONDITIONER requires smearing.")
         END IF
         IF (dft_control%restricted) THEN
            CPABORT("Restricted open-shell K-point OT is not implemented.")
         END IF
      END IF

      IF (scf_env%outer_scf%iter_count > 0) THEN
         IF (scf_env%qs_ot_env(1)%settings%ot_state == 1) THEN
            scf_control%max_scf = MAX(scf_env%qs_ot_env(1)%settings%max_scf_diis, &
                                      scf_control%max_scf)
         END IF
      END IF

      ! Keep a note that we are restricted.
      IF (dft_control%restricted) THEN
         scf_env%qs_ot_env(:)%restricted = .TRUE.
         ! requires rotation
         IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) THEN
            CALL cp_abort(__LOCATION__, &
                          "Restricted calculation with OT requires orbital rotation. Please "// &
                          "activate the OT%ROTATION keyword!")
         END IF
      ELSE
         scf_env%qs_ot_env(:)%restricted = .FALSE.
      END IF

      ! This will rotate the MOs to be eigen states, which is not compatible with rotation.
      ! e.g. mo_derivs here do not yet include potentially different occupation numbers.
      do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
      ! Only full all needs rotation.
      is_full_all = scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all
      IF (do_rotation .AND. is_full_all .AND. .NOT. do_kpoints) THEN
         CPABORT('PRECONDITIONER FULL_ALL is not compatible with ROTATION.')
      END IF

      DO ikpoint = 1, nkpoint
         local_kpoint = 0
         IF (do_kpoints) THEN
            IF (ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2)) THEN
               local_kpoint = ikpoint - kp_range(1) + 1
            END IF
         END IF
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            IF (do_kpoints) THEN
               kpoint_weight = 1.0_dp
               IF (ASSOCIATED(wkp)) kpoint_weight = wkp(ikpoint)
               CALL qs_ot_set_context(scf_env%qs_ot_env(ot_channel), &
                                      spin_index=ispin, &
                                      kpoint_index=ikpoint, &
                                      local_kpoint_index=local_kpoint, &
                                      kpoint_weight=kpoint_weight)
            ELSE
               CALL qs_ot_set_context(scf_env%qs_ot_env(ot_channel), spin_index=ispin)
            END IF
         END DO
      END DO

      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, dft_control%nspins, &
                                       nkpoint=nkpoint, &
                                       restricted=dft_control%restricted, &
                                       require_kpoint=do_kpoints, &
                                       kp_range=kp_range, &
                                       wkp=wkp)

   END SUBROUTINE allocate_qs_ot_envs

! **************************************************************************************************
!> \brief checks whether every local K-point channel has nonzero orbitals
!> \param kpoints K-point environment
!> \return true if every channel has a nonzero real or imaginary MO coefficient
! **************************************************************************************************
   FUNCTION qs_ot_kpoint_mos_initialized(kpoints) RESULT(initialized)
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL                                            :: initialized

      INTEGER                                            :: icomponent, ikpoint_local, ispin, nocc
      REAL(KIND=dp)                                      :: channel_norm, coefficient_norm
      TYPE(cp_fm_struct_type), POINTER                   :: active_mo_struct
      TYPE(cp_fm_type)                                   :: active_mo_coeff
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (active_mo_struct, kp, mo_coeff, para_env)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))

      coefficient_norm = HUGE(1.0_dp)
      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, SIZE(kp%mos, 2)
            CALL get_mo_set(kp%mos(1, ispin), homo=nocc, mo_coeff=mo_coeff)
            CPASSERT(ASSOCIATED(mo_coeff))
            IF (nocc < 1) THEN
               initialized = .FALSE.
               RETURN
            END IF
            CALL cp_fm_struct_create(active_mo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                                     ncol_global=nocc)
            channel_norm = 0.0_dp
            DO icomponent = 1, SIZE(kp%mos, 1)
               CALL get_mo_set(kp%mos(icomponent, ispin), mo_coeff=mo_coeff)
               CPASSERT(ASSOCIATED(mo_coeff))
               CALL cp_fm_create(active_mo_coeff, active_mo_struct)
               CALL cp_fm_to_fm(mo_coeff, active_mo_coeff, nocc)
               channel_norm = MAX(channel_norm, cp_fm_norm(active_mo_coeff, "M"))
               CALL cp_fm_release(active_mo_coeff)
            END DO
            CALL cp_fm_struct_release(active_mo_struct)
            coefficient_norm = MIN(coefficient_norm, channel_norm)
         END DO
      END DO

      CALL get_kpoint_info(kpoints, para_env=para_env)
      CPASSERT(ASSOCIATED(para_env))
      CALL para_env%min(coefficient_norm)
      initialized = coefficient_norm > 100.0_dp*EPSILON(1.0_dp)

   END FUNCTION qs_ot_kpoint_mos_initialized

! **************************************************************************************************
!> \brief copy the spin-independent core-Hamiltonian eigenvectors into all UKS start channels
!> \param kpoints K-point environment
!> \param nspin number of spin channels
! **************************************************************************************************
   SUBROUTINE qs_ot_copy_kpoint_spin_mos(kpoints, nspin)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nspin

      INTEGER                                            :: icomponent, ikpoint_local, ispin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_source, eigenvalues_target
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_source, mo_coeff_target
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (eigenvalues_source, eigenvalues_target, kp, mo_coeff_source, mo_coeff_target)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(nspin >= 1)

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(SIZE(kp%mos, 2) >= nspin)
         DO icomponent = 1, SIZE(kp%mos, 1)
            CALL get_mo_set(kp%mos(icomponent, 1), mo_coeff=mo_coeff_source, &
                            eigenvalues=eigenvalues_source)
            CPASSERT(ASSOCIATED(mo_coeff_source))
            CPASSERT(ASSOCIATED(eigenvalues_source))
            DO ispin = 2, nspin
               CALL get_mo_set(kp%mos(icomponent, ispin), mo_coeff=mo_coeff_target, &
                               eigenvalues=eigenvalues_target)
               CPASSERT(ASSOCIATED(mo_coeff_target))
               CPASSERT(ASSOCIATED(eigenvalues_target))
               CALL cp_fm_to_fm(mo_coeff_source, mo_coeff_target)
               eigenvalues_target(:) = eigenvalues_source(:)
            END DO
         END DO
      END DO

   END SUBROUTINE qs_ot_copy_kpoint_spin_mos

! **************************************************************************************************
!> \brief restore the fixed occupied rank required by K-point OT after global Aufbau filling
!> \param kpoints K-point environment
!> \param nspin number of spin channels
! **************************************************************************************************
   SUBROUTINE qs_ot_set_fixed_kpoint_occupations(kpoints, nspin)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nspin

      INTEGER                                            :: icomponent, ikpoint_local, ispin, &
                                                            nelectron, nmo, nocc
      REAL(KIND=dp)                                      :: maxocc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (kp, occupation_numbers)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(nspin >= 1)

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(SIZE(kp%mos, 2) >= nspin)
         DO ispin = 1, nspin
            DO icomponent = 1, SIZE(kp%mos, 1)
               CALL get_mo_set(kp%mos(icomponent, ispin), homo=nocc, maxocc=maxocc, &
                               nelectron=nelectron, nmo=nmo, occupation_numbers=occupation_numbers)
               CPASSERT(ASSOCIATED(occupation_numbers))
               IF (nocc < 1 .OR. nocc > nmo) THEN
                  CPABORT("Fixed-occupation K-point OT requires a nonempty occupied subspace.")
               END IF
               IF (ABS(maxocc*REAL(nocc, KIND=dp) - REAL(nelectron, KIND=dp)) > &
                   100.0_dp*EPSILON(1.0_dp)) THEN
                  CPABORT("Fixed-occupation K-point OT requires uniform integer occupations.")
               END IF
               occupation_numbers(:) = 0.0_dp
               occupation_numbers(1:nocc) = maxocc
               kp%mos(icomponent, ispin)%kTS = 0.0_dp
            END DO
         END DO
      END DO

   END SUBROUTINE qs_ot_set_fixed_kpoint_occupations

! **************************************************************************************************
!> \brief allocate local k-point OT minimizer state for labelled spin/k-point channels
!> \param qs_env ...
!> \param scf_env ...
!> \param dft_control ...
!> \param kpoints ...
!> \param matrix_ks_kp ...
!> \param matrix_s_kp ...
!> \param matrix_t_kp ...
!> \param nspin_ot ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_kpoint_state(qs_env, scf_env, dft_control, &
                                          kpoints, matrix_ks_kp, matrix_s_kp, matrix_t_kp, nspin_ot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, matrix_t_kp
      INTEGER, INTENT(IN)                                :: nspin_ot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_qs_ot_kpoint_state'

      INTEGER                                            :: handle, ikpoint, ikpoint_local, ispin, &
                                                            nao, nmo, nocc, ot_channel
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_struct_type), POINTER                   :: active_mo_struct
      TYPE(cp_fm_type)                                   :: active_mo_coeff, active_mo_coeff_im
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_im
      TYPE(dbcsr_type), POINTER                          :: matrix_sk_im, matrix_sk_re
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_im

      CALL timeset(routineN, handle)

      NULLIFY (active_mo_struct, eigenvalues, kp, matrix_sk_im, matrix_sk_re, &
               mo_coeff, mo_coeff_im, mo_set, mo_set_im)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(ASSOCIATED(matrix_ks_kp))
      CPASSERT(ASSOCIATED(matrix_s_kp))
      CPASSERT(ASSOCIATED(matrix_s_kp(1, 1)%matrix))
      CPASSERT(.NOT. dft_control%restricted .OR. nspin_ot == 1)

      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)
      IF (use_real_wfn) THEN
         CPABORT("K-point OT currently requires complex K-point wavefunctions.")
      END IF

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(nspin_ot <= SIZE(kp%mos, 2))
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= 1)
         CPASSERT(SIZE(kp%ot_smat) >= 2)
         CALL qs_scf_get_ot_kpoint_operator(matrix_s_kp, kpoints, kp, 1, &
                                            kp%ot_smat(1), kp%ot_smat(2), &
                                            matrix_sk_re, matrix_sk_im)
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            CPASSERT(.NOT. scf_env%qs_ot_env(ot_channel)%state_allocated)

            NULLIFY (mo_coeff, mo_coeff_im, mo_set_im)
            mo_set => kp%mos(1, ispin)
            CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                            homo=nocc, nao=nao, nmo=nmo)
            CPASSERT(ASSOCIATED(mo_coeff))
            CPASSERT(ASSOCIATED(eigenvalues))
            CPASSERT(nao > 0)
            CPASSERT(nmo > 0)
            IF (scf_env%qs_ot_env(ot_channel)%settings%do_ener) nocc = nmo
            IF (nocc < 1 .OR. nocc > nmo) THEN
               CALL cp_abort(__LOCATION__, &
                             "K-point OT requires at least one occupied orbital in every spin channel.")
            END IF

            CALL cp_fm_struct_create(active_mo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                                     ncol_global=nocc)
            CALL cp_fm_create(active_mo_coeff, active_mo_struct)
            CALL cp_fm_to_fm(mo_coeff, active_mo_coeff, nocc)

            CALL qs_ot_allocate(scf_env%qs_ot_env(ot_channel), &
                                matrix_sk_re, &
                                active_mo_struct)

            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_c0, 0.0_dp)
            CALL copy_fm_to_dbcsr(active_mo_coeff, scf_env%qs_ot_env(ot_channel)%matrix_c0)

            IF (.NOT. use_real_wfn) THEN
               CPASSERT(SIZE(kp%mos, 1) >= 2)
               mo_set_im => kp%mos(2, ispin)
               CALL get_mo_set(mo_set=mo_set_im, mo_coeff=mo_coeff_im)
               CPASSERT(ASSOCIATED(mo_coeff_im))
               CALL cp_fm_create(active_mo_coeff_im, active_mo_struct)
               CALL cp_fm_to_fm(mo_coeff_im, active_mo_coeff_im, nocc)
               CALL allocate_qs_ot_kpoint_complex_state(scf_env%qs_ot_env(ot_channel), &
                                                        matrix_sk_re)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_c0_im, 0.0_dp)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sc0_im, 0.0_dp)
               CALL copy_fm_to_dbcsr(active_mo_coeff_im, scf_env%qs_ot_env(ot_channel)%matrix_c0_im)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_re, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_c0, 0.0_dp, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_sc0)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_im, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_c0_im, 0.0_dp, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_x)
               CALL dbcsr_add(scf_env%qs_ot_env(ot_channel)%matrix_sc0, &
                              scf_env%qs_ot_env(ot_channel)%matrix_x, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

               CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_re, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_c0_im, 0.0_dp, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_sc0_im)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sk_im, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_c0, 0.0_dp, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_x)
               CALL dbcsr_add(scf_env%qs_ot_env(ot_channel)%matrix_sc0_im, &
                              scf_env%qs_ot_env(ot_channel)%matrix_x, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END IF

            CALL qs_ot_init(scf_env%qs_ot_env(ot_channel))
            IF (scf_env%qs_ot_env(ot_channel)%settings%do_ener) THEN
               scf_env%qs_ot_env(ot_channel)%ener_x(:) = eigenvalues(1:nocc)
            END IF
            scf_env%qs_ot_env(ot_channel)%broyden_adaptive_sigma = qs_env%broyden_adaptive_sigma

            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_x, 0.0_dp)
            CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sx, 0.0_dp)
            IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_x_im, 0.0_dp)
               CALL dbcsr_set(scf_env%qs_ot_env(ot_channel)%matrix_sx_im, 0.0_dp)
            END IF

            SELECT CASE (scf_env%qs_ot_env(ot_channel)%settings%ot_algorithm)
            CASE ("TOD")
               IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
                  CALL qs_ot_get_p_complex(scf_env%qs_ot_env(ot_channel)%matrix_x, &
                                           scf_env%qs_ot_env(ot_channel)%matrix_x_im, &
                                           scf_env%qs_ot_env(ot_channel)%matrix_sx, &
                                           scf_env%qs_ot_env(ot_channel)%matrix_sx_im, &
                                           scf_env%qs_ot_env(ot_channel))
               ELSE
                  CALL qs_ot_get_p(scf_env%qs_ot_env(ot_channel)%matrix_x, &
                                   scf_env%qs_ot_env(ot_channel)%matrix_sx, &
                                   scf_env%qs_ot_env(ot_channel))
               END IF
            CASE ("REF")
               CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_x, &
                               scf_env%qs_ot_env(ot_channel)%matrix_c0)
               CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_sx, &
                               scf_env%qs_ot_env(ot_channel)%matrix_sc0)
               IF (scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_x_im, &
                                  scf_env%qs_ot_env(ot_channel)%matrix_c0_im)
                  CALL dbcsr_copy(scf_env%qs_ot_env(ot_channel)%matrix_sx_im, &
                                  scf_env%qs_ot_env(ot_channel)%matrix_sc0_im)
                  CALL qs_ot_get_orbitals_ref_complex(scf_env%qs_ot_env(ot_channel)%matrix_c0, &
                                                      scf_env%qs_ot_env(ot_channel)%matrix_c0_im, &
                                                      matrix_sk_re, matrix_sk_im, &
                                                      scf_env%qs_ot_env(ot_channel))
               END IF
            CASE DEFAULT
               CPABORT("Algorithm not yet implemented")
            END SELECT

            CPASSERT(scf_env%qs_ot_env(ot_channel)%state_allocated)
            IF (.NOT. use_real_wfn) THEN
               CPASSERT(scf_env%qs_ot_env(ot_channel)%has_complex_kpoint_state)
            END IF
            CALL cp_fm_release(active_mo_coeff)
            IF (.NOT. use_real_wfn) THEN
               CALL cp_fm_release(active_mo_coeff_im)
            END IF
            CALL cp_fm_struct_release(active_mo_struct)
         END DO
         CALL dbcsr_release_p(matrix_sk_re)
         CALL dbcsr_release_p(matrix_sk_im)
      END DO

      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, dft_control%nspins, &
                                       nkpoint=kpoints%nkp, &
                                       restricted=ANY(scf_env%qs_ot_env(:)%restricted), &
                                       require_kpoint=.TRUE., &
                                       kp_range=kpoints%kp_range, &
                                       require_local_state=.TRUE., &
                                       require_complex_state=.NOT. use_real_wfn)

      scf_env%qs_ot_env(:)%initial_density_sync_pending = .TRUE.

      CALL prepare_qs_ot_kpoint_preconditioners(scf_env, kpoints, matrix_ks_kp, &
                                                matrix_s_kp, matrix_t_kp, nspin_ot)

      CALL timestop(handle)

   END SUBROUTINE allocate_qs_ot_kpoint_state

! **************************************************************************************************
!> \brief prepare and attach orbital preconditioners for complex k-point OT channels
!> \param scf_env ...
!> \param kpoints ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_t ...
!> \param nspin_ot ...
! **************************************************************************************************
   SUBROUTINE prepare_qs_ot_kpoint_preconditioners(scf_env, kpoints, matrix_ks, matrix_s, matrix_t, nspin_ot)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_t
      INTEGER, INTENT(IN)                                :: nspin_ot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_qs_ot_kpoint_preconditioners'

      INTEGER                                            :: handle, ikpoint, ispin, ks_spin, &
                                                            local_kpoint, n_ot_channels, &
                                                            ot_channel, prec_type, source_channel
      TYPE(dbcsr_type), POINTER                          :: matrix_h_im, matrix_h_re, matrix_s_im, &
                                                            matrix_s_re, matrix_t_im, matrix_t_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL timeset(routineN, handle)

      NULLIFY (kp, matrix_h_im, matrix_h_re, matrix_s_im, matrix_s_re, matrix_t_im, matrix_t_re)
      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))

      prec_type = scf_env%qs_ot_env(1)%settings%preconditioner_type
      IF (prec_type == ot_precond_none) THEN
         CALL timestop(handle)
         RETURN
      END IF
      SELECT CASE (prec_type)
      CASE (ot_precond_full_all, ot_precond_full_single_inverse, &
            ot_precond_full_kinetic, ot_precond_s_inverse)
      CASE DEFAULT
         CPABORT("Unsupported complex K-point OT preconditioner")
      END SELECT
      n_ot_channels = SIZE(scf_env%qs_ot_env)
      IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
         DO ot_channel = 1, SIZE(scf_env%ot_preconditioner)
            IF (ASSOCIATED(scf_env%ot_preconditioner(ot_channel)%preconditioner)) THEN
               CALL destroy_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner)
               DEALLOCATE (scf_env%ot_preconditioner(ot_channel)%preconditioner)
            END IF
         END DO
         DEALLOCATE (scf_env%ot_preconditioner)
         NULLIFY (scf_env%ot_preconditioner)
      END IF
      ALLOCATE (scf_env%ot_preconditioner(n_ot_channels))
      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            CPASSERT(scf_env%qs_ot_env(ot_channel)%state_allocated)
            ALLOCATE (scf_env%ot_preconditioner(ot_channel)%preconditioner)
            CALL init_preconditioner(scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                                     para_env=scf_env%qs_ot_env(ot_channel)%para_env, &
                                     blacs_env=scf_env%qs_ot_env(ot_channel)%blacs_env)
         END DO
      END DO

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(SIZE(kp%ot_smat) >= 2)
         ikpoint = kp%nkpoint
         CALL qs_scf_get_ot_kpoint_operator(matrix_s, kpoints, kp, 1, &
                                            kp%ot_smat(1), kp%ot_smat(2), &
                                            matrix_s_re, matrix_s_im)
         IF (prec_type == ot_precond_full_kinetic) THEN
            CPASSERT(ASSOCIATED(matrix_t))
            CPASSERT(ASSOCIATED(kp%ot_tmat))
            CPASSERT(SIZE(kp%ot_tmat) >= 2)
            CALL qs_scf_get_ot_kpoint_operator(matrix_t, kpoints, kp, 1, &
                                               kp%ot_tmat(1), kp%ot_tmat(2), &
                                               matrix_t_re, matrix_t_im)
         END IF
         DO ispin = 1, nspin_ot
            ot_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot)
            IF (ispin > 1 .AND. &
                (prec_type == ot_precond_full_kinetic .OR. prec_type == ot_precond_s_inverse)) THEN
               source_channel = qs_ot_channel_index(1, ikpoint, nspin_ot)
               CPASSERT(ASSOCIATED(scf_env%ot_preconditioner(source_channel)%preconditioner%complex_fm))
               scf_env%ot_preconditioner(ot_channel)%preconditioner%complex_fm => &
                  scf_env%ot_preconditioner(source_channel)%preconditioner%complex_fm
               scf_env%ot_preconditioner(ot_channel)%preconditioner%owns_complex_fm = .FALSE.
               scf_env%ot_preconditioner(ot_channel)%preconditioner%energy_gap = &
                  scf_env%ot_preconditioner(source_channel)%preconditioner%energy_gap
               scf_env%ot_preconditioner(ot_channel)%preconditioner%in_use = &
                  scf_env%ot_preconditioner(source_channel)%preconditioner%in_use
               scf_env%ot_preconditioner(ot_channel)%preconditioner%solver = &
                  scf_env%ot_preconditioner(source_channel)%preconditioner%solver
               CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ot_channel), &
                                             scf_env%ot_preconditioner(ot_channel)%preconditioner)
               CYCLE
            END IF
            IF (prec_type == ot_precond_full_all .OR. &
                prec_type == ot_precond_full_single_inverse) THEN
               ks_spin = MIN(ispin, SIZE(kp%ot_hmat, 2))
               CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
               CALL qs_scf_get_ot_kpoint_operator(matrix_ks, kpoints, kp, ks_spin, &
                                                  kp%ot_hmat(1, ks_spin), kp%ot_hmat(2, ks_spin), &
                                                  matrix_h_re, matrix_h_im)
            END IF
            SELECT CASE (prec_type)
            CASE (ot_precond_full_all)
               CALL make_preconditioner_complex_full_all( &
                  scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                  scf_env%qs_ot_env(ot_channel)%matrix_c0, &
                  scf_env%qs_ot_env(ot_channel)%matrix_c0_im, &
                  matrix_h_re, matrix_h_im, matrix_s_re, matrix_s_im, &
                  kp%mos(1, ispin), &
                  scf_env%qs_ot_env(ot_channel)%settings%energy_gap, &
                  scf_env%qs_ot_env(ot_channel)%settings%precond_solver_type)
            CASE (ot_precond_full_single_inverse)
               CALL make_preconditioner_complex_full_single_inverse( &
                  scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                  scf_env%qs_ot_env(ot_channel)%matrix_c0, &
                  scf_env%qs_ot_env(ot_channel)%matrix_c0_im, &
                  matrix_h_re, matrix_h_im, matrix_s_re, matrix_s_im, &
                  scf_env%qs_ot_env(ot_channel)%settings%energy_gap, &
                  scf_env%qs_ot_env(ot_channel)%settings%precond_solver_type)
            CASE (ot_precond_full_kinetic)
               CALL make_preconditioner_complex_full_kinetic( &
                  scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                  matrix_t_re, matrix_t_im, matrix_s_re, matrix_s_im, &
                  scf_env%qs_ot_env(ot_channel)%settings%energy_gap, &
                  scf_env%qs_ot_env(ot_channel)%settings%precond_solver_type)
            CASE (ot_precond_s_inverse)
               CALL make_preconditioner_complex_full_s_inverse( &
                  scf_env%ot_preconditioner(ot_channel)%preconditioner, &
                  matrix_s_re, matrix_s_im, &
                  scf_env%qs_ot_env(ot_channel)%settings%precond_solver_type)
            END SELECT
            CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ot_channel), &
                                          scf_env%ot_preconditioner(ot_channel)%preconditioner)
            IF (ASSOCIATED(matrix_h_re)) CALL dbcsr_release_p(matrix_h_re)
            IF (ASSOCIATED(matrix_h_im)) CALL dbcsr_release_p(matrix_h_im)
         END DO
         CALL dbcsr_release_p(matrix_s_re)
         CALL dbcsr_release_p(matrix_s_im)
         IF (ASSOCIATED(matrix_t_re)) CALL dbcsr_release_p(matrix_t_re)
         IF (ASSOCIATED(matrix_t_im)) CALL dbcsr_release_p(matrix_t_im)
      END DO

      CALL timestop(handle)

   END SUBROUTINE prepare_qs_ot_kpoint_preconditioners

! **************************************************************************************************
!> \brief allocate imaginary orbital images for complex k-point OT channel state
!> \param qs_ot_env ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE allocate_qs_ot_kpoint_complex_state(qs_ot_env, matrix_s)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      INTEGER                                            :: i, ncoef, nmo, nmo_ortho

      CPASSERT(qs_ot_env%state_allocated)
      CPASSERT(.NOT. qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_s))

      CALL dbcsr_get_info(qs_ot_env%matrix_c0, nfullcols_total=nmo)
      CALL dbcsr_get_info(qs_ot_env%matrix_sc0, nfullcols_total=nmo_ortho)

      CALL dbcsr_init_p(qs_ot_env%matrix_c0_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_c0_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_sc0_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sc0_im, template=matrix_s, n=nmo_ortho, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_x_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_x_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_sx_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_sx_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(qs_ot_env%matrix_gx_im)
      CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx_im, template=matrix_s, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)
      IF (qs_ot_env%settings%occupation_preconditioner) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_preconditioned_gx_im)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_preconditioned_gx_im, &
                                                template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
         IF (qs_ot_env%use_dx) THEN
            CALL dbcsr_init_p(qs_ot_env%matrix_response_gx)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_response_gx, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
            CALL dbcsr_init_p(qs_ot_env%matrix_response_gx_im)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_response_gx_im, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
            CALL dbcsr_init_p(qs_ot_env%matrix_mermin_g0)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_mermin_g0, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
            CALL dbcsr_init_p(qs_ot_env%matrix_mermin_g0_im)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_mermin_g0_im, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
            IF (qs_ot_env%settings%do_ener) ALLOCATE (qs_ot_env%ener_mermin_g0(nmo))
         END IF
      END IF
      SELECT CASE (qs_ot_env%settings%ot_algorithm)
      CASE ("TOD")
         CPASSERT(nmo_ortho == nmo)
         CALL allocate_complex_copy(qs_ot_env%matrix_p_im, qs_ot_env%matrix_p, "matrix_p_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_r_im, qs_ot_env%matrix_r, "matrix_r_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_sinp_im, qs_ot_env%matrix_sinp, "matrix_sinp_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_cosp_im, qs_ot_env%matrix_cosp, "matrix_cosp_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_buf1, "matrix_buf1_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf2_im, qs_ot_env%matrix_buf2, "matrix_buf2_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf3_im, qs_ot_env%matrix_buf3, "matrix_buf3_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf4_im, qs_ot_env%matrix_buf4, "matrix_buf4_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_os_im, qs_ot_env%matrix_os, "matrix_os_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf1_ortho_im, qs_ot_env%matrix_buf1_ortho, &
                                    "matrix_buf1_ortho_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf2_ortho_im, qs_ot_env%matrix_buf2_ortho, &
                                    "matrix_buf2_ortho_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_tmp_ortho, qs_ot_env%matrix_buf1_ortho, &
                                    "matrix_tmp_ortho")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_x, "matrix_buf_nk")
         CALL allocate_complex_copy(qs_ot_env%matrix_buf_nk_im, qs_ot_env%matrix_x, "matrix_buf_nk_im")
         CALL allocate_complex_copy(qs_ot_env%matrix_tmp_nk, qs_ot_env%matrix_x, "matrix_tmp_nk")
      CASE ("REF")
         CALL dbcsr_init_p(qs_ot_env%matrix_ref_inv_sqrt)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_ref_inv_sqrt, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%matrix_ref_inv_sqrt_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%matrix_ref_inv_sqrt_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%buf1_n_k_dp)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%buf1_n_k_dp, template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
      CASE DEFAULT
         CPABORT("Complex K-point OT state requires ALGORITHM STRICT or IRAC")
      END SELECT
      IF (qs_ot_env%use_dx) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_dx_im)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_dx_im, template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
      END IF
      IF (qs_ot_env%use_gx_old) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_gx_old_im)
         CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_gx_old_im, template=matrix_s, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)
      END IF
      IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
          qs_ot_env%settings%ot_method == "BROY" .OR. &
          qs_ot_env%settings%ot_method == "LBFG") THEN
         ncoef = qs_ot_env%settings%diis_m
         IF (qs_ot_env%settings%ot_method == "LBFG") ncoef = ncoef + 1
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_e_im, ncoef)
         CALL dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_x_im, ncoef)
         DO i = 1, ncoef
            CALL dbcsr_init_p(qs_ot_env%matrix_h_e_im(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_e_im(i)%matrix, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
            CALL dbcsr_init_p(qs_ot_env%matrix_h_x_im(i)%matrix)
            CALL cp_dbcsr_m_by_n_from_row_template(qs_ot_env%matrix_h_x_im(i)%matrix, &
                                                   template=matrix_s, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
         END DO
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_init_p(qs_ot_env%rot_mat_u_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_u_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_x_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_x_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_dedu_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dedu_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_chc_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_chc_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         CALL dbcsr_init_p(qs_ot_env%rot_mat_gx_im)
         CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx_im, template=matrix_s, m=nmo, n=nmo, &
                                            sym=dbcsr_type_no_symmetry)
         IF (qs_ot_env%settings%occupation_preconditioner) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_response_gx_im)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_response_gx_im, &
                                               template=matrix_s, m=nmo, n=nmo, &
                                               sym=dbcsr_type_no_symmetry)
            IF (qs_ot_env%use_dx) THEN
               CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_g0)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_g0, &
                                                  template=matrix_s, m=nmo, n=nmo, &
                                                  sym=dbcsr_type_no_symmetry)
               CALL dbcsr_init_p(qs_ot_env%rot_mat_mermin_g0_im)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_mermin_g0_im, &
                                                  template=matrix_s, m=nmo, n=nmo, &
                                                  sym=dbcsr_type_no_symmetry)
            END IF
         END IF
         IF (qs_ot_env%use_dx) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_dx_im)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_dx_im, template=matrix_s, m=nmo, n=nmo, &
                                               sym=dbcsr_type_no_symmetry)
         END IF
         IF (qs_ot_env%use_gx_old) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_gx_old_im)
            CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_gx_old_im, template=matrix_s, m=nmo, n=nmo, &
                                               sym=dbcsr_type_no_symmetry)
         END IF
         IF (qs_ot_env%settings%ot_method == "DIIS" .OR. &
             qs_ot_env%settings%ot_method == "BROY" .OR. &
             qs_ot_env%settings%ot_method == "LBFG") THEN
            ncoef = qs_ot_env%settings%diis_m
            IF (qs_ot_env%settings%ot_method == "LBFG") ncoef = ncoef + 1
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_e_im, ncoef)
            CALL dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_x_im, ncoef)
            DO i = 1, ncoef
               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_e_im(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_e_im(i)%matrix, &
                                                  template=matrix_s, m=nmo, n=nmo, &
                                                  sym=dbcsr_type_no_symmetry)
               CALL dbcsr_init_p(qs_ot_env%rot_mat_h_x_im(i)%matrix)
               CALL cp_dbcsr_m_by_n_from_template(qs_ot_env%rot_mat_h_x_im(i)%matrix, &
                                                  template=matrix_s, m=nmo, n=nmo, &
                                                  sym=dbcsr_type_no_symmetry)
            END DO
         END IF
      END IF

      qs_ot_env%has_complex_kpoint_state = .TRUE.

   CONTAINS

! **************************************************************************************************
!> \brief allocate a zeroed matrix with the distribution and sparsity of a real companion
!> \param matrix output matrix
!> \param template real companion matrix
!> \param name DBCSR matrix name
! **************************************************************************************************
      SUBROUTINE allocate_complex_copy(matrix, template, name)
      TYPE(dbcsr_type), POINTER                          :: matrix, template
      CHARACTER(LEN=*), INTENT(IN)                       :: name

         CALL dbcsr_init_p(matrix)
         CALL dbcsr_copy(matrix, template, name=name)
         CALL dbcsr_set(matrix, 0.0_dp)
      END SUBROUTINE allocate_complex_copy

   END SUBROUTINE allocate_qs_ot_kpoint_complex_state

! **************************************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param scf_env ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_cleanup(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scf_env_cleanup'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Release SCF work storage
      CALL cp_fm_release(scf_env%scf_work1)

      IF (ASSOCIATED(scf_env%scf_work1_red)) THEN
         CALL cp_fm_release(scf_env%scf_work1_red)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2)) THEN
         CALL cp_fm_release(scf_env%scf_work2)
         DEALLOCATE (scf_env%scf_work2)
         NULLIFY (scf_env%scf_work2)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2_red)) THEN
         CALL cp_fm_release(scf_env%scf_work2_red)
         DEALLOCATE (scf_env%scf_work2_red)
         NULLIFY (scf_env%scf_work2_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho)) THEN
         CALL cp_fm_release(scf_env%ortho)
         DEALLOCATE (scf_env%ortho)
         NULLIFY (scf_env%ortho)
      END IF
      IF (ASSOCIATED(scf_env%ortho_red)) THEN
         CALL cp_fm_release(scf_env%ortho_red)
         DEALLOCATE (scf_env%ortho_red)
         NULLIFY (scf_env%ortho_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1)) THEN
         CALL cp_fm_release(scf_env%ortho_m1)
         DEALLOCATE (scf_env%ortho_m1)
         NULLIFY (scf_env%ortho_m1)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1_red)) THEN
         CALL cp_fm_release(scf_env%ortho_m1_red)
         DEALLOCATE (scf_env%ortho_m1_red)
         NULLIFY (scf_env%ortho_m1_red)
      END IF

      IF (ASSOCIATED(scf_env%ortho_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%ortho_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf1_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf1_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf2_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf2_dbcsr)
      END IF

      IF (ASSOCIATED(scf_env%p_mix_new)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_mix_new)
      END IF

      IF (ASSOCIATED(scf_env%p_delta)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_delta)
      END IF

      ! Method dependent cleanup
      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         !
      CASE (ot_diag_method_nr)
         !
      CASE (general_diag_method_nr)
         !
      CASE (special_diag_method_nr)
         !
      CASE (block_krylov_diag_method_nr)
      CASE (block_davidson_diag_method_nr)
         CALL block_davidson_deallocate(scf_env%block_davidson_env)
      CASE (filter_matrix_diag_method_nr)
         !
      CASE (smeagol_method_nr)
         !
      CASE DEFAULT
         CPABORT("unknown scf method method:"//cp_to_string(scf_env%method))
      END SELECT

      IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
         DEALLOCATE (scf_env%outer_scf%variables)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
         DEALLOCATE (scf_env%outer_scf%count)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
         DEALLOCATE (scf_env%outer_scf%gradient)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
         DEALLOCATE (scf_env%outer_scf%energy)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. &
          scf_env%outer_scf%deallocate_jacobian) THEN
         DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      END IF

      CALL timestop(handle)

   END SUBROUTINE scf_env_cleanup

! **************************************************************************************************
!> \brief perform a CDFT scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param should_stop flag determining if calculation should stop
!> \par History
!>      12.2015 Created
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_scf(qs_env, should_stop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cdft_scf'

      INTEGER                                            :: handle, iatom, iimage, ispin, ivar, nmo, &
                                                            nvar, output_unit, tsteps
      LOGICAL                                            :: cdft_loop_converged, converged, &
                                                            exit_cdft_loop, first_iteration, &
                                                            my_uocc, uniform_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gapw_cdft_values
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_occupations
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gapw_wmat, matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(local_rho_type), POINTER                      :: gapw_operator_rho
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (atomic_kind_set, gapw_operator_rho, gapw_wmat, para_env, qs_kind_set, &
               scf_env, ks_env, energy, rho, matrix_s, rho_ao, cdft_control, logger, &
               dft_control, pw_env, auxbas_pw_pool, energy, ks_env, scf_env, dft_section, &
               input, scf_section, scf_control, mos, mo_occupations)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, energy=energy, &
                      dft_control=dft_control, scf_control=scf_control, &
                      ks_env=ks_env, input=input)

      CALL timeset(routineN//"_loop", handle)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      first_iteration = .TRUE.

      cdft_control => dft_control%qs_control%cdft_control

      scf_env%outer_scf%iter_count = 0
      cdft_control%total_steps = 0

      ! Write some info about the CDFT calculation
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "CDFT EXTERNAL SCF WAVEFUNCTION OPTIMIZATION"
         CALL qs_scf_cdft_initial_info(output_unit, cdft_control)
      END IF
      IF (cdft_control%reuse_precond) THEN
         reuse_precond = .FALSE.
         cdft_control%nreused = 0
      END IF
      cdft_outer_loop: DO
         ! Change outer_scf settings to OT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
         ! Solve electronic structure with fixed value of constraint
         CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                             converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ! Decide whether to reuse the preconditioner on the next iteration
         IF (cdft_control%reuse_precond) THEN
            ! For convergence in exactly one step, the preconditioner is always reused (assuming max_reuse > 0)
            ! usually this means that the electronic structure has already converged to the correct state
            ! but the constraint optimizer keeps jumping over the optimal solution
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count == 1 &
                .AND. cdft_control%total_steps /= 1) THEN
               cdft_control%nreused = cdft_control%nreused - 1
            END IF
            ! SCF converged in less than precond_freq steps
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count <= cdft_control%precond_freq .AND. &
                cdft_control%total_steps /= 1 .AND. cdft_control%nreused < cdft_control%max_reuse) THEN
               reuse_precond = .TRUE.
               cdft_control%nreused = cdft_control%nreused + 1
            ELSE
               reuse_precond = .FALSE.
               cdft_control%nreused = 0
            END IF
         END IF
         ! Update history purging counters
         IF (first_iteration .AND. cdft_control%purge_history) THEN
            cdft_control%istep = cdft_control%istep + 1
            IF (scf_env%outer_scf%iter_count > 1) THEN
               cdft_control%nbad_conv = cdft_control%nbad_conv + 1
               IF (cdft_control%nbad_conv >= cdft_control%purge_freq .AND. &
                   cdft_control%istep >= cdft_control%purge_offset) THEN
                  cdft_control%nbad_conv = 0
                  cdft_control%istep = 0
                  cdft_control%should_purge = .TRUE.
               END IF
            END IF
         END IF
         first_iteration = .FALSE.
         ! Change outer_scf settings to CDFT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
         CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                      cdft_loop_converged, exit_cdft_loop)
         CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                               energy, cdft_control%total_steps, &
                               should_stop, cdft_loop_converged, cdft_loop=.TRUE.)
         IF (exit_cdft_loop) EXIT cdft_outer_loop
         ! Check if the inverse Jacobian needs to be calculated
         CALL qs_calculate_inverse_jacobian(qs_env)
         ! Check if a line search should be performed to find an optimal step size for the optimizer
         CALL qs_cdft_line_search(qs_env)
         ! Optimize constraint
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
      END DO cdft_outer_loop

      cdft_control%ienergy = cdft_control%ienergy + 1

      ! Store needed arrays for ET coupling calculation
      IF (cdft_control%do_et) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos)
         nvar = SIZE(cdft_control%target)
         IF (dft_control%qs_control%gapw) THEN
            IF (dft_control%nimages /= 1) THEN
               CALL cp_abort(__LOCATION__, &
                             "GAPW CDFT-CI currently requires a Gamma-point calculation.")
            END IF
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, para_env=para_env, &
                            qs_kind_set=qs_kind_set, rho=rho)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            CALL local_rho_set_create(gapw_operator_rho)
            CALL allocate_rho_atom_internals(gapw_operator_rho%rho_atom_set, atomic_kind_set, &
                                             qs_kind_set, dft_control, para_env)
            ALLOCATE (gapw_cdft_values(nvar), gapw_wmat(dft_control%nspins*dft_control%nimages))
            DO iimage = 1, dft_control%nimages
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_init_p(gapw_wmat(dft_control%nspins*(iimage - 1) + ispin)%matrix)
                  CALL dbcsr_copy(gapw_wmat(dft_control%nspins*(iimage - 1) + ispin)%matrix, &
                                  matrix_s(iimage)%matrix, name="GAPW CDFT WEIGHT MATRIX")
               END DO
            END DO
         END IF
         ! Matrix representation of weight function
         ALLOCATE (cdft_control%wmat(nvar))
         DO ivar = 1, nvar
            CALL dbcsr_init_p(cdft_control%wmat(ivar)%matrix)
            CALL dbcsr_copy(cdft_control%wmat(ivar)%matrix, matrix_s(1)%matrix, &
                            name="ET_RESTRAINT_MATRIX")
            CALL dbcsr_set(cdft_control%wmat(ivar)%matrix, 0.0_dp)
            CALL integrate_v_rspace(cdft_control%group(ivar)%weight, &
                                    hmat=cdft_control%wmat(ivar), qs_env=qs_env, &
                                    calculate_forces=.FALSE., &
                                    gapw=dft_control%qs_control%gapw)
            IF (dft_control%qs_control%gapw) THEN
               DO ispin = 1, SIZE(gapw_wmat)
                  CALL dbcsr_set(gapw_wmat(ispin)%matrix, 0.0_dp)
               END DO
               CALL zero_rho_atom_integrals(gapw_operator_rho%rho_atom_set)
               CALL gapw_cdft_one_center(qs_env, energy_only=.FALSE., calculate_forces=.FALSE., &
                                         values=gapw_cdft_values, operator_group=ivar, &
                                         rho_atom_operator_set=gapw_operator_rho%rho_atom_set)
               CALL update_ks_atom(qs_env, gapw_wmat, rho_ao, forces=.FALSE., &
                                   rho_atom_external=gapw_operator_rho%rho_atom_set)
               CALL dbcsr_add(cdft_control%wmat(ivar)%matrix, gapw_wmat(1)%matrix, 1.0_dp, 1.0_dp)
            END IF
         END DO
         IF (dft_control%qs_control%gapw) THEN
            CALL dbcsr_deallocate_matrix_set(gapw_wmat)
            CALL local_rho_set_release(gapw_operator_rho)
            DEALLOCATE (gapw_cdft_values)
         END IF
         ! Overlap matrix
         CALL dbcsr_init_p(cdft_control%matrix_s%matrix)
         CALL dbcsr_copy(cdft_control%matrix_s%matrix, matrix_s(1)%matrix, &
                         name="OVERLAP")
         ! Molecular orbital coefficients
         NULLIFY (cdft_control%mo_coeff)
         ALLOCATE (cdft_control%mo_coeff(dft_control%nspins))
         DO ispin = 1, dft_control%nspins
            CALL cp_fm_create(matrix=cdft_control%mo_coeff(ispin), &
                              matrix_struct=qs_env%mos(ispin)%mo_coeff%matrix_struct, &
                              name="MO_COEFF_A"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
            CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_coeff, &
                             cdft_control%mo_coeff(ispin))
         END DO
         ! Density matrix
         IF (cdft_control%calculate_metric) THEN
            CALL get_qs_env(qs_env, rho=rho)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            ALLOCATE (cdft_control%matrix_p(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               NULLIFY (cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_init_p(cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_copy(cdft_control%matrix_p(ispin)%matrix, rho_ao(ispin)%matrix, &
                               name="DENSITY MATRIX")
            END DO
         END IF
         ! Copy occupation numbers if non-uniform occupation
         uniform_occupation = .TRUE.
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin), uniform_occupation=my_uocc)
            uniform_occupation = uniform_occupation .AND. my_uocc
         END DO
         IF (.NOT. uniform_occupation) THEN
            ALLOCATE (cdft_control%occupations(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin), &
                               nmo=nmo, &
                               occupation_numbers=mo_occupations)
               ALLOCATE (cdft_control%occupations(ispin)%array(nmo))
               cdft_control%occupations(ispin)%array(1:nmo) = mo_occupations(1:nmo)
            END DO
         END IF
      END IF

      ! Deallocate constraint storage if forces are not needed
      ! In case of a simulation with multiple force_evals,
      ! deallocate only if weight function should not be copied to different force_evals
      IF (.NOT. (cdft_control%save_pot .OR. cdft_control%transfer_pot)) THEN
         CALL get_qs_env(qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         DO iatom = 1, SIZE(cdft_control%group)
            CALL auxbas_pw_pool%give_back_pw(cdft_control%group(iatom)%weight)
            DEALLOCATE (cdft_control%group(iatom)%weight)
         END DO
         IF (cdft_control%atomic_charges) THEN
            DO iatom = 1, cdft_control%natoms
               CALL auxbas_pw_pool%give_back_pw(cdft_control%charge(iatom))
            END DO
            DEALLOCATE (cdft_control%charge)
         END IF
         IF (cdft_control%type == outer_scf_becke_constraint .AND. &
             cdft_control%becke_control%cavity_confine) THEN
            IF (.NOT. ASSOCIATED(cdft_control%becke_control%cavity_mat)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%becke_control%cavity)
            ELSE
               DEALLOCATE (cdft_control%becke_control%cavity_mat)
            END IF
         ELSE IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
            IF (ASSOCIATED(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)
            END IF
         END IF
         IF (ASSOCIATED(cdft_control%charges_fragment)) DEALLOCATE (cdft_control%charges_fragment)
         cdft_control%need_pot = .TRUE.
         cdft_control%external_control = .FALSE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE cdft_scf

! **************************************************************************************************
!> \brief perform cleanup operations for cdft_control
!> \param cdft_control container for the external CDFT SCF loop variables
!> \par History
!>      12.2015 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_control_cleanup(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      IF (ASSOCIATED(cdft_control%constraint%variables)) THEN
         DEALLOCATE (cdft_control%constraint%variables)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%count)) THEN
         DEALLOCATE (cdft_control%constraint%count)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%gradient)) THEN
         DEALLOCATE (cdft_control%constraint%gradient)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%energy)) THEN
         DEALLOCATE (cdft_control%constraint%energy)
      END IF
      IF (ASSOCIATED(cdft_control%constraint%inv_jacobian) .AND. &
          cdft_control%constraint%deallocate_jacobian) THEN
         DEALLOCATE (cdft_control%constraint%inv_jacobian)
      END IF

   END SUBROUTINE cdft_control_cleanup

! **************************************************************************************************
!> \brief Calculates the finite difference inverse Jacobian
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \par History
!>      01.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_calculate_inverse_jacobian(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_calculate_inverse_jacobian'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: counter, handle, i, ispin, iter_count, &
                                                            iwork, j, max_scf, nspins, nsteps, &
                                                            nvar, nwork, output_unit, pwork, &
                                                            tsteps, twork
      LOGICAL                                            :: converged, explicit_jacobian, &
                                                            should_build, should_stop, &
                                                            use_md_history
      REAL(KIND=dp)                                      :: inv_error, step_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, dh, step_multiplier
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: mos_stashed
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, cdft_control, &
               inv_jacobian, para_env, tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      explicit_jacobian = .FALSE.
      should_build = .FALSE.
      use_md_history = .FALSE.
      iter_count = scf_env%outer_scf%iter_count
      ! Quick exit if optimizer does not require Jacobian
      IF (.NOT. ASSOCIATED(scf_control%outer_scf%cdft_opt_control)) RETURN
      ! Check if Jacobian should be calculated and initialize
      CALL timeset(routineN, handle)
      CALL initialize_inverse_jacobian(scf_control, scf_env, explicit_jacobian, should_build, used_history)
      IF (scf_control%outer_scf%cdft_opt_control%jacobian_restart) THEN
         ! Restart from previously calculated inverse Jacobian
         should_build = .FALSE.
         CALL restart_inverse_jacobian(qs_env)
      END IF
      IF (should_build) THEN
         scf_env%outer_scf%deallocate_jacobian = .FALSE.
         ! Actually need to (re)build the Jacobian
         IF (explicit_jacobian) THEN
            ! Build Jacobian with finite differences
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that need the explicit Jacobian can"// &
                             " only be used together with a valid CDFT constraint.")
            END IF
            ! Redirect output from Jacobian calculation to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-JacobianInfo.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            ALLOCATE (jacobian(nvar, nvar))
            jacobian = 0.0_dp
            nsteps = cdft_control%total_steps
            ! Setup finite difference scheme
            CALL prepare_jacobian_stencil(qs_env, output_unit, nwork, pwork, coeff, step_multiplier, dh)
            twork = pwork - nwork
            DO i = 1, nvar
               jacobian(i, :) = coeff(0)*scf_env%outer_scf%gradient(i, iter_count)
            END DO
            ! Calculate the Jacobian by perturbing each Lagrangian and recalculating the energy self-consistently
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, nvar
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I3,A,I3,A)') &
                     " ###  Constraint        ", i, " of ", nvar, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               counter = 0
               DO iwork = nwork, pwork
                  IF (iwork == 0) CYCLE
                  counter = counter + 1
                  IF (output_unit > 0) THEN
                     WRITE (output_unit, FMT="(A)") " #####################################"
                     WRITE (output_unit, '(A,I3,A,I3,A)') &
                        " ###  Energy evaluation ", counter, " of ", twork, " ###"
                     WRITE (output_unit, FMT="(A)") " #####################################"
                  END IF
                  IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) == 1) THEN
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(1)
                  ELSE
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(i)
                  END IF
                  scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count)
                  scf_env%outer_scf%variables(i, iter_count + 1) = scf_env%outer_scf%variables(i, iter_count) + &
                                                                   step_multiplier(iwork)*step_size
                  CALL outer_loop_update_qs_env(qs_env, scf_env)
                  CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
                  CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                      converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
                  ! Update (iter_count + 1) element of gradient and print constraint info
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
                  CALL outer_loop_gradient(qs_env, scf_env)
                  CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                        energy_qs, cdft_control%total_steps, &
                                        should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
                  ! Update Jacobian
                  DO j = 1, nvar
                     jacobian(j, i) = jacobian(j, i) + coeff(iwork)*scf_env%outer_scf%gradient(j, iter_count + 1)
                  END DO
                  ! Reset everything to last converged state
                  scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
                  scf_env%outer_scf%gradient = gradient
                  scf_env%outer_scf%energy = energy
                  cdft_control%total_steps = nsteps
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
               END DO
            END DO
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Finalize and invert Jacobian
            DO j = 1, nvar
               DO i = 1, nvar
                  jacobian(i, j) = jacobian(i, j)/dh(j)
               END DO
            END DO
            IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
            END IF
            inv_jacobian => scf_env%outer_scf%inv_jacobian
            CALL invert_matrix(jacobian, inv_jacobian, inv_error)
            scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, jacobian, gradient, energy, coeff, step_multiplier, dh)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== JACOBIAN CALCULATED =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         ELSE
            ! Build a strictly diagonal Jacobian from history and invert it
            CALL build_diagonal_jacobian(qs_env, used_history)
         END IF
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. para_env%is_source()) THEN
         ! Write restart file for inverse Jacobian
         CALL print_inverse_jacobian(logger, scf_env%outer_scf%inv_jacobian, iter_count)
      END IF
      ! Update counter
      scf_control%outer_scf%cdft_opt_control%ijacobian(1) = scf_control%outer_scf%cdft_opt_control%ijacobian(1) + 1
      CALL timestop(handle)

   END SUBROUTINE qs_calculate_inverse_jacobian

! **************************************************************************************************
!> \brief Perform backtracking line search to find the optimal step size for the CDFT constraint
!>        optimizer. Assumes that the CDFT gradient function is a smooth function of the constraint
!>        variables.
!> \param qs_env the qs_environment_type where to perform the line search
!> \par History
!>      02.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_cdft_line_search(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_cdft_line_search'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: handle, i, ispin, iter_count, &
                                                            max_linesearch, max_scf, nspins, &
                                                            nsteps, nvar, output_unit, tsteps
      LOGICAL :: continue_ls, continue_ls_exit, converged, do_linesearch, found_solution, &
         reached_maxls, should_exit, should_stop, sign_changed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: positive_sign
      REAL(KIND=dp)                                      :: alpha, alpha_ls, factor, norm_ls
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      REAL(KIND=dp), EXTERNAL                            :: dnrm2
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, &
               cdft_control, inv_jacobian, para_env, &
               tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      do_linesearch = .FALSE.
      SELECT CASE (scf_control%outer_scf%optimizer)
      CASE DEFAULT
         do_linesearch = .FALSE.
      CASE (outer_scf_optimizer_newton_ls)
         do_linesearch = .TRUE.
      CASE (outer_scf_optimizer_broyden)
         SELECT CASE (scf_control%outer_scf%cdft_opt_control%broyden_type)
         CASE (broyden_type_1, broyden_type_2, broyden_type_1_explicit, broyden_type_2_explicit)
            do_linesearch = .FALSE.
         CASE (broyden_type_1_ls, broyden_type_1_explicit_ls, broyden_type_2_ls, broyden_type_2_explicit_ls)
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            END IF
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               do_linesearch = .TRUE.
            END IF
         END SELECT
      END SELECT
      IF (do_linesearch) THEN
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), ALLOCATABLE :: mos_ls, mos_stashed
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            END IF
            CPASSERT(ASSOCIATED(scf_env%outer_scf%inv_jacobian))
            CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control))
            alpha = scf_control%outer_scf%cdft_opt_control%newton_step_save
            iter_count = scf_env%outer_scf%iter_count
            ! Redirect output from line search procedure to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-LineSearch.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            nsteps = cdft_control%total_steps
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            max_linesearch = scf_control%outer_scf%cdft_opt_control%max_ls
            continue_ls = scf_control%outer_scf%cdft_opt_control%continue_ls
            factor = scf_control%outer_scf%cdft_opt_control%factor_ls
            continue_ls_exit = .FALSE.
            found_solution = .FALSE.
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            reached_maxls = .FALSE.
            ! Broyden optimizers: perform update of inv_jacobian if necessary
            IF (scf_control%outer_scf%cdft_opt_control%broyden_update) THEN
               CALL outer_loop_optimize(scf_env, scf_control)
               ! Reset the variables and prevent a reupdate of inv_jacobian
               scf_env%outer_scf%variables(:, iter_count + 1) = 0
               scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            END IF
            ! Print some info
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH STARTED  =================================="
               WRITE (output_unit, FMT="(A,I5,A)") &
                  " Evaluating optimal step size for optimizer using a maximum of", max_linesearch, " steps"
               IF (continue_ls) THEN
                  WRITE (output_unit, FMT="(A)") &
                     " Line search continues until best step size is found or max steps are reached"
               END IF
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Initial step size: ", alpha
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Step size update factor: ", factor
               WRITE (output_unit, '(/,A,I10,A,I10)') &
                  " Energy evaluation: ", cdft_control%ienergy, ", CDFT SCF iteration: ", iter_count
            END IF
            ! Perform backtracking line search
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, max_linesearch
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I10,A)') &
                     " ###  Line search step: ", i, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               ! Newton update of CDFT variables with a step size of alpha
               scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count) - alpha* &
                                                                MATMUL(inv_jacobian, scf_env%outer_scf%gradient(:, iter_count))
               ! With updated CDFT variables, perform SCF
               CALL outer_loop_update_qs_env(qs_env, scf_env)
               CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
               CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                   converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
               ! Update (iter_count + 1) element of gradient and print constraint info
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
               CALL outer_loop_gradient(qs_env, scf_env)
               CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                     energy_qs, cdft_control%total_steps, &
                                     should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
               ! Store sign of initial gradient for each variable for continue_ls
               IF (continue_ls .AND. .NOT. ALLOCATED(positive_sign)) THEN
                  ALLOCATE (positive_sign(nvar))
                  DO ispin = 1, nvar
                     positive_sign(ispin) = scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp
                  END DO
               END IF
               ! Check if the L2 norm of the gradient decreased
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) < &
                   dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count), 1)) THEN
                  ! Optimal step size found
                  IF (.NOT. continue_ls) THEN
                     should_exit = .TRUE.
                  ELSE
                     ! But line search continues for at least one more iteration in an attempt to find a better solution
                     ! if max number of steps is not exceeded
                     IF (found_solution) THEN
                        ! Check if the norm also decreased w.r.t. to previously found solution
                        IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) > norm_ls) THEN
                           ! Norm increased => accept previous solution and exit
                           continue_ls_exit = .TRUE.
                        END IF
                     END IF
                     ! Store current state and the value of alpha
                     IF (.NOT. continue_ls_exit) THEN
                        should_exit = .FALSE.
                        alpha_ls = alpha
                        found_solution = .TRUE.
                        norm_ls = dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1)
                        ! Check if the sign of the gradient has changed for all variables (w.r.t initial gradient)
                        ! In this case we should exit because further line search steps will just increase the norm
                        sign_changed = .TRUE.
                        DO ispin = 1, nvar
                           sign_changed = sign_changed .AND. (positive_sign(ispin) .NEQV. &
                                                              scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp)
                        END DO
                        IF (.NOT. ALLOCATED(mos_ls)) THEN
                           ALLOCATE (mos_ls(nspins))
                        ELSE
                           DO ispin = 1, nspins
                              CALL deallocate_mo_set(mos_ls(ispin))
                           END DO
                        END IF
                        DO ispin = 1, nspins
                           CALL duplicate_mo_set(mos_ls(ispin), mos(ispin))
                        END DO
                        alpha = alpha*factor
                        ! Exit on last iteration
                        IF (i == max_linesearch) continue_ls_exit = .TRUE.
                        ! Exit if constraint target is satisfied to requested tolerance
                        IF (SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count + 1)**2)) < &
                            scf_control%outer_scf%eps_scf) THEN
                           continue_ls_exit = .TRUE.
                        END IF
                        ! Exit if line search jumped over the optimal step length
                        IF (sign_changed) continue_ls_exit = .TRUE.
                     END IF
                  END IF
               ELSE
                  ! Gradient increased => alpha is too large (if the gradient function is smooth)
                  should_exit = .FALSE.
                  ! Update alpha using Armijo's scheme
                  alpha = alpha*factor
               END IF
               IF (continue_ls_exit) THEN
                  ! Continuation of line search did not yield a better alpha, use previously located solution and exit
                  alpha = alpha_ls
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_ls(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                     CALL deallocate_mo_set(mos_ls(ispin))
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
                  DEALLOCATE (mos_ls)
                  should_exit = .TRUE.
               END IF
               ! Reached max steps and SCF converged: continue with last iterated step size
               IF (.NOT. should_exit .AND. &
                   (i == max_linesearch .AND. converged .AND. .NOT. found_solution)) THEN
                  should_exit = .TRUE.
                  reached_maxls = .TRUE.
                  alpha = alpha*(1.0_dp/factor)
               END IF
               ! Reset outer SCF environment to last converged state
               scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
               scf_env%outer_scf%gradient = gradient
               scf_env%outer_scf%energy = energy
               ! Exit line search if a suitable step size was found
               IF (should_exit) EXIT
               ! Reset the electronic structure
               cdft_control%total_steps = nsteps
               DO ispin = 1, nspins
                  CALL deallocate_mo_set(mos(ispin))
                  CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                  CALL calculate_density_matrix(mos(ispin), &
                                                p_rmpv(ispin)%matrix)
               END DO
               CALL qs_rho_update_rho(rho, qs_env=qs_env)
               CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
            END DO
            scf_control%outer_scf%cdft_opt_control%newton_step = alpha
            IF (.NOT. should_exit) THEN
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with fixed step size.")
               scf_control%outer_scf%cdft_opt_control%newton_step = scf_control%outer_scf%cdft_opt_control%newton_step_save
            END IF
            IF (reached_maxls) THEN
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with lasted iterated step size.")
            END IF
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, gradient, energy)
            IF (ALLOCATED(positive_sign)) DEALLOCATE (positive_sign)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH COMPLETE =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         END BLOCK
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_cdft_line_search

END MODULE qs_scf

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
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
   USE cp_control_types,                ONLY: cdft_control_type,&
                                              dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_constants,                 ONLY: &
        cdft2ot, cdft_combined_acceptor, cdft_combined_all, cdft_combined_constraint, &
        cdft_combined_donor, history_guess, jacobian_fd1, jacobian_fd1_backward, &
        jacobian_fd1_central, jacobian_fd2, jacobian_fd2_backward, ot2cdft, ot_precond_full_all, &
        ot_precond_full_single, ot_precond_full_single_inverse, ot_precond_none, &
        ot_precond_s_inverse, outer_scf_hirshfeld_constraint, outer_scf_optimizer_broyden, &
        outer_scf_optimizer_newton
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_io,                       ONLY: write_kpoints_restart
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: invert_matrix
   USE message_passing,                 ONLY: mp_send
   USE particle_types,                  ONLY: particle_type
   USE preconditioner,                  ONLY: prepare_preconditioner,&
                                              restart_preconditioner
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_density_mixing_types,         ONLY: gspace_mixing_nr
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_io,                        ONLY: write_mo_set
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type
   USE qs_ot,                           ONLY: qs_ot_new_preconditioner
   USE qs_ot_scf,                       ONLY: ot_scf_init,&
                                              ot_scf_read_input
   USE qs_outer_scf,                    ONLY: outer_loop_gradient,&
                                              outer_loop_optimize,&
                                              outer_loop_purge_history,&
                                              outer_loop_switch,&
                                              outer_loop_update_qs_env
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
   USE qs_scf_loop_utils,               ONLY: qs_scf_check_inner_exit,&
                                              qs_scf_check_outer_exit,&
                                              qs_scf_density_mixing,&
                                              qs_scf_inner_finalize,&
                                              qs_scf_new_mos,&
                                              qs_scf_new_mos_kp,&
                                              qs_scf_rho_update,&
                                              qs_scf_set_loop_flags
   USE qs_scf_output,                   ONLY: qs_scf_cdft_info,&
                                              qs_scf_cdft_initial_info,&
                                              qs_scf_loop_info,&
                                              qs_scf_loop_print,&
                                              qs_scf_outer_loop_info,&
                                              qs_scf_write_mos
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_scf_env_type, scf_env_release, &
        special_diag_method_nr
   USE qs_wf_history_methods,           ONLY: wfi_purge_history,&
                                              wfi_update
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'
   LOGICAL, PRIVATE                     :: reuse_precond = .FALSE.

   PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf, cdft_scf

CONTAINS

! **************************************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
!> \note
! **************************************************************************************************
   SUBROUTINE scf(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'scf', routineP = moduleN//':'//routineN

      INTEGER                                            :: max_scf_tmp
      LOGICAL                                            :: converged, outer_scf_loop, should_stop
      LOGICAL, SAVE                                      :: first_step_flag = .TRUE.
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, input=input, &
                      dft_control=dft_control, scf_control=scf_control)
      IF (scf_control%max_scf > 0) THEN

         dft_section => section_vals_get_subs_vals(input, "DFT")
         scf_section => section_vals_get_subs_vals(dft_section, "SCF")

         IF (.NOT. ASSOCIATED(scf_env)) THEN
            CALL qs_scf_env_initialize(qs_env, scf_env)
            ! Moved here from qs_scf_env_initialize to be able to have more scf_env
            CALL set_qs_env(qs_env, scf_env=scf_env)
            CALL scf_env_release(scf_env)
            CALL get_qs_env(qs_env=qs_env, scf_env=scf_env)
         ELSE
            CALL qs_scf_env_initialize(qs_env, scf_env)
         ENDIF

         IF ((scf_control%density_guess .EQ. history_guess) .AND. (first_step_flag)) THEN
            max_scf_tmp = scf_control%max_scf
            scf_control%max_scf = 1
            outer_scf_loop = scf_control%outer_scf%have_scf
            scf_control%outer_scf%have_scf = .FALSE.
         END IF

         IF (.NOT. dft_control%qs_control%cdft) THEN
            CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                converged=converged, should_stop=should_stop)
         ELSE
            ! Third SCF loop needed for CDFT with OT to properly restart OT inner loop
            CALL cdft_scf(qs_env=qs_env, should_stop=should_stop)
         END IF

         ! If SCF has not converged, then we should not start MP2
         IF (ASSOCIATED(qs_env%mp2_env)) qs_env%mp2_env%hf_fail = .NOT. converged

         !   *** add the converged wavefunction to the wavefunction history
         IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
             ((scf_control%density_guess .NE. history_guess) .OR. &
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
         ELSE IF ((scf_control%density_guess .EQ. history_guess) .AND. &
                  (first_step_flag)) THEN
            scf_control%max_scf = max_scf_tmp
            scf_control%outer_scf%have_scf = outer_scf_loop
            first_step_flag = .FALSE.
         END IF

         ! *** compute properties that depend on the converged wavefunction
         IF (.NOT. (should_stop)) CALL qs_scf_compute_properties(qs_env, dft_section)

         ! *** cleanup
         CALL scf_env_cleanup(scf_env)
         IF (dft_control%qs_control%cdft) &
            CALL cdft_control_cleanup(dft_control%qs_control%cdft_control)

      END IF

   END SUBROUTINE scf

! **************************************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param should_stop ...
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!> \author Matthias Krack
!> \note
! **************************************************************************************************
   SUBROUTINE scf_env_do_scf(scf_env, scf_control, qs_env, converged, should_stop)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: converged, should_stop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_do_scf', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: description, name
      INTEGER :: ext_master_id, external_comm, handle, handle2, i_tmp, iatom, ic, ispin, &
         iter_count, output_unit, scf_energy_message_tag, total_steps
      LOGICAL :: diis_step, do_kpoints, energy_only, exit_inner_loop, exit_outer_loop, &
         inner_loop_converged, just_energy, outer_loop_converged
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), DIMENSION(3)                        :: res_val_3
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, rho, energy, &
               logger, qs_charges, ks_env, mos, atomic_kind_set, qs_kind_set, &
               particle_set, dft_section, input, &
               scf_section, para_env, results, kpoints, pw_env, rho_ao_kp)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)

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
                      input=input, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      results=results, &
                      pw_env=pw_env, &
                      para_env=para_env)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
         "SCF WAVEFUNCTION OPTIMIZATION"

      IF ((output_unit > 0) .AND. (.NOT. scf_control%use_ot)) THEN
         WRITE (UNIT=output_unit, &
                FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)") &
            "Step", "Update method", "Time", "Convergence", "Total energy", "Change", &
            REPEAT("-", 78)
      END IF
      CALL cp_add_iter_level(logger%iter_info, "QS_SCF")

      ! check for external communicator and if the inermediate energy should be sended
      res_val_3(:) = -1.0_dp
      description = "[EXT_SCF_ENER_COMM]"
      IF (test_for_result(results, description=description)) THEN
         CALL get_results(results, description=description, &
                          values=res_val_3, n_entries=i_tmp)
         CPASSERT(i_tmp .EQ. 3)
         IF (ALL(res_val_3(:) .LE. 0.0)) &
            CALL cp_abort(__LOCATION__, &
                          " Trying to access result ("//TRIM(description)// &
                          ") which is not correctly stored.")
      END IF
      external_comm = NINT(res_val_3(1))
      ext_master_id = NINT(res_val_3(2))
      scf_energy_message_tag = NINT(res_val_3(3))

      ! *** outer loop of the scf, can treat other variables,
      ! *** such as lagrangian multipliers
      scf_env%outer_scf%iter_count = 0
      iter_count = 0
      total_steps = 0
      energy%tot_old = 0.0_dp

      scf_outer_loop: DO

         CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
                            scf_section=scf_section)

         CALL qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)

         scf_loop: DO

            CALL timeset(routineN//"_inner_loop", handle2)

            scf_env%iter_count = scf_env%iter_count+1
            iter_count = iter_count+1
            CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter_count)

            IF (output_unit > 0) CALL m_flush(output_unit)

            total_steps = total_steps+1
            just_energy = energy_only

            CALL qs_ks_update_qs_env(qs_env, just_energy=just_energy, &
                                     calculate_forces=.FALSE.)

            ! print 'heavy weight' or relatively expensive quantities
            CALL qs_scf_loop_print(qs_env, scf_env, para_env)

            IF (do_kpoints) THEN
               ! kpoints
               CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step)
            ELSE
               ! Gamma points only
               CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only)
            END IF

            ! another heavy weight print object, print controlled by dft_section
            CALL qs_scf_write_mos(mos, atomic_kind_set, qs_kind_set, particle_set, dft_section)

            CALL qs_scf_density_mixing(scf_env, rho, para_env, diis_step)

            t2 = m_walltime()

            CALL qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy)

            IF (.NOT. just_energy) energy%tot_old = energy%total

            ! check for external communicator and if the intermediate energy should be sent
            IF (scf_energy_message_tag .GT. 0) THEN
               CALL mp_send(energy%total, ext_master_id, scf_energy_message_tag, external_comm)
            END IF

            CALL qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, exit_inner_loop, &
                                         inner_loop_converged, output_unit)

            ! In case we decide to exit we perform few more check to see if this one
            ! is really the last SCF step
            IF (exit_inner_loop) THEN

               CALL qs_scf_inner_finalize(scf_env, qs_env, diis_step, output_unit)

               CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                            outer_loop_converged, exit_outer_loop)

               ! Let's tag the last SCF cycle so we can print informations only of the last step
               IF (exit_outer_loop) CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter_count)

            END IF

            IF (do_kpoints) THEN
               CALL write_kpoints_restart(rho_ao_kp, kpoints, scf_env, dft_section, particle_set, qs_kind_set)
            ELSE
               ! Write Wavefunction restart file
               CALL write_mo_set(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set)
            END IF

            ! Exit if we have finished with the SCF inner loop
            IF (exit_inner_loop) THEN
               CALL timestop(handle2)
               EXIT scf_loop
            END IF

            IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       scf_section, "PRINT%ITERATION_INFO/TIME_CUMUL"), cp_p_file)) &
               t1 = m_walltime()

            ! mixing methods have the new density matrix in p_mix_new
            IF (scf_env%mixing_method > 0) THEN
               DO ic = 1, SIZE(rho_ao_kp, 2)
                  DO ispin = 1, dft_control%nspins
                     CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                     CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
                  END DO
               END DO
            END IF

            CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, &
                                   mix_rho=scf_env%mixing_method >= gspace_mixing_nr)

            CALL timestop(handle2)

         END DO scf_loop

         IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

         ! In case we use the OUTER SCF loop let's print some info..
         CALL qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)

         IF (exit_outer_loop) EXIT scf_outer_loop

         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)

      END DO scf_outer_loop

      converged = inner_loop_converged .AND. outer_loop_converged

      IF (dft_control%qs_control%cdft) &
         dft_control%qs_control%cdft_control%total_steps = &
         dft_control%qs_control%cdft_control%total_steps+total_steps

      IF (.NOT. converged) CPWARN("SCF run NOT converged")

      ! if needed copy mo_coeff dbcsr->fm for later use in post_scf!fm->dbcsr
      DO ispin = 1, SIZE(mos) !fm -> dbcsr
         IF (mos(ispin)%mo_set%use_mo_coeff_b) THEN !fm->dbcsr
            IF (.NOT. ASSOCIATED(mos(ispin)%mo_set%mo_coeff_b)) & !fm->dbcsr
               CPABORT("mo_coeff_b is not allocated") !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b, & !fm->dbcsr
                                  mos(ispin)%mo_set%mo_coeff) !fm -> dbcsr
         ENDIF !fm->dbcsr
      ENDDO !fm -> dbcsr

      IF (dft_control%qs_control%becke_restraint) THEN
         ! Check if Becke potential is not needed for force evaluation or CDFT and deallocate
         IF (.NOT. dft_control%qs_control%becke_control%save_pot .AND. &
             .NOT. dft_control%qs_control%cdft) THEN
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                      dft_control%qs_control%becke_control%becke_pot%pw)
            IF (dft_control%qs_control%becke_control%atomic_charges) THEN
               DO iatom = 1, dft_control%qs_control%becke_control%natoms
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%charge(iatom)%pw)
               END DO
               DEALLOCATE (dft_control%qs_control%becke_control%charge)
            END IF
            IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               IF (.NOT. ASSOCIATED(dft_control%qs_control%becke_control%cavity_mat)) THEN
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%cavity%pw)
               ELSE
                  DEALLOCATE (dft_control%qs_control%becke_control%cavity_mat)
               END IF
            END IF
            IF (ASSOCIATED(dft_control%qs_control%becke_control%charges_fragment)) &
               DEALLOCATE (dft_control%qs_control%becke_control%charges_fragment)
            IF (dft_control%qs_control%becke_control%constraint_type == cdft_combined_constraint) THEN
               SELECT CASE (dft_control%qs_control%becke_control%combined_type)
               CASE (cdft_combined_all)
                  ! Do nothing
               CASE (cdft_combined_acceptor, cdft_combined_donor)
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%combined_weight%pw)
               END SELECT
            END IF
            dft_control%qs_control%becke_control%need_pot = .TRUE.
            dft_control%qs_control%becke_control%external_control = .FALSE.
         END IF
      END IF

      CALL cp_rm_iter_level(logger%iter_info, level_name="QS_SCF")
      CALL timestop(handle)

   END SUBROUTINE scf_env_do_scf

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_loop', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nmo, number_of_OT_envs
      LOGICAL                                            :: do_rotation, has_unit_metric, is_full_all
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (scf_control, matrix_s, matrix_ks, dft_control, mos, mo_coeff)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)

      CALL get_qs_env(qs_env=qs_env, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      mos=mos)

      ! if using mo_coeff_b then copy to fm
      DO ispin = 1, SIZE(mos) !fm->dbcsr
         IF (mos(1)%mo_set%use_mo_coeff_b) THEN !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b, mos(ispin)%mo_set%mo_coeff) !fm->dbcsr
         ENDIF !fm->dbcsr
      ENDDO !fm->dbcsr

      ! this just guarantees that all mo_occupations match the eigenvalues, if smear
      DO ispin = 1, dft_control%nspins
         ! do not reset mo_occupations if the maximum overlap method is in use
         IF (.NOT. scf_control%diagonalization%mom) &
            CALL set_mo_occupation(mo_set=mos(ispin)%mo_set, &
                                   smear=scf_control%smear)
      ENDDO

      SELECT CASE (scf_env%method)
      CASE DEFAULT

         CPABORT("unknown scf method method:"//cp_to_string(scf_env%method))

      CASE (filter_matrix_diag_method_nr)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

      CASE (general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr)
         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

      CASE (ot_diag_method_nr)
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

         ! disable DFTB and SE for now
         IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
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
                                     scf_control%diagonalization%ot_settings%energy_gap, dft_control%nspins)

      CASE (block_davidson_diag_method_nr)
         ! Preconditioner initialized within the loop, when required
      CASE (ot_method_nr)
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
                     CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
                     IF (has_unit_metric) THEN
                        CALL make_basis_simple(mo_coeff, nmo)
                     ELSE
                        CALL make_basis_sm(mo_coeff, nmo, matrix_s(1)%matrix)
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ELSE
            ! dont need any dirty trick for the numerically stable irac algorithm.
         ENDIF

         IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN

            ! restricted calculations require just one set of OT orbitals
            number_of_OT_envs = dft_control%nspins
            IF (dft_control%restricted) number_of_OT_envs = 1

            ALLOCATE (scf_env%qs_ot_env(number_of_OT_envs))

            ! XXX Joost XXX should disentangle reading input from this part
            CALL ot_scf_read_input(scf_env%qs_ot_env, scf_section)

            ! keep a note that we are restricted
            IF (dft_control%restricted) THEN
               scf_env%qs_ot_env(1)%restricted = .TRUE.
               ! requires rotation
               IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) &
                  CALL cp_abort(__LOCATION__, &
                                "Restricted calculation with OT requires orbital rotation. Please "// &
                                "activate the OT%ROTATION keyword!")
            ELSE
               scf_env%qs_ot_env(:)%restricted = .FALSE.
            ENDIF

            ! might need the KS matrix to init properly
            CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                                     calculate_forces=.FALSE.)

            ! if an old preconditioner is still around (i.e. outer SCF is active),
            ! remove it if this could be worthwhile
            IF (.NOT. reuse_precond) &
               CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           dft_control%nspins)

            !
            ! preconditioning still needs to be done correctly with has_unit_metric
            ! notice that a big part of the preconditioning (S^-1) is fine anyhow
            !
            IF (has_unit_metric) THEN
               NULLIFY (orthogonality_metric)
            ELSE
               orthogonality_metric => matrix_s(1)%matrix
            ENDIF

            IF (.NOT. reuse_precond) &
               CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                           scf_env%qs_ot_env(1)%settings%energy_gap, dft_control%nspins, &
                                           has_unit_metric=has_unit_metric, &
                                           chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type)
            IF (reuse_precond) reuse_precond = .FALSE.

            CALL ot_scf_init(mo_array=mos, matrix_s=orthogonality_metric, &
                             broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma, &
                             qs_ot_env=scf_env%qs_ot_env, matrix_ks=matrix_ks(1)%matrix)

            SELECT CASE (scf_env%qs_ot_env (1)%settings%preconditioner_type)
            CASE (ot_precond_none)
            CASE (ot_precond_full_all, ot_precond_full_single_inverse)
               ! this will rotate the MOs to be eigen states, which is not compatible with rotation
               ! e.g. mo_derivs here do not yet include potentially different occupations numbers
               do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
               ! only full all needs rotation
               is_full_all = scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all
               CPASSERT(.NOT. (do_rotation .AND. is_full_all))
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(ispin)%preconditioner)
               ENDDO
            CASE (ot_precond_s_inverse, ot_precond_full_single)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               ENDDO
            CASE DEFAULT
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               ENDDO
            END SELECT
         ENDIF

         ! if we have non-uniform occupations we should be using rotation
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         DO ispin = 1, SIZE(mos)
            IF (.NOT. mos(ispin)%mo_set%uniform_occupation) THEN
               CPASSERT(do_rotation)
            ENDIF
         ENDDO
      END SELECT

      ! another safety check
      IF (dft_control%low_spin_roks) THEN
         CPASSERT(scf_env%method == ot_method_nr)
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         CPASSERT(do_rotation)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE init_scf_loop

! **************************************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param scf_env ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_cleanup(scf_env)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)

!   *** Release SCF work storage ***

      IF (ASSOCIATED(scf_env%scf_work1)) THEN
         DO ispin = 1, SIZE(scf_env%scf_work1)
            CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix)
         ENDDO
         DEALLOCATE (scf_env%scf_work1)
      ENDIF
      IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2)
      IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho)
      IF (ASSOCIATED(scf_env%ortho_m1)) CALL cp_fm_release(scf_env%ortho_m1)

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

! *** method dependent cleanup
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
      CASE DEFAULT
         CPABORT("unknown scf method method:"//cp_to_string(scf_env%method))
      END SELECT

      IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
         DEALLOCATE (scf_env%outer_scf%variables)
      ENDIF
      IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
         DEALLOCATE (scf_env%outer_scf%count)
      ENDIF
      IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
         DEALLOCATE (scf_env%outer_scf%gradient)
      ENDIF
      IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
         DEALLOCATE (scf_env%outer_scf%energy)
      ENDIF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
         DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE scf_env_cleanup

! **************************************************************************************************
!> \brief perform a CDFT scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param should_stop flag determing if calculation should stop
!> \par History
!>      12.2015 Created
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_scf(qs_env, should_stop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_scf', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ispin, output_unit
      LOGICAL                                            :: cdft_loop_converged, converged, &
                                                            exit_cdft_loop, first_iteration
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env, ks_env, energy, rho, matrix_s, rho_ao, cdft_control, logger, &
               dft_control, pw_env, auxbas_pw_pool, energy, ks_env, scf_env, dft_section, &
               input, scf_section, scf_control)
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
         CALL qs_scf_cdft_initial_info(output_unit, cdft_control, dft_control)
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
                             converged=converged, should_stop=should_stop)
         ! Decide whether to reuse the preconditioner on the next iteration
         IF (cdft_control%reuse_precond) THEN
            ! For convergence in exactly one step, the preconditioner is always reused (assuming max_reuse > 0)
            ! usually this means that the electronic structure has already converged to the correct state
            ! but the constraint optimizer keeps jumping over the optimal solution
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count == 1 &
                .AND. cdft_control%total_steps /= 1) &
               cdft_control%nreused = cdft_control%nreused-1
            ! SCF converged in less than precond_freq steps
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count .LE. cdft_control%precond_freq .AND. &
                cdft_control%total_steps /= 1 .AND. cdft_control%nreused .LT. cdft_control%max_reuse) THEN
               reuse_precond = .TRUE.
               cdft_control%nreused = cdft_control%nreused+1
            ELSE
               reuse_precond = .FALSE.
               cdft_control%nreused = 0
            END IF
         END IF
         ! Update history purging counters
         IF (first_iteration .AND. cdft_control%purge_history) THEN
            cdft_control%istep = cdft_control%istep+1
            IF (scf_env%outer_scf%iter_count .GT. 1) THEN
               cdft_control%nbad_conv = cdft_control%nbad_conv+1
               IF (cdft_control%nbad_conv .GE. cdft_control%purge_freq .AND. &
                   cdft_control%istep .GE. cdft_control%purge_offset) THEN
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
         ! Check if inverse Jacobian needs to be calculated
         CALL outer_loop_calculate_inverse_jacobian(qs_env)
         ! Optimize constraint
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
      END DO cdft_outer_loop

      IF (dft_control%qs_control%becke_restraint) THEN
         ! Store needed arrays for ET coupling calculation
         IF (cdft_control%do_et) THEN
            CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
            CALL dbcsr_init_p(cdft_control%wmat%matrix)
            CALL dbcsr_copy(cdft_control%wmat%matrix, matrix_s(1)%matrix, &
                            name="ET_RESTRAINT_MATRIX")
            CALL dbcsr_set(cdft_control%wmat%matrix, 0.0_dp)
            CALL integrate_v_rspace(dft_control%qs_control%becke_control%becke_pot, &
                                    hmat=cdft_control%wmat, qs_env=qs_env, &
                                    calculate_forces=.FALSE., &
                                    gapw=dft_control%qs_control%gapw)
            CALL dbcsr_init_p(cdft_control%matrix_s%matrix)
            CALL dbcsr_copy(cdft_control%matrix_s%matrix, matrix_s(1)%matrix, &
                            name="OVERLAP")
            NULLIFY (cdft_control%mo_coeff)
            ALLOCATE (cdft_control%mo_coeff(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               NULLIFY (cdft_control%mo_coeff(ispin)%matrix)
               CALL cp_fm_create(matrix=cdft_control%mo_coeff(ispin)%matrix, &
                                 matrix_struct=qs_env%mos(ispin)%mo_set%mo_coeff%matrix_struct, &
                                 name="MO_COEFF_A"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
               CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_set%mo_coeff, &
                                cdft_control%mo_coeff(ispin)%matrix)
            END DO
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
         END IF
         ! Forces are not needed => time to deallocate Becke
         IF (.NOT. dft_control%qs_control%becke_control%save_pot) THEN
            CALL get_qs_env(qs_env, pw_env=pw_env)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                      dft_control%qs_control%becke_control%becke_pot%pw)
            IF (dft_control%qs_control%becke_control%atomic_charges) THEN
               DO iatom = 1, dft_control%qs_control%becke_control%natoms
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%charge(iatom)%pw)
               END DO
               DEALLOCATE (dft_control%qs_control%becke_control%charge)
            END IF
            IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               IF (.NOT. ASSOCIATED(dft_control%qs_control%becke_control%cavity_mat)) THEN
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%cavity%pw)
               ELSE
                  DEALLOCATE (dft_control%qs_control%becke_control%cavity_mat)
               END IF
            END IF
            IF (ASSOCIATED(dft_control%qs_control%becke_control%charges_fragment)) &
               DEALLOCATE (dft_control%qs_control%becke_control%charges_fragment)
            IF (dft_control%qs_control%becke_control%constraint_type == cdft_combined_constraint) THEN
               SELECT CASE (dft_control%qs_control%becke_control%combined_type)
               CASE (cdft_combined_all)
                  ! Do nothing
               CASE (cdft_combined_acceptor, cdft_combined_donor)
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                            dft_control%qs_control%becke_control%combined_weight%pw)
               END SELECT
            END IF
            dft_control%qs_control%becke_control%need_pot = .TRUE.
            dft_control%qs_control%becke_control%external_control = .FALSE.
         END IF
      END IF

      IF (dft_control%qs_control%cdft) THEN
         IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
            IF (.NOT. cdft_control%save_pot) THEN
               CALL get_qs_env(qs_env, pw_env=pw_env)
               CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, cdft_control%weight%pw)
               cdft_control%need_pot = .TRUE.
            END IF
         END IF
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

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_cleanup', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(cdft_control%constraint%variables)) &
         DEALLOCATE (cdft_control%constraint%variables)
      IF (ASSOCIATED(cdft_control%constraint%count)) &
         DEALLOCATE (cdft_control%constraint%count)
      IF (ASSOCIATED(cdft_control%constraint%gradient)) &
         DEALLOCATE (cdft_control%constraint%gradient)
      IF (ASSOCIATED(cdft_control%constraint%energy)) &
         DEALLOCATE (cdft_control%constraint%energy)

   END SUBROUTINE cdft_control_cleanup

! **************************************************************************************************
!> \brief Calculates the finite difference inverse Jacobian
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \par History
!>      01.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE outer_loop_calculate_inverse_jacobian(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_calculate_inverse_jacobian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, iter_count, iwork, j, &
                                                            max_scf, nspins, nvar, nwork, &
                                                            output_unit, pwork
      LOGICAL                                            :: converged, should_stop
      REAL(KIND=dp)                                      :: dh, inv_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, step_multiplier
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_stashed
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input, scf_section

      CALL timeset(routineN, handle)

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               mos_stashed, ks_env, scf_env, scf_control, dft_control, &
               cdft_control, input, scf_section, inv_jacobian)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, input=input)

      SELECT CASE (scf_control%outer_scf%optimizer)
      CASE DEFAULT
         scf_control%outer_scf%build_jacobian = .FALSE.
      CASE (outer_scf_optimizer_newton)
         scf_control%outer_scf%build_jacobian = .TRUE.
      CASE (outer_scf_optimizer_broyden)
         ! TODO: Add Broyden variants where the explicit Jacobian is built
         scf_control%outer_scf%build_jacobian = .FALSE.
      END SELECT

      IF (scf_control%outer_scf%build_jacobian) THEN
         ! Get output_unit
         scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                            extension=".scfLog")
         ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
         nspins = dft_control%nspins
         ALLOCATE (mos_stashed(nspins))
         DO ispin = 1, nspins
            CALL duplicate_mo_set(mos_stashed(ispin)%mo_set, mos(ispin)%mo_set)
         END DO
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         p_rmpv => rho_ao_kp(:, 1)
         cdft_control => dft_control%qs_control%cdft_control
         IF (.NOT. ASSOCIATED(cdft_control)) &
            CALL cp_abort(__LOCATION__, "Optimizers that need the explicit Jacobian can"// &
                          " only be used together with a valid CDFT constraint")
         ! Allocate work
         nvar = SIZE(scf_env%outer_scf%variables, 1)
         max_scf = scf_control%outer_scf%max_scf+1
         iter_count = scf_env%outer_scf%iter_count
         ALLOCATE (gradient(nvar, max_scf))
         gradient = scf_env%outer_scf%gradient
         ALLOCATE (energy(max_scf))
         energy = scf_env%outer_scf%energy
         ALLOCATE (jacobian(nvar, nvar))
         jacobian = 0.0_dp
         ! Setup finite difference scheme
         SELECT CASE (scf_control%outer_scf%jacobian_type)
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, "Unknown Jacobian type: "// &
                          cp_to_string(scf_control%outer_scf%jacobian_type))
         CASE (jacobian_fd1)
            nwork = 0
            pwork = 1
            ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
            coeff(nwork) = -1.0_dp
            coeff(pwork) = 1.0_dp
            step_multiplier = 1.0_dp
            dh = scf_control%outer_scf%jacobian_step
         CASE (jacobian_fd1_backward)
            nwork = -1
            pwork = 0
            ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
            coeff(nwork) = -1.0_dp
            coeff(pwork) = 1.0_dp
            step_multiplier = -1.0_dp
            dh = scf_control%outer_scf%jacobian_step
         CASE (jacobian_fd2)
            nwork = 0
            pwork = 2
            ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
            coeff(0) = -3.0_dp/2.0_dp
            coeff(1) = 2.0_dp
            coeff(2) = -1.0_dp/2.0_dp
            step_multiplier = 1.0_dp
            step_multiplier(2) = 2.0_dp
            dh = 2.0_dp*scf_control%outer_scf%jacobian_step
         CASE (jacobian_fd2_backward)
            nwork = -2
            pwork = 0
            ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
            coeff(0) = 3.0_dp/2.0_dp
            coeff(-1) = -2.0_dp
            coeff(-2) = 1.0_dp/2.0_dp
            step_multiplier = -1.0_dp
            step_multiplier(-2) = -2.0_dp
            dh = 2.0_dp*scf_control%outer_scf%jacobian_step
         CASE (jacobian_fd1_central)
            nwork = -1
            pwork = 1
            ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
            coeff(nwork) = -1.0_dp/2.0_dp
            coeff(pwork) = 1.0_dp/2.0_dp
            step_multiplier = 1.0_dp
            step_multiplier(nwork) = -1.0_dp
            dh = 2.0_dp*scf_control%outer_scf%jacobian_step
         END SELECT
         DO i = 1, nvar
            jacobian(i, :) = coeff(0)*scf_env%outer_scf%gradient(i, iter_count)
         END DO
         ! Print some info
         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(/,A)") &
               " ============================== CDFT SCF LOOP =================================="
            WRITE (output_unit, FMT="(A)") &
               " Evaluating inverse Jacobian using finite differences"
            SELECT CASE (scf_control%outer_scf%jacobian_type)
            CASE (jacobian_fd1)
               WRITE (output_unit, '(A)') " Type               : First order forward difference"
            CASE (jacobian_fd1_backward)
               WRITE (output_unit, '(A)') " Type               : First order backward difference"
            CASE (jacobian_fd2)
               WRITE (output_unit, '(A)') " Type               : Second order forward difference"
            CASE (jacobian_fd2_backward)
               WRITE (output_unit, '(A)') " Type               : Second order backward difference"
            CASE (jacobian_fd1_central)
               WRITE (output_unit, '(A)') " Type               : First order central difference"
            CASE DEFAULT
               CALL cp_abort(__LOCATION__, "Unknown Jacobian type: "// &
                             cp_to_string(scf_control%outer_scf%jacobian_type))
            END SELECT
            WRITE (output_unit, '(A,ES12.4)') " Step size          : ", scf_control%outer_scf%jacobian_step
         END IF
         ! Calculate the Jacobian by perturbing each Lagrangian and recalculating the energy self-consistently
         DO i = 1, nvar
            DO iwork = nwork, pwork
               IF (iwork == 0) CYCLE
               scf_env%outer_scf%variables(i, iter_count+1) = scf_env%outer_scf%variables(i, iter_count)+ &
                                                              step_multiplier(iwork)* &
                                                              scf_control%outer_scf%jacobian_step
               CALL outer_loop_update_qs_env(qs_env, scf_env)
               CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
               CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                   converged=converged, should_stop=should_stop)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
               ! Update (iter_count + 1) element of gradient
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count+1
               CALL outer_loop_gradient(qs_env, scf_env)
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count-1
               ! Update Jacobian
               DO j = 1, nvar
                  jacobian(j, i) = jacobian(j, i)+coeff(iwork)*scf_env%outer_scf%gradient(j, iter_count+1)
               END DO
               ! Reset everything to last converged state
               scf_env%outer_scf%variables(i, iter_count+1) = 0.0_dp
               scf_env%outer_scf%gradient = gradient
               scf_env%outer_scf%energy = energy
               DO ispin = 1, nspins
                  CALL duplicate_mo_set(mos(ispin)%mo_set, mos_stashed(ispin)%mo_set)
                  IF (mos(ispin)%mo_set%use_mo_coeff_b) THEN
                     CALL copy_fm_to_dbcsr(mos(ispin)%mo_set%mo_coeff, &
                                           mos(ispin)%mo_set%mo_coeff_b)
                  ENDIF
                  CALL calculate_density_matrix(mos(ispin)%mo_set, &
                                                p_rmpv(ispin)%matrix)
               END DO
               CALL qs_rho_update_rho(rho, qs_env=qs_env)
               CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
            END DO
         END DO
         ! Finalize and invert Jacobian
         jacobian = jacobian/dh
         IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
            ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
         inv_jacobian => scf_env%outer_scf%inv_jacobian
         CALL invert_matrix(jacobian, inv_jacobian, inv_error)
         ! Release temporary storage
         DO ispin = 1, nspins
            CALL deallocate_mo_set(mos_stashed(ispin)%mo_set)
         END DO
         DEALLOCATE (mos_stashed, jacobian, gradient, energy, coeff, step_multiplier)
         IF (output_unit > 0) &
            WRITE (output_unit, FMT="(/,A)") &
            " ============================= JACOBIAN CALCULATED ============================="
      END IF
      CALL timestop(handle)

   END SUBROUTINE outer_loop_calculate_inverse_jacobian

END MODULE qs_scf

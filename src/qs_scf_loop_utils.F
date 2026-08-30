!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Utility routines for qs_scf
! **************************************************************************************************
MODULE qs_scf_loop_utils
   USE cp_control_types,                ONLY: dft_control_type,&
                                              hairy_probes_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_release_p,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_get_diag,&
                                              dbcsr_scale_by_vector
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE input_constants,                 ONLY: diag_update_method_adiis,&
                                              smear_fermi_dirac
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_set_mo_occupation
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_density_mixing_types,         ONLY: broyden_mixing_nr,&
                                              direct_mixing_nr,&
                                              gspace_mixing_nr,&
                                              modified_broyden_mixing_nr,&
                                              multisecant_mixing_nr,&
                                              new_pulay_mixing_nr,&
                                              no_mixing_nr,&
                                              pulay_mixing_nr
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fb_env_methods,               ONLY: fb_env_do_diag
   USE qs_gspace_mixing,                ONLY: gspace_mixing
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mixing_utils,                 ONLY: self_consistency_check
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_mom_methods,                  ONLY: do_mom_diag
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_ot,                           ONLY: qs_ot_fixed_n_energy_gradient,&
                                              qs_ot_get_orbitals_complex,&
                                              qs_ot_get_orbitals_ref_complex
   USE qs_ot_minimizer,                 ONLY: ot_mini
   USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                              ot_scf_mini
   USE qs_ot_types,                     ONLY: qs_ot_channel_index,&
                                              qs_ot_check_channel_context,&
                                              qs_ot_type
   USE qs_outer_scf,                    ONLY: outer_loop_gradient
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_diagonalization,          ONLY: &
        do_block_davidson_diag, do_block_krylov_diag, do_general_diag, do_general_diag_kp, &
        do_ot_diag, do_roks_diag, do_scf_diag_subspace, do_special_diag, &
        qs_scf_get_ot_kpoint_operator
   USE qs_scf_methods,                  ONLY: cp_sm_mix,&
                                              scf_env_density_mixing
   USE qs_scf_output,                   ONLY: qs_scf_print_summary
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_scf_env_type, &
        smeagol_method_nr, special_diag_method_nr
   USE scf_control_types,               ONLY: scf_control_type,&
                                              smear_type
   USE smeagol_interface,               ONLY: run_smeagol_emtrans
   USE smearing_utils,                  ONLY: Smearkp,&
                                              Smearkp2
   USE tblite_interface,                ONLY: tb_native_scc_mixer_active
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_loop_utils'

   PUBLIC :: qs_scf_set_loop_flags, &
             qs_scf_new_mos, qs_scf_new_mos_kp, &
             qs_scf_kp_ref_refresh_ready, &
             qs_scf_candidate_density_delta, qs_scf_commit_density_candidate, &
             qs_scf_density_mixing, qs_scf_check_inner_exit, &
             qs_scf_check_outer_exit, qs_scf_inner_finalize, qs_scf_rho_update

CONTAINS

! **************************************************************************************************
!> \brief decide whether a mature variable-occupation REF chart should be rebuilt
!> \param accepted_searches accepted search directions in the current REF chart
!> \param residual current OT residual
!> \param eps_scf requested SCF threshold
!> \param base_state whether the current point can start a new search direction
!> \param allow_exit_window whether a prior physical REF rebuild permits an exit-only rebuild
!> \return true inside the asymptotic residual window after sufficient accepted history
! **************************************************************************************************
   PURE FUNCTION qs_scf_kp_ref_refresh_ready( &
      accepted_searches, residual, eps_scf, base_state, allow_exit_window) &
      RESULT(refresh)

      INTEGER, INTENT(IN)                                :: accepted_searches
      REAL(KIND=dp), INTENT(IN)                          :: residual, eps_scf
      LOGICAL, INTENT(IN)                                :: base_state
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_exit_window
      LOGICAL                                            :: refresh

      INTEGER, PARAMETER                                 :: min_accepted_searches = 12

      LOGICAL                                            :: exit_window
      REAL(KIND=dp)                                      :: eps

      eps = MAX(eps_scf, EPSILON(1.0_dp))
      exit_window = .FALSE.
      IF (PRESENT(allow_exit_window)) exit_window = allow_exit_window
      refresh = base_state .AND. accepted_searches >= min_accepted_searches .AND. &
                ((residual > 10.0_dp*eps .AND. residual <= 0.2_dp*SQRT(eps)) .OR. &
                 (exit_window .AND. residual > 0.25_dp*eps .AND. residual <= eps))

   END FUNCTION qs_scf_kp_ref_refresh_ready

! **************************************************************************************************
!> \brief computes properties for a given hamiltonian using the current wfn
!> \param scf_env ...
!> \param diis_step ...
!> \param energy_only ...
!> \param just_energy ...
!> \param exit_inner_loop ...
! **************************************************************************************************
   SUBROUTINE qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      LOGICAL                                            :: diis_step, energy_only, just_energy, &
                                                            exit_inner_loop

! Some flags needed to be set at the beginning of the loop

      diis_step = .FALSE.
      energy_only = .FALSE.
      just_energy = .FALSE.

      ! SCF loop, optimisation of the wfn coefficients
      ! qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date here

      scf_env%iter_count = 0
      scf_env%adiis_check_next = .FALSE.
      scf_env%adiis_validated = .FALSE.
      scf_env%step_norm = 0.0_dp
      scf_env%raw_map_delta = 0.0_dp
      scf_env%oda_lambda = 0.0_dp
      scf_env%oda_energy = 0.0_dp
      scf_env%oda_gradient0 = 0.0_dp
      scf_env%oda_gradient1 = 0.0_dp
      scf_env%oda_evaluations = 0
      scf_env%oda_status = ""
      scf_env%raw_map_delta_valid = .FALSE.
      exit_inner_loop = .FALSE.

   END SUBROUTINE qs_scf_set_loop_flags

! **************************************************************************************************
!> \brief takes known energy and derivatives and produces new wfns
!>        and or density matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param energy_only ...
!> \param probe ...
! **************************************************************************************************
   SUBROUTINE qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, &
                             energy_only, probe)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: scf_section
      LOGICAL                                            :: diis_step, energy_only
      TYPE(hairy_probes_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: probe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_new_mos'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: disable_diis, has_unit_metric, &
                                                            skip_diag_sub
      REAL(KIND=dp)                                      :: saved_eps_diis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (energy, ks_env, matrix_ks, matrix_s, rho, mos, dft_control)

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, energy=energy, &
                      ks_env=ks_env, &
                      matrix_ks=matrix_ks, rho=rho, mos=mos, &
                      dft_control=dft_control, &
                      has_unit_metric=has_unit_metric)
      scf_env%iter_param = 0.0_dp

      disable_diis = dft_control%qs_control%xtb_control%do_tblite .AND. &
                     tb_native_scc_mixer_active(dft_control)
      IF (disable_diis) THEN
         saved_eps_diis = scf_control%eps_diis
         scf_control%eps_diis = 0.0_dp
      END IF

      ! transfer total_zeff_corr from qs_env to scf_env only if
      ! correct_el_density_dip is switched on [SGh]
      IF (dft_control%correct_el_density_dip) THEN
         scf_env%sum_zeff_corr = qs_env%total_zeff_corr
         IF (ABS(qs_env%total_zeff_corr) > 0.0_dp) THEN
            IF (scf_env%method /= general_diag_method_nr) THEN
               CALL cp_abort(__LOCATION__, &
                             "Please use ALGORITHM STANDARD in "// &
                             "SCF%DIAGONALIZATION if "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE ")
            ELSE IF (dft_control%roks) THEN
               CALL cp_abort(__LOCATION__, &
                             "Combination of "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE "// &
                             "is not implemented with ROKS")
            ELSE IF (scf_control%diagonalization%mom) THEN
               CALL cp_abort(__LOCATION__, &
                             "Combination of "// &
                             "CORE_CORRECTION /= 0.0 and "// &
                             "SURFACE_DIPOLE_CORRECTION TRUE "// &
                             "is not implemented with SCF%MOM")
            END IF
         END IF
      END IF

      SELECT CASE (scf_env%method)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf method: "// &
                       cp_to_string(scf_env%method))

         ! *************************************************************************
         ! Filter matrix diagonalisation: ugly implementation at this point of time
         ! *************************************************************************
      CASE (filter_matrix_diag_method_nr)

         IF (ABS(qs_env%total_zeff_corr) > 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "CORE_CORRECTION /= 0.0 plus SURFACE_DIPOLE_CORRECTION TRUE "// &
                          "requires SCF%DIAGONALIZATION: ALGORITHM STANDARD")
         END IF
         CALL fb_env_do_diag(scf_env%filter_matrix_env, qs_env, &
                             matrix_ks, matrix_s, scf_section, diis_step)

         ! Diagonlization in non orthonormal case
      CASE (general_diag_method_nr)
         IF (dft_control%roks) THEN
            CALL do_roks_diag(scf_env, mos, matrix_ks, matrix_s, &
                              scf_control, scf_section, diis_step, &
                              has_unit_metric)
         ELSE
            IF (scf_control%diagonalization%mom) THEN
               CALL do_mom_diag(scf_env, mos, matrix_ks, &
                                matrix_s, scf_control, scf_section, &
                                diis_step)
            ELSE
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  CALL do_general_diag(scf_env, mos, matrix_ks, &
                                       matrix_s, scf_control, scf_section, &
                                       diis_step, &
                                       probe)
               ELSE
                  CALL do_general_diag(scf_env, mos, matrix_ks, &
                                       matrix_s, scf_control, scf_section, &
                                       diis_step)
               END IF
            END IF
            IF (scf_control%do_diag_sub) THEN
               skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                               (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
               IF (.NOT. skip_diag_sub) THEN
                  CALL do_scf_diag_subspace(qs_env, scf_env, scf_env%subspace_env, mos, rho, &
                                            ks_env, scf_section, scf_control)
               END IF
            END IF
         END IF
         ! Diagonlization in orthonormal case
      CASE (special_diag_method_nr)
         IF (dft_control%roks) THEN
            CALL do_roks_diag(scf_env, mos, matrix_ks, matrix_s, &
                              scf_control, scf_section, diis_step, &
                              has_unit_metric)
         ELSE
            CALL do_special_diag(scf_env, mos, matrix_ks, &
                                 scf_control, scf_section, &
                                 diis_step)
         END IF
         ! OT diagonalization
      CASE (ot_diag_method_nr)
         CALL do_ot_diag(scf_env, mos, matrix_ks, matrix_s, &
                         scf_control, scf_section, diis_step)
         ! Block Krylov diagonlization
      CASE (block_krylov_diag_method_nr)
         IF ((scf_env%krylov_space%eps_std_diag > 0.0_dp) .AND. &
             (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%krylov_space%eps_std_diag)) THEN
            IF (scf_env%krylov_space%always_check_conv) THEN
               CALL do_block_krylov_diag(scf_env, mos, matrix_ks, &
                                         scf_control, scf_section, check_moconv_only=.TRUE.)
            END IF
            CALL do_general_diag(scf_env, mos, matrix_ks, &
                                 matrix_s, scf_control, scf_section, diis_step)
         ELSE
            CALL do_block_krylov_diag(scf_env, mos, matrix_ks, &
                                      scf_control, scf_section)
         END IF
         IF (scf_control%do_diag_sub) THEN
            skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                            (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
            IF (.NOT. skip_diag_sub) THEN
               CALL do_scf_diag_subspace(qs_env, scf_env, scf_env%subspace_env, mos, rho, &
                                         ks_env, scf_section, scf_control)
            END IF
         END IF
         ! Block Davidson diagonlization
      CASE (block_davidson_diag_method_nr)
         CALL do_block_davidson_diag(qs_env, scf_env, mos, matrix_ks, matrix_s, scf_control, &
                                     scf_section, .FALSE.)
         ! OT without diagonlization. Needs special treatment for SCP runs
      CASE (ot_method_nr)
         CALL qs_scf_loop_do_ot(qs_env, scf_env, scf_control%smear, mos, rho, &
                                qs_env%mo_derivs, energy%total, &
                                matrix_s, energy_only=energy_only, has_unit_metric=has_unit_metric)
      END SELECT
      IF (disable_diis) scf_control%eps_diis = saved_eps_diis

      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, SIZE(mos)
         energy%kTS = energy%kTS + mos(ispin)%kTS
         energy%efermi = energy%efermi + mos(ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos), KIND=dp)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_new_mos

! **************************************************************************************************
!> \brief Updates MOs and density matrix using diagonalization
!>        Kpoint code
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param diis_step ...
!> \param probe ...
!> \param ot_kp_subspace_refresh ...
!> \param allow_ot_kp_subspace_refresh ...
!> \param allow_ot_kp_exit_refresh ...
!> \param accepted_ot_kp_searches ...
!> \param added_mos_auto_grow ...
!> \param energy_only ...
! **************************************************************************************************
   SUBROUTINE qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, probe, &
                                ot_kp_subspace_refresh, allow_ot_kp_subspace_refresh, &
                                allow_ot_kp_exit_refresh, accepted_ot_kp_searches, &
                                added_mos_auto_grow, energy_only)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: diis_step
      TYPE(hairy_probes_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: probe
      LOGICAL, INTENT(OUT), OPTIONAL                     :: ot_kp_subspace_refresh
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_ot_kp_subspace_refresh, &
                                                            allow_ot_kp_exit_refresh
      INTEGER, INTENT(IN), OPTIONAL                      :: accepted_ot_kp_searches
      LOGICAL, INTENT(OUT), OPTIONAL                     :: added_mos_auto_grow, energy_only

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_new_mos_kp'

      INTEGER                                            :: accepted_searches, handle, ispin
      LOGICAL :: allow_exit_refresh, allow_refresh, base_state, disable_diis, has_unit_metric, &
         my_added_mos_auto_grow, ot_energy_only, refresh
      REAL(dp)                                           :: diis_error, residual, saved_eps_diis
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)

      allow_refresh = .FALSE.
      allow_exit_refresh = .FALSE.
      accepted_searches = 0
      my_added_mos_auto_grow = .FALSE.
      refresh = .FALSE.
      IF (PRESENT(allow_ot_kp_subspace_refresh)) allow_refresh = allow_ot_kp_subspace_refresh
      IF (PRESENT(allow_ot_kp_exit_refresh)) allow_exit_refresh = allow_ot_kp_exit_refresh
      IF (PRESENT(accepted_ot_kp_searches)) accepted_searches = accepted_ot_kp_searches
      IF (PRESENT(ot_kp_subspace_refresh)) ot_kp_subspace_refresh = .FALSE.
      IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .FALSE.
      IF (PRESENT(energy_only)) energy_only = .FALSE.

      NULLIFY (dft_control, energy, kpoints, matrix_ks, matrix_s)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy, kpoints=kpoints)
      scf_env%iter_param = 0.0_dp
      disable_diis = dft_control%qs_control%xtb_control%do_tblite .AND. &
                     tb_native_scc_mixer_active(dft_control)
      IF (disable_diis) THEN
         saved_eps_diis = scf_control%eps_diis
         scf_control%eps_diis = 0.0_dp
      END IF

      IF (dft_control%roks) THEN
         CPABORT("KP code: ROKS method not available: ")
      END IF

      SELECT CASE (scf_env%method)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "KP code: Unknown scf method: "// &
                       cp_to_string(scf_env%method))
      CASE (general_diag_method_nr)
         ! Diagonlization in non orthonormal case
         CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s)
         IF (scf_control%diagonalization%update_method == diag_update_method_adiis) THEN
            CPASSERT(ASSOCIATED(scf_env%scf_subspace_buffer))
            diis_error = scf_env%scf_subspace_buffer%last_old_fock_weight
         END IF
         IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
            scf_control%smear%do_smear = .FALSE.
            CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .TRUE., &
                                    diis_step, diis_error, qs_env, probe)
         ELSE
            CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .TRUE., &
                                    diis_step, diis_error, qs_env, &
                                    added_mos_auto_grow=my_added_mos_auto_grow)
         END IF
         IF (my_added_mos_auto_grow) THEN
            IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
            IF (disable_diis) scf_control%eps_diis = saved_eps_diis
            CALL timestop(handle)
            RETURN
         END IF
         IF (scf_control%diagonalization%update_method == diag_update_method_adiis .AND. &
             .NOT. diis_step) THEN
            scf_env%iter_param = diis_error
            IF (scf_env%scf_subspace_buffer%last_restart) THEN
               scf_env%iter_method = "ADIIS/Rst."
            ELSE
               scf_env%iter_method = "ADIIS/Diag."
            END IF
         ELSE IF (diis_step) THEN
            scf_env%iter_param = diis_error
            scf_env%iter_method = "DIIS/Diag."
         ELSE
            IF (scf_env%mixing_method == 0) THEN
               scf_env%iter_method = "NoMix/Diag."
            ELSE IF (scf_env%mixing_method == 1) THEN
               scf_env%iter_param = scf_env%p_mix_alpha
               scf_env%iter_method = "P_Mix/Diag."
            ELSE IF (scf_env%mixing_method > 1) THEN
               scf_env%iter_param = scf_env%mixing_store%alpha
               scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
            END IF
         END IF
      CASE (special_diag_method_nr)
         CALL get_qs_env(qs_env=qs_env, has_unit_metric=has_unit_metric)
         CPASSERT(has_unit_metric)
         ! Diagonlization in orthonormal case
         CALL cp_abort(__LOCATION__, &
                       "KP code: Scf method not available: "// &
                       cp_to_string(scf_env%method))
      CASE (ot_diag_method_nr, &
            block_krylov_diag_method_nr, &
            block_davidson_diag_method_nr)
         CALL cp_abort(__LOCATION__, &
                       "KP code: Scf method not available: "// &
                       cp_to_string(scf_env%method))
      CASE (ot_method_nr)
         CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s)
         base_state = .NOT. scf_env%qs_ot_env(1)%energy_only .AND. &
                      (scf_env%qs_ot_env(1)%line_search_count == 0 .OR. &
                       scf_env%qs_ot_env(1)%line_search_might_be_done)
         residual = MAX(scf_env%iter_delta, scf_env%qs_ot_env(1)%delta)
         refresh = allow_refresh .AND. scf_env%qs_ot_env(1)%settings%do_ener .AND. &
                   scf_env%qs_ot_env(1)%settings%occupation_preconditioner .AND. &
                   (dft_control%smear .OR. scf_control%smear%do_smear) .AND. &
                   qs_scf_kp_ref_refresh_ready(accepted_searches, residual, &
                                               scf_control%eps_scf, base_state, allow_exit_refresh)
         IF (refresh) THEN
            IF (PRESENT(ot_kp_subspace_refresh)) ot_kp_subspace_refresh = .TRUE.
         ELSE
            CALL qs_scf_loop_do_ot_kp(qs_env, scf_env, matrix_ks, matrix_s, &
                                      my_added_mos_auto_grow, ot_energy_only)
            IF (PRESENT(energy_only)) energy_only = ot_energy_only
            IF (my_added_mos_auto_grow) THEN
               IF (PRESENT(added_mos_auto_grow)) added_mos_auto_grow = .TRUE.
               IF (disable_diis) scf_control%eps_diis = saved_eps_diis
               CALL timestop(handle)
               RETURN
            END IF
            IF (allow_exit_refresh .AND. scf_env%iter_delta <= scf_control%eps_scf) THEN
               base_state = TRIM(scf_env%iter_method) == "OT CG"
               refresh = qs_scf_kp_ref_refresh_ready( &
                         accepted_searches, scf_env%iter_delta, scf_control%eps_scf, &
                         base_state, allow_exit_window=.TRUE.)
               IF (refresh .AND. PRESENT(ot_kp_subspace_refresh)) THEN
                  ot_kp_subspace_refresh = .TRUE.
               END IF
            END IF
         END IF
      CASE (smeagol_method_nr)
         ! SMEAGOL interface
         diis_step = .FALSE.
         IF (scf_env%mixing_method == 0) THEN
            scf_env%iter_method = "NoMix/SMGL"
         ELSE IF (scf_env%mixing_method == 1) THEN
            scf_env%iter_param = scf_env%p_mix_alpha
            scf_env%iter_method = "P_Mix/SMGL"
         ELSE IF (scf_env%mixing_method > 1) THEN
            scf_env%iter_param = scf_env%mixing_store%alpha
            scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/SMGL"
         END IF
         CALL run_smeagol_emtrans(qs_env, last=.FALSE., iter=scf_env%iter_count, rho_ao_kp=scf_env%p_mix_new)
      END SELECT
      IF (disable_diis) scf_control%eps_diis = saved_eps_diis

      CALL get_qs_env(qs_env=qs_env, energy=energy)
      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      mos => kpoints%kp_env(1)%kpoint_env%mos
      DO ispin = 1, SIZE(mos, 2)
         energy%kTS = energy%kTS + mos(1, ispin)%kTS
         energy%efermi = energy%efermi + mos(1, ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos, 2), KIND=dp)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_new_mos_kp

! **************************************************************************************************
!> \brief K-point OT driver for local complex spin/k-point channels
!> \param qs_env ...
!> \param scf_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param added_mos_auto_grow ...
!> \param energy_only ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_do_ot_kp(qs_env, scf_env, matrix_ks, matrix_s, added_mos_auto_grow, energy_only)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      LOGICAL, INTENT(OUT)                               :: added_mos_auto_grow, energy_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_do_ot_kp'

      INTEGER :: first_channel, handle, homo, ikpoint, ispin, local_channel, local_kpoint, nao, &
         nkp_groups, nkpoint, nlocal_channels, nmo, nspin_ot
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: transfer_step, use_real_wfn
      REAL(KIND=dp)                                      :: trial_kTS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: orbital_scaling_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trial_eigenvalues, trial_occupation
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      TYPE(cp_fm_struct_type), POINTER                   :: active_mo_struct, chc_struct
      TYPE(cp_fm_type)                                   :: active_mo_coeff, active_mo_coeff_im, &
                                                            chc_re, hc_im, hc_re
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_im
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hc, matrix_hc_im, &
                                                            matrix_hc_physical, &
                                                            matrix_hc_physical_im
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_target, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: matrix_k_im, matrix_k_re, matrix_s_im, &
                                                            matrix_s_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)
      added_mos_auto_grow = .FALSE.

      NULLIFY (active_mo_struct, chc_struct, density_target, dft_control, energy, fmwork, kpoints, &
               local_ot_env, matrix_hc, matrix_hc_im, matrix_k_im, matrix_k_re, matrix_s_im, &
               matrix_hc_physical, matrix_hc_physical_im, matrix_s_re, mo_coeff, mo_coeff_im, &
               occupation_numbers, para_env_inter_kp, rho, &
               rho_ao_kp, sab_nl, scf_control, wkp)

      CPASSERT(ASSOCIATED(scf_env%qs_ot_env))
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy, kpoints=kpoints, rho=rho, &
                      scf_control=scf_control)
      CALL get_kpoint_info(kpoints, nkp=nkpoint, wkp=wkp, use_real_wfn=use_real_wfn, &
                           kp_range=kp_range, nkp_groups=nkp_groups, para_env_inter_kp=para_env_inter_kp, &
                           sab_nl=sab_nl)
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))

      IF (use_real_wfn) THEN
         CALL cp_abort(__LOCATION__, &
                       "K-point OT driver currently supports complex k-point wavefunctions only.")
      END IF
      IF (scf_env%qs_ot_env(1)%settings%ot_algorithm /= "TOD" .AND. &
          scf_env%qs_ot_env(1)%settings%ot_algorithm /= "REF") THEN
         CALL cp_abort(__LOCATION__, &
                       "K-point OT driver supports OT%ALGORITHM STRICT or IRAC.")
      END IF
      kp => kpoints%kp_env(1)%kpoint_env
      nspin_ot = SIZE(kp%mos, 2)
      CALL qs_ot_check_channel_context(scf_env%qs_ot_env, nspin_ot, &
                                       nkpoint=nkpoint, &
                                       restricted=ANY(scf_env%qs_ot_env(:)%restricted), &
                                       require_kpoint=.TRUE., &
                                       kp_range=kp_range, &
                                       wkp=wkp, &
                                       require_local_state=.TRUE., &
                                       require_complex_state=.TRUE.)

      IF (nkp_groups > 1) THEN
         transfer_step = .FALSE.
         CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints, scf_env, scf_control, .FALSE., &
                                 transfer_step, store_ot_matrices=.TRUE., transfer_only=.TRUE.)
      END IF

      first_channel = qs_ot_channel_index(1, kp_range(1), nspin_ot)
      nlocal_channels = SIZE(kpoints%kp_env)*nspin_ot
      local_ot_env => scf_env%qs_ot_env(first_channel:first_channel + nlocal_channels - 1)
      local_ot_env(1)%etotal = energy%total

      ALLOCATE (matrix_hc(nlocal_channels))
      ALLOCATE (matrix_hc_im(nlocal_channels))
      IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
         ALLOCATE (matrix_hc_physical(nlocal_channels))
         ALLOCATE (matrix_hc_physical_im(nlocal_channels))
      END IF
      DO local_channel = 1, nlocal_channels
         NULLIFY (matrix_hc(local_channel)%matrix)
         NULLIFY (matrix_hc_im(local_channel)%matrix)
         ALLOCATE (matrix_hc(local_channel)%matrix)
         ALLOCATE (matrix_hc_im(local_channel)%matrix)
         CALL dbcsr_copy(matrix_hc(local_channel)%matrix, local_ot_env(local_channel)%matrix_x, name="matrix_hc")
         CALL dbcsr_copy(matrix_hc_im(local_channel)%matrix, local_ot_env(local_channel)%matrix_x_im, &
                         name="matrix_hc_im")
         CALL dbcsr_set(matrix_hc(local_channel)%matrix, 0.0_dp)
         CALL dbcsr_set(matrix_hc_im(local_channel)%matrix, 0.0_dp)
         IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
            NULLIFY (matrix_hc_physical(local_channel)%matrix)
            NULLIFY (matrix_hc_physical_im(local_channel)%matrix)
            ALLOCATE (matrix_hc_physical(local_channel)%matrix)
            ALLOCATE (matrix_hc_physical_im(local_channel)%matrix)
            CALL dbcsr_copy(matrix_hc_physical(local_channel)%matrix, &
                            local_ot_env(local_channel)%matrix_x, name="matrix_hc_physical")
            CALL dbcsr_copy(matrix_hc_physical_im(local_channel)%matrix, &
                            local_ot_env(local_channel)%matrix_x_im, name="matrix_hc_physical_im")
            CALL dbcsr_set(matrix_hc_physical(local_channel)%matrix, 0.0_dp)
            CALL dbcsr_set(matrix_hc_physical_im(local_channel)%matrix, 0.0_dp)
         END IF
      END DO

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL qs_scf_get_ot_kpoint_operator( &
               matrix_ks, kpoints, kp, MIN(ispin, SIZE(matrix_ks, 1)), &
               kp%ot_hmat(1, MIN(ispin, SIZE(kp%ot_hmat, 2))), &
               kp%ot_hmat(2, MIN(ispin, SIZE(kp%ot_hmat, 2))), matrix_k_re, matrix_k_im)

            CALL get_mo_set(kp%mos(1, ispin), homo=homo, mo_coeff=mo_coeff, nao=nao, nmo=nmo, &
                            occupation_numbers=occupation_numbers)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            IF (local_ot_env(local_channel)%settings%do_ener) homo = nmo
            CPASSERT(homo >= 1 .AND. homo <= nmo)
            IF (.NOT. local_ot_env(local_channel)%settings%do_ener) THEN
               IF (MAXVAL(ABS(occupation_numbers(1:homo) - occupation_numbers(1))) > &
                   100.0_dp*EPSILON(1.0_dp)) THEN
                  CALL cp_abort(__LOCATION__, &
                                "Fixed-occupation K-point OT requires uniform occupations within each channel.")
               END IF
               IF (occupation_numbers(1) <= 100.0_dp*EPSILON(1.0_dp)) THEN
                  CALL cp_abort(__LOCATION__, &
                                "Fixed-occupation K-point OT encountered an empty active orbital block.")
               END IF
               IF (homo < nmo) THEN
                  IF (MAXVAL(ABS(occupation_numbers(homo + 1:nmo))) > 100.0_dp*EPSILON(1.0_dp)) THEN
                     CALL cp_abort(__LOCATION__, &
                                   "Fixed-occupation K-point OT cannot optimize occupied padded MO columns.")
                  END IF
               END IF
            END IF
            CALL cp_fm_struct_create(active_mo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                                     ncol_global=homo)
            CALL cp_fm_create(active_mo_coeff, active_mo_struct)
            CALL cp_fm_create(active_mo_coeff_im, active_mo_struct)
            CALL cp_fm_to_fm(mo_coeff, active_mo_coeff, homo)
            CALL cp_fm_to_fm(mo_coeff_im, active_mo_coeff_im, homo)
            CALL cp_fm_create(hc_re, active_mo_struct)
            CALL cp_fm_create(hc_im, active_mo_struct)

            CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, active_mo_coeff, hc_re, homo)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, active_mo_coeff_im, hc_re, homo, &
                                         alpha=-1.0_dp, beta=1.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_re, active_mo_coeff_im, hc_im, homo)
            CALL cp_dbcsr_sm_fm_multiply(matrix_k_im, active_mo_coeff, hc_im, homo, &
                                         alpha=1.0_dp, beta=1.0_dp)

            IF (local_ot_env(local_channel)%settings%do_ener) THEN
               CALL cp_fm_struct_create(chc_struct, template_fmstruct=active_mo_struct, &
                                        nrow_global=homo, ncol_global=homo)
               CALL cp_fm_create(chc_re, chc_struct)
               CALL cp_fm_gemm("T", "N", homo, homo, nao, 1.0_dp, &
                               active_mo_coeff, hc_re, 0.0_dp, chc_re)
               CALL cp_fm_gemm("T", "N", homo, homo, nao, 1.0_dp, &
                               active_mo_coeff_im, hc_im, 1.0_dp, chc_re)
               CALL copy_fm_to_dbcsr(chc_re, local_ot_env(local_channel)%rot_mat_chc)
               CALL dbcsr_get_diag(local_ot_env(local_channel)%rot_mat_chc, &
                                   local_ot_env(local_channel)%ener_rayleigh)
               CPASSERT(ASSOCIATED(local_ot_env(local_channel)%para_env))
               CALL local_ot_env(local_channel)%para_env%sum( &
                  local_ot_env(local_channel)%ener_rayleigh)
               CALL cp_fm_release(chc_re)
               CALL cp_fm_struct_release(chc_struct)
            END IF

            CALL copy_fm_to_dbcsr(hc_re, matrix_hc(local_channel)%matrix)
            CALL copy_fm_to_dbcsr(hc_im, matrix_hc_im(local_channel)%matrix)

            CALL cp_fm_release(active_mo_coeff)
            CALL cp_fm_release(active_mo_coeff_im)
            CALL cp_fm_release(hc_re)
            CALL cp_fm_release(hc_im)
            CALL cp_fm_struct_release(active_mo_struct)
            CALL dbcsr_release_p(matrix_k_re)
            CALL dbcsr_release_p(matrix_k_im)
         END DO
      END DO

      IF (local_ot_env(1)%settings%do_ener) THEN
         CALL qs_scf_loop_collect_ot_kp_trial_smearing( &
            kpoints, local_ot_env, first_channel, nspin_ot, kp_range, wkp, para_env_inter_kp, &
            scf_control%smear, trial_eigenvalues, trial_occupation, trial_kTS)
      END IF

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), homo=homo, nmo=nmo, &
                            occupation_numbers=occupation_numbers)
            IF (local_ot_env(local_channel)%settings%do_ener) homo = nmo
            IF (local_ot_env(local_channel)%settings%occupation_preconditioner) THEN
               CALL dbcsr_copy(matrix_hc_physical(local_channel)%matrix, &
                               matrix_hc(local_channel)%matrix)
               CALL dbcsr_copy(matrix_hc_physical_im(local_channel)%matrix, &
                               matrix_hc_im(local_channel)%matrix)
               ALLOCATE (orbital_scaling_factor(homo))
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)* &
                                           trial_occupation(1:homo, ikpoint, ispin)
               CALL dbcsr_scale_by_vector(matrix_hc_physical(local_channel)%matrix, &
                                          alpha=orbital_scaling_factor, side='right')
               CALL dbcsr_scale_by_vector(matrix_hc_physical_im(local_channel)%matrix, &
                                          alpha=orbital_scaling_factor, side='right')
               DEALLOCATE (orbital_scaling_factor)
            END IF
            ALLOCATE (orbital_scaling_factor(homo))
            IF (local_ot_env(local_channel)%settings%occupation_preconditioner) THEN
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)
            ELSE IF (local_ot_env(local_channel)%settings%do_ener) THEN
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)* &
                                           trial_occupation(1:homo, ikpoint, ispin)
            ELSE
               orbital_scaling_factor(:) = 2.0_dp*wkp(ikpoint)*occupation_numbers(1:homo)
            END IF
            CALL dbcsr_scale_by_vector(matrix_hc(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            CALL dbcsr_scale_by_vector(matrix_hc_im(local_channel)%matrix, &
                                       alpha=orbital_scaling_factor, side='right')
            DEALLOCATE (orbital_scaling_factor)
         END DO
      END DO

      IF (local_ot_env(1)%initial_density_sync_pending) THEN
         ! The initial Hamiltonian can come from a density that does not match the
         ! fixed-rank orbital projector. Rebuild that density before starting a line search.
         IF (local_ot_env(1)%settings%ot_algorithm == "REF") THEN
            DO local_channel = 1, nlocal_channels
               CALL dbcsr_copy(matrix_hc(local_channel)%matrix, local_ot_env(local_channel)%matrix_x)
               CALL dbcsr_copy(matrix_hc_im(local_channel)%matrix, local_ot_env(local_channel)%matrix_x_im)
            END DO
         END IF
         scf_env%iter_delta = 2.0_dp*scf_control%eps_scf
         local_ot_env(1)%OT_METHOD_FULL = "OT Sync"
         local_ot_env(:)%initial_density_sync_pending = .FALSE.
      ELSE
         IF (local_ot_env(1)%settings%do_ener) THEN
            local_ot_env(1)%etotal = energy%total - energy%kTS + trial_kTS
            CALL qs_scf_loop_prepare_ot_kp_energy_gradient( &
               kpoints, local_ot_env, first_channel, nspin_ot, wkp, para_env_inter_kp, &
               scf_control%smear, trial_eigenvalues, trial_occupation)
         END IF
         IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
            CALL ot_mini(local_ot_env, matrix_hc, matrix_hc_im=matrix_hc_im, &
                         matrix_hc_physical=matrix_hc_physical, &
                         matrix_hc_physical_im=matrix_hc_physical_im, &
                         para_env_inter_kp=para_env_inter_kp)
         ELSE
            CALL ot_mini(local_ot_env, matrix_hc, matrix_hc_im=matrix_hc_im, &
                         para_env_inter_kp=para_env_inter_kp)
         END IF
         scf_env%iter_delta = local_ot_env(1)%delta
      END IF
      energy_only = local_ot_env(1)%energy_only

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(SIZE(kp%ot_smat) >= 2)
         CALL qs_scf_get_ot_kpoint_operator(matrix_s, kpoints, kp, 1, &
                                            kp%ot_smat(1), kp%ot_smat(2), &
                                            matrix_s_re, matrix_s_im)
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), homo=homo, mo_coeff=mo_coeff, eigenvalues=eigenvalues, &
                            nao=nao, nmo=nmo)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            IF (local_ot_env(local_channel)%settings%do_ener) homo = nmo
            SELECT CASE (local_ot_env(local_channel)%settings%ot_algorithm)
            CASE ("TOD")
               CALL qs_ot_get_orbitals_complex(matrix_hc(local_channel)%matrix, &
                                               matrix_hc_im(local_channel)%matrix, &
                                               matrix_s_re, matrix_s_im, &
                                               local_ot_env(local_channel))
            CASE ("REF")
               CALL qs_ot_get_orbitals_ref_complex(matrix_hc(local_channel)%matrix, &
                                                   matrix_hc_im(local_channel)%matrix, &
                                                   matrix_s_re, matrix_s_im, &
                                                   local_ot_env(local_channel), local_ot_env(1))
            CASE DEFAULT
               CPABORT("K-point OT orbital update requires ALGORITHM STRICT or IRAC")
            END SELECT
            CALL cp_fm_struct_create(active_mo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                                     ncol_global=homo)
            CALL cp_fm_create(active_mo_coeff, active_mo_struct)
            CALL cp_fm_create(active_mo_coeff_im, active_mo_struct)
            CALL copy_dbcsr_to_fm(matrix_hc(local_channel)%matrix, active_mo_coeff)
            CALL copy_dbcsr_to_fm(matrix_hc_im(local_channel)%matrix, active_mo_coeff_im)
            CALL cp_fm_to_fm(active_mo_coeff, mo_coeff, homo)
            CALL cp_fm_to_fm(active_mo_coeff_im, mo_coeff_im, homo)
            IF (local_ot_env(local_channel)%settings%do_ener) THEN
               eigenvalues(1:homo) = local_ot_env(local_channel)%ener_x(:)
               kp%mos(2, ispin)%eigenvalues(1:homo) = local_ot_env(local_channel)%ener_x(:)
            END IF
            CALL cp_fm_release(active_mo_coeff)
            CALL cp_fm_release(active_mo_coeff_im)
            CALL cp_fm_struct_release(active_mo_struct)
         END DO
         CALL dbcsr_release_p(matrix_s_re)
         CALL dbcsr_release_p(matrix_s_im)
      END DO

      IF (local_ot_env(1)%settings%do_ener) THEN
         CALL kpoint_set_mo_occupation(kpoints, scf_control%smear, &
                                       added_mos_auto=scf_control%added_mos_auto, &
                                       added_mos_auto_grow=added_mos_auto_grow)
      END IF

      IF (ALLOCATED(trial_eigenvalues)) DEALLOCATE (trial_eigenvalues)
      IF (ALLOCATED(trial_occupation)) DEALLOCATE (trial_occupation)

      IF (.NOT. added_mos_auto_grow) THEN
         fmwork => scf_env%scf_work1
         IF (ASSOCIATED(scf_env%p_mix_new)) THEN
            density_target => scf_env%p_mix_new
         ELSE
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            density_target => rho_ao_kp
         END IF
         CALL kpoint_density_matrices(kpoints)
         CALL kpoint_density_transform(kpoints, density_target, .FALSE., &
                                       matrix_s(1, 1)%matrix, sab_nl, fmwork, &
                                       overlap_rs=matrix_s)
      END IF
      scf_env%iter_method = scf_env%qs_ot_env(first_channel)%OT_METHOD_FULL
      scf_env%iter_param = scf_env%qs_ot_env(first_channel)%ds_min
      qs_env%broyden_adaptive_sigma = scf_env%qs_ot_env(first_channel)%broyden_adaptive_sigma

      DO local_channel = 1, nlocal_channels
         CALL dbcsr_release(matrix_hc(local_channel)%matrix)
         CALL dbcsr_release(matrix_hc_im(local_channel)%matrix)
         DEALLOCATE (matrix_hc(local_channel)%matrix)
         DEALLOCATE (matrix_hc_im(local_channel)%matrix)
         IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
            CALL dbcsr_release(matrix_hc_physical(local_channel)%matrix)
            CALL dbcsr_release(matrix_hc_physical_im(local_channel)%matrix)
            DEALLOCATE (matrix_hc_physical(local_channel)%matrix)
            DEALLOCATE (matrix_hc_physical_im(local_channel)%matrix)
         END IF
      END DO
      DEALLOCATE (matrix_hc, matrix_hc_im)
      IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
         DEALLOCATE (matrix_hc_physical, matrix_hc_physical_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_do_ot_kp

! **************************************************************************************************
!> \brief collect the global fixed-N trial smearing state used by k-point OT
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param kp_range ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
!> \param trial_eigenvalues ...
!> \param trial_occupation ...
!> \param trial_kTS ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_collect_ot_kp_trial_smearing( &
      kpoints, local_ot_env, first_channel, nspin_ot, kp_range, wkp, para_env_inter_kp, smear, &
      trial_eigenvalues, trial_occupation, trial_kTS)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      INTEGER, DIMENSION(2), INTENT(IN)                  :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: trial_eigenvalues, trial_occupation
      REAL(KIND=dp), INTENT(OUT)                         :: trial_kTS

      INTEGER                                            :: ikpoint, ispin, local_channel, &
                                                            local_kpoint, ne_a, ne_b, nelectron, &
                                                            nkpoint, nmo
      REAL(KIND=dp)                                      :: kTS, mu, nel
      REAL(KIND=dp), DIMENSION(2)                        :: mus
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(smear%do_smear)
      CPASSERT(smear%method == smear_fermi_dirac)
      CPASSERT(ASSOCIATED(para_env_inter_kp))
      CALL get_kpoint_info(kpoints, nkp=nkpoint)
      kp => kpoints%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nelectron=nelectron, nmo=nmo)
      ne_a = nelectron
      ne_b = 0
      IF (nspin_ot == 2) CALL get_mo_set(kp%mos(1, 2), nelectron=ne_b)

      ALLOCATE (trial_eigenvalues(nmo, nkpoint, nspin_ot))
      ALLOCATE (trial_occupation(nmo, nkpoint, nspin_ot))
      trial_eigenvalues(:, :, :) = 0.0_dp
      trial_occupation(:, :, :) = 0.0_dp
      mus(:) = 0.0_dp

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         CPASSERT(ikpoint >= kp_range(1) .AND. ikpoint <= kp_range(2))
         DO ispin = 1, nspin_ot
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            trial_eigenvalues(1:nmo, ikpoint, ispin) = &
               local_ot_env(local_channel)%ener_x(1:nmo)
         END DO
      END DO
      CALL para_env_inter_kp%sum(trial_eigenvalues)

      IF (nspin_ot == 1) THEN
         nel = REAL(nelectron, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 1), mus(1), trial_kTS, &
                      trial_eigenvalues(:, :, 1), nel, wkp, smear%electronic_temperature, &
                      2.0_dp, smear_fermi_dirac)
      ELSE IF (smear%fixed_mag_mom > 0.0_dp) THEN
         nel = REAL(ne_a, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 1), mus(1), kTS, &
                      trial_eigenvalues(:, :, 1), nel, wkp, smear%electronic_temperature, &
                      1.0_dp, smear_fermi_dirac)
         trial_kTS = kTS
         nel = REAL(ne_b, KIND=dp)
         CALL Smearkp(trial_occupation(:, :, 2), mus(2), kTS, &
                      trial_eigenvalues(:, :, 2), nel, wkp, smear%electronic_temperature, &
                      1.0_dp, smear_fermi_dirac)
         trial_kTS = trial_kTS + kTS
      ELSE
         nel = REAL(ne_a + ne_b, KIND=dp)
         CALL Smearkp2(trial_occupation, mu, trial_kTS, trial_eigenvalues, nel, wkp, &
                       smear%electronic_temperature, smear_fermi_dirac)
         mus(1:2) = mu
      END IF
   END SUBROUTINE qs_scf_loop_collect_ot_kp_trial_smearing

! **************************************************************************************************
!> \brief construct auxiliary-energy gradients on the global fixed-N Mermin surface
!> \param kpoints ...
!> \param local_ot_env ...
!> \param first_channel ...
!> \param nspin_ot ...
!> \param wkp ...
!> \param para_env_inter_kp ...
!> \param smear ...
!> \param trial_eigenvalues ...
!> \param trial_occupation ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_gradient( &
      kpoints, local_ot_env, first_channel, nspin_ot, wkp, para_env_inter_kp, smear, &
      trial_eigenvalues, trial_occupation)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: local_ot_env
      INTEGER, INTENT(IN)                                :: first_channel, nspin_ot
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(mp_para_env_type), POINTER                    :: para_env_inter_kp
      TYPE(smear_type), INTENT(IN)                       :: smear
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: trial_eigenvalues, trial_occupation

      INTEGER                                            :: group, ikpoint, imo, ispin, &
                                                            local_channel, local_kpoint, ngroups, &
                                                            nmo
      REAL(KIND=dp)                                      :: fixed_n_mean, maxocc, sigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fixed_n_residual, fixed_n_weight, &
                                                            preconditioner_count, &
                                                            preconditioner_shift, response_weight
      TYPE(kpoint_env_type), POINTER                     :: kp

      IF (.NOT. local_ot_env(1)%settings%do_ener) RETURN

      CPASSERT(smear%do_smear)
      CPASSERT(smear%method == smear_fermi_dirac)
      CPASSERT(ASSOCIATED(para_env_inter_kp))
      sigma = smear%electronic_temperature
      CPASSERT(sigma > EPSILON(sigma))

      ngroups = 1
      IF (nspin_ot == 2 .AND. smear%fixed_mag_mom > 0.0_dp) ngroups = 2
      ALLOCATE (fixed_n_residual(ngroups), fixed_n_weight(ngroups), &
                preconditioner_count(ngroups), preconditioner_shift(ngroups))
      fixed_n_residual(:) = 0.0_dp
      fixed_n_weight(:) = 0.0_dp
      preconditioner_count(:) = 0.0_dp
      preconditioner_shift(:) = 0.0_dp

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
            CPASSERT(maxocc > EPSILON(maxocc))
            CPASSERT(SIZE(local_ot_env(local_channel)%ener_x) == nmo)
            CPASSERT(SIZE(local_ot_env(local_channel)%ener_rayleigh) == nmo)
            DO imo = 1, nmo
               fixed_n_weight(group) = fixed_n_weight(group) + &
                                       wkp(ikpoint)*MAX(0.0_dp, trial_occupation(imo, ikpoint, ispin)* &
                                                        (maxocc - trial_occupation(imo, ikpoint, ispin))/(maxocc*sigma))
               fixed_n_residual(group) = fixed_n_residual(group) + &
                                         wkp(ikpoint)*MAX(0.0_dp, trial_occupation(imo, ikpoint, ispin)* &
                                                          (maxocc - trial_occupation(imo, ikpoint, ispin))/(maxocc*sigma))* &
                                         (local_ot_env(local_channel)%ener_rayleigh(imo) - &
                                          trial_eigenvalues(imo, ikpoint, ispin))
            END DO
         END DO
      END DO
      CALL para_env_inter_kp%sum(fixed_n_weight)
      CALL para_env_inter_kp%sum(fixed_n_residual)

      IF (local_ot_env(1)%settings%occupation_preconditioner) THEN
         DO local_kpoint = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(local_kpoint)%kpoint_env
            ikpoint = kp%nkpoint
            DO ispin = 1, nspin_ot
               group = MERGE(ispin, 1, ngroups == nspin_ot)
               local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
               CALL get_mo_set(kp%mos(1, ispin), nmo=nmo)
               fixed_n_mean = 0.0_dp
               IF (fixed_n_weight(group) > EPSILON(fixed_n_weight(group))) THEN
                  fixed_n_mean = fixed_n_residual(group)/fixed_n_weight(group)
               END IF
               local_ot_env(local_channel)%ener_preconditioned_gx(:) = fixed_n_mean - &
                                                                       (local_ot_env(local_channel)%ener_rayleigh(:) - &
                                                                        trial_eigenvalues(:, ikpoint, ispin))
               preconditioner_shift(group) = preconditioner_shift(group) + &
                                             SUM(local_ot_env(local_channel)%ener_preconditioned_gx)
               preconditioner_count(group) = preconditioner_count(group) + REAL(nmo, KIND=dp)
            END DO
         END DO
         CALL para_env_inter_kp%sum(preconditioner_shift)
         CALL para_env_inter_kp%sum(preconditioner_count)
      END IF

      DO local_kpoint = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(local_kpoint)%kpoint_env
         ikpoint = kp%nkpoint
         DO ispin = 1, nspin_ot
            group = MERGE(ispin, 1, ngroups == nspin_ot)
            local_channel = qs_ot_channel_index(ispin, ikpoint, nspin_ot) - first_channel + 1
            CALL get_mo_set(kp%mos(1, ispin), maxocc=maxocc, nmo=nmo)
            IF (local_ot_env(local_channel)%settings%occupation_preconditioner .AND. &
                preconditioner_count(group) > 0.0_dp) THEN
               local_ot_env(local_channel)%ener_preconditioned_gx(:) = &
                  local_ot_env(local_channel)%ener_preconditioned_gx(:) - &
                  preconditioner_shift(group)/preconditioner_count(group)
            END IF
            ALLOCATE (response_weight(nmo))
            response_weight(:) = wkp(ikpoint)*MAX(0.0_dp, trial_occupation(1:nmo, ikpoint, ispin)* &
                                                  (maxocc - trial_occupation(1:nmo, ikpoint, ispin))/ &
                                                  (maxocc*sigma))
            CALL qs_ot_fixed_n_energy_gradient( &
               local_ot_env(local_channel)%ener_rayleigh, trial_eigenvalues(:, ikpoint, ispin), &
               response_weight, fixed_n_weight(group), fixed_n_residual(group), &
               local_ot_env(local_channel)%ener_gx)
            DEALLOCATE (response_weight)
         END DO
      END DO

      DEALLOCATE (fixed_n_residual, fixed_n_weight, preconditioner_count, preconditioner_shift)

   END SUBROUTINE qs_scf_loop_prepare_ot_kp_energy_gradient

! **************************************************************************************************
!> \brief the inner loop of scf, specific to using to the orbital transformation method
!>       basically, in goes the ks matrix out goes a new p matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param smear ...
!> \param mos ...
!> \param rho ...
!> \param mo_derivs ...
!> \param total_energy ...
!> \param matrix_s ...
!> \param energy_only ...
!> \param has_unit_metric ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      2013    moved from qs_scf [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_do_ot(qs_env, scf_env, smear, mos, rho, mo_derivs, total_energy, &
                                matrix_s, energy_only, has_unit_metric)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(smear_type), POINTER                          :: smear
      TYPE(mo_set_type), DIMENSION(:), INTENT(INOUT)     :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
      REAL(KIND=dp), INTENT(IN)                          :: total_energy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      LOGICAL, INTENT(INOUT)                             :: energy_only
      LOGICAL, INTENT(IN)                                :: has_unit_metric

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_loop_do_ot'

      INTEGER                                            :: handle, ispin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric

      CALL timeset(routineN, handle)
      NULLIFY (rho_ao)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (has_unit_metric) THEN
         NULLIFY (orthogonality_metric)
      ELSE
         orthogonality_metric => matrix_s(1)%matrix
      END IF

      ! in case of LSD the first spin qs_ot_env will drive the minimization
      ! in the case of a restricted calculation, it will make sure the spin orbitals are equal

      CALL ot_scf_mini(mos, mo_derivs, smear, orthogonality_metric, &
                       total_energy, energy_only, scf_env%iter_delta, &
                       scf_env%qs_ot_env)

      DO ispin = 1, SIZE(mos)
         CALL set_mo_occupation(mo_set=mos(ispin), smear=smear)
      END DO

      DO ispin = 1, SIZE(mos)
         CALL calculate_density_matrix(mos(ispin), &
                                       rho_ao(ispin)%matrix, &
                                       use_dbcsr=.TRUE.)
      END DO

      scf_env%iter_method = scf_env%qs_ot_env(1)%OT_METHOD_FULL
      scf_env%iter_param = scf_env%qs_ot_env(1)%ds_min
      qs_env%broyden_adaptive_sigma = scf_env%qs_ot_env(1)%broyden_adaptive_sigma

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_do_ot

! **************************************************************************************************
!> \brief Measure the distance between the diagonalized candidate density and the current density.
!> \param scf_env Holds SCF environment information.
!> \param rho Current electron density.
!> \param para_env Parallel environment.
!> \param delta Maximum matrix-element difference.
! **************************************************************************************************
   SUBROUTINE qs_scf_candidate_density_delta(scf_env, rho, para_env, delta)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(OUT)                         :: delta

      INTEGER                                            :: ic, ispin
      REAL(KIND=dp)                                      :: tmp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      CPASSERT(ASSOCIATED(scf_env%p_mix_new))
      delta = 0.0_dp
      DO ic = 1, SIZE(scf_env%p_mix_new, 2)
         DO ispin = 1, SIZE(scf_env%p_mix_new, 1)
            ! p_mix=1 leaves the candidate unchanged while cp_sm_mix evaluates
            ! the same max-element density difference used by direct mixing.
            CALL cp_sm_mix(m1=scf_env%p_mix_new(ispin, ic)%matrix, &
                           m2=rho_ao_kp(ispin, ic)%matrix, p_mix=1.0_dp, &
                           delta=tmp, para_env=para_env)
            delta = MAX(delta, tmp)
         END DO
      END DO

   END SUBROUTINE qs_scf_candidate_density_delta

! **************************************************************************************************
!> \brief Commit the diagonalized candidate density without numerical mixing.
!> \param scf_env Holds SCF environment information.
!> \param rho Electron density to update.
! **************************************************************************************************
   SUBROUTINE qs_scf_commit_density_candidate(scf_env, rho)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: ic, ispin
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      CPASSERT(ASSOCIATED(scf_env%p_mix_new))
      DO ic = 1, SIZE(scf_env%p_mix_new, 2)
         DO ispin = 1, SIZE(scf_env%p_mix_new, 1)
            CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name)
            CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
         END DO
      END DO

   END SUBROUTINE qs_scf_commit_density_candidate

! **************************************************************************************************
!> \brief Performs the requested density mixing if any needed
!> \param scf_env   Holds SCF environment information
!> \param rho       All data for the electron density
!> \param para_env  Parallel environment
!> \param diis_step Did we do a DIIS step?
! **************************************************************************************************
   SUBROUTINE qs_scf_density_mixing(scf_env, rho, para_env, diis_step)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: diis_step

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      SELECT CASE (scf_env%mixing_method)
      CASE (direct_mixing_nr)
         CALL scf_env_density_mixing(scf_env%p_mix_new, &
                                     scf_env%mixing_store, rho_ao_kp, para_env, scf_env%iter_delta, scf_env%iter_count, &
                                     diis=diis_step)
      CASE (gspace_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, modified_broyden_mixing_nr, &
            multisecant_mixing_nr, new_pulay_mixing_nr)
         ! Compute the difference p_out-p_in
         CALL self_consistency_check(rho_ao_kp, scf_env%p_delta, para_env, scf_env%p_mix_new, &
                                     delta=scf_env%iter_delta)
      CASE (no_mixing_nr)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf mixing method: "// &
                       cp_to_string(scf_env%mixing_method))
      END SELECT

   END SUBROUTINE qs_scf_density_mixing

! **************************************************************************************************
!> \brief checks whether exit conditions for outer loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param outer_loop_converged ...
!> \param exit_outer_loop ...
! **************************************************************************************************
   SUBROUTINE qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                      outer_loop_converged, exit_outer_loop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: should_stop, outer_loop_converged, &
                                                            exit_outer_loop

      REAL(KIND=dp)                                      :: outer_loop_eps

      outer_loop_converged = .TRUE.
      IF (scf_control%outer_scf%have_scf) THEN
         ! We have an outer SCF loop...
         scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
         outer_loop_converged = .FALSE.

         CALL outer_loop_gradient(qs_env, scf_env)
         ! Multiple constraints: get largest deviation
         outer_loop_eps = SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count)**2))

         IF (outer_loop_eps < scf_control%outer_scf%eps_scf) outer_loop_converged = .TRUE.
      END IF

      exit_outer_loop = should_stop .OR. outer_loop_converged .OR. &
                        scf_env%outer_scf%iter_count > scf_control%outer_scf%max_scf

   END SUBROUTINE qs_scf_check_outer_exit

! **************************************************************************************************
!> \brief checks whether exit conditions for inner loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param just_energy ...
!> \param exit_inner_loop ...
!> \param inner_loop_converged ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, just_energy, &
                                      exit_inner_loop, inner_loop_converged, output_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: should_stop, just_energy, &
                                                            exit_inner_loop, inner_loop_converged
      INTEGER                                            :: output_unit

      LOGICAL                                            :: update_validated

      inner_loop_converged = .FALSE.
      exit_inner_loop = .FALSE.
      update_validated = scf_control%diagonalization%update_method /= diag_update_method_adiis
      IF (.NOT. update_validated) THEN
         CPASSERT(ASSOCIATED(scf_env%scf_subspace_buffer))
         update_validated = scf_env%adiis_validated .OR. &
                            scf_env%scf_subspace_buffer%diis_weight >= 1.0_dp
      END IF

      CALL external_control(should_stop, "SCF", target_time=qs_env%target_time, &
                            start_time=qs_env%start_time)
      IF (scf_env%iter_delta < scf_control%eps_scf .AND. update_validated) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,I5,A/)") &
               "*** SCF run converged in ", scf_env%iter_count, " steps ***"
         END IF
         inner_loop_converged = .TRUE.
         exit_inner_loop = .TRUE.
      ELSE IF (should_stop .OR. scf_env%iter_count >= scf_control%max_scf) THEN
         inner_loop_converged = .FALSE.
         IF (just_energy) THEN
            exit_inner_loop = .FALSE.
         ELSE
            exit_inner_loop = .TRUE.
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T3,A,I5,A/)") &
                  "Leaving inner SCF loop after reaching ", scf_env%iter_count, " steps."
            END IF
         END IF
      END IF

   END SUBROUTINE qs_scf_check_inner_exit

! **************************************************************************************************
!> \brief undoing density mixing. Important upon convergence
!> \param scf_env ...
!> \param rho ...
!> \param dft_control ...
!> \param para_env ...
!> \param diis_step ...
! **************************************************************************************************
   SUBROUTINE qs_scf_undo_mixing(scf_env, rho, dft_control, para_env, diis_step)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: diis_step

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: ic, ispin, nc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp

      NULLIFY (rho_ao_kp)

      IF (scf_env%mixing_method > 0) THEN
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         nc = SIZE(scf_env%p_mix_new, 2)
         SELECT CASE (scf_env%mixing_method)
         CASE (direct_mixing_nr)
            CALL scf_env_density_mixing(scf_env%p_mix_new, scf_env%mixing_store, &
                                        rho_ao_kp, para_env, scf_env%iter_delta, &
                                        scf_env%iter_count, diis=diis_step, &
                                        invert=.TRUE.)
            DO ic = 1, nc
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                  CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
               END DO
            END DO
         CASE (gspace_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, modified_broyden_mixing_nr, &
               multisecant_mixing_nr, new_pulay_mixing_nr)
            DO ic = 1, nc
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                  CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
               END DO
            END DO
         END SELECT
      END IF
   END SUBROUTINE qs_scf_undo_mixing

! **************************************************************************************************
!> \brief Performs the updates rho (takes care of mixing as well)
!> \param rho ...
!> \param qs_env ...
!> \param scf_env ...
!> \param ks_env ...
!> \param mix_rho ...
! **************************************************************************************************
   SUBROUTINE qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      LOGICAL, INTENT(IN)                                :: mix_rho

      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      ! ** update qs_env%rho
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      ! ** Density mixing through density matrix or on the reciprocal space grid (exclusive)
      IF (mix_rho) THEN
         CALL gspace_mixing(qs_env, scf_env%mixing_method, scf_env%mixing_store, rho, &
                            para_env, scf_env%iter_count)

      END IF
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

   END SUBROUTINE qs_scf_rho_update

! **************************************************************************************************
!> \brief Performs the necessary steps before leaving innner scf loop
!> \param scf_env ...
!> \param qs_env ...
!> \param diis_step ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_inner_finalize(scf_env, qs_env, diis_step, output_unit)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: diis_step
      INTEGER, INTENT(IN)                                :: output_unit

      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (energy, rho, dft_control, ks_env)

      CALL get_qs_env(qs_env=qs_env, energy=energy, ks_env=ks_env, &
                      rho=rho, dft_control=dft_control, para_env=para_env, &
                      do_kpoints=do_kpoints)

      CALL cleanup_scf_loop(scf_env)

      ! now, print out energies and charges corresponding to the obtained wfn
      ! (this actually is not 100% consistent at this point)!
      CALL qs_scf_print_summary(output_unit, qs_env)

      CALL qs_scf_undo_mixing(scf_env, rho, dft_control, para_env, diis_step)

      !   *** update rspace rho since the mo changed
      !   *** this might not always be needed (i.e. no post calculation / no forces )
      !   *** but guarantees that rho and wfn are consistent at this point
      CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho=.FALSE.)

   END SUBROUTINE qs_scf_inner_finalize

! **************************************************************************************************
!> \brief perform cleanup operations at the end of an scf loop
!> \param scf_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE cleanup_scf_loop(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cleanup_scf_loop'

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
      CASE (ot_diag_method_nr)
         !
      CASE (general_diag_method_nr)
         !
      CASE (special_diag_method_nr)
         !
      CASE (block_krylov_diag_method_nr, block_davidson_diag_method_nr)
         !
      CASE (filter_matrix_diag_method_nr)
         !
      CASE (smeagol_method_nr)
         !
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "unknown scf method method:"// &
                       cp_to_string(scf_env%method))
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE cleanup_scf_loop

END MODULE qs_scf_loop_utils

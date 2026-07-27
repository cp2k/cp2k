!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2024 created [Maximilian Graml]
! **************************************************************************************************

MODULE bse_main

   USE bse_full_diag,                   ONLY: create_A,&
                                              create_B,&
                                              create_hermitian_form_of_ABBA,&
                                              diagonalize_A,&
                                              diagonalize_C
   USE bse_iterative,                   ONLY: do_subspace_iterations,&
                                              fill_local_3c_arrays
   USE bse_print,                       ONLY: print_BSE_start_flag
   USE bse_util,                        ONLY: adapt_BSE_input_params,&
                                              deallocate_matrices_bse,&
                                              determine_bse_combined_window,&
                                              estimate_BSE_resources,&
                                              get_bse_spin_block_layout,&
                                              mult_B_with_W,&
                                              truncate_BSE_matrices
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE group_dist_types,                ONLY: group_dist_d1_type
   USE input_constants,                 ONLY: bse_abba,&
                                              bse_both,&
                                              bse_fulldiag,&
                                              bse_iterdiag,&
                                              bse_screening_alpha,&
                                              bse_screening_tdhf,&
                                              bse_tda
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_main'

   PUBLIC :: start_bse_calculation

CONTAINS

! **************************************************************************************************
!> \brief Main subroutine managing BSE calculations
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param bse_lev_virt ...
!> \param gd_array ...
!> \param color_sub ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE start_bse_calculation(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_Q_static_bse_gemm, &
                                    Eigenval, Eigenval_scf, &
                                    homo, virtual, dimen_RI, dimen_RI_red, bse_lev_virt, &
                                    gd_array, color_sub, mp2_env, qs_env, mo_coeff, unit_nr)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S_ia_bse, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: Eigenval, Eigenval_scf
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      INTEGER, DIMENSION(:), INTENT(IN)                  :: bse_lev_virt
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'start_bse_calculation'

      INTEGER                                            :: first_active_mo, handle, ispin, &
                                                            last_active_mo, n_ov_joint, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo_red_arr, n_ov_arr, offsets_arr, &
                                                            virt_red_arr
      LOGICAL                                            :: my_do_abba, my_do_fulldiag, &
                                                            my_do_iterat_diag, my_do_tda
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced, Eigenval_reduced_1, &
                                                            Eigenval_reduced_2, &
                                                            Eigenval_reduced_joint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_iaQ_bse_local
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                            fm_inv_sqrt_A_minus_B, fm_Q_copy, &
                                                            fm_sqrt_A_minus_B
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fm_mat_S_ab_trunc_arr, &
         fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_mat_S_ia_trunc_arr, &
         fm_mat_S_ij_trunc_arr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      para_env => fm_mat_S_ia_bse(1)%matrix_struct%para_env

      my_do_fulldiag = .FALSE.
      my_do_iterat_diag = .FALSE.
      my_do_tda = .FALSE.
      my_do_abba = .FALSE.
      !Method: Iterative or full diagonalization
      SELECT CASE (mp2_env%bse%bse_diag_method)
      CASE (bse_iterdiag)
         my_do_iterat_diag = .TRUE.
         !MG: Basics of the Davidson solver are implemented, but not rigorously checked.
         CPABORT("Iterative BSE not yet implemented")
      CASE (bse_fulldiag)
         my_do_fulldiag = .TRUE.
      END SELECT
      !Approximation: TDA and/or full ABBA matrix
      SELECT CASE (mp2_env%bse%flag_tda)
      CASE (bse_tda)
         my_do_tda = .TRUE.
      CASE (bse_abba)
         my_do_abba = .TRUE.
      CASE (bse_both)
         my_do_tda = .TRUE.
         my_do_abba = .TRUE.
      END SELECT

      CALL print_BSE_start_flag(my_do_tda, my_do_abba, unit_nr)

      ! Link BSE debug flag against debug print level
      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level == debug_print_level) THEN
         mp2_env%bse%bse_debug_print = .TRUE.
      END IF

      CALL fm_mat_S_ia_bse(1)%matrix_struct%para_env%sync()

      ALLOCATE (homo_red_arr(nspins), virt_red_arr(nspins))
      ALLOCATE (fm_mat_S_ia_trunc_arr(nspins), fm_mat_S_ij_trunc_arr(nspins), fm_mat_S_ab_trunc_arr(nspins))
      ALLOCATE (fm_mat_S_bar_ia_bse_arr(nspins), fm_mat_S_bar_ij_bse_arr(nspins))

      IF (nspins > 1) THEN
         CALL cp_warn(__LOCATION__, &
                      "Open-shell (UKS/LSD) BSE is a recent addition and has not been "// &
                      "extensively validated. Verify results carefully before using them "// &
                      "for production calculations.")
         ! === Open-shell path: combined-window cutoff, per-spin truncation, joint A (+B for ABBA) ===
         ! Determine union of per-spin active MO windows
         CALL determine_bse_combined_window(Eigenval_scf(:, 1, :), homo, virtual, &
                                            mp2_env%bse%bse_cutoff_occ, &
                                            mp2_env%bse%bse_cutoff_empty, &
                                            first_active_mo, last_active_mo)
         DO ispin = 1, nspins
            CALL truncate_BSE_matrices(fm_mat_S_ia_bse(ispin), fm_mat_S_ij_bse(ispin), &
                                       fm_mat_S_ab_bse(ispin), &
                                       fm_mat_S_ia_trunc_arr(ispin), fm_mat_S_ij_trunc_arr(ispin), &
                                       fm_mat_S_ab_trunc_arr(ispin), &
                                       Eigenval_scf(:, 1, ispin), Eigenval(:, 1, ispin), &
                                       Eigenval_reduced, homo(ispin), virtual(ispin), dimen_RI, &
                                       unit_nr, bse_lev_virt(ispin), homo_red_arr(ispin), &
                                       virt_red_arr(ispin), mp2_env, &
                                       homo_incl_in=first_active_mo, &
                                       virt_incl_in=last_active_mo - homo(ispin))
            IF (ispin == 1) THEN
               ALLOCATE (Eigenval_reduced_1(SIZE(Eigenval_reduced)))
               Eigenval_reduced_1(:) = Eigenval_reduced(:)
            ELSE
               ALLOCATE (Eigenval_reduced_2(SIZE(Eigenval_reduced)))
               Eigenval_reduced_2(:) = Eigenval_reduced(:)
            END IF
            DEALLOCATE (Eigenval_reduced)
         END DO
         ! Flat eigenvalue layout: [sigma=1 levels, sigma=2 levels]
         ALLOCATE (Eigenval_reduced_joint(SIZE(Eigenval_reduced_1) + SIZE(Eigenval_reduced_2)))
         Eigenval_reduced_joint(1:SIZE(Eigenval_reduced_1)) = Eigenval_reduced_1
         Eigenval_reduced_joint(SIZE(Eigenval_reduced_1) + 1:) = Eigenval_reduced_2
         DEALLOCATE (Eigenval_reduced_1, Eigenval_reduced_2)

         ALLOCATE (n_ov_arr(nspins), offsets_arr(nspins))
         CALL get_bse_spin_block_layout(homo_red_arr, virt_red_arr, n_ov_arr, offsets_arr, n_ov_joint)

         CALL adapt_BSE_input_params(n_ov_joint, 1, unit_nr, mp2_env, qs_env)

         ! W: mult_B_with_W modifies Q in-place (Cholesky); copy the original Q for each spin call
         DO ispin = 1, nspins
            CALL cp_fm_create(fm_Q_copy, fm_mat_Q_static_bse_gemm%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_static_bse_gemm, fm_Q_copy)
            CALL mult_B_with_W(fm_mat_S_ij_trunc_arr(ispin), fm_mat_S_ia_trunc_arr(ispin), &
                               fm_mat_S_bar_ia_bse_arr(ispin), fm_mat_S_bar_ij_bse_arr(ispin), &
                               fm_Q_copy, dimen_RI_red, homo_red_arr(ispin), virt_red_arr(ispin))
            CALL cp_fm_release(fm_Q_copy)
         END DO
         CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

         IF (my_do_fulldiag) THEN
            CALL estimate_BSE_resources(n_ov_joint, unit_nr, my_do_abba, para_env, diag_runtime_est)
            IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
                mp2_env%bse%screening_method == bse_screening_alpha) THEN
               CALL create_A(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, &
                             fm_A_BSE, Eigenval_reduced_joint, unit_nr, &
                             homo_red_arr, virt_red_arr, dimen_RI, mp2_env, para_env, qs_env)
            ELSE
               CALL create_A(fm_mat_S_ia_trunc_arr, fm_mat_S_bar_ij_bse_arr, fm_mat_S_ab_trunc_arr, &
                             fm_A_BSE, Eigenval_reduced_joint, unit_nr, &
                             homo_red_arr, virt_red_arr, dimen_RI, mp2_env, para_env, qs_env)
            END IF
            IF (my_do_abba) THEN
               IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
                   mp2_env%bse%screening_method == bse_screening_alpha) THEN
                  CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ia_trunc_arr, fm_B_BSE, &
                                homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
               ELSE
                  CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_bar_ia_bse_arr, fm_B_BSE, &
                                homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
               END IF
               CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                                  unit_nr, mp2_env, diag_runtime_est)
               CALL cp_fm_release(fm_B_BSE)
            END IF
            NULLIFY (dft_control, tddfpt_control)
            CALL get_qs_env(qs_env, dft_control=dft_control)
            tddfpt_control => dft_control%tddfpt2_control
            ! 4th arg (homo_irred) = full per-spin occupied counts: per-spin absolute-MO labels for
            ! the amplitude table, and N_e = SUM(homo_irred) for the TRK print.
            IF (my_do_tda .AND. (.NOT. tddfpt_control%do_bse)) THEN
               CALL diagonalize_A(fm_A_BSE, homo_red_arr, virt_red_arr, homo, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_A_BSE)
            IF (my_do_abba) THEN
               CALL diagonalize_C(fm_C_BSE, homo_red_arr, virt_red_arr, homo, &
                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
               CALL cp_fm_release(fm_C_BSE)
            END IF
         END IF

         DO ispin = 1, nspins
            CALL cp_fm_release(fm_mat_S_bar_ia_bse_arr(ispin))
            CALL cp_fm_release(fm_mat_S_bar_ij_bse_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ia_trunc_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ij_trunc_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ab_trunc_arr(ispin))
         END DO
         IF (mp2_env%bse%do_nto_analysis) DEALLOCATE (mp2_env%bse%bse_nto_state_list_final)
         DEALLOCATE (Eigenval_reduced_joint, n_ov_arr, offsets_arr)

      ELSE
         ! === Closed-shell n_spin=1 path (bit-identical) ===
         CALL truncate_BSE_matrices(fm_mat_S_ia_bse(1), fm_mat_S_ij_bse(1), fm_mat_S_ab_bse(1), &
                                    fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                                    fm_mat_S_ab_trunc_arr(1), &
                                    Eigenval_scf(:, 1, 1), Eigenval(:, 1, 1), Eigenval_reduced, &
                                    homo(1), virtual(1), dimen_RI, unit_nr, &
                                    bse_lev_virt(1), homo_red_arr(1), virt_red_arr(1), mp2_env)
         CALL mult_B_with_W(fm_mat_S_ij_trunc_arr(1), fm_mat_S_ia_trunc_arr(1), &
                            fm_mat_S_bar_ia_bse_arr(1), fm_mat_S_bar_ij_bse_arr(1), &
                            fm_mat_Q_static_bse_gemm, dimen_RI_red, homo_red_arr(1), virt_red_arr(1))

         IF (my_do_iterat_diag) THEN
            CALL fill_local_3c_arrays(fm_mat_S_ab_trunc_arr(1), fm_mat_S_ia_trunc_arr(1), &
                                      fm_mat_S_bar_ia_bse_arr(1), fm_mat_S_bar_ij_bse_arr(1), &
                                      B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                      B_iaQ_bse_local, dimen_RI_red, homo_red_arr(1), &
                                      virt_red_arr(1), gd_array, color_sub, para_env)
         END IF

         CALL adapt_BSE_input_params(homo_red_arr(1), virt_red_arr(1), unit_nr, mp2_env, qs_env)

         IF (my_do_fulldiag) THEN
            n_ov_joint = homo_red_arr(1)*virt_red_arr(1)
            CALL estimate_BSE_resources(n_ov_joint, unit_nr, my_do_abba, para_env, diag_runtime_est)
            IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
                mp2_env%bse%screening_method == bse_screening_alpha) THEN
               CALL create_A(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, &
                             fm_A_BSE, Eigenval_reduced, unit_nr, &
                             homo_red_arr, virt_red_arr, dimen_RI, mp2_env, para_env, qs_env)
            ELSE
               CALL create_A(fm_mat_S_ia_trunc_arr, fm_mat_S_bar_ij_bse_arr, fm_mat_S_ab_trunc_arr, &
                             fm_A_BSE, Eigenval_reduced, unit_nr, &
                             homo_red_arr, virt_red_arr, dimen_RI, mp2_env, para_env, qs_env)
            END IF
            IF (my_do_abba) THEN
               IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
                   mp2_env%bse%screening_method == bse_screening_alpha) THEN
                  CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ia_trunc_arr, fm_B_BSE, &
                                homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
               ELSE
                  CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_bar_ia_bse_arr, fm_B_BSE, &
                                homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
               END IF
               CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                                  unit_nr, mp2_env, diag_runtime_est)
            END IF
            CALL cp_fm_release(fm_B_BSE)

            NULLIFY (dft_control, tddfpt_control)
            CALL get_qs_env(qs_env, dft_control=dft_control)
            tddfpt_control => dft_control%tddfpt2_control
            IF (my_do_tda .AND. (.NOT. tddfpt_control%do_bse)) THEN
               CALL diagonalize_A(fm_A_BSE, homo_red_arr, virt_red_arr, homo, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_A_BSE)
            IF (my_do_abba) THEN
               CALL diagonalize_C(fm_C_BSE, homo_red_arr, virt_red_arr, homo, &
                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_C_BSE)
         END IF

         CALL deallocate_matrices_bse(fm_mat_S_bar_ia_bse_arr(1), fm_mat_S_bar_ij_bse_arr(1), &
                                      fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                                      fm_mat_S_ab_trunc_arr(1), fm_mat_Q_static_bse_gemm, mp2_env)
         DEALLOCATE (Eigenval_reduced)
         IF (my_do_iterat_diag) THEN
            CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                        B_iaQ_bse_local, homo(1), virtual(1), &
                                        mp2_env%bse%bse_spin_config, unit_nr, &
                                        Eigenval(:, 1, 1), para_env, mp2_env)
            DEALLOCATE (B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, B_iaQ_bse_local)
         END IF

      END IF

      DEALLOCATE (homo_red_arr, virt_red_arr)
      DEALLOCATE (fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr)
      DEALLOCATE (fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A53)') 'BSE|', 'The BSE was successfully calculated. Have a nice day!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE start_bse_calculation

END MODULE bse_main

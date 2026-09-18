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

   USE bse_davidson,                    ONLY: abba_indefinite,&
                                              bse_davidson_abba_mk,&
                                              bse_davidson_abba_os,&
                                              bse_davidson_refcheck_abba,&
                                              bse_davidson_refcheck_tda,&
                                              bse_davidson_tda
   USE bse_full_diag,                   ONLY: create_A,&
                                              create_B,&
                                              create_hermitian_form_of_ABBA,&
                                              diagonalize_A,&
                                              diagonalize_C,&
                                              postprocess_bse
   USE bse_matvec,                      ONLY: bse_matvec_create,&
                                              bse_matvec_env_type,&
                                              bse_matvec_release,&
                                              bse_matvec_selfcheck
   USE bse_print,                       ONLY: print_BSE_start_flag
   USE bse_util,                        ONLY: adapt_BSE_input_params,&
                                              deallocate_matrices_bse,&
                                              estimate_BSE_resources,&
                                              get_bse_spin_block_layout,&
                                              mult_B_with_W,&
                                              truncate_BSE_matrices
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE group_dist_types,                ONLY: group_dist_d1_type
   USE input_constants,                 ONLY: &
        bse_abba, bse_abba_mk_davidson, bse_abba_olsen_stratmann, bse_both, bse_fulldiag, &
        bse_iterdiag, bse_screening_alpha, bse_screening_rpa, bse_screening_tdhf, bse_singlet, &
        bse_tda, bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mo_window,                       ONLY: combine_mo_windows,&
                                              determine_mo_window,&
                                              mo_window_type
   USE mp2_types,                       ONLY: mp2_type
   USE physcon,                         ONLY: evolt
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

      INTEGER                                            :: handle, ispin, n_ov_joint, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo_red_arr, n_ov_arr, offsets_arr, &
                                                            virt_red_arr
      LOGICAL                                            :: my_do_abba, my_do_fulldiag, &
                                                            my_do_iterat_diag, my_do_tda, &
                                                            spins_differ
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced, Eigenval_reduced_1, &
                                                            Eigenval_reduced_2, &
                                                            Eigenval_reduced_joint
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                            fm_inv_sqrt_A_minus_B, fm_Q_copy, &
                                                            fm_sqrt_A_minus_B
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fm_mat_S_ab_trunc_arr, &
         fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_mat_S_ia_trunc_arr, &
         fm_mat_S_ij_trunc_arr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_window_type)                               :: combined_window
      TYPE(mo_window_type), ALLOCATABLE, DIMENSION(:)    :: spin_windows
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      MARK_USED(gd_array)
      MARK_USED(color_sub)

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
      ALLOCATE (spin_windows(nspins))
      ALLOCATE (fm_mat_S_ia_trunc_arr(nspins), fm_mat_S_ij_trunc_arr(nspins), fm_mat_S_ab_trunc_arr(nspins))
      ALLOCATE (fm_mat_S_bar_ia_bse_arr(nspins), fm_mat_S_bar_ij_bse_arr(nspins))

      DO ispin = 1, nspins
         CALL determine_mo_window(Eigenval_scf(:, 1, ispin), homo(ispin) + virtual(ispin), homo(ispin), &
                                  mp2_env%bse%bse_cutoff_occ, mp2_env%bse%bse_cutoff_empty, &
                                  spin_windows(ispin))
      END DO

      IF (nspins > 1) THEN
         IF (my_do_iterat_diag) THEN
            CPABORT("Iterative BSE is not yet available for open-shell references")
         END IF
         CALL cp_warn(__LOCATION__, &
                      "Open-shell (UKS/LSD) BSE is a recent addition and has not been "// &
                      "extensively validated. Verify results carefully before using them "// &
                      "for production calculations.")
         ! === Open-shell path: combined-window cutoff, per-spin truncation, joint A (+B for ABBA) ===
         CALL combine_mo_windows(spin_windows, combined_window, spins_differ)
         IF (spins_differ) THEN
            CPWARN("BSE: spin-resolved active MO cutoff candidates differ; using combined window.")
         END IF
         DO ispin = 1, nspins
            CALL truncate_BSE_matrices(fm_mat_S_ia_bse(ispin), fm_mat_S_ij_bse(ispin), &
                                       fm_mat_S_ab_bse(ispin), &
                                       fm_mat_S_ia_trunc_arr(ispin), fm_mat_S_ij_trunc_arr(ispin), &
                                       fm_mat_S_ab_trunc_arr(ispin), &
                                       Eigenval_scf(:, 1, ispin), Eigenval(:, 1, ispin), &
                                       Eigenval_reduced, homo(ispin), virtual(ispin), dimen_RI, &
                                       unit_nr, bse_lev_virt(ispin), homo_red_arr(ispin), &
                                       virt_red_arr(ispin), mp2_env, combined_window, .FALSE.)
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
                                    bse_lev_virt(1), homo_red_arr(1), virt_red_arr(1), mp2_env, &
                                    spin_windows(1), .TRUE.)
         CALL mult_B_with_W(fm_mat_S_ij_trunc_arr(1), fm_mat_S_ia_trunc_arr(1), &
                            fm_mat_S_bar_ia_bse_arr(1), fm_mat_S_bar_ij_bse_arr(1), &
                            fm_mat_Q_static_bse_gemm, dimen_RI_red, homo_red_arr(1), virt_red_arr(1))

         CALL adapt_BSE_input_params(homo_red_arr(1), virt_red_arr(1), unit_nr, mp2_env, qs_env)

         IF (my_do_iterat_diag) THEN
            CALL solve_bse_iteratively(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, &
                                       fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, &
                                       fm_mat_S_bar_ij_bse_arr, Eigenval_reduced, &
                                       homo_red_arr, virt_red_arr, homo, dimen_RI, &
                                       my_do_tda, my_do_abba, mp2_env, qs_env, mo_coeff, unit_nr)
         END IF

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

      END IF

      DEALLOCATE (homo_red_arr, virt_red_arr)
      DEALLOCATE (spin_windows)
      DEALLOCATE (fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr)
      DEALLOCATE (fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A53)') 'BSE|', 'The BSE was successfully calculated. Have a nice day!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE start_bse_calculation

! **************************************************************************************************
!> \brief Lowest excitations of a closed-shell reference from the block Davidson solvers on top of
!>        the matrix-free application of A and B, post-processed as in the full diagonalization
!> \param fm_mat_S_ia_trunc_arr ...
!> \param fm_mat_S_ij_trunc_arr ...
!> \param fm_mat_S_ab_trunc_arr ...
!> \param fm_mat_S_bar_ia_bse_arr ...
!> \param fm_mat_S_bar_ij_bse_arr ...
!> \param Eigenval_reduced ...
!> \param homo_red_arr ...
!> \param virt_red_arr ...
!> \param homo ...
!> \param dimen_RI ...
!> \param do_tda ...
!> \param do_abba ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE solve_bse_iteratively(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, &
                                    fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, &
                                    fm_mat_S_bar_ij_bse_arr, Eigenval_reduced, &
                                    homo_red_arr, virt_red_arr, homo, dimen_RI, &
                                    do_tda, do_abba, mp2_env, qs_env, mo_coeff, unit_nr)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, &
         fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval_reduced
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo_red_arr, virt_red_arr, homo
      INTEGER, INTENT(IN)                                :: dimen_RI
      LOGICAL, INTENT(IN)                                :: do_tda, do_abba
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_bse_iteratively'
      INTEGER, PARAMETER                                 :: n_ov_max_debug = 20000

      CHARACTER(LEN=16)                                  :: margin_str, solver_str
      INTEGER                                            :: abba_status, handle, n_ov
      LOGICAL                                            :: do_debug
      REAL(KIND=dp)                                      :: ab_margin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exc_ens, exc_ens_abba
      TYPE(bse_matvec_env_type)                          :: mv_env
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_X, fm_X_abba, &
                                                            fm_X_ia, fm_Y_abba, fm_Y_ia
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%tddfpt2_control%do_bse .OR. dft_control%tddfpt2_control%do_bse_w_only .OR. &
          dft_control%tddfpt2_control%do_bse_gw_only) THEN
         CALL cp_abort(__LOCATION__, &
                       "BSE_DIAG_METHOD ITERDIAG cannot provide the BSE kernel of TDDFPT "// &
                       "(DO_BSE, DO_BSE_W_ONLY, DO_BSE_GW_ONLY). Use BSE_DIAG_METHOD FULLDIAG.")
      END IF

      para_env => fm_mat_S_ia_trunc_arr(1)%matrix_struct%para_env
      n_ov = homo_red_arr(1)*virt_red_arr(1)
      ! the debug checks hold the explicit N_ov x N_ov matrices
      do_debug = mp2_env%bse%bse_debug_print .AND. n_ov <= n_ov_max_debug
      IF (mp2_env%bse%bse_debug_print .AND. .NOT. do_debug .AND. unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A10,T13,A,T71,I10)') 'BSE|DEBUG|', &
            'Checks against the explicit A and B skipped, N_ov above', n_ov_max_debug
      END IF

      CALL create_matvec_env(mv_env, fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                             fm_mat_S_ab_trunc_arr(1), fm_mat_S_bar_ia_bse_arr(1), &
                             fm_mat_S_bar_ij_bse_arr(1), Eigenval_reduced, &
                             homo_red_arr(1), virt_red_arr(1), do_abba, mp2_env, unit_nr)

      IF (do_debug) THEN
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
         IF (do_abba) THEN
            IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
                mp2_env%bse%screening_method == bse_screening_alpha) THEN
               CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ia_trunc_arr, fm_B_BSE, &
                             homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
            ELSE
               CALL create_B(fm_mat_S_ia_trunc_arr, fm_mat_S_bar_ia_bse_arr, fm_B_BSE, &
                             homo_red_arr, virt_red_arr, dimen_RI, unit_nr, mp2_env)
            END IF
            CALL bse_matvec_selfcheck(mv_env, fm_A_BSE, unit_nr, fm_B_BSE)
         ELSE
            CALL bse_matvec_selfcheck(mv_env, fm_A_BSE, unit_nr)
         END IF
      END IF

      IF (do_tda) THEN
         CALL bse_davidson_tda(mv_env, mp2_env, unit_nr, exc_ens, fm_X)
         IF (do_debug) CALL bse_davidson_refcheck_tda(fm_A_BSE, exc_ens, fm_X, unit_nr)

         CALL copy_to_ia_context(fm_X, fm_mat_S_ia_trunc_arr(1), fm_X_ia)
         CALL postprocess_bse(exc_ens, fm_X_ia, mp2_env, qs_env, mo_coeff, &
                              homo_red_arr(1), virt_red_arr(1), homo(1), unit_nr, .TRUE.)
         CALL cp_fm_release(fm_X_ia)
         DEALLOCATE (exc_ens)
      END IF

      IF (do_abba) THEN
         ! the TDA vectors are the initial guess of a combined run
         SELECT CASE (mp2_env%bse%abba_solver)
         CASE (bse_abba_mk_davidson)
            solver_str = 'MK_DAVIDSON'
            IF (do_tda) THEN
               CALL bse_davidson_abba_mk(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                         abba_status, ab_margin, fm_X)
            ELSE
               CALL bse_davidson_abba_mk(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                         abba_status, ab_margin)
            END IF
         CASE (bse_abba_olsen_stratmann)
            solver_str = 'OLSEN_STRATMANN'
            IF (do_tda) THEN
               CALL bse_davidson_abba_os(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                         abba_status, ab_margin, fm_X)
            ELSE
               CALL bse_davidson_abba_os(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                         abba_status, ab_margin)
            END IF
         CASE DEFAULT
            CPABORT("Unknown BSE_ITERAT%ABBA_SOLVER")
         END SELECT
         IF (abba_status == abba_indefinite) THEN
            WRITE (margin_str, '(F16.6)') ab_margin*evolt
            IF (.NOT. do_tda) THEN
               CALL cp_abort(__LOCATION__, &
                             "Matrix (A-B) is not positive definite on the trial space of the "// &
                             TRIM(solver_str)//" solver (eigenvalue "//TRIM(ADJUSTL(margin_str))//" eV). "// &
                             "The full ABBA problem is ill-defined for this reference.")
            END IF
            ! the TDA results of a combined run stand on their own
            CALL cp_warn(__LOCATION__, &
                         "BSE: A-B is not positive definite on the trial space (margin "// &
                         TRIM(ADJUSTL(margin_str))//" eV); the ABBA calculation is skipped, "// &
                         "TDA results above are complete.")
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A4)') 'BSE|'
               WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.6)') 'BSE|', &
                  'A-B is not positive definite on the trial space (eV)', ab_margin*evolt
               WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
                  'The ABBA calculation is skipped, the TDA results above are complete'
            END IF
         ELSE
            IF (do_debug) THEN
               CALL bse_davidson_refcheck_abba(fm_A_BSE, fm_B_BSE, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                               mp2_env, unit_nr)
            END IF

            CALL copy_to_ia_context(fm_X_abba, fm_mat_S_ia_trunc_arr(1), fm_X_ia)
            CALL copy_to_ia_context(fm_Y_abba, fm_mat_S_ia_trunc_arr(1), fm_Y_ia)
            CALL cp_fm_release(fm_X_abba)
            CALL cp_fm_release(fm_Y_abba)
            CALL postprocess_bse(exc_ens_abba, fm_X_ia, mp2_env, qs_env, mo_coeff, &
                                 homo_red_arr(1), virt_red_arr(1), homo(1), unit_nr, .FALSE., fm_Y_ia)
            CALL cp_fm_release(fm_X_ia)
            CALL cp_fm_release(fm_Y_ia)
            DEALLOCATE (exc_ens_abba)
         END IF
      END IF

      IF (do_tda) CALL cp_fm_release(fm_X)
      IF (do_debug) THEN
         CALL cp_fm_release(fm_A_BSE)
         IF (do_abba) CALL cp_fm_release(fm_B_BSE)
      END IF
      CALL bse_matvec_release(mv_env)

      CALL timestop(handle)

   END SUBROUTINE solve_bse_iteratively

! **************************************************************************************************
!> \brief Copies eigenvectors from the process grid of the solver to the one of the RI slabs, from
!>        which postprocess_bse takes its grid
!> \param fm_vec ...
!> \param fm_S_ia slab that provides the target context
!> \param fm_vec_ia ...
! **************************************************************************************************
   SUBROUTINE copy_to_ia_context(fm_vec, fm_S_ia, fm_vec_ia)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_vec, fm_S_ia
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_vec_ia

      INTEGER                                            :: n_exc, n_ov
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL cp_fm_get_info(fm_vec, nrow_global=n_ov, ncol_global=n_exc)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=fm_S_ia%matrix_struct%para_env, &
                               context=fm_S_ia%matrix_struct%context, &
                               nrow_global=n_ov, ncol_global=n_exc)
      CALL cp_fm_create(fm_vec_ia, fm_struct, name="fm_davidson_ia")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm_submat_general(fm_vec, fm_vec_ia, n_ov, n_exc, 1, 1, 1, 1, &
                                      fm_S_ia%matrix_struct%context)

   END SUBROUTINE copy_to_ia_context

! **************************************************************************************************
!> \brief Sets up the matrix-free application of A and B with the prefactors and the screened
!>        slabs that create_A and create_B use for the chosen SPIN_CONFIG and SCREENING_IN_W
!> \param mv_env ...
!> \param fm_S_ia ...
!> \param fm_S_ij ...
!> \param fm_S_ab ...
!> \param fm_S_bar_ia ...
!> \param fm_S_bar_ij ...
!> \param Eigenval_reduced ...
!> \param homo_red ...
!> \param virt_red ...
!> \param do_abba ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_matvec_env(mv_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_S_bar_ia, fm_S_bar_ij, &
                                Eigenval_reduced, homo_red, virt_red, do_abba, mp2_env, unit_nr)

      TYPE(bse_matvec_env_type), INTENT(OUT)             :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_S_ia, fm_S_ij, fm_S_ab, fm_S_bar_ia, &
                                                            fm_S_bar_ij
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval_reduced
      INTEGER, INTENT(IN)                                :: homo_red, virt_red
      LOGICAL, INTENT(IN)                                :: do_abba
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      REAL(KIND=dp)                                      :: alpha, w_fac

      alpha = 0.0_dp
      SELECT CASE (mp2_env%bse%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      SELECT CASE (mp2_env%bse%screening_method)
      CASE (bse_screening_alpha)
         w_fac = mp2_env%bse%screening_factor
      CASE (bse_screening_rpa)
         w_fac = 0.0_dp
      CASE DEFAULT
         w_fac = 1.0_dp
      END SELECT

      ! TDHF and ALPHA screening contract the bare slabs, all others the W-multiplied ones
      IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
          mp2_env%bse%screening_method == bse_screening_alpha) THEN
         CALL bse_matvec_create(mv_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_S_ia, Eigenval_reduced, &
                                homo_red, virt_red, alpha, w_fac, do_abba, unit_nr)
      ELSE
         CALL bse_matvec_create(mv_env, fm_S_ia, fm_S_bar_ij, fm_S_ab, fm_S_bar_ia, Eigenval_reduced, &
                                homo_red, virt_red, alpha, w_fac, do_abba, unit_nr)
      END IF

   END SUBROUTINE create_matvec_env

END MODULE bse_main

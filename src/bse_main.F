!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
                                              diagonalize_C,&
                                              truncate_BSE_matrices
   USE bse_iterative,                   ONLY: do_subspace_iterations,&
                                              fill_local_3c_arrays
   USE bse_util,                        ONLY: deallocate_matrices_bse,&
                                              estimate_BSE_resources,&
                                              mult_B_with_W,&
                                              print_BSE_start_flag
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: group_dist_d1_type
   USE input_constants,                 ONLY: bse_fulldiag_abba,&
                                              bse_fulldiag_both,&
                                              bse_fulldiag_tda,&
                                              bse_iter_abba,&
                                              bse_iter_both,&
                                              bse_iter_tda
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_main'

   PUBLIC :: start_bse_calculation

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param gd_array ...
!> \param color_sub ...
!> \param para_env ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE start_bse_calculation(fm_mat_S, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
                                    Eigenval, homo, virtual, dimen_RI, dimen_RI_red, &
                                    gd_array, color_sub, &
                                    para_env, mp2_env, unit_nr)
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_mat_S
      TYPE(cp_fm_type)                                   :: fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                                            fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: Eigenval
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp2_type)                                     :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'start_bse_calculation'

      INTEGER                                            :: handle, homo_red, virtual_red
      LOGICAL                                            :: bse_abba, bse_full_diag, &
                                                            bse_iterat_diag, bse_tda
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_iaQ_bse_local
      TYPE(cp_fm_type) :: fm_mat_A_BSE, fm_mat_B_BSE, fm_mat_C_BSE, fm_mat_S_ab_trunc, &
         fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_ij_trunc, fm_mat_S_trunc, &
         fm_mat_work_diff_raised_neg_BSE, fm_mat_work_diff_raised_pos_BSE

      CALL timeset(routineN, handle)

      bse_full_diag = .FALSE.
      bse_iterat_diag = .FALSE.
      bse_tda = .FALSE.
      bse_abba = .FALSE.
      !Method: Iterative or full diagonalization
      SELECT CASE (mp2_env%ri_g0w0%bse_method)
      CASE (bse_fulldiag_tda, bse_fulldiag_abba, bse_fulldiag_both)
         bse_full_diag = .TRUE.
      CASE (bse_iter_tda, bse_iter_abba, bse_iter_both)
         bse_iterat_diag = .TRUE.
         !MG: Basics of the Davidson solver are implemented, but not rigorously checked.
         CPABORT("Iterative BSE not yet implemented")
      END SELECT
      !Approximation: TDA and/or full ABBA matrix
      SELECT CASE (mp2_env%ri_g0w0%bse_method)
      CASE (bse_fulldiag_tda, bse_iter_tda)
         bse_tda = .TRUE.
      CASE (bse_fulldiag_abba, bse_iter_abba)
         bse_abba = .TRUE.
      CASE (bse_iter_both, bse_fulldiag_both)
         bse_tda = .TRUE.
         bse_abba = .TRUE.
      END SELECT

      CALL print_BSE_start_flag(bse_tda, bse_abba, mp2_env, unit_nr)

      CALL fm_mat_S(1)%matrix_struct%para_env%sync()

      !Reduce matrices in case of energy cutoff for occupied and unoccupied in A/B-BSE-matrices
      CALL truncate_BSE_matrices(fm_mat_S(1), fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                 fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                 Eigenval(:, 1, 1), Eigenval_reduced, &
                                 homo(1), virtual(1), dimen_RI, unit_nr, &
                                 homo_red, virtual_red, &
                                 mp2_env)

      CALL mult_B_with_W(fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, fm_mat_S_trunc, &
                         fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                         fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
                         dimen_RI_red, homo_red, virtual_red)

      IF (bse_iterat_diag) THEN
         CALL fill_local_3c_arrays(fm_mat_S_ab_trunc, fm_mat_S_trunc, &
                                   fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                   B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                   B_iaQ_bse_local, dimen_RI_red, homo_red, virtual_red, &
                                   gd_array, color_sub, para_env)
      END IF

      IF (bse_full_diag) THEN
         CALL estimate_BSE_resources(homo_red, virtual_red, unit_nr, bse_abba, &
                                     para_env, diag_runtime_est)

         CALL create_A(fm_mat_S_trunc, fm_mat_S_bar_ij_bse, fm_mat_S_ab_trunc, &
                       fm_mat_A_BSE, Eigenval_reduced, unit_nr, &
                       homo_red, virtual_red, dimen_RI, mp2_env, &
                       para_env)
         IF (bse_abba) THEN
            CALL create_B(fm_mat_S_trunc, fm_mat_S_bar_ia_bse, fm_mat_B_BSE, &
                          homo_red, virtual_red, dimen_RI, unit_nr, mp2_env)
            CALL create_hermitian_form_of_ABBA(fm_mat_A_BSE, fm_mat_B_BSE, fm_mat_C_BSE, &
                                               fm_mat_work_diff_raised_pos_BSE, fm_mat_work_diff_raised_neg_BSE, &
                                               homo_red, virtual_red, unit_nr, mp2_env, diag_runtime_est)
         END IF
         CALL cp_fm_release(fm_mat_B_BSE)
         IF (bse_tda) THEN
            CALL diagonalize_A(fm_mat_A_BSE, homo_red, virtual_red, homo(1), virtual(1), &
                               unit_nr, diag_runtime_est, mp2_env)
         END IF
         ! Release to avoid faulty use of changed A matrix
         CALL cp_fm_release(fm_mat_A_BSE)
         IF (bse_abba) THEN
            CALL diagonalize_C(fm_mat_C_BSE, homo_red, virtual_red, homo(1), virtual(1), &
                               fm_mat_work_diff_raised_pos_BSE, fm_mat_work_diff_raised_neg_BSE, &
                               unit_nr, diag_runtime_est, mp2_env)
         END IF
         ! Release to avoid faulty use of changed C matrix
         CALL cp_fm_release(fm_mat_C_BSE)
      END IF

      CALL deallocate_matrices_bse(fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                   fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                   fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)
      IF (bse_iterat_diag) THEN
         CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo(1), virtual(1), mp2_env%ri_g0w0%bse_spin_config, unit_nr, &
                                     Eigenval(:, 1, 1), para_env, mp2_env)
         !Deallocate "three center" integrals
         DEALLOCATE (B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, B_iaQ_bse_local)
      END IF

      CALL timestop(handle)

   END SUBROUTINE start_bse_calculation

END MODULE bse_main

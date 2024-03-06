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
                                              diagonalize_C
   USE bse_iterative,                   ONLY: do_subspace_iterations,&
                                              fill_local_3c_arrays
   USE bse_util,                        ONLY: deallocate_matrices_bse,&
                                              estimate_BSE_resources,&
                                              mult_B_with_W,&
                                              print_BSE_start_flag,&
                                              truncate_BSE_matrices
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: group_dist_d1_type
   USE input_constants,                 ONLY: bse_abba,&
                                              bse_both,&
                                              bse_fulldiag,&
                                              bse_iterdiag,&
                                              bse_tda
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
!> \brief Main subroutine managing BSE calculations
!> \param fm_mat_S_ia_bse ...
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
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE start_bse_calculation(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
                                    Eigenval, homo, virtual, dimen_RI, dimen_RI_red, &
                                    gd_array, color_sub, mp2_env, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ia_bse, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: Eigenval
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(mp2_type)                                     :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'start_bse_calculation'

      INTEGER                                            :: handle, homo_red, virtual_red
      LOGICAL                                            :: my_do_abba, my_do_fulldiag, &
                                                            my_do_iterat_diag, my_do_tda
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_iaQ_bse_local
      TYPE(cp_fm_type) :: fm_A_BSE, fm_B_BSE, fm_C_BSE, fm_inv_sqrt_A_minus_B, fm_mat_S_ab_trunc, &
         fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_ia_trunc, fm_mat_S_ij_trunc, &
         fm_sqrt_A_minus_B
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      para_env => fm_mat_S_ia_bse%matrix_struct%para_env

      my_do_fulldiag = .FALSE.
      my_do_iterat_diag = .FALSE.
      my_do_tda = .FALSE.
      my_do_abba = .FALSE.
      !Method: Iterative or full diagonalization
      SELECT CASE (mp2_env%ri_g0w0%bse_diag_method)
      CASE (bse_iterdiag)
         my_do_iterat_diag = .TRUE.
         !MG: Basics of the Davidson solver are implemented, but not rigorously checked.
         CPABORT("Iterative BSE not yet implemented")
      CASE (bse_fulldiag)
         my_do_fulldiag = .TRUE.
      END SELECT
      !Approximation: TDA and/or full ABBA matrix
      SELECT CASE (mp2_env%ri_g0w0%bse_approx)
      CASE (bse_tda)
         my_do_tda = .TRUE.
      CASE (bse_abba)
         my_do_abba = .TRUE.
      CASE (bse_both)
         my_do_tda = .TRUE.
         my_do_abba = .TRUE.
      END SELECT

      CALL print_BSE_start_flag(my_do_tda, my_do_abba, mp2_env, unit_nr)

      CALL fm_mat_S_ia_bse%matrix_struct%para_env%sync()

      ! Reduce matrices in case of energy cutoff for occupied and unoccupied in A/B-BSE-matrices
      CALL truncate_BSE_matrices(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                 fm_mat_S_ia_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                 Eigenval(:, 1, 1), Eigenval_reduced, &
                                 homo(1), virtual(1), dimen_RI, unit_nr, &
                                 homo_red, virtual_red, &
                                 mp2_env)
      ! \bar{B}^P_rs = \sum_R W_PR B^R_rs where B^R_rs = \sum_T [1/sqrt(v)]_RT (T|rs)
      ! r,s: MO-index, P,R,T: RI-index
      ! B: fm_mat_S_..., W: fm_mat_Q_...
      CALL mult_B_with_W(fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, fm_mat_S_ia_trunc, &
                         fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                         fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
                         dimen_RI_red, homo_red, virtual_red)

      IF (my_do_iterat_diag) THEN
         CALL fill_local_3c_arrays(fm_mat_S_ab_trunc, fm_mat_S_ia_trunc, &
                                   fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                   B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                   B_iaQ_bse_local, dimen_RI_red, homo_red, virtual_red, &
                                   gd_array, color_sub, para_env)
      END IF

      IF (my_do_fulldiag) THEN
         ! Quick estimate of memory consumption and runtime of diagonalizations
         CALL estimate_BSE_resources(homo_red, virtual_red, unit_nr, my_do_abba, &
                                     para_env, diag_runtime_est)
         ! Matrix A constructed from GW energies and 3c-B-matrices (cf. subroutine mult_B_with_W)
         ! A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab
         ! ε_a, ε_i are GW singleparticle energies from Eigenval_reduced
         ! α is a spin-dependent factor
         ! v_ia,jb = \sum_P B^P_ia B^P_jb (unscreened Coulomb interaction)
         ! W_ij,ab = \sum_P \bar{B}^P_ij B^P_ab (screened Coulomb interaction)
         CALL create_A(fm_mat_S_ia_trunc, fm_mat_S_bar_ij_bse, fm_mat_S_ab_trunc, &
                       fm_A_BSE, Eigenval_reduced, unit_nr, &
                       homo_red, virtual_red, dimen_RI, mp2_env, &
                       para_env)
         IF (my_do_abba) THEN
            ! Matrix B constructed from 3c-B-matrices (cf. subroutine mult_B_with_W)
            ! B_ia,jb = α * v_ia,jb - W_ib,aj
            ! α is a spin-dependent factor
            ! v_ia,jb = \sum_P B^P_ia B^P_jb (unscreened Coulomb interaction)
            ! W_ib,aj = \sum_P \bar{B}^P_ib B^P_aj (screened Coulomb interaction)
            CALL create_B(fm_mat_S_ia_trunc, fm_mat_S_bar_ia_bse, fm_B_BSE, &
                          homo_red, virtual_red, dimen_RI, unit_nr, mp2_env)
            ! Construct Matrix C=(A-B)^0.5 (A+B) (A-B)^0.5 to solve full BSE matrix as a hermitian problem
            ! (cf. Eq. (A7) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001)).
            ! We keep fm_sqrt_A_minus_B and fm_inv_sqrt_A_minus_B for print of singleparticle transitions
            ! of ABBA as described in Eq. (A10) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001).
            CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                               fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                               homo_red, virtual_red, unit_nr, mp2_env, diag_runtime_est)
         END IF
         CALL cp_fm_release(fm_B_BSE)
         IF (my_do_tda) THEN
            ! Solving the hermitian eigenvalue equation A X^n = Ω^n X^n
            CALL diagonalize_A(fm_A_BSE, homo_red, virtual_red, homo(1), virtual(1), &
                               unit_nr, diag_runtime_est, mp2_env)
         END IF
         ! Release to avoid faulty use of changed A matrix
         CALL cp_fm_release(fm_A_BSE)
         IF (my_do_abba) THEN
            ! Solving eigenvalue equation C Z^n = (Ω^n)^2 Z^n .
            ! Here, the eigenvectors Z^n relate to X^n via
            ! Eq. (A10) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001).
            CALL diagonalize_C(fm_C_BSE, homo_red, virtual_red, homo(1), virtual(1), &
                               fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                               unit_nr, diag_runtime_est, mp2_env)
         END IF
         ! Release to avoid faulty use of changed C matrix
         CALL cp_fm_release(fm_C_BSE)
      END IF

      CALL deallocate_matrices_bse(fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                   fm_mat_S_ia_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                   fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)
      IF (my_do_iterat_diag) THEN
         ! Contains untested Block-Davidson algorithm
         CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo(1), virtual(1), mp2_env%ri_g0w0%bse_spin_config, unit_nr, &
                                     Eigenval(:, 1, 1), para_env, mp2_env)
         ! Deallocate local 3c-B-matrices
         DEALLOCATE (B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, B_iaQ_bse_local)
      END IF

      CALL timestop(handle)

   END SUBROUTINE start_bse_calculation

END MODULE bse_main

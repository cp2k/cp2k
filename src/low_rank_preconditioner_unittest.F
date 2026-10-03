!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM low_rank_preconditioner_unittest
   USE kinds,                           ONLY: dp
   USE low_rank_preconditioner_model,   ONLY: apply_low_rank_dense,&
                                              low_rank_inverse_weight,&
                                              low_rank_select_rank

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: nao = 5, nocc = 3, rank = 2
   REAL(KIND=dp), PARAMETER                            :: tolerance = 2.0E-14_dp
   INTEGER                                             :: i, selected_rank
   REAL(KIND=dp)                                       :: angle, energy_gap, mu, window
   REAL(KIND=dp), DIMENSION(7)                         :: full_eigenvalues
   REAL(KIND=dp), DIMENSION(rank)                      :: eigenvalues
   REAL(KIND=dp), DIMENSION(nao, nao)                  :: identity, overlap_inverse
   REAL(KIND=dp), DIMENSION(nao, rank)                 :: vectors
   REAL(KIND=dp), DIMENSION(nao, nocc)                 :: gradient, gradient_rotated, RESULT, &
                                                          result_rotated
   REAL(KIND=dp), DIMENSION(nocc, nocc)                :: rotation

   identity = 0.0_dp
   DO i = 1, nao
      identity(i, i) = 1.0_dp
   END DO

   overlap_inverse = identity
   overlap_inverse(1, 1) = 1.4_dp
   overlap_inverse(2, 2) = 0.8_dp
   overlap_inverse(1, 2) = 0.1_dp
   overlap_inverse(2, 1) = 0.1_dp

   vectors = 0.0_dp
   vectors(3, 1) = 1.0_dp
   vectors(4, 2) = 0.8_dp
   vectors(5, 2) = 0.6_dp
   eigenvalues = [0.35_dp, 0.8_dp]
   mu = -0.2_dp
   energy_gap = 0.08_dp
   window = 1.0_dp

   gradient = RESHAPE([0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.2_dp, &
                       0.6_dp, 0.5_dp, -0.3_dp, 0.1_dp, 0.7_dp, &
                       -0.4_dp, 0.8_dp, 0.2_dp, -0.5_dp, 0.9_dp], [nao, nocc])
   angle = 0.37_dp
   rotation = 0.0_dp
   rotation(1, 1) = COS(angle)
   rotation(1, 2) = -SIN(angle)
   rotation(2, 1) = SIN(angle)
   rotation(2, 2) = COS(angle)
   rotation(3, 3) = 1.0_dp
   gradient_rotated = MATMUL(gradient, rotation)

   CALL apply_low_rank_dense(overlap_inverse, vectors, eigenvalues, mu, energy_gap, window, &
                             gradient, RESULT)
   CALL apply_low_rank_dense(overlap_inverse, vectors, eigenvalues, mu, energy_gap, window, &
                             gradient_rotated, result_rotated)
   IF (MAXVAL(ABS(result_rotated - MATMUL(RESULT, rotation))) > tolerance) THEN
      ERROR STOP "Low-rank preconditioner is not right-covariant under occupied rotations"
   END IF

   full_eigenvalues = [-0.8_dp, -0.2_dp, 0.2_dp, 0.4_dp, 0.7_dp, 0.7_dp, 1.4_dp]
   selected_rank = low_rank_select_rank(full_eigenvalues, 2, 3, 1.0E-8_dp)
   IF (selected_rank /= 2) ERROR STOP "Rank cap split a degenerate complementary manifold"
   selected_rank = low_rank_select_rank(full_eigenvalues, 2, 5, 1.0E-8_dp)
   IF (selected_rank /= 5) ERROR STOP "Full complementary space did not retain its rank"
   IF (ABS(low_rank_inverse_weight(1.5_dp, 0.0_dp, 0.08_dp) - 2.0_dp/3.0_dp) > tolerance) THEN
      ERROR STOP "Common-reference inverse weight differs from the spectral model"
   END IF

END PROGRAM low_rank_preconditioner_unittest

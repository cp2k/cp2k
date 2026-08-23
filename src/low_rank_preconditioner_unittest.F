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
   REAL(KIND=dp)                                       :: angle, common_shift, energy_gap, mu, window
   REAL(KIND=dp), DIMENSION(7)                         :: full_eigenvalues
   REAL(KIND=dp), DIMENSION(rank)                      :: eigenvalues
   REAL(KIND=dp), DIMENSION(nao, nao)                  :: complement_operator, &
                                                          complement_operator_rotated, hamiltonian, &
                                                          identity, overlap_inverse, projector, &
                                                          projector_rotated
   REAL(KIND=dp), DIMENSION(nao, 2)                    :: occupied, occupied_rotated, vectors
   REAL(KIND=dp), DIMENSION(nao, nocc)                 :: gradient, gradient_rotated, RESULT, &
                                                          result_rotated
   REAL(KIND=dp), DIMENSION(nocc, nocc)                :: rotation
   REAL(KIND=dp), DIMENSION(2, 2)                      :: occupied_rotation

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

   ! Verify that the projected complementary operator depends only on the occupied subspace.
   occupied = 0.0_dp
   occupied(1, 1) = 1.0_dp/SQRT(2.0_dp)
   occupied(2, 1) = 1.0_dp/SQRT(2.0_dp)
   occupied(3, 2) = 1.0_dp/SQRT(2.0_dp)
   occupied(4, 2) = 1.0_dp/SQRT(2.0_dp)
   occupied_rotation = RESHAPE([COS(angle), SIN(angle), -SIN(angle), COS(angle)], [2, 2])
   occupied_rotated = MATMUL(occupied, occupied_rotation)

   hamiltonian = RESHAPE([1.0_dp, 0.2_dp, 0.1_dp, 0.0_dp, 0.3_dp, &
                          0.2_dp, 1.4_dp, 0.0_dp, 0.2_dp, 0.1_dp, &
                          0.1_dp, 0.0_dp, 2.0_dp, 0.4_dp, 0.2_dp, &
                          0.0_dp, 0.2_dp, 0.4_dp, 2.3_dp, 0.1_dp, &
                          0.3_dp, 0.1_dp, 0.2_dp, 0.1_dp, 3.0_dp], [nao, nao])
   projector = identity - MATMUL(occupied, TRANSPOSE(occupied))
   projector_rotated = identity - MATMUL(occupied_rotated, TRANSPOSE(occupied_rotated))
   common_shift = -10.2_dp
   complement_operator = MATMUL(TRANSPOSE(projector), MATMUL(hamiltonian, projector)) + &
                         common_shift*MATMUL(occupied, TRANSPOSE(occupied))
   complement_operator_rotated = &
      MATMUL(TRANSPOSE(projector_rotated), MATMUL(hamiltonian, projector_rotated)) + &
      common_shift*MATMUL(occupied_rotated, TRANSPOSE(occupied_rotated))
   IF (MAXVAL(ABS(complement_operator_rotated - complement_operator)) > tolerance) THEN
      ERROR STOP "Complementary-state construction depends on the occupied orbital gauge"
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

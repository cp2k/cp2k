!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM optbas_frontier_utils_unittest
   USE kinds,                           ONLY: dp
   USE optbas_frontier_utils,           ONLY: closed_shell_nocc,&
                                              compute_frontier_losses,&
                                              frontier_result_type,&
                                              smooth_frontier_weight

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: tolerance = 1.0E-12_dp
   REAL(KIND=dp), DIMENSION(4) :: candidate_eigenvalues, reference_eigenvalues
   REAL(KIND=dp), DIMENSION(2) :: coefficients, coefficients_initial
   REAL(KIND=dp), DIMENSION(4, 4) :: overlap
   TYPE(frontier_result_type) :: result, result_rotated

   reference_eigenvalues = [-1.0_dp, -0.5_dp, 0.1_dp, 0.2_dp]
   candidate_eigenvalues = reference_eigenvalues
   coefficients_initial = [1.0_dp, 0.0_dp]
   coefficients = coefficients_initial

   IF (closed_shell_nocc(56) /= 28) ERROR STOP "restricted occupied-orbital count is incorrect"

   ! Identical spaces have zero loss. A permutation/sign change inside equally weighted subspaces
   ! leaves the trace objectives invariant.
   overlap = 0.0_dp
   overlap(1, 1) = 1.0_dp
   overlap(2, 2) = 1.0_dp
   overlap(3, 3) = 1.0_dp
   overlap(4, 4) = 1.0_dp
   CALL compute_frontier_losses(overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
                                10.0_dp, 0.1_dp, 0.2_dp, coefficients, coefficients_initial, RESULT)
   CALL assert_close(result%loss_occupied, 0.0_dp, "identical occupied subspace")
   CALL assert_close(result%loss_virtual, 0.0_dp, "identical virtual subspace")
   CALL assert_close(result%loss_empty_overlap, 0.0_dp, "identical empty completeness")
   CALL assert_close(result%virtual_overlap_percent, 100.0_dp, "identical empty percentage")

   overlap = 0.0_dp
   overlap(1, 2) = -1.0_dp
   overlap(2, 1) = 1.0_dp
   overlap(3, 4) = -1.0_dp
   overlap(4, 3) = 1.0_dp
   CALL compute_frontier_losses(overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
                                10.0_dp, 0.1_dp, 0.2_dp, coefficients, coefficients_initial, result_rotated)
   CALL assert_close(result_rotated%loss_occupied, result%loss_occupied, "occupied permutation/sign invariance")
   CALL assert_close(result_rotated%loss_virtual, result%loss_virtual, "virtual permutation/sign invariance")
   CALL assert_close(result_rotated%loss_empty_overlap, result%loss_empty_overlap, &
                     "empty completeness permutation/sign invariance")

   ! Orthogonal equal-rank/equally weighted spaces give unit normalized losses.
   overlap = 0.0_dp
   candidate_eigenvalues = reference_eigenvalues
   CALL compute_frontier_losses(overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
                                10.0_dp, 0.1_dp, 0.2_dp, coefficients, coefficients_initial, RESULT)
   CALL assert_close(result%loss_occupied, 1.0_dp, "orthogonal occupied normalization")
   CALL assert_close(result%loss_virtual, 1.0_dp, "orthogonal virtual normalization")
   CALL assert_close(result%loss_empty_overlap, 1.0_dp, "orthogonal empty completeness")

   ! A gap error of E_scale and a coefficient displacement of ||c0|| both give unit loss.
   overlap = 0.0_dp
   candidate_eigenvalues(3:4) = candidate_eigenvalues(3:4) + 0.2_dp
   coefficients = [2.0_dp, 0.0_dp]
   CALL compute_frontier_losses(overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
                                10.0_dp, 0.1_dp, 0.2_dp, coefficients, coefficients_initial, RESULT)
   CALL assert_close(result%loss_gap, 1.0_dp, "gap normalization")
   CALL assert_close(result%loss_coefficients, 1.0_dp, "coefficient normalization")

   ! The frontier weight is continuous and monotonic across the cutoff.
   IF (.NOT. smooth_frontier_weight(0.999_dp, 1.0_dp, 0.1_dp) > &
       smooth_frontier_weight(1.001_dp, 1.0_dp, 0.1_dp)) ERROR STOP "frontier weight is not monotonic"
   CALL assert_close(smooth_frontier_weight(1.0_dp, 1.0_dp, 0.1_dp), 0.5_dp, "cutoff midpoint")

CONTAINS

   SUBROUTINE assert_close(actual, expected, description)
      REAL(KIND=dp), INTENT(IN)                          :: actual, expected
      CHARACTER(LEN=*), INTENT(IN)                       :: description

      IF (ABS(actual - expected) > tolerance) THEN
         WRITE (*, '(A,2ES24.14)') TRIM(description)//": ", actual, expected
         ERROR STOP "frontier objective unit test failed"
      END IF

   END SUBROUTINE assert_close

END PROGRAM optbas_frontier_utils_unittest

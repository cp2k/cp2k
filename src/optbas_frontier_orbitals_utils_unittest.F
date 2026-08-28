!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM optbas_frontier_orbitals_utils_unittest
   USE kinds,                           ONLY: dp
   USE optbas_frontier_orbitals_utils,  ONLY: angular_momentum_diagnostics_type,&
                                              compute_angular_momentum_diagnostics,&
                                              compute_frontier_orbital_loss_function,&
                                              frontier_orbitals_result_type,&
                                              smooth_frontier_orbital_weight

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: tolerance = 1.0E-12_dp
   REAL(KIND=dp), DIMENSION(4) :: candidate_eigenvalues, reference_eigenvalues
   REAL(KIND=dp), DIMENSION(4, 4) :: orbital_overlap
   REAL(KIND=dp), DIMENSION(2, 2) :: lowdin_coefficients, orthogonal_ao_overlap
   INTEGER, DIMENSION(2) :: angular_momentum
   TYPE(angular_momentum_diagnostics_type) :: angular_diagnostics
   TYPE(frontier_orbitals_result_type) :: result, result_rotated

   reference_eigenvalues = [-1.0_dp, -0.5_dp, 0.1_dp, 0.2_dp]
   candidate_eigenvalues = reference_eigenvalues

   ! Identical spaces have zero loss. A permutation/sign change inside equally weighted subspaces
   ! leaves the trace objectives invariant.
   orbital_overlap = 0.0_dp
   orbital_overlap(1, 1) = 1.0_dp
   orbital_overlap(2, 2) = 1.0_dp
   orbital_overlap(3, 3) = 1.0_dp
   orbital_overlap(4, 4) = 1.0_dp
   CALL compute_frontier_orbital_loss_function( &
      orbital_overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
      10.0_dp, 0.1_dp, 0.2_dp, RESULT)
   CALL assert_close(result%loss_occupied_subspace, 0.0_dp, "identical occupied subspace")
   CALL assert_close(result%loss_virtual_subspace, 0.0_dp, "identical virtual subspace")
   CALL assert_close(result%empty_subspace_overlap, 1.0_dp, "identical empty subspace")
   CALL assert_close(result%loss_empty_subspace, 0.0_dp, "identical empty-space loss")
   CALL assert_close(result%virtual_subspace_similarity, 1.0_dp, &
                     "identical virtual-space similarity")

   orbital_overlap = 0.0_dp
   orbital_overlap(1, 2) = -1.0_dp
   orbital_overlap(2, 1) = 1.0_dp
   orbital_overlap(3, 4) = -1.0_dp
   orbital_overlap(4, 3) = 1.0_dp
   CALL compute_frontier_orbital_loss_function( &
      orbital_overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
      10.0_dp, 0.1_dp, 0.2_dp, result_rotated)
   CALL assert_close(result_rotated%loss_occupied_subspace, result%loss_occupied_subspace, &
                     "occupied permutation/sign invariance")
   CALL assert_close(result_rotated%loss_virtual_subspace, result%loss_virtual_subspace, &
                     "virtual permutation/sign invariance")

   ! Orthogonal equal-rank and equally weighted spaces give unit normalized losses.
   orbital_overlap = 0.0_dp
   candidate_eigenvalues = reference_eigenvalues
   CALL compute_frontier_orbital_loss_function( &
      orbital_overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
      10.0_dp, 0.1_dp, 0.2_dp, RESULT)
   CALL assert_close(result%loss_occupied_subspace, 1.0_dp, "orthogonal occupied normalization")
   CALL assert_close(result%loss_virtual_subspace, 1.0_dp, "orthogonal virtual normalization")
   CALL assert_close(result%empty_subspace_overlap, 0.0_dp, "orthogonal empty subspace")
   CALL assert_close(result%loss_empty_subspace, 1.0_dp, "orthogonal empty-space loss")
   CALL assert_close(result%virtual_subspace_similarity, 0.0_dp, &
                     "orthogonal virtual-space similarity")

   ! A gap error equal to E_scale gives a unit gap loss.
   candidate_eigenvalues(3:4) = candidate_eigenvalues(3:4) + 0.2_dp
   CALL compute_frontier_orbital_loss_function( &
      orbital_overlap, reference_eigenvalues, candidate_eigenvalues, 2, &
      10.0_dp, 0.1_dp, 0.2_dp, RESULT)
   CALL assert_close(result%loss_gap, 1.0_dp, "gap normalization")

   ! The virtual-orbital weight is continuous and monotonic across the cutoff.
   IF (.NOT. smooth_frontier_orbital_weight(0.999_dp, 1.0_dp, 0.1_dp) > &
       smooth_frontier_orbital_weight(1.001_dp, 1.0_dp, 0.1_dp)) THEN
      ERROR STOP "virtual-orbital weight is not monotonic"
   END IF
   CALL assert_close(smooth_frontier_orbital_weight(1.0_dp, 1.0_dp, 0.1_dp), &
                     0.5_dp, "cutoff midpoint")

   ! With orthonormal and identical AO bases, the occupied s and virtual p components
   ! both have unit reference contribution and unit corresponding-channel coverage.
   lowdin_coefficients = 0.0_dp
   lowdin_coefficients(1, 1) = 1.0_dp
   lowdin_coefficients(2, 2) = 1.0_dp
   orthogonal_ao_overlap = lowdin_coefficients
   angular_momentum = [0, 1]
   CALL compute_angular_momentum_diagnostics( &
      lowdin_coefficients, orthogonal_ao_overlap, angular_momentum, angular_momentum, &
      [-1.0_dp, 0.1_dp], 1, 10.0_dp, 0.1_dp, angular_diagnostics)
   CALL assert_close(angular_diagnostics%occupied_contribution(0), 1.0_dp, &
                     "occupied s contribution")
   CALL assert_close(angular_diagnostics%virtual_contribution(1), 1.0_dp, &
                     "virtual p contribution")
   CALL assert_close(angular_diagnostics%occupied_coverage(0), 1.0_dp, &
                     "occupied s coverage")
   CALL assert_close(angular_diagnostics%virtual_coverage(1), 1.0_dp, &
                     "virtual p coverage")

CONTAINS

   SUBROUTINE assert_close(actual, expected, description)
      REAL(KIND=dp), INTENT(IN)                          :: actual, expected
      CHARACTER(LEN=*), INTENT(IN)                       :: description

      IF (ABS(actual - expected) > tolerance) THEN
         WRITE (*, '(A,2ES24.14)') TRIM(description)//": ", actual, expected
         ERROR STOP "frontier-orbital objective unit test failed"
      END IF

   END SUBROUTINE assert_close

END PROGRAM optbas_frontier_orbitals_utils_unittest

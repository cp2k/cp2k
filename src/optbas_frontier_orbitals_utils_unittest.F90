!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Unit tests for the frontier-orbital objective
!> \author Jan Wilhelm
! **************************************************************************************************
PROGRAM optbas_frontier_orbitals_utils_unittest
   USE kinds,                           ONLY: dp
   USE optbas_frontier_orbitals_utils,  ONLY: compute_frontier_orbital_loss_function,&
                                              frontier_orbitals_result_type,&
                                              smooth_frontier_orbital_weight

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: tolerance = 1.0E-12_dp
   REAL(KIND=dp), PARAMETER :: eigenvalues(4) = [-1.0_dp, -0.5_dp, 0.1_dp, 0.2_dp]
   INTEGER :: i, j, orientation
   REAL(KIND=dp) :: candidate(4), overlap(4, 4), scale, losses(6), expected(6)
   TYPE(frontier_orbitals_result_type) :: result

   ! Equal virtual weights give L_occ=L_vir=L_empty=1-scale**2.
   DO i = 0, 2
      scale = 0.5_dp*REAL(i, KIND=dp)
      candidate = eigenvalues
      candidate(3:4) = candidate(3:4) + 0.2_dp*scale
      overlap = 0.0_dp
      DO j = 1, 4
         overlap(j, j) = scale
      END DO
      expected(1:3) = 1.0_dp - scale**2
      expected(4:6) = scale**2
      DO orientation = 1, 2
         CALL compute_frontier_orbital_loss_function(overlap, eigenvalues, candidate, 2, &
                                                     10.0_dp, 0.1_dp, 0.2_dp, RESULT)
         losses = [result%loss_occupied_subspace, result%loss_virtual_subspace, &
                   result%loss_empty_subspace, result%empty_subspace_overlap, &
                   result%virtual_subspace_similarity, result%loss_gap]
         IF (.NOT. ALL(ABS(losses - expected) <= tolerance)) ERROR STOP "Incorrect frontier-orbital objective"
         ! Permutations and signs within each subspace leave the objective unchanged.
         overlap = overlap(:, [2, 1, 4, 3])
         overlap(:, 1) = -overlap(:, 1)
      END DO
   END DO
   IF (.NOT. (smooth_frontier_orbital_weight(0.999_dp, 1.0_dp, 0.1_dp) > &
              smooth_frontier_orbital_weight(1.001_dp, 1.0_dp, 0.1_dp))) ERROR STOP "Nonmonotonic virtual weight"
   IF (.NOT. (ABS(smooth_frontier_orbital_weight(1.0_dp, 1.0_dp, 0.1_dp) - 0.5_dp) <= tolerance)) THEN
      ERROR STOP "Incorrect virtual-weight cutoff midpoint"
   END IF

END PROGRAM optbas_frontier_orbitals_utils_unittest

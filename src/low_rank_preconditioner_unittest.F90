!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM low_rank_preconditioner_unittest
   USE kinds,                           ONLY: dp
   USE low_rank_preconditioner_model,   ONLY: low_rank_inverse_weight,&
                                              low_rank_select_rank

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER                            :: tolerance = 2.0E-14_dp
   INTEGER                                             :: selected_rank
   REAL(KIND=dp), DIMENSION(7)                         :: full_eigenvalues

   full_eigenvalues = [-0.8_dp, -0.2_dp, 0.2_dp, 0.4_dp, 0.7_dp, 0.7_dp, 1.4_dp]
   selected_rank = low_rank_select_rank(full_eigenvalues, 2, 3, 1.0E-8_dp)
   IF (selected_rank /= 2) ERROR STOP "Rank cap split a degenerate complementary manifold"
   selected_rank = low_rank_select_rank(full_eigenvalues, 2, 5, 1.0E-8_dp)
   IF (selected_rank /= 5) ERROR STOP "Full complementary space did not retain its rank"
   IF (ABS(low_rank_inverse_weight(1.5_dp, 0.0_dp, 0.08_dp) - 2.0_dp/3.0_dp) > tolerance) THEN
      ERROR STOP "Common-reference inverse weight differs from the spectral model"
   END IF

END PROGRAM low_rank_preconditioner_unittest

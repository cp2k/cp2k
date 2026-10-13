!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM hfx_pair_list_methods_unittest
   USE hfx_pair_list_methods,           ONLY: min_image_gap_sq
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat

   ! Orthorhombic cell with two reversed cell vectors (positive determinant):
   ! moving r = (4,0,0) by the image (-10,0,0) gives |r'|**2 = 36, i.e. a gap of 20.
   hmat = 0.0_dp
   hmat(1, 1) = -10.0_dp
   hmat(2, 2) = -12.0_dp
   hmat(3, 3) = 14.0_dp
   IF (ABS(min_image_gap_sq([4.0_dp, 0.0_dp, 0.0_dp], hmat) - 20.0_dp) > 1.0E-12_dp) THEN
      ERROR STOP "Wrong minimum-image gap for reversed cell vectors"
   END IF

   CALL check_against_all_images([10.0_dp, 12.0_dp, 14.0_dp])
   CALL check_against_all_images([-10.0_dp, -12.0_dp, 14.0_dp])
   CALL check_against_all_images([-10.0_dp, 12.0_dp, 14.0_dp])
   CALL check_against_all_images([10.0_dp, -12.0_dp, -14.0_dp])

CONTAINS

! **************************************************************************************************
!> \brief Compare min_image_gap_sq with an explicit search over neighbouring images.
!> \param diag diagonal entries of an orthorhombic cell matrix
! **************************************************************************************************
   SUBROUTINE check_against_all_images(diag)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: diag

      INTEGER                                            :: i, j, k, n1, n2, n3
      REAL(KIND=dp)                                      :: exact, gap
      REAL(KIND=dp), DIMENSION(3)                        :: r, shifted
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell

      cell = 0.0_dp
      cell(1, 1) = diag(1)
      cell(2, 2) = diag(2)
      cell(3, 3) = diag(3)

      ! minimum-image vectors on a grid of fractional coordinates in [-1/2, 1/2]
      DO i = -8, 8
         DO j = -8, 8
            DO k = -8, 8
               r = [REAL(i, dp)*ABS(diag(1)), REAL(j, dp)*ABS(diag(2)), REAL(k, dp)*ABS(diag(3))]/16.0_dp
               gap = min_image_gap_sq(r, cell)
               exact = HUGE(1.0_dp)
               DO n1 = -2, 2
                  DO n2 = -2, 2
                     DO n3 = -2, 2
                        IF (n1 == 0 .AND. n2 == 0 .AND. n3 == 0) CYCLE
                        shifted = r - MATMUL(cell, REAL([n1, n2, n3], dp))
                        exact = MIN(exact, SUM(shifted**2) - SUM(r**2))
                     END DO
                  END DO
               END DO
               IF (ABS(gap - exact) > 1.0E-10_dp) THEN
                  ERROR STOP "min_image_gap_sq differs from the explicit image search"
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE check_against_all_images

END PROGRAM hfx_pair_list_methods_unittest

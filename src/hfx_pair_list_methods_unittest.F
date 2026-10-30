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

   REAL(KIND=dp), PARAMETER :: cells(3, 4) = RESHAPE([ &
                                                     10.0_dp, 12.0_dp, 14.0_dp, -10.0_dp, -12.0_dp, 14.0_dp, &
                                                     -10.0_dp, 12.0_dp, 14.0_dp, 10.0_dp, -12.0_dp, -14.0_dp], [3, 4])
   ! Centre, axial/mixed displacements, opposite corners and a point just outside the cell.
   REAL(KIND=dp), PARAMETER :: points(3, 8) = RESHAPE([ &
                                                      0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, &
                                                      0.0_dp, 0.0_dp, -6.0_dp, 4.0_dp, -5.0_dp, 6.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, &
                                                      -5.0_dp, -6.0_dp, -7.0_dp, 5.1_dp, 0.0_dp, 0.0_dp], [3, 8])
   INTEGER :: axis, icell, ipoint, ix, iy, iz
   REAL(KIND=dp) :: exact, hmat(3, 3), r(3), shifted(3)

   DO icell = 1, SIZE(cells, 2)
      hmat = 0.0_dp
      DO axis = 1, 3
         hmat(axis, axis) = cells(axis, icell)
      END DO
      DO ipoint = 1, SIZE(points, 2)
         r = points(:, ipoint)
         exact = HUGE(1.0_dp)
         DO ix = -1, 1
            DO iy = -1, 1
               DO iz = -1, 1
                  IF (ix == 0 .AND. iy == 0 .AND. iz == 0) CYCLE
                  shifted = r - MATMUL(hmat, REAL([ix, iy, iz], KIND=dp))
                  exact = MIN(exact, SUM(shifted**2) - SUM(r**2))
               END DO
            END DO
         END DO
         IF (.NOT. (ABS(min_image_gap_sq(r, hmat) - exact) <= 1.0E-10_dp)) THEN
            ERROR STOP "Minimum-image gap differs from the explicit image search"
         END IF
      END DO
   END DO

END PROGRAM hfx_pair_list_methods_unittest

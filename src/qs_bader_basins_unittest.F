!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_bader_basins_test
   USE kinds,                           ONLY: dp
   USE qs_bader_basins,                 ONLY: bader_integrate_grid

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: tol = 1.0E-12_dp, SPACING(3) = 1.0_dp
   INTEGER :: n_other
   LOGICAL :: assigned(2)
   REAL(KIND=dp) :: atomic(2, 2), other(2), vacuum(2)

   ! Point four flows into both basins, in the ratio 1:3.
   CALL bader_integrate_grid([0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 4.0_dp, 0.0_dp], &
                             [6, 1, 1], spacing, [.FALSE., .FALSE., .FALSE.], 0.0_dp, [3, 5], &
                             atomic, assigned, other, vacuum, n_other)
   IF (.NOT. ALL(assigned) .OR. n_other /= 0) ERROR STOP "Incorrect basin assignments"
   IF (MAXVAL(ABS(atomic(:, 1) - [3.25_dp, 2.25_dp])) > tol) ERROR STOP "Incorrect left basin"
   IF (MAXVAL(ABS(atomic(:, 2) - [4.75_dp, 1.75_dp])) > tol) ERROR STOP "Incorrect right basin"
   IF (MAXVAL(ABS(vacuum - [0.0_dp, 2.0_dp])) > tol) ERROR STOP "Incorrect vacuum population"

   ! When the last point wraps to the first, the two seed points share a single basin.
   CALL bader_integrate_grid([4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other)
   IF (ANY(assigned) .OR. n_other /= 1) ERROR STOP "A shared basin was not marked as ambiguous"
   IF (MAXVAL(ABS(other - [10.0_dp, 4.0_dp])) > tol) ERROR STOP "Incorrect periodic basin population"

   ! In nearest-atom mode, a maximum equidistant across the periodic boundary is shared.
   CALL bader_integrate_grid([4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other, &
                             atom_grid=RESHAPE([0.5_dp, 0.0_dp, 0.0_dp, 3.5_dp, 0.0_dp, 0.0_dp], [3, 2]))
   IF (.NOT. ALL(assigned) .OR. n_other /= 0) ERROR STOP "An equidistant basin was not assigned"
   IF (MAXVAL(ABS(atomic(:, 1) - [5.0_dp, 2.0_dp])) > tol) ERROR STOP "Incorrect shared basin"
   IF (MAXVAL(ABS(atomic(:, 2) - atomic(:, 1))) > tol) ERROR STOP "Basin sharing broke symmetry"
END PROGRAM qs_bader_basins_test

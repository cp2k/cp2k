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
   INTEGER :: labels(6), n_other, n_extra
   LOGICAL :: assigned(2)
   REAL(KIND=dp) :: atomic(2, 2), extra(2), other(2), vacuum(2)

   ! Point four flows into both basins, in the ratio 1:3.
   CALL bader_integrate_grid([0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 4.0_dp, 0.0_dp], &
                             [6, 1, 1], spacing, [.FALSE., .FALSE., .FALSE.], 0.0_dp, [3, 5], &
                             atomic, assigned, other, vacuum, n_other, labels=labels)
   IF (.NOT. ALL(assigned) .OR. n_other /= 0) ERROR STOP "Incorrect basin assignments"
   IF (MAXVAL(ABS(atomic(:, 1) - [3.25_dp, 2.25_dp])) > tol) ERROR STOP "Incorrect left basin"
   IF (MAXVAL(ABS(atomic(:, 2) - [4.75_dp, 1.75_dp])) > tol) ERROR STOP "Incorrect right basin"
   IF (MAXVAL(ABS(vacuum - [0.0_dp, 2.0_dp])) > tol) ERROR STOP "Incorrect vacuum population"
   IF (ANY(labels /= [0, 1, 1, 2, 2, 0])) ERROR STOP "Incorrect strongest-uphill basin labels"

   ! When the last point wraps to the first, the two seed points share a single basin.
   CALL bader_integrate_grid([4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other, labels=labels(:4))
   IF (ANY(assigned) .OR. n_other /= 1) ERROR STOP "A shared basin was not marked as ambiguous"
   IF (MAXVAL(ABS(other - [10.0_dp, 4.0_dp])) > tol) ERROR STOP "Incorrect periodic basin population"
   IF (ANY(labels(:4) /= -1)) ERROR STOP "A basin with two seeds has no unique label"

   ! In nearest-atom mode, a maximum equidistant across the periodic boundary is shared.
   CALL bader_integrate_grid([4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other, &
                             atom_grid=RESHAPE([0.5_dp, 0.0_dp, 0.0_dp, 3.5_dp, 0.0_dp, 0.0_dp], [3, 2]), &
                             n_extra=n_extra, extra=extra, labels=labels(:4))
   IF (.NOT. ALL(assigned) .OR. n_other /= 0) ERROR STOP "An equidistant basin was not assigned"
   IF (MAXVAL(ABS(atomic(:, 1) - [5.0_dp, 2.0_dp])) > tol) ERROR STOP "Incorrect shared basin"
   IF (MAXVAL(ABS(atomic(:, 2) - atomic(:, 1))) > tol) ERROR STOP "Basin sharing broke symmetry"
   IF (n_extra /= 1) ERROR STOP "Grouping hid the maximum without a unique nuclear seed"
   IF (MAXVAL(ABS(extra - [10.0_dp, 4.0_dp])) > tol) ERROR STOP "Incorrect population before grouping"
   IF (ANY(labels(:4) /= -1)) ERROR STOP "A shared nearest-atom basin has no unique label"

   ! A maximum plateau straddling a periodic boundary belongs to one seed.
   CALL bader_integrate_grid([4.0_dp, 1.0_dp, 0.0_dp, 4.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [4], &
                             atomic(:, :1), assigned(:1), other, vacuum, n_other, labels=labels(:4))
   IF (.NOT. assigned(1) .OR. n_other /= 0) ERROR STOP "Periodic maximum plateau was split"
   IF (MAXVAL(ABS(atomic(:, 1) - [9.0_dp, 3.0_dp])) > tol) ERROR STOP "Wrong maximum plateau population"
   IF (ANY(labels(:4) /= [1, 1, 0, 1])) ERROR STOP "Periodic plateau labels are inconsistent"

   CALL bader_integrate_grid([4.0_dp, 1.0_dp, 0.0_dp, 4.0_dp], &
                             [4, 1, 1], spacing, [.TRUE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other, &
                             atom_grid=RESHAPE([0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp], [3, 2]))
   IF (.NOT. ALL(assigned)) ERROR STOP "Only one end of the maximum plateau was considered"
   IF (MAXVAL(ABS(atomic(:, 2) - atomic(:, 1))) > tol) ERROR STOP "Plateau assignment depends on the representative"

   ! A nonterminal plateau drains collectively through both of its higher faces.
   CALL bader_integrate_grid([4.0_dp, 2.0_dp, 2.0_dp, 5.0_dp], &
                             [4, 1, 1], spacing, [.FALSE., .FALSE., .FALSE.], 0.0_dp, [1, 4], &
                             atomic, assigned, other, vacuum, n_other)
   IF (.NOT. ALL(assigned) .OR. n_other /= 0) ERROR STOP "Nonterminal plateau has no basin"
   IF (MAXVAL(ABS(atomic(:, 1) - [5.6_dp, 1.8_dp])) > tol) ERROR STOP "Wrong plateau flux to the left"
   IF (MAXVAL(ABS(atomic(:, 2) - [7.4_dp, 2.2_dp])) > tol) ERROR STOP "Wrong plateau flux to the right"
END PROGRAM qs_bader_basins_test

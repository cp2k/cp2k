!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local Wilson-loop linear algebra, independent of the electronic-structure representation.
!>        WCC are arg(eigenvalue)/(2*pi), as in Z2Pack. The occupied subspace must be isolated.
! **************************************************************************************************
MODULE topology_wilson
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wilson_step, wilson_spectrum, z2_from_wcc, wcc_distance, surface_resolved
   REAL(KIND=dp), PARAMETER :: two_pi = 6.283185307179586476925286766559_dp
   ! Reduced-coordinate tolerance for an ambiguous largest-gap crossing.
   REAL(KIND=dp), PARAMETER :: crossing_tol = 1.e-12_dp
   ! Fraction of the adjacent largest gap allowed for motion and gap separation.
   REAL(KIND=dp), PARAMETER :: gap_fraction = 0.3_dp
CONTAINS

! **************************************************************************************************
!> \brief Conservative neighbouring-line movement and largest-gap separation checks.
!> \param wcc centres, band index first and surface-line index second
!> \return whether the sampled surface is locally resolved (not a proof between samples)
! **************************************************************************************************
   FUNCTION surface_resolved(wcc) RESULT(resolved)
      REAL(KIND=dp), INTENT(IN)                          :: wcc(:, :)
      LOGICAL                                            :: resolved

      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: delta, left_gap, left_size, right_gap, &
                                                            right_size
      REAL(KIND=dp), ALLOCATABLE                         :: sorted(:, :)

      resolved = .FALSE.
      n = SIZE(wcc, 1)
      IF (n < 1 .OR. SIZE(wcc, 2) < 2) RETURN
      ALLOCATE (sorted(n, SIZE(wcc, 2)))
      sorted(:, :) = MODULO(wcc, 1.0_dp)
      DO i = 1, SIZE(wcc, 2)
         CALL sort_wcc(sorted(:, i))
      END DO
      DO i = 2, SIZE(wcc, 2)
         CALL largest_gap(sorted(:, i - 1), left_gap, left_size)
         CALL largest_gap(sorted(:, i), right_gap, right_size)
         IF (wcc_distance(sorted(:, i - 1), sorted(:, i)) >= gap_fraction*MIN(left_size, right_size)) RETURN
         DO j = 1, n
            delta = ABS(sorted(j, i) - left_gap)
            IF (MIN(delta, 1.0_dp - delta) <= gap_fraction*left_size) RETURN
            delta = ABS(sorted(j, i - 1) - right_gap)
            IF (MIN(delta, 1.0_dp - delta) <= gap_fraction*right_size) RETURN
         END DO
      END DO
      resolved = .TRUE.
   END FUNCTION surface_resolved

! **************************************************************************************************
!> \brief Multiply by the unitary polar factor of an overlap. Reject rank-deficient links.
!> \param product running Wilson matrix
!> \param overlap overlap between adjacent occupied subspaces
!> \param minimum_sv smallest singular value of this link
!> \param status zero on success; negative for invalid/singular input; positive LAPACK failure
!> \param sv_tol rank tolerance
! **************************************************************************************************
   SUBROUTINE wilson_step(product, overlap, minimum_sv, status, sv_tol)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: product(:, :)
      COMPLEX(KIND=dp), INTENT(IN)                       :: overlap(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: minimum_sv
      INTEGER, INTENT(OUT)                               :: status
      REAL(KIND=dp), INTENT(IN)                          :: sv_tol

      COMPLEX(KIND=dp), ALLOCATABLE                      :: a(:, :), u(:, :), vh(:, :), work(:)
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:), sv(:)

      n = SIZE(overlap, 1)
      status = -1
      minimum_sv = 0.0_dp
      IF (n < 1 .OR. SIZE(overlap, 2) /= n .OR. ANY(SHAPE(product) /= [n, n])) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(overlap, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(overlap)))) RETURN
      ALLOCATE (a(n, n), u(n, n), vh(n, n), sv(n), rwork(5*n), work(MAX(1, 4*n)))
      a(:, :) = overlap
      CALL zgesvd('A', 'A', n, n, a, n, sv, u, n, vh, n, work, SIZE(work), rwork, status)
      IF (status /= 0) RETURN
      minimum_sv = MINVAL(sv)
      IF (minimum_sv <= sv_tol) THEN
         status = -2
         RETURN
      END IF
      product = MATMUL(product, MATMUL(u, vh))
   END SUBROUTINE wilson_step

! **************************************************************************************************
!> \brief Sorted Wilson eigenphases in reduced units and total Berry phase in radians.
!> \param product Wilson matrix
!> \param wcc sorted centres in [0,1)
!> \param berry Berry phase in [-pi,pi)
!> \param status zero on success
! **************************************************************************************************
   SUBROUTINE wilson_spectrum(product, wcc, berry, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: product(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: wcc(:), berry
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(KIND=dp)                                   :: dummy(1, 1)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: a(:, :), eig(:), work(:)
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:)

      n = SIZE(product, 1)
      status = -1
      berry = 0.0_dp
      IF (n < 1 .OR. SIZE(product, 2) /= n .OR. SIZE(wcc) /= n) RETURN
      ALLOCATE (a(n, n), eig(n), work(MAX(1, 4*n)), rwork(2*n))
      a(:, :) = product
      CALL zgeev('N', 'N', n, a, n, eig, dummy, 1, dummy, 1, work, SIZE(work), rwork, status)
      IF (status /= 0) RETURN
      wcc = MODULO(ATAN2(AIMAG(eig), REAL(eig, dp))/two_pi, 1.0_dp)
      CALL sort_wcc(wcc)
      berry = two_pi*(MODULO(SUM(wcc) + 0.5_dp, 1.0_dp) - 0.5_dp)
   END SUBROUTINE wilson_spectrum

! **************************************************************************************************
!> \brief Largest-gap crossing parity for an ordered time-reversal half-surface.
!>        Endpoint degeneracy is necessary but does not establish time-reversal symmetry.
!>        Sampling convergence and an isolated fixed-rank subspace must be checked by the caller.
!> \param wcc Wilson centres, band index first and surface-line index second
!> \param invariant Z2 parity, or -1 if checks fail
!> \param status 0 success, -1 invalid input, -2 missing Kramers pairs, -3 ambiguous crossing
!> \param pair_tol tolerance on boundary Kramers degeneracy, in reduced units
! **************************************************************************************************
   SUBROUTINE z2_from_wcc(wcc, invariant, status, pair_tol)
      REAL(KIND=dp), INTENT(IN)                          :: wcc(:, :)
      INTEGER, INTENT(OUT)                               :: invariant, status
      REAL(KIND=dp), INTENT(IN)                          :: pair_tol

      INTEGER                                            :: crossings, i, j, lines, n
      REAL(KIND=dp)                                      :: gap_size, hi, lo
      REAL(KIND=dp), ALLOCATABLE                         :: gaps(:), sorted(:, :)

      invariant = -1
      status = -1
      n = SIZE(wcc, 1)
      lines = SIZE(wcc, 2)
      IF (n < 2 .OR. MOD(n, 2) /= 0 .OR. lines < 2) RETURN
      IF (.NOT. ALL(ieee_is_finite(wcc))) RETURN
      ALLOCATE (sorted(n, lines), gaps(lines))
      sorted(:, :) = MODULO(wcc, 1.0_dp)
      DO i = 1, lines
         CALL sort_wcc(sorted(:, i))
         CALL largest_gap(sorted(:, i), gaps(i), gap_size)
      END DO
      status = -2
      IF (.NOT. kramers_pairs(sorted(:, 1), pair_tol)) RETURN
      IF (.NOT. kramers_pairs(sorted(:, lines), pair_tol)) RETURN
      crossings = 0
      DO i = 2, lines
         lo = MIN(gaps(i - 1), gaps(i))
         hi = MAX(gaps(i - 1), gaps(i))
         DO j = 1, n
            IF (ABS(sorted(j, i) - gaps(i - 1)) < crossing_tol) THEN
               status = -3
               RETURN
            END IF
            IF (sorted(j, i) > lo .AND. sorted(j, i) < hi) crossings = crossings + 1
         END DO
      END DO
      invariant = MOD(crossings, 2)
      status = 0
   END SUBROUTINE z2_from_wcc

! **************************************************************************************************
!> \brief Minimum maximal cyclic matching distance between two WCC sets.
!> \param a first WCC set
!> \param b second WCC set
!> \return distance, or huge if sizes differ
! **************************************************************************************************
   FUNCTION wcc_distance(a, b) RESULT(distance)
      REAL(KIND=dp), INTENT(IN)                          :: a(:), b(:)
      REAL(KIND=dp)                                      :: distance

      INTEGER                                            :: i, j, n, shift
      REAL(KIND=dp)                                      :: delta, error
      REAL(KIND=dp), ALLOCATABLE                         :: aa(:), bb(:)

      distance = HUGE(1.0_dp)
      n = SIZE(a)
      IF (n /= SIZE(b) .OR. n < 1) RETURN
      ALLOCATE (aa(n), bb(n))
      aa(:) = MODULO(a, 1.0_dp)
      bb(:) = MODULO(b, 1.0_dp)
      CALL sort_wcc(aa)
      CALL sort_wcc(bb)
      DO shift = 0, n - 1
         error = 0.0_dp
         DO i = 1, n
            j = MOD(i - 1 + shift, n) + 1
            delta = ABS(aa(i) - bb(j))
            error = MAX(error, MIN(delta, 1.0_dp - delta))
         END DO
         distance = MIN(distance, error)
      END DO
   END FUNCTION wcc_distance

! **************************************************************************************************
!> \brief Test circularly adjacent Kramers pairs (including pairs across the branch cut).
!> \param wcc sorted centres
!> \param tolerance degeneracy tolerance
!> \return whether one of the two cyclic pairings succeeds
! **************************************************************************************************
   FUNCTION kramers_pairs(wcc, tolerance) RESULT(paired)
      REAL(KIND=dp), INTENT(IN)                          :: wcc(:), tolerance
      LOGICAL                                            :: paired

      INTEGER                                            :: i, j, k, n, offset
      LOGICAL                                            :: candidate
      REAL(KIND=dp)                                      :: delta

      paired = .FALSE.
      n = SIZE(wcc)
      DO offset = 0, 1
         candidate = .TRUE.
         DO i = 1, n, 2
            j = MOD(i - 1 + offset, n) + 1
            k = MOD(i + offset, n) + 1
            delta = ABS(wcc(j) - wcc(k))
            candidate = candidate .AND. MIN(delta, 1.0_dp - delta) <= tolerance
         END DO
         paired = paired .OR. candidate
      END DO
   END FUNCTION kramers_pairs

! **************************************************************************************************
!> \brief Find the midpoint of the largest gap between sorted centres.
!> \param wcc sorted centres
!> \param centre gap midpoint
!> \param width gap width
! **************************************************************************************************
   SUBROUTINE largest_gap(wcc, centre, width)
      REAL(KIND=dp), INTENT(IN)                          :: wcc(:)
      REAL(KIND=dp), INTENT(OUT)                         :: centre, width

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: delta

      n = SIZE(wcc)
      width = -1.0_dp
      DO i = 1, n
         IF (i < n) THEN
            delta = wcc(i + 1) - wcc(i)
         ELSE
            delta = wcc(1) + 1.0_dp - wcc(n)
         END IF
         IF (delta > width) THEN
            width = delta
            centre = MODULO(wcc(i) + delta/2.0_dp, 1.0_dp)
         END IF
      END DO
   END SUBROUTINE largest_gap

! **************************************************************************************************
!> \brief In-place insertion sort of Wilson centres.
!> \param values values to sort
! **************************************************************************************************
   SUBROUTINE sort_wcc(values)
      REAL(KIND=dp), INTENT(INOUT)                       :: values(:)

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: value

      DO i = 2, SIZE(values)
         value = values(i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (values(j) <= value) EXIT
            values(j + 1) = values(j)
            j = j - 1
         END DO
         values(j + 1) = value
      END DO
   END SUBROUTINE sort_wcc
END MODULE topology_wilson

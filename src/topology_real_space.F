!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite and periodic physical-window Chern markers in nonorthogonal AO representations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_real_space
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   USE kubo_projected, ONLY: projected_eigensystem
   USE mathconstants, ONLY: twopi

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: finite_chern_marker, periodic_chern_marker
CONTAINS

! **************************************************************************************************
!> \brief Evaluate -2*pi*i/A Tr(Wmixed P [[Xmixed,P],[Ymixed,P]]) without forming S^-1.
!> \param occupied Metric-orthonormal coefficients of an isolated occupied subspace, without spin weights
!> \param metric Positive definite covariant overlap, including spin blocks when applicable
!> \param positions Two covariant position matrices in the ordered oriented plane
!> \param window Covariant matrix of a physical window 0 <= w(r) <= 1, not an AO-index selection
!> \param area Integral of the physical window in the oriented plane, in matching squared-length units
!> \param marker Real finite-window value, not rounded to an integer
!> \param status 0 success; -1 invalid input; -2 bad metric; -3 occupied normalization; -4 window range;
!>               positive LAPACK failure
!> \note Only open finite systems. No wrapping of AO indices or reuse as a periodic Wilson/Z2 invariant.
!>       For Berry connection i<u|grad u>, the plane ordering fixes the sign. Reversing it reverses C.
!>       Projected AO positions need not commute, so a full-window trace need not vanish.
! **************************************************************************************************
   SUBROUTINE finite_chern_marker(occupied, metric, positions, window, area, marker, status)
      COMPLEX(dp), INTENT(IN)                            :: occupied(:, :), metric(:, :), &
                                                            positions(:, :, :), window(:, :)
      REAL(dp), INTENT(IN)                               :: area
      REAL(dp), INTENT(OUT)                              :: marker
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: trace
      COMPLEX(dp), ALLOCATABLE                           :: chol(:, :), cx(:, :), cy(:, :), &
                                                            operators(:, :), PRODUCT(:, :), &
                                                            projector(:, :)
      INTEGER                                            :: axis, i, info, n, nocc
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      marker = 0.0_dp
      status = -1
      n = SIZE(metric, 1)
      nocc = SIZE(occupied, 2)
      IF (n < 1 .OR. SIZE(metric, 2) /= n .OR. SIZE(occupied, 1) /= n .OR. nocc > n) RETURN
      IF (ANY(SHAPE(positions) /= [n, n, 2]) .OR. ANY(SHAPE(window) /= [n, n])) RETURN
      IF (.NOT. ieee_is_finite(area)) RETURN
      IF (area <= 0.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(occupied, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(occupied)))) RETURN
      DO axis = 1, 2
         IF (.NOT. ALL(ieee_is_finite(REAL(positions(:, :, axis), dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(positions(:, :, axis))))) RETURN
         IF (MAXVAL(ABS(positions(:, :, axis) - CONJG(TRANSPOSE(positions(:, :, axis))))) > &
             1.e-10_dp*MAX(TINY(1.0_dp), MAXVAL(ABS(positions(:, :, axis))))) RETURN
      END DO

      CALL validate_window_projector(occupied, metric, window, status)
      IF (status /= 0) RETURN
      IF (nocc == 0 .OR. nocc == n) THEN
         status = 0
         RETURN
      END IF

      ALLOCATE (chol(n, n), operators(n, 3*n), projector(n, n), cx(n, n), cy(n, n), PRODUCT(n, n))
      chol(:, :) = metric
      operators(:, :n) = positions(:, :, 1)
      operators(:, n + 1:2*n) = positions(:, :, 2)
      operators(:, 2*n + 1:) = window
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zposv('U', n, 3*n, chol, n, operators, n, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 1
      IF (info /= 0) RETURN
      projector(:, :) = MATMUL(occupied, MATMUL(CONJG(TRANSPOSE(occupied)), metric))
      cx(:, :) = MATMUL(operators(:, :n), projector) - MATMUL(projector, operators(:, :n))
      cy(:, :) = MATMUL(operators(:, n + 1:2*n), projector) - MATMUL(projector, operators(:, n + 1:2*n))
      PRODUCT(:, :) = MATMUL(operators(:, 2*n + 1:), &
                             MATMUL(projector, MATMUL(cx, cy) - MATMUL(cy, cx)))
      trace = 0.0_dp
      DO i = 1, n
         trace = trace + PRODUCT(i, i)
      END DO
      status = -1
      IF (.NOT. ieee_is_finite(REAL(trace, dp)) .OR. .NOT. ieee_is_finite(AIMAG(trace))) RETURN
      IF (ABS(REAL(trace, dp)) > 1.e-9_dp*MAX(1.0_dp, ABS(trace))) RETURN
      marker = twopi*(AIMAG(trace)/area)
      IF (.NOT. ieee_is_finite(marker)) THEN
         marker = 0.0_dp
         RETURN
      END IF
      status = 0
   END SUBROUTINE finite_chern_marker

! **************************************************************************************************
!> \brief Symmetric periodic Chern marker from metric-dual states at +/- reciprocal torus vectors.
!> \param occupied Metric-orthonormal occupied coefficients on the complete Born-von-Karman torus
!> \param metric Positive definite covariant AO overlap
!> \param twists Covariant integrals of exp(-i*b_j.r), with reciprocal torus vectors in plane order
!> \param window Covariant physical periodic window, 0 <= W <= S
!> \param cell_area Physical area of the complete torus, not of the primitive cell
!> \param window_area Integral of the physical window on that torus
!> \param eps_dual Minimum singular value of occupied twist overlaps, between zero and one
!> \param marker Unrounded window average in the same Berry/orientation convention as the finite marker
!> \param minimum_sv Minimum occupied twist-overlap singular value
!> \param status Finite-marker statuses, or -6 for an unresolved dual-state overlap
!> \note Implements Bau/Marrazzo (2024), symmetric Eq. (15). The dual-state outer products are
!>       not orthogonal projectors. Do not orthonormalize them or replace them by polar factors.
! **************************************************************************************************
   SUBROUTINE periodic_chern_marker(occupied, metric, twists, window, cell_area, window_area, eps_dual, &
                                    marker, minimum_sv, status)
      COMPLEX(dp), INTENT(IN)                            :: occupied(:, :), metric(:, :), &
                                                            twists(:, :, :), window(:, :)
      REAL(dp), INTENT(IN)                               :: cell_area, window_area, eps_dual
      REAL(dp), INTENT(OUT)                              :: marker, minimum_sv
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: dummy(1, 1), trace
      COMPLEX(dp), ALLOCATABLE :: applied(:, :), chol(:, :), dual(:, :), factor(:, :), &
                                  overlap(:, :), PRODUCT(:, :), projector(:, :), q(:, :, :), rhs(:, :), weighted(:, :), &
                                  work(:)
      INTEGER                                            :: axis, branch, first, i, index, info, m, n
      INTEGER, ALLOCATABLE                               :: pivots(:)
      REAL(dp)                                           :: sign
      REAL(dp), ALLOCATABLE                              :: rwork(:), singular(:)
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      marker = 0.0_dp
      minimum_sv = 1.0_dp
      status = -1
      n = SIZE(metric, 1)
      m = SIZE(occupied, 2)
      IF (n < 1 .OR. ANY(SHAPE(twists) /= [n, n, 2])) RETURN
      IF (.NOT. ALL(ieee_is_finite([cell_area, window_area, eps_dual]))) RETURN
      IF (cell_area <= 0.0_dp .OR. window_area <= 0.0_dp) RETURN
      IF (window_area > cell_area*(1.0_dp + 1.e-10_dp)) RETURN
      IF (eps_dual <= 0.0_dp .OR. eps_dual >= 1.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(twists, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(twists)))) RETURN
      CALL validate_window_projector(occupied, metric, window, status)
      IF (status /= 0 .OR. m == 0 .OR. m == n) RETURN
      ALLOCATE (chol(n, n), applied(n, 4*m), factor(m, m), overlap(m, m), rhs(m, n), &
                dual(n, m), q(n, n, 2), projector(n, n), weighted(n, n), PRODUCT(n, n), &
                pivots(m), singular(m), rwork(5*m), work(5*m))
      DO axis = 1, 2
         first = 2*(axis - 1)*m
         applied(:, first + 1:first + m) = MATMUL(twists(:, :, axis), occupied)
         applied(:, first + m + 1:first + 2*m) = MATMUL(CONJG(TRANSPOSE(twists(:, :, axis))), occupied)
      END DO
      chol(:, :) = metric
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zposv('U', n, 4*m, chol, n, applied, n, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 1
      IF (info /= 0) RETURN
      q(:, :, :) = 0.0_dp
      DO axis = 1, 2
         DO branch = 1, 2
            index = 2*(axis - 1) + branch - 1
            first = index*m
            dual(:, :) = applied(:, first + 1:first + m)
            overlap(:, :) = MATMUL(CONJG(TRANSPOSE(occupied)), MATMUL(metric, dual))
            factor(:, :) = overlap
            ! LAPACK's values-only SVD probes IEEE arithmetic internally.
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
            CALL zgesvd('N', 'N', m, m, factor, m, singular, dummy, 1, dummy, 1, &
                        work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
            status = ABS(info) + 1
            IF (info /= 0) RETURN
            minimum_sv = MIN(minimum_sv, MINVAL(singular))
            status = -6
            IF (.NOT. ALL(ieee_is_finite(singular)) .OR. minimum_sv <= eps_dual) RETURN
            ! Right solve D*M= S^-1*E*C. No explicit inverse, and no change of the dual gauge.
            factor(:, :) = TRANSPOSE(overlap)
            rhs(:, :) = TRANSPOSE(dual)
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
            CALL zgesv(m, n, factor, m, pivots, rhs, m, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
            status = ABS(info) + 1
            IF (info /= 0) RETURN
            dual(:, :) = TRANSPOSE(rhs)
            sign = REAL(3 - 2*branch, dp)
            q(:, :, axis) = q(:, :, axis) + sign*MATMUL(dual, MATMUL(CONJG(TRANSPOSE(dual)), metric))
         END DO
      END DO
      projector(:, :) = MATMUL(occupied, MATMUL(CONJG(TRANSPOSE(occupied)), metric))
      chol(:, :) = metric
      weighted(:, :) = window
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zposv('U', n, n, chol, n, weighted, n, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 1
      IF (info /= 0) RETURN
      PRODUCT(:, :) = MATMUL(weighted, MATMUL(MATMUL(q(:, :, 1), q(:, :, 2)) - &
                                              MATMUL(q(:, :, 2), q(:, :, 1)), projector))
      trace = 0.0_dp
      DO i = 1, n
         trace = trace + PRODUCT(i, i)
      END DO
      status = -1
      IF (.NOT. ieee_is_finite(REAL(trace, dp)) .OR. .NOT. ieee_is_finite(AIMAG(trace))) RETURN
      marker = -(cell_area/window_area)*AIMAG(trace)/(4*twopi)
      IF (.NOT. ieee_is_finite(marker)) THEN
         marker = 0.0_dp
         RETURN
      END IF
      status = 0
   END SUBROUTINE periodic_chern_marker

! **************************************************************************************************
!> \brief Common shape, positive-metric, physical-window and occupied normalization checks.
!> \param occupied Metric-orthonormal occupied coefficients
!> \param metric Covariant positive overlap
!> \param window Covariant physical window
!> \param status Zero for success or shared marker error status
! **************************************************************************************************
   SUBROUTINE validate_window_projector(occupied, metric, window, status)
      COMPLEX(dp), INTENT(IN)                            :: occupied(:, :), metric(:, :), &
                                                            window(:, :)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: frame(:, :), gram(:, :)
      INTEGER                                            :: i, m, n
      REAL(dp), ALLOCATABLE                              :: spectrum(:)

      status = -1
      n = SIZE(metric, 1)
      m = SIZE(occupied, 2)
      IF (n < 1 .OR. SIZE(metric, 2) /= n .OR. SIZE(occupied, 1) /= n .OR. m > n) RETURN
      IF (ANY(SHAPE(window) /= [n, n])) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(occupied, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(occupied)))) RETURN
      ALLOCATE (frame(n, n), spectrum(n), gram(m, m))
      CALL projected_eigensystem(window, metric, spectrum, frame, status)
      IF (status /= 0) RETURN
      status = -4
      IF (MINVAL(spectrum) < -1.e-10_dp .OR. MAXVAL(spectrum) > 1.0_dp + 1.e-10_dp) RETURN
      gram(:, :) = MATMUL(CONJG(TRANSPOSE(occupied)), MATMUL(metric, occupied))
      DO i = 1, m
         gram(i, i) = gram(i, i) - 1.0_dp
      END DO
      status = -3
      IF (m > 0) THEN
         IF (MAXVAL(ABS(gram)) > 1.e-9_dp) RETURN
      END IF
      status = 0
   END SUBROUTINE validate_window_projector
END MODULE topology_real_space

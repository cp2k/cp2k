!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite-window Chern markers in nonorthogonal AO representations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE topology_real_space
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE kubo_projected,                  ONLY: projected_eigensystem
   USE mathconstants,                   ONLY: twopi

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: finite_chern_marker
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
                                                            frame(:, :), gram(:, :), &
                                                            operators(:, :), PRODUCT(:, :), &
                                                            projector(:, :)
      INTEGER                                            :: axis, i, info, n, nocc
      REAL(dp), ALLOCATABLE                              :: spectrum(:)

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

      ALLOCATE (frame(n, n), spectrum(n))
      ! The shared generalized solver checks both Hermiticity and metric conditioning.
      CALL projected_eigensystem(window, metric, spectrum, frame, status)
      IF (status /= 0) RETURN
      status = -4
      IF (MINVAL(spectrum) < -1.e-10_dp .OR. MAXVAL(spectrum) > 1.0_dp + 1.e-10_dp) RETURN
      DEALLOCATE (frame, spectrum)

      ALLOCATE (gram(nocc, nocc))
      gram(:, :) = MATMUL(CONJG(TRANSPOSE(occupied)), MATMUL(metric, occupied))
      DO i = 1, nocc
         gram(i, i) = gram(i, i) - 1.0_dp
      END DO
      status = -3
      IF (nocc > 0) THEN
         IF (MAXVAL(ABS(gram)) > 1.e-9_dp) RETURN
      END IF
      IF (nocc == 0 .OR. nocc == n) THEN
         status = 0
         RETURN
      END IF

      ALLOCATE (chol(n, n), operators(n, 3*n), projector(n, n), cx(n, n), cy(n, n), PRODUCT(n, n))
      chol(:, :) = metric
      operators(:, :n) = positions(:, :, 1)
      operators(:, n + 1:2*n) = positions(:, :, 2)
      operators(:, 2*n + 1:) = window
      CALL zposv('U', n, 3*n, chol, n, operators, n, info)
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
END MODULE topology_real_space

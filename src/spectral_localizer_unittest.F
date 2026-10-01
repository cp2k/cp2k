!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Independent finite-model, metric, orientation and Pfaffian tests.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM spectral_localizer_unittest
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode
#endif
   USE kinds, ONLY: dp
   USE spectral_localizer, ONLY: localizer_aii_skew, &
                                 localizer_evaluate, &
                                 localizer_inertia, &
                                 localizer_periodic_evaluate, &
                                 localizer_pfaffian_sign

   IMPLICIT NONE
#if defined(__HAS_IEEE_EXCEPTIONS)
   LOGICAL :: original_halt(SIZE(IEEE_ALL)), current_halt(SIZE(IEEE_ALL))
#endif
   COMPLEX(KIND=dp), ALLOCATABLE :: h(:, :), s(:, :), r(:, :, :), transform(:, :), &
                                    hc(:, :), sc(:, :), rc(:, :, :), hs(:, :), ss(:, :), rs(:, :, :)
   INTEGER :: axis, i, index, reference, status, n, pf, permutation(4)
   REAL(KIND=dp) :: gap, gap_reference, ratio, skew(4, 4), exact

#if defined(__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, original_halt)
#endif
   CALL inertia_tests()
   CALL real_inertia_tests()
   CALL gap_tests()
   CALL pfaffian_orientation_tests()
   CALL aii_transform_tests()
   CALL aii_roundoff_tests()
   CALL periodic_tests()
   ! Independent Pfaffian polynomial for size four, including odd permutations.
   skew = 0.0_dp
   skew(1, 2) = 0.3_dp
   skew(1, 3) = 2.0_dp
   skew(1, 4) = -1.0_dp
   skew(2, 3) = 0.7_dp
   skew(2, 4) = 3.0_dp
   skew(3, 4) = 0.8_dp
   skew = skew - TRANSPOSE(skew)
   exact = skew(1, 2)*skew(3, 4) - skew(1, 3)*skew(2, 4) + skew(1, 4)*skew(2, 3)
   CALL localizer_pfaffian_sign(skew, pf, status)
   IF (status /= 0 .OR. pf /= INT(SIGN(1.0_dp, exact))) ERROR STOP 'Pfaffian polynomial failed'
   permutation = [2, 1, 3, 4]
   CALL localizer_pfaffian_sign(skew(permutation, permutation), pf, status)
   IF (status /= 0 .OR. pf /= -INT(SIGN(1.0_dp, exact))) ERROR STOP 'Pfaffian orientation failed'
   skew = 0.0_dp
   CALL localizer_pfaffian_sign(skew, pf, status)
   IF (status == 0 .OR. pf /= 0) ERROR STOP 'Singular Pfaffian accepted'

   ALLOCATE (h(1, 1), s(1, 1), r(1, 1, 2))
   h = 4.0_dp
   s = 2.0_dp
   r(1, 1, :) = [6.0_dp, 8.0_dp]
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp], 1.0_dp, 0.2_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - SQRT(2.0_dp)) > 1.e-12_dp) THEN
      ERROR STOP 'Analytic generalized localizer failed'
   END IF
   CALL evaluate(h, s, r, [3.0_dp, 4.0_dp], 2.0_dp, 0.2_dp, .FALSE., index, gap, status)
   IF (status /= 1) ERROR STOP 'Zero gap was assigned an index'
   DEALLOCATE (h, s, r)

   ALLOCATE (h(2, 2), s(2, 2), r(2, 2, 2))
   h = 0.0_dp
   s = 0.0_dp
   r = 0.0_dp
   s(1, 1) = 1.e-9_dp
   s(2, 2) = 1.0_dp
   h(1, 1) = 1.e-16_dp
   h(2, 2) = 1.0_dp
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.2_dp, .FALSE., index, gap, status)
   IF (status /= 1) ERROR STOP 'Metric-amplified roundoff ignored'
   DEALLOCATE (h, s, r)

   CALL qwz_model(7, -1.0_dp, h, s, r)
   n = SIZE(h, 1)
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., reference, gap_reference, status)
   IF (status /= 0 .OR. reference /= -1 .OR. ABS(gap_reference - 0.740327333453389_dp) > 1.e-11_dp) THEN
      ERROR STOP 'Nontrivial open QWZ model failed'
   END IF
   CALL evaluate(h, s, r(:, :, [2, 1]), [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= -reference .OR. ABS(gap - gap_reference) > 1.e-11_dp) THEN
      ERROR STOP 'Oriented-plane reversal failed'
   END IF
   CALL evaluate(CONJG(h), s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= -reference) ERROR STOP 'Time reversed Chern sign failed'
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.15_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= reference) ERROR STOP 'Kappa plateau failed'
   CALL evaluate(h, s, r, [20.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= 0) ERROR STOP 'Outside-vacuum limit failed'

   ALLOCATE (transform(n, n), hc(n, n), sc(n, n), rc(n, n, 2))
   transform = 0.0_dp
   DO i = 1, n
      transform(i, i) = (1.0_dp + REAL(MOD(i, 5), dp)/5.0_dp)*EXP(CMPLX(0.0_dp, 0.17_dp*i, dp))
      IF (i < n) transform(i, i + 1) = CMPLX(0.05_dp, 0.03_dp, dp)
   END DO
   hc(:, :) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(h, transform))
   sc(:, :) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(s, transform))
   DO axis = 1, 2
      rc(:, :, axis) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(r(:, :, axis), transform))
   END DO
   CALL evaluate(hc, sc, rc, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= reference .OR. ABS(gap - gap_reference) > 1.e-11_dp) THEN
      ERROR STOP 'Complex nonorthogonal basis invariance failed'
   END IF
   rc(:, :, 1) = rc(:, :, 1) + 13.0_dp*sc
   rc(:, :, 2) = rc(:, :, 2) - 7.0_dp*sc
   hc(:, :) = hc + 11.0_dp*sc
   CALL evaluate(hc, sc, rc, [13.0_dp, -7.0_dp], 11.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= reference .OR. ABS(gap - gap_reference) > 1.e-11_dp) THEN
      ERROR STOP 'Origin/energy shift covariance failed'
   END IF
   sc = 0.0_dp
   DO i = 1, n
      sc(i, i) = 1.0_dp
   END DO
   sc(1, 1) = 1.e-14_dp
   CALL evaluate(hc, sc, rc, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= -2) ERROR STOP 'Ill-conditioned metric accepted'
   sc(1, 1) = -1.0_dp
   CALL evaluate(hc, sc, rc, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= -2) ERROR STOP 'Indefinite metric accepted'
   hc(1, 2) = hc(1, 2) + 1.0_dp
   CALL evaluate(hc, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= -1) ERROR STOP 'Non-Hermitian Hamiltonian accepted'

   ! Two time-reversed Chern copies, with TR-preserving, spin-nonconserving mixing.
   ALLOCATE (hs(2*n, 2*n), ss(2*n, 2*n), rs(2*n, 2*n, 2))
   hs = 0.0_dp
   ss = 0.0_dp
   rs = 0.0_dp
   hs(:n, :n) = h
   hs(n + 1:, n + 1:) = CONJG(h)
   ss(:n, :n) = s
   ss(n + 1:, n + 1:) = s
   rs(:n, :n, :) = r
   rs(n + 1:, n + 1:, :) = r
   DO i = 1, n, 2
      hs(i, n + i + 1) = 0.1_dp
      hs(i + 1, n + i) = -0.1_dp
   END DO
   hs(n + 1:, :n) = CONJG(TRANSPOSE(hs(:n, n + 1:)))
   CALL evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., index, gap, status)
   IF (status /= 0 .OR. index /= 1) ERROR STOP 'Spin-mixed nontrivial AII model failed'
   gap_reference = gap
   ! Real, nonorthogonal AO change applied to both spin components preserves fixed Theta.
   transform = REAL(transform, dp)
   DO i = 1, n
      transform(i, i) = 1.0_dp + REAL(MOD(i, 5), dp)/5.0_dp
   END DO
   DO axis = 1, 2
      DO i = 1, 2
         hs((i - 1)*n + 1:i*n, (axis - 1)*n + 1:axis*n) = &
            MATMUL(TRANSPOSE(transform), MATMUL(hs((i - 1)*n + 1:i*n, (axis - 1)*n + 1:axis*n), transform))
      END DO
   END DO
   sc(:, :) = MATMUL(TRANSPOSE(transform), MATMUL(s, transform))
   ss(:n, :n) = sc
   ss(n + 1:, n + 1:) = sc
   DO axis = 1, 2
      rc(:, :, axis) = MATMUL(TRANSPOSE(transform), MATMUL(r(:, :, axis), transform))
      rs(:n, :n, axis) = rc(:, :, axis)
      rs(n + 1:, n + 1:, axis) = rc(:, :, axis)
   END DO
   CALL evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., index, gap, status)
   IF (status /= 0 .OR. index /= 1 .OR. ABS(gap - gap_reference) > 1.e-11_dp) THEN
      ERROR STOP 'AII metric covariance failed'
   END IF
   hs(1, 1) = hs(1, 1) + 1.e-12_dp
   CALL evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., index, gap, status)
   IF (status /= 0 .OR. index /= 1) ERROR STOP 'Roundoff-level AII residual failed'
   hs(1, 1) = hs(1, 1) + 1.0_dp
   CALL evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., index, gap, status)
   IF (status /= -3) ERROR STOP 'Broken physical time reversal accepted'
   hs(:, :) = ss
   CALL evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., index, gap, status)
   IF (status /= 0 .OR. index /= 0) ERROR STOP 'Atomic AII calibration failed'

   DEALLOCATE (h, s, r)
   CALL qwz_model(7, -3.0_dp, h, s, r)
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., index, gap, status)
   IF (status /= 0 .OR. index /= 0) ERROR STOP 'Trivial QWZ model failed'
#if defined(__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, current_halt)
   IF (ANY(current_halt .NEQV. original_halt)) ERROR STOP 'Caller floating-point traps changed'
#endif

CONTAINS

! **************************************************************************************************
!> \brief Signed-permutation AII conversion against the explicit dense unitary transform.
! **************************************************************************************************
   SUBROUTINE aii_transform_tests()
      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :), q(:, :), reference(:, :)
      INTEGER                                            :: i, j, m, n, status
      REAL(dp)                                           :: error, largest, residual
      REAL(dp), ALLOCATABLE                              :: actual(:, :), expected(:, :)

      largest = 0.0_dp
      DO n = 1, 13
         m = 4*n
         ALLOCATE (q(m, m), matrix(m, m), reference(m, m), expected(m, m), actual(m, m))
         q = 0.0_dp
         DO i = 1, m
            q(i, i) = 1.0_dp
            DO j = 1, m
               expected(i, j) = SIN(REAL(i + 2*j, dp)) + COS(REAL(3*i - j, dp))
            END DO
         END DO
         expected = expected - TRANSPOSE(expected)
         DO i = 1, n
            q(i, 3*n + i) = CMPLX(0.0_dp, -1.0_dp, dp)
            q(n + i, 2*n + i) = CMPLX(0.0_dp, 1.0_dp, dp)
            q(2*n + i, n + i) = CMPLX(0.0_dp, 1.0_dp, dp)
            q(3*n + i, i) = CMPLX(0.0_dp, -1.0_dp, dp)
         END DO
         q = q/SQRT(2.0_dp)
         matrix = CMPLX(0.0_dp, -1.0_dp, dp)*MATMUL(q, MATMUL(expected, CONJG(TRANSPOSE(q))))
         reference = MATMUL(CONJG(TRANSPOSE(q)), MATMUL(matrix, q))
         CALL localizer_aii_skew(matrix, actual, 1.e-10_dp, residual, status)
         error = MAX(MAXVAL(ABS(actual - expected)), MAXVAL(ABS(actual + AIMAG(reference))))
         IF (status /= 0 .OR. error > 1.e-13_dp) ERROR STOP 'AII sparse-unitary transform failed'
         largest = MAX(largest, error)
         DO i = 1, m
            matrix(i, i) = matrix(i, i) - CMPLX(0.0_dp, 5.e-12_dp, dp)
         END DO
         CALL localizer_aii_skew(matrix, actual, 1.e-10_dp, residual, status)
         IF (status /= 0 .OR. MAXVAL(ABS(actual - expected)) > 1.e-13_dp .OR. residual < 4.9e-12_dp) THEN
            ERROR STOP 'AII skew projection or discarded-error bound failed'
         END IF
         IF (MAXVAL(ABS(actual + TRANSPOSE(actual))) /= 0.0_dp) ERROR STOP 'AII projection is not exactly skew'
         DO i = 1, m
            matrix(i, i) = matrix(i, i) + CMPLX(0.0_dp, 5.e-12_dp, dp)
         END DO
         matrix(1, 1) = matrix(1, 1) + 0.1_dp
         CALL localizer_aii_skew(matrix, actual, 1.e-10_dp, residual, status)
         IF (status == 0) ERROR STOP 'Broken AII transform accepted'
         CALL localizer_aii_skew(matrix(:m - 1, :m - 1), actual, 1.e-10_dp, residual, status)
         IF (status == 0) ERROR STOP 'Invalid AII dimensions accepted'
         DEALLOCATE (q, matrix, reference, expected, actual)
      END DO
      WRITE (*, '(A,ES14.5)') 'AII signed-permutation maximum error: ', largest
   END SUBROUTINE aii_transform_tests

! **************************************************************************************************
!> \brief Accepted Hermiticity roundoff must not fail the stricter Pfaffian skew check.
! **************************************************************************************************
   SUBROUTINE aii_roundoff_tests()
      COMPLEX(dp)                                        :: h0(4, 4), r0(4, 4, 2), s0(4, 4)
      INTEGER                                            :: i, index, status
      REAL(dp)                                           :: gap

      h0 = 0.0_dp
      s0 = 0.0_dp
      r0 = 0.0_dp
      DO i = 1, 4
         s0(i, i) = 1.0_dp
         h0(i, i) = REAL(3*MOD(i + 1, 2) - 1, dp)
      END DO
      h0(1, 2) = 5.e-12_dp
      h0(3, 4) = 5.e-12_dp
      CALL evaluate(h0, s0, r0, [0.0_dp, 0.0_dp], 0.0_dp, 0.2_dp, .TRUE., index, gap, status)
      IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - 1.0_dp) > 1.e-12_dp) ERROR STOP 'AII roundoff rejected'
      h0(1, 1) = -1.e-10_dp
      h0(3, 3) = -1.e-10_dp
      CALL evaluate(h0, s0, r0, [0.0_dp, 0.0_dp], 0.0_dp, 0.2_dp, .TRUE., index, gap, status)
      IF (status /= 1) ERROR STOP 'AII unresolved gap assigned an index'
      h0(1, 2) = 1.e-5_dp
      CALL evaluate(h0, s0, r0, [0.0_dp, 0.0_dp], 0.0_dp, 0.2_dp, .TRUE., index, gap, status)
      IF (status /= -1) ERROR STOP 'AII non-Hermitian input projected away'
   END SUBROUTINE aii_roundoff_tests

! **************************************************************************************************
!> \brief Torus models with known Chern/Z2 phases, origin periodicity and metric covariance.
! **************************************************************************************************
   SUBROUTINE periodic_tests()
      INTEGER, PARAMETER                                 :: length = 8

      COMPLEX(dp), ALLOCATABLE :: co(:, :, :), ct(:, :, :), h0(:, :), ht(:, :), pos(:, :, :), &
         s0(:, :), si(:, :, :), sn(:, :, :), spin_c(:, :, :), spin_h(:, :), spin_s(:, :), &
         spin_sin(:, :, :), st(:, :), t(:, :)
      INTEGER                                            :: col, expected, i, idx, info, j, kx, ky, &
                                                            n, orbital, phase, row, x, y
      REAL(dp)                                           :: cond, g, g0, mass, two_pi, wave(2)

      two_pi = 8.0_dp*ATAN(1.0_dp)
      DO phase = 1, 3
         mass = -1.0_dp
         expected = -1
         IF (phase == 2) THEN
            mass = 1.0_dp
            expected = 1
         ELSE IF (phase == 3) THEN
            mass = -3.0_dp
            expected = 0
         END IF
         CALL qwz_model(length, mass, h0, s0, pos, torus=.TRUE.)
         n = SIZE(h0, 1)
         ALLOCATE (co(n, n, 2), si(n, n, 2))
         co = 0.0_dp
         si = 0.0_dp
         DO j = 1, 2
            DO i = 1, n
               co(i, i, j) = COS(two_pi*REAL(pos(i, i, j), dp)/length)
               si(i, i, j) = SIN(two_pi*REAL(pos(i, i, j), dp)/length)
            END DO
         END DO
         CALL localizer_periodic_evaluate(h0, s0, co, si, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g0, cond, info)
         IF (info /= 0 .OR. idx /= expected) ERROR STOP 'Periodic QWZ Chern phase failed'
         CALL localizer_periodic_evaluate(h0, s0, co, si, [two_pi, -2*two_pi], 0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 0 .OR. idx /= expected .OR. ABS(g - g0) > 1.e-11_dp) THEN
            ERROR STOP 'Torus origin periodicity failed'
         END IF
         CALL localizer_periodic_evaluate(h0, s0, co, si, two_pi*[1.0_dp, 2.0_dp]/length, &
                                          0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 0 .OR. idx /= expected .OR. ABS(g - g0) > 1.e-11_dp) THEN
            ERROR STOP 'Periodic lattice translation covariance failed'
         END IF
         CALL localizer_periodic_evaluate(h0, s0, co(:, :, [2, 1]), si(:, :, [2, 1]), &
                                          [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 0 .OR. idx /= -expected .OR. ABS(g - g0) > 1.e-11_dp) THEN
            ERROR STOP 'Periodic orientation reversal failed'
         END IF

         ALLOCATE (spin_h(2*n, 2*n), spin_s(2*n, 2*n), spin_c(2*n, 2*n, 2), spin_sin(2*n, 2*n, 2))
         spin_h = 0.0_dp
         spin_s = 0.0_dp
         spin_c = 0.0_dp
         spin_sin = 0.0_dp
         spin_h(:n, :n) = h0
         spin_h(n + 1:, n + 1:) = CONJG(h0)
         DO i = 1, n, 2
            spin_h(i, n + i + 1) = 0.1_dp
            spin_h(i + 1, n + i) = -0.1_dp
         END DO
         spin_h(n + 1:, :n) = CONJG(TRANSPOSE(spin_h(:n, n + 1:)))
         spin_s(:n, :n) = s0
         spin_s(n + 1:, n + 1:) = s0
         spin_c(:n, :n, :) = co
         spin_c(n + 1:, n + 1:, :) = co
         spin_sin(:n, :n, :) = si
         spin_sin(n + 1:, n + 1:, :) = si
         CALL localizer_periodic_evaluate(spin_h, spin_s, spin_c, spin_sin, [0.0_dp, 0.0_dp], &
                                          0.0_dp, 1.0_dp, .TRUE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 0 .OR. idx /= ABS(expected)) ERROR STOP 'Periodic spin-mixed Z2 phase failed'
         spin_h = spin_s
         CALL localizer_periodic_evaluate(spin_h, spin_s, spin_c, spin_sin, [0.0_dp, 0.0_dp], &
                                          0.0_dp, 1.0_dp, .TRUE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 0 .OR. idx /= 0) ERROR STOP 'Periodic atomic Pfaffian calibration failed'
         DEALLOCATE (spin_h, spin_s, spin_c, spin_sin)

         IF (phase == 1) THEN
            ALLOCATE (t(n, n), ht(n, n), st(n, n), ct(n, n, 2), sn(n, n, 2))
            ! Full Fourier transform: H becomes block diagonal, but torus coordinates do not.
            t = 0.0_dp
            DO ky = 0, length - 1
               DO kx = 0, length - 1
                  col = 2*(ky*length + kx) + 1
                  DO y = 0, length - 1
                     DO x = 0, length - 1
                        row = 2*(y*length + x) + 1
                        DO orbital = 0, 1
                           t(row + orbital, col + orbital) = &
                              EXP(CMPLX(0.0_dp, two_pi*(x*kx + y*ky)/length, dp))/length
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            ht = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(h0, t))
            st = s0
            DO ky = 0, length - 1
               DO kx = 0, length - 1
                  col = 2*(ky*length + kx) + 1
                  wave = two_pi*[REAL(kx, dp), REAL(ky, dp)]/length
                  IF (ABS(ht(col, col) - (mass + SUM(COS(wave)))) > 1.e-12_dp) THEN
                     ERROR STOP 'Periodic real-space/k-space Hamiltonian mismatch'
                  END IF
                  IF (ABS(ht(col, col + 1) - CMPLX(SIN(wave(1)), -SIN(wave(2)), dp)) > 1.e-12_dp) THEN
                     ERROR STOP 'Periodic Bloch Hamiltonian offdiagonal mismatch'
                  END IF
               END DO
            END DO
            DO j = 1, 2
               ct(:, :, j) = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(co(:, :, j), t))
               sn(:, :, j) = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(si(:, :, j), t))
            END DO
            CALL localizer_periodic_evaluate(ht, st, ct, sn, [0.0_dp, 0.0_dp], &
                                             0.0_dp, 1.0_dp, .FALSE., &
                                             1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
            IF (info /= 0 .OR. idx /= expected .OR. ABS(g - g0) > 1.e-11_dp) THEN
               ERROR STOP 'Full k-mesh/real-space periodic localizer mismatch'
            END IF
            t = 0.0_dp
            DO i = 1, n
               t(i, i) = (1.0_dp + MOD(i, 5)/5.0_dp)*EXP(CMPLX(0.0_dp, 0.17_dp*i, dp))
               IF (i < n) t(i, i + 1) = CMPLX(0.07_dp, -0.04_dp, dp)
            END DO
            ht = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(h0, t))
            st = MATMUL(CONJG(TRANSPOSE(t)), t)
            DO j = 1, 2
               ct(:, :, j) = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(co(:, :, j), t))
               sn(:, :, j) = MATMUL(CONJG(TRANSPOSE(t)), MATMUL(si(:, :, j), t))
            END DO
            CALL localizer_periodic_evaluate(ht + 7.0_dp*st, st, ct, sn, [0.0_dp, 0.0_dp], &
                                             7.0_dp, 1.0_dp, .FALSE., &
                                             1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
            IF (info /= 0 .OR. idx /= expected .OR. ABS(g - g0) > 1.e-11_dp) THEN
               ERROR STOP 'Periodic metric/energy covariance failed'
            END IF
            CALL localizer_periodic_evaluate(2*ht, st, ct, sn, [0.0_dp, 0.0_dp], &
                                             0.0_dp, 2.0_dp, .FALSE., &
                                             1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
            IF (info /= 0 .OR. idx /= expected .OR. ABS(g - 2*g0) > 1.e-11_dp) THEN
               ERROR STOP 'Periodic energy-scale covariance failed'
            END IF
            DEALLOCATE (t, ht, st, ct, sn)
         END IF
         CALL localizer_periodic_evaluate(h0, s0, co, si, [0.0_dp, 0.0_dp], 0.0_dp, 0.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= -1) ERROR STOP 'Zero torus energy scale accepted'
         co(1, 2, 1) = co(1, 2, 1) + 0.1_dp
         CALL localizer_periodic_evaluate(h0, s0, co, si, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= -1) ERROR STOP 'Non-Hermitian torus cosine accepted'
         h0 = 0.0_dp
         si = 0.0_dp
         co(:, :, 1) = s0
         co(:, :, 2) = s0
         CALL localizer_periodic_evaluate(h0, s0, co, si, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                          1.e-10_dp, 1.e-10_dp, 1.e-9_dp, idx, g, cond, info)
         IF (info /= 1) ERROR STOP 'Closed periodic localizer assigned an index'
         DEALLOCATE (h0, s0, pos, co, si)
      END DO
   END SUBROUTINE periodic_tests

! **************************************************************************************************
!> \brief Exercise 2x2 pivots, complex congruences, scale changes and singularity.
! **************************************************************************************************
   SUBROUTINE inertia_tests()
      COMPLEX(KIND=dp)                                   :: change(4, 4), matrix(4, 4), &
                                                            transformed(4, 4)
      INTEGER                                            :: exponent, i, j, m, p, RESULT

      matrix = 0.0_dp
      matrix(1, 2) = CMPLX(0.0_dp, 2.0_dp, dp)
      matrix(2, 1) = CONJG(matrix(1, 2))
      matrix(3, 3) = 1.e-12_dp
      matrix(4, 4) = 3.0_dp
      change = 0.0_dp
      DO i = 1, 4
         change(i, i) = CMPLX(1.0_dp + 0.2_dp*i, 0.1_dp*i, dp)
         DO j = i + 1, 4
            change(i, j) = CMPLX(0.03_dp*j, -0.07_dp*i, dp)
         END DO
      END DO
      transformed = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(matrix, change))
      DO exponent = -150, 150, 150
         CALL localizer_inertia((10.0_dp**exponent)*transformed, p, m, RESULT)
         IF (RESULT /= 0 .OR. p /= 3 .OR. m /= 1) ERROR STOP 'LDL block inertia or congruence failed'
      END DO
      matrix(3, 3) = -1.e-12_dp
      CALL localizer_inertia(matrix, p, m, RESULT)
      IF (RESULT /= 0 .OR. p /= 2 .OR. m /= 2) ERROR STOP 'LDL gap crossing failed'
      matrix(3, 3) = 0.0_dp
      CALL localizer_inertia(matrix, p, m, RESULT)
      IF (RESULT /= 1) ERROR STOP 'Singular LDL accepted'
      matrix(2, 1) = matrix(1, 2)
      CALL localizer_inertia(matrix, p, m, RESULT)
      IF (RESULT /= -1) ERROR STOP 'Complex symmetric matrix accepted as Hermitian'
      matrix = 0.0_dp
      matrix(1, 2) = CMPLX(0.0_dp, 1.e-12_dp, dp)
      matrix(2, 1) = CONJG(matrix(1, 2))
      matrix(3, 3) = 1.0_dp
      matrix(4, 4) = -1.0_dp
      CALL localizer_inertia(matrix, p, m, RESULT)
      IF (RESULT /= 0 .OR. p /= 2 .OR. m /= 2) ERROR STOP 'Small imaginary coupling discarded'
   END SUBROUTINE inertia_tests

! **************************************************************************************************
!> \brief Cross-check real and complex LDL under congruences, all signatures and permutations.
! **************************************************************************************************
   SUBROUTINE real_inertia_tests()
      COMPLEX(dp), ALLOCATABLE                           :: cmatrix(:, :), phases(:)
      INTEGER                                            :: exponent, i, info, j, m, n, negative, p
      REAL(dp)                                           :: pair(2, 2), scale
      REAL(dp), ALLOCATABLE                              :: change(:, :), diagonal(:, :), &
                                                            matrix(:, :)

      DO n = 1, 12
         ALLOCATE (diagonal(n, n), change(n, n), matrix(n, n), phases(n), cmatrix(n, n))
         change = 0.0_dp
         DO i = 1, n
            change(i, i) = 0.5_dp + REAL(i, dp)/n
            phases(i) = EXP(CMPLX(0.0_dp, 0.37_dp*i, dp))
            DO j = i + 1, n
               change(i, j) = 0.1_dp*SIN(REAL(i + 3*j, dp))
            END DO
         END DO
         DO negative = 0, n
            diagonal = 0.0_dp
            DO i = 1, n
               diagonal(i, i) = 1.0_dp + REAL(i, dp)/n
               IF (i <= negative) diagonal(i, i) = -diagonal(i, i)
            END DO
            matrix = MATMUL(TRANSPOSE(change), MATMUL(diagonal, change))
            DO i = 1, n
               DO j = 1, n
                  cmatrix(i, j) = CONJG(phases(i))*matrix(i, j)*phases(j)
               END DO
            END DO
            DO exponent = -150, 150, 150
               scale = 10.0_dp**exponent
               CALL localizer_inertia(scale*matrix, p, m, info)
               IF (info /= 0 .OR. p /= n - negative .OR. m /= negative) ERROR STOP 'Real LDL signature failed'
               CALL localizer_inertia(scale*cmatrix, p, m, info)
               IF (info /= 0 .OR. p /= n - negative .OR. m /= negative) ERROR STOP 'Complex LDL signature failed'
               CALL localizer_inertia(CMPLX(scale*matrix, 0.0_dp, dp), p, m, info)
               IF (info /= 0 .OR. p /= n - negative .OR. m /= negative) ERROR STOP 'Real LDL dispatch failed'
            END DO
         END DO
         DEALLOCATE (diagonal, change, matrix, phases, cmatrix)
      END DO
      pair = RESHAPE([0.0_dp, 2.0_dp, 2.0_dp, 0.0_dp], [2, 2])
      CALL localizer_inertia(pair, p, m, info)
      IF (info /= 0 .OR. p /= 1 .OR. m /= 1) ERROR STOP 'Zero-diagonal real 2x2 pivot failed'
      pair(1, 2) = 0.0_dp
      CALL localizer_inertia(pair, p, m, info)
      IF (info /= -1) ERROR STOP 'Nonsymmetric real matrix accepted'
      pair = 0.0_dp
      CALL localizer_inertia(pair, p, m, info)
      IF (info /= 1) ERROR STOP 'Singular real matrix accepted'
   END SUBROUTINE real_inertia_tests

! **************************************************************************************************
!> \brief Generalized gaps have an analytic reference independent of LDL pivot magnitudes.
! **************************************************************************************************
   SUBROUTINE gap_tests()
      COMPLEX(dp)                                        :: change(2, 2), h(2, 2), hc(2, 2), &
                                                            r(2, 2, 2), s(2, 2), sc(2, 2)
      INTEGER                                            :: i, index, status
      REAL(dp)                                           :: crossing(5), expected, gap, ratio, trace

      h = RESHAPE([1.0_dp, 20.0_dp, 20.0_dp, 399.0_dp], [2, 2])
      s = 0.0_dp
      s(1, 1) = 2.0_dp
      s(2, 2) = 2.0_dp
      r = 0.0_dp
      trace = REAL(h(1, 1) + h(2, 2), dp)
      expected = 1.0_dp/(trace + SQRT(trace**2 + 4.0_dp))
      CALL localizer_evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - expected) > 1.e-12_dp) ERROR STOP 'Analytic LDL/gap separation failed'
      change = RESHAPE([0.5_dp, 0.0_dp, 0.2_dp, 3.0_dp], [2, 2])
      hc = MATMUL(TRANSPOSE(change), MATMUL(h, change))
      sc = MATMUL(TRANSPOSE(change), MATMUL(s, change))
      CALL localizer_evaluate(hc, sc, r, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. ABS(gap - expected) > 1.e-12_dp) ERROR STOP 'Metric LDL/gap separation failed'
      crossing = [-1.e-6_dp, -1.e-12_dp, 0.0_dp, 1.e-12_dp, 1.e-6_dp]
      h = 0.0_dp
      h(2, 2) = 5.0_dp
      DO i = 1, SIZE(crossing)
         h(1, 1) = crossing(i)
         CALL localizer_evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                 1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
         IF (i == 1 .OR. i == 5) THEN
            IF (status /= 0 .OR. ABS(gap - 5.e-7_dp) > 1.e-14_dp) ERROR STOP 'Resolved gap crossing failed'
         ELSE
            IF (status /= 1) ERROR STOP 'Unresolved gap assigned an index'
         END IF
      END DO
   END SUBROUTINE gap_tests

! **************************************************************************************************
!> \brief Independent Pfaffian polynomials distinguish signs despite identical determinants.
! **************************************************************************************************
   SUBROUTINE pfaffian_orientation_tests()
      INTEGER                                            :: exponent, i, j, n, pf_sign, sample, &
                                                            status
      INTEGER, ALLOCATABLE                               :: permutation(:)
      REAL(dp)                                           :: exact, scale
      REAL(dp), ALLOCATABLE                              :: a(:, :), change(:, :), transformed(:, :)

      DO n = 2, 8, 2
         ALLOCATE (a(n, n), change(n, n), transformed(n, n), permutation(n))
         DO sample = 1, 8
            a = 0.0_dp
            DO j = 2, n
               DO i = 1, j - 1
                  a(i, j) = SIN(REAL(11*i + 17*j + sample, dp))
               END DO
            END DO
            DO i = 1, n, 2
               a(i, i + 1) = a(i, i + 1) + 3.0_dp
            END DO
            a = a - TRANSPOSE(a)
            exact = pfaffian_polynomial(a)
            IF (ABS(exact) < 1.e-6_dp) ERROR STOP 'Ill-conditioned Pfaffian fixture'
            permutation = [(i, i=1, n)]
            permutation(:2) = [2, 1]
            DO exponent = -150, 150, 150
               scale = 10.0_dp**exponent
               CALL localizer_pfaffian_sign(scale*a, pf_sign, status)
               IF (status /= 0 .OR. pf_sign /= INT(SIGN(1.0_dp, exact))) ERROR STOP 'Pfaffian polynomial sign failed'
               CALL localizer_pfaffian_sign(-scale*a, pf_sign, status)
               IF (status /= 0 .OR. pf_sign /= (-1)**(n/2)*INT(SIGN(1.0_dp, exact))) ERROR STOP 'Negative Pfaffian scale failed'
               CALL localizer_pfaffian_sign(scale*a(permutation, permutation), pf_sign, status)
               IF (status /= 0 .OR. pf_sign /= -INT(SIGN(1.0_dp, exact))) ERROR STOP 'Equal determinant/opposite Pfaffian failed'
            END DO
            change = 0.0_dp
            DO i = 1, n
               change(i, i) = 0.8_dp + 0.1_dp*i
               IF (i < n) change(i, i + 1) = 0.1_dp
            END DO
            change(1, :) = -change(1, :)
            transformed = MATMUL(change, MATMUL(a, TRANSPOSE(change)))
            CALL localizer_pfaffian_sign(transformed, pf_sign, status)
            IF (status /= 0 .OR. pf_sign /= -INT(SIGN(1.0_dp, exact))) ERROR STOP 'Pfaffian congruence orientation failed'
         END DO
         DEALLOCATE (a, change, transformed, permutation)
      END DO
      ALLOCATE (a(4, 4))
      a = 0.0_dp
      a(1, 2) = 1.0_dp
      a(1, 3) = 1.0_dp
      a(2, 4) = 1.0_dp
      a(3, 4) = 1.0_dp
      a = a - TRANSPOSE(a)
      CALL localizer_pfaffian_sign(a, pf_sign, status)
      IF (status /= -2 .OR. pf_sign /= 0) ERROR STOP 'Pfaffian polynomial cancellation ignored'
   END SUBROUTINE pfaffian_orientation_tests

! **************************************************************************************************
!> \brief Small independent recursive Pfaffian oracle, without pivoting or elimination.
!> \param a Even-order skew matrix, at most eight rows in these tests
!> \return Pfaffian from the expansion along the first row
! **************************************************************************************************
   RECURSIVE FUNCTION pfaffian_polynomial(a) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a(:, :)
      REAL(dp)                                           :: value

      INTEGER                                            :: i, j, n, keep(SIZE(a, 1) - 2)

      n = SIZE(a, 1)
      IF (n == 2) THEN
         value = a(1, 2)
         RETURN
      END IF
      value = 0.0_dp
      DO j = 2, n
         keep = PACK([(i, i=2, n)], [(i /= j, i=2, n)])
         value = value + REAL((-1)**j, dp)*a(1, j)*pfaffian_polynomial(a(keep, keep))
      END DO
   END FUNCTION pfaffian_polynomial

! **************************************************************************************************
!> \brief Apply uniform numerical tolerances in tests.
!> \param h ...
!> \param s ...
!> \param r ...
!> \param point ...
!> \param energy ...
!> \param kappa ...
!> \param aii ...
!> \param invariant ...
!> \param gap ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE evaluate(h, s, r, point, energy, kappa, aii, invariant, gap, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: h(:, :), s(:, :), r(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: point(2), energy, kappa
      LOGICAL, INTENT(IN)                                :: aii
      INTEGER, INTENT(OUT)                               :: invariant
      REAL(KIND=dp), INTENT(OUT)                         :: gap
      INTEGER, INTENT(OUT)                               :: status

      CALL localizer_evaluate(h, s, r, point, energy, kappa, aii, 1.e-10_dp, 1.e-10_dp, &
                              1.e-9_dp, invariant, gap, ratio, status)
   END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief Open square Qi-Wu-Zhang model: sin(kx)sx+sin(ky)sy+(m+cos(kx)+cos(ky))sz.
!> \param length Number of sites along each direction
!> \param mass Onsite mass; -1 nontrivial, -3 trivial in the infinite model
!> \param h Hamiltonian in site-major, two-orbital order
!> \param s Identity metric
!> \param r Commuting diagonal positions centered on the flake
!> \param torus Optional periodic boundary hoppings
! **************************************************************************************************
   SUBROUTINE qwz_model(length, mass, h, s, r, torus)
      INTEGER, INTENT(IN)                                :: length
      REAL(KIND=dp), INTENT(IN)                          :: mass
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: h(:, :), s(:, :), r(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: torus

      COMPLEX(KIND=dp)                                   :: tx(2, 2), ty(2, 2)
      INTEGER                                            :: i, j, n, orbital, x, y

      n = 2*length**2
      ALLOCATE (h(n, n), s(n, n), r(n, n, 2))
      h = 0.0_dp
      s = 0.0_dp
      r = 0.0_dp
      tx = 0.0_dp
      tx(1, 1) = 0.5_dp
      tx(2, 2) = -0.5_dp
      tx(1, 2) = CMPLX(0.0_dp, -0.5_dp, dp)
      tx(2, 1) = tx(1, 2)
      ty = 0.0_dp
      ty(1, 1) = 0.5_dp
      ty(2, 2) = -0.5_dp
      ty(1, 2) = -0.5_dp
      ty(2, 1) = 0.5_dp
      DO y = 1, length
         DO x = 1, length
            i = 2*((y - 1)*length + x) - 1
            h(i, i) = mass
            h(i + 1, i + 1) = -mass
            DO orbital = 0, 1
               s(i + orbital, i + orbital) = 1.0_dp
               r(i + orbital, i + orbital, 1) = REAL(x, dp) - REAL(length + 1, dp)/2.0_dp
               r(i + orbital, i + orbital, 2) = REAL(y, dp) - REAL(length + 1, dp)/2.0_dp
            END DO
            IF (x < length) THEN
               j = i + 2
               h(i:i + 1, j:j + 1) = tx
               h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(tx))
            END IF
            IF (y < length) THEN
               j = i + 2*length
               h(i:i + 1, j:j + 1) = ty
               h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(ty))
            END IF
            IF (PRESENT(torus)) THEN
               IF (torus .AND. x == length) THEN
                  j = 2*(y - 1)*length + 1
                  h(i:i + 1, j:j + 1) = tx
                  h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(tx))
               END IF
               IF (torus .AND. y == length) THEN
                  j = 2*x - 1
                  h(i:i + 1, j:j + 1) = ty
                  h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(ty))
               END IF
            END IF
         END DO
      END DO
   END SUBROUTINE qwz_model

END PROGRAM spectral_localizer_unittest

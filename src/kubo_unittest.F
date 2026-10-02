!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE kubo_test_cases
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: check_projected_current, check_bloch_current, check_hall_response, check_real_space_chern
CONTAINS

! **************************************************************************************************
!> \brief Validate projected currents and their AO-metric covariance.
! **************************************************************************************************
   SUBROUTINE check_projected_current()
      USE kinds, ONLY: dp
      USE kubo_projected, ONLY: projected_current
      COMPLEX(dp) :: b(4, 4), h(4, 4), h2(4, 4), j(4, 4, 3), j2(4, 4, 3), norm(3, 3), norm2(3, 3), &
         s(4, 4), s2(4, 4), x(4, 4, 3), x2(4, 4, 3)
      INTEGER                                            :: a, c, i, status
      REAL(dp)                                           :: e(4), e2(4)

      h = 0.0_dp
      s = 0.0_dp
      x = 0.0_dp
      DO i = 1, 4
         s(i, i) = 1.0_dp
      END DO
      ! Kramers-degenerate spinor model with imaginary SOC matrix elements.
      h(2, 2) = 0.8_dp
      h(4, 4) = 0.8_dp
      h(1, 2) = CMPLX(0.0_dp, 0.3_dp, dp)
      h(2, 1) = CONJG(h(1, 2))
      h(3, 4) = CONJG(h(1, 2))
      h(4, 3) = h(1, 2)
      x(1, 2, 1) = 1.0_dp
      x(2, 1, 1) = 1.0_dp
      x(3, 4, 1) = 1.0_dp
      x(4, 3, 1) = 1.0_dp
      x(1, 1, 2) = -0.5_dp
      x(2, 2, 2) = 0.5_dp
      x(3, 3, 2) = -0.5_dp
      x(4, 4, 2) = 0.5_dp
      x(:, :, 3) = 2*x(:, :, 1) - x(:, :, 2)
      CALL projected_current(h, s, x, e, j, status)
      IF (status /= 0) ERROR STOP 'Spinor current failed'
      IF (MAXVAL(ABS(e - [-0.1_dp, -0.1_dp, 0.9_dp, 0.9_dp])) > 1.e-12_dp) ERROR STOP 'Spinor spectrum failed'
      DO a = 1, 3
         IF (MAXVAL(ABS(j(:, :, a) + CONJG(TRANSPOSE(j(:, :, a))))) > 1.e-12_dp) ERROR STOP 'Current is not anti-Hermitian'
         DO c = 1, 3
            norm(a, c) = SUM(CONJG(j(:, :, a))*j(:, :, c))
         END DO
      END DO
      IF (ABS(norm(1, 1) - 4.0_dp) > 1.e-12_dp) ERROR STOP 'Analytic current strength failed'
      ! Complex nonorthogonal basis mixing spin and orbitals, including a degenerate subspace.
      b = s
      b(1, 1) = 1.7_dp
      b(2, 2) = 0.8_dp
      b(3, 1) = CMPLX(0.2_dp, 0.4_dp, dp)
      b(1, 4) = CMPLX(-0.3_dp, 0.1_dp, dp)
      b(2, 3) = 0.4_dp
      h2 = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(h, b))
      s2 = MATMUL(CONJG(TRANSPOSE(b)), b)
      DO a = 1, 3
         x2(:, :, a) = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(x(:, :, a), b)) + REAL(a, dp)*s2
      END DO
      CALL projected_current(h2, s2, x2, e2, j2, status)
      IF (status /= 0 .OR. MAXVAL(ABS(e - e2)) > 1.e-12_dp) ERROR STOP 'Nonorthogonal current spectrum failed'
      DO a = 1, 3
         DO c = 1, 3
            norm2(a, c) = SUM(CONJG(j2(:, :, a))*j2(:, :, c))
         END DO
      END DO
      IF (MAXVAL(ABS(norm - norm2)) > 1.e-11_dp) ERROR STOP 'Basis/origin/gauge invariance failed'
      s(1, 1) = -1.0_dp
      CALL projected_current(h, s, x, e, j, status)
      IF (status /= -2) ERROR STOP 'Indefinite overlap accepted'
      s(1, 1) = 1.0_dp
      x(1, 2, 1) = 2.0_dp
      CALL projected_current(h, s, x, e, j, status)
      IF (status /= -1) ERROR STOP 'Non-Hermitian position accepted'
      PRINT *, 'Projected spinor current, metric and gauge tests passed.'
   END SUBROUTINE check_projected_current

! **************************************************************************************************
!> \brief Validate Bloch currents, degenerate blocks and k-dependent AO gauges.
! **************************************************************************************************
   SUBROUTINE check_bloch_current()
      USE kinds, ONLY: dp
      USE kubo_bloch, ONLY: bloch_current, &
                            bloch_current_from_frame, &
                            fermi_slope
      USE kubo_projected, ONLY: projected_current
      COMPLEX(dp), PARAMETER                             :: iu = (0.0_dp, 1.0_dp)

      COMPLEX(dp) :: b(2, 2), d(2, 2, 3), d2(2, 2, 3), db(2, 2), dh(2, 2, 3), dh2(2, 2, 3), &
         ds(2, 2, 3), ds2(2, 2, 3), expected(2, 2), frame(2, 2), h(2, 2), h2(2, 2), j(2, 2, 3), &
         j2(2, 2, 3), s(2, 2), s2(2, 2)
      INTEGER                                            :: i, ik, status
      REAL(dp)                                           :: e(2), e2(2), k, residuals(2), root, &
                                                            shift, speed, strength, strength2

      s = 0.0_dp
      s(1, 1) = 1.0_dp
      s(2, 2) = 1.0_dp
      DO ik = -7, 7
         k = REAL(ik, dp)/5
         h = k*k*s
         h(1, 1) = h(1, 1) + 0.3_dp
         h(2, 2) = h(2, 2) - 0.3_dp
         h(1, 2) = 0.4_dp*SIN(k) - iu*0.2_dp*COS(k)
         h(2, 1) = CONJG(h(1, 2))
         dh = 0.0_dp
         dh(:, :, 1) = 2*k*s
         dh(1, 2, 1) = 0.4_dp*COS(k) + iu*0.2_dp*SIN(k)
         dh(2, 1, 1) = CONJG(dh(1, 2, 1))
         ds = 0.0_dp
         d = 0.0_dp
         CALL bloch_current(h, s, dh, ds, d, e, j, status, frame)
         IF (status /= 0) ERROR STOP 'Bloch current failed'
         CALL bloch_current_from_frame(h, s, dh, ds, d, e, frame, j2, 1.e-9_dp, residuals, status)
         IF (status /= 0 .OR. MAXVAL(ABS(j2 - j)) > 1.e-12_dp) ERROR STOP 'Reused eigenframe failed'
         CALL bloch_current_from_frame(h, s, dh, ds, d, e + 0.1_dp, frame, j2, 1.e-9_dp, residuals, status)
         IF (status /= -4) ERROR STOP 'Wrong frame eigenvalues accepted'
         CALL bloch_current_from_frame(h, s, dh, ds, d, e, 1.1_dp*frame, j2, 1.e-9_dp, residuals, status)
         IF (status /= -4) ERROR STOP 'Nonorthonormal frame accepted'
         root = SQRT(0.09_dp + 0.16_dp*SIN(k)**2 + 0.04_dp*COS(k)**2)
         IF (MAXVAL(ABS(e - [k*k - root, k*k + root])) > 1.e-12_dp) ERROR STOP 'Analytic bands failed'
         speed = 0.12_dp*SIN(k)*COS(k)/root
         IF (ABS(iu*j(1, 1, 1) - (2*k - speed)) > 1.e-12_dp) ERROR STOP 'Lower band velocity failed'
         IF (ABS(iu*j(2, 2, 1) - (2*k + speed)) > 1.e-12_dp) ERROR STOP 'Upper band velocity failed'
         IF (MAXVAL(ABS(j(:, :, 1) + CONJG(TRANSPOSE(j(:, :, 1))))) > 1.e-12_dp) THEN
            ERROR STOP 'Bloch current not anti-Hermitian'
         END IF

         ! A k-dependent nonunitary AO change requires both overlap and connection derivatives.
         db = 0.0_dp
         db(1, 1) = 0.13_dp
         db(2, 1) = CMPLX(0.17_dp, 0.23_dp, dp)
         db(1, 2) = CMPLX(-0.09_dp, 0.11_dp, dp)
         b = s + k*db
         h2 = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(h, b))
         s2 = MATMUL(CONJG(TRANSPOSE(b)), b)
         dh2 = 0.0_dp
         ds2 = 0.0_dp
         d2 = 0.0_dp
         dh2(:, :, 1) = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(dh(:, :, 1), b)) + &
                        MATMUL(CONJG(TRANSPOSE(db)), MATMUL(h, b)) + MATMUL(CONJG(TRANSPOSE(b)), MATMUL(h, db))
         ds2(:, :, 1) = MATMUL(CONJG(TRANSPOSE(db)), b) + MATMUL(CONJG(TRANSPOSE(b)), db)
         d2(:, :, 1) = -iu*MATMUL(CONJG(TRANSPOSE(db)), b) + 7*s2
         CALL bloch_current(h2, s2, dh2, ds2, d2, e2, j2, status)
         IF (status /= 0 .OR. MAXVAL(ABS(e - e2)) > 1.e-12_dp) ERROR STOP 'Bloch basis covariance failed'
         DO i = 1, 2
            IF (ABS(j(i, i, 1) - j2(i, i, 1)) > 1.e-11_dp) ERROR STOP 'Band velocity basis dependence'
         END DO
         IF (MAXVAL(ABS(ABS(j) - ABS(j2))) > 1.e-11_dp) ERROR STOP 'Interband velocity gauge dependence'
         ! An energy-zero shift must include the k derivative of the metric.
         shift = 4.2_dp
         CALL bloch_current(h2 + shift*s2, s2, dh2 + shift*ds2, ds2, d2, e2, j2, status)
         IF (status /= 0 .OR. MAXVAL(ABS(e + shift - e2)) > 1.e-11_dp) ERROR STOP 'Energy shift failed'
         IF (MAXVAL(ABS(ABS(j) - ABS(j2))) > 1.e-10_dp) ERROR STOP 'Velocity energy-zero dependence'
         d2(:, :, 1) = 0.0_dp
         CALL bloch_current(h2, s2, dh2, ds2, d2, e2, j2, status)
         IF (status /= -3) ERROR STOP 'Inconsistent connection accepted'
      END DO

      dh = 0.0_dp
      ds = 0.0_dp
      d = 0.0_dp
      d(1, 2, 1) = iu
      d(2, 1, 1) = -iu
      CALL bloch_current(h, s, dh, ds, d, e, j, status)
      IF (status /= 0) ERROR STOP 'Finite Bloch limit failed'
      CALL projected_current(h, s, d, e2, j2, status)
      IF (status /= 0 .OR. MAXVAL(ABS(j - j2)) > 1.e-12_dp) ERROR STOP 'Finite commutator limit failed'

      ! Every matrix element within an exactly degenerate block has the same Fermi weight.
      h = 0.0_dp
      dh = 0.0_dp
      d = 0.0_dp
      dh(:, :, 1) = s
      dh(1, 2, 1) = iu
      dh(2, 1, 1) = -iu
      CALL bloch_current(h, s, dh, ds, d, e, j, status)
      IF (status /= 0) ERROR STOP 'Degenerate Bloch block failed'
      strength = SUM(ABS(j(:, :, 1))**2)*(-fermi_slope(0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp))
      IF (ABS(strength - 10.0_dp) > 1.e-12_dp) ERROR STOP 'Intraband Drude weight failed'
      b = RESHAPE([iu, (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), iu], [2, 2])/SQRT(2.0_dp)
      expected = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(dh(:, :, 1), b))
      dh(:, :, 1) = expected
      CALL bloch_current(h, s, dh, ds, d, e, j, status)
      strength2 = SUM(ABS(j(:, :, 1))**2)*(-fermi_slope(0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp))
      IF (status /= 0 .OR. ABS(strength - strength2) > 1.e-12_dp) ERROR STOP 'Degenerate response gauge failed'
      IF (ABS(fermi_slope(1.e-9_dp, -1.e-9_dp, 0.0_dp, 0.1_dp) + 2.5_dp) > 1.e-12_dp) THEN
         ERROR STOP 'Fermi near-degeneracy limit failed'
      END IF
      IF (ABS(fermi_slope(0.1_dp, -0.1_dp, 0.0_dp, 0.1_dp) + TANH(0.5_dp)/0.2_dp) > 1.e-12_dp) THEN
         ERROR STOP 'Fermi divided difference failed'
      END IF
      PRINT *, 'Bloch metric/connection, band velocity, gauge, degeneracy and Drude tests passed.'
   END SUBROUTINE check_bloch_current

! **************************************************************************************************
!> \brief Validate Hall signs against Liouvillian and Berry-curvature references.
! **************************************************************************************************
   SUBROUTINE check_hall_response()
      USE ieee_arithmetic, ONLY: ieee_quiet_nan, &
                                 ieee_value
      USE kinds, ONLY: dp
      USE kubo_bloch, ONLY: bloch_current, &
                            fermi_slope
      USE kubo_hall, ONLY: hall_response

      IMPLICIT NONE
      COMPLEX(dp), PARAMETER :: iu = (0.0_dp, 1.0_dp)
      REAL(dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp
      REAL(dp) :: energies(4), hall(3, 3), check(3, 3), reference(3, 3), rotation(3, 3)
      REAL(dp) :: analytic, error, gauge_error, mass, kx, ky, radius, omega, integrated, curvature, gamma
      REAL(dp) :: e2(2), masses(3), expected(3), unit(3)
      COMPLEX(dp) :: current(4, 4, 3), changed(4, 4, 3), u(4, 4)
      COMPLEX(dp) :: h2(2, 2), s2(2, 2), dh2(2, 2, 3), ds2(2, 2, 3), d2(2, 2, 3), j2(2, 2, 3)
      INTEGER :: a, b, i, j, imass, ix, iy, status
      INTEGER, PARAMETER :: mesh = 48

      ! A two-level response fixes the sign, occupation factor and relaxation denominator.
      j2 = 0.0_dp
      j2(1, 2, 1) = -iu
      j2(2, 1, 1) = -iu
      j2(1, 2, 2) = -1.0_dp
      j2(2, 1, 2) = 1.0_dp
      CALL hall_response([-0.7_dp, 0.7_dp], j2, 0.0_dp, 0.2_dp, 0.3_dp, hall, status)
      analytic = 2.0_dp*TANH(0.7_dp/0.4_dp)/(1.4_dp**2 + 0.3_dp**2)
      IF (status /= 0 .OR. ABS(hall(1, 2) - analytic) > 1.e-13_dp) ERROR STOP 'Two-level Hall sign or scale'
      IF (MAXVAL(ABS(hall + TRANSPOSE(hall))) /= 0.0_dp) ERROR STOP 'Hall not antisymmetric'
      CALL hall_response([0.0_dp, 0.0_dp], j2, 0.0_dp, 0.2_dp, 0.3_dp, hall, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall)) /= 0.0_dp) ERROR STOP 'Degenerate Hall block'
      DO i = 2, 14, 2
         radius = 10.0_dp**(-i)
         CALL hall_response([-radius, radius], j2, 0.0_dp, 0.2_dp, 0.3_dp, hall, status)
         analytic = 2.0_dp*TANH(radius/0.4_dp)/(4*radius**2 + 0.3_dp**2)
         IF (status /= 0 .OR. ABS(hall(1, 2)/analytic - 1.0_dp) > 1.e-12_dp) THEN
            ERROR STOP 'Near-degenerate Hall limit'
         END IF
      END DO

      energies = [-0.9_dp, -0.9_dp, 0.8_dp, 0.8_dp]
      DO a = 1, 3
         DO j = 1, 4
            DO i = 1, 4
               current(i, j, a) = CMPLX(SIN(REAL(i + 3*j + a, dp)), COS(REAL(2*i - j + a, dp)), dp)
            END DO
         END DO
         current(:, :, a) = current(:, :, a) - CONJG(TRANSPOSE(current(:, :, a)))
      END DO
      CALL hall_response(energies, current, 0.1_dp, 0.15_dp, 0.07_dp, hall, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall)) < 0.1_dp) ERROR STOP 'Nontrivial four-band response'
      CALL liouville_reference(energies, current, 0.1_dp, 0.15_dp, 0.07_dp, reference)
      error = MAXVAL(ABS(hall - 0.5_dp*(reference - TRANSPOSE(reference))))
      IF (error > 2.e-12_dp) ERROR STOP 'Full Liouvillian reference mismatch'
      PRINT '(A,ES14.6)', 'Hall Liouvillian residual: ', error
      DO a = 1, 3
         unit = 0.0_dp
         unit(a) = 1.0_dp
         IF (DOT_PRODUCT(unit, MATMUL(reference, unit)) < -1.e-12_dp) ERROR STOP 'Dissipation sign'
      END DO

      u = 0.0_dp
      u(1:2, 1:2) = RESHAPE([iu, (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), iu], [2, 2])/SQRT(2.0_dp)
      u(3:4, 3:4) = RESHAPE([(0.6_dp, 0.0_dp), (0.0_dp, 0.8_dp), (0.0_dp, 0.8_dp), &
                             (0.6_dp, 0.0_dp)], [2, 2])
      DO a = 1, 3
         changed(:, :, a) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(current(:, :, a), u))
      END DO
      CALL hall_response(energies, changed, 0.1_dp, 0.15_dp, 0.07_dp, check, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall - check)) > 2.e-12_dp) ERROR STOP 'Degenerate gauge covariance'
      CALL hall_response(energies + 2.3_dp, current, 2.4_dp, 0.15_dp, 0.07_dp, check, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall - check)) > 2.e-12_dp) ERROR STOP 'Energy zero covariance'
      changed = current([3, 1, 4, 2], [3, 1, 4, 2], :)
      CALL hall_response(energies([3, 1, 4, 2]), changed, 0.1_dp, 0.15_dp, 0.07_dp, check, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall - check)) > 2.e-12_dp) ERROR STOP 'Band permutation covariance'
      CALL hall_response(energies, CONJG(current), 0.1_dp, 0.15_dp, 0.07_dp, check, status)
      IF (status /= 0 .OR. MAXVAL(ABS(hall + check)) > 2.e-12_dp) ERROR STOP 'Time-reversed Hall sign'
      rotation = RESHAPE([0.6_dp, 0.8_dp, 0.0_dp, -0.8_dp, 0.6_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -1.0_dp], [3, 3])
      changed = 0.0_dp
      DO a = 1, 3
         DO b = 1, 3
            changed(:, :, a) = changed(:, :, a) + rotation(a, b)*current(:, :, b)
         END DO
      END DO
      CALL hall_response(energies, changed, 0.1_dp, 0.15_dp, 0.07_dp, check, status)
      reference = MATMUL(rotation, MATMUL(hall, TRANSPOSE(rotation)))
      IF (status /= 0 .OR. MAXVAL(ABS(reference - check)) > 2.e-12_dp) ERROR STOP 'Cartesian tensor covariance'

      CALL hall_response(energies, current, 0.0_dp, 0.0_dp, 0.1_dp, check, status)
      IF (status == 0) ERROR STOP 'Zero temperature accepted'
      CALL hall_response(energies, current, 0.0_dp, 0.1_dp, 0.0_dp, check, status)
      IF (status == 0) ERROR STOP 'Zero relaxation accepted'
      CALL hall_response(energies, current, ieee_value(0.0_dp, ieee_quiet_nan), 0.1_dp, 0.1_dp, check, status)
      IF (status == 0) ERROR STOP 'NaN chemical potential accepted'
      changed = current
      changed(1, 2, 1) = changed(1, 2, 1) + 0.1_dp
      CALL hall_response(energies, changed, 0.0_dp, 0.1_dp, 0.1_dp, check, status)
      IF (status == 0) ERROR STOP 'Non-anti-Hermitian current accepted'

      ! Two-band lattice model: compare every point against analytic Berry curvature,
      ! then integrate its two opposite Chern phases and a trivial phase.
      masses = [-1.0_dp, 1.0_dp, 3.0_dp]
      expected = [-1.0_dp, 1.0_dp, 0.0_dp]
      s2 = 0.0_dp
      s2(1, 1) = 1.0_dp
      s2(2, 2) = 1.0_dp
      ds2 = 0.0_dp
      d2 = 0.0_dp
      gamma = 1.e-3_dp
      DO imass = 1, 3
         mass = masses(imass)
         integrated = 0.0_dp
         curvature = 0.0_dp
         error = 0.0_dp
         gauge_error = 0.0_dp
         DO iy = 0, mesh - 1
            ky = 2*pi*REAL(iy, dp)/mesh
            DO ix = 0, mesh - 1
               kx = 2*pi*REAL(ix, dp)/mesh
               h2(1, 1) = mass + COS(kx) + COS(ky)
               h2(2, 2) = -h2(1, 1)
               h2(1, 2) = SIN(kx) - iu*SIN(ky)
               h2(2, 1) = CONJG(h2(1, 2))
               dh2 = 0.0_dp
               dh2(1, 1, 1) = -SIN(kx)
               dh2(2, 2, 1) = SIN(kx)
               dh2(1, 2, 1) = COS(kx)
               dh2(2, 1, 1) = COS(kx)
               dh2(1, 1, 2) = -SIN(ky)
               dh2(2, 2, 2) = SIN(ky)
               dh2(1, 2, 2) = -iu*COS(ky)
               dh2(2, 1, 2) = iu*COS(ky)
               CALL bloch_current(h2, s2, dh2, ds2, d2, e2, j2, status)
               IF (status /= 0) ERROR STOP 'Chern-model Bloch current'
               CALL hall_response(e2, j2, 0.0_dp, 0.01_dp, gamma, check, status)
               IF (status /= 0) ERROR STOP 'Chern-model Hall kernel'
               radius = SQRT(SIN(kx)**2 + SIN(ky)**2 + (mass + COS(kx) + COS(ky))**2)
               omega = (COS(kx) + COS(ky) + mass*COS(kx)*COS(ky))/(2*radius**3)
               analytic = -omega*4*radius**2/(4*radius**2 + gamma**2)*TANH(radius/0.02_dp)
               error = MAX(error, ABS(check(1, 2) - analytic))
               CALL ao_gauge_reference(kx, ky, h2, dh2, gamma, reference)
               gauge_error = MAX(gauge_error, MAXVAL(ABS(check - reference)))
               integrated = integrated + check(1, 2)*2*pi/mesh**2
               curvature = curvature + omega*2*pi/mesh**2
            END DO
         END DO
         IF (error > 1.e-12_dp) ERROR STOP 'Analytic Berry-curvature response mismatch'
         IF (ABS(integrated - expected(imass)) > 1.e-6_dp) ERROR STOP 'Chern Hall limit'
         IF (ABS(curvature + expected(imass)) > 1.e-10_dp) ERROR STOP 'Independent Chern integral'
         IF (gauge_error > 1.e-12_dp) ERROR STOP 'Nonorthogonal AO Hall covariance'
         PRINT '(A,F5.1,3ES18.10)', 'Hall Chern model mass, sigma, Chern, residual: ', mass, integrated, curvature, error
         PRINT '(A,ES14.6)', 'Hall k-dependent nonorthogonal AO residual: ', gauge_error
      END DO
      PRINT *, 'Hall sign, Liouvillian, gauge, time reversal, degeneracy and Chern tests passed.'
   CONTAINS
! **************************************************************************************************
!> \brief Transform the lattice model to a k-dependent nonorthogonal AO frame.
!> \param kx First wave-vector component
!> \param ky Second wave-vector component
!> \param h Orthogonal Hamiltonian
!> \param dh Its Cartesian derivatives
!> \param gamma Relaxation energy
!> \param hall Hall response in the transformed AO frame
! **************************************************************************************************
      SUBROUTINE ao_gauge_reference(kx, ky, h, dh, gamma, hall)
      REAL(dp), INTENT(IN)                               :: kx, ky
      COMPLEX(dp), INTENT(IN)                            :: h(2, 2), dh(2, 2, 3)
      REAL(dp), INTENT(IN)                               :: gamma
      REAL(dp), INTENT(OUT)                              :: hall(3, 3)

      COMPLEX(dp)                                        :: basis(2, 2), db(2, 2, 3), dhx(2, 2, 3), &
                                                            dsx(2, 2, 3), dx(2, 2, 3), hx(2, 2), &
                                                            jx(2, 2, 3), sx(2, 2)
      INTEGER                                            :: a, status
      REAL(dp)                                           :: e(2)

         basis = 0.0_dp
         basis(1, 1) = 1.1_dp + 0.1_dp*COS(kx)
         basis(2, 2) = 0.9_dp + 0.1_dp*SIN(ky)
         basis(1, 2) = 0.13_dp*EXP(iu*(kx + ky))
         basis(2, 1) = 0.07_dp*EXP(iu*(2*kx - ky))
         db = 0.0_dp
         db(1, 1, 1) = -0.1_dp*SIN(kx)
         db(2, 2, 2) = 0.1_dp*COS(ky)
         db(1, 2, 1:2) = iu*basis(1, 2)
         db(2, 1, 1) = 2*iu*basis(2, 1)
         db(2, 1, 2) = -iu*basis(2, 1)
         hx = MATMUL(CONJG(TRANSPOSE(basis)), MATMUL(h, basis))
         sx = MATMUL(CONJG(TRANSPOSE(basis)), basis)
         DO a = 1, 3
            dhx(:, :, a) = MATMUL(CONJG(TRANSPOSE(basis)), MATMUL(dh(:, :, a), basis)) + &
                           MATMUL(CONJG(TRANSPOSE(db(:, :, a))), MATMUL(h, basis)) + &
                           MATMUL(CONJG(TRANSPOSE(basis)), MATMUL(h, db(:, :, a)))
            dsx(:, :, a) = MATMUL(CONJG(TRANSPOSE(db(:, :, a))), basis) + &
                           MATMUL(CONJG(TRANSPOSE(basis)), db(:, :, a))
            dx(:, :, a) = -iu*MATMUL(CONJG(TRANSPOSE(db(:, :, a))), basis) + REAL(a, dp)*sx
         END DO
         CALL bloch_current(hx, sx, dhx, dsx, dx, e, jx, status)
         IF (status /= 0) ERROR STOP 'Nonorthogonal AO Chern model'
         CALL hall_response(e, jx, 0.0_dp, 0.01_dp, gamma, hall, status)
         IF (status /= 0) ERROR STOP 'Nonorthogonal AO Hall kernel'
      END SUBROUTINE ao_gauge_reference

! **************************************************************************************************
!> \brief Invert the full complex Liouvillian in a nondiagonal four-level frame.
!> \param e Eigenvalues
!> \param currents Anti-Hermitian currents in the diagonal frame
!> \param mu Chemical potential
!> \param kt Temperature
!> \param gamma Relaxation rate
!> \param sigma Complete real DC tensor without physical prefactors
! **************************************************************************************************
      SUBROUTINE liouville_reference(e, currents, mu, kt, gamma, sigma)
      REAL(dp), INTENT(IN)                               :: e(4)
      COMPLEX(dp), INTENT(IN)                            :: currents(4, 4, 3)
      REAL(dp), INTENT(IN)                               :: mu, kt, gamma
      REAL(dp), INTENT(OUT)                              :: sigma(3, 3)

      COMPLEX(dp)                                        :: basis(4, 4), drive(4, 4), h(4, 4), &
                                                            liouvillian(16, 16), response(4, 4), &
                                                            rhs(16, 3), trace, v(4, 4, 3)
      INTEGER                                            :: a, b, column, i, info, j, k, l, &
                                                            pivot(16), row

         DO j = 1, 4
            DO i = 1, 4
               basis(i, j) = EXP(iu*2*pi*REAL((i - 1)*(j - 1), dp)/4)/2
            END DO
         END DO
         h = 0.0_dp
         DO i = 1, 4
            h(i, i) = e(i)
         END DO
         h = MATMUL(basis, MATMUL(h, CONJG(TRANSPOSE(basis))))
         DO a = 1, 3
            v(:, :, a) = MATMUL(basis, MATMUL(iu*currents(:, :, a), CONJG(TRANSPOSE(basis))))
            DO j = 1, 4
               DO i = 1, 4
                  drive(i, j) = fermi_slope(e(i), e(j), mu, kt)*iu*currents(i, j, a)
               END DO
            END DO
            drive = MATMUL(basis, MATMUL(drive, CONJG(TRANSPOSE(basis))))
            rhs(:, a) = RESHAPE(drive, [16])
         END DO
         liouvillian = 0.0_dp
         DO j = 1, 4
            DO i = 1, 4
               row = i + 4*(j - 1)
               liouvillian(row, row) = gamma
               DO l = 1, 4
                  DO k = 1, 4
                     column = k + 4*(l - 1)
                     IF (j == l) liouvillian(row, column) = liouvillian(row, column) + iu*h(i, k)
                     IF (i == k) liouvillian(row, column) = liouvillian(row, column) - iu*h(l, j)
                  END DO
               END DO
            END DO
         END DO
         CALL zgesv(16, 3, liouvillian, 16, pivot, rhs, 16, info)
         IF (info /= 0) ERROR STOP 'Liouvillian solve failed'
         DO b = 1, 3
            response = RESHAPE(rhs(:, b), [4, 4])
            DO a = 1, 3
               trace = SUM(v(:, :, a)*TRANSPOSE(response))
               sigma(a, b) = -REAL(trace, dp)
            END DO
         END DO
      END SUBROUTINE liouville_reference
   END SUBROUTINE check_hall_response

! **************************************************************************************************
!> \brief Crosscheck the finite-window marker with existing native Wilson analysis and AO gauges.
! **************************************************************************************************
   SUBROUTINE check_real_space_chern()
      USE ieee_arithmetic, ONLY: ieee_quiet_nan, ieee_value
      USE kinds, ONLY: dp
      USE kubo_projected, ONLY: projected_eigensystem
      USE topology_real_space, ONLY: finite_chern_marker
      USE topology_wilson, ONLY: wilson_step, wilson_spectrum, chern_from_wcc
      COMPLEX(dp), PARAMETER :: iu = (0.0_dp, 1.0_dp)
      INTEGER, PARAMETER :: length = 9, n = 2*length**2, mesh = 32
      REAL(dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp
      COMPLEX(dp) :: h(n, n), s(n, n), r(n, n, 2), w(n, n), c(n, n), b(n, n), cc(n, n), &
                     sc(n, n), rc(n, n, 2), wc(n, n), factor(n, n), tx(2, 2), ty(2, 2), &
                     h2(2, 2), s2(2, 2), c2(2, 2), states(2, mesh), PRODUCT(1, 1), overlap(1, 1)
      REAL(dp) :: e(n), e2(2), masses(3), expected(3), value, changed, full_trace, &
                  wcc(1, mesh + 1), centre(1), winding, berry, minimum_sv, kx, ky, residual
      INTEGER :: phase, x, y, i, j, orbital, axis, status, invariant, pivots(n)

      masses = [-1.0_dp, 1.0_dp, 3.0_dp]
      expected = [1.0_dp, -1.0_dp, 0.0_dp]
      tx = RESHAPE([CMPLX(0.5_dp, 0.0_dp, dp), -0.5_dp*iu, -0.5_dp*iu, &
                    CMPLX(-0.5_dp, 0.0_dp, dp)], [2, 2])
      ty = RESHAPE([CMPLX(0.5_dp, 0.0_dp, dp), CMPLX(0.5_dp, 0.0_dp, dp), &
                    CMPLX(-0.5_dp, 0.0_dp, dp), CMPLX(-0.5_dp, 0.0_dp, dp)], [2, 2])
      s2 = 0.0_dp
      s2(1, 1) = 1.0_dp
      s2(2, 2) = 1.0_dp
      DO phase = 1, 3
         h = 0.0_dp
         s = 0.0_dp
         r = 0.0_dp
         w = 0.0_dp
         DO y = 1, length
            DO x = 1, length
               i = 2*((y - 1)*length + x) - 1
               h(i, i) = masses(phase)
               h(i + 1, i + 1) = -masses(phase)
               DO orbital = 0, 1
                  s(i + orbital, i + orbital) = 1.0_dp
                  r(i + orbital, i + orbital, 1) = REAL(x - (length + 1)/2, dp)
                  r(i + orbital, i + orbital, 2) = REAL(y - (length + 1)/2, dp)
                  IF (ABS(x - (length + 1)/2) <= 1 .AND. ABS(y - (length + 1)/2) <= 1) THEN
                     w(i + orbital, i + orbital) = 1.0_dp
                  END IF
               END DO
               IF (x < length) THEN
                  h(i:i + 1, i + 2:i + 3) = tx
                  h(i + 2:i + 3, i:i + 1) = CONJG(TRANSPOSE(tx))
               END IF
               IF (y < length) THEN
                  j = i + 2*length
                  h(i:i + 1, j:j + 1) = ty
                  h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(ty))
               END IF
            END DO
         END DO
         CALL projected_eigensystem(h, s, e, c, status)
         IF (status /= 0) ERROR STOP 'Open Chern model eigensystem'
         CALL finite_chern_marker(c(:, :n/2), s, r, w, 9.0_dp, value, status)
         IF (status /= 0) ERROR STOP 'Finite-window Chern marker'
         CALL finite_chern_marker(c(:, :n/2), s, r, s, REAL(length**2, dp), full_trace, status)
         IF (status /= 0 .OR. ABS(full_trace) > 1.e-11_dp) ERROR STOP 'Open lattice boundary cancellation'

         ! Use the existing native Wilson implementation, not a second Chern integrator.
         DO y = 0, mesh
            ky = 2*pi*REAL(y, dp)/mesh
            DO x = 0, mesh - 1
               kx = 2*pi*REAL(x, dp)/mesh
               h2(1, 1) = masses(phase) + COS(kx) + COS(ky)
               h2(2, 2) = -h2(1, 1)
               h2(1, 2) = SIN(kx) - iu*SIN(ky)
               h2(2, 1) = CONJG(h2(1, 2))
               CALL projected_eigensystem(h2, s2, e2, c2, status)
               IF (status /= 0) ERROR STOP 'Wilson reference eigensystem'
               states(:, x + 1) = c2(:, 1)
            END DO
            product = 1.0_dp
            DO x = 1, mesh
               overlap(1, 1) = DOT_PRODUCT(states(:, x), states(:, MODULO(x, mesh) + 1))
               CALL wilson_step(product, overlap, minimum_sv, status, 1.e-10_dp)
               IF (status /= 0) ERROR STOP 'Native Wilson reference link'
            END DO
            CALL wilson_spectrum(product, centre, berry, status)
            IF (status /= 0) ERROR STOP 'Native Wilson reference spectrum'
            wcc(:, y + 1) = centre
         END DO
         CALL chern_from_wcc(wcc, invariant, winding, status, 1.e-9_dp)
         IF (status /= 0 .OR. ABS(REAL(invariant, dp) - expected(phase)) > 1.e-12_dp) THEN
            ERROR STOP 'Native Wilson reference orientation'
         END IF
         IF (ABS(value - winding) > 0.01_dp) ERROR STOP 'Finite marker versus native Wilson'

         b = 0.0_dp
         DO i = 1, n
            b(i, i) = (1.0_dp + REAL(MODULO(i, 5), dp)/5)*EXP(iu*0.17_dp*i)
            IF (i < n) b(i, i + 1) = CMPLX(0.05_dp, 0.03_dp, dp)
         END DO
         sc = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(s, b))
         wc = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(w, b))
         DO axis = 1, 2
            rc(:, :, axis) = MATMUL(CONJG(TRANSPOSE(b)), MATMUL(r(:, :, axis), b))
         END DO
         factor = b
         cc = c
         CALL zgesv(n, n/2, factor, n, pivots, cc, n, status)
         IF (status /= 0) ERROR STOP 'Occupied AO basis change'
         CALL finite_chern_marker(cc(:, :n/2), sc, rc, wc, 9.0_dp, changed, status)
         residual = ABS(value - changed)
         IF (status /= 0 .OR. residual > 1.e-11_dp) ERROR STOP 'Metric/window basis covariance'
         rc(:, :, 1) = rc(:, :, 1) + 2.7_dp*sc
         rc(:, :, 2) = rc(:, :, 2) - 1.2_dp*sc
         CALL finite_chern_marker(cc(:, :n/2), sc, rc, wc, 9.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value - changed) > 1.e-11_dp) ERROR STOP 'Marker origin covariance'
         cc(:, 1) = (c(:, 1) + iu*c(:, 2))/SQRT(2.0_dp)
         cc(:, 2) = (iu*c(:, 1) + c(:, 2))/SQRT(2.0_dp)
         cc(:, 3:) = c(:, 3:)
         CALL finite_chern_marker(cc(:, :n/2), s, r, w, 9.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value - changed) > 1.e-11_dp) ERROR STOP 'Occupied subspace gauge'
         CALL finite_chern_marker(CONJG(c(:, :n/2)), s, r, w, 9.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value + changed) > 1.e-11_dp) ERROR STOP 'Marker time reversal'
         CALL finite_chern_marker(c(:, :n/2), s, r(:, :, [2, 1]), w, 9.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value + changed) > 1.e-11_dp) ERROR STOP 'Marker plane orientation'
         CALL finite_chern_marker(c(:, :n/2), s, 2*r, w, 36.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value - changed) > 1.e-11_dp) ERROR STOP 'Marker physical length units'
         CALL finite_chern_marker(c(:, :n/2), s, r, w, 18.0_dp, changed, status)
         IF (status /= 0 .OR. ABS(value - 2*changed) > 1.e-11_dp) ERROR STOP 'Marker area normalization'
         PRINT '(A,F5.1,4ES18.9)', 'Real-space Chern mass / marker / Wilson / full trace / AO residual:', &
            masses(phase), value, winding, full_trace, residual
      END DO

      CALL finite_chern_marker(c(:, :0), s, r, w, 9.0_dp, changed, status)
      IF (status /= 0 .OR. changed /= 0.0_dp) ERROR STOP 'Empty occupied marker'
      CALL finite_chern_marker(c, s, r, w, 9.0_dp, changed, status)
      IF (status /= 0 .OR. changed /= 0.0_dp) ERROR STOP 'Complete occupied marker'
      CALL finite_chern_marker(1.1_dp*c(:, :n/2), s, r, w, 9.0_dp, changed, status)
      IF (status /= -3) ERROR STOP 'Unnormalized occupied marker accepted'
      CALL finite_chern_marker(c(:, :n/2), s, r, -w, 9.0_dp, changed, status)
      IF (status /= -4) ERROR STOP 'Negative physical window accepted'
      CALL finite_chern_marker(c(:, :n/2), s, r, 2*w, 9.0_dp, changed, status)
      IF (status /= -4) ERROR STOP 'Window greater than identity accepted'
      CALL finite_chern_marker(c(:, :n/2), s, r, w, 0.0_dp, changed, status)
      IF (status /= -1) ERROR STOP 'Zero window area accepted'
      CALL finite_chern_marker(c(:, :n/2), s, r, w, ieee_value(0.0_dp, ieee_quiet_nan), changed, status)
      IF (status /= -1) ERROR STOP 'NaN window area accepted'
      s(1, 1) = -1.0_dp
      CALL finite_chern_marker(c(:, :n/2), s, r, w, 9.0_dp, changed, status)
      IF (status /= -2) ERROR STOP 'Indefinite marker metric accepted'
      PRINT *, 'Finite-window Chern, native Wilson, metric, gauge, orientation and boundary tests passed.'
   END SUBROUTINE check_real_space_chern

END MODULE kubo_test_cases

PROGRAM kubo_unittest
   USE kubo_test_cases, ONLY: check_projected_current, check_bloch_current, check_hall_response, check_real_space_chern
   IMPLICIT NONE
   CALL check_projected_current()
   CALL check_bloch_current()
   CALL check_hall_response()
   CALL check_real_space_chern()
END PROGRAM kubo_unittest

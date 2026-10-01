!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Independent chiral singular values, strong/weak/trivial Wilson-Dirac models and AO covariance.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM spectral_localizer_3d_unittest
   USE ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   USE spectral_localizer, ONLY: localizer_3d_evaluate, localizer_3d_matrix_evaluate
   IMPLICIT NONE
   COMPLEX(dp), PARAMETER :: ii = CMPLX(0.0_dp, 1.0_dp, dp)
   INTEGER, PARAMETER :: side = 3, sites = side**3, n = 4*sites
   COMPLEX(dp) :: h(n, n), s(n, n), r(n, n, 3), u(n, n), hc(n, n), sc(n, n), rc(n, n, 3)
   COMPLEX(dp) :: pauli(2, 2, 3), a(2*n, 2*n), dummy(1, 1), work(10*n)
   REAL(dp) :: singular(2*n), rwork(10*n), gap, reference, ratio, mass, shift(3), defect, lower, upper, midpoint
   INTEGER :: i, j, axis, model, index, expected, status, sign_product, kx, ky, kz
#if defined(__HAS_IEEE_EXCEPTIONS)
   LOGICAL :: halt(SIZE(IEEE_ALL)), current_halt(SIZE(IEEE_ALL))
   CALL ieee_get_halting_mode(IEEE_ALL, halt)
#endif

   CALL determinant_controls()
   pauli = 0.0_dp
   pauli(1, 2, 1) = 1.0_dp
   pauli(2, 1, 1) = 1.0_dp
   pauli(1, 2, 2) = -ii
   pauli(2, 1, 2) = ii
   pauli(1, 1, 3) = 1.0_dp
   pauli(2, 2, 3) = -1.0_dp
   DO model = 1, 3
      mass = REAL(2*(model - 1), dp)
      CALL dirac_model(mass)
      ! Bulk inversion eigenvalues at the eight TRIM: strong, not weak, parity.
      sign_product = 1
      DO kx = -1, 1, 2
         DO ky = -1, 1, 2
            DO kz = -1, 1, 2
               sign_product = sign_product*INT(SIGN(1.0_dp, mass + kx + ky + kz))
            END DO
         END DO
      END DO
      expected = (1 - sign_product)/2
      CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
      IF (status /= 0 .OR. index /= expected) ERROR STOP '3D strong/trivial parity mismatch'
      reference = gap
      a = 0.0_dp
      a(:n, :n) = h
      a(n + 1:, n + 1:) = h
      DO axis = 1, 3
         DO j = 1, 2
            DO i = 1, 2
               a((i - 1)*n + 1:i*n, (j - 1)*n + 1:j*n) = &
                  a((i - 1)*n + 1:i*n, (j - 1)*n + 1:j*n) - 0.5_dp*ii*pauli(i, j, axis)*r(:, :, axis)
            END DO
         END DO
      END DO
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zgesvd('N', 'N', 2*n, 2*n, a, 2*n, singular, dummy, 1, dummy, 1, work, SIZE(work), rwork, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0 .OR. ABS(reference - MINVAL(singular)) > 2.e-12_dp) ERROR STOP '3D singular gap mismatch'
      WRITE (*, '(A,F4.1,I3,2ES18.9)') '3D Dirac mass, index, gap, SVD error: ', &
         mass, expected, reference, ABS(reference - MINVAL(singular))
   END DO

   CALL dirac_model(2.0_dp)
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   reference = gap
   u = 0.0_dp
   DO i = 1, n/2
      u(i, i) = CMPLX(0.8_dp + 0.004_dp*i, 0.07_dp*SIN(REAL(i, dp)), dp)
      u(i + n/2, i + n/2) = CONJG(u(i, i))
      u(i, i + n/2) = CMPLX(0.08_dp, 0.03_dp, dp)
      u(i + n/2, i) = -CONJG(u(i, i + n/2))
      IF (i < n/2) THEN
         u(i, i + 1) = CMPLX(0.02_dp, 0.01_dp, dp)
         u(i + n/2, i + n/2 + 1) = CONJG(u(i, i + 1))
      END IF
   END DO
   hc = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(h, u))
   sc = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(s, u))
   DO axis = 1, 3
      rc(:, :, axis) = MATMUL(CONJG(TRANSPOSE(u)), MATMUL(r(:, :, axis), u))
   END DO
   CALL evaluate(hc, sc, rc, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   defect = ABS(gap - reference)
   IF (status /= 0 .OR. index /= 1 .OR. defect > 3.e-12_dp) ERROR STOP '3D nonorthogonal spin-mixing basis'
   WRITE (*, '(A,ES18.9)') '3D nonorthogonal spin-mixing gap residual: ', defect
   shift = [1.2_dp, -0.8_dp, 0.4_dp]
   DO axis = 1, 3
      rc(:, :, axis) = r(:, :, axis) + shift(axis)*s
   END DO
   CALL evaluate(h + 3.0_dp*s, s, rc, shift, 3.0_dp)
   IF (status /= 0 .OR. index /= 1 .OR. ABS(gap - reference) > 3.e-12_dp) ERROR STOP '3D origin shifts'
   rc(:, :, 1) = 0.6_dp*r(:, :, 1) + 0.8_dp*r(:, :, 2)
   rc(:, :, 2) = -0.8_dp*r(:, :, 1) + 0.6_dp*r(:, :, 2)
   rc(:, :, 3) = -r(:, :, 3)
   CALL evaluate(h, s, rc, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   IF (status /= 0 .OR. index /= 1 .OR. ABS(gap - reference) > 3.e-12_dp) ERROR STOP '3D reflection/rotation'
   ! A too-small finite box/scale is not automatically classified by the bulk parity.
   CALL localizer_3d_evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 0.3_dp, &
                              1.e-10_dp, 1.e-10_dp, 1.e-10_dp, index, gap, ratio, status)
   IF (status /= 0 .OR. index /= 0) ERROR STOP '3D finite-size counterexample lost'
   lower = 0.3_dp
   upper = 0.5_dp
   DO i = 1, 32
      midpoint = 0.5_dp*(lower + upper)
      CALL localizer_3d_evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, midpoint, &
                                 1.e-10_dp, 1.e-10_dp, 0.0_dp, index, gap, ratio, status)
      IF (status == 1) EXIT
      IF (status /= 0) ERROR STOP '3D transition factorization failed'
      IF (index == 0) THEN
         lower = midpoint
      ELSE
         upper = midpoint
      END IF
   END DO
   CALL localizer_3d_evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, midpoint, &
                              1.e-10_dp, 1.e-10_dp, 1.e-8_dp, index, gap, ratio, status)
   IF (status /= 1) ERROR STOP '3D closing gap assigned an index'
   WRITE (*, '(A,2ES18.9)') '3D localizer transition scale and unresolved gap: ', midpoint, gap
   CALL localizer_3d_evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 0.5_dp, &
                              1.e-10_dp, -1.0_dp, 1.e-10_dp, index, gap, ratio, status)
   IF (status /= -1) ERROR STOP '3D invalid metric tolerance accepted'
   CALL localizer_3d_evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 0.5_dp, &
                              1.e-10_dp, 1.e-10_dp, ieee_value(0.0_dp, ieee_quiet_nan), index, gap, ratio, status)
   IF (status /= -1) ERROR STOP '3D nonfinite gap tolerance accepted'
   h(1, 1) = h(1, 1) + 0.1_dp
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   IF (status /= -3) ERROR STOP '3D broken time reversal accepted'
   h = 0.0_dp
   r = 0.0_dp
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   IF (status /= 1 .OR. index /= 0) ERROR STOP '3D zero gap accepted'
   s(1, 1) = -1.0_dp
   s(n/2 + 1, n/2 + 1) = -1.0_dp
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   IF (status /= -2) ERROR STOP '3D indefinite metric accepted'
   h(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
   CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
   IF (status /= -1) ERROR STOP '3D nonfinite operator accepted'
#if defined(__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, current_halt)
   IF (ANY(halt .NEQV. current_halt)) ERROR STOP '3D localizer changed IEEE exception modes'
#endif
   PRINT *, '3D AII determinant, generalized gap, basis, origin, TR and finite-size tests passed.'
CONTAINS
! **************************************************************************************************
!> \brief Chiral-block determinant controls independent of the Wilson-Dirac model.
!> \note Wilkinson's matrix has exponential LU pivot growth despite a large singular gap.
! **************************************************************************************************
   SUBROUTINE determinant_controls()
      INTEGER, PARAMETER                                 :: m = 64, quarter = m/4

      COMPLEX(dp)                                        :: block(m, m), metric(m/2, m/2), &
                                                            unused(1, 1), workspace(5*m)
      INTEGER                                            :: col, col_sign, info, orientation, &
                                                            paired_col, paired_row, RESULT, row, &
                                                            row_sign
      REAL(dp)                                           :: matrix(m, m), measured_gap, &
                                                            metric_condition, real_work(5*m), &
                                                            svd_gap, values(m)

      metric = 0.0_dp
      DO row = 1, m/2
         metric(row, row) = 1.0_dp
      END DO
      matrix = 0.0_dp
      DO row = 1, m
         matrix(row, row) = 1.0_dp
         matrix(row, :row - 1) = -1.0_dp
         matrix(row, m) = 1.0_dp + 0.25_dp*MODULO(row, 3)
      END DO
      DO orientation = 0, 1
         ! Positive last-column entries give a positive determinant. One row reversal changes its sign.
         IF (orientation == 1) matrix(1, :) = -matrix(1, :)
         ! Exact dyadic entries in Q W Q^dagger retain the pivot-growth fixture after conversion.
         DO col = 1, m
            paired_col = (3 - (col - 1)/quarter)*quarter + MODULO(col - 1, quarter) + 1
            col_sign = MERGE(1, -1, col <= quarter .OR. col > 3*quarter)
            DO row = 1, m
               paired_row = (3 - (row - 1)/quarter)*quarter + MODULO(row - 1, quarter) + 1
               row_sign = MERGE(1, -1, row <= quarter .OR. row > 3*quarter)
               block(row, col) = 0.5_dp*(matrix(row, col) + row_sign*col_sign*matrix(paired_row, paired_col) + &
                                         ii*(col_sign*matrix(row, paired_col) - row_sign*matrix(paired_row, col)))
            END DO
         END DO
         CALL localizer_3d_matrix_evaluate(block, metric, 1.e-10_dp, 1.e-10_dp, 1.e-10_dp, &
                                           RESULT, measured_gap, metric_condition, info)
         IF (info /= 0 .OR. RESULT /= orientation) ERROR STOP '3D high-growth determinant sign'
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
         CALL zgesvd('N', 'N', m, m, block, m, values, unused, 1, unused, 1, workspace, &
                     SIZE(workspace), real_work, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
         svd_gap = MINVAL(values)
         IF (info /= 0 .OR. ABS(measured_gap - svd_gap) > 1.e-12_dp) ERROR STOP '3D high-growth singular gap'
         WRITE (*, '(A,I2,2ES18.9)') '3D high-growth determinant, gap, SVD error: ', &
            RESULT, measured_gap, ABS(measured_gap - svd_gap)
      END DO
      ! Independent atomic limit for both signs of a spin-degenerate energy.
      DO col = -1, 1, 2
         block = 0.0_dp
         DO row = 1, m
            block(row, row) = REAL(col, dp)
         END DO
         CALL localizer_3d_matrix_evaluate(block, metric, 1.e-10_dp, 1.e-10_dp, 1.e-10_dp, &
                                           RESULT, measured_gap, metric_condition, info)
         IF (info /= 0 .OR. RESULT /= 0 .OR. ABS(measured_gap - 1.0_dp) > 1.e-12_dp) THEN
            ERROR STOP '3D atomic reference orientation'
         END IF
      END DO
   END SUBROUTINE determinant_controls

! **************************************************************************************************
!> \brief Fixed tolerances and localizer scale for the small Wilson-Dirac control.
!> \param hh Physical Hamiltonian
!> \param ss Physical overlap
!> \param rr Three Cartesian moments
!> \param point Query position
!> \param energy Query energy
! **************************************************************************************************
   SUBROUTINE evaluate(hh, ss, rr, point, energy)
      COMPLEX(dp), INTENT(IN)                            :: hh(:, :), ss(:, :), rr(:, :, :)
      REAL(dp), INTENT(IN)                               :: point(3), energy

      CALL localizer_3d_evaluate(hh, ss, rr, point, energy, 0.5_dp, 1.e-10_dp, 1.e-10_dp, 1.e-10_dp, &
                                 index, gap, ratio, status)
   END SUBROUTINE evaluate
! **************************************************************************************************
!> \brief Open 3x3x3 Wilson-Dirac cube in fixed physical up/down ordering.
!> \param m Onsite mass; periodic Bloch mass is m+cos(kx)+cos(ky)+cos(kz)
! **************************************************************************************************
   SUBROUTINE dirac_model(m)
      REAL(dp), INTENT(IN)                               :: m

      INTEGER                                            :: col, direction, neighbor, orbital, &
                                                            other, pos(3), row, site, spin, x, y, z

      h = 0.0_dp
      s = 0.0_dp
      r = 0.0_dp
      DO z = 0, side - 1
         DO y = 0, side - 1
            DO x = 0, side - 1
               pos = [x, y, z]
               site = 1 + x + side*y + side**2*z
               DO spin = 1, 2
                  DO orbital = 1, 2
                     row = site + sites*(orbital - 1) + 2*sites*(spin - 1)
                     s(row, row) = 1.0_dp
                     h(row, row) = m*REAL(3 - 2*orbital, dp)
                     r(row, row, :) = REAL(pos - 1, dp)
                     DO direction = 1, 3
                        IF (pos(direction) == side - 1) CYCLE
                        neighbor = site + side**(direction - 1)
                        col = neighbor + sites*(orbital - 1) + 2*sites*(spin - 1)
                        h(row, col) = 0.5_dp*REAL(3 - 2*orbital, dp)
                        h(col, row) = CONJG(h(row, col))
                        DO other = 1, 2
                           col = neighbor + sites*(2 - orbital) + 2*sites*(other - 1)
                           h(row, col) = -0.5_dp*ii*pauli(spin, other, direction)
                           h(col, row) = CONJG(h(row, col))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE dirac_model
END PROGRAM spectral_localizer_3d_unittest

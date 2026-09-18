!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM topology_wilson_unittest
   USE iso_fortran_env,                 ONLY: output_unit
   USE kinds,                           ONLY: dp
   USE topology_wilson,                 ONLY: chern_from_wcc,&
                                              surface_resolved,&
                                              wcc_distance,&
                                              wilson_spectrum,&
                                              wilson_step,&
                                              z2_from_wcc

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp, tolerance = 1.e-10_dp

   INTEGER :: io_unit

   io_unit = output_unit
   CALL test_model(-1.0_dp, 1)
   CALL test_model(-3.0_dp, 0)
   CALL test_model(1.0_dp, 1)
   CALL test_model(3.0_dp, 0)
   CALL test_checks()
   CALL test_chern_model(-1.0_dp, 1)
   CALL test_chern_model(1.0_dp, -1)
   CALL test_chern_model(-3.0_dp, 0)
   CALL test_chern_checks()
   WRITE (io_unit, *) "Wilson-loop, Z2 and Chern unit tests passed."

CONTAINS

! **************************************************************************************************
!> \brief Known first Chern numbers of the occupied two-band lattice Dirac model.
!> \param mass Dirac mass in h = (sin kx, sin ky, m+cos kx+cos ky).sigma
!> \param expected Chern number in the native/Z2Pack Wilson-winding convention
! **************************************************************************************************
   SUBROUTINE test_chern_model(mass, expected)
      REAL(KIND=dp), INTENT(IN)                          :: mass
      INTEGER, INTENT(IN)                                :: expected

      INTEGER, PARAMETER                                 :: nline = 65, npoint = 41

      COMPLEX(KIND=dp)                                   :: h(2, 2), overlap(1, 1), PRODUCT(1, 1), &
                                                            states(2, npoint), work(8)
      INTEGER                                            :: i, invariant, j, line, status
      REAL(KIND=dp)                                      :: berry, d, doubled(2, nline), &
                                                            eigenvalues(2), kx, ky, minimum_sv, &
                                                            rwork(4), wcc(1, nline), winding

      DO line = 1, nline
         ky = 2.0_dp*pi*REAL(line - 1, dp)/REAL(nline - 1, dp)
         DO i = 1, npoint
            kx = 2.0_dp*pi*REAL(i - 1, dp)/REAL(npoint, dp)
            d = mass + COS(kx) + COS(ky)
            h(1, 1) = d
            h(2, 2) = -d
            h(1, 2) = CMPLX(SIN(kx), -SIN(ky), dp)
            h(2, 1) = CONJG(h(1, 2))
            CALL zheev('V', 'U', 2, h, 2, eigenvalues, work, SIZE(work), rwork, status)
            IF (status /= 0) ERROR STOP "Chern model diagonalization failed"
            states(:, i) = h(:, 1)
         END DO
         PRODUCT(:, :) = 1.0_dp
         DO i = 1, npoint
            j = MOD(i, npoint) + 1
            overlap(1, 1) = DOT_PRODUCT(states(:, i), states(:, j))
            CALL wilson_step(product, overlap, minimum_sv, status, tolerance)
            IF (status /= 0) ERROR STOP "Chern model overlap is singular"
         END DO
         CALL wilson_spectrum(product, wcc(:, line), berry, status)
         IF (status /= 0) ERROR STOP "Chern model spectrum failed"
      END DO
      CALL chern_from_wcc(wcc, invariant, winding, status, tolerance)
      IF (status /= 0 .OR. invariant /= expected) ERROR STOP "Incorrect model Chern number"
      IF (ABS(winding - REAL(expected, dp)) > tolerance) ERROR STOP "Incorrect Chern winding"
      CALL chern_from_wcc(wcc(:, nline:1:-1), invariant, winding, status, tolerance)
      IF (status /= 0 .OR. invariant /= -expected) ERROR STOP "Chern orientation reversal failed"
      doubled(1, :) = wcc(1, :)
      doubled(2, :) = wcc(1, :)
      CALL chern_from_wcc(doubled, invariant, winding, status, tolerance)
      IF (status /= 0 .OR. invariant /= 2*expected) ERROR STOP "Direct-sum Chern additivity failed"
   END SUBROUTINE test_chern_model

! **************************************************************************************************
!> \brief Reject unclosed surfaces and unresolved determinant-phase steps.
! **************************************************************************************************
   SUBROUTINE test_chern_checks()
      INTEGER                                            :: invariant, status
      REAL(KIND=dp)                                      :: surface(1, 5), winding

      surface(1, :) = [0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
      CALL chern_from_wcc(surface, invariant, winding, status, tolerance)
      IF (status /= -2) ERROR STOP "Unclosed Chern surface was accepted"
      surface(1, :) = [0.0_dp, 0.3_dp, 0.6_dp, 0.9_dp, 1.0_dp]
      CALL chern_from_wcc(surface, invariant, winding, status, tolerance)
      IF (status /= -3) ERROR STOP "Unresolved Chern phase step was accepted"
      surface(:, :) = 0.1_dp
      CALL chern_from_wcc(surface, invariant, winding, status, tolerance)
      IF (status /= 0 .OR. invariant /= 0) ERROR STOP "Constant Chern surface was rejected"
   END SUBROUTINE test_chern_checks

! **************************************************************************************************
!> \brief Known trivial/nontrivial lattice BHZ models, non-Abelian gauge changes and loop reversal.
!> \param mass Dirac mass; |mass| < 2 (excluding zero) is topological, |mass| > 2 is trivial
!> \param expected expected Z2 invariant
! **************************************************************************************************
   SUBROUTINE test_model(mass, expected)
      REAL(KIND=dp), INTENT(IN)                          :: mass
      INTEGER, INTENT(IN)                                :: expected

      INTEGER, PARAMETER                                 :: nline = 31, npoint = 41

      COMPLEX(KIND=dp) :: gauged(2, 2), gauges(2, 2, npoint), links(2, 2, npoint), overlap(2, 2), &
         PRODUCT(2, 2), reversed(2, 2), states(4, 2, npoint)
      INTEGER                                            :: i, invariant, j, line, status
      REAL(KIND=dp)                                      :: berry, gauge_wcc(2, nline), kx, ky, &
                                                            minimum_sv, reverse_wcc(2, nline), &
                                                            wcc(2, nline)

      DO line = 1, nline
         ky = pi*REAL(line - 1, dp)/REAL(nline - 1, dp)
         DO i = 1, npoint
            kx = 2.0_dp*pi*REAL(i - 1, dp)/REAL(npoint, dp)
            CALL occupied_states(kx, ky, mass, states(:, :, i))
            CALL gauge_matrix(kx, ky, gauges(:, :, i))
         END DO
         PRODUCT(:, :) = 0.0_dp
         PRODUCT(1, 1) = 1.0_dp
         PRODUCT(2, 2) = 1.0_dp
         gauged(:, :) = product
         reversed(:, :) = product
         DO i = 1, npoint
            j = MOD(i, npoint) + 1
            links(:, :, i) = MATMUL(CONJG(TRANSPOSE(states(:, :, i))), states(:, :, j))
            CALL wilson_step(product, links(:, :, i), minimum_sv, status, tolerance)
            IF (status /= 0) ERROR STOP "BHZ overlap is singular"
            overlap(:, :) = MATMUL(CONJG(TRANSPOSE(gauges(:, :, i))), &
                                   MATMUL(links(:, :, i), gauges(:, :, j)))
            CALL wilson_step(gauged, overlap, minimum_sv, status, tolerance)
            IF (status /= 0) ERROR STOP "Gauge-transformed BHZ overlap is singular"
         END DO
         DO i = npoint, 1, -1
            overlap(:, :) = CONJG(TRANSPOSE(links(:, :, i)))
            CALL wilson_step(reversed, overlap, minimum_sv, status, tolerance)
            IF (status /= 0) ERROR STOP "Reversed BHZ overlap is singular"
         END DO
         CALL wilson_spectrum(product, wcc(:, line), berry, status)
         IF (status /= 0) ERROR STOP "BHZ spectrum failed"
         CALL wilson_spectrum(gauged, gauge_wcc(:, line), berry, status)
         IF (status /= 0) ERROR STOP "Gauge-transformed spectrum failed"
         CALL wilson_spectrum(reversed, reverse_wcc(:, line), berry, status)
         IF (status /= 0) ERROR STOP "Reversed spectrum failed"
         IF (wcc_distance(wcc(:, line), gauge_wcc(:, line)) > tolerance) THEN
            ERROR STOP "Wilson spectrum depends on occupied-space gauge"
         END IF
         IF (wcc_distance(wcc(:, line), -reverse_wcc(:, line)) > tolerance) THEN
            ERROR STOP "Loop reversal did not conjugate Wilson eigenvalues"
         END IF
      END DO
      CALL z2_from_wcc(wcc, invariant, status, tolerance)
      IF (status /= 0 .OR. invariant /= expected) ERROR STOP "Incorrect BHZ Z2 invariant"
      CALL z2_from_wcc(gauge_wcc, invariant, status, tolerance)
      IF (status /= 0 .OR. invariant /= expected) ERROR STOP "Gauge-dependent Z2 invariant"
      CALL z2_from_wcc(reverse_wcc, invariant, status, tolerance)
      IF (status /= 0 .OR. invariant /= expected) ERROR STOP "Reversal-dependent Z2 invariant"
   END SUBROUTINE test_model

! **************************************************************************************************
!> \brief Occupied eigenvectors of diag(h(k), h*(-k)), with h = (sin kx, sin ky, m+cos kx+cos ky).sigma.
!> \param kx first reduced angle
!> \param ky second reduced angle
!> \param mass Dirac mass
!> \param states two occupied eigenvectors
! **************************************************************************************************
   SUBROUTINE occupied_states(kx, ky, mass, states)
      REAL(KIND=dp), INTENT(IN)                          :: kx, ky, mass
      COMPLEX(KIND=dp), INTENT(OUT)                      :: states(4, 2)

      COMPLEX(KIND=dp)                                   :: h(4, 4), work(16)
      INTEGER                                            :: status
      REAL(KIND=dp)                                      :: d, eigenvalues(4), rwork(10)

      d = mass + COS(kx) + COS(ky)
      h(:, :) = 0.0_dp
      h(1, 1) = d
      h(2, 2) = -d
      h(1, 2) = CMPLX(SIN(kx), -SIN(ky), dp)
      h(2, 1) = CONJG(h(1, 2))
      h(3, 3) = d
      h(4, 4) = -d
      h(3, 4) = CMPLX(-SIN(kx), -SIN(ky), dp)
      h(4, 3) = CONJG(h(3, 4))
      CALL zheev('V', 'U', 4, h, 4, eigenvalues, work, SIZE(work), rwork, status)
      IF (status /= 0) ERROR STOP "BHZ diagonalization failed"
      states(:, :) = h(:, 1:2)
   END SUBROUTINE occupied_states

! **************************************************************************************************
!> \brief Deterministic k-dependent U(2) rotations, mixing the two occupied states.
!> \param kx first reduced angle
!> \param ky second reduced angle
!> \param matrix gauge matrix
! **************************************************************************************************
   SUBROUTINE gauge_matrix(kx, ky, matrix)
      REAL(KIND=dp), INTENT(IN)                          :: kx, ky
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrix(2, 2)

      COMPLEX(KIND=dp)                                   :: determinant_phase, phase
      REAL(KIND=dp)                                      :: angle

      angle = 0.37_dp + 0.61_dp*SIN(3.0_dp*kx + ky)
      phase = EXP(CMPLX(0.0_dp, 0.73_dp*COS(kx - 2.0_dp*ky), dp))
      determinant_phase = EXP(CMPLX(0.0_dp, 0.29_dp*SIN(kx + ky), dp))
      matrix(1, 1) = COS(angle)
      matrix(1, 2) = SIN(angle)*phase
      matrix(2, 1) = -SIN(angle)*CONJG(phase)
      matrix(2, 2) = COS(angle)
      matrix(:, 1) = matrix(:, 1)*determinant_phase
   END SUBROUTINE gauge_matrix

! **************************************************************************************************
!> \brief Analytic polar/Berry phases and singular-link, Kramers-pair and surface-resolution checks.
! **************************************************************************************************
   SUBROUTINE test_checks()
      COMPLEX(KIND=dp)                                   :: overlap(2, 2), PRODUCT(2, 2)
      INTEGER                                            :: invariant, status
      REAL(KIND=dp)                                      :: berry, minimum_sv, surface(2, 3), wcc(2)

      PRODUCT(:, :) = 0.0_dp
      PRODUCT(1, 1) = 1.0_dp
      PRODUCT(2, 2) = 1.0_dp
      overlap(:, :) = 0.0_dp
      CALL wilson_step(product, overlap, minimum_sv, status, tolerance)
      IF (status /= -2) ERROR STOP "Singular overlap was accepted"

      overlap(1, 1) = 0.8_dp*EXP(CMPLX(0.0_dp, 0.2_dp*pi, dp))
      overlap(2, 2) = 0.6_dp*EXP(CMPLX(0.0_dp, 0.6_dp*pi, dp))
      CALL wilson_step(product, overlap, minimum_sv, status, tolerance)
      IF (status /= 0 .OR. ABS(minimum_sv - 0.6_dp) > tolerance) THEN
         ERROR STOP "Polar factor singular value is incorrect"
      END IF
      CALL wilson_spectrum(product, wcc, berry, status)
      IF (status /= 0) ERROR STOP "Analytic spectrum failed"
      IF (wcc_distance(wcc, [0.1_dp, 0.3_dp]) > tolerance) ERROR STOP "Incorrect analytic WCC"
      IF (ABS(berry - 0.8_dp*pi) > tolerance) ERROR STOP "Incorrect analytic Berry phase"

      surface(:, 1) = wcc
      surface(:, 2) = wcc
      surface(:, 3) = wcc
      CALL z2_from_wcc(surface, invariant, status, tolerance)
      IF (status /= -2 .OR. invariant /= -1) ERROR STOP "Missing Kramers pairs were accepted"
      surface(:, :) = 0.1_dp
      surface(:, 2) = 0.6_dp
      IF (surface_resolved(surface)) ERROR STOP "Unresolved surface was accepted"
      surface(:, :) = 0.1_dp
      IF (.NOT. surface_resolved(surface)) ERROR STOP "Constant surface was rejected"
      CALL z2_from_wcc(surface, invariant, status, tolerance)
      IF (status /= 0 .OR. invariant /= 0) ERROR STOP "Incorrect constant-surface parity"
      IF (wcc_distance([0.99_dp, 0.2_dp], [0.2_dp, -0.01_dp]) > tolerance) THEN
         ERROR STOP "Circular WCC matching failed"
      END IF
   END SUBROUTINE test_checks

END PROGRAM topology_wilson_unittest

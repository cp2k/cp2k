!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Antisymmetric DC charge response with a scalar relaxation rate.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE kubo_hall
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE kubo_bloch,                      ONLY: fermi_slope

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hall_response

   REAL(dp), PARAMETER :: antihermiticity_tolerance = 1.e-9_dp
CONTAINS
! **************************************************************************************************
!> \brief Hall tensor, with first index the current and second index the electric field.
!> \param energies Complete eigenvalue list for one spin channel and k point
!> \param current Anti-Hermitian -i*hbar*v matrix triplet, in exactly the same eigenframe
!> \param mu Chemical potential
!> \param kt Positive temperature in energy units
!> \param gamma Positive relaxation rate in energy units
!> \param hall Antisymmetric tensor before e^2/hbar, volume, spin and k weights
!> \param status Zero on success; negative for invalid input or nonfinite arithmetic
!> \note The trace uses (gamma+i[H,.])^-1. Pairing n,m gives
!>       2*f[en,em]*(en-em)*Im[J_a(m,n)*J_b(n,m)]/(gamma^2+(en-em)^2).
!>       Equal-energy blocks contribute zero without a band-gap division.
! **************************************************************************************************
   SUBROUTINE hall_response(energies, current, mu, kt, gamma, hall, status)
      REAL(dp), INTENT(IN)                               :: energies(:)
      COMPLEX(dp), INTENT(IN)                            :: current(:, :, :)
      REAL(dp), INTENT(IN)                               :: mu, kt, gamma
      REAL(dp), INTENT(OUT)                              :: hall(3, 3)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: a, b, m, n, nstate
      REAL(dp)                                           :: delta, factor, scale, total(3, 3)

      hall = 0.0_dp
      status = -1
      nstate = SIZE(energies)
      IF (nstate < 1 .OR. ANY(SHAPE(current) /= [nstate, nstate, 3])) RETURN
      IF (.NOT. ALL(ieee_is_finite([mu, kt, gamma]))) RETURN
      IF (kt <= 0.0_dp .OR. gamma <= 0.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(energies))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(current, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(current)))) RETURN
      DO a = 1, 3
         scale = MAX(1.0_dp, MAXVAL(ABS(current(:, :, a))))
         IF (MAXVAL(ABS(current(:, :, a) + CONJG(TRANSPOSE(current(:, :, a))))) > antihermiticity_tolerance*scale) RETURN
      END DO
      total = 0.0_dp
      DO n = 1, nstate
         DO m = n + 1, nstate
            delta = energies(n) - energies(m)
            IF (.NOT. ieee_is_finite(delta)) RETURN
            IF (delta == 0.0_dp) CYCLE
            scale = MAX(gamma, ABS(delta))
            factor = 2.0_dp*fermi_slope(energies(n), energies(m), mu, kt)* &
                     (delta/scale)/scale/((gamma/scale)**2 + (delta/scale)**2)
            IF (.NOT. ieee_is_finite(factor)) RETURN
            DO b = 2, 3
               DO a = 1, b - 1
                  total(a, b) = total(a, b) + factor*AIMAG(current(m, n, a)*current(n, m, b))
               END DO
            END DO
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(total))) RETURN
      hall = total - TRANSPOSE(total)
      status = 0
   END SUBROUTINE hall_response
END MODULE kubo_hall

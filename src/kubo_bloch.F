!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Metric-aware Bloch velocities without eigenvector derivatives or band-gap divisions.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE kubo_bloch
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE kubo_projected,                  ONLY: projected_eigensystem

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: bloch_current, bloch_current_from_frame, fermi_slope
CONTAINS
! **************************************************************************************************
!> \brief Evaluate v_nm = C^H(dH - E_n*dS)C + i*(E_n-E_m)*C^H*D*C.
!> \param h Full covariant Hamiltonian, including SOC when applicable
!> \param s Overlap in the same Bloch convention
!> \param dh Cartesian k derivatives of h, in Ha*bohr
!> \param ds Cartesian k derivatives of s, in bohr
!> \param d Cell-gauge image first moments, satisfying D-D^H=-i*dS, in bohr
!> \param energies Generalized spectrum
!> \param current Anti-Hermitian current -i*v in the eigenbasis, including diagonal velocities
!> \param status 0 success; -1 invalid data; -2 metric failure; -3 connection failure; positive LAPACK failure
!> \param coefficients Optional complete eigenframe for symmetry reconstruction
!> \note For a general k-dependent basis D transforms as B^H*D*B-i*(dB)^H*S*B.
!>       The connection terms are essential even when the spectrum is accurately reproduced.
! **************************************************************************************************
   SUBROUTINE bloch_current(h, s, dh, ds, d, energies, current, status, coefficients)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :), dh(:, :, :), &
                                                            ds(:, :, :), d(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energies(:)
      COMPLEX(dp), INTENT(OUT)                           :: current(:, :, :)
      INTEGER, INTENT(OUT)                               :: status
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: coefficients(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: c(:, :)
      INTEGER                                            :: n
      REAL(dp)                                           :: residuals(2)

      n = SIZE(h, 1)
      status = -1
      energies = 0.0_dp
      current = 0.0_dp
      IF (PRESENT(coefficients)) THEN
         coefficients = 0.0_dp
         IF (ANY(SHAPE(coefficients) /= [n, n])) RETURN
      END IF
      ALLOCATE (c(n, n))
      CALL projected_eigensystem(h, s, energies, c, status)
      IF (status /= 0) RETURN
      CALL bloch_current_from_frame(h, s, dh, ds, d, energies, c, current, 1.e-8_dp, residuals, status)
      IF (status == 0 .AND. PRESENT(coefficients)) coefficients = c
   END SUBROUTINE bloch_current
! **************************************************************************************************
!> \brief Evaluate the current in a validated reconstructed eigenframe, without rediagonalization.
!> \param h Target Hamiltonian
!> \param s Target AO metric
!> \param dh Cartesian Hamiltonian derivatives
!> \param ds Cartesian overlap derivatives
!> \param d Cell-gauge image first moments
!> \param energies Source spectrum, including complete degenerate blocks
!> \param c Candidate target eigenframe
!> \param current Target current in exactly the supplied gauge
!> \param tolerance Maximum orthonormality and relative projected-eigenproblem errors
!> \param residuals Metric error and energy-scaled projected-Hamiltonian error
!> \param status 0 success; -1 invalid data; -3 connection failure; -4 invalid eigenframe
! **************************************************************************************************
   SUBROUTINE bloch_current_from_frame(h, s, dh, ds, d, energies, c, current, tolerance, residuals, status)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :), dh(:, :, :), &
                                                            ds(:, :, :), d(:, :, :)
      REAL(dp), INTENT(IN)                               :: energies(:)
      COMPLEX(dp), INTENT(IN)                            :: c(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: current(:, :, :)
      REAL(dp), INTENT(IN)                               :: tolerance
      REAL(dp), INTENT(OUT)                              :: residuals(2)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), PARAMETER                             :: iu = (0.0_dp, 1.0_dp)

      COMPLEX(dp), ALLOCATABLE                           :: dc(:, :), hc(:, :), sc(:, :), v(:, :)
      INTEGER                                            :: axis, i, j, n
      REAL(dp)                                           :: scale

      status = -1
      current = 0.0_dp
      residuals = HUGE(1.0_dp)
      n = SIZE(h, 1)
      IF (n < 1 .OR. ANY(SHAPE(h) /= [n, n]) .OR. ANY(SHAPE(s) /= [n, n])) RETURN
      IF (SIZE(energies) /= n .OR. ANY(SHAPE(current) /= [n, n, 3])) RETURN
      IF (ANY(SHAPE(dh) /= [n, n, 3]) .OR. ANY(SHAPE(ds) /= [n, n, 3])) RETURN
      IF (ANY(SHAPE(d) /= [n, n, 3])) RETURN
      IF (ANY(SHAPE(c) /= [n, n])) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(energies))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(c, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(c)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(h, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(h)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(s, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(s)))) RETURN
      IF (MAXVAL(ABS(h - CONJG(TRANSPOSE(h)))) > 1.e-9_dp*MAX(1.0_dp, MAXVAL(ABS(h)))) RETURN
      IF (MAXVAL(ABS(s - CONJG(TRANSPOSE(s)))) > 1.e-9_dp*MAX(1.0_dp, MAXVAL(ABS(s)))) RETURN
      IF (.NOT. finite(dh) .OR. .NOT. finite(ds) .OR. .NOT. finite(d)) RETURN
      DO axis = 1, 3
         scale = MAX(TINY(1.0_dp), MAXVAL(ABS(dh(:, :, axis))))
         IF (MAXVAL(ABS(dh(:, :, axis) - CONJG(TRANSPOSE(dh(:, :, axis))))) > 1.e-9_dp*scale) RETURN
         scale = MAX(TINY(1.0_dp), MAXVAL(ABS(ds(:, :, axis))))
         IF (MAXVAL(ABS(ds(:, :, axis) - CONJG(TRANSPOSE(ds(:, :, axis))))) > 1.e-9_dp*scale) RETURN
         status = -3
         scale = MAX(1.0_dp, MAXVAL(ABS(ds(:, :, axis))), MAXVAL(ABS(d(:, :, axis))))
         IF (MAXVAL(ABS(d(:, :, axis) - CONJG(TRANSPOSE(d(:, :, axis))) + iu*ds(:, :, axis))) &
             > 1.e-9_dp*scale) RETURN
         status = -1
      END DO
      ALLOCATE (hc(n, n), sc(n, n), dc(n, n), v(n, n))
      sc = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(s, c))
      hc = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(h, c))
      DO i = 1, n
         sc(i, i) = sc(i, i) - 1.0_dp
         hc(i, i) = hc(i, i) - energies(i)
      END DO
      residuals = [MAXVAL(ABS(sc)), MAXVAL(ABS(hc))/MAX(1.0_dp, MAXVAL(ABS(energies)))]
      status = -4
      IF (.NOT. ALL(ieee_is_finite(residuals)) .OR. MAXVAL(residuals) > tolerance) RETURN
      status = -1
      DO axis = 1, 3
         hc = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(dh(:, :, axis), c))
         sc = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(ds(:, :, axis), c))
         dc = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(d(:, :, axis), c))
         DO j = 1, n
            DO i = 1, n
               v(i, j) = hc(i, j) - energies(i)*sc(i, j) + iu*(energies(i) - energies(j))*dc(i, j)
            END DO
         END DO
         current(:, :, axis) = -iu*v
      END DO
      IF (.NOT. finite(current)) RETURN
      status = 0
   END SUBROUTINE bloch_current_from_frame
! **************************************************************************************************
!> \brief Finite check on all components, including connection matrices which are not Hermitian.
!> \param a Matrix triplet
!> \return Whether all entries are finite
! **************************************************************************************************
   LOGICAL FUNCTION finite(a)
      COMPLEX(dp), INTENT(IN)                            :: a(:, :, :)

      finite = ALL(ieee_is_finite(REAL(a, dp))) .AND. ALL(ieee_is_finite(AIMAG(a)))
   END FUNCTION finite
! **************************************************************************************************
!> \brief Stable Fermi divided difference including the intraband and degenerate limits.
!> \param ea First energy
!> \param eb Second energy
!> \param mu Chemical potential
!> \param kt Positive temperature in energy units
!> \return (f(ea)-f(eb))/(ea-eb) for unit occupation, nonpositive
! **************************************************************************************************
   PURE REAL(dp) FUNCTION fermi_slope(ea, eb, mu, kt) RESULT(slope)
      REAL(dp), INTENT(IN) :: ea, eb, mu, kt
      REAL(dp) :: delta, f

      delta = ea - eb
      IF (ABS(delta) < 1.e-4_dp*kt) THEN
         f = occupation((ea + eb)*0.5_dp)
         slope = -f*(1.0_dp - f)/kt*(1.0_dp + (delta/kt)**2*(1.0_dp - 6*f + 6*f*f)/24.0_dp)
      ELSE
         slope = (occupation(ea) - occupation(eb))/delta
      END IF
   CONTAINS
! **************************************************************************************************
!> \brief Stable logistic occupation.
!> \param e Energy
!> \return Unit-occupation Fermi function
! **************************************************************************************************
      PURE REAL(dp) FUNCTION occupation(e) RESULT(f)
      REAL(dp), INTENT(IN)                               :: e

      REAL(dp)                                           :: x

         x = EXP(-ABS((e - mu)/kt))
         f = x/(1 + x)
         IF (e < mu) f = 1.0_dp - f
      END FUNCTION occupation
   END FUNCTION fermi_slope
END MODULE kubo_bloch

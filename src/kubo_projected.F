!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite projected-Hamiltonian current in a nonorthogonal, possibly spinor AO basis.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE kubo_projected
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: projected_current, projected_eigensystem
CONTAINS
! **************************************************************************************************
!> \brief Compute C^dagger (H S^-1 X - X S^-1 H) C without forming S^-1.
!> \param h Complete scalar/spinor covariant Hamiltonian, not an occupied-state projection
!> \param s Positive definite overlap, with explicit spin blocks when applicable
!> \param positions Three analytic covariant Cartesian first moments
!> \param energies Eigenvalues of H C = S C E
!> \param current Anti-Hermitian commutators in the eigenbasis; velocity is i*current in atomic units
!> \param status 0 success; -1 invalid matrices; -2 ill-conditioned metric; positive LAPACK failure
!> \note This is the commutator of projected operators. It is not the projection of the
!>       continuum commutator in an incomplete basis, nor a periodic Berry velocity.
! **************************************************************************************************
   SUBROUTINE projected_current(h, s, positions, energies, current, status)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :), positions(:, :, :)
      REAL(dp), INTENT(OUT)                              :: energies(:)
      COMPLEX(dp), INTENT(OUT)                           :: current(:, :, :)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: c(:, :), moment(:, :)
      INTEGER                                            :: axis, i, j, n

      status = -1
      energies = 0.0_dp
      current = 0.0_dp
      n = SIZE(h, 1)
      IF (n < 1 .OR. SIZE(h, 2) /= n .OR. ANY(SHAPE(s) /= [n, n])) RETURN
      IF (SIZE(energies) /= n .OR. ANY(SHAPE(current) /= [n, n, 3])) RETURN
      IF (ANY(SHAPE(positions) /= [n, n, 3])) RETURN
      IF (.NOT. valid_hermitian(h) .OR. .NOT. valid_hermitian(s)) RETURN
      DO axis = 1, 3
         IF (.NOT. valid_hermitian(positions(:, :, axis))) RETURN
      END DO
      ALLOCATE (c(n, n), moment(n, n))
      CALL projected_eigensystem(h, s, energies, c, status)
      IF (status /= 0) RETURN
      status = -1
      DO axis = 1, 3
         moment = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(positions(:, :, axis), c))
         DO j = 1, n
            DO i = 1, n
               current(i, j, axis) = (energies(i) - energies(j))*moment(i, j)
            END DO
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(current, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(current)))) RETURN
      status = 0
   END SUBROUTINE projected_current
! **************************************************************************************************
!> \brief Complete generalized eigensystem with a positive, well-conditioned AO metric.
!> \param h Hermitian covariant Hamiltonian
!> \param s Hermitian positive definite overlap
!> \param energies Complete spectrum
!> \param c Metric-orthonormal eigenvectors
!> \param status 0 success; -1 invalid matrices; -2 metric failure; positive LAPACK failure
! **************************************************************************************************
   SUBROUTINE projected_eigensystem(h, s, energies, c, status)
      COMPLEX(dp), INTENT(IN) :: h(:, :), s(:, :)
      REAL(dp), INTENT(OUT) :: energies(:)
      COMPLEX(dp), INTENT(OUT) :: c(:, :)
      INTEGER, INTENT(OUT) :: status

      COMPLEX(dp), ALLOCATABLE :: b(:, :), work(:)
      REAL(dp), ALLOCATABLE :: rwork(:)
      INTEGER :: n, info
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      status = -1
      energies = 0.0_dp
      c = 0.0_dp
      n = SIZE(h, 1)
      IF (n < 1 .OR. SIZE(h, 2) /= n .OR. ANY(SHAPE(s) /= [n, n])) RETURN
      IF (SIZE(energies) /= n .OR. ANY(SHAPE(c) /= [n, n])) RETURN
      IF (.NOT. valid_hermitian(h) .OR. .NOT. valid_hermitian(s)) RETURN
      ALLOCATE (b(n, n), work(MAX(1, 2*n)), rwork(MAX(1, 3*n)))
      b = s
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zheev('N', 'U', n, b, n, energies, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 1
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(energies))) RETURN
      status = -2
      IF (MAXVAL(energies) <= 0.0_dp) RETURN
      IF (MINVAL(energies) <= 1.e-10_dp*MAXVAL(energies)) RETURN
      c = h
      b = s
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zhegv(1, 'V', 'U', n, c, n, b, n, energies, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 1
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(energies))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(c, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(c)))) RETURN
      status = 0
   END SUBROUTINE projected_eigensystem
! **************************************************************************************************
!> \brief Reject nonfinite or non-Hermitian data instead of discarding a triangle.
!> \param a Covariant matrix
!> \return Whether both components are finite and Hermitian within tolerance
! **************************************************************************************************
   LOGICAL FUNCTION valid_hermitian(a) RESULT(valid)
      COMPLEX(dp), INTENT(IN)                            :: a(:, :)

      valid = ALL(ieee_is_finite(REAL(a, dp))) .AND. ALL(ieee_is_finite(AIMAG(a)))
      IF (.NOT. valid) RETURN
      valid = MAXVAL(ABS(a - CONJG(TRANSPOSE(a)))) <= 1.e-10_dp*MAX(TINY(1.0_dp), MAXVAL(ABS(a)))
   END FUNCTION valid_hermitian
END MODULE kubo_projected

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gauge-covariant C1/C2 plaquette kernels, independent of the state representation.
!> \note C2 is the second Chern character pairing, not generally the second Chern class.
! **************************************************************************************************
MODULE topology_curvature
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE topology_wilson,                 ONLY: wilson_step

   IMPLICIT NONE
   PRIVATE
   REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp, unitary_tol = 1.e-8_dp
   PUBLIC :: polar_link, link_plaquette, curvature_density
   INTERFACE
      SUBROUTINE zgees(jobvs, sort, SELECT, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info)
         IMPORT :: dp
         CHARACTER(LEN=1), INTENT(IN) :: jobvs, sort
         INTEGER, INTENT(IN) :: n, lda, ldvs, lwork
         INTEGER, INTENT(OUT) :: sdim, info
         COMPLEX(KIND=dp), INTENT(INOUT) :: a(lda, *)
         COMPLEX(KIND=dp), INTENT(OUT) :: w(*), vs(ldvs, *), work(*)
         REAL(KIND=dp), INTENT(OUT) :: rwork(*)
         LOGICAL, INTENT(OUT) :: bwork(*)
         INTERFACE
            LOGICAL FUNCTION select(eigenvalue)
               IMPORT :: dp
               COMPLEX(KIND=dp), INTENT(IN) :: eigenvalue
            END FUNCTION select
         END INTERFACE
      END SUBROUTINE zgees
   END INTERFACE
CONTAINS

! **************************************************************************************************
!> \brief Polar-unitarize a physical link, retaining its minimum singular value.
!> \param overlap Raw selected-state overlap matrix
!> \param link Unitary polar factor
!> \param minimum Smallest singular value of the raw overlap
!> \param status Zero on success; nonzero for an invalid or singular link
!> \param tolerance Minimum allowed singular value
! **************************************************************************************************
   SUBROUTINE polar_link(overlap, link, minimum, status, tolerance)
      COMPLEX(KIND=dp), INTENT(IN)                       :: overlap(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: link(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: minimum
      INTEGER, INTENT(OUT)                               :: status
      REAL(KIND=dp), INTENT(IN)                          :: tolerance

      INTEGER                                            :: i

      link(:, :) = 0.0_dp
      DO i = 1, MIN(SIZE(link, 1), SIZE(link, 2))
         link(i, i) = 1.0_dp
      END DO
      CALL wilson_step(link, overlap, minimum, status, tolerance)
   END SUBROUTINE polar_link

! **************************************************************************************************
!> \brief Oriented plaquette U_mu(x) U_nu(x+mu) U_mu(x+nu)^dagger U_nu(x)^dagger.
!> \param u_mu Link in direction mu based at x
!> \param u_nu_at_mu Link in direction nu based at x+mu
!> \param u_mu_at_nu Link in direction mu based at x+nu
!> \param u_nu Link in direction nu based at x
!> \param p Oriented unitary plaquette based at x
! **************************************************************************************************
   SUBROUTINE link_plaquette(u_mu, u_nu_at_mu, u_mu_at_nu, u_nu, p)
      COMPLEX(KIND=dp), INTENT(IN)                       :: u_mu(:, :), u_nu_at_mu(:, :), &
                                                            u_mu_at_nu(:, :), u_nu(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: p(:, :)

      p(:, :) = MATMUL(MATMUL(u_mu, u_nu_at_mu), MATMUL(CONJG(TRANSPOSE(u_mu_at_nu)), CONJG(TRANSPOSE(u_nu))))
   END SUBROUTINE link_plaquette

! **************************************************************************************************
!> \brief Local contribution to C1 or C2. Report unrounded values and reject unresolved phases.
!> \param plaquettes Order 12 for C1; 12,13,14,23,24,34 for C2; all based at the same vertex
!> \param ndim Parameter dimension, 2 or 4
!> \param phase_limit Maximum allowed principal plaquette phase, strictly between zero and pi
!> \param value Unrounded local contribution
!> \param maximum Maximum absolute plaquette phase
!> \param status Zero on success, negative invalid/unresolved input, positive LAPACK failure
! **************************************************************************************************
   SUBROUTINE curvature_density(plaquettes, ndim, phase_limit, value, maximum, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: plaquettes(:, :, :)
      INTEGER, INTENT(IN)                                :: ndim
      REAL(KIND=dp), INTENT(IN)                          :: phase_limit
      REAL(KIND=dp), INTENT(OUT)                         :: value, maximum
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(KIND=dp), ALLOCATABLE                      :: a(:, :), e(:), f(:, :, :), term(:, :), &
                                                            v(:, :), work(:)
      INTEGER                                            :: i, j, n, np, sdim
      LOGICAL, ALLOCATABLE                               :: bwork(:)
      REAL(KIND=dp)                                      :: determinant_phase
      REAL(KIND=dp), ALLOCATABLE                         :: phase(:), rwork(:)

      status = -1
      value = 0.0_dp
      maximum = 0.0_dp
      IF (ndim /= 2 .AND. ndim /= 4) RETURN
      IF (.NOT. ieee_is_finite(phase_limit)) RETURN
      IF (phase_limit <= 0.0_dp .OR. phase_limit >= pi) RETURN
      n = SIZE(plaquettes, 1)
      np = ndim*(ndim - 1)/2
      IF (n < 1 .OR. SIZE(plaquettes, 2) /= n .OR. SIZE(plaquettes, 3) /= np) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(plaquettes, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(plaquettes)))) RETURN
      ALLOCATE (a(n, n), v(n, n), e(n), work(4*n), rwork(n), phase(n), bwork(n), f(n, n, np), term(n, n))
      DO j = 1, np
         a(:, :) = plaquettes(:, :, j)
         term(:, :) = MATMUL(CONJG(TRANSPOSE(a)), a)
         DO i = 1, n
            term(i, i) = term(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(term)) > unitary_tol) THEN
            status = -2
            RETURN
         END IF
         CALL zgees('V', 'N', select_none, n, a, n, sdim, e, v, n, work, SIZE(work), rwork, bwork, status)
         IF (status /= 0) RETURN
         phase(:) = ATAN2(AIMAG(e), REAL(e, dp))
         IF (ndim == 2) THEN
            determinant_phase = ATAN2(SIN(SUM(phase)), COS(SUM(phase)))
            maximum = ABS(determinant_phase)
            value = -determinant_phase/(2.0_dp*pi)
         ELSE
            maximum = MAX(maximum, MAXVAL(ABS(phase)))
            a(:, :) = v
            DO i = 1, n
               a(:, i) = CMPLX(0.0_dp, phase(i), dp)*v(:, i)
            END DO
            f(:, :, j) = MATMUL(a, CONJG(TRANSPOSE(v)))
         END IF
         IF (maximum >= phase_limit) THEN
            status = -3
            RETURN
         END IF
      END DO
      IF (ndim == 4) THEN
         term(:, :) = MATMUL(f(:, :, 1), f(:, :, 6)) - MATMUL(f(:, :, 2), f(:, :, 5)) + MATMUL(f(:, :, 3), f(:, :, 4))
         DO i = 1, n
            value = value - REAL(term(i, i), dp)/(4.0_dp*pi*pi)
         END DO
      END IF
      status = 0
   END SUBROUTINE curvature_density

! **************************************************************************************************
!> \brief Unused eigenvalue-selection callback for an unsorted complex Schur decomposition.
!> \param z Eigenvalue supplied by LAPACK
!> \return Always false; no eigenvalue reordering is requested
! **************************************************************************************************
   LOGICAL FUNCTION select_none(z)
      COMPLEX(KIND=dp), INTENT(IN)                       :: z

      select_none = .FALSE.
      IF (.NOT. ieee_is_finite(REAL(z, dp))) select_none = .FALSE.
   END FUNCTION select_none
END MODULE topology_curvature

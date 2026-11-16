!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inversion action on CP2K Gaussian Bloch coefficients, including lattice phases.
! **************************************************************************************************
MODULE topology_inversion
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gaussian_inversion_action
CONTAINS

! **************************************************************************************************
!> \brief Construct sparse AO inversion at a TRIM for the operation r -> 2*center-r.
!> \param qs_env Electronic environment
!> \param k Fractional reciprocal coordinate at a TRIM
!> \param center Fractional inversion center in the user's cell setting
!> \param mapping Target AO row for every source row
!> \param phase AO parity times exp(+2*pi*i*k.L), with 2c-r_a=r_b+L
!> \param tolerance Geometry tolerance in Bohr
!> \param status Zero on success, negative if geometry/basis/periodicity is incompatible
! **************************************************************************************************
   SUBROUTINE gaussian_inversion_action(qs_env, k, center, mapping, phase, tolerance, status)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: k(3), center(3)
      INTEGER, INTENT(OUT)                               :: mapping(:)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: phase(:)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: a, b, found, i, iset, l, n, shell, total
      INTEGER, ALLOCATABLE                               :: kind(:), offset(:), TARGET(:)
      REAL(KIND=dp)                                      :: angle, delta(3), image(3)
      REAL(KIND=dp), ALLOCATABLE                         :: position(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds_set(:)

      status = -1
      mapping = 0
      phase = CMPLX(0, 0, dp)
      IF (tolerance <= 0.0_dp .OR. SIZE(phase) /= SIZE(mapping)) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. .NOT. ALL(ieee_is_finite(center))) RETURN
      IF (.NOT. ALL(ieee_is_finite(k))) RETURN
      IF (MAXVAL(ABS(2*k - NINT(2*k))) > 1.e-8_dp) RETURN
      CALL get_qs_env(qs_env, cell=cell, particle_set=particles, qs_kind_set=kinds_set)
      IF (ANY(cell%perd /= 1)) RETURN
      n = SIZE(particles)
      ALLOCATE (KIND(n), offset(n), TARGET(n), position(3, n))
      total = 0
      DO a = 1, n
         CALL get_atomic_kind(particles(a)%atomic_kind, kind_number=KIND(a))
         CALL get_qs_kind(kinds_set(KIND(a)), basis_set=basis)
         IF (.NOT. ASSOCIATED(basis)) RETURN
         offset(a) = total
         total = total + basis%nsgf
         position(:, a) = MATMUL(cell%h_inv, pbc(particles(a)%r, cell))
      END DO
      IF (total /= SIZE(mapping)) RETURN
      DO a = 1, n
         found = 0
         DO b = 1, n
            IF (KIND(a) /= KIND(b)) CYCLE
            delta = 2*center - position(:, a) - position(:, b)
            IF (SQRT(SUM(MATMUL(cell%hmat, delta - NINT(delta))**2)) > tolerance) CYCLE
            found = found + 1
            TARGET(a) = b
         END DO
         IF (found /= 1) RETURN
      END DO
      DO a = 1, n
         b = TARGET(a)
         IF (TARGET(b) /= a) RETURN
         image = REAL(NINT(2*center - position(:, a) - position(:, b)), dp)
         angle = twopi*DOT_PRODUCT(k, image)
         CALL get_qs_kind(kinds_set(KIND(a)), basis_set=basis)
         DO iset = 1, basis%nset
            DO shell = 1, basis%nshell(iset)
               l = basis%l(shell, iset)
               DO i = basis%first_sgf(shell, iset), basis%last_sgf(shell, iset)
                  mapping(offset(a) + i) = offset(b) + i
                  phase(offset(a) + i) = (-1)**l*CMPLX(COS(angle), SIN(angle), dp)
               END DO
            END DO
         END DO
      END DO
      IF (ANY(mapping == 0)) RETURN
      status = 0
   END SUBROUTINE gaussian_inversion_action
END MODULE topology_inversion

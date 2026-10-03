!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Spherical USPP augmentation for atomic initial guesses and atomic SCF.
! **************************************************************************************************
MODULE atom_uspp
   USE atom_types,                      ONLY: atom_basis_gridrep,&
                                              atom_basis_type,&
                                              release_atom_basis
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE kinds,                           ONLY: dp
   USE splines,                         ONLY: spline3ders
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_uspp'
   PUBLIC :: atom_uspp_setup, atom_uspp_density, atom_uspp_integrate

CONTAINS

! **************************************************************************************************
!> \brief Prepare spherical projector products; add local augmentation and metric correction.
!> \param pot Scalar USPP in the original projector convention.
!> \param basis Atomic basis and radial density quadrature.
!> \param core Bare local matrix, augmented in place by the local ionic potential.
!> \param metric Gaussian overlap on entry; generalized overlap on return.
! **************************************************************************************************
   SUBROUTINE atom_uspp_setup(pot, basis, core, metric)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      REAL(dp), INTENT(INOUT)                            :: core(:, :, 0:), metric(:, :, 0:)

      INTEGER                                            :: first, i, j, k, l, m, n, nr, pair
      REAL(dp), ALLOCATABLE                              :: deriv(:), potential(:), qr(:), values(:)
      TYPE(atom_basis_type)                              :: source_basis
      TYPE(upf_projector_type)                           :: projectors

      CPASSERT(pot%is_ultrasoft .AND. .NOT. pot%is_paw .AND. .NOT. pot%has_so)
      CPASSERT(ALLOCATED(pot%qijl))
      CPASSERT(ALL(basis%grid%rad > 0.0_dp))
      CPASSERT(SIZE(pot%r) >= 4)
      IF (ALLOCATED(basis%uspp)) DEALLOCATE (basis%uspp)
      ALLOCATE (basis%uspp)
      m = pot%number_of_proj
      n = MAXVAL(basis%nbas)
      nr = basis%grid%nr
      ALLOCATE (basis%uspp%lbeta(m), basis%uspp%beta(n, m), &
                basis%uspp%radial(nr, m*(m + 1)/2, 0:1))
      basis%uspp%lbeta(:) = pot%lbeta
      basis%uspp%beta = 0.0_dp
      basis%uspp%radial = 0.0_dp
      CALL upf_projector_init(pot, 0, projectors)
      CALL atom_basis_gridrep(basis, source_basis, pot%r, pot%rab)
      DO i = 1, m
         l = pot%lbeta(i)
         CPASSERT(l <= UBOUND(basis%nbas, 1))
         DO j = 1, basis%nbas(l)
            basis%uspp%beta(j, i) = DOT_PRODUCT(projectors%weighted_beta(:, i), source_basis%bf(:, j, l))
         END DO
      END DO
      CALL release_atom_basis(source_basis)

      ALLOCATE (qr(SIZE(pot%r)), values(nr), deriv(nr), potential(nr))
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      DO j = 1, m
         DO i = 1, j
            IF (pot%lbeta(i) /= pot%lbeta(j)) CYCLE
            pair = j*(j - 1)/2 + i
            qr(first:) = pot%qijl(first:, pair, 0)/pot%r(first:)**2
            ! L=0 is smooth and even at the origin; infer its finite value from r**2.
            IF (first == 2) THEN
               qr(1) = (qr(2)*pot%r(3)**2 - qr(3)*pot%r(2)**2)/(pot%r(3)**2 - pot%r(2)**2)
            END IF
            CALL interpolate_radial(pot%r, qr, basis%grid%rad, values, deriv)
            WHERE (basis%grid%rad > pot%r(SIZE(pot%r)))
               values = 0.0_dp
               deriv = 0.0_dp
            END WHERE
            basis%uspp%radial(:, pair, 0) = values
            basis%uspp%radial(:, pair, 1) = deriv
         END DO
      END DO
      CALL interpolate_radial(pot%r, pot%vlocal, basis%grid%rad, potential, deriv)
      DO k = 1, nr
         IF (basis%grid%rad(k) > pot%r(SIZE(pot%r))) potential(k) = -pot%zion/basis%grid%rad(k)
      END DO
      CALL atom_uspp_integrate(core, potential, basis, 0)
      values = 1.0_dp
      ! The atomic metric uses the same interpolant/quadrature as its density.
      CALL atom_uspp_integrate(metric, values, basis, 0)
   END SUBROUTINE atom_uspp_setup

! **************************************************************************************************
!> \brief Interpolate a radial table and its derivative on either ordered atomic grid.
!> \param x Ascending tabulation coordinates.
!> \param y Tabulated values.
!> \param r Ascending or descending evaluation grid.
!> \param values Interpolated values.
!> \param deriv Derivative of the same interpolant.
! **************************************************************************************************
   SUBROUTINE interpolate_radial(x, y, r, values, deriv)
      REAL(dp), INTENT(IN)                               :: x(:), y(:), r(:)
      REAL(dp), INTENT(OUT)                              :: values(:), deriv(:)

      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: d(:), sorted(:), v(:)

      n = SIZE(r)
      ALLOCATE (sorted(n), v(n), d(n))
      IF (r(1) > r(n)) THEN
         sorted(:) = r(n:1:-1)
      ELSE
         sorted(:) = r
      END IF
      CALL spline3ders(x, y, sorted, ynew=v, dynew=d)
      IF (r(1) > r(n)) THEN
         values(:) = v(n:1:-1)
         deriv(:) = d(n:1:-1)
      ELSE
         values(:) = v
         deriv(:) = d
      END IF
   END SUBROUTINE interpolate_radial

! **************************************************************************************************
!> \brief Add the spherical augmentation to 4*pi times density or radial derivative.
!> \param density Radial quantity, incremented in place.
!> \param pmat Atomic density matrices, with occupations summed over magnetic states.
!> \param basis Atomic basis and prepared augmentation.
!> \param maxl Largest occupied angular momentum.
!> \param typ RHO or DER.
! **************************************************************************************************
   SUBROUTINE atom_uspp_density(density, pmat, basis, maxl, typ)
      REAL(dp), INTENT(INOUT)                            :: density(:)
      REAL(dp), INTENT(IN)                               :: pmat(:, :, 0:)
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      INTEGER, INTENT(IN)                                :: maxl
      CHARACTER(LEN=*), INTENT(IN)                       :: typ

      INTEGER                                            :: i, j, l, n, order, pair
      REAL(dp)                                           :: coefficient

      IF (.NOT. ALLOCATED(basis%uspp)) RETURN
      SELECT CASE (typ)
      CASE ('RHO')
         order = 0
      CASE ('DER')
         order = 1
      CASE DEFAULT
         CPABORT("Atomic USPP kinetic/Laplacian augmentation is not implemented")
      END SELECT
      DO j = 1, SIZE(basis%uspp%lbeta)
         l = basis%uspp%lbeta(j)
         IF (l > maxl) CYCLE
         n = basis%nbas(l)
         DO i = 1, j
            IF (basis%uspp%lbeta(i) /= l) CYCLE
            pair = j*(j - 1)/2 + i
            coefficient = DOT_PRODUCT(basis%uspp%beta(:n, i), &
                                      MATMUL(pmat(:n, :n, l), basis%uspp%beta(:n, j)))
            IF (i /= j) coefficient = 2.0_dp*coefficient
            density(:) = density + coefficient*basis%uspp%radial(:, pair, order)
         END DO
      END DO
   END SUBROUTINE atom_uspp_density

! **************************************************************************************************
!> \brief Add the adjoint spherical augmentation potential to an atomic AO matrix.
!> \param matrix Atomic matrix, incremented in place.
!> \param potential Potential on the same radial grid as the augmentation.
!> \param basis Atomic basis and prepared augmentation.
!> \param order Zero for density, one for its radial derivative.
! **************************************************************************************************
   SUBROUTINE atom_uspp_integrate(matrix, potential, basis, order)
      REAL(dp), INTENT(INOUT)                            :: matrix(:, :, 0:)
      REAL(dp), INTENT(IN)                               :: potential(:)
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      INTEGER, INTENT(IN)                                :: order

      INTEGER                                            :: i, j, k, l, n, pair
      REAL(dp)                                           :: coefficient

      IF (.NOT. ALLOCATED(basis%uspp)) RETURN
      CPASSERT(order == 0 .OR. order == 1)
      DO j = 1, SIZE(basis%uspp%lbeta)
         l = basis%uspp%lbeta(j)
         n = basis%nbas(l)
         DO i = 1, SIZE(basis%uspp%lbeta)
            IF (basis%uspp%lbeta(i) /= l) CYCLE
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            coefficient = SUM(potential*basis%grid%wr*basis%uspp%radial(:, pair, order))
            DO k = 1, n
               matrix(:n, k, l) = matrix(:n, k, l) + &
                                  coefficient*basis%uspp%beta(:n, i)*basis%uspp%beta(k, j)
            END DO
         END DO
      END DO
   END SUBROUTINE atom_uspp_integrate
END MODULE atom_uspp

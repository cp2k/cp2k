!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Half-numerical overlaps with the original scalar UPF radial projectors.
!>        Cartesian polynomials and angular products are integrated analytically.
!>        Only a radial quadrature remains; no fit or rotation of projector identities
!>        is performed. The UPF reader preserves the tabulated beta amplitude and
!>        includes the Ry-to-Hartree factor in the coupling matrix D.
! **************************************************************************************************
MODULE upf_projector_integrals
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT
   USE kinds,                           ONLY: dp
   USE libgrpp,                         ONLY: libgrpp_init
   USE mathconstants,                   ONLY: dfac,&
                                              fourpi,&
                                              maxfac
   USE mathlib,                         ONLY: binomial
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE upf_types,                       ONLY: atom_upfpot_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'upf_projector_integrals'

   TYPE upf_projector_type
      INTEGER :: degree = -1, lmax = -1, harmonic_max = -1, nproj = 0
      REAL(KIND=dp) :: radius = 0.0_dp
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:)
      REAL(KIND=dp), ALLOCATABLE :: r(:), weighted_beta(:, :), d(:, :), gaunt(:, :, :)
   END TYPE upf_projector_type

   PUBLIC :: upf_projector_type, upf_projector_init, upf_projector_overlap

   INTERFACE
      FUNCTION scaled_bessel(n, x) RESULT(value) BIND(C, name="libgrpp_modified_bessel_scaled")
         IMPORT :: C_DOUBLE, C_INT
         INTEGER(C_INT), VALUE :: n
         REAL(C_DOUBLE), VALUE :: x
         REAL(C_DOUBLE) :: value
      END FUNCTION scaled_bessel
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Prepare immutable radial and angular data outside parallel regions.
!> \param pot Original UPF data, including beta(r), not r*beta(r).
!> \param degree Largest Cartesian degree needed (AO degree + derivatives/moments).
!> \param projectors Prepared data; reinitializing also releases the previous allocation.
! **************************************************************************************************
   SUBROUTINE upf_projector_init(pot, degree, projectors)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_projector_type), INTENT(OUT)              :: projectors

      INTEGER                                            :: a, b, i, j, k, l, lambda, m, n, &
                                                            nsimpson, p, powers(3), q, s
      REAL(KIND=dp), ALLOCATABLE                         :: moments(:, :), weights(:)

      CPASSERT(degree >= 0)
      IF (pot%has_so) THEN
         CPABORT("Spin-orbit UPF projectors are not supported")
      END IF
      CPASSERT(pot%number_of_proj > 0)
      CPASSERT(ALL(pot%lbeta >= 0))
      projectors%degree = degree
      projectors%lmax = MAXVAL(pot%lbeta)
      projectors%harmonic_max = degree + projectors%lmax
      ! Both the factorial and bundled Bessel tables cover this range.
      CPASSERT(projectors%harmonic_max <= MIN(20, maxfac/2))
      CALL init_orbital_pointers(projectors%harmonic_max)
      CALL init_spherical_harmonics(projectors%harmonic_max, -1)
      CALL libgrpp_init()

      n = SIZE(pot%r)
      CPASSERT(n >= 3)
      CPASSERT(SIZE(pot%rab) == n)
      CPASSERT(SIZE(pot%beta, 1) == n)
      CPASSERT(SIZE(pot%beta, 2) == pot%number_of_proj)
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)))
      CPASSERT(ALL(pot%rab > 0.0_dp))
      ALLOCATE (projectors%r(n), projectors%lbeta(pot%number_of_proj))
      projectors%r(:) = pot%r
      projectors%lbeta(:) = pot%lbeta
      ALLOCATE (weights(n))
      weights = 0.0_dp
      ! Composite Simpson in the mesh index; PP_RAB is dr/d(index).
      ! An even mesh has one final trapezoidal interval.
      nsimpson = n - 1 + MODULO(n, 2)
      weights(1:nsimpson:2) = 2.0_dp/3.0_dp
      weights(2:nsimpson - 1:2) = 4.0_dp/3.0_dp
      weights(1) = 1.0_dp/3.0_dp
      weights(nsimpson) = 1.0_dp/3.0_dp
      IF (nsimpson < n) weights(n - 1:n) = weights(n - 1:n) + 0.5_dp
      weights(:) = weights*pot%rab*pot%r**2
      ALLOCATE (projectors%weighted_beta(n, pot%number_of_proj))
      DO i = 1, pot%number_of_proj
         projectors%weighted_beta(:, i) = weights*pot%beta(:, i)
         DO j = n, 1, -1
            IF (pot%beta(j, i) /= 0.0_dp) THEN
               projectors%radius = MAX(projectors%radius, pot%r(j))
               EXIT
            END IF
         END DO
      END DO

      projectors%nproj = SUM(2*pot%lbeta + 1)
      ALLOCATE (projectors%radial_index(projectors%nproj), projectors%angular_index(projectors%nproj))
      p = 0
      DO i = 1, pot%number_of_proj
         l = pot%lbeta(i)
         DO m = 1, 2*l + 1
            p = p + 1
            projectors%radial_index(p) = i
            projectors%angular_index(p) = l*l + m
         END DO
      END DO
      ALLOCATE (projectors%d(projectors%nproj, projectors%nproj))
      projectors%d = 0.0_dp
      DO p = 1, projectors%nproj
         DO q = 1, projectors%nproj
            i = projectors%radial_index(p)
            j = projectors%radial_index(q)
            IF (projectors%angular_index(p) == projectors%angular_index(q)) THEN
               projectors%d(p, q) = pot%dion(i, j)
            ELSE IF (pot%lbeta(i) /= pot%lbeta(j)) THEN
               CPASSERT(ABS(pot%dion(i, j)) < 1.e-12_dp)
            END IF
         END DO
      END DO

      ! G(s,lm,lambda*nu) = integral Omega**s Y_lm(Omega) Y_lambda,nu(Omega).
      ALLOCATE (projectors%gaunt(ncoset(degree), (projectors%lmax + 1)**2, &
                                 (projectors%harmonic_max + 1)**2))
      projectors%gaunt = 0.0_dp
      DO s = 1, ncoset(degree)
         k = SUM(indco(:, s))
         DO l = 0, projectors%lmax
            IF (.NOT. ANY(pot%lbeta == l)) CYCLE
            DO lambda = MAX(0, l - k), l + k
               IF (MODULO(k + l + lambda, 2) /= 0) CYCLE
               ALLOCATE (moments(nco(l), nco(lambda)))
               DO b = 1, nco(lambda)
                  DO a = 1, nco(l)
                     powers = indco(:, s) + indco(:, ncoset(l - 1) + a) + indco(:, ncoset(lambda - 1) + b)
                     moments(a, b) = 0.0_dp
                     IF (ALL(MODULO(powers, 2) == 0)) THEN
                        moments(a, b) = fourpi*PRODUCT(dfac(powers - 1))/dfac(SUM(powers) + 1)
                     END IF
                  END DO
               END DO
               projectors%gaunt(s, l*l + 1:(l + 1)**2, lambda*lambda + 1:(lambda + 1)**2) = &
                  MATMUL(MATMUL(orbtramat(l)%slm, moments), TRANSPOSE(orbtramat(lambda)%slm))
               DEALLOCATE (moments)
            END DO
         END DO
      END DO
   END SUBROUTINE upf_projector_init

! **************************************************************************************************
!> \brief Evaluate <Cartesian Gaussian|beta_i Y_lm> and AO-center derivatives.
!> \param projectors Immutable prepared data.
!> \param la Largest AO Cartesian degree (all degrees 0..la are returned).
!> \param zeta Primitive Gaussian exponents.
!> \param rac Projector center minus AO center, as in ai_overlap.
!> \param nder Maximum derivative degree, or moment degree in moment mode.
!> \param integrals Primitive-major AO order, original radial/m order, coset derivative order.
!> \param moment_origin Optional AO center minus moment reference; requests moments
!>        instead of derivatives, in the same coset ordering including the overlap.
! **************************************************************************************************
   SUBROUTINE upf_projector_overlap(projectors, la, zeta, rac, nder, integrals, moment_origin)
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: zeta(:), rac(3)
      INTEGER, INTENT(IN)                                :: nder
      REAL(KIND=dp), INTENT(OUT)                         :: integrals(:, :, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: moment_origin(3)

      INTEGER                                            :: a(3), axis, b(3), degree, deriv(3), i, &
                                                            id, ip, j, k, l, lambda, nb, nc, nh, &
                                                            nr, p, power, s, shifted(3), step, t
      REAL(KIND=dp)                                      :: alpha, direction(3), distance, factor
      REAL(KIND=dp), ALLOCATABLE :: angular(:, :, :), cartesian(:, :), centered(:, :), &
         coeff(:, :), envelope(:), harmonics(:), monomials(:), next_coeff(:), radial(:, :, :), &
         radial_power(:), radial_work(:), shift(:, :)

      degree = la + nder
      CPASSERT(la >= 0 .AND. nder >= 0)
      CPASSERT(degree <= projectors%degree)
      CPASSERT(ALL(zeta > 0.0_dp))
      CPASSERT(SIZE(integrals, 1) == SIZE(zeta)*ncoset(la))
      CPASSERT(SIZE(integrals, 2) == projectors%nproj)
      CPASSERT(SIZE(integrals, 3) == ncoset(nder))
      nc = ncoset(degree)
      nh = degree + projectors%lmax
      nr = SIZE(projectors%r)
      nb = SIZE(projectors%weighted_beta, 2)
      ALLOCATE (harmonics((nh + 1)**2), monomials(ncoset(nh)))
      distance = NORM2(rac)
      direction = [0.0_dp, 0.0_dp, 1.0_dp]
      IF (distance > 0.0_dp) direction = -rac/distance
      DO i = 1, ncoset(nh)
         monomials(i) = PRODUCT(direction**indco(:, i))
      END DO
      DO l = 0, nh
         harmonics(l*l + 1:(l + 1)**2) = &
            MATMUL(orbtramat(l)%slm, monomials(ncoset(l - 1) + 1:ncoset(l)))
      END DO
      ALLOCATE (angular(nc, (projectors%lmax + 1)**2, 0:nh))
      DO lambda = 0, nh
         DO i = 1, SIZE(angular, 2)
            angular(:, i, lambda) = MATMUL( &
                                    projectors%gaunt(:nc, i, lambda*lambda + 1:(lambda + 1)**2), &
                                    harmonics(lambda*lambda + 1:(lambda + 1)**2))
         END DO
      END DO
      ALLOCATE (shift(nc, nc))
      shift = 0.0_dp
      DO t = 1, nc
         a = indco(:, t)
         DO s = 1, nc
            b = indco(:, s)
            IF (ANY(b > a)) CYCLE
            factor = 1.0_dp
            DO axis = 1, 3
               factor = factor*binomial(a(axis), b(axis))*rac(axis)**(a(axis) - b(axis))
            END DO
            shift(t, s) = factor
         END DO
      END DO
      ALLOCATE (radial(0:degree, nb, 0:nh), radial_work(nr), radial_power(nr), envelope(nr))
      ALLOCATE (centered(nc, projectors%nproj), cartesian(nc, projectors%nproj))
      ALLOCATE (coeff(0:degree, 3), next_coeff(0:degree))
      integrals = 0.0_dp
      DO ip = 1, SIZE(zeta)
         alpha = zeta(ip)
         envelope(:) = EXP(-alpha*(projectors%r - distance)**2)
         DO lambda = 0, nh
            DO i = 1, nr
               radial_work(i) = envelope(i)*scaled_bessel(INT(lambda, C_INT), &
                                                          REAL(2.0_dp*alpha*projectors%r(i)*distance, C_DOUBLE))
            END DO
            radial_power(:) = radial_work
            DO k = 0, degree
               radial(k, :, lambda) = MATMUL(radial_power, projectors%weighted_beta)
               radial_power(:) = radial_power*projectors%r
            END DO
         END DO
         centered = 0.0_dp
         DO p = 1, projectors%nproj
            i = projectors%radial_index(p)
            j = projectors%angular_index(p)
            DO s = 1, nc
               k = SUM(indco(:, s))
               centered(s, p) = fourpi*DOT_PRODUCT(angular(s, j, :), radial(k, i, :))
            END DO
         END DO
         cartesian(:, :) = MATMUL(shift, centered)
         DO id = 1, ncoset(nder)
            deriv = indco(:, id)
            DO t = 1, ncoset(la)
               a = indco(:, t)
               coeff = 0.0_dp
               DO axis = 1, 3
                  coeff(a(axis), axis) = 1.0_dp
                  IF (PRESENT(moment_origin)) THEN
                     coeff(:, axis) = 0.0_dp
                     DO power = 0, deriv(axis)
                        coeff(a(axis) + power, axis) = binomial(deriv(axis), power)* &
                                                       moment_origin(axis)**(deriv(axis) - power)
                     END DO
                  ELSE
                     DO step = 1, deriv(axis)
                        next_coeff = 0.0_dp
                        DO power = 0, a(axis) + step - 1
                           next_coeff(power + 1) = next_coeff(power + 1) + 2.0_dp*alpha*coeff(power, axis)
                           IF (power > 0) next_coeff(power - 1) = next_coeff(power - 1) - power*coeff(power, axis)
                        END DO
                        coeff(:, axis) = next_coeff
                     END DO
                  END IF
               END DO
               DO s = 1, nc
                  shifted = indco(:, s)
                  factor = coeff(shifted(1), 1)*coeff(shifted(2), 2)*coeff(shifted(3), 3)
                  IF (factor == 0.0_dp) CYCLE
                  integrals((ip - 1)*ncoset(la) + t, :, id) = integrals((ip - 1)*ncoset(la) + t, :, id) + &
                                                              factor*cartesian(s, :)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE upf_projector_overlap

END MODULE upf_projector_integrals

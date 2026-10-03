!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Coulomb coupling of UPF augmentation to compensated Gaussian GAPW multipoles.
!>        The source radial quadrature is identical to the augmentation/overlap quadrature.
!>        These geometry-dependent integrals are building blocks, not an SCF driver.
! **************************************************************************************************
MODULE gapw_coulomb_integrals
   USE gamma,                           ONLY: fgamma_0,&
                                              fgamma_ref,&
                                              init_md_ftable
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE spherical_harmonics,             ONLY: y_lm
   USE uspp_augmentation,               ONLY: uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gapw_coulomb_integrals'
   PUBLIC :: gapw_coulomb_init, gapw_compensated_potential, gapw_augmentation_integrals
   PUBLIC :: gapw_periodic_augmentation_integrals

CONTAINS

! **************************************************************************************************
!> \brief Initialize shared angular/Boys tables before entering parallel evaluation.
!> \param maxl Largest Gaussian multipole angular momentum.
!> \param maxpower Largest additional power of r**2.
! **************************************************************************************************
   SUBROUTINE gapw_coulomb_init(maxl, maxpower)
      INTEGER, INTENT(IN)                                :: maxl, maxpower

      CPASSERT(MIN(maxl, maxpower) >= 0)
      CALL init_orbital_pointers(maxl)
      CALL init_spherical_harmonics(maxl, -1)
      CALL init_md_ftable(maxl + maxpower + 2)
   END SUBROUTINE gapw_coulomb_init

! **************************************************************************************************
!> \brief Potential of [r**(2*power) exp(-alpha*r**2) - c exp(-eta*r**2)] r**l Y_lm.
!>        c cancels the multipole exactly. Both density and potential use normalized
!>        real spherical harmonics in CP2K ordering. No normalization of the Gaussian
!>        itself is included. A complementary-tail expression avoids subtracting
!>        equal long-range multipole potentials. Values/derivatives are regular at r=0.
!> \param alpha Gaussian exponent.
!> \param power Additional nonnegative radial power of r**2.
!> \param l Angular momentum.
!> \param eta Compensation Gaussian exponent.
!> \param displacement Evaluation point minus Gaussian center, in bohr.
!> \param potential Values for m=-l..l.
!> \param gradient Cartesian derivatives with respect to displacement.
!> \param hessian Optional second Cartesian derivatives.
! **************************************************************************************************
   SUBROUTINE gapw_compensated_potential(alpha, power, l, eta, displacement, potential, gradient, hessian)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: eta, displacement(3)
      REAL(dp), INTENT(OUT)                              :: potential(2*l + 1), gradient(3, 2*l + 1)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: hessian(3, 3, 2*l + 1)

      INTEGER                                            :: a, axis, axis2, j, m, powers(3), &
                                                            reduced(3)
      LOGICAL                                            :: tail
      REAL(dp)                                           :: angular, c, dangular(3), &
                                                            ddangular(3, 3), ddh, ddh0, dh, dh0, &
                                                            factor, h, h0, s, term

      CPASSERT(alpha > 0.0_dp .AND. eta > 0.0_dp .AND. MIN(l, power) >= 0)
      potential = 0.0_dp
      gradient = 0.0_dp
      IF (PRESENT(hessian)) hessian = 0.0_dp
      IF (alpha == eta .AND. power == 0) RETURN
      c = (eta/alpha)**(REAL(l, dp) + 1.5_dp)
      DO j = 1, power
         c = c*(REAL(l + j, dp) + 0.5_dp)/alpha
      END DO
      s = SUM(displacement**2)
      ! Switch beyond the radial moment's peak. Low multipoles require the
      ! tail form earlier, before their residual is lost to cancellation.
      tail = MIN(alpha, eta)*s > REAL(l + power + 2, dp)
      CALL gaussian_coulomb_radial(alpha, power, l, s, tail, h, dh, ddh)
      CALL gaussian_coulomb_radial(eta, 0, l, s, tail, h0, dh0, ddh0)
      h = h - c*h0
      dh = dh - c*dh0
      ddh = ddh - c*ddh0
      DO m = 1, 2*l + 1
         angular = 0.0_dp
         dangular = 0.0_dp
         ddangular = 0.0_dp
         DO a = 1, nco(l)
            powers = indco(:, ncoset(l - 1) + a)
            factor = orbtramat(l)%slm(m, a)
            IF (factor == 0.0_dp) CYCLE
            angular = angular + factor*PRODUCT(displacement**powers)
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               term = factor*REAL(powers(axis), dp)
               dangular(axis) = dangular(axis) + term*PRODUCT(displacement**reduced)
               IF (PRESENT(hessian)) THEN
                  DO axis2 = 1, 3
                     IF (reduced(axis2) == 0) CYCLE
                     reduced(axis2) = reduced(axis2) - 1
                     ddangular(axis, axis2) = ddangular(axis, axis2) + &
                                              term*REAL(reduced(axis2) + 1, dp)*PRODUCT(displacement**reduced)
                     reduced(axis2) = reduced(axis2) + 1
                  END DO
               END IF
            END DO
         END DO
         potential(m) = h*angular
         gradient(:, m) = 2.0_dp*displacement*dh*angular + h*dangular
         IF (PRESENT(hessian)) THEN
            DO axis2 = 1, 3
               DO axis = 1, 3
                  hessian(axis, axis2, m) = 4.0_dp*displacement(axis)*displacement(axis2)*ddh*angular + &
                                            2.0_dp*dh*(displacement(axis)*dangular(axis2) + &
                                                       displacement(axis2)*dangular(axis)) + &
                                            h*ddangular(axis, axis2)
               END DO
               hessian(axis2, axis2, m) = hessian(axis2, axis2, m) + 2.0_dp*dh*angular
            END DO
         END IF
      END DO
   END SUBROUTINE gapw_compensated_potential

! **************************************************************************************************
!> \brief Regular radial potential coefficient and its first two r**2 derivatives.
!> \param alpha Gaussian exponent.
!> \param power Additional power of r**2.
!> \param l Angular momentum.
!> \param s Squared radius.
!> \param tail Subtract the asymptotic multipole; valid only away from the origin.
!> \param h Coefficient multiplying the solid harmonic.
!> \param dh First derivative with respect to s.
!> \param ddh Second derivative with respect to s.
! **************************************************************************************************
   SUBROUTINE gaussian_coulomb_radial(alpha, power, l, s, tail, h, dh, ddh)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: s
      LOGICAL, INTENT(IN)                                :: tail
      REAL(dp), INTENT(OUT)                              :: h, dh, ddh

      INTEGER                                            :: j, n
      REAL(dp)                                           :: density, ex, f(0:l + power + 2), outer, &
                                                            upper, x

      x = alpha*s
      ex = EXP(-x)
      outer = 0.5_dp*ex/alpha
      DO j = 1, power
         outer = (REAL(j, dp)*outer + 0.5_dp*s**j*ex)/alpha
      END DO
      n = l + power + 1
      IF (tail) THEN
         ! Half-integer upper incomplete gamma: only exponentially small terms remain.
         upper = 0.5_dp*SQRT(pi)*ERFC(SQRT(x))
         DO j = 1, n
            upper = (REAL(j, dp) - 0.5_dp)*upper + 0.5_dp*x**(REAL(j, dp) - 0.5_dp)*ex
         END DO
         upper = upper/alpha**(REAL(n, dp) + 0.5_dp)
         density = s**power*ex
         h = fourpi/REAL(2*l + 1, dp)*(outer - upper/s**(REAL(l, dp) + 0.5_dp))
         dh = 2.0_dp*pi*upper/s**(REAL(l, dp) + 1.5_dp)
         ddh = -pi*REAL(2*l + 3, dp)*upper/s**(REAL(l, dp) + 2.5_dp) - pi*density/s
      ELSE
         IF (x <= 12.0_dp) THEN
            CALL fgamma_0(n + 1, x, f)
         ELSE IF (x >= REAL(n + 1, dp)) THEN
            ! The large-x polynomial approximation in fgamma_0 is unsuitable for
            ! differentiating a small difference of compensated potentials.
            f(0) = 0.5_dp*SQRT(pi/x)*ERF(SQRT(x))
            DO j = 1, n + 1
               f(j) = (REAL(2*j - 1, dp)*f(j - 1) - ex)/(2.0_dp*x)
            END DO
         ELSE
            f = fgamma_ref(n + 1, x)
         END IF
         h = fourpi/REAL(2*l + 1, dp)*(s**(power + 1)*f(n) + outer)
         dh = -2.0_dp*pi*s**power*f(n)
         ddh = 2.0_dp*pi*alpha*s**power*f(n + 1)
         IF (power > 0) ddh = ddh - 2.0_dp*pi*REAL(power, dp)*s**(power - 1)*f(n)
      END IF
   END SUBROUTINE gaussian_coulomb_radial

! **************************************************************************************************
!> \brief Integrate original UPF augmentation against compensated Gaussian multipole potentials.
!>        A caller can contract these bilinear integrals with projector occupations and
!>        Gaussian product coefficients. The same tensor supplies both density adjoints.
!>        The displacement derivative gives opposite forces on the two centers.
!>        Periodic images and Gaussian product/Gaunt coefficients belong to the caller.
!> \param augmentation Original radial weights and angular projector data.
!> \param displacement Augmentation center minus Gaussian center.
!> \param alpha Gaussian exponent.
!> \param power Additional Gaussian radial power.
!> \param l Gaussian multipole angular momentum.
!> \param eta Compensation exponent.
!> \param directions Unit vectors of the angular quadrature.
!> \param weights Angular weights, summing to 4*pi.
!> \param integral Projector-pair/multipole integrals.
!> \param derivative Optional derivatives with respect to displacement.
! **************************************************************************************************
   SUBROUTINE gapw_augmentation_integrals(augmentation, displacement, alpha, power, l, eta, &
                                          directions, weights, integral, derivative)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3), alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: eta, directions(:, :), weights(:)
      REAL(dp), INTENT(OUT)                              :: integral(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative(:, :, :, :)

      INTEGER                                            :: axis, i, ia, ir, j, la, lm, m, n, p, &
                                                            pair, q
      REAL(dp)                                           :: gradient(3, 2*l + 1), point(3), &
                                                            potential(2*l + 1), radial
      REAL(dp), ALLOCATABLE                              :: dprojected(:, :, :), harmonics(:, :), &
                                                            projected(:, :), values(:)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(integral) == [n, n, 2*l + 1]))
      CPASSERT(SIZE(directions, 1) == 3 .AND. SIZE(directions, 2) == SIZE(weights))
      CPASSERT(ABS(SUM(weights) - fourpi) < 1.e-10_dp)
      CPASSERT(MAXVAL(ABS(SUM(directions**2, DIM=1) - 1.0_dp)) < 1.e-12_dp)
      IF (PRESENT(derivative)) THEN
         CPASSERT(ALL(SHAPE(derivative) == [n, n, 2*l + 1, 3]))
      END IF
      ALLOCATE (harmonics(SIZE(weights), (augmentation%lmax + 1)**2), values(SIZE(weights)), &
                projected((augmentation%lmax + 1)**2, 2*l + 1), &
                dprojected((augmentation%lmax + 1)**2, 2*l + 1, 3))
      DO la = 0, augmentation%lmax
         DO m = -la, la
            lm = la*la + la + m + 1
            CALL y_lm(directions, values, la, m)
            harmonics(:, lm) = weights*values
         END DO
      END DO
      integral = 0.0_dp
      IF (PRESENT(derivative)) derivative = 0.0_dp
      DO ir = 1, SIZE(augmentation%weighted_q, 1)
         IF (ALL(augmentation%weighted_q(ir, :, :) == 0.0_dp)) CYCLE
         projected = 0.0_dp
         dprojected = 0.0_dp
         DO ia = 1, SIZE(weights)
            point = displacement + augmentation%r(ir)*directions(:, ia)
            CALL gapw_compensated_potential(alpha, power, l, eta, point, potential, gradient)
            DO m = 1, 2*l + 1
               projected(:, m) = projected(:, m) + harmonics(ia, :)*potential(m)
               IF (PRESENT(derivative)) THEN
                  DO axis = 1, 3
                     dprojected(:, m, axis) = dprojected(:, m, axis) + harmonics(ia, :)*gradient(axis, m)
                  END DO
               END IF
            END DO
         END DO
         DO q = 1, n
            j = augmentation%radial_index(q)
            DO p = 1, q
               i = augmentation%radial_index(p)
               pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
               DO la = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
                  augmentation%lbeta(i) + augmentation%lbeta(j), 2
                  radial = augmentation%weighted_q(ir, pair, la)
                  IF (radial == 0.0_dp) CYCLE
                  DO lm = la*la + 1, (la + 1)**2
                     IF (augmentation%gaunt(p, q, lm) == 0.0_dp) CYCLE
                     integral(p, q, :) = integral(p, q, :) + &
                                         radial*augmentation%gaunt(p, q, lm)*projected(lm, :)
                     IF (PRESENT(derivative)) THEN
                        derivative(p, q, :, :) = derivative(p, q, :, :) + &
                                                 radial*augmentation%gaunt(p, q, lm)*dprojected(lm, :, :)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO q = 1, n
         DO p = 1, q - 1
            integral(q, p, :) = integral(p, q, :)
            IF (PRESENT(derivative)) derivative(q, p, :, :) = derivative(p, q, :, :)
         END DO
      END DO
   END SUBROUTINE gapw_augmentation_integrals

! **************************************************************************************************
!> \brief Periodic compensated integrals and consistent displacement/strain derivatives.
!>        Gaussian tails are screened using an absolute primitive-potential envelope, not the
!>        current density coefficients. The original UPF radial mesh is never shortened.
!> \param augmentation Original radial and angular augmentation data.
!> \param displacement Augmentation center minus Gaussian center.
!> \param alpha Gaussian exponent.
!> \param power Additional radial power of r**2.
!> \param l Gaussian multipole angular momentum.
!> \param eta Compensation exponent.
!> \param hmat Periodic cell in bohr.
!> \param h_inv Inverse cell matrix.
!> \param directions Unit vectors of the angular rule.
!> \param weights Angular quadrature weights, summing to 4*pi.
!> \param tolerance Absolute tolerance for screening compensated primitive tails.
!> \param tensor Potential, three displacement derivatives and nine strain derivatives.
!> \param periodicity Periodic lattice directions (zero or one); defaults to three dimensions.
! **************************************************************************************************
   SUBROUTINE gapw_periodic_augmentation_integrals(augmentation, displacement, alpha, power, l, eta, &
                                                   hmat, h_inv, directions, weights, tolerance, tensor, periodicity)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(dp), INTENT(IN)                               :: displacement(3), alpha
      INTEGER, INTENT(IN)                                :: power, l
      REAL(dp), INTENT(IN)                               :: eta, hmat(3, 3), h_inv(3, 3), &
                                                            directions(:, :), weights(:), tolerance
      REAL(dp), INTENT(OUT)                              :: tensor(:, :, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: periodicity(3)

      INTEGER                                            :: axis, axis2, i, ir, ix, iy, iz, &
                                                            lower(3), n, periodic(3), upper(3)
      REAL(dp)                                           :: c, ddh, dh, dh0, envelope, extent, &
                                                            fractional(3), h, h0, padding, &
                                                            point(3), qnorm, radius, rmax, s, &
                                                            shell_factor
      REAL(dp), ALLOCATABLE                              :: derivative(:, :, :, :), integral(:, :, :)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(tensor) == [n, n, 2*l + 1, 13]))
      CPASSERT(tolerance > 0.0_dp)
      periodic = 1
      IF (PRESENT(periodicity)) periodic = periodicity
      CPASSERT(ALL(periodic == 0 .OR. periodic == 1))
      tensor = 0.0_dp
      IF (alpha == eta .AND. power == 0) RETURN
      rmax = 0.0_dp
      DO ir = 1, SIZE(augmentation%weighted_q, 1)
         IF (ANY(augmentation%weighted_q(ir, :, :) /= 0.0_dp)) rmax = MAX(rmax, augmentation%r(ir))
      END DO
      IF (rmax == 0.0_dp) RETURN
      c = (eta/alpha)**(REAL(l, dp) + 1.5_dp)
      DO i = 1, power
         c = c*(REAL(l + i, dp) + 0.5_dp)/alpha
      END DO
      qnorm = MAX(1.0_dp, SUM(ABS(augmentation%weighted_q))*SUM(ABS(augmentation%gaunt))* &
                  SQRT(fourpi*REAL(2*augmentation%lmax + 1, dp)))
      ! Start beyond all polynomial peaks; use separate absolute Gaussian tails so
      ! an accidental cancellation of two exponents cannot underestimate the support.
      radius = SQRT(8.0_dp*REAL(l + power + 2, dp)/MIN(alpha, eta))
      DO
         s = radius**2
         CALL gaussian_coulomb_radial(alpha, power, l, s, .TRUE., h, dh, ddh)
         CALL gaussian_coulomb_radial(eta, 0, l, s, .TRUE., h0, dh0, ddh)
         envelope = SQRT(REAL(2*l + 1, dp)/fourpi)*radius**l* &
                    (REAL(2*l + 2, dp)*(ABS(h) + c*ABS(h0)) + 2.0_dp*s*(ABS(dh) + c*ABS(dh0)))
         shell_factor = (1.0_dp + 2.0_dp*(radius + rmax)*MAXVAL(SQRT(SUM(h_inv**2, DIM=2))))**3
         IF (qnorm*shell_factor*envelope < tolerance) EXIT
         radius = 1.2_dp*radius
      END DO
      extent = radius + rmax
      fractional = MATMUL(h_inv, displacement)
      fractional = fractional - REAL(periodic, dp)*ANINT(fractional)
      lower = 0
      upper = 0
      DO axis = 1, 3
         IF (periodic(axis) == 0) CYCLE
         padding = extent*SQRT(SUM(h_inv(axis, :)**2))
         lower(axis) = CEILING(-fractional(axis) - padding)
         upper(axis) = FLOOR(-fractional(axis) + padding)
      END DO
      ALLOCATE (integral(n, n, 2*l + 1), derivative(n, n, 2*l + 1, 3))
      DO iz = lower(3), upper(3)
         DO iy = lower(2), upper(2)
            DO ix = lower(1), upper(1)
               point = MATMUL(hmat, fractional + REAL([ix, iy, iz], dp))
               IF (SUM(point**2) > extent**2) CYCLE
               CALL gapw_augmentation_integrals(augmentation, point, alpha, power, l, eta, &
                                                directions, weights, integral, derivative)
               tensor(:, :, :, 1) = tensor(:, :, :, 1) + integral
               DO axis = 1, 3
                  tensor(:, :, :, 1 + axis) = tensor(:, :, :, 1 + axis) + derivative(:, :, :, axis)
                  DO axis2 = 1, 3
                     i = 4 + axis + 3*(axis2 - 1)
                     tensor(:, :, :, i) = tensor(:, :, :, i) + derivative(:, :, :, axis)*point(axis2)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE gapw_periodic_augmentation_integrals

END MODULE gapw_coulomb_integrals

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM gapw_coulomb_integrals_unittest
   USE gapw_coulomb_integrals,          ONLY: gapw_augmentation_integrals,&
                                              gapw_compensated_potential,&
                                              gapw_coulomb_init,&
                                              gapw_periodic_augmentation_integrals
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE spherical_harmonics,             ONLY: y_lm
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   REAL(dp) :: max_error

   CALL cp2k_init()
   CALL gapw_coulomb_init(6, 4)
   max_error = 0.0_dp
   CALL check_multipoles()
   CALL check_augmentation()
   CALL check_periodic()
   PRINT '(A,ES16.7)', 'GAPW compensated Coulomb maximum scaled error: ', max_error
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_multipoles()
      REAL(dp), PARAMETER                                :: alpha = 1.3_dp, eta = 0.55_dp, &
                                                            step = 1.e-3_dp

      INTEGER                                            :: axis, j, l, m, n, power, sample
      REAL(dp) :: angular(1), c, density, direction(3, 1), gm(3, 13), gmm(3, 13), gp(3, 13), &
         gpp(3, 13), grad(3, 13), hess(3, 3, 13), higher(13), point(3), points(3, 7), r, &
         reference, s, shifted(3), tail_difference, value(13), vm(13), vmm(13), vp(13), vpp(13)

      points(:, 1) = 0.0_dp
      points(:, 2) = [1.e-8_dp, -2.e-8_dp, 3.e-8_dp]
      points(:, 3) = [0.31_dp, -0.52_dp, 0.73_dp]
      points(:, 4) = [0.0_dp, 0.0_dp, 1.2_dp]
      points(:, 7) = [8.0_dp, -3.0_dp, 2.0_dp]
      DO l = 0, 6
         n = 2*l + 1
         DO power = 0, 3
            points(:, 5) = [SQRT(REAL(l + power + 2, dp)/eta)*(1.0_dp - 1.e-8_dp), 0.0_dp, 0.0_dp]
            points(:, 6) = [SQRT(REAL(l + power + 2, dp)/eta)*(1.0_dp + 1.e-8_dp), 0.0_dp, 0.0_dp]
            c = (eta/alpha)**(REAL(l, dp) + 1.5_dp)
            DO j = 1, power
               c = c*(REAL(l + j, dp) + 0.5_dp)/alpha
            END DO
            DO sample = 1, SIZE(points, 2)
               point = points(:, sample)
               s = SUM(point**2)
               r = SQRT(s)
               CALL gapw_compensated_potential(alpha, power, l, eta, point, value(:n), grad(:, :n), hess(:, :, :n))
               DO m = -l, l
                  angular = 0.0_dp
                  IF (r > 0.0_dp) THEN
                     direction(:, 1) = point/r
                     CALL y_lm(direction, angular, l, m)
                     angular = angular*r**l
                  ELSE IF (l == 0) THEN
                     angular = 1.0_dp/SQRT(fourpi)
                  END IF
                  density = (s**power*EXP(-alpha*s) - c*EXP(-eta*s))*angular(1)
                  reference = -fourpi*density
                  CALL check_error(hess(1, 1, m + l + 1) + hess(2, 2, m + l + 1) + &
                                   hess(3, 3, m + l + 1), reference, 2.e-11_dp)
               END DO
               DO axis = 1, 3
                  ! Fourth-order differences permit a larger step without magnifying
                  ! roundoff where the direct and complementary-tail formulas join.
                  shifted = point
                  shifted(axis) = shifted(axis) + step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vp(:n), gp(:, :n))
                  shifted(axis) = shifted(axis) - 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vm(:n), gm(:, :n))
                  shifted(axis) = point(axis) + 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vpp(:n), gpp(:, :n))
                  shifted(axis) = point(axis) - 2.0_dp*step
                  CALL gapw_compensated_potential(alpha, power, l, eta, shifted, vmm(:n), gmm(:, :n))
                  DO m = 1, n
                     CALL check_error((8.0_dp*(vp(m) - vm(m)) - vpp(m) + vmm(m))/(12.0_dp*step), &
                                      grad(axis, m), 2.e-8_dp)
                     DO j = 1, 3
                        CALL check_error((8.0_dp*(gp(j, m) - gm(j, m)) - gpp(j, m) + gmm(j, m))/(12.0_dp*step), &
                                         hess(j, axis, m), 2.e-8_dp)
                     END DO
                  END DO
               END DO
               ! -d/d(alpha) raises the radial power, including the compensation moment.
               CALL gapw_compensated_potential(alpha + step, power, l, eta, point, vp(:n), gp(:, :n))
               CALL gapw_compensated_potential(alpha - step, power, l, eta, point, vm(:n), gm(:, :n))
               CALL gapw_compensated_potential(alpha + 2.0_dp*step, power, l, eta, point, vpp(:n), gp(:, :n))
               CALL gapw_compensated_potential(alpha - 2.0_dp*step, power, l, eta, point, vmm(:n), gm(:, :n))
               CALL gapw_compensated_potential(alpha, power + 1, l, eta, point, higher(:n), gp(:, :n))
               DO m = 1, n
                  CALL check_error(-(8.0_dp*(vp(m) - vm(m)) - vpp(m) + vmm(m))/(12.0_dp*step), &
                                   higher(m), 2.e-8_dp)
               END DO
            END DO
         END DO
      END DO
      ! Independent monopole solution, including the exponentially small far field.
      DO sample = 1, SIZE(points, 2)
         point = points(:, sample)
         r = SQRT(SUM(point**2))
         IF (r == 0.0_dp) THEN
            reference = SQRT(pi)/alpha*(1.0_dp - SQRT(eta/alpha))
         ELSE
            tail_difference = ERFC(SQRT(eta)*r) - ERFC(SQRT(alpha)*r)
            reference = pi/(2.0_dp*alpha**1.5_dp)*tail_difference/r
         END IF
         CALL gapw_compensated_potential(alpha, 0, 0, eta, point, value(:1), grad(:, :1))
         IF (r > 1.e-6_dp) THEN
            CALL check_error(value(1)/reference, 1.0_dp, 2.e-10_dp)
         ELSE
            ! Avoid cancellation in the independent erfc difference near r=0.
            reference = SQRT(pi)/alpha*(1.0_dp - SQRT(eta/alpha))
            CALL check_error(value(1), reference, 2.e-12_dp)
         END IF
      END DO
      CALL gapw_compensated_potential(alpha, 0, 6, alpha, point, value, grad, hess)
      IF (ANY(value /= 0.0_dp) .OR. ANY(grad /= 0.0_dp) .OR. ANY(hess /= 0.0_dp)) THEN
         ERROR STOP 'Identical compensated Gaussian did not cancel exactly'
      END IF
   END SUBROUTINE check_multipoles

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_augmentation()
      INTEGER, PARAMETER                                 :: nr = 801
      REAL(dp), PARAMETER                                :: a = 1.1_dp, b = 1.3_dp, dr = 0.01_dp, &
                                                            eta = 0.55_dp, step = 2.e-5_dp

      INTEGER                                            :: axis, ir, lg, sample
      REAL(dp) :: derivative(2, 2, 1, 3), dref(3), integral(2, 2, 1), local_error, minus(2, 2, 1), &
         plus(2, 2, 1), point(3), radial_weights(nr), reference, shifted(3)
      TYPE(uspp_augmentation_type)                       :: aug

      aug%nproj = 2
      aug%lmax = 2
      ALLOCATE (aug%r(nr), aug%weighted_q(nr, 3, 0:2), aug%radial_index(2), aug%lbeta(2), aug%gaunt(2, 2, 9))
      aug%radial_index(:) = [1, 2]
      aug%lbeta(:) = [0, 1]
      aug%gaunt = 0.0_dp
      aug%gaunt(1, 1, 1) = 1.0_dp/SQRT(fourpi)
      aug%gaunt(1, 2, 3) = 1.0_dp/SQRT(fourpi)
      aug%gaunt(2, 1, 3) = aug%gaunt(1, 2, 3)
      DO ir = 1, nr
         aug%r(ir) = REAL(ir - 1, dp)*dr
      END DO
      radial_weights = 2.0_dp*dr/3.0_dp
      radial_weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
      radial_weights(1) = dr/3.0_dp
      radial_weights(nr) = dr/3.0_dp
      aug%weighted_q = 0.0_dp
      aug%weighted_q(:, 1, 0) = radial_weights*aug%r**2*EXP(-a*aug%r**2)
      aug%weighted_q(:, 2, 1) = aug%r*aug%weighted_q(:, 1, 0)
      CALL init_lebedev_grids()
      lg = get_number_of_lebedev_grid(n=194)
      DO sample = 1, 3
         point = REAL(sample - 1, dp)*[0.31_dp, -0.52_dp, 0.73_dp]
         CALL gapw_augmentation_integrals(aug, point, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                          fourpi*lebedev_grid(lg)%w, integral, derivative)
         CALL gaussian_reference(point, a, b, eta, reference, dref)
         CALL check_error(integral(1, 1, 1), reference, 2.e-10_dp)
         ! s-p augmentation is sqrt(3)/(2*a) times the center derivative of s-s.
         CALL check_error(integral(1, 2, 1), SQRT(3.0_dp)/(2.0_dp*a)*dref(3), 2.e-10_dp)
         IF (ANY(integral(:, :, 1) /= TRANSPOSE(integral(:, :, 1)))) ERROR STOP 'Asymmetric Coulomb tensor'
         DO axis = 1, 3
            CALL check_error(derivative(1, 1, 1, axis), dref(axis), 2.e-10_dp)
            IF (sample /= 2) CYCLE
            shifted = point
            shifted(axis) = shifted(axis) + step
            CALL gapw_augmentation_integrals(aug, shifted, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                             fourpi*lebedev_grid(lg)%w, plus)
            shifted(axis) = shifted(axis) - 2.0_dp*step
            CALL gapw_augmentation_integrals(aug, shifted, b, 0, 0, eta, lebedev_grid(lg)%r, &
                                             fourpi*lebedev_grid(lg)%w, minus)
            local_error = MAXVAL(ABS((plus(:, :, 1) - minus(:, :, 1))/(2.0_dp*step) - derivative(:, :, 1, axis)))
            CALL check_error(local_error, 0.0_dp, 2.e-9_dp)
         END DO
      END DO
      CALL deallocate_lebedev_grids()
   END SUBROUTINE check_augmentation

! **************************************************************************************************
!> \brief ...
!> \param point ...
!> \param a ...
!> \param b ...
!> \param eta ...
!> \param value ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE gaussian_reference(point, a, b, eta, value, gradient)
      REAL(dp), INTENT(IN)                               :: point(3), a, b, eta
      REAL(dp), INTENT(OUT)                              :: value, gradient(3)

      REAL(dp)                                           :: c, delta, r, radial, u, v

      r = SQRT(SUM(point**2))
      u = SQRT(a*b/(a + b))
      v = SQRT(a*eta/(a + eta))
      c = (pi*pi/(a*b))**1.5_dp/(fourpi*SQRT(fourpi))
      gradient = 0.0_dp
      IF (r == 0.0_dp) THEN
         value = c*2.0_dp/SQRT(pi)*(u - v)
      ELSE
         delta = ERF(u*r) - ERF(v*r)
         value = c*delta/r
         radial = c*(2.0_dp/SQRT(pi)*(u*EXP(-(u*r)**2) - v*EXP(-(v*r)**2))/r - delta/r**2)
         gradient = radial*point/r
      END IF
   END SUBROUTINE gaussian_reference

! **************************************************************************************************
!> \brief Triclinic lattice, wrapped representatives and strain against analytic Gaussian pairs.
! **************************************************************************************************
   SUBROUTINE check_periodic()
      INTEGER, PARAMETER                                 :: nr = 121
      REAL(dp), PARAMETER                                :: a = 2.0_dp, b = 1.3_dp, dr = 0.05_dp, &
                                                            eta = 0.55_dp, step = 2.e-5_dp

      INTEGER                                            :: axis, axis2, boundary, ir, ix, iy, iz, &
                                                            lg, periodicity(3), sample
      REAL(dp) :: boundary_point(3), changed(3, 3), dref(3), fractional(3), h_inv(3, 3), &
         hmat(3, 3), image_point(3), inverse(3, 3), minus(1, 1, 1, 13), plus(1, 1, 1, 13), &
         point(3), radial_weights(nr), reference(13), shifted(1, 1, 1, 13), tensor(1, 1, 1, 13), &
         value
      TYPE(uspp_augmentation_type)                       :: aug, compact

      aug%nproj = 1
      aug%lmax = 0
      ALLOCATE (aug%r(nr), aug%weighted_q(nr, 1, 0:0), aug%radial_index(1), aug%lbeta(1), aug%gaunt(1, 1, 1))
      aug%radial_index = 1
      aug%lbeta = 0
      aug%gaunt = 1.0_dp/SQRT(fourpi)
      DO ir = 1, nr
         aug%r(ir) = REAL(ir - 1, dp)*dr
      END DO
      radial_weights = 2.0_dp*dr/3.0_dp
      radial_weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
      radial_weights([1, nr]) = dr/3.0_dp
      aug%weighted_q(:, 1, 0) = radial_weights*aug%r**2*EXP(-a*aug%r**2)
      hmat = RESHAPE([6.2_dp, 0.0_dp, 0.0_dp, 1.3_dp, 6.7_dp, 0.0_dp, -0.8_dp, 0.6_dp, 7.1_dp], [3, 3])
      CALL triangular_inverse(hmat, h_inv)
      point = [0.31_dp, -0.52_dp, 0.73_dp]
      fractional = MATMUL(h_inv, point)
      reference = 0.0_dp
      DO iz = -4, 4
         DO iy = -4, 4
            DO ix = -4, 4
               image_point = point + MATMUL(hmat, REAL([ix, iy, iz], dp))
               CALL gaussian_reference(image_point, a, b, eta, value, dref)
               reference(1) = reference(1) + value
               reference(2:4) = reference(2:4) + dref
               DO axis = 1, 3
                  DO axis2 = 1, 3
                     ir = 4 + axis + 3*(axis2 - 1)
                     reference(ir) = reference(ir) + dref(axis)*image_point(axis2)
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL init_lebedev_grids()
      lg = get_number_of_lebedev_grid(n=194)
      CALL gapw_periodic_augmentation_integrals(aug, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor)
      DO ir = 1, 13
         CALL check_error(tensor(1, 1, 1, ir), reference(ir), 3.e-9_dp)
      END DO
      CALL gapw_periodic_augmentation_integrals(aug, point + MATMUL(hmat, REAL([2, -3, 1], dp)), b, 0, 0, &
                                                eta, hmat, h_inv, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-16_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 2.e-12_dp)
      ! A diagonal strain and a shear, with fixed fractional centers and rigid radial functions.
      DO sample = 1, 2
         axis = 1
         axis2 = sample
         changed = hmat
         changed(axis, :) = changed(axis, :) + step*hmat(axis2, :)
         CALL triangular_inverse(changed, inverse)
         CALL gapw_periodic_augmentation_integrals(aug, MATMUL(changed, fractional), b, 0, 0, eta, changed, &
                                                   inverse, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, plus)
         changed = hmat
         changed(axis, :) = changed(axis, :) - step*hmat(axis2, :)
         CALL triangular_inverse(changed, inverse)
         CALL gapw_periodic_augmentation_integrals(aug, MATMUL(changed, fractional), b, 0, 0, eta, changed, &
                                                   inverse, lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, minus)
         CALL check_error((plus(1, 1, 1, 1) - minus(1, 1, 1, 1))/(2.0_dp*step), &
                          tensor(1, 1, 1, 4 + axis + 3*(axis2 - 1)), 3.e-9_dp)
      END DO
      ! Isolated, wire and slab boundaries, including nonconsecutive periodic axes.
      ! Move beyond a nonperiodic face to detect accidental minimum-image wrapping.
      DO boundary = 0, 2
         SELECT CASE (boundary)
         CASE (0)
            periodicity = [0, 0, 0]
         CASE (1)
            periodicity = [0, 1, 0]
         CASE (2)
            periodicity = [1, 0, 1]
         END SELECT
         DO axis = 1, 3
            IF (periodicity(axis) == 0) EXIT
         END DO
         boundary_point = point + 0.55_dp*hmat(:, axis)
         reference = 0.0_dp
         DO iz = -4*periodicity(3), 4*periodicity(3)
            DO iy = -4*periodicity(2), 4*periodicity(2)
               DO ix = -4*periodicity(1), 4*periodicity(1)
                  image_point = boundary_point + MATMUL(hmat, REAL([ix, iy, iz], dp))
                  CALL gaussian_reference(image_point, a, b, eta, value, dref)
                  reference(1) = reference(1) + value
                  reference(2:4) = reference(2:4) + dref
                  DO axis = 1, 3
                     DO axis2 = 1, 3
                        ir = 4 + axis + 3*(axis2 - 1)
                        reference(ir) = reference(ir) + dref(axis)*image_point(axis2)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL gapw_periodic_augmentation_integrals(aug, boundary_point, b, 0, 0, eta, hmat, h_inv, &
                                                   lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor, &
                                                   periodicity=periodicity)
         DO ir = 1, 13
            CALL check_error(tensor(1, 1, 1, ir), reference(ir), 3.e-9_dp)
         END DO
         boundary_point = boundary_point + MATMUL(hmat, REAL(periodicity*[2, -3, 1], dp))
         CALL gapw_periodic_augmentation_integrals(aug, boundary_point, b, 0, 0, eta, hmat, h_inv, &
                                                   lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted, &
                                                   periodicity=periodicity)
         CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 2.e-12_dp)
      END DO
      ! A zero radial suffix may be omitted from weighted_q while r retains the
      ! original mesh for spline construction. Check all thirteen tensor entries.
      aug%weighted_q(nr/2 + 1:, :, :) = 0.0_dp
      CALL gapw_periodic_augmentation_integrals(aug, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, tensor)
      compact = aug
      DEALLOCATE (compact%weighted_q)
      ALLOCATE (compact%weighted_q(nr/2, 1, 0:0))
      compact%weighted_q(:, :, :) = aug%weighted_q(:nr/2, :, :)
      CALL gapw_periodic_augmentation_integrals(compact, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted - tensor)), 0.0_dp, 1.e-14_dp)
      DEALLOCATE (compact%weighted_q)
      ALLOCATE (compact%weighted_q(0, 1, 0:0))
      CALL gapw_periodic_augmentation_integrals(compact, point, b, 0, 0, eta, hmat, h_inv, &
                                                lebedev_grid(lg)%r, fourpi*lebedev_grid(lg)%w, 1.e-14_dp, shifted)
      CALL check_error(MAXVAL(ABS(shifted)), 0.0_dp, 0.0_dp)
      CALL deallocate_lebedev_grids()
   END SUBROUTINE check_periodic

! **************************************************************************************************
!> \brief Inverse of the upper-triangular test cell.
!> \param hmat Test cell.
!> \param inverse Inverse cell.
! **************************************************************************************************
   SUBROUTINE triangular_inverse(hmat, inverse)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3)
      REAL(dp), INTENT(OUT)                              :: inverse(3, 3)

      INTEGER                                            :: i

      inverse = 0.0_dp
      DO i = 1, 3
         inverse(i, i) = 1.0_dp/hmat(i, i)
      END DO
      inverse(1, 2) = -hmat(1, 2)/(hmat(1, 1)*hmat(2, 2))
      inverse(2, 3) = -hmat(2, 3)/(hmat(2, 2)*hmat(3, 3))
      inverse(1, 3) = (hmat(1, 2)*hmat(2, 3)/hmat(2, 2) - hmat(1, 3))/(hmat(1, 1)*hmat(3, 3))
   END SUBROUTINE triangular_inverse

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param reference ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_error(value, reference, tolerance)
      REAL(dp), INTENT(IN)                               :: value, reference, tolerance

      REAL(dp)                                           :: error

      error = ABS(value - reference)/MAX(1.0_dp, ABS(reference))
      max_error = MAX(max_error, error)
      IF (.NOT. (error <= tolerance)) THEN
         PRINT *, 'Coulomb check: value, reference, error, tolerance:', value, reference, error, tolerance
         ERROR STOP 'Compensated Coulomb integral check failed'
      END IF
   END SUBROUTINE check_error

END PROGRAM gapw_coulomb_integrals_unittest

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original scalar radial fields in the native real-harmonic projector order.
!>        Regular radial factors are splined in r**2, with Cartesian solid harmonics.
!>        An exact zero suffix supplies a C2 closure; no numerical cutoff is introduced.
! **************************************************************************************************
MODULE upf_radial_fields
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE splines,                         ONLY: spline3ders
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE upf_radial_fields_type
      INTEGER :: nproj = 0, lmax = -1
      LOGICAL :: finite_table = .FALSE.
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:), radial_end(:)
      REAL(dp) :: radius = 0.0_dp
      REAL(dp), ALLOCATABLE :: x(:), coefficients(:, :, :)
   END TYPE upf_radial_fields_type

   PUBLIC :: upf_radial_fields_type, upf_radial_fields_init, upf_radial_fields_evaluate
   PUBLIC :: upf_radial_fields_periodic

CONTAINS

! **************************************************************************************************
!> \brief Cache original f_l(r)Y_lm, without weights, normalization changes or channel rotations.
!>        f_l/r**l is continued from the first positive-radius cubic to the origin.
!>        The final interval into the exact zero suffix matches value and two derivatives.
!>        Noncompact tables require an explicit boundary policy and are rejected here.
!> \param r Original strictly increasing nonnegative radial mesh.
!> \param radial Original f_l(r), including its full exact zero suffix.
!> \param lbeta Angular momentum of each original radial channel.
!> \param fields Immutable cache, initialized outside threaded regions.
!> \param finite_table Retain a noncompact reference table. Queries beyond its last node are invalid.
! **************************************************************************************************
   SUBROUTINE upf_radial_fields_init(r, radial, lbeta, fields, finite_table)
      REAL(dp), INTENT(IN)                               :: r(:), radial(:, :)
      INTEGER, INTENT(IN)                                :: lbeta(:)
      TYPE(upf_radial_fields_type), INTENT(OUT)          :: fields
      LOGICAL, INTENT(IN), OPTIONAL                      :: finite_table

      INTEGER                                            :: first, i, j, l, last, m, n, p
      REAL(dp)                                           :: h
      REAL(dp), ALLOCATABLE                              :: slope(:), value(:), y(:)

      n = SIZE(r)
      CPASSERT(n >= 4 .AND. SIZE(lbeta) > 0)
      CPASSERT(ALL(SHAPE(radial) == [n, SIZE(lbeta)]))
      CPASSERT(ALL(ieee_is_finite(r)) .AND. ALL(ieee_is_finite(radial)))
      CPASSERT(ALL(r >= 0.0_dp) .AND. ALL(r(2:) > r(:n - 1)))
      CPASSERT(ALL(lbeta >= 0) .AND. ALL(lbeta <= 20))
      first = 1
      IF (r(1) == 0.0_dp) first = 2
      n = SIZE(r) - first + 1
      CPASSERT(n >= 4)
      fields%lmax = MAXVAL(lbeta)
      fields%nproj = SUM(2*lbeta + 1)
      IF (PRESENT(finite_table)) fields%finite_table = finite_table
      IF (fields%finite_table) fields%radius = r(SIZE(r))
      CALL init_orbital_pointers(fields%lmax)
      CALL init_spherical_harmonics(fields%lmax, -1)
      ALLOCATE (fields%lbeta, SOURCE=lbeta)
      ALLOCATE (fields%x(n + 1), fields%coefficients(0:5, n, SIZE(lbeta)), &
                fields%radial_end(SIZE(lbeta)), y(n), value(n + 1), slope(n + 1))
      fields%x(1) = 0.0_dp
      fields%x(2:) = r(first:)**2
      fields%coefficients = 0.0_dp
      fields%radial_end = 1
      DO i = 1, SIZE(lbeta)
         IF (ALL(radial(:, i) == 0.0_dp)) CYCLE
         last = SIZE(r)
         IF (fields%finite_table) THEN
            fields%radial_end(i) = SIZE(fields%x)
         ELSE
            DO WHILE (radial(last, i) == 0.0_dp)
               last = last - 1
            END DO
            IF (last == SIZE(r)) THEN
               CPABORT("Original radial field requires an exact zero suffix")
            END IF
            fields%radial_end(i) = last - first + 3
            fields%radius = MAX(fields%radius, r(last + 1))
         END IF
         CPASSERT(last >= first)
         y(:) = radial(first:, i)/r(first:)**lbeta(i)
         CPASSERT(ALL(ieee_is_finite(y)))
         CALL spline3ders(fields%x(2:), y, fields%x, ynew=value, dynew=slope)
         DO j = 1, fields%radial_end(i) - 1
            h = fields%x(j + 1) - fields%x(j)
            ASSOCIATE (c => fields%coefficients(:, j, i))
               c(1) = value(j)
               c(2) = h*slope(j)
               c(3) = 3.0_dp*(value(j + 1) - c(1)) - 2.0_dp*c(2) - h*slope(j + 1)
               c(4) = 2.0_dp*(c(1) - value(j + 1)) + c(2) + h*slope(j + 1)
               IF (j + 1 == fields%radial_end(i) .AND. .NOT. fields%finite_table) THEN
                  c(4) = -10.0_dp*c(1) - 6.0_dp*c(2) - 3.0_dp*c(3)
                  c(5) = 15.0_dp*c(1) + 8.0_dp*c(2) + 3.0_dp*c(3)
                  c(6) = -6.0_dp*c(1) - 3.0_dp*c(2) - c(3)
               END IF
            END ASSOCIATE
         END DO
      END DO
      ALLOCATE (fields%radial_index(fields%nproj), fields%angular_index(fields%nproj))
      p = 0
      DO i = 1, SIZE(lbeta)
         l = lbeta(i)
         DO m = 1, 2*l + 1
            p = p + 1
            fields%radial_index(p) = i
            fields%angular_index(p) = l*l + m
         END DO
      END DO
   END SUBROUTINE upf_radial_fields_init

! **************************************************************************************************
!> \brief Evaluate values and Cartesian derivatives from the same radial/harmonic polynomial.
!> \param fields Initialized original radial fields.
!> \param displacement Point minus field center, in bohr.
!> \param value Values in the native magnetic-projector order.
!> \param gradient Optional Cartesian point gradients; center gradients have opposite sign.
!> \param hessian Optional Cartesian Hessians, including origin and coordinate axes.
! **************************************************************************************************
   SUBROUTINE upf_radial_fields_evaluate(fields, displacement, value, gradient, hessian)
      TYPE(upf_radial_fields_type), INTENT(IN)           :: fields
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: value(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :), hessian(:, :, :)

      INTEGER                                            :: a, axis, axis2, hi, i, l, lo, m, mid, p, &
                                                            powers(3), reduced(3)
      REAL(dp)                                           :: angular, da(3), dda(3, 3), factor, h, &
                                                            radial, rd, rdd, s, t, term

      CPASSERT(ALLOCATED(fields%x) .AND. SIZE(value) == fields%nproj)
      value = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [3, fields%nproj]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(hessian)) THEN
         CPASSERT(ALL(SHAPE(hessian) == [3, 3, fields%nproj]))
         hessian = 0.0_dp
      END IF
      s = SUM(displacement**2)
      IF (fields%finite_table) THEN
         CPASSERT(s <= fields%radius**2)
      ELSE IF (s >= fields%radius**2) THEN
         RETURN
      END IF
      lo = 1
      hi = SIZE(fields%x)
      DO WHILE (hi - lo > 1)
         mid = (lo + hi)/2
         IF (s < fields%x(mid)) THEN
            hi = mid
         ELSE
            lo = mid
         END IF
      END DO
      h = fields%x(hi) - fields%x(lo)
      t = (s - fields%x(lo))/h
      DO p = 1, fields%nproj
         i = fields%radial_index(p)
         IF (lo >= fields%radial_end(i)) CYCLE
         l = fields%lbeta(i)
         m = fields%angular_index(p) - l*l
         ASSOCIATE (c => fields%coefficients(:, lo, i))
            radial = c(1) + t*(c(2) + t*(c(3) + t*(c(4) + t*(c(5) + t*c(6)))))
            rd = (c(2) + t*(2.0_dp*c(3) + t*(3.0_dp*c(4) + t*(4.0_dp*c(5) + t*5.0_dp*c(6)))))/h
            rdd = (2.0_dp*c(3) + t*(6.0_dp*c(4) + t*(12.0_dp*c(5) + t*20.0_dp*c(6))))/h**2
         END ASSOCIATE
         angular = 0.0_dp
         da = 0.0_dp
         dda = 0.0_dp
         DO a = 1, nco(l)
            powers = indco(:, ncoset(l - 1) + a)
            factor = orbtramat(l)%slm(m, a)
            angular = angular + factor*PRODUCT(displacement**powers)
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               term = factor*REAL(powers(axis), dp)
               da(axis) = da(axis) + term*PRODUCT(displacement**reduced)
               DO axis2 = 1, 3
                  IF (reduced(axis2) == 0) CYCLE
                  reduced(axis2) = reduced(axis2) - 1
                  dda(axis, axis2) = dda(axis, axis2) + &
                                     term*REAL(reduced(axis2) + 1, dp)*PRODUCT(displacement**reduced)
                  reduced(axis2) = reduced(axis2) + 1
               END DO
            END DO
         END DO
         value(p) = radial*angular
         IF (PRESENT(gradient)) gradient(:, p) = 2.0_dp*displacement*rd*angular + radial*da
         IF (PRESENT(hessian)) THEN
            DO axis2 = 1, 3
               DO axis = 1, 3
                  hessian(axis, axis2, p) = 4.0_dp*displacement(axis)*displacement(axis2)*rdd*angular + &
                                            2.0_dp*rd*(displacement(axis)*da(axis2) + displacement(axis2)*da(axis)) + &
                                            radial*dda(axis, axis2)
               END DO
               hessian(axis2, axis2, p) = hessian(axis2, axis2, p) + 2.0_dp*rd*angular
            END DO
         END IF
      END DO
   END SUBROUTINE upf_radial_fields_evaluate

! **************************************************************************************************
!> \brief Sum supplied Bloch images on common Cartesian quadrature points, including responses.
!>        Caller supplies complete image bounds; no minimum-image approximation is made.
!>        Cell responses use Cartesian cell columns, with reduced k fixed unless differentiated.
!> \param fields Original radial cache.
!> \param points Cartesian quadrature points, one per column.
!> \param center Cartesian original center.
!> \param cell Cartesian cell vectors, one per column.
!> \param images Integer image vectors, one per column.
!> \param reduced_k Reduced Bloch wavevector, in reciprocal lattice units.
!> \param value Bloch values, point-major and native projector-minor.
!> \param gradient Optional spatial gradients, Cartesian axis last.
!> \param response Optional geometry/Bloch response, including image phase derivatives.
!> \param gradient_response Optional joint spatial/geometry response, requiring Hessians.
!> \param points_response Optional Cartesian point displacements.
!> \param center_response Optional Cartesian center displacement.
!> \param cell_response Optional cell-vector response.
!> \param reduced_k_response Optional reduced Bloch-wavevector response.
! **************************************************************************************************
   SUBROUTINE upf_radial_fields_periodic(fields, points, center, cell, images, reduced_k, value, &
                                         gradient, response, gradient_response, points_response, &
                                         center_response, cell_response, reduced_k_response)
      TYPE(upf_radial_fields_type), INTENT(IN)           :: fields
      REAL(dp), INTENT(IN)                               :: points(:, :), center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: value(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: points_response(:, :), &
                                                            center_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp)                                        :: dphase, phase
      INTEGER                                            :: axis, image, point
      REAL(dp)                                           :: displacement(3), dr(3), dshift(3), &
                                                            shift(3), v(fields%nproj), &
                                                            hh(3, 3, fields%nproj), &
                                                            g(3, fields%nproj)

      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(value) == [SIZE(points, 2), fields%nproj]))
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [SIZE(points, 2), fields%nproj, 3]))
         gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(value)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(gradient_response)) THEN
         IF (.NOT. (ALL(SHAPE(gradient_response) == [SIZE(points, 2), fields%nproj, 3]))) THEN
            CPABORT("Invalid radial field response dimensions")
         END IF
         gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(points_response)) THEN
         CPASSERT(ALL(SHAPE(points_response) == SHAPE(points)))
      END IF
      DO image = 1, SIZE(images, 2)
         shift = center + MATMUL(cell, REAL(images(:, image), dp))
         dshift = 0.0_dp
         IF (PRESENT(center_response)) dshift = center_response
         IF (PRESENT(cell_response)) dshift = dshift + MATMUL(cell_response, REAL(images(:, image), dp))
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, image), dp)), dp))
         dphase = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (PRESENT(reduced_k_response)) THEN
            dphase = phase*CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k_response, REAL(images(:, image), dp)), dp)
         END IF
         DO point = 1, SIZE(points, 2)
            displacement = points(:, point) - shift
            dr = -dshift
            IF (PRESENT(points_response)) dr = dr + points_response(:, point)
            IF (SUM(displacement**2) >= fields%radius**2) CYCLE
            IF (PRESENT(gradient_response)) THEN
               CALL upf_radial_fields_evaluate(fields, displacement, v, g, hh)
            ELSE IF (PRESENT(gradient) .OR. PRESENT(response)) THEN
               CALL upf_radial_fields_evaluate(fields, displacement, v, g)
            ELSE
               CALL upf_radial_fields_evaluate(fields, displacement, v)
            END IF
            value(point, :) = value(point, :) + phase*v
            IF (PRESENT(response)) response(point, :) = response(point, :) + &
                                                        phase*MATMUL(dr, g) + dphase*v
            DO axis = 1, 3
               IF (PRESENT(gradient)) gradient(point, :, axis) = gradient(point, :, axis) + phase*g(axis, :)
               IF (PRESENT(gradient_response)) THEN
                  gradient_response(point, :, axis) = gradient_response(point, :, axis) + &
                                                      phase*MATMUL(dr, hh(axis, :, :)) + dphase*g(axis, :)
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE upf_radial_fields_periodic

END MODULE upf_radial_fields

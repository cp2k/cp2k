!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native contracted Gaussian fields and their consistent Cartesian/Bloch responses.
!>        Uses the stored sphi coefficients and native spherical AO order without renormalization.
! **************************************************************************************************
MODULE gaussian_basis_fields
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gaussian_basis_fields_evaluate, gaussian_basis_fields_periodic

CONTAINS

! **************************************************************************************************
!> \brief Evaluate a native contracted basis at explicit point-minus-center displacements.
!>        No minimum-image operation is applied; an explicit cutoff can skip excluded points.
!> \param basis Initialized native basis, including its original normalization in sphi.
!> \param displacements Cartesian displacements, indexed (axis,point).
!> \param value New values, indexed (point,AO).
!> \param gradient Optional Cartesian derivatives, indexed (point,AO,axis).
!> \param hessian Optional Cartesian Hessians, indexed (point,AO,axis,axis).
!> \param cutoff_squared Optional squared cutoff radius; fields beyond it are zero.
! **************************************************************************************************
   SUBROUTINE gaussian_basis_fields_evaluate(basis, displacements, value, gradient, hessian, cutoff_squared)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: displacements(:, :)
      REAL(dp), INTENT(OUT)                              :: value(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :, :), hessian(:, :, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: cutoff_squared

      INTEGER                                            :: axis, axis2, first, ico, ipgf, ipoint, &
                                                            irow, iset, isgf, ishell, l, last, &
                                                            ncart, powers(3), reduced(3)
      REAL(dp) :: alpha, exponential, polynomial, polynomial_gradient(3), &
         polynomial_hessian(3, 3), primitive_gradient(3), primitive_hessian(3, 3), r(3), &
         r_squared, weight, coordinate_power(3, 0:MAXVAL(basis%lmax)), &
         cartesian_value(ncoset(MAXVAL(basis%lmax))), &
         cartesian_hessian(3, 3, ncoset(MAXVAL(basis%lmax))), &
         cartesian_gradient(3, ncoset(MAXVAL(basis%lmax)))

      CPASSERT(basis%nset > 0 .AND. basis%nsgf > 0)
      CPASSERT(ASSOCIATED(basis%sphi) .AND. ASSOCIATED(basis%zet))
      CPASSERT(UBOUND(ncoset, 1) >= MAXVAL(basis%lmax))
      CPASSERT(SIZE(displacements, 1) == 3)
      CPASSERT(ALL(SHAPE(value) == [SIZE(displacements, 2), basis%nsgf]))
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [SIZE(value, 1), SIZE(value, 2), 3]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(hessian)) THEN
         CPASSERT(ALL(SHAPE(hessian) == [SIZE(value, 1), SIZE(value, 2), 3, 3]))
         hessian = 0.0_dp
      END IF
      value = 0.0_dp
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(basis, displacements, value, gradient, hessian, cutoff_squared, ncoset, indco) &
      !$OMP PRIVATE(ipoint, r, iset, ncart, ishell, l, first, last, ipgf, alpha, exponential, &
      !$OMP         isgf, ico, irow, weight, powers, reduced, polynomial, polynomial_gradient, &
      !$OMP         polynomial_hessian, primitive_gradient, primitive_hessian, cartesian_value, &
      !$OMP         cartesian_gradient, cartesian_hessian, coordinate_power, r_squared, axis, axis2)
      DO ipoint = 1, SIZE(displacements, 2)
         r = displacements(:, ipoint)
         r_squared = SUM(r**2)
         IF (PRESENT(cutoff_squared)) THEN
            IF (r_squared > cutoff_squared) CYCLE
         END IF
         coordinate_power(:, 0) = 1.0_dp
         DO l = 1, UBOUND(coordinate_power, 2)
            coordinate_power(:, l) = coordinate_power(:, l - 1)*r
         END DO
         ! Cartesian monomials and their derivatives do not depend on the primitive or contraction.
         DO ico = 1, SIZE(cartesian_value)
            powers = indco(:, ico)
            cartesian_value(ico) = coordinate_power(1, powers(1))*coordinate_power(2, powers(2))* &
                                   coordinate_power(3, powers(3))
            IF (.NOT. PRESENT(gradient) .AND. .NOT. PRESENT(hessian)) CYCLE
            cartesian_gradient(:, ico) = 0.0_dp
            cartesian_hessian(:, :, ico) = 0.0_dp
            DO axis = 1, 3
               IF (powers(axis) == 0) CYCLE
               reduced = powers
               reduced(axis) = reduced(axis) - 1
               cartesian_gradient(axis, ico) = REAL(powers(axis), dp)*coordinate_power(1, reduced(1))* &
                                               coordinate_power(2, reduced(2))*coordinate_power(3, reduced(3))
               IF (.NOT. PRESENT(hessian)) CYCLE
               DO axis2 = 1, 3
                  IF (reduced(axis2) == 0) CYCLE
                  reduced(axis2) = reduced(axis2) - 1
                  cartesian_hessian(axis, axis2, ico) = REAL(powers(axis)*(reduced(axis2) + 1), dp)* &
                                                        coordinate_power(1, reduced(1))*coordinate_power(2, reduced(2))* &
                                                        coordinate_power(3, reduced(3))
                  reduced(axis2) = reduced(axis2) + 1
               END DO
            END DO
         END DO
         DO iset = 1, basis%nset
            ncart = ncoset(basis%lmax(iset))
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               exponential = EXP(-alpha*r_squared)
               DO ishell = 1, basis%nshell(iset)
                  l = basis%l(ishell, iset)
                  first = basis%first_sgf(ishell, iset)
                  last = basis%last_sgf(ishell, iset)
                  ! Primitive Cartesian derivatives are shared by every spherical contraction.
                  DO ico = ncoset(l - 1) + 1, ncoset(l)
                     irow = (ipgf - 1)*ncart + ico
                     polynomial = cartesian_value(ico)
                     IF (PRESENT(gradient) .OR. PRESENT(hessian)) THEN
                        polynomial_gradient = cartesian_gradient(:, ico)
                        polynomial_hessian = cartesian_hessian(:, :, ico)
                        DO axis = 1, 3
                           primitive_gradient(axis) = polynomial_gradient(axis) - 2.0_dp*alpha*r(axis)*polynomial
                           IF (.NOT. PRESENT(hessian)) CYCLE
                           DO axis2 = 1, 3
                              primitive_hessian(axis, axis2) = polynomial_hessian(axis, axis2) - 2.0_dp*alpha*( &
                                                        r(axis)*polynomial_gradient(axis2) + r(axis2)*polynomial_gradient(axis)) + &
                                                               4.0_dp*alpha**2*r(axis)*r(axis2)*polynomial
                              IF (axis == axis2) THEN
                                 primitive_hessian(axis, axis2) = primitive_hessian(axis, axis2) - 2.0_dp*alpha*polynomial
                              END IF
                           END DO
                        END DO
                     END IF
                     DO isgf = first, last
                        weight = basis%sphi(irow, isgf)*exponential
                        value(ipoint, isgf) = value(ipoint, isgf) + weight*polynomial
                        IF (.NOT. PRESENT(gradient) .AND. .NOT. PRESENT(hessian)) CYCLE
                        DO axis = 1, 3
                           IF (PRESENT(gradient)) THEN
                              gradient(ipoint, isgf, axis) = gradient(ipoint, isgf, axis) + weight*primitive_gradient(axis)
                           END IF
                           IF (.NOT. PRESENT(hessian)) CYCLE
                           DO axis2 = 1, 3
                              hessian(ipoint, isgf, axis, axis2) = hessian(ipoint, isgf, axis, axis2) + &
                                                                   weight*primitive_hessian(axis, axis2)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE gaussian_basis_fields_evaluate

! **************************************************************************************************
!> \brief Sum explicit Gaussian images with physical Bloch gradients and joint geometry response.
!>        Caller owns the finite image set and its convergence; the set is held fixed in a response.
!> \param basis ...
!> \param points Cartesian points, indexed (axis,point).
!> \param center Physical basis center, without implicit folding.
!> \param cell Cartesian lattice vectors, one vector per column.
!> \param images Explicit integer lattice images, indexed (axis,image).
!> \param reduced_k Reduced Bloch vector; phase is exp(i*2*pi*k.image).
!> \param value New Bloch values, indexed (point,AO).
!> \param gradient Optional physical Cartesian Bloch gradients, axis last.
!> \param response Optional joint first response of the values.
!> \param gradient_response Optional response of the physical spatial gradients.
!> \param point_response Optional Cartesian point response, default zero.
!> \param center_response Optional Cartesian center response, default zero.
!> \param cell_response Optional lattice response, default zero.
!> \param reduced_k_response Optional reduced Bloch-vector response, default zero.
!> \param cutoff_squared Optional squared spatial radius, fixed during geometry responses.
! **************************************************************************************************
   SUBROUTINE gaussian_basis_fields_periodic(basis, points, center, cell, images, reduced_k, value, gradient, &
                                             response, gradient_response, point_response, center_response, &
                                             cell_response, reduced_k_response, cutoff_squared)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(dp), INTENT(IN)                               :: points(:, :), center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: value(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: gradient(:, :, :), response(:, :), &
                                                            gradient_response(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_response(:, :), &
                                                            center_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3), cutoff_squared

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: axis, axis2, image, npoint
      LOGICAL                                            :: need_gradient, need_hessian
      REAL(dp)                                           :: dc(3, 3), dcenter(3), dk(3), lattice(3)
      REAL(dp), ALLOCATABLE                              :: ddisplacement(:, :), displacement(:, :), &
                                                            field(:, :), g(:, :, :), h(:, :, :, :)

      npoint = SIZE(points, 2)
      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(value) == [npoint, basis%nsgf]))
      need_gradient = PRESENT(gradient) .OR. PRESENT(response) .OR. PRESENT(gradient_response)
      need_hessian = PRESENT(gradient_response)
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [npoint, basis%nsgf, 3]))
         gradient = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(value)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (PRESENT(gradient_response)) THEN
         CPASSERT(ALL(SHAPE(gradient_response) == [npoint, basis%nsgf, 3]))
         gradient_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      dc = 0.0_dp
      dcenter = 0.0_dp
      dk = 0.0_dp
      IF (PRESENT(cell_response)) dc = cell_response
      IF (PRESENT(center_response)) dcenter = center_response
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      ALLOCATE (displacement(3, npoint), ddisplacement(3, npoint), field(npoint, basis%nsgf))
      IF (need_gradient) ALLOCATE (g(npoint, basis%nsgf, 3))
      IF (need_hessian) ALLOCATE (h(npoint, basis%nsgf, 3, 3))
      value = CMPLX(0.0_dp, 0.0_dp, dp)
      DO image = 1, SIZE(images, 2)
         lattice = REAL(images(:, image), dp)
         displacement(:, :) = points - SPREAD(center + MATMUL(cell, lattice), 2, npoint)
         ddisplacement(:, :) = -SPREAD(dcenter + MATMUL(dc, lattice), 2, npoint)
         IF (PRESENT(point_response)) THEN
            CPASSERT(ALL(SHAPE(point_response) == SHAPE(points)))
            ddisplacement(:, :) = ddisplacement + point_response
         END IF
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, lattice), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, lattice), dp)*phase
         IF (need_hessian) THEN
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, g, h, cutoff_squared)
         ELSE IF (need_gradient) THEN
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, g, cutoff_squared=cutoff_squared)
         ELSE
            CALL gaussian_basis_fields_evaluate(basis, displacement, field, cutoff_squared=cutoff_squared)
         END IF
         value = value + phase*field
         IF (PRESENT(gradient)) gradient = gradient + phase*g
         IF (PRESENT(response)) THEN
            response = response + phase_response*field
            DO axis = 1, 3
               response = response + phase*SPREAD(ddisplacement(axis, :), 2, basis%nsgf)*g(:, :, axis)
            END DO
         END IF
         IF (PRESENT(gradient_response)) THEN
            gradient_response = gradient_response + phase_response*g
            DO axis = 1, 3
               DO axis2 = 1, 3
                  gradient_response(:, :, axis) = gradient_response(:, :, axis) + &
                                                  phase*SPREAD(ddisplacement(axis2, :), 2, basis%nsgf)*h(:, :, axis, axis2)
               END DO
            END DO
         END IF
      END DO
   END SUBROUTINE gaussian_basis_fields_periodic

END MODULE gaussian_basis_fields

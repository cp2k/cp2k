!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common point quadratures and semilocal adjoints for inverse UPF basis blocks.
!>        Point rows may be MPI-distributed; each supplied column block is replicated.
!>        Original projector augmentation remains a separate, explicitly supplied operator.
! **************************************************************************************************
MODULE upf_inverse_quadrature
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: upf_inverse_pair, upf_inverse_pair_response
   PUBLIC :: upf_inverse_semilocal_pair, upf_inverse_semilocal_response, upf_inverse_smooth_density

CONTAINS

! **************************************************************************************************
!> \brief Integrate left^dagger * weights * right, without discarding intercenter terms.
!>        Use the same rule for beta/Delta, beta/Gaussian and all operator blocks.
!> \param left Local point rows and selected left columns.
!> \param right Matching local point rows and selected right columns.
!> \param weights Original quadrature measure for the local points.
!> \param integral New rectangular block, globally reduced if comm is supplied.
!> \param comm Optional communicator partitioning point rows; empty ranks participate.
! **************************************************************************************************
   SUBROUTINE upf_inverse_pair(left, right, weights, integral, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      CALL check_pair(left, right, weights, SHAPE(integral))
      integral = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL add_pair(left, right, weights, integral)
      IF (PRESENT(comm)) CALL comm%sum(integral)
   END SUBROUTINE upf_inverse_pair

! **************************************************************************************************
!> \brief Differentiate both field legs and the original quadrature measure.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param left_response Response of local left fields.
!> \param right_response Response of local right fields.
!> \param weights_response Response of the same local measure, including strain Jacobians.
!> \param response New rectangular response block.
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_pair_response(left, right, weights, left_response, right_response, &
                                        weights_response, response, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(IN)                            :: left_response(:, :), right_response(:, :)
      REAL(dp), INTENT(IN)                               :: weights_response(:)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      CALL check_pair(left, right, weights, SHAPE(response))
      CPASSERT(ALL(SHAPE(left_response) == SHAPE(left)))
      CPASSERT(ALL(SHAPE(right_response) == SHAPE(right)))
      CPASSERT(SIZE(weights_response) == SIZE(weights))
      response = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL add_pair(left_response, right, weights, response)
      CALL add_pair(left, right_response, weights, response)
      CALL add_pair(left, right, weights_response, response)
      IF (PRESENT(comm)) CALL comm%sum(response)
   END SUBROUTINE upf_inverse_pair_response

! **************************************************************************************************
!> \brief Semilocal rectangular block for density, Cartesian density gradient and tau potentials.
!>        Integral of vn L*R + vg.(grad(L*)R+L*grad(R)) + vt grad(L*).grad(R)/2.
!>        vgradient is dE/d(grad n), not dE/dsigma. Cartesian gradients are physical Bloch
!>        gradients; supplying gradients of periodic factors alone would omit k terms.
!>        Apply to GG/GD/DG/DD blocks and retain the original augmentation pullback separately.
!> \param left Local original Gaussian or difference-wave values.
!> \param right Matching right values.
!> \param left_gradient Cartesian spatial derivatives, axis last.
!> \param right_gradient ...
!> \param weights Original local quadrature measure.
!> \param vdensity Potential conjugate to smooth density.
!> \param vgradient Cartesian potential conjugate to its gradient, axis last.
!> \param vtau Potential conjugate to tau=sum(f*|grad psi|**2)/2.
!> \param integral New rectangular operator block.
!> \param comm Optional point-row communicator; reduction occurs once per block.
! **************************************************************************************************
   SUBROUTINE upf_inverse_semilocal_pair(left, right, left_gradient, right_gradient, weights, &
                                         vdensity, vgradient, vtau, integral, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :), &
                                                            left_gradient(:, :, :), &
                                                            right_gradient(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights(:), vdensity(:), &
                                                            vgradient(:, :), vtau(:)
      COMPLEX(dp), INTENT(OUT)                           :: integral(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)
      INTEGER                                            :: axis, nright
      LOGICAL                                            :: gradient_active(3), tau_active

      CALL check_semilocal(left, right, left_gradient, right_gradient, weights, &
                           vdensity, vgradient, vtau, SHAPE(integral))
      integral = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (.NOT. ANY(weights /= 0.0_dp) .OR. &
          .NOT. (ANY(vdensity /= 0.0_dp) .OR. ANY(vgradient /= 0.0_dp) .OR. ANY(vtau /= 0.0_dp))) THEN
         IF (PRESENT(comm)) CALL comm%sum(integral)
         RETURN
      END IF
      DO axis = 1, 3
         gradient_active(axis) = ANY(vgradient(:, axis) /= 0.0_dp)
      END DO
      tau_active = ANY(vtau /= 0.0_dp)
      nright = SIZE(right, 2)
      ALLOCATE (work(SIZE(right, 1), nright))
      work(:, :) = SPREAD(vdensity, 2, nright)*right
      DO axis = 1, 3
         IF (.NOT. gradient_active(axis)) CYCLE
         work(:, :) = work + SPREAD(vgradient(:, axis), 2, nright)*right_gradient(:, :, axis)
      END DO
      CALL add_pair(left, work, weights, integral)
      DO axis = 1, 3
         IF (.NOT. gradient_active(axis) .AND. .NOT. tau_active) CYCLE
         work(:, :) = SPREAD(vgradient(:, axis), 2, nright)*right
         IF (tau_active) work(:, :) = work + 0.5_dp*SPREAD(vtau, 2, nright)*right_gradient(:, :, axis)
         CALL add_pair(left_gradient(:, :, axis), work, weights, integral)
      END DO
      IF (PRESENT(comm)) CALL comm%sum(integral)
   END SUBROUTINE upf_inverse_semilocal_pair

! **************************************************************************************************
!> \brief Response of the same semilocal block, including moving gradients, potentials and measure.
!>        Geometry responses of spatial gradients must include the inverse and field Hessians.
!> \param left ...
!> \param right ...
!> \param left_gradient ...
!> \param right_gradient ...
!> \param weights ...
!> \param vdensity ...
!> \param vgradient ...
!> \param vtau ...
!> \param left_response ...
!> \param right_response ...
!> \param left_gradient_response ...
!> \param right_gradient_response ...
!> \param weights_response ...
!> \param vdensity_response Zero for a frozen conjugate potential.
!> \param vgradient_response ...
!> \param vtau_response ...
!> \param response New rectangular response block.
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_semilocal_response(left, right, left_gradient, right_gradient, weights, &
                                             vdensity, vgradient, vtau, left_response, right_response, &
                                             left_gradient_response, right_gradient_response, &
                                             weights_response, vdensity_response, vgradient_response, &
                                             vtau_response, response, comm)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :), &
                                                            left_gradient(:, :, :), &
                                                            right_gradient(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights(:), vdensity(:), &
                                                            vgradient(:, :), vtau(:)
      COMPLEX(dp), INTENT(IN) :: left_response(:, :), right_response(:, :), &
         left_gradient_response(:, :, :), right_gradient_response(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights_response(:), &
                                                            vdensity_response(:), &
                                                            vgradient_response(:, :), &
                                                            vtau_response(:)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)
      TYPE(mp_comm_type), INTENT(IN), OPTIONAL           :: comm

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)

      CALL check_semilocal(left, right, left_gradient, right_gradient, weights, &
                           vdensity, vgradient, vtau, SHAPE(response))
      CALL check_semilocal(left_response, right_response, left_gradient_response, right_gradient_response, &
                           weights_response, vdensity_response, vgradient_response, vtau_response, SHAPE(response))
      CPASSERT(ALL(SHAPE(left_response) == SHAPE(left)))
      CPASSERT(ALL(SHAPE(right_response) == SHAPE(right)))
      ALLOCATE (work(SIZE(response, 1), SIZE(response, 2)))
      ! Linearity in either field leg and in each potential/measure gives the complete response.
      CALL upf_inverse_semilocal_pair(left_response, right, left_gradient_response, right_gradient, &
                                      weights, vdensity, vgradient, vtau, response)
      CALL upf_inverse_semilocal_pair(left, right_response, left_gradient, right_gradient_response, &
                                      weights, vdensity, vgradient, vtau, work)
      response = response + work
      CALL upf_inverse_semilocal_pair(left, right, left_gradient, right_gradient, &
                                      weights_response, vdensity, vgradient, vtau, work)
      response = response + work
      CALL upf_inverse_semilocal_pair(left, right, left_gradient, right_gradient, &
                                      weights, vdensity_response, vgradient_response, vtau_response, work)
      response = response + work
      IF (PRESENT(comm)) CALL comm%sum(response)
   END SUBROUTINE upf_inverse_semilocal_response

! **************************************************************************************************
!> \brief Smooth density descriptors whose adjoint is the semilocal pair operator.
!>        Per-spin AO density includes occupations and k weights. Values/gradients are U*chi.
!>        Caller must add original charge, gradient and independent kinetic augmentation.
!>        A charge kernel does not determine a missing kinetic augmentation kernel.
!> \param soft Local inverse basis values.
!> \param soft_gradient Physical Cartesian Bloch gradients, axis last.
!> \param ao_density Hermitian density matrix in reconstructed Gaussian variables.
!> \param density Smooth local point density.
!> \param gradient Smooth Cartesian density gradient, axis last.
!> \param tau Smooth kinetic-energy density, with factor one half.
! **************************************************************************************************
   SUBROUTINE upf_inverse_smooth_density(soft, soft_gradient, ao_density, density, gradient, tau)
      COMPLEX(dp), INTENT(IN)                            :: soft(:, :), soft_gradient(:, :, :), &
                                                            ao_density(:, :)
      REAL(dp), INTENT(OUT)                              :: density(:), gradient(:, :), tau(:)

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :), work(:, :)
      INTEGER                                            :: axis, n, npoint
      REAL(dp)                                           :: tolerance

      n = SIZE(soft, 2)
      npoint = SIZE(soft, 1)
      CPASSERT(ALL(SHAPE(ao_density) == [n, n]))
      CPASSERT(ALL(SHAPE(soft_gradient) == [npoint, n, 3]))
      CPASSERT(SIZE(density) == npoint .AND. SIZE(tau) == npoint)
      CPASSERT(ALL(SHAPE(gradient) == [npoint, 3]))
      IF (n > 0) THEN
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(ao_density)))
         IF (.NOT. (MAXVAL(ABS(ao_density - CONJG(TRANSPOSE(ao_density)))) <= tolerance)) THEN
            CPABORT("Inconsistent inverse quadrature dimensions")
         END IF
      END IF
      ALLOCATE (work(npoint, n), derivative(npoint, n))
      work(:, :) = MATMUL(soft, ao_density)
      density = REAL(SUM(CONJG(soft)*work, DIM=2), dp)
      tau = 0.0_dp
      DO axis = 1, 3
         gradient(:, axis) = 2.0_dp*REAL(SUM(CONJG(soft_gradient(:, :, axis))*work, DIM=2), dp)
         derivative(:, :) = MATMUL(soft_gradient(:, :, axis), ao_density)
         tau = tau + 0.5_dp*REAL(SUM(CONJG(soft_gradient(:, :, axis))*derivative, DIM=2), dp)
      END DO
   END SUBROUTINE upf_inverse_smooth_density

! **************************************************************************************************
!> \brief Accumulate a local rectangular contraction, allowing zero local points/columns.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param integral ...
! **************************************************************************************************
   SUBROUTINE add_pair(left, right, weights, integral)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      COMPLEX(dp), INTENT(INOUT)                         :: integral(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)

      IF (SIZE(weights) == 0 .OR. SIZE(left, 2) == 0 .OR. SIZE(right, 2) == 0) RETURN
      IF (.NOT. ANY(weights /= 0.0_dp)) RETURN
      ALLOCATE (work(SIZE(right, 1), SIZE(right, 2)))
      work(:, :) = SPREAD(weights, 2, SIZE(right, 2))*right
      integral = integral + MATMUL(CONJG(TRANSPOSE(left)), work)
   END SUBROUTINE add_pair

! **************************************************************************************************
!> \brief Validate local point dimensions and finite original measure.
!> \param left ...
!> \param right ...
!> \param weights ...
!> \param output_shape ...
! **************************************************************************************************
   SUBROUTINE check_pair(left, right, weights, output_shape)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :)
      REAL(dp), INTENT(IN)                               :: weights(:)
      INTEGER, INTENT(IN)                                :: output_shape(2)

      IF (.NOT. (SIZE(left, 1) == SIZE(weights) .AND. SIZE(right, 1) == SIZE(weights))) THEN
         CPABORT("Inconsistent inverse quadrature dimensions")
      END IF
      CPASSERT(ALL(output_shape == [SIZE(left, 2), SIZE(right, 2)]))
      CPASSERT(ALL(ieee_is_finite(weights)))
   END SUBROUTINE check_pair

! **************************************************************************************************
!> \brief ...
!> \param left ...
!> \param right ...
!> \param left_gradient ...
!> \param right_gradient ...
!> \param weights ...
!> \param vdensity ...
!> \param vgradient ...
!> \param vtau ...
!> \param output_shape ...
! **************************************************************************************************
   SUBROUTINE check_semilocal(left, right, left_gradient, right_gradient, weights, vdensity, vgradient, vtau, output_shape)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :), right(:, :), &
                                                            left_gradient(:, :, :), &
                                                            right_gradient(:, :, :)
      REAL(dp), INTENT(IN)                               :: weights(:), vdensity(:), &
                                                            vgradient(:, :), vtau(:)
      INTEGER, INTENT(IN)                                :: output_shape(2)

      CALL check_pair(left, right, weights, output_shape)
      IF (.NOT. (ALL(SHAPE(left_gradient) == [SIZE(left, 1), SIZE(left, 2), 3]))) THEN
         CPABORT("Inconsistent inverse quadrature dimensions")
      END IF
      IF (.NOT. (ALL(SHAPE(right_gradient) == [SIZE(right, 1), SIZE(right, 2), 3]))) THEN
         CPABORT("Inconsistent inverse quadrature dimensions")
      END IF
      IF (.NOT. (SIZE(vdensity) == SIZE(weights) .AND. SIZE(vtau) == SIZE(weights))) THEN
         CPABORT("Inconsistent inverse quadrature dimensions")
      END IF
      CPASSERT(ALL(SHAPE(vgradient) == [SIZE(weights), 3]))
   END SUBROUTINE check_semilocal

END MODULE upf_inverse_quadrature

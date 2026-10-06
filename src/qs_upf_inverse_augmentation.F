!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Original scalar USPP charge augmentation in the inverse Quickstep projector frame.
!>        Charge kernels are periodic, without a Bloch phase; W carries the k dependence.
!>        Kinetic augmentation is independent and is not inferred from the charge kernels.
! **************************************************************************************************
MODULE qs_upf_inverse_augmentation
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_derive,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate,&
                                              uspp_augmentation_collocate_atom,&
                                              uspp_augmentation_init,&
                                              uspp_augmentation_integrate,&
                                              uspp_augmentation_integrate_atom,&
                                              uspp_augmentation_prepare_fourier,&
                                              uspp_augmentation_prepare_realspace,&
                                              uspp_augmentation_realspace
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_augmentation_type, qs_upf_inverse_augmentation_prepare
   PUBLIC :: qs_upf_inverse_augmentation_collocate, qs_upf_inverse_augmentation_integrate
   PUBLIC :: qs_upf_inverse_augmentation_response, qs_upf_inverse_augmentation_release
   PUBLIC :: qs_upf_inverse_augmentation_fourier_collocate, qs_upf_inverse_augmentation_fourier_integrate
   PUBLIC :: qs_upf_inverse_augmentation_fourier_response

CONTAINS

! **************************************************************************************************
!> \brief Own original Q caches in the native magnetic-projector order, without changing QS state.
!> \param ks_env Live native KS environment used for frame.
!> \param frame Prepared inverse frame with selected nonghost USPP centers.
!> \param augmentation New owned original radial/angular caches.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_prepare(ks_env, frame, augmentation)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(OUT)                                     :: augmentation

      INTEGER                                            :: center, first, ikind, last
      REAL(dp)                                           :: metric_difference
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CPASSERT((frame%atom_local .OR. frame%inverse%ready))
      CALL get_ks_env(ks_env, qs_kind_set=kinds)
      CPASSERT(ASSOCIATED(kinds))
      ALLOCATE (augmentation%kinds(SIZE(kinds)))
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         IF (augmentation%kinds(ikind)%nproj /= 0) CYCLE
         CALL uspp_augmentation_init(kinds(ikind)%sgp_potential%upf_data, &
                                     kinds(ikind)%sgp_potential%upf_projectors, augmentation%kinds(ikind))
         CALL uspp_augmentation_prepare_realspace(kinds(ikind)%sgp_potential%upf_data, augmentation%kinds(ikind))
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         CPASSERT(augmentation%kinds(ikind)%nproj == last - first + 1)
         metric_difference = MAXVAL(ABS(augmentation%kinds(ikind)%q - REAL(frame%radial_q(first:last, first:last), dp)))
         CPASSERT(metric_difference == 0.0_dp)
      END DO
   END SUBROUTINE qs_upf_inverse_augmentation_prepare

! **************************************************************************************************
!> \brief Add original charge augmentation and its Cartesian gradient to matching local point rows.
!>        Only same-center projector blocks enter Q; inverse intercenter coupling remains in W.
!> \param frame ...
!> \param augmentation ...
!> \param projector_density Full Hermitian W P W^dagger, including occupations and k weights.
!> \param density Local density, incremented without weights.
!> \param gradient Local physical Cartesian density gradient, axis last, incremented.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_collocate(frame, augmentation, projector_density, density, gradient)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: projector_density(:, :)
      REAL(dp), INTENT(INOUT)                            :: density(:), gradient(:, :)

      INTEGER                                            :: atom, center, first, ikind, image, last, &
                                                            n
      REAL(dp)                                           :: hermiticity, position(3), tolerance
      REAL(dp), ALLOCATABLE                              :: displacements(:, :), g(:, :)

      CPASSERT((frame%atom_local .OR. frame%inverse%ready))
      n = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_density) == [n, n]))
      CPASSERT(SIZE(density) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(gradient) == [SIZE(density), 3]))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(projector_density)))
      hermiticity = MAXVAL(ABS(projector_density - CONJG(TRANSPOSE(projector_density))))
      CPASSERT(hermiticity <= tolerance)
      ALLOCATE (displacements(3, SIZE(density)), g(3, SIZE(density)))
      g = 0.0_dp
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            displacements(:, :) = frame%points - SPREAD(position, 2, SIZE(density))
            CALL uspp_augmentation_collocate_atom(augmentation%kinds(ikind), displacements, &
                                                  REAL(projector_density(first:last, first:last), dp), density, g)
         END DO
      END DO
      gradient = gradient + TRANSPOSE(g)
   END SUBROUTINE qs_upf_inverse_augmentation_collocate

! **************************************************************************************************
!> \brief Integrate the adjoint of the same original charge/gradient image sum, with one reduction.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity Unweighted potential conjugate to the augmented density.
!> \param vgradient Unweighted potential conjugate to its Cartesian gradient, axis last.
!> \param projector_operator New full projector operator; original Q contributes only same-center blocks.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_integrate(frame, augmentation, vdensity, vgradient, projector_operator)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN)                               :: vdensity(:), vgradient(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: projector_operator(:, :)

      INTEGER                                            :: atom, center, first, ikind, image, last
      REAL(dp)                                           :: position(3)
      REAL(dp), ALLOCATABLE                              :: displacements(:, :), matrix(:, :)

      CALL check_operator(frame, vdensity, vgradient, projector_operator)
      projector_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (displacements(3, SIZE(frame%weights)))
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (matrix(last - first + 1, last - first + 1))
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            displacements(:, :) = frame%points - SPREAD(position, 2, SIZE(frame%weights))
            CALL uspp_augmentation_integrate_atom(augmentation%kinds(ikind), displacements, frame%weights, &
                                                  vdensity, TRANSPOSE(vgradient), matrix)
            projector_operator(first:last, first:last) = projector_operator(first:last, first:last) + matrix
         END DO
         DEALLOCATE (matrix)
      END DO
      CALL frame%point_comm%sum(projector_operator)
   END SUBROUTINE qs_upf_inverse_augmentation_integrate

! **************************************************************************************************
!> \brief Add periodic original Q through native FFTs; its zero mode retains the reference radial q.
!>        Charge has no Bloch phase. The real-grid projection precedes each spectral gradient.
!> \param frame Prepared inverse frame on this physical-cell PW grid.
!> \param augmentation Original kernels; only their bounded Fourier caches are updated.
!> \param projector_density Full Hermitian W P W^dagger, with occupations and k weights.
!> \param density Matching native real density, incremented.
!> \param gradient Matching physical Cartesian density gradient, incremented.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_fourier_collocate(frame, augmentation, projector_density, density, gradient)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: projector_density(:, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, gradient(3)

      INTEGER                                            :: atom, axis, center, first, ikind, last, &
                                                            orders(3)
      REAL(dp)                                           :: hermiticity, tolerance
      TYPE(pw_c1d_gs_type)                               :: charge_g, derivative_g
      TYPE(pw_r3d_rs_type)                               :: work_r

      CALL check_fourier_grid(frame, density%pw_grid)
      CPASSERT(ALL(SHAPE(projector_density) == SIZE(frame%coefficients, 1)))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(projector_density)))
      hermiticity = MAXVAL(ABS(projector_density - CONJG(TRANSPOSE(projector_density))))
      CPASSERT(hermiticity <= tolerance)
      CALL charge_g%create(density%pw_grid)
      CALL derivative_g%create(density%pw_grid)
      CALL work_r%create(density%pw_grid)
      charge_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         CALL uspp_augmentation_prepare_fourier(augmentation%kinds(ikind), density%pw_grid%g)
         CALL uspp_augmentation_collocate(augmentation%kinds(ikind), density%pw_grid%g, &
                                          frame%positions(:, atom), density%pw_grid%vol, &
                                          REAL(projector_density(first:last, first:last), dp), charge_g%array)
      END DO
      CALL pw_transfer(charge_g, work_r)
      density%array = density%array + work_r%array
      ! Retain the same real function for gradients, including the FULLSPACE Nyquist projection.
      CALL pw_transfer(work_r, charge_g)
      DO axis = 1, 3
         CPASSERT(ASSOCIATED(gradient(axis)%pw_grid, density%pw_grid))
         orders = 0
         orders(axis) = 1
         derivative_g%array = charge_g%array
         CALL pw_derive(derivative_g, orders)
         CALL pw_transfer(derivative_g, work_r)
         gradient(axis)%array = gradient(axis)%array + work_r%array
      END DO
      CALL charge_g%release()
      CALL derivative_g%release()
      CALL work_r%release()
   END SUBROUTINE qs_upf_inverse_augmentation_fourier_collocate

! **************************************************************************************************
!> \brief Adjoint of the native real-grid Fourier charge/gradient map.
!>        Real-grid projection also enters -div(vgradient), including FULLSPACE Nyquist modes.
!>        The original Q/q data are unchanged; unweighted conjugate potentials are required.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity Potential conjugate to augmented density.
!> \param vgradient Potentials conjugate to physical Cartesian grad(n), not sigma.
!> \param projector_operator Full original-Q projector operator, overwritten and reduced once.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_fourier_integrate(frame, augmentation, vdensity, vgradient, projector_operator)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3)
      COMPLEX(dp), INTENT(OUT)                           :: projector_operator(:, :)

      INTEGER                                            :: atom, center, first, ikind, last
      REAL(dp), ALLOCATABLE                              :: matrix(:, :), multiplicity(:)
      TYPE(pw_c1d_gs_type)                               :: potential_g

      CALL check_fourier_grid(frame, vdensity%pw_grid)
      CPASSERT(ALL(SHAPE(projector_operator) == SIZE(frame%coefficients, 1)))
      CALL potential_g%create(vdensity%pw_grid)
      CALL fourier_potential(vdensity, vgradient, potential_g)
      ALLOCATE (multiplicity(SIZE(potential_g%array)))
      multiplicity = 1.0_dp
      IF (vdensity%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (vdensity%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      projector_operator = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (matrix(last - first + 1, last - first + 1))
         CALL uspp_augmentation_prepare_fourier(augmentation%kinds(ikind), vdensity%pw_grid%g)
         CALL uspp_augmentation_integrate(augmentation%kinds(ikind), vdensity%pw_grid%g, &
                                          frame%positions(:, atom), vdensity%pw_grid%vol, potential_g%array, multiplicity, matrix)
         projector_operator(first:last, first:last) = matrix
         DEALLOCATE (matrix)
      END DO
      IF (vdensity%pw_grid%para%mode /= PW_MODE_LOCAL) CALL frame%point_comm%sum(projector_operator)
      CALL potential_g%release()
   END SUBROUTINE qs_upf_inverse_augmentation_fourier_integrate

! **************************************************************************************************
!> \brief Geometry/potential response of the same native Fourier charge/gradient adjoint.
!>        Reciprocal indices and real-grid rows are fixed; points follow the cell affinely.
!>        Includes moving centers, Cartesian reciprocal vectors, volume measure and gradient adjoint.
!>        W and its response enter the subsequent AO pullback, not this original-Q operator.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity Unweighted conjugate potential on fixed fractional grid rows.
!> \param vgradient Potentials conjugate to the physical Cartesian density gradient.
!> \param atom_response Cartesian responses of all native atoms.
!> \param cell_response Cartesian lattice-vector response.
!> \param vdensity_response Response of the unweighted potential at each moving grid row.
!> \param vgradient_response Response of the Cartesian conjugate potentials at each grid row.
!> \param response Original-Q projector operator response, overwritten and reduced once.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_fourier_response(frame, augmentation, vdensity, vgradient, &
                                                           atom_response, cell_response, vdensity_response, &
                                                           vgradient_response, response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity_response, vgradient_response(3)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)

      INTEGER                                            :: atom, axis, center, first, ikind, info, &
                                                            last, second
      REAL(dp)                                           :: cell_inverse(3, 3), displacement(3), &
                                                            strain(3, 3), volume_response
      REAL(dp), ALLOCATABLE :: derivative(:, :, :), g_response(:, :), matrix(:, :), &
         multiplicity(:), potential_matrix(:, :), strain_matrix(:, :, :, :)
      TYPE(pw_c1d_gs_type)                               :: derivative_g, potential_g, response_g
      TYPE(pw_r3d_rs_type)                               :: work_r

      CALL check_fourier_grid(frame, vdensity%pw_grid)
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(ALL(SHAPE(response) == SIZE(frame%coefficients, 1)))
      CPASSERT(ASSOCIATED(vdensity_response%pw_grid, vdensity%pw_grid))
      cell_inverse = frame%cell
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      ALLOCATE (g_response(3, SIZE(vdensity%pw_grid%g, 2)))
      g_response(:, :) = -MATMUL(TRANSPOSE(strain), vdensity%pw_grid%g)
      CALL potential_g%create(vdensity%pw_grid)
      CALL response_g%create(vdensity%pw_grid)
      CALL derivative_g%create(vdensity%pw_grid)
      CALL work_r%create(vdensity%pw_grid)
      CALL fourier_potential(vdensity, vgradient, potential_g)
      CALL fourier_potential(vdensity_response, vgradient_response, response_g)
      ! Differentiate -div(vgradient) itself. Its real-grid projection is the same as in the adjoint.
      DO axis = 1, 3
         CALL pw_transfer(vgradient(axis), derivative_g)
         derivative_g%array = derivative_g%array*CMPLX(0.0_dp, -g_response(axis, :), dp)
         CALL pw_transfer(derivative_g, work_r)
         CALL pw_transfer(work_r, derivative_g)
         response_g%array = response_g%array + derivative_g%array
      END DO
      ALLOCATE (multiplicity(SIZE(potential_g%array)))
      multiplicity = 1.0_dp
      IF (vdensity%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (vdensity%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      response = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         ALLOCATE (matrix(last - first + 1, last - first + 1), &
                   potential_matrix(last - first + 1, last - first + 1), &
                   derivative(last - first + 1, last - first + 1, 3), &
                   strain_matrix(last - first + 1, last - first + 1, 3, 3))
         CALL uspp_augmentation_prepare_fourier(augmentation%kinds(ikind), vdensity%pw_grid%g)
         CALL uspp_augmentation_integrate(augmentation%kinds(ikind), vdensity%pw_grid%g, &
                                          frame%positions(:, atom), vdensity%pw_grid%vol, potential_g%array, &
                                          multiplicity, matrix, derivative, strain_matrix)
         CALL uspp_augmentation_integrate(augmentation%kinds(ikind), vdensity%pw_grid%g, &
                                          frame%positions(:, atom), vdensity%pw_grid%vol, response_g%array, &
                                          multiplicity, potential_matrix)
         displacement = atom_response(:, atom) - MATMUL(strain, frame%positions(:, atom))
         matrix(:, :) = volume_response*matrix + potential_matrix
         DO axis = 1, 3
            matrix(:, :) = matrix + displacement(axis)*derivative(:, :, axis)
            DO second = 1, 3
               matrix(:, :) = matrix + strain(axis, second)*strain_matrix(:, :, axis, second)
            END DO
         END DO
         response(first:last, first:last) = matrix
         DEALLOCATE (matrix, potential_matrix, derivative, strain_matrix)
      END DO
      IF (vdensity%pw_grid%para%mode /= PW_MODE_LOCAL) CALL frame%point_comm%sum(response)
      CALL potential_g%release()
      CALL response_g%release()
      CALL derivative_g%release()
      CALL work_r%release()
   END SUBROUTINE qs_upf_inverse_augmentation_fourier_response

! **************************************************************************************************
!> \brief Fold the projected spectral gradient adjoint into an unweighted scalar potential.
!> \param vdensity ...
!> \param vgradient ...
!> \param potential_g Allocated matching Fourier field, overwritten.
! **************************************************************************************************
   SUBROUTINE fourier_potential(vdensity, vgradient, potential_g)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: potential_g

      INTEGER                                            :: axis, orders(3)
      TYPE(pw_c1d_gs_type)                               :: derivative_g
      TYPE(pw_r3d_rs_type)                               :: work_r

      CPASSERT(ASSOCIATED(potential_g%pw_grid, vdensity%pw_grid))
      CALL derivative_g%create(vdensity%pw_grid)
      CALL work_r%create(vdensity%pw_grid)
      CALL pw_transfer(vdensity, potential_g)
      DO axis = 1, 3
         CPASSERT(ASSOCIATED(vgradient(axis)%pw_grid, vdensity%pw_grid))
         orders = 0
         orders(axis) = 1
         CALL pw_transfer(vgradient(axis), derivative_g)
         CALL pw_derive(derivative_g, orders)
         CALL pw_transfer(derivative_g, work_r)
         CALL pw_transfer(work_r, derivative_g)
         potential_g%array = potential_g%array - derivative_g%array
      END DO
      CALL derivative_g%release()
      CALL work_r%release()
   END SUBROUTINE fourier_potential

! **************************************************************************************************
!> \brief Require the native real/reciprocal FFT layout of the same physical-cell point frame.
!> \param frame ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE check_fourier_grid(frame, grid)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid

      INTEGER                                            :: axis
      REAL(dp)                                           :: difference, tolerance

      CPASSERT((frame%atom_local .OR. frame%inverse%ready) .AND. ASSOCIATED(grid))
      CPASSERT(ASSOCIATED(grid%g) .AND. grid%vol > 0.0_dp)
      CPASSERT(SIZE(frame%weights) == grid%ngpts_local)
      CPASSERT(ALL(frame%weights == grid%dvol))
      CPASSERT(grid%grid_span == FULLSPACE .OR. grid%grid_span == HALFSPACE)
      IF (grid%grid_span == HALFSPACE) THEN
         CPASSERT(ALL(MODULO(grid%npts, 2) == 1))
      END IF
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(frame%cell)))
      DO axis = 1, 3
         difference = MAXVAL(ABS(grid%dh(:, axis)*grid%npts(axis) - frame%cell(:, axis)))
         CPASSERT(difference <= tolerance)
      END DO
      IF (grid%para%mode /= PW_MODE_LOCAL) THEN
         CPASSERT(grid%para%group%mp_comm_type == frame%point_comm)
      ELSE
         CPASSERT(frame%point_comm%num_pe == 1)
      END IF
   END SUBROUTINE check_fourier_grid

! **************************************************************************************************
!> \brief Differentiate the same Q image sum, including its spatial gradient and original measure.
!>        Image membership is fixed. Inverse-coefficient responses belong to the later W pullback.
!> \param frame ...
!> \param augmentation ...
!> \param vdensity ...
!> \param vgradient ...
!> \param point_response Local Cartesian point responses, axis first.
!> \param atom_response Responses of all native atom positions, axis first.
!> \param cell_response Physical-cell response for the same fixed images.
!> \param weights_response Response of the same original local quadrature measure.
!> \param vdensity_response Explicit response of the conjugate potential.
!> \param vgradient_response ...
!> \param response New full projector-operator response, reduced once over point rows.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_response(frame, augmentation, vdensity, vgradient, &
                                                   point_response, atom_response, cell_response, weights_response, &
                                                   vdensity_response, vgradient_response, response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(qs_upf_inverse_augmentation_type), INTENT(IN) :: augmentation
      REAL(dp), INTENT(IN) :: vdensity(:), vgradient(:, :), point_response(:, :), &
         atom_response(:, :), cell_response(3, 3), weights_response(:), vdensity_response(:), &
         vgradient_response(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :)

      INTEGER                                            :: atom, axis, center, first, ikind, image, &
                                                            last, n, point
      REAL(dp)                                           :: dr(3), position(3), r(3), radius2
      REAL(dp), ALLOCATABLE                              :: dq(:, :, :), dterm(:, :), &
                                                            hq(:, :, :, :), q(:, :), term(:, :)

      CALL check_operator(frame, vdensity, vgradient, response)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      CPASSERT(SIZE(vdensity_response) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(vgradient_response) == SHAPE(vgradient)))
      response = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         n = last - first + 1
         radius2 = augmentation%kinds(ikind)%realspace_radius**2
         ALLOCATE (q(n, n), dq(n, n, 3), hq(n, n, 3, 3), term(n, n), dterm(n, n))
         DO image = 1, SIZE(frame%images, 2)
            position = frame%positions(:, atom) + MATMUL(frame%cell, REAL(frame%images(:, image), dp))
            DO point = 1, SIZE(frame%weights)
               r = frame%points(:, point) - position
               IF (DOT_PRODUCT(r, r) > radius2) CYCLE
               dr = point_response(:, point) - atom_response(:, atom) - MATMUL(cell_response, REAL(frame%images(:, image), dp))
               CALL uspp_augmentation_realspace(augmentation%kinds(ikind), r, q, dq, hq)
               term(:, :) = vdensity(point)*q
               dterm(:, :) = vdensity_response(point)*q
               DO axis = 1, 3
                  term(:, :) = term + vgradient(point, axis)*dq(:, :, axis)
                  dterm(:, :) = dterm + vdensity(point)*dq(:, :, axis)*dr(axis) + &
                                vgradient_response(point, axis)*dq(:, :, axis)
                  dterm(:, :) = dterm + vgradient(point, axis)*(hq(:, :, axis, 1)*dr(1) + &
                                                                hq(:, :, axis, 2)*dr(2) + hq(:, :, axis, 3)*dr(3))
               END DO
               response(first:last, first:last) = response(first:last, first:last) + &
                                                  frame%weights(point)*dterm + weights_response(point)*term
            END DO
         END DO
         DEALLOCATE (q, dq, hq, term, dterm)
      END DO
      CALL frame%point_comm%sum(response)
   END SUBROUTINE qs_upf_inverse_augmentation_response

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param vdensity ...
!> \param vgradient ...
!> \param projector_operator ...
! **************************************************************************************************
   SUBROUTINE check_operator(frame, vdensity, vgradient, projector_operator)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN)                               :: vdensity(:), vgradient(:, :)
      COMPLEX(dp), INTENT(IN)                            :: projector_operator(:, :)

      CPASSERT((frame%atom_local .OR. frame%inverse%ready))
      CPASSERT(SIZE(vdensity) == SIZE(frame%weights))
      CPASSERT(ALL(SHAPE(vgradient) == [SIZE(frame%weights), 3]))
      CPASSERT(ALL(SHAPE(projector_operator) == SIZE(frame%coefficients, 1)))
   END SUBROUTINE check_operator

! **************************************************************************************************
!> \brief Release owned radial/angular caches.
!> \param augmentation ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_augmentation_release(augmentation)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation

      TYPE(qs_upf_inverse_augmentation_type)             :: empty

      augmentation = empty
   END SUBROUTINE qs_upf_inverse_augmentation_release

END MODULE qs_upf_inverse_augmentation

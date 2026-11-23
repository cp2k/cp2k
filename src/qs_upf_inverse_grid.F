!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse USPP fields on native PW rows and real/imaginary Quickstep AO matrices.
!>        This development adapter replicates AO columns and dense AO contractions.
!>        It does not change the reference SCF representation or supply missing augmentation.
! **************************************************************************************************
MODULE qs_upf_inverse_grid
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_desymmetrize, dbcsr_finalize, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, &
        dbcsr_release, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED,&
                                              PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_from_ks,&
                                              qs_upf_inverse_release,&
                                              qs_upf_inverse_response,&
                                              qs_upf_inverse_type,&
                                              split_reconstruction
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_collocate,&
                                              qs_upf_inverse_augmentation_fourier_collocate,&
                                              qs_upf_inverse_augmentation_fourier_integrate,&
                                              qs_upf_inverse_augmentation_fourier_response,&
                                              qs_upf_inverse_augmentation_integrate,&
                                              qs_upf_inverse_augmentation_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_grid_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_inverse,                     ONLY: upf_inverse_basis,&
                                              upf_inverse_fields_init
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response,&
                                              upf_inverse_semilocal_pair,&
                                              upf_inverse_semilocal_response,&
                                              upf_inverse_smooth_density
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Composite midpoint quadrature of the compact reconstruction products.
   INTEGER, PARAMETER :: subdivisions = 2, interpolation_order = 8

   PUBLIC :: qs_upf_inverse_grid_type, qs_upf_inverse_grid_prepare, qs_upf_inverse_grid_release
   PUBLIC :: qs_upf_inverse_grid_collocate, qs_upf_inverse_grid_integrate
   PUBLIC :: qs_upf_inverse_grid_metric
   PUBLIC :: qs_upf_inverse_grid_collocate_dense, qs_upf_inverse_grid_integrate_dense, qs_upf_inverse_grid_metric_dense
   PUBLIC :: qs_upf_inverse_matrix_gather, qs_upf_inverse_matrix_scatter
   PUBLIC :: qs_upf_inverse_grid_response

CONTAINS

! **************************************************************************************************
!> \brief Prepare original inverse fields with the native grid coordinates and volume element.
!>        A local PW grid is complete on each rank and uses SELF for point integrals.
!>        A distributed PW grid uses its own communicator, with each point counted once.
!> \param ks_env Prepared native KS environment; its ORB bases remain borrowed.
!> \param grid Native PW grid, including its auxiliary-cell geometry when applicable.
!> \param images Explicit images of the physical KS cell, held fixed for responses.
!> \param reduced_k Reduced Bloch vector of the physical KS cell.
!> \param frame New owned adapter; neither the grid nor the KS environment is modified.
!> \param info Inverse-factorization status.
!> \param atom_local Optional nonoverlapping one-center inverse.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_prepare(ks_env, grid, images, reduced_k, frame, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(qs_upf_inverse_grid_type), INTENT(OUT)        :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      INTEGER                                            :: axis, ix, iy, iz, n, point
      INTEGER, ALLOCATABLE                               :: patch_center(:)
      LOGICAL                                            :: refine
      REAL(dp), ALLOCATABLE                              :: metric_weights(:), points(:, :), &
                                                            weights(:)
      TYPE(mp_comm_type)                                 :: comm

      frame%bounds = grid%bounds
      frame%bounds_local = grid%bounds_local
      frame%npts = grid%npts
      frame%dh = grid%dh
      frame%dvol = grid%dvol
      frame%mode = grid%para%mode
      CPASSERT(frame%dvol > 0.0_dp .AND. ALL(frame%npts > 0))
      n = PRODUCT(MAX(0, frame%bounds_local(2, :) - frame%bounds_local(1, :) + 1))
      ! Reciprocal blocks need not have the size of the local real-space rectangle.
      CPASSERT(n == PRODUCT(MAX(0, grid%npts_local)))
      IF (frame%mode == PW_MODE_LOCAL) THEN
         CPASSERT(ALL(frame%bounds_local == frame%bounds))
         comm = mp_comm_self
      ELSE
         CPASSERT(frame%mode == PW_MODE_DISTRIBUTED)
         comm = grid%para%group%mp_comm_type
      END IF
      ALLOCATE (points(3, n), weights(n))
      weights = frame%dvol
      point = 0
      DO iz = frame%bounds_local(1, 3), frame%bounds_local(2, 3)
         DO iy = frame%bounds_local(1, 2), frame%bounds_local(2, 2)
            DO ix = frame%bounds_local(1, 1), frame%bounds_local(2, 1)
               point = point + 1
               points(:, point) = MATMUL(frame%dh, REAL([ix, iy, iz] - frame%bounds(1, :), dp))
            END DO
         END DO
      END DO
      refine = .TRUE.
      IF (PRESENT(direct)) refine = .NOT. direct
      IF (refine) THEN
         CALL reconstruction_quadrature(frame, ks_env, points, weights, metric_weights, patch_center)
      END IF
      CALL qs_upf_inverse_from_ks(ks_env, points, weights, images, reduced_k, comm, frame%ao, info, atom_local, direct)
      IF (info /= 0 .OR. .NOT. refine) RETURN
      frame%ao%coarse_points = n
      frame%ao%patch_points = SIZE(patch_center)
      CALL MOVE_ALLOC(patch_center, frame%ao%patch_center)
      CALL MOVE_ALLOC(metric_weights, frame%ao%metric_weights)
      CALL split_reconstruction(frame%ao, frame%ao%delta)
      CALL upf_inverse_basis(frame%ao%gaussian, frame%ao%delta, frame%ao%coefficients, frame%ao%soft)
      DO axis = 1, 3
         CALL split_reconstruction(frame%ao, frame%ao%delta_gradient(:, :, axis))
         CALL upf_inverse_basis(frame%ao%gaussian_gradient(:, :, axis), frame%ao%delta_gradient(:, :, axis), &
                                frame%ao%coefficients, frame%ao%soft_gradient(:, :, axis))
      END DO
   END SUBROUTINE qs_upf_inverse_grid_prepare

! **************************************************************************************************
!> \brief Sample each compact density correction on a translated fine lattice centered on its atom.
!>        The transpose of the potential interpolation deposits the density on the GPW grid.
!> \param frame Native grid geometry and local indices of the interpolation stencil.
!> \param ks_env Atomic positions and original UPF radial fields.
!> \param points Global grid samples followed by the two copies of each atom-centered sample.
!> \param weights Signed volume weights used for the density and potential.
!> \param metric_weights Weights with unique ownership across overlapping MPI halos.
!> \param patch_center Reconstruction center of each additional quadrature point.
! **************************************************************************************************
   SUBROUTINE reconstruction_quadrature(frame, ks_env, points, weights, metric_weights, patch_center)
      TYPE(qs_upf_inverse_grid_type), INTENT(INOUT)      :: frame
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: points(:, :), weights(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: metric_weights(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: patch_center(:)

      INTEGER :: axis, base(3), box(3), center, count, extent(3), &
         grid_indices(interpolation_order, 3), i, info, j, kind, n, nc, owner(3), pass, q(3), &
         step(3), x, y, z
      REAL(dp)                                           :: coord(3), fraction, inverse_dh(3, 3), &
                                                            offset(3), origin(3), radius, width(3)
      REAL(dp), ALLOCATABLE                              :: radii(:), refined(:, :)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_radial_fields_type)                       :: beta, delta

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, uspp=uspp)
      CPASSERT(ASSOCIATED(kinds) .AND. ASSOCIATED(atoms) .AND. ASSOCIATED(uspp))
      nc = SIZE(uspp%center_atoms)
      ALLOCATE (radii(SIZE(kinds)), SOURCE=-1.0_dp)
      DO center = 1, nc
         kind = uspp%center_kinds(center)
         IF (radii(kind) >= 0.0_dp) CYCLE
         CALL upf_inverse_fields_init(kinds(kind)%sgp_potential%upf_data, beta, delta)
         radii(kind) = delta%radius
      END DO
      n = SIZE(weights)
      inverse_dh = frame%dh
      CALL invmat(inverse_dh, info)
      CPASSERT(info == 0)
      width = SQRT(SUM(inverse_dh**2, DIM=2))
      extent = subdivisions*frame%npts
      ALLOCATE (frame%stencil(interpolation_order, subdivisions, 3, nc), &
                frame%dstencil(interpolation_order, subdivisions, 3, nc))
      DO pass = 1, 2
         count = 0
         DO center = 1, nc
            radius = radii(uspp%center_kinds(center))
            IF (radius <= 0.0_dp) CYCLE
            origin = MATMUL(inverse_dh, atoms(uspp%center_atoms(center))%r)
            box = CEILING(radius*width*subdivisions + 0.5_dp)
            IF (pass == 1) THEN
               DO axis = 1, 3
                  DO j = 1, subdivisions
                     fraction = MODULO(origin(axis) + (j - 0.5_dp)/subdivisions, 1.0_dp)
                     CALL cardinal_stencil(fraction, frame%stencil(:, j, axis, center), frame%dstencil(:, j, axis, center))
                  END DO
               END DO
            END IF
            DO z = -box(3), box(3)
               DO y = -box(2), box(2)
                  DO x = -box(1), box(1)
                     q = [x, y, z]
                     offset = (REAL(q, dp) + 0.5_dp)/subdivisions
                     IF (SUM(MATMUL(frame%dh, offset)**2) >= radius**2) CYCLE
                     IF (ANY(2*box + 1 > extent)) THEN
                        IF (previous_image(q, box, extent, frame%dh, radius)) CYCLE
                     END IF
                     coord = origin + offset
                     step = MODULO(q, subdivisions) + 1
                     ! Use the same fractional argument as the stencil, also at a grid plane.
                     base = FLOOR(origin + (REAL(step, dp) - 0.5_dp)/subdivisions) + (q - (step - 1))/subdivisions
                     DO axis = 1, 3
                        DO j = 1, interpolation_order
                           i = MODULO(base(axis) + j - interpolation_order/2, frame%npts(axis)) + frame%bounds(1, axis)
                           IF (i < frame%bounds_local(1, axis) .OR. i > frame%bounds_local(2, axis)) i = HUGE(0)
                           grid_indices(j, axis) = i
                        END DO
                     END DO
                     IF (ANY(ALL(grid_indices == HUGE(0), DIM=1))) CYCLE
                     count = count + 1
                     IF (pass == 1) CYCLE
                     frame%patch_index(:, :, count) = grid_indices
                     frame%patch_offset(:, count) = step
                     patch_center(count) = center
                     refined(:, n + count) = MATMUL(frame%dh, MODULO(coord, REAL(frame%npts, dp)))
                     ! A unique coarse node owns each metric contribution across MPI quadrature halos.
                     owner = MODULO(base, frame%npts) + frame%bounds(1, :)
                     IF (ALL(owner >= frame%bounds_local(1, :)) .AND. ALL(owner <= frame%bounds_local(2, :))) THEN
                        metric_weights(n + count) = frame%dvol/REAL(subdivisions**3, dp)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (pass == 1) THEN
            ALLOCATE (frame%patch_index(interpolation_order, 3, count), frame%patch_offset(3, count), patch_center(count), &
                      refined(3, n + 2*count), metric_weights(n + 2*count))
            refined(:, :n) = points
            metric_weights = 0.0_dp
         END IF
      END DO
      refined(:, n + count + 1:) = refined(:, n + 1:n + count)
      CALL MOVE_ALLOC(refined, points)
      DEALLOCATE (weights)
      ALLOCATE (weights(n + 2*count))
      weights(:n) = frame%dvol
      weights(n + 1:n + count) = frame%dvol/REAL(subdivisions**3, dp)
      weights(n + count + 1:) = -weights(n + 1:n + count)
      metric_weights(:n) = weights(:n)
      metric_weights(n + count + 1:) = -metric_weights(n + 1:n + count)
   END SUBROUTINE reconstruction_quadrature

! **************************************************************************************************
!> \brief Keep each periodic fine node once when a reconstruction sphere intersects its own images.
!> \param q Fine-grid offset within the bounding box of the sphere.
!> \param box Maximum offsets in each direction.
!> \param extent Fine-grid periods.
!> \param dh Coarse-grid basis vectors.
!> \param radius Support of the radial reconstruction.
!> \retval duplicate An earlier image contains the same quadrature point.
! **************************************************************************************************
   LOGICAL FUNCTION previous_image(q, box, extent, dh, radius) RESULT(duplicate)
      INTEGER, INTENT(IN)                                :: q(3), box(3), extent(3)
      REAL(dp), INTENT(IN)                               :: dh(3, 3), radius

      INTEGER                                            :: first(3), last(3), q2(3), x, y, z

      duplicate = .FALSE.
      first = -((box + q)/extent)
      last = (box - q)/extent
      DO z = first(3), MIN(0, last(3))
         DO y = first(2), last(2)
            DO x = first(1), last(1)
               IF (z == 0 .AND. (y > 0 .OR. (y == 0 .AND. x >= 0))) CYCLE
               q2 = q + extent*[x, y, z]
               IF (SUM(MATMUL(dh, (REAL(q2, dp) + 0.5_dp)/subdivisions)**2) < radius**2) THEN
                  duplicate = .TRUE.
                  RETURN
               END IF
            END DO
         END DO
      END DO
   END FUNCTION previous_image

! **************************************************************************************************
!> \brief Interpolate the resolved density spectrum with a moment-corrected cardinal B-spline.
!>        Values and two derivatives agree as a moving atom crosses a lattice plane.
!> \param t Fractional position between neighboring grid nodes.
!> \param weights Interpolation weights on the eight nearest nodes.
!> \param derivatives Derivative of the weights with respect to t.
! **************************************************************************************************
   SUBROUTINE cardinal_stencil(t, weights, derivatives)
      REAL(dp), INTENT(IN)                               :: t
      REAL(dp), INTENT(OUT) :: weights(interpolation_order), derivatives(interpolation_order)

      REAL(dp), PARAMETER :: second_moment = REAL(interpolation_order, dp)/24.0_dp

      INTEGER                                            :: cell, degree, j
      REAL(dp) :: b(-interpolation_order:3, 0:interpolation_order - 1), curvature, third

      ! Cardinal B-splines with their second moment removed. Low spatial frequencies
      ! are reproduced through cubic order while unresolved frequencies are attenuated.
      degree = interpolation_order - 1
      b = 0.0_dp
      b(0, 0) = 1.0_dp
      DO j = 1, degree
         DO cell = -j, 0
            b(cell, j) = ((t - cell)*b(cell, j - 1) + (cell + j + 1 - t)*b(cell + 1, j - 1))/j
         END DO
      END DO
      DO j = 1, interpolation_order
         cell = j - interpolation_order
         curvature = b(cell, degree - 2) - 2*b(cell + 1, degree - 2) + b(cell + 2, degree - 2)
         third = b(cell, degree - 3) - 3*b(cell + 1, degree - 3) + 3*b(cell + 2, degree - 3) - b(cell + 3, degree - 3)
         weights(j) = b(cell, degree) - second_moment*curvature
         derivatives(j) = b(cell, degree - 1) - b(cell + 1, degree - 1) - second_moment*third
      END DO
   END SUBROUTINE cardinal_stencil

! **************************************************************************************************
!> \brief Collocate one spin/k contribution and its projector density in the same inverse frame.
!>        The caller supplies occupations/k weights in P and optional original charge augmentation.
!>        Outputs are overwritten; this routine neither sums k points nor updates QS validity flags.
!> \param frame Prepared native-grid inverse frame.
!> \param density_real Native real AO density component in reconstructed Gaussian variables.
!> \param density Density on an allocated native PW grid; includes charge augmentation if supplied.
!> \param gradient Its three physical Cartesian gradient components.
!> \param tau Smooth kinetic density, including the factor one half.
!> \param density_imag Optional native imaginary AO density component, with its actual sign.
!> \param projector_density Optional W P W^dagger for remaining original augmentation.
!> \param augmentation Optional owned original scalar USPP charge/gradient kernels.
!>        Tau remains smooth; supplying charge augmentation does not supply kinetic augmentation.
!> \param fourier_charge Optional native Fourier Q map, instead of point Q; default is false.
!>        Its zero mode equals the reference radial q; a supplied-q metric retains the explicit defect.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_collocate(frame, density_real, density, gradient, tau, density_imag, &
                                            projector_density, augmentation, fourier_charge)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: density_real
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, gradient(3), tau
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: density_imag
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: projector_density(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge

      COMPLEX(dp), ALLOCATABLE                           :: p(:, :)

      CALL qs_upf_inverse_matrix_gather(frame%ao, density_real, p, density_imag)
      CALL qs_upf_inverse_grid_collocate_dense(frame, p, density, gradient, tau, projector_density, augmentation, fourier_charge)
   END SUBROUTINE qs_upf_inverse_grid_collocate

! **************************************************************************************************
!> \brief The same native density map for a replicated complex AO density matrix.
!>        The caller supplies occupations and any k weight exactly once.
!> \param frame ...
!> \param p Hermitian AO density in the inverse variables.
!> \param density ...
!> \param gradient ...
!> \param tau Smooth kinetic density only.
!> \param projector_density ...
!> \param augmentation ...
!> \param fourier_charge ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_collocate_dense(frame, p, density, gradient, tau, &
                                                  projector_density, augmentation, fourier_charge)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      COMPLEX(dp), INTENT(IN)                            :: p(:, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, gradient(3), tau
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: projector_density(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :)
      INTEGER                                            :: axis, n, nproj
      LOGICAL                                            :: use_fourier
      REAL(dp), ALLOCATABLE                              :: dn(:, :), rho(:), t(:)

      CPASSERT((frame%ao%atom_local .OR. frame%ao%inverse%ready))
      use_fourier = .FALSE.
      IF (PRESENT(fourier_charge)) use_fourier = fourier_charge
      IF (use_fourier) THEN
         CPASSERT(PRESENT(augmentation))
      END IF
      CPASSERT(ALL(SHAPE(p) == SIZE(frame%ao%soft, 2)))
      n = SIZE(frame%ao%weights)
      ALLOCATE (rho(n), dn(n, 3), t(n))
      IF (frame%native_gaussian) THEN
         rho = 0.0_dp
         dn = 0.0_dp
         t = 0.0_dp
      ELSE
         CALL upf_inverse_smooth_density(frame%ao%soft, frame%ao%soft_gradient, p, rho, dn, t)
      END IF
      IF (PRESENT(projector_density) .OR. PRESENT(augmentation)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         projected(:, :) = MATMUL(MATMUL(frame%ao%coefficients, p), CONJG(TRANSPOSE(frame%ao%coefficients)))
         IF (PRESENT(augmentation) .AND. .NOT. use_fourier) THEN
            CALL qs_upf_inverse_augmentation_collocate(frame%ao, augmentation, projected, rho, dn)
         END IF
         IF (PRESENT(projector_density)) THEN
            CPASSERT(ALL(SHAPE(projector_density) == [nproj, nproj]))
            projector_density = projected
         END IF
      END IF
      CALL unpack_grid(frame, rho, density)
      CALL unpack_grid(frame, t, tau)
      DO axis = 1, 3
         CALL unpack_grid(frame, dn(:, axis), gradient(axis))
      END DO
      IF (use_fourier) THEN
         CALL qs_upf_inverse_augmentation_fourier_collocate(frame%ao, augmentation, projected, density, gradient)
      END IF
   END SUBROUTINE qs_upf_inverse_grid_collocate_dense

! **************************************************************************************************
!> \brief Integrate the density/gradient/tau adjoint into owned native real S / imaginary A matrices.
!>        Inverse-mediated AO fill is retained, independently of the template's stored blocks.
!> \param frame ...
!> \param vdensity Unweighted potential conjugate to smooth density.
!> \param vgradient Three unweighted potentials conjugate to physical grad(n), not sigma.
!> \param vtau Unweighted potential conjugate to smooth tau.
!> \param template Native AO block layout/distribution; its values and sparsity are not copied.
!> \param operator_real New owned real symmetric matrix; release it after use.
!> \param operator_imag New owned imaginary antisymmetric matrix; release it after use.
!> \param projector_operator Optional explicitly supplied remaining projector operator.
!>        Missing charge or kinetic augmentation must not be inferred from the other kernel.
!> \param augmentation Optional original charge/gradient adjoint on the same native point rows.
!> \param fourier_charge Select the native Fourier adjoint of the same density map; default is false.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_integrate(frame, vdensity, vgradient, vtau, template, &
                                            operator_real, operator_imag, projector_operator, augmentation, fourier_charge)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: operator_real, operator_imag
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge

      COMPLEX(dp), ALLOCATABLE                           :: h(:, :)

      CALL qs_upf_inverse_grid_integrate_dense(frame, vdensity, vgradient, vtau, h, &
                                               projector_operator, augmentation, fourier_charge)
      CALL qs_upf_inverse_matrix_scatter(frame%ao, h, template, operator_real, operator_imag)
   END SUBROUTINE qs_upf_inverse_grid_integrate

! **************************************************************************************************
!> \brief Complex AO adjoint, before its native matrix distribution is selected.
!> \param frame ...
!> \param vdensity ...
!> \param vgradient ...
!> \param vtau ...
!> \param h New replicated Hermitian operator, without occupation or k weights.
!> \param projector_operator ...
!> \param augmentation ...
!> \param fourier_charge ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_integrate_dense(frame, vdensity, vgradient, vtau, h, &
                                                  projector_operator, augmentation, fourier_charge)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: h(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT), OPTIONAL                         :: augmentation
      LOGICAL, INTENT(IN), OPTIONAL                      :: fourier_charge

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :)
      INTEGER                                            :: axis, n, nao, nproj
      LOGICAL                                            :: use_fourier
      REAL(dp), ALLOCATABLE                              :: dg(:, :), dn(:), dt(:)

      CPASSERT((frame%ao%atom_local .OR. frame%ao%inverse%ready))
      use_fourier = .FALSE.
      IF (PRESENT(fourier_charge)) use_fourier = fourier_charge
      IF (use_fourier) THEN
         CPASSERT(PRESENT(augmentation))
      END IF
      n = SIZE(frame%ao%weights)
      nao = SIZE(frame%ao%soft, 2)
      ALLOCATE (dn(n), dg(n, 3), dt(n), h(nao, nao))
      CALL pack_grid(frame, vdensity, dn)
      CALL pack_grid(frame, vtau, dt)
      DO axis = 1, 3
         CALL pack_grid(frame, vgradient(axis), dg(:, axis))
      END DO
      IF (frame%native_gaussian) THEN
         h = 0.0_dp
      ELSE
         CALL upf_inverse_semilocal_pair(frame%ao%soft, frame%ao%soft, frame%ao%soft_gradient, frame%ao%soft_gradient, &
                                         frame%ao%weights, dn, dg, dt, h, frame%ao%point_comm)
      END IF
      IF (PRESENT(augmentation)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         ALLOCATE (projected(nproj, nproj))
         IF (use_fourier) THEN
            CALL qs_upf_inverse_augmentation_fourier_integrate(frame%ao, augmentation, vdensity, vgradient, projected)
         ELSE
            CALL qs_upf_inverse_augmentation_integrate(frame%ao, augmentation, dn, dg, projected)
         END IF
         h(:, :) = h + MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projected, frame%ao%coefficients))
      END IF
      IF (PRESENT(projector_operator)) THEN
         nproj = SIZE(frame%ao%coefficients, 1)
         CPASSERT(ALL(SHAPE(projector_operator) == [nproj, nproj]))
         h(:, :) = h + MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projector_operator, frame%ao%coefficients))
      END IF
   END SUBROUTINE qs_upf_inverse_grid_integrate_dense

! **************************************************************************************************
!> \brief Assemble the original generalized metric in the inverse AO variables.
!>        The point rule is shared with the density/operator map; no metric residual is removed.
!> \param frame ...
!> \param projector_metric Explicit original projector metric (supplied q or reference radial q).
!>        The caller selects the comparison functional; this routine does not substitute either.
!> \param template ...
!> \param metric_real New owned real symmetric metric.
!> \param metric_imag New owned imaginary antisymmetric metric.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_metric(frame, projector_metric, template, metric_real, metric_imag)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: metric_real, metric_imag

      COMPLEX(dp), ALLOCATABLE                           :: metric(:, :)

      CALL qs_upf_inverse_grid_metric_dense(frame, projector_metric, metric)
      CALL qs_upf_inverse_matrix_scatter(frame%ao, metric, template, metric_real, metric_imag)
   END SUBROUTINE qs_upf_inverse_grid_metric

! **************************************************************************************************
!> \brief Replicated complex inverse metric with the caller's explicit projector metric.
!> \param frame ...
!> \param projector_metric ...
!> \param metric New owned dense metric; radial Q and supplied PP_Q remain distinct.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_metric_dense(frame, projector_metric, metric)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: metric(:, :)

      INTEGER                                            :: nao, nproj

      CPASSERT((frame%ao%atom_local .OR. frame%ao%inverse%ready))
      nao = SIZE(frame%ao%soft, 2)
      nproj = SIZE(frame%ao%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [nproj, nproj]))
      ALLOCATE (metric(nao, nao))
      IF (frame%native_gaussian) THEN
         metric = 0.0_dp
      ELSE
         IF (ALLOCATED(frame%ao%metric_weights)) THEN
            CALL upf_inverse_pair(frame%ao%soft, frame%ao%soft, frame%ao%metric_weights, metric, frame%ao%point_comm)
         ELSE
            CALL upf_inverse_pair(frame%ao%soft, frame%ao%soft, frame%ao%weights, metric, frame%ao%point_comm)
         END IF
      END IF
      metric(:, :) = metric + MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), &
                                     MATMUL(projector_metric, frame%ao%coefficients))
   END SUBROUTINE qs_upf_inverse_grid_metric_dense

! **************************************************************************************************
!> \brief Joint native-grid geometry response of the inverse metric and original-Q semilocal operator.
!>        Fixed fractional grid rows, image membership and reciprocal indices define the derivative.
!>        Both moving inverse legs, their Cartesian gradients and the cell measure are differentiated.
!> \param frame ...
!> \param augmentation Original Fourier charge kernels; does not supply kinetic augmentation.
!> \param projector_metric Explicit original supplied-q or radial-Q metric, unchanged by geometry.
!> \param vdensity Unweighted conjugate potential.
!> \param vgradient Potentials conjugate to physical Cartesian grad(n).
!> \param vtau Potential conjugate to smooth tau, with the factor one half.
!> \param atom_response Cartesian response of every native atom.
!> \param cell_response Lattice-vector response.
!> \param reduced_k_response Response of the reduced Bloch vector, zero for physical cell stress.
!> \param vdensity_response Potential response at each moving grid row.
!> \param vgradient_response ...
!> \param vtau_response ...
!> \param operator_response New replicated complex operator response, without k/spin weights.
!> \param metric_response New metric response in the same inverse AO variables.
!> \param projector_operator Optional original geometry-independent operator, e.g. Dion.
!> \param dw Optional geometry derivative of the inverse projector coefficients used in the operator response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_response(frame, augmentation, projector_metric, vdensity, vgradient, vtau, &
                                           atom_response, cell_response, reduced_k_response, vdensity_response, &
                                           vgradient_response, vtau_response, operator_response, metric_response, &
                                           projector_operator, dw)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity_response, &
                                                            vgradient_response(3), vtau_response
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operator_response(:, :), &
                                                            metric_response(:, :)
      COMPLEX(dp), INTENT(IN), OPTIONAL                  :: projector_operator(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: dw(:, :)

      COMPLEX(dp), ALLOCATABLE :: coefficients_response(:, :), gradient_response(:, :, :), &
         projected(:, :), projected_response(:, :), soft_response(:, :)
      INTEGER                                            :: atom, axis, center, info, n, nao, np, &
                                                            nproj, p, row
      REAL(dp)                                           :: cell_inverse(3, 3), inverse_dh(3, 3), &
                                                            relative(3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: dg(:, :), dn(:), dt(:), dvg(:, :), &
                                                            dvn(:), dvt(:), motion(:, :), &
                                                            points_response(:, :), &
                                                            weights_response(:)

      CPASSERT((frame%ao%atom_local .OR. frame%ao%inverse%ready))
      n = SIZE(frame%ao%weights)
      nao = SIZE(frame%ao%soft, 2)
      nproj = SIZE(frame%ao%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [nproj, nproj]))
      ALLOCATE (dn(n), dg(n, 3), dt(n), dvn(n), dvg(n, 3), dvt(n), points_response(3, n), weights_response(n), &
                coefficients_response(nproj, nao), soft_response(n, nao), gradient_response(n, nao, 3), &
                projected(nproj, nproj), projected_response(nproj, nproj), &
                operator_response(nao, nao), metric_response(nao, nao))
      cell_inverse = frame%ao%cell
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      points_response(:, :) = MATMUL(strain, frame%ao%points)
      weights_response(:) = volume_response*frame%ao%weights
      CALL pack_grid(frame, vdensity, dn)
      CALL pack_grid(frame, vtau, dt)
      CALL pack_grid(frame, vdensity_response, dvn)
      CALL pack_grid(frame, vtau_response, dvt)
      DO axis = 1, 3
         CALL pack_grid(frame, vgradient(axis), dg(:, axis))
         CALL pack_grid(frame, vgradient_response(axis), dvg(:, axis))
      END DO
      IF (ALLOCATED(frame%ao%patch_center)) THEN
         inverse_dh = frame%dh
         CALL invmat(inverse_dh, info)
         CPASSERT(info == 0)
         ALLOCATE (motion(3, SIZE(frame%ao%center_atoms)))
         DO center = 1, SIZE(frame%ao%center_atoms)
            atom = frame%ao%center_atoms(center)
            relative = atom_response(:, atom) - MATMUL(strain, frame%ao%positions(:, atom))
            motion(:, center) = MATMUL(inverse_dh, relative)
         END DO
         np = frame%ao%patch_points
         DO p = 1, np
            center = frame%ao%patch_center(p)
            row = frame%ao%coarse_points + p
            points_response(:, row) = points_response(:, row) + MATMUL(frame%dh, motion(:, center))
            points_response(:, row + np) = points_response(:, row)
         END DO
         IF (ANY(motion /= 0.0_dp)) THEN
            CALL transfer_patch(frame, vdensity, gathered=dvn, motion=motion)
            CALL transfer_patch(frame, vtau, gathered=dvt, motion=motion)
            DO axis = 1, 3
               CALL transfer_patch(frame, vgradient(axis), gathered=dvg(:, axis), motion=motion)
            END DO
         END IF
      END IF
      IF (ANY(dg /= 0.0_dp) .OR. ANY(dt /= 0.0_dp) .OR. ANY(dvg /= 0.0_dp) .OR. ANY(dvt /= 0.0_dp)) THEN
         CALL qs_upf_inverse_response(frame%ao, points_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, coefficients_response, soft_response, gradient_response)
      ELSE
         ! Scalar potentials and the metric require field responses, but no spatial Hessians.
         gradient_response(:, :, :) = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL qs_upf_inverse_response(frame%ao, points_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, coefficients_response, soft_response)
      END IF
      IF (frame%native_gaussian) THEN
         operator_response = 0.0_dp
      ELSE
         CALL upf_inverse_semilocal_response(frame%ao%soft, frame%ao%soft, &
                                             frame%ao%soft_gradient, frame%ao%soft_gradient, frame%ao%weights, dn, dg, dt, &
                                             soft_response, soft_response, gradient_response, gradient_response, &
                                             weights_response, dvn, dvg, dvt, operator_response, frame%ao%point_comm)
      END IF
      CALL qs_upf_inverse_augmentation_fourier_integrate(frame%ao, augmentation, vdensity, vgradient, projected)
      CALL qs_upf_inverse_augmentation_fourier_response(frame%ao, augmentation, vdensity, vgradient, atom_response, &
                                                        cell_response, vdensity_response, vgradient_response, &
                                                        projected_response)
      IF (PRESENT(projector_operator)) THEN
         CPASSERT(ALL(SHAPE(projector_operator) == [nproj, nproj]))
         projected(:, :) = projected + projector_operator
      END IF
      operator_response(:, :) = operator_response + MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), &
                                                           MATMUL(projected_response, frame%ao%coefficients))
      operator_response(:, :) = operator_response + MATMUL(CONJG(TRANSPOSE(coefficients_response)), &
                                                           MATMUL(projected, frame%ao%coefficients)) + &
                                MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projected, coefficients_response))
      IF (frame%native_gaussian) THEN
         metric_response = 0.0_dp
      ELSE
         IF (ALLOCATED(frame%ao%metric_weights)) THEN
            CALL upf_inverse_pair_response(frame%ao%soft, frame%ao%soft, frame%ao%metric_weights, soft_response, soft_response, &
                                           volume_response*frame%ao%metric_weights, metric_response, frame%ao%point_comm)
         ELSE
            CALL upf_inverse_pair_response(frame%ao%soft, frame%ao%soft, frame%ao%weights, soft_response, soft_response, &
                                           weights_response, metric_response, frame%ao%point_comm)
         END IF
      END IF
      metric_response(:, :) = metric_response + MATMUL(CONJG(TRANSPOSE(coefficients_response)), &
                                                       MATMUL(projector_metric, frame%ao%coefficients)) + &
                              MATMUL(CONJG(TRANSPOSE(frame%ao%coefficients)), MATMUL(projector_metric, coefficients_response))
      IF (PRESENT(dw)) CALL MOVE_ALLOC(coefficients_response, dw)
   END SUBROUTINE qs_upf_inverse_grid_response

! **************************************************************************************************
!> \brief Gather native AO components collectively, using native desymmetrization for S/A signs.
!>        All matrix-communicator ranks participate; the dense result is replicated.
!> \param frame Native atom/AO order used to validate both components.
!> \param matrix_real Native real AO component.
!> \param dense New owned complex dense matrix, real + i*imaginary.
!> \param matrix_imag Optional native imaginary AO component.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_matrix_gather(frame, matrix_real, dense, matrix_imag)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_real
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: dense(:, :)
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix_imag

      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: component(:, :)

      n = frame%ao_last(SIZE(frame%ao_last))
      ALLOCATE (dense(n, n), component(n, n))
      CALL gather_component(frame, matrix_real, component)
      dense(:, :) = CMPLX(component, 0.0_dp, dp)
      IF (PRESENT(matrix_imag)) THEN
         CALL gather_component(frame, matrix_imag, component)
         dense(:, :) = dense + CMPLX(0.0_dp, component, dp)
      END IF
   END SUBROUTINE qs_upf_inverse_matrix_gather

! **************************************************************************************************
!> \brief Store every nonzero AO block on its native owner, including fill absent from the template.
!> \param frame ...
!> \param dense Hermitian complex operator in the native AO order.
!> \param template ...
!> \param matrix_real New owned symmetric real component.
!> \param matrix_imag New owned antisymmetric imaginary component.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_matrix_scatter(frame, dense, template, matrix_real, matrix_imag)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_real, matrix_imag

      INTEGER                                            :: c, i, j, nc, nr, owner, r
      REAL(dp)                                           :: tolerance
      TYPE(mp_comm_type)                                 :: comm

      CALL check_matrix(frame, template)
      nr = frame%ao_last(SIZE(frame%ao_last))
      CPASSERT(ALL(SHAPE(dense) == [nr, nr]))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(dense)))
      CPASSERT(MAXVAL(ABS(dense - CONJG(TRANSPOSE(dense)))) <= tolerance)
      CALL dbcsr_get_info(template, group=comm)
      CALL dbcsr_create(matrix_real, name="Inverse UPF real", template=template, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(matrix_imag, name="Inverse UPF imaginary", template=template, matrix_type=dbcsr_type_antisymmetric)
      DO j = 1, SIZE(frame%ao_first)
         c = frame%ao_first(j)
         nc = frame%ao_last(j) - c + 1
         DO i = 1, j
            r = frame%ao_first(i)
            nr = frame%ao_last(i) - r + 1
            CALL dbcsr_get_stored_coordinates(matrix_real, i, j, owner)
            IF (owner /= comm%mepos) CYCLE
            IF (ANY(REAL(dense(r:r + nr - 1, c:c + nc - 1), dp) /= 0.0_dp)) THEN
               CALL dbcsr_put_block(matrix_real, i, j, REAL(dense(r:r + nr - 1, c:c + nc - 1), dp))
            END IF
            IF (ANY(AIMAG(dense(r:r + nr - 1, c:c + nc - 1)) /= 0.0_dp)) THEN
               CALL dbcsr_put_block(matrix_imag, i, j, AIMAG(dense(r:r + nr - 1, c:c + nc - 1)))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix_real)
      CALL dbcsr_finalize(matrix_imag)
   END SUBROUTINE qs_upf_inverse_matrix_scatter

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param matrix ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE gather_component(frame, matrix, dense)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp), INTENT(OUT)                              :: dense(:, :)

      INTEGER                                            :: c, nc, nr, r
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: full
      TYPE(mp_comm_type)                                 :: comm

      CALL check_matrix(frame, matrix)
      CALL dbcsr_get_info(matrix, group=comm)
      CALL dbcsr_create(full, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix, full)
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iterator, full)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=r, col_offset=c, transposed=transposed)
         nr = SIZE(block, 1)
         nc = SIZE(block, 2)
         IF (transposed) THEN
            dense(r:r + nc - 1, c:c + nr - 1) = TRANSPOSE(block)
         ELSE
            dense(r:r + nr - 1, c:c + nc - 1) = block
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_release(full)
      CALL comm%sum(dense)
   END SUBROUTINE gather_component

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE check_matrix(frame, matrix)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix

      INTEGER, POINTER                                   :: cols(:), rows(:)

      CALL dbcsr_get_info(matrix, row_blk_size=rows, col_blk_size=cols)
      CPASSERT(SIZE(rows) == SIZE(frame%ao_first) .AND. SIZE(cols) == SIZE(frame%ao_first))
      CPASSERT(ALL(rows == frame%ao_last - frame%ao_first + 1))
      CPASSERT(ALL(cols == rows))
   END SUBROUTINE check_matrix

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param field ...
! **************************************************************************************************
   SUBROUTINE check_grid(frame, field)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field

      INTEGER                                            :: extents(3)

      CPASSERT(ASSOCIATED(field%array) .AND. ASSOCIATED(field%pw_grid))
      CPASSERT(ALL(field%pw_grid%bounds == frame%bounds))
      CPASSERT(ALL(field%pw_grid%bounds_local == frame%bounds_local))
      CPASSERT(ALL(field%pw_grid%npts == frame%npts))
      CPASSERT(field%pw_grid%para%mode == frame%mode)
      CPASSERT(ALL(field%pw_grid%dh == frame%dh) .AND. field%pw_grid%dvol == frame%dvol)
      extents = MAX(0, frame%bounds_local(2, :) - frame%bounds_local(1, :) + 1)
      CPASSERT(ALL(SHAPE(field%array) == extents))
      IF (SIZE(field%array) > 0) THEN
         CPASSERT(ALL(LBOUND(field%array) == frame%bounds_local(1, :)))
         CPASSERT(ALL(UBOUND(field%array) == frame%bounds_local(2, :)))
      END IF
      IF (frame%mode == PW_MODE_DISTRIBUTED) THEN
         CPASSERT(field%pw_grid%para%group%mp_comm_type == frame%ao%point_comm)
      END IF
   END SUBROUTINE check_grid

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param field ...
!> \param points ...
! **************************************************************************************************
   SUBROUTINE pack_grid(frame, field, points)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      REAL(dp), INTENT(OUT)                              :: points(:)

      INTEGER                                            :: n

      CALL check_grid(frame, field)
      n = SIZE(field%array)
      points(:n) = RESHAPE(field%array, [n])
      points(n + 1:) = 0.0_dp
      IF (ALLOCATED(frame%ao%patch_center)) CALL transfer_patch(frame, field, gathered=points)
   END SUBROUTINE pack_grid

! **************************************************************************************************
!> \brief ...
!> \param frame ...
!> \param points ...
!> \param field ...
! **************************************************************************************************
   SUBROUTINE unpack_grid(frame, points, field)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      REAL(dp), INTENT(IN)                               :: points(:)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: field

      INTEGER                                            :: n

      CALL check_grid(frame, field)
      n = SIZE(field%array)
      field%array = RESHAPE(points(:n), SHAPE(field%array))
      IF (ALLOCATED(frame%ao%patch_center)) CALL transfer_patch(frame, field, scattered=points)
   END SUBROUTINE unpack_grid

! **************************************************************************************************
!> \brief Add the interpolation or its transpose, including only grid nodes owned by this rank.
!> \param frame Atom-centered quadrature and interpolation weights.
!> \param field Global-grid field. Scattering accumulates in its local array.
!> \param gathered Local samples to which interpolated values are added.
!> \param scattered Paired samples to deposit with opposite quadrature weights.
!> \param motion Optional motion of each center in grid coordinates. Gather the interpolation derivative.
! **************************************************************************************************
   SUBROUTINE transfer_patch(frame, field, gathered, scattered, motion)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: gathered(:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: scattered(:), motion(:, :)

      INTEGER                                            :: a, axis, b, c, center, i, j, k, n, np, &
                                                            offset, p
      REAL(dp)                                           :: dw(interpolation_order, 3), v, &
                                                            w(interpolation_order, 3), &
                                                            wx(interpolation_order), yz

      CPASSERT(PRESENT(gathered) .NEQV. PRESENT(scattered))
      IF (PRESENT(gathered)) THEN
         IF (ALL(field%array == 0.0_dp)) RETURN
      END IF
      n = frame%ao%coarse_points
      np = frame%ao%patch_points
      DO p = 1, np
         center = frame%ao%patch_center(p)
         IF (PRESENT(motion)) THEN
            IF (ALL(motion(:, center) == 0.0_dp)) CYCLE
         END IF
         v = 0.0_dp
         IF (PRESENT(scattered)) THEN
            v = (scattered(n + p) - scattered(n + np + p))/REAL(subdivisions**3, dp)
            IF (v == 0.0_dp) CYCLE
         END IF
         DO axis = 1, 3
            offset = frame%patch_offset(axis, p)
            w(:, axis) = frame%stencil(:, offset, axis, center)
            IF (PRESENT(motion)) dw(:, axis) = frame%dstencil(:, offset, axis, center)*motion(axis, center)
         END DO
         DO k = 1, interpolation_order
            c = frame%patch_index(k, 3, p)
            IF (c == HUGE(0)) CYCLE
            DO j = 1, interpolation_order
               b = frame%patch_index(j, 2, p)
               IF (b == HUGE(0)) CYCLE
               yz = w(j, 2)*w(k, 3)
               IF (PRESENT(motion)) THEN
                  wx = dw(:, 1)*yz + w(:, 1)*(dw(j, 2)*w(k, 3) + w(j, 2)*dw(k, 3))
               ELSE
                  wx = w(:, 1)*yz
               END IF
               IF (PRESENT(scattered)) THEN
                  DO i = 1, interpolation_order
                     a = frame%patch_index(i, 1, p)
                     IF (a /= HUGE(0)) field%array(a, b, c) = field%array(a, b, c) + v*wx(i)
                  END DO
               ELSE
                  DO i = 1, interpolation_order
                     a = frame%patch_index(i, 1, p)
                     IF (a /= HUGE(0)) v = v + field%array(a, b, c)*wx(i)
                  END DO
               END IF
            END DO
         END DO
         IF (PRESENT(gathered)) THEN
            gathered(n + p) = gathered(n + p) + v
            gathered(n + np + p) = gathered(n + np + p) + v
         END IF
      END DO
   END SUBROUTINE transfer_patch

! **************************************************************************************************
!> \brief Release owned inverse fields before their borrowed KS basis owner is destroyed.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_grid_release(frame)
      TYPE(qs_upf_inverse_grid_type), INTENT(INOUT)      :: frame

      TYPE(qs_upf_inverse_grid_type)                     :: empty

      CALL qs_upf_inverse_release(frame%ao)
      frame = empty
   END SUBROUTINE qs_upf_inverse_grid_release

END MODULE qs_upf_inverse_grid

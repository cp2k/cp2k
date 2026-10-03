!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Direct reciprocal-space reference for scalar UPF augmentation.
!>        Uses the complete tabulated radial mesh and CP2K real spherical harmonics.
!>        The overlap correction is the zero-frequency limit of this same transform.
! **************************************************************************************************
MODULE uspp_augmentation
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE bessel_lib,                      ONLY: bessel0
   USE cp_array_sort,                   ONLY: cp_1d_r_sort
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE spherical_harmonics,             ONLY: y_lm
   USE splines,                         ONLY: spline3ders
   USE upf_projector_integrals,         ONLY: upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uspp_augmentation'

   TYPE uspp_augmentation_type
      INTEGER :: nproj = 0, lmax = -1
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:), radial_end(:, :)
      REAL(dp) :: realspace_radius = 0.0_dp
      ! weighted_q can omit an exact zero suffix; r retains the full mesh for real-space splines.
      REAL(KIND=dp), ALLOCATABLE :: r(:), weighted_q(:, :, :), gaunt(:, :, :), q(:, :)
      REAL(KIND=dp), ALLOCATABLE :: radial_x(:), radial_value(:, :, :), radial_slope(:, :, :)
      INTEGER(int_8) :: fourier_budget = -1_int_8
      INTEGER, ALLOCATABLE :: fourier_map(:)
      REAL(dp), ALLOCATABLE :: fourier_norm(:), fourier_radial(:, :, :), fourier_next(:, :, :)
   END TYPE uspp_augmentation_type

   PUBLIC :: uspp_augmentation_type, uspp_augmentation_init, uspp_augmentation_fourier, &
             uspp_augmentation_collocate, uspp_augmentation_integrate, &
             uspp_augmentation_prepare_realspace, uspp_augmentation_realspace, &
             uspp_augmentation_prepare_fourier, &
             uspp_augmentation_collocate_atom, uspp_augmentation_integrate_atom

CONTAINS

! **************************************************************************************************
!> \brief Prepare radial quadrature and Gaunt coefficients in the native beta convention.
!> \param pot Original scalar USPP with angular-resolved augmentation.
!> \param projectors Prepared native projectors, with Cartesian degree at least max(lbeta).
!> \param augmentation Immutable augmentation data, in the same magnetic-projector ordering.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_init(pot, projectors, augmentation)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(uspp_augmentation_type), INTENT(OUT)          :: augmentation

      INTEGER                                            :: i, j, l, m, n, nc, p, pair, q, s
      REAL(KIND=dp), ALLOCATABLE                         :: weights(:)

      CPASSERT(.NOT. pot%has_so .AND. .NOT. pot%is_paw)
      CPASSERT(ALLOCATED(pot%qijl))
      CPASSERT(projectors%degree >= MAXVAL(pot%lbeta))
      CPASSERT(ALL(projectors%lbeta == pot%lbeta))
      n = SIZE(pot%r)
      CPASSERT(n >= 3 .AND. SIZE(pot%rab) == n)
      CPASSERT(ALL(pot%r >= 0.0_dp) .AND. ALL(pot%rab > 0.0_dp))
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)))
      nc = pot%number_of_proj*(pot%number_of_proj + 1)/2
      augmentation%nproj = projectors%nproj
      augmentation%lmax = 2*MAXVAL(pot%lbeta)
      CALL init_orbital_pointers(augmentation%lmax)
      CALL init_spherical_harmonics(augmentation%lmax, -1)
      CPASSERT(ALL(SHAPE(pot%qijl) == [n, nc, augmentation%lmax + 1]))
      CPASSERT(LBOUND(pot%qijl, 3) == 0)
      ALLOCATE (augmentation%r(n), weights(n))
      augmentation%r(:) = pot%r
      weights(:) = upf_radial_weights(pot)
      ! Keep the quadrature weights of the original mesh. Only an exact zero
      ! suffix can be omitted; shortening the mesh before forming its weights
      ! would change the endpoint rule and hence the augmentation metric.
      DO WHILE (n > 0)
         IF (ANY(pot%qijl(n, :, :) /= 0.0_dp)) EXIT
         n = n - 1
      END DO
      ALLOCATE (augmentation%weighted_q(n, nc, 0:augmentation%lmax))
      DO l = 0, augmentation%lmax
         DO pair = 1, nc
            augmentation%weighted_q(:, pair, l) = weights(:n)*pot%qijl(:n, pair, l)
         END DO
      END DO
      n = augmentation%nproj
      ALLOCATE (augmentation%radial_index(n), augmentation%angular_index(n), &
                augmentation%lbeta(pot%number_of_proj))
      augmentation%radial_index(:) = projectors%radial_index
      augmentation%angular_index(:) = projectors%angular_index
      augmentation%lbeta(:) = pot%lbeta
      ALLOCATE (augmentation%gaunt(n, n, (augmentation%lmax + 1)**2), augmentation%q(n, n))
      augmentation%gaunt = 0.0_dp
      augmentation%q = 0.0_dp
      DO q = 1, n
         j = projectors%radial_index(q)
         l = pot%lbeta(j)
         m = projectors%angular_index(q) - l*l
         DO p = 1, n
            i = projectors%radial_index(p)
            ! Contract the existing analytic monomial/Y/Y integrals with Y_j's polynomial.
            DO s = 1, nco(l)
               augmentation%gaunt(p, q, :) = augmentation%gaunt(p, q, :) + &
                                             orbtramat(l)%slm(m, s)*projectors%gaunt(ncoset(l - 1) + s, &
                                                                           projectors%angular_index(p), :(augmentation%lmax + 1)**2)
            END DO
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            IF (projectors%angular_index(p) == projectors%angular_index(q)) THEN
               augmentation%q(p, q) = SUM(augmentation%weighted_q(:, pair, 0))
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_init

! **************************************************************************************************
!> \brief Prepare regular radial splines for atom-grid augmentation quadrature.
!>        UPF stores r**2 Q_ij^L(r). Interpolate Q_ij^L(r)/r**L in r**2 and multiply
!>        by solid harmonics, so values and Cartesian derivatives remain regular at the origin.
!>        The existing reciprocal-space quadrature and its overlap correction are unchanged.
!> \param pot Original scalar USPP used to initialize augmentation.
!> \param augmentation Prepared angular data; owns the optional real-space cache on return.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_prepare_realspace(pot, augmentation)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(uspp_augmentation_type), INTENT(INOUT)        :: augmentation

      INTEGER                                            :: first, l, last, n, nc, pair
      REAL(KIND=dp), ALLOCATABLE                         :: values(:), x(:)

      CPASSERT(ALLOCATED(pot%qijl) .AND. ALLOCATED(augmentation%r))
      CPASSERT(SIZE(pot%r) == SIZE(augmentation%r))
      CPASSERT(ALL(pot%r == augmentation%r))
      CPASSERT(ALL(pot%lbeta == augmentation%lbeta))
      IF (ALLOCATED(augmentation%radial_x)) RETURN
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      n = SIZE(pot%r) - first + 1
      CPASSERT(n >= 4)
      nc = SIZE(pot%qijl, 2)
      ALLOCATE (x(n), values(n), augmentation%radial_x(n + 1), &
                augmentation%radial_value(n + 1, nc, 0:augmentation%lmax), &
                augmentation%radial_slope(n + 1, nc, 0:augmentation%lmax), &
                augmentation%radial_end(nc, 0:augmentation%lmax))
      x(:) = pot%r(first:)**2
      augmentation%radial_x(1) = 0.0_dp
      augmentation%radial_x(2:) = x
      augmentation%radial_value = 0.0_dp
      augmentation%radial_slope = 0.0_dp
      augmentation%radial_end = 1
      augmentation%realspace_radius = 0.0_dp
      DO l = 0, augmentation%lmax
         DO pair = 1, nc
            IF (ALL(pot%qijl(:, pair, l) == 0.0_dp)) CYCLE
            ! Preserve the exact zero suffix of each supplied augmentation channel.
            ! No numerical threshold is used, and a nonzero last node keeps the full mesh.
            last = SIZE(pot%r)
            DO WHILE (pot%qijl(last, pair, l) == 0.0_dp)
               last = last - 1
            END DO
            augmentation%radial_end(pair, l) = last - first + 3
            augmentation%realspace_radius = MAX(augmentation%realspace_radius, pot%r(MIN(last + 1, SIZE(pot%r))))
            values(:) = pot%qijl(first:, pair, l)/pot%r(first:)**(l + 2)
            ! The first interval continues the endpoint cubic to x=0. A tabulated
            ! r**2 Q value at r=0 alone does not determine the regular radial limit.
            CALL spline3ders(x, values, augmentation%radial_x, &
                             ynew=augmentation%radial_value(:, pair, l), &
                             dynew=augmentation%radial_slope(:, pair, l))
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_prepare_realspace

! **************************************************************************************************
!> \brief Evaluate the augmentation matrix and derivatives at a Cartesian displacement.
!>        Uses cubic splines with a C2 quintic closure into an exact tabulated zero suffix.
!>        Values, gradients and Hessians share the same polynomial. No angular division
!>        by radius is needed, including at the center and on coordinate axes.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacement Evaluation point minus augmentation center, in bohr.
!> \param qr Symmetric Q_ij(r), without radial Jacobian or quadrature weights.
!> \param gradient Optional Cartesian gradient of qr with respect to the evaluation point.
!> \param hessian Optional Cartesian Hessian of qr, also needed for GGA geometry derivatives.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_realspace(augmentation, displacement, qr, gradient, hessian)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacement(3)
      REAL(KIND=dp), INTENT(OUT)                         :: qr(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: gradient(:, :, :), hessian(:, :, :, :)

      REAL(dp), ALLOCATABLE                              :: ddharmonics(:, :, :), dharmonics(:, :), &
                                                            harmonics(:)

      ALLOCATE (harmonics((augmentation%lmax + 1)**2), &
                dharmonics(3, (augmentation%lmax + 1)**2), &
                ddharmonics(3, 3, (augmentation%lmax + 1)**2))
      CALL realspace_with_workspace(augmentation, displacement, qr, gradient, hessian, &
                                    harmonics, dharmonics, ddharmonics)
   END SUBROUTINE uspp_augmentation_realspace

! **************************************************************************************************
!> \brief Point kernel with reusable angular scratch space for atom-grid batches.
!> \param augmentation Initialized radial and angular data
!> \param displacement Point relative to augmentation center
!> \param qr Augmentation matrix
!> \param gradient Optional Cartesian gradient
!> \param hessian Optional Cartesian Hessian
!> \param harmonics Scratch values of solid harmonics
!> \param dharmonics Scratch gradients of solid harmonics
!> \param ddharmonics Scratch Hessians of solid harmonics
! **************************************************************************************************
   SUBROUTINE realspace_with_workspace(augmentation, displacement, qr, gradient, hessian, &
                                       harmonics, dharmonics, ddharmonics)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacement(3)
      REAL(KIND=dp), INTENT(OUT)                         :: qr(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: gradient(:, :, :), hessian(:, :, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: harmonics(:), dharmonics(:, :), &
                                                            ddharmonics(:, :, :)

      INTEGER                                            :: a, axis, axis2, hi, i, j, l, lm, lo, m, &
                                                            mid, p, pair, powers(3), q, reduced(3)
      REAL(KIND=dp) :: angular, angular_gradient(3), angular_hessian(3, 3), c0, c1, c2, c3, c4, &
         c5, factor, h, radial, radial_d, radial_dd, s, t, term, value, value_gradient(3), &
         value_hessian(3, 3)

      CPASSERT(ALLOCATED(augmentation%radial_x))
      CPASSERT(ALL(SHAPE(qr) == augmentation%nproj))
      qr = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [augmentation%nproj, augmentation%nproj, 3]))
         gradient = 0.0_dp
      END IF
      IF (PRESENT(hessian)) THEN
         CPASSERT(ALL(SHAPE(hessian) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         hessian = 0.0_dp
      END IF
      s = SUM(displacement**2)
      hi = SIZE(augmentation%radial_x)
      IF (s > augmentation%realspace_radius**2) RETURN
      lo = 1
      DO WHILE (hi - lo > 1)
         mid = (lo + hi)/2
         IF (s < augmentation%radial_x(mid)) THEN
            hi = mid
         ELSE
            lo = mid
         END IF
      END DO
      h = augmentation%radial_x(hi) - augmentation%radial_x(lo)
      t = (s - augmentation%radial_x(lo))/h
      harmonics = 0.0_dp
      dharmonics = 0.0_dp
      ddharmonics = 0.0_dp
      DO l = 0, augmentation%lmax
         DO m = -l, l
            lm = l*l + l + m + 1
            DO a = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + a)
               factor = orbtramat(l)%slm(l + m + 1, a)
               harmonics(lm) = harmonics(lm) + factor*PRODUCT(displacement**powers)
               DO axis = 1, 3
                  IF (powers(axis) == 0) CYCLE
                  reduced = powers
                  reduced(axis) = reduced(axis) - 1
                  term = factor*REAL(powers(axis), dp)
                  dharmonics(axis, lm) = dharmonics(axis, lm) + term*PRODUCT(displacement**reduced)
                  DO axis2 = 1, 3
                     IF (reduced(axis2) == 0) CYCLE
                     reduced(axis2) = reduced(axis2) - 1
                     ddharmonics(axis, axis2, lm) = ddharmonics(axis, axis2, lm) + &
                                                    term*REAL(reduced(axis2) + 1, dp)*PRODUCT(displacement**reduced)
                     reduced(axis2) = reduced(axis2) + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      DO q = 1, augmentation%nproj
         j = augmentation%radial_index(q)
         DO p = 1, q
            i = augmentation%radial_index(p)
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
               augmentation%lbeta(i) + augmentation%lbeta(j), 2
               IF (lo >= augmentation%radial_end(pair, l)) CYCLE
               c0 = augmentation%radial_value(lo, pair, l)
               c1 = h*augmentation%radial_slope(lo, pair, l)
               c2 = 3.0_dp*(augmentation%radial_value(hi, pair, l) - c0) - &
                    2.0_dp*c1 - h*augmentation%radial_slope(hi, pair, l)
               c3 = 2.0_dp*(c0 - augmentation%radial_value(hi, pair, l)) + &
                    c1 + h*augmentation%radial_slope(hi, pair, l)
               c4 = 0.0_dp
               c5 = 0.0_dp
               IF (hi == augmentation%radial_end(pair, l)) THEN
                  ! Match value and first/second x derivatives at the last nonzero node;
                  ! all three vanish at the first node of the exact zero suffix.
                  c3 = -10.0_dp*c0 - 6.0_dp*c1 - 3.0_dp*c2
                  c4 = 15.0_dp*c0 + 8.0_dp*c1 + 3.0_dp*c2
                  c5 = -6.0_dp*c0 - 3.0_dp*c1 - c2
               END IF
               radial = c0 + t*(c1 + t*(c2 + t*(c3 + t*(c4 + t*c5))))
               radial_d = (c1 + t*(2.0_dp*c2 + t*(3.0_dp*c3 + t*(4.0_dp*c4 + t*5.0_dp*c5))))/h
               radial_dd = (2.0_dp*c2 + t*(6.0_dp*c3 + t*(12.0_dp*c4 + t*20.0_dp*c5)))/h**2
               angular = 0.0_dp
               angular_gradient = 0.0_dp
               angular_hessian = 0.0_dp
               DO lm = l*l + 1, (l + 1)**2
                  factor = augmentation%gaunt(p, q, lm)
                  angular = angular + factor*harmonics(lm)
                  angular_gradient = angular_gradient + factor*dharmonics(:, lm)
                  angular_hessian = angular_hessian + factor*ddharmonics(:, :, lm)
               END DO
               value = radial*angular
               qr(p, q) = qr(p, q) + value
               IF (PRESENT(gradient)) THEN
                  value_gradient = 2.0_dp*displacement*radial_d*angular + radial*angular_gradient
                  gradient(p, q, :) = gradient(p, q, :) + value_gradient
               END IF
               IF (PRESENT(hessian)) THEN
                  DO axis2 = 1, 3
                     DO axis = 1, 3
                        value_hessian(axis, axis2) = &
                           4.0_dp*displacement(axis)*displacement(axis2)*radial_dd*angular + &
                           2.0_dp*radial_d*(displacement(axis)*angular_gradient(axis2) + &
                                            displacement(axis2)*angular_gradient(axis)) + &
                           radial*angular_hessian(axis, axis2)
                     END DO
                     value_hessian(axis2, axis2) = value_hessian(axis2, axis2) + 2.0_dp*radial_d*angular
                  END DO
                  hessian(p, q, :, :) = hessian(p, q, :, :) + value_hessian
               END IF
            END DO
            qr(q, p) = qr(p, q)
            IF (PRESENT(gradient)) gradient(q, p, :) = gradient(p, q, :)
            IF (PRESENT(hessian)) hessian(q, p, :, :) = hessian(p, q, :, :)
         END DO
      END DO
   END SUBROUTINE realspace_with_workspace

! **************************************************************************************************
!> \brief Add one center's augmentation density and gradient to arbitrary atom-grid points.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacements Cartesian quadrature points relative to this center, in bohr.
!> \param density Real symmetric on-site density matrix, including off-diagonal entries.
!> \param rho Accumulated density values; no quadrature weights included.
!> \param gradient Accumulated Cartesian density gradients.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_collocate_atom(augmentation, displacements, density, rho, gradient)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacements(:, :), density(:, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: rho(:), gradient(:, :)

      INTEGER                                            :: axis, n, point
      REAL(dp), ALLOCATABLE                              :: ddharmonics(:, :, :), dharmonics(:, :), &
                                                            harmonics(:)
      REAL(KIND=dp), ALLOCATABLE                         :: dqr(:, :, :), qr(:, :)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(density) == n))
      CPASSERT(ALL(SHAPE(displacements) == [3, SIZE(rho)]))
      CPASSERT(ALL(SHAPE(gradient) == SHAPE(displacements)))
      ALLOCATE (qr(n, n), dqr(n, n, 3))
      ALLOCATE (harmonics((augmentation%lmax + 1)**2), &
                dharmonics(3, (augmentation%lmax + 1)**2), &
                ddharmonics(3, 3, (augmentation%lmax + 1)**2))
      DO point = 1, SIZE(rho)
         CALL realspace_with_workspace(augmentation, displacements(:, point), qr, dqr, &
                                       harmonics=harmonics, dharmonics=dharmonics, ddharmonics=ddharmonics)
         rho(point) = rho(point) + SUM(density*qr)
         DO axis = 1, 3
            gradient(axis, point) = gradient(axis, point) + SUM(density*dqr(:, :, axis))
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_collocate_atom

! **************************************************************************************************
!> \brief Adjoint of the atom-grid density map, including the GGA gradient term.
!> \param augmentation Data with an initialized real-space cache.
!> \param displacements Cartesian quadrature points relative to this center, in bohr.
!> \param weights Quadrature weights, optionally signed for a hard/soft energy difference.
!> \param v_rho Derivative of the energy density with respect to density, without weights.
!> \param v_gradient Derivative with respect to the Cartesian density gradient, without weights.
!> \param matrix Resulting symmetric on-site operator.
!> \param density Optional on-site density, required for center_gradient.
!> \param center_gradient Optional energy gradient from translating this center at fixed grid points.
!>        Projector, moving-grid and partition-weight derivatives are separate contributions.
!> \param point_gradient Optional energy gradient of each quadrature point at fixed center and weights.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_integrate_atom(augmentation, displacements, weights, v_rho, v_gradient, &
                                               matrix, density, center_gradient, point_gradient)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: displacements(:, :), weights(:), &
                                                            v_rho(:), v_gradient(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: matrix(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: density(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: center_gradient(3), point_gradient(:, :)

      INTEGER                                            :: axis, axis2, n, point
      REAL(dp), ALLOCATABLE                              :: ddharmonics(:, :, :), dharmonics(:, :), &
                                                            harmonics(:)
      REAL(KIND=dp)                                      :: derivative
      REAL(KIND=dp), ALLOCATABLE                         :: ddqr(:, :, :, :), dqr(:, :, :), qr(:, :)

      n = augmentation%nproj
      CPASSERT(ALL(SHAPE(matrix) == n))
      CPASSERT(ALL(SHAPE(displacements) == [3, SIZE(weights)]))
      CPASSERT(SIZE(v_rho) == SIZE(weights))
      CPASSERT(ALL(SHAPE(v_gradient) == SHAPE(displacements)))
      ALLOCATE (qr(n, n), dqr(n, n, 3))
      ALLOCATE (harmonics((augmentation%lmax + 1)**2), &
                dharmonics(3, (augmentation%lmax + 1)**2), &
                ddharmonics(3, 3, (augmentation%lmax + 1)**2))
      matrix = 0.0_dp
      IF (PRESENT(center_gradient) .OR. PRESENT(point_gradient)) THEN
         CPASSERT(PRESENT(density))
         CPASSERT(ALL(SHAPE(density) == n))
         ALLOCATE (ddqr(n, n, 3, 3))
      END IF
      IF (PRESENT(center_gradient)) center_gradient = 0.0_dp
      IF (PRESENT(point_gradient)) THEN
         CPASSERT(ALL(SHAPE(point_gradient) == SHAPE(displacements)))
         point_gradient = 0.0_dp
      END IF
      DO point = 1, SIZE(weights)
         IF (weights(point) == 0.0_dp) CYCLE
         IF (v_rho(point) == 0.0_dp .AND. ALL(v_gradient(:, point) == 0.0_dp)) CYCLE
         IF (ALLOCATED(ddqr)) THEN
            CALL realspace_with_workspace(augmentation, displacements(:, point), qr, dqr, ddqr, &
                                          harmonics, dharmonics, ddharmonics)
         ELSE
            CALL realspace_with_workspace(augmentation, displacements(:, point), qr, dqr, &
                                          harmonics=harmonics, dharmonics=dharmonics, ddharmonics=ddharmonics)
         END IF
         matrix(:, :) = matrix + weights(point)*v_rho(point)*qr
         DO axis = 1, 3
            matrix(:, :) = matrix + weights(point)*v_gradient(axis, point)*dqr(:, :, axis)
         END DO
         IF (ALLOCATED(ddqr)) THEN
            DO axis = 1, 3
               derivative = v_rho(point)*SUM(density*dqr(:, :, axis))
               DO axis2 = 1, 3
                  derivative = derivative + v_gradient(axis2, point)*SUM(density*ddqr(:, :, axis2, axis))
               END DO
               IF (PRESENT(center_gradient)) center_gradient(axis) = center_gradient(axis) - weights(point)*derivative
               IF (PRESENT(point_gradient)) point_gradient(axis, point) = weights(point)*derivative
            END DO
         END IF
      END DO
   END SUBROUTINE uspp_augmentation_integrate_atom

! **************************************************************************************************
!> \brief Cache exact radial transforms for repeated reciprocal-vector lengths.
!>        The retained arrays are bounded; uncached lengths use direct quadrature.
!>        Initialization of a new augmentation object also releases this cache.
!> \param augmentation Prepared radial data, with the cache updated in place.
!> \param gvecs Local Cartesian reciprocal vectors, in their current ordering.
!> \param max_bytes Optional retained-memory limit per kind and rank (default 64 MiB).
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_prepare_fourier(augmentation, gvecs, max_bytes)
      TYPE(uspp_augmentation_type), INTENT(INOUT)        :: augmentation
      REAL(dp), INTENT(IN)                               :: gvecs(:, :)
      INTEGER(int_8), INTENT(IN), OPTIONAL               :: max_bytes

      INTEGER                                            :: capacity, ig, index, lmax, nc, ng, unique
      INTEGER(int_8)                                     :: budget, entry_bytes, index_bytes
      INTEGER, ALLOCATABLE                               :: order(:)
      REAL(dp)                                           :: previous
      REAL(dp), ALLOCATABLE                              :: norms(:), sorted(:)

      CPASSERT(SIZE(gvecs, 1) == 3)
      budget = 64_int_8*1024*1024
      IF (PRESENT(max_bytes)) budget = max_bytes
      CPASSERT(budget >= 0_int_8)
      ng = SIZE(gvecs, 2)
      nc = SIZE(augmentation%weighted_q, 2)
      lmax = augmentation%lmax
      CPASSERT(nc > 0 .AND. lmax >= 0)
      index_bytes = INT(ng, int_8)*(STORAGE_SIZE(0.0_dp) + STORAGE_SIZE(0))/8
      entry_bytes = 2_int_8*nc*(lmax + 1)*STORAGE_SIZE(0.0_dp)/8
      capacity = INT(MIN(INT(ng, int_8), MAX(0_int_8, (budget - index_bytes)/entry_bytes)))
      IF (capacity > 0) THEN
         ALLOCATE (norms(ng))
         DO ig = 1, ng
            norms(ig) = SQRT(SUM(gvecs(:, ig)*gvecs(:, ig)))
         END DO
         IF (ALLOCATED(augmentation%fourier_norm) .AND. augmentation%fourier_budget == budget) THEN
            IF (SIZE(augmentation%fourier_norm) == ng) THEN
               IF (ALL(augmentation%fourier_norm == norms)) RETURN
            END IF
         END IF
      END IF
      IF (ALLOCATED(augmentation%fourier_norm)) THEN
         DEALLOCATE (augmentation%fourier_norm, augmentation%fourier_map, &
                     augmentation%fourier_radial, augmentation%fourier_next)
      END IF
      augmentation%fourier_budget = budget
      IF (capacity == 0) RETURN
      ALLOCATE (augmentation%fourier_norm, SOURCE=norms)
      ALLOCATE (augmentation%fourier_map(ng), order(ng), sorted(ng))
      augmentation%fourier_map = 0
      sorted(:) = norms
      CALL cp_1d_r_sort(sorted, ng, order)
      unique = 0
      previous = -1.0_dp
      DO ig = 1, ng
         IF (sorted(ig) == 0.0_dp) CYCLE
         IF (sorted(ig) /= previous) unique = unique + 1
         previous = sorted(ig)
      END DO
      capacity = MIN(capacity, unique)
      ALLOCATE (augmentation%fourier_radial(nc, 0:lmax, capacity), &
                augmentation%fourier_next(nc, 0:lmax, capacity))
      index = 0
      previous = -1.0_dp
      DO ig = 1, ng
         IF (sorted(ig) == 0.0_dp) CYCLE
         IF (sorted(ig) /= previous) THEN
            IF (index == capacity) EXIT
            index = index + 1
            previous = sorted(ig)
            CALL augmentation_radial(augmentation%r, augmentation%weighted_q, previous, &
                                     augmentation%fourier_radial(:, :, index), augmentation%fourier_next(:, :, index))
         END IF
         augmentation%fourier_map(order(ig)) = index
      END DO
   END SUBROUTINE uspp_augmentation_prepare_fourier

! **************************************************************************************************
!> \brief Direct radial quadrature shared by cached and uncached Fourier evaluations.
!> \param r Original radial mesh.
!> \param weighted_q Weighted radial augmentation.
!> \param norm_g Reciprocal-vector length.
!> \param radial Integrals against j_l(Gr).
!> \param radial_next Optional integrals against r*j_(l+1)(Gr), for strain derivatives.
! **************************************************************************************************
   SUBROUTINE augmentation_radial(r, weighted_q, norm_g, radial, radial_next)
      REAL(dp), INTENT(IN)                               :: r(:), weighted_q(:, :, 0:), norm_g
      REAL(dp), INTENT(OUT)                              :: radial(:, 0:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: radial_next(:, 0:)

      INTEGER                                            :: l, nr
      REAL(dp), ALLOCATABLE                              :: work(:)

      nr = SIZE(weighted_q, 1)
      ALLOCATE (work(nr))
      DO l = 0, UBOUND(weighted_q, 3)
         IF (PRESENT(radial_next)) THEN
            work(:) = r(:nr)*bessel0(norm_g*r(:nr), l + 1)
            radial_next(:, l) = MATMUL(work, weighted_q(:, :, l))
         END IF
         work(:) = bessel0(norm_g*r(:nr), l)
         radial(:, l) = MATMUL(work, weighted_q(:, :, l))
      END DO
   END SUBROUTINE augmentation_radial

! **************************************************************************************************
!> \brief Evaluate Q_ij(G) = integral Q_ij(r-R) exp(-i G.r) dr / volume.
!>        Direct radial quadrature is a reference without reciprocal interpolation.
!> \param augmentation Prepared data.
!> \param g Cartesian reciprocal vector in inverse bohr.
!> \param position Augmentation center in bohr.
!> \param volume Cell volume in bohr cubed.
!> \param qg Complex symmetric form-factor matrix. Q(-G) = conjugate(Q(G)).
!> \param strain Optional dQ/d epsilon_ab at fixed fractional center and reciprocal index.
!>        Includes the inverse-volume derivative; the center phase is strain invariant.
!> \param grid_index Optional index in the prepared Fourier cache; absent entries use direct quadrature.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_fourier(augmentation, g, position, volume, qg, strain, grid_index)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: g(3), position(3), volume
      COMPLEX(KIND=dp), INTENT(OUT)                      :: qg(:, :)
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: strain(:, :, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: grid_index

      COMPLEX(KIND=dp)                                   :: phase, phase_l
      INTEGER                                            :: a, axis, cached, i, j, l, lm, m, p, &
                                                            pair, powers(3), q
      REAL(KIND=dp)                                      :: angular, angular_derivative(3), &
                                                            direction(3, 1), monomial, norm_g, &
                                                            value(1)
      REAL(KIND=dp), ALLOCATABLE                         :: harmonics(:), radial(:, :), &
                                                            radial_next(:, :), solid_gradient(:, :)

      CPASSERT(volume > 0.0_dp)
      CPASSERT(ALL(SHAPE(qg) == augmentation%nproj))
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         strain = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      norm_g = SQRT(SUM(g*g))
      IF (norm_g == 0.0_dp) THEN
         qg(:, :) = CMPLX(augmentation%q/volume, 0.0_dp, KIND=dp)
         IF (PRESENT(strain)) THEN
            DO axis = 1, 3
               strain(:, :, axis, axis) = -qg
            END DO
         END IF
         RETURN
      END IF
      direction(:, 1) = g/norm_g
      ALLOCATE (harmonics((augmentation%lmax + 1)**2))
      ALLOCATE (radial(SIZE(augmentation%weighted_q, 2), 0:augmentation%lmax))
      IF (PRESENT(strain)) THEN
         ALLOCATE (solid_gradient(3, (augmentation%lmax + 1)**2), &
                   radial_next(SIZE(radial, 1), 0:augmentation%lmax))
         solid_gradient = 0.0_dp
      END IF
      cached = 0
      IF (PRESENT(grid_index) .AND. ALLOCATED(augmentation%fourier_map)) THEN
         ! Independent callers can supply a changed or reordered list without preparing it.
         ! Reuse only a matching length; the direct path remains valid for every other vector.
         IF (grid_index >= 1 .AND. grid_index <= SIZE(augmentation%fourier_map)) THEN
            IF (augmentation%fourier_norm(grid_index) == norm_g) cached = augmentation%fourier_map(grid_index)
         END IF
      END IF
      IF (cached > 0) THEN
         radial(:, :) = augmentation%fourier_radial(:, :, cached)
         IF (PRESENT(strain)) radial_next(:, :) = augmentation%fourier_next(:, :, cached)
      ELSE
         CALL augmentation_radial(augmentation%r, augmentation%weighted_q, norm_g, radial, radial_next)
      END IF
      DO l = 0, augmentation%lmax
         DO m = -l, l
            lm = l*l + l + m + 1
            CALL y_lm(direction, value, l, m)
            harmonics(lm) = value(1)
            IF (PRESENT(strain)) THEN
               ! Cartesian solid-harmonic derivatives avoid polar-axis singularities.
               DO a = 1, nco(l)
                  DO axis = 1, 3
                     powers = indco(:, ncoset(l - 1) + a)
                     IF (powers(axis) == 0) CYCLE
                     monomial = REAL(powers(axis), dp)
                     powers(axis) = powers(axis) - 1
                     solid_gradient(axis, lm) = solid_gradient(axis, lm) + &
                                                orbtramat(l)%slm(l + m + 1, a)*monomial*PRODUCT(direction(:, 1)**powers)
                  END DO
               END DO
            END IF
         END DO
      END DO
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, position), KIND=dp))*fourpi/volume
      qg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO q = 1, augmentation%nproj
         j = augmentation%radial_index(q)
         DO p = 1, q
            i = augmentation%radial_index(p)
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
               augmentation%lbeta(i) + augmentation%lbeta(j), 2
               angular = DOT_PRODUCT(augmentation%gaunt(p, q, l*l + 1:(l + 1)**2), &
                                     harmonics(l*l + 1:(l + 1)**2))
               phase_l = CMPLX(0.0_dp, -1.0_dp, KIND=dp)**l
               qg(p, q) = qg(p, q) + phase_l*angular*radial(pair, l)
               IF (PRESENT(strain)) THEN
                  angular_derivative = MATMUL(solid_gradient(:, l*l + 1:(l + 1)**2), &
                                              augmentation%gaunt(p, q, l*l + 1:(l + 1)**2))
                  ! d[j_l(Gr)Y_l(G)]/dG combines the l/G terms analytically.
                  angular_derivative = radial(pair, l)/norm_g*angular_derivative &
                                       - radial_next(pair, l)*angular*direction(:, 1)
                  DO axis = 1, 3
                     strain(p, q, :, axis) = strain(p, q, :, axis) - phase*phase_l*g*angular_derivative(axis)
                  END DO
               END IF
            END DO
            qg(p, q) = phase*qg(p, q)
            qg(q, p) = qg(p, q)
            IF (PRESENT(strain)) THEN
               DO axis = 1, 3
                  strain(p, q, axis, axis) = strain(p, q, axis, axis) - qg(p, q)
               END DO
               strain(q, p, :, :) = strain(p, q, :, :)
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_fourier

! **************************************************************************************************
!> \brief Add the augmentation density on a list of reciprocal vectors.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param density_matrix Real symmetric projector density matrix, with occupations included.
!> \param rho Fourier density coefficients (integral divided by volume), accumulated in place.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_collocate(augmentation, gvecs, position, volume, density_matrix, rho)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume, &
                                                            density_matrix(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: qg(:, :)
      INTEGER                                            :: ig

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(rho))
      CPASSERT(ALL(SHAPE(density_matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      DO ig = 1, SIZE(rho)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg, grid_index=ig)
         rho(ig) = rho(ig) + SUM(density_matrix*qg)
      END DO
   END SUBROUTINE uspp_augmentation_collocate

! **************************************************************************************************
!> \brief Adjoint of augmentation collocation: integrate a real local potential against Q.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param vg Fourier potential coefficients (integral divided by volume).
!> \param multiplicity One for a complete reciprocal list; two for a retained +/-G pair.
!> \param matrix Real symmetric screened projector contribution, overwritten.
!> \param derivative Optional derivatives of matrix with respect to the augmentation center.
!> \param strain Optional density strain contribution at fixed Fourier potential and volume weight.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_integrate(augmentation, gvecs, position, volume, vg, multiplicity, matrix, derivative, strain)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume
      COMPLEX(KIND=dp), INTENT(IN)                       :: vg(:)
      REAL(KIND=dp), INTENT(IN)                          :: multiplicity(:)
      REAL(KIND=dp), INTENT(OUT)                         :: matrix(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: derivative(:, :, :), strain(:, :, :, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: dq(:, :, :, :), qg(:, :)
      INTEGER                                            :: idir, ig, jdir

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(vg))
      CPASSERT(SIZE(multiplicity) == SIZE(vg))
      CPASSERT(ALL(SHAPE(matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      matrix = 0.0_dp
      IF (PRESENT(derivative)) THEN
         CPASSERT(ALL(SHAPE(derivative) == [augmentation%nproj, augmentation%nproj, 3]))
         derivative = 0.0_dp
      END IF
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         ALLOCATE (dq(augmentation%nproj, augmentation%nproj, 3, 3))
         strain = 0.0_dp
      END IF
      DO ig = 1, SIZE(vg)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg, dq, grid_index=ig)
         IF (PRESENT(strain)) THEN
            DO jdir = 1, 3
               DO idir = 1, 3
                  strain(:, :, idir, jdir) = strain(:, :, idir, jdir) + volume*multiplicity(ig)* &
                                             REAL(CONJG(dq(:, :, idir, jdir))*vg(ig), KIND=dp)
               END DO
            END DO
         END IF
         matrix(:, :) = matrix + volume*multiplicity(ig)*REAL(CONJG(qg)*vg(ig), KIND=dp)
         IF (PRESENT(derivative)) THEN
            DO idir = 1, 3
               derivative(:, :, idir) = derivative(:, :, idir) + volume*multiplicity(ig)* &
                                        REAL(CMPLX(0.0_dp, gvecs(idir, ig), dp)*CONJG(qg)*vg(ig), KIND=dp)
            END DO
         END IF
      END DO
   END SUBROUTINE uspp_augmentation_integrate

END MODULE uspp_augmentation

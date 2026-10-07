!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fixed periodic atom quadrature with independent energy and descriptor domains.
! **************************************************************************************************
MODULE skala_atom_grid
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_atom_grid'

   TYPE seed_environment_type
      REAL(dp), ALLOCATABLE :: displacement(:, :), distance(:), membership(:)
      INTEGER, ALLOCATABLE :: atom(:), image(:, :)
   END TYPE seed_environment_type

   TYPE skala_partition_type
      PRIVATE
      INTEGER :: shell = 1
      REAL(dp) :: lattice(3, 3) = 0.0_dp, inverse(3, 3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: centres(:, :)
      TYPE(seed_environment_type), ALLOCATABLE :: neighbours(:)
      TYPE(seed_environment_type) :: self_neighbours
   END TYPE skala_partition_type

   PUBLIC :: skala_partition_type, skala_partition_init, skala_partition_weights, &
             skala_partition_pullback, skala_atom_grid_build, skala_atom_grid_pullback, &
             skala_atom_grid_direction, skala_radial_grid

CONTAINS

! **************************************************************************************************
!> \brief Cache translated Becke seed environments. Lattice vectors are columns, in bohr.
!> \param partition Cached periodic partition geometry.
!> \param lattice Cell vectors as columns, in bohr.
!> \param centres Cartesian atom centres, in bohr.
!> \param shell Retained lattice-image shell.
! **************************************************************************************************
   SUBROUTINE skala_partition_init(partition, lattice, centres, shell)
      TYPE(skala_partition_type), INTENT(OUT)            :: partition
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: shell

      INTEGER                                            :: a

      CPASSERT(shell >= 0)
      CPASSERT(SIZE(centres, 1) == 3 .AND. SIZE(centres, 2) > 0)
      CPASSERT(ALL(ieee_is_finite(lattice)) .AND. ALL(ieee_is_finite(centres)))
      CPASSERT(ABS(det_3x3(lattice)) > 1.0e-12_dp)
      partition%shell = shell
      partition%lattice = lattice
      partition%inverse = inv_3x3(lattice)
      ALLOCATE (partition%centres(3, SIZE(centres, 2)), partition%neighbours(SIZE(centres, 2)))
      partition%centres(:, :) = centres
      DO a = 1, SIZE(centres, 2)
         CALL build_environment(partition, a, .FALSE., partition%neighbours(a))
      END DO
      CALL build_environment(partition, 1, .TRUE., partition%self_neighbours)
   END SUBROUTINE skala_partition_init

! **************************************************************************************************
!> \brief Each image has the same translated environment, not a finite-cluster edge environment.
!> \param partition Cached periodic partition geometry.
!> \param atom Target atom index.
!> \param self_only Use only translated images of the target atom.
!> \param environment Translated seed environment.
! **************************************************************************************************
   SUBROUTINE build_environment(partition, atom, self_only, environment)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      LOGICAL, INTENT(IN)                                :: self_only
      TYPE(seed_environment_type), INTENT(OUT)           :: environment

      INTEGER                                            :: b, capacity, count, hi(3), i, j, k, lo(3)
      INTEGER, ALLOCATABLE                               :: atoms(:), images(:, :)
      REAL(dp)                                           :: capacity_real, delta(3), distance, &
                                                            fractional(3), weight
      REAL(dp), ALLOCATABLE                              :: distances(:), offsets(:, :), weights(:)

      capacity_real = REAL(SIZE(partition%centres, 2), dp)*(2.0_dp*partition%shell + 3)**3
      CPASSERT(capacity_real < REAL(HUGE(capacity), dp))
      capacity = SIZE(partition%centres, 2)*(2*partition%shell + 3)**3
      ALLOCATE (offsets(3, capacity), weights(capacity), distances(capacity))
      ALLOCATE (atoms(capacity), images(3, capacity))
      count = 0
      DO b = 1, SIZE(partition%centres, 2)
         IF (self_only .AND. b /= atom) CYCLE
         delta = MATMUL(partition%inverse, partition%centres(:, b) - partition%centres(:, atom))
         lo = CEILING(-delta - (partition%shell + 1))
         hi = FLOOR(-delta + (partition%shell + 1))
         DO k = lo(3), hi(3)
            DO j = lo(2), hi(2)
               DO i = lo(1), hi(1)
                  IF (b == atom .AND. ALL([i, j, k] == 0)) CYCLE
                  fractional = delta + [i, j, k]
                  weight = compact_window(fractional, partition%shell)
                  IF (weight == 0.0_dp) CYCLE
                  distance = NORM2(MATMUL(partition%lattice, fractional))
                  IF (distance < 1.0e-10_dp) THEN
                     CPABORT("Coincident atoms in the periodic Skala quadrature")
                  END IF
                  count = count + 1
                  offsets(:, count) = MATMUL(partition%lattice, fractional)
                  distances(count) = distance
                  weights(count) = weight
                  atoms(count) = b
                  images(:, count) = [i, j, k]
               END DO
            END DO
         END DO
      END DO
      ALLOCATE (environment%displacement(3, count), environment%distance(count), environment%membership(count))
      ALLOCATE (environment%atom(count), environment%image(3, count))
      environment%displacement(:, :) = offsets(:, :count)
      environment%distance(:) = distances(:count)
      environment%membership(:) = weights(:count)
      environment%atom(:) = atoms(:count)
      environment%image(:, :) = images(:, :count)
   END SUBROUTINE build_environment

! **************************************************************************************************
!> \brief Energy partition of atom at zero image and descriptor turn-on of its self-image partition.
!> \param partition Cached fixed geometry.
!> \param atom Target atom block.
!> \param offset Cartesian displacement from the target atom, in bohr.
!> \param energy All-atom, all-image partition weight.
!> \param descriptor Self-image domain window, independent of all other atoms.
! **************************************************************************************************
   SUBROUTINE skala_partition_weights(partition, atom, offset, energy, descriptor)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3)
      REAL(dp), INTENT(OUT)                              :: energy, descriptor

      REAL(dp)                                           :: d, z

      CPASSERT(atom >= 1 .AND. atom <= SIZE(partition%centres, 2))
      CPASSERT(ALL(ieee_is_finite(offset)))
      energy = 0.0_dp
      descriptor = 0.0_dp
      IF (compact_window(MATMUL(partition%inverse, offset), partition%shell) == 0.0_dp) RETURN
      d = normalized_seed(partition, atom, offset, .TRUE.)
      IF (SIZE(partition%centres, 2) == 1) THEN
         energy = d
      ELSE
         energy = normalized_seed(partition, atom, offset, .FALSE.)
      END IF
      IF (d >= 1.0e-11_dp) THEN
         descriptor = 1.0_dp
      ELSE IF (d > 1.0e-12_dp) THEN
         z = (d - 1.0e-12_dp)/9.0e-12_dp
         descriptor = z**3*(10.0_dp - 15.0_dp*z + 6.0_dp*z*z)
      END IF
   END SUBROUTINE skala_partition_weights

! **************************************************************************************************
!> \brief Normalize all nonzero translated seeds using a scaled log sum.
!> \param partition Cached periodic partition geometry.
!> \param atom Target atom index.
!> \param offset Displacement from the target atom, in bohr.
!> \param self_only Use only translated images of the target atom.
!> \param log_max Largest log seed used in the normalization.
!> \param scaled_sum Sum of seeds relative to log_max.
!> \return Normalized periodic seed weight.
! **************************************************************************************************
   FUNCTION normalized_seed(partition, atom, offset, self_only, log_max, scaled_sum) RESULT(weight)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3)
      LOGICAL, INTENT(IN)                                :: self_only
      REAL(dp), INTENT(OUT), OPTIONAL                    :: log_max, scaled_sum
      REAL(dp)                                           :: weight

      INTEGER                                            :: b, hi(3), i, j, k, lo(3)
      REAL(dp)                                           :: fractional(3), largest, relative(3), &
                                                            TARGET, total, value, x(3)

      ! A finite sentinel avoids infinities in builds trapping floating-point exceptions.
      TARGET = -HUGE(1.0_dp)
      largest = TARGET
      total = 0.0_dp
      DO b = 1, SIZE(partition%centres, 2)
         IF (self_only .AND. b /= atom) CYCLE
         relative = offset + partition%centres(:, atom) - partition%centres(:, b)
         fractional = MATMUL(partition%inverse, relative)
         lo = CEILING(fractional - (partition%shell + 1))
         hi = FLOOR(fractional + (partition%shell + 1))
         DO k = lo(3), hi(3)
            DO j = lo(2), hi(2)
               DO i = lo(1), hi(1)
                  x = relative - MATMUL(partition%lattice, REAL([i, j, k], dp))
                  IF (self_only) THEN
                     value = log_seed(partition, partition%self_neighbours, x)
                  ELSE
                     value = log_seed(partition, partition%neighbours(b), x)
                  END IF
                  IF (b == atom .AND. ALL([i, j, k] == 0)) TARGET = value
                  IF (value == -HUGE(1.0_dp)) CYCLE
                  IF (value > largest) THEN
                     total = total*EXP(largest - value) + 1.0_dp
                     largest = value
                  ELSE
                     total = total + EXP(value - largest)
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF (total <= 0.0_dp) CPABORT("Empty periodic Skala partition")
      weight = 0.0_dp
      IF (TARGET /= -HUGE(1.0_dp)) weight = EXP(TARGET - largest)/total
      IF (PRESENT(log_max)) log_max = largest
      IF (PRESENT(scaled_sum)) scaled_sum = total
   END FUNCTION normalized_seed

! **************************************************************************************************
!> \brief Pull back independent energy and descriptor partition weights at fixed row topology.
!> \param partition Cached geometry.
!> \param atom Target atom block.
!> \param offset Cartesian displacement from that atom, independent of the lattice.
!> \param energy_adjoint Derivative with respect to the energy partition weight.
!> \param descriptor_adjoint Derivative with respect to the descriptor domain window.
!> \param doffset Derivative with respect to the independent Cartesian offset.
!> \param dcentres Derivatives with respect to all Cartesian atom positions at fixed offset.
!> \param dlattice Derivative with respect to lattice columns at fixed Cartesian centres and offset.
! **************************************************************************************************
   SUBROUTINE skala_partition_pullback(partition, atom, offset, energy_adjoint, descriptor_adjoint, &
                                       doffset, dcentres, dlattice)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3), energy_adjoint, &
                                                            descriptor_adjoint
      REAL(dp), INTENT(OUT)                              :: doffset(3), dcentres(:, :), &
                                                            dlattice(3, 3)

      REAL(dp)                                           :: adjoint, d, z

      CPASSERT(atom >= 1 .AND. atom <= SIZE(partition%centres, 2))
      CPASSERT(ALL(SHAPE(dcentres) == SHAPE(partition%centres)))
      CPASSERT(ALL(ieee_is_finite(offset)))
      CPASSERT(ieee_is_finite(energy_adjoint) .AND. ieee_is_finite(descriptor_adjoint))
      doffset = 0.0_dp
      dcentres = 0.0_dp
      dlattice = 0.0_dp
      IF (compact_window(MATMUL(partition%inverse, offset), partition%shell) == 0.0_dp) RETURN
      adjoint = 0.0_dp
      IF (descriptor_adjoint /= 0.0_dp) THEN
         d = normalized_seed(partition, atom, offset, .TRUE.)
         IF (d > 1.0e-12_dp .AND. d < 1.0e-11_dp) THEN
            z = (d - 1.0e-12_dp)/9.0e-12_dp
            adjoint = descriptor_adjoint*30.0_dp*z*z*(1.0_dp - z)**2/9.0e-12_dp
         END IF
      END IF
      IF (SIZE(partition%centres, 2) == 1) THEN
         adjoint = adjoint + energy_adjoint
      ELSE IF (energy_adjoint /= 0.0_dp) THEN
         CALL seed_pullback(partition, atom, offset, .FALSE., energy_adjoint, doffset, dcentres, dlattice)
      END IF
      IF (adjoint /= 0.0_dp) THEN
         CALL seed_pullback(partition, atom, offset, .TRUE., adjoint, doffset, dcentres, dlattice)
      END IF
      CPASSERT(ALL(ieee_is_finite(doffset)))
      CPASSERT(ALL(ieee_is_finite(dcentres)))
      CPASSERT(ALL(ieee_is_finite(dlattice)))
   END SUBROUTINE skala_partition_pullback

! **************************************************************************************************
!> \brief Reverse the normalized image sum without differentiating its log-sum scaling.
!> \param partition Cached periodic partition geometry.
!> \param atom Target atom index.
!> \param offset Displacement from the target atom, in bohr.
!> \param self_only Use only translated images of the target atom.
!> \param adjoint Incoming scalar covector.
!> \param doffset Displacement covector.
!> \param dcentres Atom-centre covectors.
!> \param dlattice Cell-vector covectors.
! **************************************************************************************************
   SUBROUTINE seed_pullback(partition, atom, offset, self_only, adjoint, doffset, dcentres, dlattice)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: offset(3)
      LOGICAL, INTENT(IN)                                :: self_only
      REAL(dp), INTENT(IN)                               :: adjoint
      REAL(dp), INTENT(INOUT)                            :: doffset(3), dcentres(:, :), &
                                                            dlattice(3, 3)

      INTEGER                                            :: b, hi(3), i, j, k, lo(3)
      REAL(dp)                                           :: coefficient, fractional(3), gx(3), &
                                                            largest, relative(3), total, value, &
                                                            weight, x(3)

      weight = normalized_seed(partition, atom, offset, self_only, largest, total)
      IF (weight == 0.0_dp) RETURN
      DO b = 1, SIZE(partition%centres, 2)
         IF (self_only .AND. b /= atom) CYCLE
         relative = offset + partition%centres(:, atom) - partition%centres(:, b)
         fractional = MATMUL(partition%inverse, relative)
         lo = CEILING(fractional - (partition%shell + 1))
         hi = FLOOR(fractional + (partition%shell + 1))
         DO k = lo(3), hi(3)
            DO j = lo(2), hi(2)
               DO i = lo(1), hi(1)
                  x = relative - MATMUL(partition%lattice, REAL([i, j, k], dp))
                  IF (self_only) THEN
                     value = log_seed(partition, partition%self_neighbours, x)
                  ELSE
                     value = log_seed(partition, partition%neighbours(b), x)
                  END IF
                  IF (value == -HUGE(1.0_dp)) CYCLE
                  coefficient = -EXP(value - largest)/total
                  IF (b == atom .AND. ALL([i, j, k] == 0)) coefficient = coefficient + 1.0_dp
                  coefficient = coefficient*weight*adjoint
                  IF (coefficient == 0.0_dp) CYCLE
                  IF (self_only) THEN
                     CALL log_seed_pullback(partition, partition%self_neighbours, b, .TRUE., x, coefficient, &
                                            gx, dcentres, dlattice)
                  ELSE
                     CALL log_seed_pullback(partition, partition%neighbours(b), b, .FALSE., x, coefficient, &
                                            gx, dcentres, dlattice)
                  END IF
                  doffset = doffset + gx
                  dcentres(:, atom) = dcentres(:, atom) + gx
                  dcentres(:, b) = dcentres(:, b) - gx
                  dlattice = dlattice - SPREAD(gx, 2, 3)*SPREAD(REAL([i, j, k], dp), 1, 3)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE seed_pullback

! **************************************************************************************************
!> \brief Differentiate a nonzero seed, including image distances and neighbour membership.
!> \param partition Cached periodic partition geometry.
!> \param environment Translated seed environment.
!> \param atom Target atom index.
!> \param self_only Use only translated images of the target atom.
!> \param x Cartesian displacement, in bohr.
!> \param adjoint Incoming scalar covector.
!> \param gx Cartesian displacement covector.
!> \param dcentres Atom-centre covectors.
!> \param dlattice Cell-vector covectors.
! **************************************************************************************************
   SUBROUTINE log_seed_pullback(partition, environment, atom, self_only, x, adjoint, gx, dcentres, dlattice)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      TYPE(seed_environment_type), INTENT(IN)            :: environment
      INTEGER, INTENT(IN)                                :: atom
      LOGICAL, INTENT(IN)                                :: self_only
      REAL(dp), INTENT(IN)                               :: x(3), adjoint
      REAL(dp), INTENT(OUT)                              :: gx(3)
      REAL(dp), INTENT(INOUT)                            :: dcentres(:, :), dlattice(3, 3)

      INTEGER                                            :: b, n
      REAL(dp) :: c, delta(3), distance, factor, fractional(3), gc(3), gd(3), gmu, mu0, radius, &
         slope, switch, unit(3), unit_x(3), window_grad(3)

      fractional = MATMUL(partition%inverse, x)
      c = compact_window(fractional, partition%shell, window_grad)
      CPASSERT(c > 0.0_dp)
      gx = MATMUL(TRANSPOSE(partition%inverse), window_grad)*(adjoint/c)
      dlattice = dlattice - SPREAD(gx, 2, 3)*SPREAD(fractional, 1, 3)
      radius = NORM2(x)
      unit_x = 0.0_dp
      IF (radius > 0.0_dp) unit_x = x/radius
      DO n = 1, SIZE(environment%distance)
         delta = environment%displacement(:, n)
         distance = NORM2(x - delta)
         unit = 0.0_dp
         IF (distance > 0.0_dp) unit = (x - delta)/distance
         mu0 = (radius - distance)/environment%distance(n)
         CALL becke_switch(mu0, switch, slope)
         c = environment%membership(n)
         factor = (1.0_dp - c) + c*switch
         IF (factor == 0.0_dp) THEN
            ! A tiny forward seed can round to a saturated pair switch on recomputation.
            ! The weighted seed derivative vanishes at this endpoint, as does membership's gradient.
            CPASSERT(c == 1.0_dp .AND. switch == 0.0_dp .AND. slope == 0.0_dp)
            CYCLE
         END IF
         CPASSERT(factor > 0.0_dp)
         gmu = adjoint*c*slope/factor
         gx = gx + gmu*(unit_x - unit)/environment%distance(n)
         gd = gmu*(unit - mu0*delta/environment%distance(n))/environment%distance(n)
         fractional = MATMUL(partition%inverse, delta)
         c = compact_window(fractional, partition%shell, window_grad)
         gc = MATMUL(TRANSPOSE(partition%inverse), window_grad)*(adjoint*(switch - 1.0_dp)/factor)
         gd = gd + gc
         b = environment%atom(n)
         IF (self_only) b = atom
         dcentres(:, b) = dcentres(:, b) + gd
         dcentres(:, atom) = dcentres(:, atom) - gd
         dlattice = dlattice + SPREAD(gd, 2, 3)*SPREAD(REAL(environment%image(:, n), dp), 1, 3) - &
                    SPREAD(gc, 2, 3)*SPREAD(fractional, 1, 3)
      END DO
   END SUBROUTINE log_seed_pullback

! **************************************************************************************************
!> \brief Compact Becke seed with triple-cubic pair switches and smooth neighbour membership.
!> \param partition Cached periodic partition geometry.
!> \param environment Translated seed environment.
!> \param x Cartesian displacement, in bohr.
!> \return Log seed, or a finite sentinel for zero weight.
! **************************************************************************************************
   FUNCTION log_seed(partition, environment, x) RESULT(value)
      TYPE(skala_partition_type), INTENT(IN)             :: partition
      TYPE(seed_environment_type), INTENT(IN)            :: environment
      REAL(dp), INTENT(IN)                               :: x(3)
      REAL(dp)                                           :: value

      INTEGER                                            :: n
      REAL(dp)                                           :: c, factor, mu, radius, switch

      value = -HUGE(1.0_dp)
      c = compact_window(MATMUL(partition%inverse, x), partition%shell)
      IF (c == 0.0_dp) RETURN
      value = LOG(c)
      radius = NORM2(x)
      DO n = 1, SIZE(environment%distance)
         mu = (radius - NORM2(x - environment%displacement(:, n)))/environment%distance(n)
         CALL becke_switch(mu, switch)
         c = environment%membership(n)
         factor = (1.0_dp - c) + c*switch
         IF (factor == 0.0_dp) THEN
            value = -HUGE(1.0_dp)
            RETURN
         END IF
         value = value + LOG(factor)
      END DO
   END FUNCTION log_seed

! **************************************************************************************************
!> \brief Triple-cubic Becke switch without the final subtraction of nearly equal numbers.
!> \param mu Normalized pair-distance difference.
!> \param weight Becke-switch value.
!> \param derivative Derivative of the switch with respect to mu.
! **************************************************************************************************
   PURE SUBROUTINE becke_switch(mu, weight, derivative)
      REAL(dp), INTENT(IN)                               :: mu
      REAL(dp), INTENT(OUT)                              :: weight
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative

      INTEGER                                            :: repetition
      REAL(dp)                                           :: slope

      weight = 0.5_dp*(1.0_dp - MAX(-1.0_dp, MIN(1.0_dp, mu)))
      slope = -0.5_dp
      DO repetition = 1, 3
         slope = slope*6.0_dp*weight*(1.0_dp - weight)
         weight = weight*weight*(3.0_dp - 2.0_dp*weight)
      END DO
      IF (PRESENT(derivative)) derivative = slope
   END SUBROUTINE becke_switch

! **************************************************************************************************
!> \brief C2 compact window, evaluated without cancellation at the outer edge.
!> \param fractional Fractional displacement.
!> \param shell Retained lattice-image shell.
!> \param gradient Window derivative with respect to fractional coordinates.
!> \return Compact-window weight in the interval from zero to one.
! **************************************************************************************************
   FUNCTION compact_window(fractional, shell, gradient) RESULT(weight)
      REAL(dp), INTENT(IN)                               :: fractional(3)
      INTEGER, INTENT(IN)                                :: shell
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(3)
      REAL(dp)                                           :: weight

      INTEGER                                            :: direction
      REAL(dp)                                           :: factors(3), slope(3), z

      weight = 1.0_dp
      factors = 1.0_dp
      slope = 0.0_dp
      IF (PRESENT(gradient)) gradient = 0.0_dp
      DO direction = 1, 3
         z = MAX(0.0_dp, ABS(fractional(direction)) - shell)
         IF (z >= 1.0_dp) THEN
            weight = 0.0_dp
            RETURN
         END IF
         IF (PRESENT(gradient)) THEN
            factors(direction) = (1.0_dp - z)**3*(1.0_dp + 3.0_dp*z + 6.0_dp*z*z)
            slope(direction) = -30.0_dp*z*z*(1.0_dp - z)**2*SIGN(1.0_dp, fractional(direction))
         END IF
         weight = weight*(1.0_dp - z)**3*(1.0_dp + 3.0_dp*z + 6.0_dp*z*z)
      END DO
      IF (PRESENT(gradient)) THEN
         DO direction = 1, 3
            gradient(direction) = slope(direction)*PRODUCT(factors(:direction - 1))*PRODUCT(factors(direction + 1:))
         END DO
      END IF
   END FUNCTION compact_window

! **************************************************************************************************
!> \brief Gauss-Legendre radial rule on [0,infinity), including r^2 dr and a scale in bohr.
!> \param n Number of radial nodes.
!> \param scale Radial scale, in bohr.
!> \param radius Radial nodes, in bohr.
!> \param weights Radial weights including r squared dr, in cubic bohr.
! **************************************************************************************************
   SUBROUTINE skala_radial_grid(n, scale, radius, weights)
      INTEGER, INTENT(IN)                                :: n
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), INTENT(OUT)                              :: radius(n), weights(n)

      INTEGER                                            :: i, info
      INTEGER, SAVE                                      :: cached_n = 0
      REAL(dp)                                           :: offdiag(MAX(1, n - 1)), t, &
                                                            work(MAX(1, 2*n - 2))
      REAL(dp), ALLOCATABLE                              :: eigenvectors(:, :)
      REAL(dp), ALLOCATABLE, SAVE                        :: cached_radius(:), cached_weights(:)
      REAL(dp), SAVE                                     :: cached_scale = 0.0_dp

      CPASSERT(n >= 2 .AND. ieee_is_finite(scale) .AND. scale > 0.0_dp)
      ! Geometry derivatives reuse this fixed rule, not atom partitions or model fields.
!$OMP CRITICAL(skala_radial_rule)
      IF (n /= cached_n .OR. scale /= cached_scale) THEN
         ALLOCATE (eigenvectors(n, n))
         ! Golub-Welsch with the Legendre Jacobi matrix. LAPACK supplies the nodes and weights.
         radius = 0.0_dp
         DO i = 1, n - 1
            offdiag(i) = REAL(i, dp)/SQRT(4.0_dp*REAL(i, dp)**2 - 1.0_dp)
         END DO
         CALL DSTEV('V', n, radius, offdiag, eigenvectors, n, work, info)
         CPASSERT(info == 0)
         DO i = 1, n
            t = 0.5_dp*(1.0_dp + radius(i))
            radius(i) = scale*t/(1.0_dp - t)
            weights(i) = eigenvectors(1, i)**2*radius(i)**2*scale/(1.0_dp - t)**2
         END DO
         IF (n /= cached_n) THEN
            IF (ALLOCATED(cached_radius)) DEALLOCATE (cached_radius, cached_weights)
            ALLOCATE (cached_radius(n), cached_weights(n))
         END IF
         cached_radius(:) = radius
         cached_weights(:) = weights
         cached_n = n
         cached_scale = scale
      ELSE
         radius = cached_radius
         weights = cached_weights
      END IF
!$OMP END CRITICAL(skala_radial_rule)
   END SUBROUTINE skala_radial_grid

! **************************************************************************************************
!> \brief Build Gauss-Legendre/Lebedev blocks, owned once across the communicator.
!> \param group Context communicator.
!> \param lattice Cell vectors as columns, in bohr.
!> \param centres Cartesian atom centres, in bohr.
!> \param nradial Number of radial nodes.
!> \param nangular Lebedev point count.
!> \param shell Retained lattice-image shell.
!> \param scale Radial scale, in bohr.
!> \param points Cartesian quadrature points, in bohr.
!> \param energy_weights Energy-partition weights.
!> \param descriptor_weights Nonlocal descriptor-domain weights.
!> \param block_sizes Point count of each local atom block.
!> \param block_centres Cartesian atom-block centres, in bohr.
!> \note Never drop a descriptor-only row.
! **************************************************************************************************
   SUBROUTINE skala_atom_grid_build(group, lattice, centres, nradial, nangular, shell, scale, &
                                    points, energy_weights, descriptor_weights, block_sizes, block_centres)
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: nradial, nangular, shell
      REAL(dp), INTENT(IN)                               :: scale
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: points(:, :), energy_weights(:), &
                                                            descriptor_weights(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: block_sizes(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: block_centres(:, :)

      INTEGER                                            :: a, angular_grid, b, ia, ir, max_rows, &
                                                            na, nb, row, start
      REAL(dp)                                           :: base_weight, d, offset(3), p
      REAL(dp), ALLOCATABLE                              :: descriptor(:), energy(:), &
                                                            positions(:, :), radial_weights(:), &
                                                            radii(:)
      TYPE(skala_partition_type)                         :: partition

      CPASSERT(nradial >= 2 .AND. nangular >= 6 .AND. nangular <= 974)
      CPASSERT(shell >= 0 .AND. scale > 0.0_dp .AND. ieee_is_finite(scale))
      CALL skala_partition_init(partition, lattice, centres, shell)
      angular_grid = get_number_of_lebedev_grid(n=nangular)
      na = lebedev_grid(angular_grid)%n
      IF (na /= nangular) THEN
         CPABORT("Select an available Lebedev point count for the Skala atom grid")
      END IF
      IF (ANY(lebedev_grid(angular_grid)%w < 0.0_dp)) THEN
         CPABORT("Skala descriptor weights require a positive Lebedev rule")
      END IF
      nb = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         nb = nb + 1
      END DO
      CPASSERT(REAL(nb, dp)*nradial*na < REAL(HUGE(max_rows), dp))
      max_rows = nb*nradial*na
      ALLOCATE (positions(3, max_rows), energy(max_rows), descriptor(max_rows), block_sizes(nb), block_centres(3, nb))
      ALLOCATE (radii(nradial), radial_weights(nradial))
      CALL skala_radial_grid(nradial, scale, radii, radial_weights)
      row = 0
      b = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         b = b + 1
         start = row
         block_centres(:, b) = centres(:, a)
         DO ir = 1, nradial
            DO ia = 1, na
               offset = radii(ir)*lebedev_grid(angular_grid)%r(:, ia)
               CALL skala_partition_weights(partition, a, offset, p, d)
               IF (p == 0.0_dp .AND. d == 0.0_dp) CYCLE
               base_weight = fourpi*radial_weights(ir)*lebedev_grid(angular_grid)%w(ia)
               row = row + 1
               positions(:, row) = centres(:, a) + offset
               energy(row) = base_weight*p
               descriptor(row) = base_weight*d
            END DO
         END DO
         block_sizes(b) = row - start
         IF (block_sizes(b) == 0) CPABORT("Empty Skala atom block. Refine the radial quadrature")
      END DO
      ALLOCATE (points(3, row), energy_weights(row), descriptor_weights(row))
      points(:, :) = positions(:, :row)
      energy_weights(:) = energy(:row)
      descriptor_weights(:) = descriptor(:row)
   END SUBROUTINE skala_atom_grid_build

! **************************************************************************************************
!> \brief Differentiate the local atom quadrature at fixed radial scale, angular rule and row topology.
!> \param group Context communicator.
!> \param lattice Cell vectors as columns, in bohr.
!> \param centres Cartesian atom centres, in bohr.
!> \param nradial Number of radial nodes.
!> \param nangular Lebedev point count.
!> \param shell Retained lattice-image shell.
!> \param scale Radial scale, in bohr.
!> \param block_sizes Point count of each local atom block.
!> \param atom_direction Cartesian atom displacement directions.
!> \param lattice_direction Cell-vector directions.
!> \param points Quadrature-point directions, in bohr.
!> \param energy_weights Energy-weight directions.
!> \param descriptor_weights Descriptor-weight directions.
!> \param block_centres Atom-block-centre directions, in bohr.
!> \note Cartesian atom and lattice-column directions are independent and replicated.
!>       Radial offsets do not deform affinely with the cell.
! **************************************************************************************************
   SUBROUTINE skala_atom_grid_direction(group, lattice, centres, nradial, nangular, shell, scale, block_sizes, &
                                        atom_direction, lattice_direction, points, energy_weights, &
                                        descriptor_weights, block_centres)
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: nradial, nangular, shell
      REAL(dp), INTENT(IN)                               :: scale
      INTEGER, INTENT(IN)                                :: block_sizes(:)
      REAL(dp), INTENT(IN)                               :: atom_direction(:, :), &
                                                            lattice_direction(3, 3)
      REAL(dp), INTENT(OUT)                              :: points(:, :), energy_weights(:), &
                                                            descriptor_weights(:), &
                                                            block_centres(:, :)

      INTEGER                                            :: a, angular_grid, b, ia, ir, row, start
      REAL(dp)                                           :: base_weight, d, gh(3, 3), gx(3), &
                                                            offset(3), p, radial_weights(nradial), &
                                                            radii(nradial), gr(3, SIZE(centres, 2))
      TYPE(skala_partition_type)                         :: partition

      CPASSERT(nradial >= 2 .AND. nangular >= 6 .AND. nangular <= 974)
      CPASSERT(ALL(SHAPE(atom_direction) == SHAPE(centres)))
      CPASSERT(ALL(ieee_is_finite(atom_direction)))
      CPASSERT(ALL(ieee_is_finite(lattice_direction)))
      CPASSERT(ALL(SHAPE(points) == [3, SUM(block_sizes)]))
      CPASSERT(ALL(SHAPE(block_centres) == [3, SIZE(block_sizes)]))
      CPASSERT(ALL(block_sizes > 0))
      CPASSERT(SIZE(energy_weights) == SUM(block_sizes))
      CPASSERT(SIZE(descriptor_weights) == SUM(block_sizes))
      CALL skala_partition_init(partition, lattice, centres, shell)
      angular_grid = get_number_of_lebedev_grid(n=nangular)
      CPASSERT(lebedev_grid(angular_grid)%n == nangular)
      CPASSERT(ALL(lebedev_grid(angular_grid)%w >= 0.0_dp))
      CALL skala_radial_grid(nradial, scale, radii, radial_weights)
      points = 0.0_dp
      energy_weights = 0.0_dp
      descriptor_weights = 0.0_dp
      block_centres = 0.0_dp
      row = 0
      b = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         b = b + 1
         CPASSERT(b <= SIZE(block_sizes))
         start = row
         block_centres(:, b) = atom_direction(:, a)
         DO ir = 1, nradial
            DO ia = 1, nangular
               offset = radii(ir)*lebedev_grid(angular_grid)%r(:, ia)
               CALL skala_partition_weights(partition, a, offset, p, d)
               IF (p == 0.0_dp .AND. d == 0.0_dp) CYCLE
               row = row + 1
               CPASSERT(row <= SIZE(energy_weights))
               points(:, row) = atom_direction(:, a)
               base_weight = fourpi*radial_weights(ir)*lebedev_grid(angular_grid)%w(ia)
               ! Contract the existing analytic partition covectors with the host direction.
               CALL skala_partition_pullback(partition, a, offset, 1.0_dp, 0.0_dp, gx, gr, gh)
               energy_weights(row) = base_weight*(SUM(gr*atom_direction) + SUM(gh*lattice_direction))
               CALL skala_partition_pullback(partition, a, offset, 0.0_dp, 1.0_dp, gx, gr, gh)
               descriptor_weights(row) = base_weight*(SUM(gr*atom_direction) + SUM(gh*lattice_direction))
            END DO
         END DO
         CPASSERT(row - start == block_sizes(b))
      END DO
      CPASSERT(b == SIZE(block_sizes) .AND. row == SIZE(energy_weights))
      CPASSERT(ALL(ieee_is_finite(points)) .AND. ALL(ieee_is_finite(block_centres)))
      CPASSERT(ALL(ieee_is_finite(energy_weights)))
      CPASSERT(ALL(ieee_is_finite(descriptor_weights)))
   END SUBROUTINE skala_atom_grid_direction

! **************************************************************************************************
!> \brief Return global geometry covectors for the atom grid at fixed radial scale and angular rule.
!> \param group Context communicator.
!> \param lattice Cell vectors as columns, in bohr.
!> \param centres Cartesian atom centres, in bohr.
!> \param nradial Number of radial nodes.
!> \param nangular Lebedev point count.
!> \param shell Retained lattice-image shell.
!> \param scale Radial scale, in bohr.
!> \param block_sizes Point count of each local atom block.
!> \param dpoints Quadrature-point covectors.
!> \param denergy_weights Energy-weight covectors.
!> \param ddescriptor_weights Descriptor-weight covectors.
!> \param dblock_centres Atom-block-centre covectors.
!> \param dcentres Atom-centre covectors.
!> \param dlattice Cell-vector covectors.
!> \note Inputs are local raw covectors, including both independent quadrature weights.
!>       Include the point-position dependence of reconstructed fields in dpoints.
!>       The host adds their explicit atom and cell derivatives separately.
!>       dcentres and dlattice are energy derivatives, not forces or stress. Lattice columns
!>       vary at fixed Cartesian centres, without affine deformation of the radial offsets.
! **************************************************************************************************
   SUBROUTINE skala_atom_grid_pullback(group, lattice, centres, nradial, nangular, shell, scale, block_sizes, &
                                       dpoints, denergy_weights, ddescriptor_weights, dblock_centres, dcentres, dlattice)
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), INTENT(IN)                               :: lattice(3, 3), centres(:, :)
      INTEGER, INTENT(IN)                                :: nradial, nangular, shell
      REAL(dp), INTENT(IN)                               :: scale
      INTEGER, INTENT(IN)                                :: block_sizes(:)
      REAL(dp), INTENT(IN)                               :: dpoints(:, :), denergy_weights(:), &
                                                            ddescriptor_weights(:), &
                                                            dblock_centres(:, :)
      REAL(dp), INTENT(OUT)                              :: dcentres(:, :), dlattice(3, 3)

      INTEGER                                            :: a, angular_grid, b, ia, ir, row, start
      REAL(dp)                                           :: base_weight, d, gh(3, 3), gx(3), &
                                                            offset(3), p, radial_weights(nradial), &
                                                            radii(nradial), gr(3, SIZE(centres, 2))
      TYPE(skala_partition_type)                         :: partition

      CPASSERT(nradial >= 2 .AND. nangular >= 6 .AND. nangular <= 974)
      CPASSERT(ALL(SHAPE(dcentres) == SHAPE(centres)))
      CPASSERT(ALL(SHAPE(dpoints) == [3, SUM(block_sizes)]))
      CPASSERT(ALL(SHAPE(dblock_centres) == [3, SIZE(block_sizes)]))
      CPASSERT(ALL(block_sizes > 0))
      CPASSERT(SIZE(denergy_weights) == SUM(block_sizes))
      CPASSERT(SIZE(ddescriptor_weights) == SUM(block_sizes))
      CPASSERT(ALL(ieee_is_finite(dpoints)))
      CPASSERT(ALL(ieee_is_finite(dblock_centres)))
      CPASSERT(ALL(ieee_is_finite(denergy_weights)))
      CPASSERT(ALL(ieee_is_finite(ddescriptor_weights)))
      CALL skala_partition_init(partition, lattice, centres, shell)
      angular_grid = get_number_of_lebedev_grid(n=nangular)
      CPASSERT(lebedev_grid(angular_grid)%n == nangular)
      CPASSERT(ALL(lebedev_grid(angular_grid)%w >= 0.0_dp))
      CALL skala_radial_grid(nradial, scale, radii, radial_weights)
      dcentres = 0.0_dp
      dlattice = 0.0_dp
      row = 0
      b = 0
      DO a = group%mepos + 1, SIZE(centres, 2), group%num_pe
         b = b + 1
         CPASSERT(b <= SIZE(block_sizes))
         start = row
         dcentres(:, a) = dcentres(:, a) + dblock_centres(:, b)
         DO ir = 1, nradial
            DO ia = 1, nangular
               offset = radii(ir)*lebedev_grid(angular_grid)%r(:, ia)
               CALL skala_partition_weights(partition, a, offset, p, d)
               IF (p == 0.0_dp .AND. d == 0.0_dp) CYCLE
               row = row + 1
               CPASSERT(row <= SIZE(denergy_weights))
               base_weight = fourpi*radial_weights(ir)*lebedev_grid(angular_grid)%w(ia)
               CALL skala_partition_pullback(partition, a, offset, base_weight*denergy_weights(row), &
                                             base_weight*ddescriptor_weights(row), gx, gr, gh)
               dcentres = dcentres + gr
               dcentres(:, a) = dcentres(:, a) + dpoints(:, row)
               dlattice = dlattice + gh
            END DO
         END DO
         CPASSERT(row - start == block_sizes(b))
      END DO
      CPASSERT(b == SIZE(block_sizes) .AND. row == SIZE(denergy_weights))
      CALL group%sum(dcentres)
      CALL group%sum(dlattice)
      CPASSERT(ALL(ieee_is_finite(dcentres)))
      CPASSERT(ALL(ieee_is_finite(dlattice)))
   END SUBROUTINE skala_atom_grid_pullback
END MODULE skala_atom_grid

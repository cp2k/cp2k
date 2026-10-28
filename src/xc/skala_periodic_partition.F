!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Smooth atom-image weights normalized over translated, equivalent environments.
! **************************************************************************************************
MODULE skala_periodic_partition
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   INTEGER, PARAMETER :: image_shell = 1
   TYPE, PUBLIC :: skala_periodic_partition_type
      PRIVATE
      INTEGER :: natom = 0, nneighbor = 0, periodicity(3) = 0
      REAL(dp) :: cell(3, 3) = 0.0_dp, inverse(3, 3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: atoms(:, :), offset(:, :, :), inverse_distance(:, :)
      REAL(dp), ALLOCATABLE :: cutoff(:, :), cutoff_gradient(:, :, :)
      INTEGER, ALLOCATABLE :: owner(:, :), self(:)
   END TYPE skala_periodic_partition_type
   PUBLIC :: build_periodic_atom_image_layout, periodic_atom_image_partition
   PUBLIC :: periodic_atom_image_partition_from_layout
CONTAINS

! **************************************************************************************************
!> \brief Prepare geometry shared by the points of an atom block.
!> \param atom_coords Base-atom coordinates.
!> \param cell Simulation cell.
!> \param target_atom Reference-cell atom whose weight will be evaluated.
!> \param image_periodicity Optional periodic directions, including auxiliary-grid periodicity.
!> \param layout Periodized compact Becke geometry.
!> \param target_image Reference-atom index in the layout.
! **************************************************************************************************
   SUBROUTINE build_periodic_atom_image_layout(atom_coords, cell, target_atom, image_periodicity, &
                                               layout, target_image)
      REAL(dp), INTENT(IN)                               :: atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: target_atom
      INTEGER, INTENT(IN), OPTIONAL                      :: image_periodicity(3)
      TYPE(skala_periodic_partition_type), INTENT(OUT)   :: layout
      INTEGER, INTENT(OUT)                               :: target_image

      INTEGER                                            :: a, b, base(3), extent(3), i, j, k, n, &
                                                            shift(3)
      REAL(dp)                                           :: distance

      CPASSERT(ASSOCIATED(cell))
      layout%natom = SIZE(atom_coords, 2)
      CPASSERT(SIZE(atom_coords, 1) == 3 .AND. layout%natom > 0)
      CPASSERT(target_atom >= 1 .AND. target_atom <= layout%natom)
      layout%periodicity = cell%perd
      IF (PRESENT(image_periodicity)) layout%periodicity = image_periodicity
      CPASSERT(ALL(layout%periodicity == 0 .OR. layout%periodicity == 1))
      layout%cell = cell%hmat
      layout%inverse = cell%h_inv
      ALLOCATE (layout%atoms(3, layout%natom))
      layout%atoms(:, :) = atom_coords
      extent = (image_shell + 1)*layout%periodicity
      layout%nneighbor = layout%natom*PRODUCT(2*extent + 1)
      ALLOCATE (layout%offset(3, layout%nneighbor, layout%natom), &
                layout%inverse_distance(layout%nneighbor, layout%natom), &
                layout%cutoff(layout%nneighbor, layout%natom), &
                layout%cutoff_gradient(3, layout%nneighbor, layout%natom), &
                layout%owner(layout%nneighbor, layout%natom), layout%self(layout%natom))
      DO a = 1, layout%natom
         n = 0
         DO b = 1, layout%natom
            base = NINT(MATMUL(layout%inverse, atom_coords(:, a) - atom_coords(:, b)))*layout%periodicity
            DO k = -extent(3), extent(3)
               DO j = -extent(2), extent(2)
                  DO i = -extent(1), extent(1)
                     n = n + 1
                     shift = base + [i, j, k]
                     layout%owner(n, a) = b
                     layout%offset(:, n, a) = atom_coords(:, b) - atom_coords(:, a) + &
                                              MATMUL(layout%cell, REAL(shift, dp))
                     CALL compact_window(layout, layout%offset(:, n, a), &
                                         layout%cutoff(n, a), layout%cutoff_gradient(:, n, a))
                     IF (a == b .AND. ALL(shift == 0)) THEN
                        layout%self(a) = n
                        layout%inverse_distance(n, a) = 0.0_dp
                     ELSE
                        distance = NORM2(layout%offset(:, n, a))
                        CPASSERT(distance > 1.e-10_dp)
                        layout%inverse_distance(n, a) = 1.0_dp/distance
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      target_image = target_atom
   END SUBROUTINE build_periodic_atom_image_layout

! **************************************************************************************************
!> \brief Evaluate a reference-atom weight and optional base-atom and explicit image-strain responses.
!> \param grid_point Cartesian point held fixed when taking the derivatives.
!> \param layout Precomputed atom-image geometry.
!> \param target_image Reference-atom index returned by the builder.
!> \param weight Normalized reference-image weight.
!> \param dweight_datom Optional base-atom derivative at fixed grid point and cell.
!> \param dweight_dstrain Optional cell response at fixed grid point and base-atom coordinates.
! **************************************************************************************************
   SUBROUTINE periodic_atom_image_partition_from_layout(grid_point, layout, target_image, &
                                                        weight, dweight_datom, dweight_dstrain)
      REAL(dp), INTENT(IN)                               :: grid_point(3)
      TYPE(skala_periodic_partition_type), INTENT(IN)    :: layout
      INTEGER, INTENT(IN)                                :: target_image
      REAL(dp), INTENT(OUT)                              :: weight
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dweight_datom(:, :), &
                                                            dweight_dstrain(:, :)

      INTEGER                                            :: a, base(3), extent(3), i, j, k, n, &
                                                            shift(3), TARGET
      LOGICAL                                            :: zero
      REAL(dp)                                           :: coefficient, moment(3, 3), &
                                                            q(layout%nneighbor), &
                                                            datom(3, layout%natom), &
                                                            centers(3, layout%nneighbor)
      LOGICAL                                            :: active(layout%nneighbor)
      INTEGER                                            :: owner(layout%nneighbor)

      CPASSERT(target_image >= 1 .AND. target_image <= layout%natom)
      CPASSERT(PRESENT(dweight_datom) .EQV. PRESENT(dweight_dstrain))
      IF (PRESENT(dweight_datom)) THEN
         CPASSERT(SIZE(dweight_datom, 1) == 3)
         CPASSERT(SIZE(dweight_datom, 2) == layout%natom)
         CPASSERT(ALL(SHAPE(dweight_dstrain) == [3, 3]))
         dweight_datom = 0.0_dp
         dweight_dstrain = 0.0_dp
      END IF
      weight = 0.0_dp
      CALL seed_value(layout, target_image, grid_point - layout%atoms(:, target_image), weight, zero, .FALSE.)
      IF (zero) RETURN
      extent = (image_shell + 1)*layout%periodicity
      n = 0
      TARGET = 0
      DO a = 1, layout%natom
         base = NINT(MATMUL(layout%inverse, grid_point - layout%atoms(:, a)))*layout%periodicity
         DO k = -extent(3), extent(3)
            DO j = -extent(2), extent(2)
               DO i = -extent(1), extent(1)
                  n = n + 1
                  shift = base + [i, j, k]
                  owner(n) = a
                  centers(:, n) = layout%atoms(:, a) + MATMUL(layout%cell, REAL(shift, dp))
                  CALL seed_value(layout, a, grid_point - centers(:, n), q(n), zero, .FALSE.)
                  active(n) = .NOT. zero
                  IF (a == target_image .AND. ALL(shift == 0)) TARGET = n
               END DO
            END DO
         END DO
      END DO
      CPASSERT(TARGET > 0)
      IF (MAXVAL(q) < 1.e-200_dp) THEN
         DO n = 1, layout%nneighbor
            q(n) = -HUGE(1.0_dp)
            IF (.NOT. active(n)) CYCLE
            CALL seed_value(layout, owner(n), grid_point - centers(:, n), q(n), zero, .TRUE.)
         END DO
         q = EXP(q - MAXVAL(q))
      END IF
      CPASSERT(SUM(q) > 0.0_dp)
      q = q/SUM(q)
      weight = q(TARGET)
      IF (.NOT. PRESENT(dweight_datom)) RETURN
      datom = 0.0_dp
      moment = 0.0_dp
      DO n = 1, layout%nneighbor
         coefficient = -weight*q(n)
         IF (n == TARGET) coefficient = coefficient + weight
         IF (coefficient == 0.0_dp) CYCLE
         CALL seed_reverse(layout, owner(n), centers(:, n), grid_point, coefficient, datom, moment)
      END DO
      dweight_datom = datom
      dweight_dstrain = moment
      DO a = 1, layout%natom
         DO j = 1, 3
            dweight_dstrain(:, j) = dweight_dstrain(:, j) - datom(:, a)*layout%atoms(j, a)
         END DO
      END DO
   END SUBROUTINE periodic_atom_image_partition_from_layout

! **************************************************************************************************
!> \brief Convenience evaluation without retaining a layout between grid points.
!> \param grid_point Cartesian evaluation point.
!> \param atom_coords Base-atom coordinates.
!> \param cell Simulation cell.
!> \param target_atom Reference-cell atom.
!> \param weight Normalized reference-image weight.
!> \param dweight_datom Optional base-atom derivative at fixed grid point and cell.
!> \param dweight_dstrain Optional cell response at fixed grid point and base atoms.
!> \param image_periodicity Optional periodic directions.
! **************************************************************************************************
   SUBROUTINE periodic_atom_image_partition(grid_point, atom_coords, cell, target_atom, weight, &
                                            dweight_datom, dweight_dstrain, image_periodicity)
      REAL(dp), INTENT(IN)                               :: grid_point(3), atom_coords(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: target_atom
      REAL(dp), INTENT(OUT)                              :: weight
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dweight_datom(:, :), &
                                                            dweight_dstrain(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: image_periodicity(3)

      INTEGER                                            :: target
      TYPE(skala_periodic_partition_type)                :: layout

      CALL build_periodic_atom_image_layout(atom_coords, cell, target_atom, image_periodicity, layout, TARGET)
      CALL periodic_atom_image_partition_from_layout(grid_point, layout, TARGET, weight, &
                                                     dweight_datom, dweight_dstrain)
   END SUBROUTINE periodic_atom_image_partition

! **************************************************************************************************
!> \brief Compact quintic window in fractional coordinates.
!> \param layout Cell geometry and periodic directions.
!> \param x Cartesian displacement from the window center.
!> \param value Window value.
!> \param gradient Cartesian displacement gradient.
! **************************************************************************************************
   PURE SUBROUTINE compact_window(layout, x, value, gradient)
      TYPE(skala_periodic_partition_type), INTENT(IN)    :: layout
      REAL(dp), INTENT(IN)                               :: x(3)
      REAL(dp), INTENT(OUT)                              :: value, gradient(3)

      INTEGER                                            :: d, j
      REAL(dp)                                           :: dh(3), f(3), h(3), z

      f = MATMUL(layout%inverse, x)
      h = 1.0_dp
      dh = 0.0_dp
      DO d = 1, 3
         IF (layout%periodicity(d) == 0) CYCLE
         z = ABS(f(d)) - REAL(image_shell, dp)
         IF (z >= 1.0_dp) THEN
            value = 0.0_dp
            gradient = 0.0_dp
            RETURN
         END IF
         IF (z <= 0.0_dp) CYCLE
         h(d) = (1.0_dp - z)**3*(1.0_dp + 3.0_dp*z + 6.0_dp*z**2)
         dh(d) = -30.0_dp*z**2*(1.0_dp - z)**2*SIGN(1.0_dp, f(d))
      END DO
      value = PRODUCT(h)
      DO d = 1, 3
         DO j = 1, 3
            IF (j /= d) dh(d) = dh(d)*h(j)
         END DO
      END DO
      gradient = MATMUL(TRANSPOSE(layout%inverse), dh)
   END SUBROUTINE compact_window

! **************************************************************************************************
!> \brief Threefold Becke polynomial and its derivative.
!> \param mu Normalized difference of distances to the two centers.
!> \param value Pair weight.
!> \param derivative Derivative with respect to mu.
! **************************************************************************************************
   PURE SUBROUTINE pair_switch(mu, value, derivative)
      REAL(dp), INTENT(IN)                               :: mu
      REAL(dp), INTENT(OUT)                              :: value, derivative

      INTEGER                                            :: i

      value = 0.5_dp*(1.0_dp - MAX(-1.0_dp, MIN(1.0_dp, mu)))
      derivative = -0.5_dp
      DO i = 1, 3
         derivative = derivative*6.0_dp*value*(1.0_dp - value)
         value = value**2*(3.0_dp - 2.0_dp*value)
      END DO
   END SUBROUTINE pair_switch

! **************************************************************************************************
!> \brief Evaluate one compact atom-image seed using its translated neighbor environment.
!> \param layout Cached neighbor geometry.
!> \param atom Base atom identifying the image family.
!> \param x Evaluation point relative to the image center.
!> \param q Seed value or logarithm.
!> \param zero Whether the seed is identically zero at this point.
!> \param logarithmic Return the logarithm to avoid product underflow.
! **************************************************************************************************
   PURE SUBROUTINE seed_value(layout, atom, x, q, zero, logarithmic)
      TYPE(skala_periodic_partition_type), INTENT(IN)    :: layout
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: x(3)
      REAL(dp), INTENT(OUT)                              :: q
      LOGICAL, INTENT(OUT)                               :: zero
      LOGICAL, INTENT(IN)                                :: logarithmic

      INTEGER                                            :: k
      REAL(dp)                                           :: c, dc(3), ds, mu, ri, rj, s

      CALL compact_window(layout, x, c, dc)
      zero = c <= 0.0_dp
      q = 0.0_dp
      IF (zero) RETURN
      q = c
      IF (logarithmic) q = LOG(c)
      ri = NORM2(x)
      DO k = 1, layout%nneighbor
         IF (k == layout%self(atom) .OR. layout%cutoff(k, atom) == 0.0_dp) CYCLE
         rj = NORM2(x - layout%offset(:, k, atom))
         mu = (ri - rj)*layout%inverse_distance(k, atom)
         CALL pair_switch(mu, s, ds)
         s = (1.0_dp - layout%cutoff(k, atom)) + layout%cutoff(k, atom)*s
         IF (s <= 0.0_dp) THEN
            q = 0.0_dp
            zero = .TRUE.
            RETURN
         END IF
         IF (logarithmic) THEN
            q = q + LOG(s)
         ELSE
            q = q*s
         END IF
      END DO
   END SUBROUTINE seed_value

! **************************************************************************************************
!> \brief Reverse one log seed into base-atom derivatives and the affine strain moment.
!> \param layout Cached neighbor geometry.
!> \param atom Base atom identifying the image family.
!> \param center Cartesian image center.
!> \param point Cartesian evaluation point.
!> \param coefficient Adjoint of the log seed after normalization.
!> \param datom Accumulated base-atom derivatives.
!> \param moment Accumulated center moment including explicit inverse-cell responses.
! **************************************************************************************************
   SUBROUTINE seed_reverse(layout, atom, center, point, coefficient, datom, moment)
      TYPE(skala_periodic_partition_type), INTENT(IN)    :: layout
      INTEGER, INTENT(IN)                                :: atom
      REAL(dp), INTENT(IN)                               :: center(3), point(3), coefficient
      REAL(dp), INTENT(INOUT)                            :: datom(:, :), moment(3, 3)

      INTEGER                                            :: b, d, k
      REAL(dp)                                           :: c, dc(3), ds, e(3), factor, gi(3), &
                                                            gj(3), mu, neighbor(3), pair, ri, rj, &
                                                            s, u(3), v(3), x(3), y(3)

      x = point - center
      CALL compact_window(layout, x, c, dc)
      IF (c <= 0.0_dp .OR. coefficient == 0.0_dp) RETURN
      gi = -coefficient*dc/c
      datom(:, atom) = datom(:, atom) + gi
      ! Include the explicit inverse-cell response of the compact window.
      DO d = 1, 3
         moment(:, d) = moment(:, d) + gi*point(d)
      END DO
      ri = NORM2(x)
      u = x/MAX(ri, 1.e-30_dp)
      DO k = 1, layout%nneighbor
         IF (k == layout%self(atom) .OR. layout%cutoff(k, atom) == 0.0_dp) CYCLE
         b = layout%owner(k, atom)
         y = x - layout%offset(:, k, atom)
         rj = NORM2(y)
         v = y/MAX(rj, 1.e-30_dp)
         mu = (ri - rj)*layout%inverse_distance(k, atom)
         CALL pair_switch(mu, s, ds)
         pair = (1.0_dp - layout%cutoff(k, atom)) + layout%cutoff(k, atom)*s
         IF (pair <= 0.0_dp) CYCLE
         ! The neighbor taper's affine cell response cancels its center moment.
         gi = coefficient*(s - 1.0_dp)/pair*layout%cutoff_gradient(:, k, atom)
         datom(:, atom) = datom(:, atom) - gi
         datom(:, b) = datom(:, b) + gi
         factor = coefficient*layout%cutoff(k, atom)*ds/pair*layout%inverse_distance(k, atom)
         e = -layout%offset(:, k, atom)*layout%inverse_distance(k, atom)
         gi = factor*(-u - mu*e)
         gj = factor*(v + mu*e)
         datom(:, atom) = datom(:, atom) + gi
         datom(:, b) = datom(:, b) + gj
         neighbor = center + layout%offset(:, k, atom)
         DO d = 1, 3
            moment(:, d) = moment(:, d) + gi*center(d) + gj*neighbor(d)
         END DO
      END DO
   END SUBROUTINE seed_reverse
END MODULE skala_periodic_partition

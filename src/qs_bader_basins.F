!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Weighted zero-flux basin integration on an orthogonal real-space grid.
! **************************************************************************************************
MODULE qs_bader_basins
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: bader_integrate_grid

CONTAINS

! **************************************************************************************************
!> \brief Integrate positive-density basins using the weighted grid method of Yu and Trinkle.
!>        A topological traversal distributes charge and volume along strictly uphill edges.
!>        Atomic seed points identify their attractors along the strongest uphill edge.
!> \param rho spin-summed electron density, first grid index varying fastest
!> \param n number of grid points in each direction
!> \param spacing grid spacing in bohr
!> \param periodic whether each direction wraps
!> \param cutoff density at or below which points are counted as vacuum
!> \param seeds linear grid index nearest each nucleus, or zero if outside the grid
!> \param atomic electron population and volume (bohr^3) assigned to each atom
!> \param assigned whether the atom has a basin distinct from all other atoms
!> \param other integrated electron population and volume of unassigned basins
!> \param vacuum integrated electron population and volume of vacuum grid points
!> \param n_other number of unassigned maxima
!> \param atom_grid optional atom positions in grid coordinates, for nearest-atom grouping of maxima
!> \param n_extra optional number of maxima without a unique nuclear seed
! **************************************************************************************************
   SUBROUTINE bader_integrate_grid(rho, n, spacing, periodic, cutoff, seeds, atomic, assigned, other, vacuum, &
                                   n_other, atom_grid, n_extra)
      REAL(KIND=dp), INTENT(IN)                          :: rho(:), SPACING(3), cutoff
      INTEGER, INTENT(IN)                                :: n(3), seeds(:)
      LOGICAL, INTENT(IN)                                :: periodic(3)
      REAL(KIND=dp), INTENT(OUT)                         :: atomic(:, :), other(2), vacuum(2)
      LOGICAL, INTENT(OUT)                               :: assigned(:)
      INTEGER, INTENT(OUT)                               :: n_other
      REAL(KIND=dp), OPTIONAL, INTENT(IN)                :: atom_grid(:, :)
      INTEGER, OPTIONAL, INTENT(OUT)                     :: n_extra

      INTEGER                                            :: head, tail, active, i, j, k, root, atom, natom, ng
      INTEGER                                            :: neighbor(6)
      INTEGER, ALLOCATABLE                              :: degree(:), order(:), attractor(:), owner(:)
      REAL(KIND=dp)                                      :: voxel, weight(6), wsum, best, factor(3)
      REAL(KIND=dp), ALLOCATABLE                         :: flow(:, :)

      ng = SIZE(rho)
      natom = SIZE(seeds)
      CPASSERT(ng == PRODUCT(n))
      CPASSERT(ALL(n > 0) .AND. ALL(spacing > 0.0_dp))
      CPASSERT(SIZE(assigned) == natom)
      CPASSERT(ALL(SHAPE(atomic) == [2, natom]))
      IF (PRESENT(atom_grid)) THEN
         CPASSERT(natom > 0)
         CPASSERT(ALL(SHAPE(atom_grid) == [3, natom]))
      END IF
      voxel = PRODUCT(spacing)
      factor = 1.0_dp/spacing**2
      ALLOCATE (degree(ng), order(ng), flow(2, ng), attractor(ng), owner(ng))
      degree = 0
      attractor = 0
      owner = 0
      flow = 0.0_dp
      atomic = 0.0_dp
      assigned = .FALSE.
      other = 0.0_dp
      vacuum = 0.0_dp
      n_other = 0
      IF (PRESENT(n_extra)) n_extra = 0
      active = 0
      DO i = 1, ng
         IF (rho(i) <= cutoff) THEN
            vacuum = vacuum + [rho(i)*voxel, voxel]
            CYCLE
         END IF
         active = active + 1
         flow(:, i) = [rho(i)*voxel, voxel]
         CALL uphill(i, neighbor, weight)
         DO j = 1, 6
            IF (weight(j) == 0.0_dp) CYCLE
            k = neighbor(j)
            degree(k) = degree(k) + 1
         END DO
      END DO

      tail = 0
      DO i = 1, ng
         IF (rho(i) <= cutoff .OR. degree(i) /= 0) CYCLE
         tail = tail + 1
         order(tail) = i
      END DO

      ! The queue visits every active point after all its incoming contributions.
      head = 1
      DO WHILE (head <= tail)
         i = order(head)
         head = head + 1
         CALL uphill(i, neighbor, weight)
         wsum = SUM(weight)
         IF (wsum == 0.0_dp) THEN
            attractor(i) = i
            CYCLE
         END IF
         best = 0.0_dp
         DO j = 1, 6
            IF (weight(j) == 0.0_dp) CYCLE
            k = neighbor(j)
            IF (weight(j) > best) THEN
               best = weight(j)
               attractor(i) = k
            END IF
            flow(:, k) = flow(:, k) + (weight(j)/wsum)*flow(:, i)
            degree(k) = degree(k) - 1
            IF (degree(k) /= 0) CYCLE
            tail = tail + 1
            order(tail) = k
         END DO
      END DO
      CPASSERT(tail == active)

      DO head = tail, 1, -1
         i = order(head)
         attractor(i) = attractor(attractor(i))
      END DO
      DO atom = 1, natom
         i = seeds(atom)
         IF (i < 1 .OR. i > ng) CYCLE
         root = attractor(i)
         IF (root == 0) CYCLE
         IF (owner(root) /= 0) THEN
            owner(root) = -1
         ELSE
            owner(root) = atom
         END IF
      END DO
      DO i = 1, ng
         IF (attractor(i) /= i) CYCLE
         atom = owner(i)
         IF (PRESENT(n_extra)) THEN
            IF (atom <= 0) n_extra = n_extra + 1
         END IF
         IF (PRESENT(atom_grid)) THEN
            CALL add_to_nearest(i)
         ELSE IF (atom > 0) THEN
            assigned(atom) = .TRUE.
            atomic(:, atom) = flow(:, i)
         ELSE
            other = other + flow(:, i)
            n_other = n_other + 1
         END IF
      END DO

   CONTAINS

! **************************************************************************************************
!> \brief Add a basin to its nearest atom, sharing exact distance ties without double counting.
!> \param point grid index of the basin maximum
! **************************************************************************************************
      SUBROUTINE add_to_nearest(point)
      INTEGER, INTENT(IN)                                :: point

      INTEGER                                            :: a, coordinate(3), d, ties
      REAL(KIND=dp)                                      :: delta(3), distance2(natom), nearest, tol

         coordinate(1) = MODULO(point - 1, n(1))
         coordinate(2) = MODULO((point - 1)/n(1), n(2))
         coordinate(3) = ((point - 1)/n(1))/n(2)
         DO a = 1, natom
            delta = REAL(coordinate, dp) - atom_grid(:, a)
            DO d = 1, 3
               IF (periodic(d)) delta(d) = delta(d) - n(d)*ANINT(delta(d)/n(d))
            END DO
            distance2(a) = SUM((delta*spacing)**2)
         END DO
         nearest = MINVAL(distance2)
         tol = 16.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, nearest)
         ties = COUNT(distance2 <= nearest + tol)
         DO a = 1, natom
            IF (distance2(a) > nearest + tol) CYCLE
            atomic(:, a) = atomic(:, a) + flow(:, point)/ties
            assigned(a) = .TRUE.
         END DO
      END SUBROUTINE add_to_nearest

! **************************************************************************************************
!> \brief Find neighbors and their uphill flux weights (face area divided by distance).
!> \param point source grid index
!> \param next indices of neighbors along six faces
!> \param weights positive flux weights, up to a common voxel-volume factor
! **************************************************************************************************
      SUBROUTINE uphill(point, next, weights)
      INTEGER, INTENT(IN)                                :: point
      INTEGER, INTENT(OUT)                               :: next(6)
      REAL(KIND=dp), INTENT(OUT)                         :: weights(6)

      INTEGER                                            :: along, coord(3), d, dir, idx, offset, &
                                                            stride, TARGET

         coord(1) = MODULO(point - 1, n(1))
         coord(2) = MODULO((point - 1)/n(1), n(2))
         coord(3) = ((point - 1)/n(1))/n(2)
         stride = 1
         next = point
         weights = 0.0_dp
         DO d = 1, 3
            DO dir = -1, 1, 2
               idx = 2*d - (1 - dir)/2
               along = coord(d) + dir
               IF (periodic(d)) along = MODULO(along, n(d))
               IF (along < 0 .OR. along >= n(d)) CYCLE
               offset = stride*(along - coord(d))
               TARGET = point + offset
               next(idx) = TARGET
               IF (rho(TARGET) > rho(point)) weights(idx) = factor(d)*(rho(TARGET) - rho(point))
            END DO
            stride = stride*n(d)
         END DO
      END SUBROUTINE uphill
   END SUBROUTINE bader_integrate_grid
END MODULE qs_bader_basins

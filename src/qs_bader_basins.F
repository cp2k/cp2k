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
!>        Equal-density neighbours form a component. A topological traversal distributes
!>        charge and volume along strictly uphill faces between components.
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
!> \param extra optional electron population and volume of these maxima before assignment
!> \param labels optional strongest-uphill labels (atom index, 0 for vacuum, -1 if no unique owner)
! **************************************************************************************************
   SUBROUTINE bader_integrate_grid(rho, n, spacing, periodic, cutoff, seeds, atomic, assigned, other, vacuum, &
                                   n_other, atom_grid, n_extra, extra, labels)
      REAL(KIND=dp), INTENT(IN)                          :: rho(:), SPACING(3), cutoff
      INTEGER, INTENT(IN)                                :: n(3), seeds(:)
      LOGICAL, INTENT(IN)                                :: periodic(3)
      REAL(KIND=dp), INTENT(OUT)                         :: atomic(:, :), other(2), vacuum(2)
      LOGICAL, INTENT(OUT)                               :: assigned(:)
      INTEGER, INTENT(OUT)                               :: n_other
      REAL(KIND=dp), OPTIONAL, INTENT(IN)                :: atom_grid(:, :)
      INTEGER, OPTIONAL, INTENT(OUT)                     :: n_extra
      REAL(KIND=dp), OPTIONAL, INTENT(OUT)                :: extra(2)
      INTEGER, OPTIONAL, INTENT(OUT)                      :: labels(:)

      INTEGER :: head, tail, active, i, j, k, member, root, atom, natom, ng
      INTEGER :: neighbor(6), other_neighbors(6)
      INTEGER, ALLOCATABLE                              :: degree(:), order(:), attractor(:), next_member(:), component(:)
      REAL(KIND=dp)                                      :: voxel, weight(6), other_weights(6), wsum, best, factor(3)
      REAL(KIND=dp), ALLOCATABLE                         :: flow(:, :)

      ng = SIZE(rho)
      natom = SIZE(seeds)
      CPASSERT(ng == PRODUCT(n))
      CPASSERT(ALL(n > 0) .AND. ALL(spacing > 0.0_dp))
      CPASSERT(SIZE(assigned) == natom)
      CPASSERT(ALL(SHAPE(atomic) == [2, natom]))
      IF (PRESENT(labels)) THEN
         CPASSERT(SIZE(labels) == ng)
      END IF
      IF (PRESENT(atom_grid)) THEN
         CPASSERT(natom > 0)
         CPASSERT(ALL(SHAPE(atom_grid) == [3, natom]))
      END IF
      voxel = PRODUCT(spacing)
      factor = 1.0_dp/spacing**2
      ALLOCATE (degree(ng), order(ng), flow(2, ng), attractor(ng), next_member(ng), component(ng))
      degree = 0
      attractor = 0
      next_member = 0
      component(:) = MERGE(-1, 0, rho > cutoff)
      flow = 0.0_dp
      atomic = 0.0_dp
      assigned = .FALSE.
      other = 0.0_dp
      vacuum = 0.0_dp
      n_other = 0
      IF (PRESENT(n_extra)) n_extra = 0
      IF (PRESENT(extra)) extra = 0.0_dp

      ! Only visit the forward faces when joining equal-density points. Wrapping
      ! uses exactly the same boundary conditions as the uphill flux.
      DO i = 1, ng
         IF (component(i) == 0) CYCLE
         CALL face_neighbors(i, neighbor)
         DO j = 2, 6, 2
            k = neighbor(j)
            IF (k == i .OR. rho(k) /= rho(i)) CYCLE
            CALL join_components(i, k)
         END DO
      END DO
      DO i = 1, ng
         IF (component(i) <= 0) CYCLE
         root = find_root(i)
         component(i) = root
      END DO

      active = 0
      DO i = 1, ng
         IF (component(i) == 0) THEN
            vacuum = vacuum + [rho(i)*voxel, voxel]
            CYCLE
         ELSE IF (component(i) < 0) THEN
            component(i) = i
            active = active + 1
         ELSE
            root = component(i)
            next_member(i) = next_member(root)
            next_member(root) = i
         END IF
         root = component(i)
         flow(:, root) = flow(:, root) + [rho(i)*voxel, voxel]
      END DO
      DO i = 1, ng
         IF (component(i) == 0) CYCLE
         CALL uphill(i, neighbor, weight)
         DO j = 1, 6
            IF (weight(j) == 0.0_dp) CYCLE
            k = component(neighbor(j))
            degree(k) = degree(k) + 1
         END DO
      END DO

      tail = 0
      DO i = 1, ng
         IF (component(i) /= i .OR. degree(i) /= 0) CYCLE
         tail = tail + 1
         order(tail) = i
      END DO

      ! The queue visits each component after all its incoming contributions.
      head = 1
      DO WHILE (head <= tail)
         i = order(head)
         head = head + 1
         CALL uphill(i, neighbor, weight)
         wsum = SUM(weight)
         member = next_member(i)
         DO WHILE (member /= 0)
            CALL uphill(member, other_neighbors, other_weights)
            wsum = wsum + SUM(other_weights)
            member = next_member(member)
         END DO
         IF (wsum == 0.0_dp) THEN
            attractor(i) = i
            CYCLE
         END IF
         best = 0.0_dp
         member = i
         DO
            IF (member /= i) CALL uphill(member, neighbor, weight)
            DO j = 1, 6
               IF (weight(j) == 0.0_dp) CYCLE
               k = component(neighbor(j))
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
            member = next_member(member)
            IF (member == 0) EXIT
         END DO
      END DO
      CPASSERT(tail == active)

      DO head = tail, 1, -1
         i = order(head)
         attractor(i) = attractor(attractor(i))
      END DO
      ! Reuse the empty in-degree array to record the unique nuclear owner.
      DO atom = 1, natom
         i = seeds(atom)
         IF (i < 1 .OR. i > ng) CYCLE
         IF (component(i) == 0) CYCLE
         root = attractor(component(i))
         IF (degree(root) /= 0) THEN
            degree(root) = -1
         ELSE
            degree(root) = atom
         END IF
      END DO
      DO i = 1, ng
         IF (attractor(i) /= i) CYCLE
         atom = degree(i)
         IF (atom <= 0) THEN
            IF (PRESENT(n_extra)) n_extra = n_extra + 1
            IF (PRESENT(extra)) extra = extra + flow(:, i)
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
      IF (PRESENT(labels)) THEN
         labels = 0
         DO i = 1, ng
            IF (component(i) == 0) CYCLE
            atom = degree(attractor(component(i)))
            labels(i) = MERGE(atom, -1, atom > 0)
         END DO
      END IF

   CONTAINS

! **************************************************************************************************
!> \brief Find the representative of an active point and compress its parent path.
!> \param point active grid point
!> \return component representative
! **************************************************************************************************
      FUNCTION find_root(point) RESULT(rep)
      INTEGER, INTENT(IN)                                :: point
      INTEGER                                            :: rep

      INTEGER                                            :: node, parent

         rep = point
         DO WHILE (component(rep) > 0)
            rep = component(rep)
         END DO
         node = point
         DO WHILE (node /= rep)
            parent = component(node)
            component(node) = rep
            node = parent
         END DO
      END FUNCTION find_root

! **************************************************************************************************
!> \brief Join two flat components by size.
!> \param first first active grid point
!> \param second second active grid point
! **************************************************************************************************
      SUBROUTINE join_components(first, second)
      INTEGER, INTENT(IN)                                :: first, second

      INTEGER                                            :: left, right

         left = find_root(first)
         right = find_root(second)
         IF (left == right) RETURN
         IF (component(left) <= component(right)) THEN
            component(left) = component(left) + component(right)
            component(right) = left
         ELSE
            component(right) = component(right) + component(left)
            component(left) = right
         END IF
      END SUBROUTINE join_components

! **************************************************************************************************
!> \brief Assign by shortest distance to any point of the maximum, sharing ties.
!> \param point representative of the basin maximum
! **************************************************************************************************
      SUBROUTINE add_to_nearest(point)
      INTEGER, INTENT(IN)                                :: point

      INTEGER                                            :: a, coordinate(3), d, node, ties
      REAL(KIND=dp)                                      :: delta(3), distance2(natom), nearest, tol

         distance2 = HUGE(1.0_dp)
         node = point
         DO WHILE (node /= 0)
            coordinate = grid_coordinate(node)
            DO a = 1, natom
               delta = REAL(coordinate, dp) - atom_grid(:, a)
               DO d = 1, 3
                  IF (periodic(d)) delta(d) = delta(d) - n(d)*ANINT(delta(d)/n(d))
               END DO
               distance2(a) = MIN(distance2(a), SUM((delta*spacing)**2))
            END DO
            node = next_member(node)
         END DO
         nearest = MINVAL(distance2)
         tol = 16.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, nearest)
         ties = COUNT(distance2 <= nearest + tol)
         IF (PRESENT(labels)) degree(point) = -1
         DO a = 1, natom
            IF (distance2(a) > nearest + tol) CYCLE
            atomic(:, a) = atomic(:, a) + flow(:, point)/ties
            assigned(a) = .TRUE.
            IF (PRESENT(labels) .AND. ties == 1) degree(point) = a
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

      INTEGER                                            :: j

         CALL face_neighbors(point, next)
         weights = 0.0_dp
         DO j = 1, 6
            IF (rho(next(j)) > rho(point)) weights(j) = factor((j + 1)/2)*(rho(next(j)) - rho(point))
         END DO
      END SUBROUTINE uphill

! **************************************************************************************************
!> \brief Find face neighbours, using the point itself for a missing face.
!> \param point grid point
!> \param next neighbours, ordered by axis and negative/positive direction
! **************************************************************************************************
      SUBROUTINE face_neighbors(point, next)
      INTEGER, INTENT(IN)                                :: point
      INTEGER, INTENT(OUT)                               :: next(6)

      INTEGER                                            :: along, coord(3), d, dir, idx, stride

         coord = grid_coordinate(point)
         stride = 1
         next = point
         DO d = 1, 3
            DO dir = -1, 1, 2
               idx = 2*d - (1 - dir)/2
               along = coord(d) + dir
               IF (periodic(d)) along = MODULO(along, n(d))
               IF (along < 0 .OR. along >= n(d)) CYCLE
               next(idx) = point + stride*(along - coord(d))
            END DO
            stride = stride*n(d)
         END DO
      END SUBROUTINE face_neighbors

! **************************************************************************************************
!> \brief Decode a linear index into zero-based grid coordinates.
!> \param point grid point
!> \return Cartesian grid indices
! **************************************************************************************************
      FUNCTION grid_coordinate(point) RESULT(coord)
      INTEGER, INTENT(IN)                                :: point
      INTEGER                                            :: coord(3)

         coord(1) = MODULO(point - 1, n(1))
         coord(2) = MODULO((point - 1)/n(1), n(2))
         coord(3) = ((point - 1)/n(1))/n(2)
      END FUNCTION grid_coordinate
   END SUBROUTINE bader_integrate_grid
END MODULE qs_bader_basins

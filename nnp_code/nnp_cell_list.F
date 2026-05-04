!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Linked-cell + Verlet-skin neighbour-finder for the NNP descriptor.
!>
!>        Owns the persistent cell-list cache (positions, image pool, head/next
!>        chain, bin geometry) and the routines that build it, refresh it on
!>        skin reuse, and walk it to populate per-element neighbour buffers.
!>        Pulled out of nnp_acsf.F so the descriptor module only deals with
!>        ACSF computation and the cell-list machinery is testable in isolation.
!>
!>        Cache state lives on the parent nnp_type (nnp%cell_list_cache) so
!>        independent &NNP force_evals carry independent caches and the
!>        lifetime tracks nnp_env_release.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
MODULE nnp_cell_list

   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              real_to_scaled
   USE kinds,                           ONLY: dp
   USE nnp_environment_types,           ONLY: nnp_cell_list_cache_type,&
                                              nnp_env_get,&
                                              nnp_neighbor_type,&
                                              nnp_type
   USE nnp_neighbor_interface,          ONLY: nnp_neigh_grp_grow
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Verlet skin parameters: skin = MIN(0.5, 0.1 * cutoff). Chain is rebuilt
   ! when any atom drifts more than skin/2 from its rebuild-time position.
   REAL(KIND=dp), PRIVATE, PARAMETER :: default_verlet_skin = 0.5_dp
   REAL(KIND=dp), PRIVATE, PARAMETER :: rebuild_eps = 1.0E-12_dp

   PUBLIC :: nnp_prepare_cell_list_cache, &
             nnp_compute_neighbors_cell_list, &
             nnp_compute_pbc_copies

CONTAINS

! **************************************************************************************************
!> \brief Prepare the cell-list cache for the current force evaluation.
!>
!>        Wraps positions into the primary cell, decides whether the bin
!>        geometry needs a full rebuild (cell change, num_atoms change,
!>        cutoff change, or any atom drifted > skin/2), and either
!>        rebuilds via nnp_build_cell_list_cache or rebins the chain via
!>        nnp_cell_list_rebin_chain (skin reuse).
!>
!>        Must be called once per force eval, BEFORE any per-atom
!>        nnp_compute_neighbors_cell_list call.
!> \param nnp NNP environment whose cell_list_cache is to be (re)prepared from current positions.
!> \par History
!>      2026-04 created [dsharma]
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_prepare_cell_list_cache(nnp)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp

      INTEGER                                            :: handle, i
      LOGICAL                                            :: rebuild
      REAL(KIND=dp)                                      :: dr2_max, exact_cutoff, list_cutoff, skin
      REAL(KIND=dp), DIMENSION(3)                        :: dr
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset('nnp_acsf_cell_list_prepare', handle)

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      exact_cutoff = nnp%max_cut
      IF (exact_cutoff > 0.0_dp) THEN
         skin = MIN(default_verlet_skin, 0.1_dp*exact_cutoff)
      ELSE
         skin = 0.0_dp
      END IF
      list_cutoff = exact_cutoff + skin

      CALL nnp_cell_list_ensure_coord_buffers(nnp%cell_list_cache, nnp%num_atoms)

      ASSOCIATE (cache => nnp%cell_list_cache)
         DO i = 1, nnp%num_atoms
            cache%coord_primary(:, i) = pbc(nnp%coord(:, i), cell, .TRUE.)
            CALL real_to_scaled(cache%coord_scaled(:, i), &
                                cache%coord_primary(:, i), cell)
         END DO

         rebuild = .NOT. cache%initialized
         IF (.NOT. rebuild) rebuild = (cache%num_atoms /= nnp%num_atoms)
         IF (.NOT. rebuild) rebuild = (ABS(cache%exact_cutoff - exact_cutoff) > rebuild_eps)
         IF (.NOT. rebuild) rebuild = (ABS(cache%list_cutoff - list_cutoff) > rebuild_eps)
         IF (.NOT. rebuild) rebuild = (ABS(cache%verlet_skin - skin) > rebuild_eps)
         IF (.NOT. rebuild) rebuild = ANY(cache%perd /= cell%perd)
         IF (.NOT. rebuild) rebuild = (cache%orthorhombic .NEQV. cell%orthorhombic)
         IF (.NOT. rebuild) rebuild = ANY(ABS(cache%hmat - cell%hmat) > rebuild_eps)
         IF (.NOT. rebuild) rebuild = ANY(ABS(cache%h_inv - cell%h_inv) > rebuild_eps)

         IF ((.NOT. rebuild) .AND. (skin > 0.0_dp)) THEN
            dr2_max = 0.0_dp
            DO i = 1, nnp%num_atoms
               dr(:) = cache%coord_primary(:, i) - cache%ref_coord_primary(:, i)
               dr(:) = pbc(dr, cell)
               dr2_max = MAX(dr2_max, DOT_PRODUCT(dr, dr))
            END DO
            IF (dr2_max > 0.25_dp*skin**2) rebuild = .TRUE.
         END IF

         IF (rebuild) THEN
            CALL nnp_build_cell_list_cache(cache, cell, exact_cutoff, list_cutoff, skin)
         ELSE IF (skin > 0.0_dp .AND. cache%initialized) THEN
            ! Verlet-skin reuse: bin geometry (lower/bin_width/nbin/bin_span/
            ! image_pool) stays frozen, but the head/next chain MUST be
            ! rebinned from current positions every step. Diagnosed via
            ! brute-force parity walk: stale chain misses pairs at the
            ! very next force eval after a Verlet skip, because images
            ! binned at rebuild-time positions sit in bins outside the
            ! walker's bin_span around the central atom's current bin.
            ! Cost is O(n_images) integer ops, no allocation.
            CALL nnp_cell_list_rebin_chain(cache)
         END IF
      END ASSOCIATE

      CALL timestop(handle)

   END SUBROUTINE nnp_prepare_cell_list_cache

! **************************************************************************************************
!> \brief Ensure the persistent coord buffers match num_atoms.
!>        Only (re-)allocates on first call or when num_atoms changes.
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_ensure_coord_buffers(cache, num_atoms)

      TYPE(nnp_cell_list_cache_type), INTENT(INOUT)      :: cache
      INTEGER, INTENT(IN)                                :: num_atoms

      IF (ALLOCATED(cache%coord_primary)) THEN
         IF (SIZE(cache%coord_primary, 2) == num_atoms) RETURN
         DEALLOCATE (cache%coord_primary)
      END IF
      IF (ALLOCATED(cache%coord_scaled)) DEALLOCATE (cache%coord_scaled)
      IF (ALLOCATED(cache%ref_coord_primary)) DEALLOCATE (cache%ref_coord_primary)

      ALLOCATE (cache%coord_primary(3, num_atoms))
      ALLOCATE (cache%coord_scaled(3, num_atoms))
      ALLOCATE (cache%ref_coord_primary(3, num_atoms))

      ! Resizing implies stale bin geometry; force the next prepare to rebuild.
      cache%initialized = .FALSE.
      cache%num_atoms = num_atoms

   END SUBROUTINE nnp_cell_list_ensure_coord_buffers

! **************************************************************************************************
!> \brief Build the Cartesian linked-cell cache from the current positions.
!>
!>        Positions must already live in cache%coord_primary / coord_scaled --
!>        this routine reads them in place and resets ref_coord_primary for
!>        the Verlet skin check.
! **************************************************************************************************
   SUBROUTINE nnp_build_cell_list_cache(cache, cell, exact_cutoff, list_cutoff, skin)

      TYPE(nnp_cell_list_cache_type), INTENT(INOUT)      :: cache
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: exact_cutoff, list_cutoff, skin

      INTEGER                                            :: i, img, n_images, nx, ny, nz, tx, ty, tz
      INTEGER, DIMENSION(3)                              :: bin_index, shift
      REAL(KIND=dp)                                      :: padding, safe_cutoff
      REAL(KIND=dp), DIMENSION(3)                        :: max_ref, min_ref, range_xyz, ref_pos

      CALL nnp_cell_list_release_cell_structure(cache)

      cache%exact_cutoff = exact_cutoff
      cache%list_cutoff = list_cutoff
      cache%verlet_skin = skin
      cache%perd = cell%perd
      cache%orthorhombic = cell%orthorhombic
      cache%hmat = cell%hmat
      cache%h_inv = cell%h_inv

      cache%ref_coord_primary(:, :) = cache%coord_primary(:, :)

      CALL nnp_compute_pbc_copies(cache%exact_pbc_copies, cell, exact_cutoff)
      CALL nnp_compute_pbc_copies(cache%list_pbc_copies, cell, list_cutoff)
      ! The bin stencil clips out-of-range bins instead of wrapping (see
      ! nnp_compute_neighbors_cell_list below), so periodic wrap-around
      ! neighbours are only found through pre-replicated shifted images in
      ! the pool. On axes with L > 2*cutoff (any tile beyond H2O-64)
      ! list_pbc_copies would drop to 0 and strip those replicas, so force
      ! at least one ring on every periodic axis. Using MAX rather than +
      ! keeps w=64 bit-exact -- at small boxes list_pbc_copies is already
      ! >= cell%perd, so the pool size and bin_width are unchanged.
      cache%image_copies = MAX(cache%list_pbc_copies, cell%perd)

      nx = 2*cache%image_copies(1) + 1
      ny = 2*cache%image_copies(2) + 1
      nz = 2*cache%image_copies(3) + 1
      n_images = cache%num_atoms*nx*ny*nz
      cache%n_images = n_images

      ALLOCATE (cache%image_atom(n_images))
      ALLOCATE (cache%image_shift(3, n_images))
      ALLOCATE (cache%image_translation(3, n_images))

      min_ref(:) = HUGE(1.0_dp)
      max_ref(:) = -HUGE(1.0_dp)
      img = 0
      DO i = 1, cache%num_atoms
         DO tx = -cache%image_copies(1), cache%image_copies(1)
            DO ty = -cache%image_copies(2), cache%image_copies(2)
               DO tz = -cache%image_copies(3), cache%image_copies(3)
                  img = img + 1
                  shift = [tx, ty, tz]
                  cache%image_atom(img) = i
                  cache%image_shift(:, img) = shift(:)
                  cache%image_translation(:, img) = MATMUL(cell%hmat, REAL(shift, KIND=dp))
                  ref_pos(:) = cache%coord_primary(:, i) + cache%image_translation(:, img)
                  min_ref(:) = MIN(min_ref(:), ref_pos(:))
                  max_ref(:) = MAX(max_ref(:), ref_pos(:))
               END DO
            END DO
         END DO
      END DO

      padding = 0.5_dp*skin + rebuild_eps
      cache%lower(:) = min_ref(:) - padding
      cache%upper(:) = max_ref(:) + padding
      range_xyz(:) = cache%upper(:) - cache%lower(:)

      ! Bin width = list_cutoff (coarse grid, bin_span = 1, 27-bin walk).
      ! NOTE: a finer grid (bin_width ~ list_cutoff/2) is faster but changes
      ! the order in which neighbors are appended to the list, which changes
      ! the floating-point summation order in the descriptor. That tripped
      ! the H2O-64 regtest at the ~3e-7 relative level (FP reordering, not a
      ! correctness bug). Stick with the coarse grid until we add an explicit
      ! sort of the neighbor list by (image index, distance).
      safe_cutoff = MAX(list_cutoff, rebuild_eps)
      DO i = 1, 3
         IF (range_xyz(i) > rebuild_eps) THEN
            cache%nbin(i) = MAX(1, CEILING(range_xyz(i)/safe_cutoff))
            cache%bin_width(i) = range_xyz(i)/REAL(cache%nbin(i), KIND=dp)
            cache%bin_span(i) = MIN(cache%nbin(i) - 1, &
                                    CEILING(list_cutoff/MAX(cache%bin_width(i), rebuild_eps)))
         ELSE
            cache%nbin(i) = 1
            cache%bin_width(i) = 1.0_dp
            cache%bin_span(i) = 0
         END IF
      END DO

      cache%n_cells = PRODUCT(cache%nbin)
      ALLOCATE (cache%head(cache%n_cells))
      ALLOCATE (cache%next(n_images))
      cache%head(:) = 0
      cache%next(:) = 0

      DO img = 1, n_images
         i = cache%image_atom(img)
         ref_pos(:) = cache%coord_primary(:, i) + cache%image_translation(:, img)
         CALL nnp_cell_list_bin_from_position(cache, ref_pos, bin_index)
         cache%next(img) = cache%head(nnp_cell_list_linear_index(cache, bin_index))
         cache%head(nnp_cell_list_linear_index(cache, bin_index)) = img
      END DO

      cache%initialized = .TRUE.

   END SUBROUTINE nnp_build_cell_list_cache

! **************************************************************************************************
!> \brief Rebin the head/next linked-cell chain from current coord_primary while
!>        keeping the rest of the cell list (image pool, bin geometry,
!>        ref_coord_primary) frozen. Called when the Verlet-skin displacement
!>        check decides no atom has drifted enough to warrant a full rebuild --
!>        but the chain MUST still reflect current positions, because the walker
!>        centres on the atom's current bin and relies on the chain being current
!>        to find images within bin_span. Without this rebin the stencil walker
!>        silently misses pairs that have crossed bin boundaries since the last
!>        rebuild (verified at w=256 NVE: ~120 pairs/step missed under skin reuse,
!>        producing non-conservative forces; bit-exact fixed by rebinning every step).
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_rebin_chain(cache)

      TYPE(nnp_cell_list_cache_type), INTENT(INOUT)      :: cache

      INTEGER                                            :: i, img, lin_idx
      INTEGER, DIMENSION(3)                              :: bin_index
      REAL(KIND=dp), DIMENSION(3)                        :: ref_pos

      cache%head(:) = 0
      cache%next(:) = 0

      DO img = 1, cache%n_images
         i = cache%image_atom(img)
         ref_pos(:) = cache%coord_primary(:, i) + cache%image_translation(:, img)
         CALL nnp_cell_list_bin_from_position(cache, ref_pos, bin_index)
         lin_idx = nnp_cell_list_linear_index(cache, bin_index)
         cache%next(img) = cache%head(lin_idx)
         cache%head(lin_idx) = img
      END DO

   END SUBROUTINE nnp_cell_list_rebin_chain

! **************************************************************************************************
!> \brief Free only the cell-structure pieces of the cache (bins + image metadata).
!>        Coordinate buffers are persistent and owned by nnp_cell_list_ensure_coord_buffers --
!>        they must NOT be freed here, otherwise rebuilds would throw away the current positions
!>        we're about to re-bin.
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_release_cell_structure(cache)

      TYPE(nnp_cell_list_cache_type), INTENT(INOUT)      :: cache

      IF (ALLOCATED(cache%image_atom)) DEALLOCATE (cache%image_atom)
      IF (ALLOCATED(cache%head)) DEALLOCATE (cache%head)
      IF (ALLOCATED(cache%next)) DEALLOCATE (cache%next)
      IF (ALLOCATED(cache%image_shift)) DEALLOCATE (cache%image_shift)
      IF (ALLOCATED(cache%image_translation)) DEALLOCATE (cache%image_translation)

      cache%n_images = 0
      cache%n_cells = 1
      cache%exact_pbc_copies = 0
      cache%list_pbc_copies = 0
      cache%image_copies = 0
      cache%nbin = 1
      cache%bin_span = 0
      cache%lower = 0.0_dp
      cache%upper = 0.0_dp
      cache%bin_width = 1.0_dp

   END SUBROUTINE nnp_cell_list_release_cell_structure

! **************************************************************************************************
!> \brief Populate ACSF neighbour buffers for one central atom by walking the
!>        linked-cell stencil around its bin. Pushes (j, dr, r) triplets into
!>        the per-(species-pair) (ind, dist) slabs reserved by
!>        nnp_neighbor_interface_prepare; entries below the per-pair cutoff
!>        squared make it through. Skips own-bin self-pair (i==j with shift 0).
!> \param nnp       NNP environment with cell_list_cache and neighbor_interface_state ready.
!> \param neighbor  per-atom neighbour view that will be filled in place
!>                  (n_rad/n_ang1/n_ang2 zeroed by reset_neighbor before this call).
!> \param i         central-atom index in the global ordering.
!> \par History
!>      2026-04 created [dsharma]
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_compute_neighbors_cell_list(nnp, neighbor, i)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      TYPE(nnp_neighbor_type), INTENT(INOUT)             :: neighbor
      INTEGER, INTENT(IN)                                :: i

      INTEGER                                            :: bx, by, bz, img, ind, j, neighbor_ind, pair_slot, s
      INTEGER, DIMENSION(3)                              :: bin_index, current_shift
      REAL(KIND=dp)                                      :: norm
      REAL(KIND=dp), DIMENSION(3)                        :: center, dr, image_pos

      ASSOCIATE (cache => nnp%cell_list_cache, &
                 state => nnp%neighbor_interface_state)

         ind = nnp%ele_ind(i)
         center(:) = cache%coord_primary(:, i)
         CALL nnp_cell_list_bin_from_position(cache, center, bin_index)

         DO bx = MAX(1, bin_index(1) - cache%bin_span(1)), &
                 MIN(cache%nbin(1), bin_index(1) + cache%bin_span(1))
            DO by = MAX(1, bin_index(2) - cache%bin_span(2)), &
                    MIN(cache%nbin(2), bin_index(2) + cache%bin_span(2))
               DO bz = MAX(1, bin_index(3) - cache%bin_span(3)), &
                       MIN(cache%nbin(3), bin_index(3) + cache%bin_span(3))
                  img = cache%head(nnp_cell_list_linear_index(cache, [bx, by, bz]))
                  DO WHILE (img > 0)
                     j = cache%image_atom(img)
                     current_shift(:) = cache%image_shift(:, img)
                     IF (j == i .AND. ALL(current_shift == 0)) THEN
                        img = cache%next(img)
                        CYCLE
                     END IF
                     IF (.NOT. nnp_cell_list_exact_image_match(cache, i, j, current_shift)) THEN
                        img = cache%next(img)
                        CYCLE
                     END IF

                     image_pos(:) = cache%coord_primary(:, j) + cache%image_translation(:, img)
                     dr(:) = center(:) - image_pos(:)
                     norm = NORM2(dr(:))
                     neighbor_ind = nnp%ele_ind(j)

                     IF (norm < state%pair_map(ind, neighbor_ind)%max_relevant_cutoff) THEN
                        DO pair_slot = 1, state%pair_map(ind, neighbor_ind)%n_rad
                           s = state%pair_map(ind, neighbor_ind)%rad_groups(pair_slot)
                           IF (norm < nnp%rad(ind)%symfgrp(s)%cutoff) THEN
                              neighbor%n_rad(s) = neighbor%n_rad(s) + 1
                              IF (neighbor%n_rad(s) > neighbor%rad(s)%cap) &
                                 CALL nnp_neigh_grp_grow(neighbor%rad(s), neighbor%n_rad(s))
                              neighbor%rad(s)%ind(neighbor%n_rad(s)) = j
                              neighbor%rad(s)%dist(1, neighbor%n_rad(s)) = dr(1)
                              neighbor%rad(s)%dist(2, neighbor%n_rad(s)) = dr(2)
                              neighbor%rad(s)%dist(3, neighbor%n_rad(s)) = dr(3)
                              neighbor%rad(s)%dist(4, neighbor%n_rad(s)) = norm
                           END IF
                        END DO

                        DO pair_slot = 1, state%pair_map(ind, neighbor_ind)%n_ang1
                           s = state%pair_map(ind, neighbor_ind)%ang1_groups(pair_slot)
                           IF (norm < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                              neighbor%n_ang1(s) = neighbor%n_ang1(s) + 1
                              IF (neighbor%n_ang1(s) > neighbor%ang1(s)%cap) &
                                 CALL nnp_neigh_grp_grow(neighbor%ang1(s), neighbor%n_ang1(s))
                              neighbor%ang1(s)%ind(neighbor%n_ang1(s)) = j
                              neighbor%ang1(s)%dist(1, neighbor%n_ang1(s)) = dr(1)
                              neighbor%ang1(s)%dist(2, neighbor%n_ang1(s)) = dr(2)
                              neighbor%ang1(s)%dist(3, neighbor%n_ang1(s)) = dr(3)
                              neighbor%ang1(s)%dist(4, neighbor%n_ang1(s)) = norm
                           END IF
                        END DO

                        DO pair_slot = 1, state%pair_map(ind, neighbor_ind)%n_ang2
                           s = state%pair_map(ind, neighbor_ind)%ang2_groups(pair_slot)
                           IF (norm < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                              neighbor%n_ang2(s) = neighbor%n_ang2(s) + 1
                              IF (neighbor%n_ang2(s) > neighbor%ang2(s)%cap) &
                                 CALL nnp_neigh_grp_grow(neighbor%ang2(s), neighbor%n_ang2(s))
                              neighbor%ang2(s)%ind(neighbor%n_ang2(s)) = j
                              neighbor%ang2(s)%dist(1, neighbor%n_ang2(s)) = dr(1)
                              neighbor%ang2(s)%dist(2, neighbor%n_ang2(s)) = dr(2)
                              neighbor%ang2(s)%dist(3, neighbor%n_ang2(s)) = dr(3)
                              neighbor%ang2(s)%dist(4, neighbor%n_ang2(s)) = norm
                           END IF
                        END DO
                     END IF

                     img = cache%next(img)
                  END DO
               END DO
            END DO
         END DO

      END ASSOCIATE

   END SUBROUTINE nnp_compute_neighbors_cell_list

! **************************************************************************************************
!> \brief Check whether a ghost image corresponds to one of the exact pbc copies.
! **************************************************************************************************
   LOGICAL FUNCTION nnp_cell_list_exact_image_match(cache, i, j, image_shift)

      TYPE(nnp_cell_list_cache_type), INTENT(IN)         :: cache
      INTEGER, INTENT(IN)                                :: i, j
      INTEGER, DIMENSION(3), INTENT(IN)                  :: image_shift

      INTEGER                                            :: d, minimal_shift

      nnp_cell_list_exact_image_match = .TRUE.
      DO d = 1, 3
         IF (cache%perd(d) == 1) THEN
            minimal_shift = NINT(cache%coord_scaled(d, i) - cache%coord_scaled(d, j))
            IF (ABS(minimal_shift - image_shift(d)) > cache%exact_pbc_copies(d)) THEN
               nnp_cell_list_exact_image_match = .FALSE.
               RETURN
            END IF
         ELSE IF (image_shift(d) /= 0) THEN
            nnp_cell_list_exact_image_match = .FALSE.
            RETURN
         END IF
      END DO

   END FUNCTION nnp_cell_list_exact_image_match

! **************************************************************************************************
!> \brief Convert Cartesian position to linked-cell bin coordinates.
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_bin_from_position(cache, pos, bin_index)

      TYPE(nnp_cell_list_cache_type), INTENT(IN)         :: cache
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: pos
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: bin_index

      INTEGER                                            :: d
      REAL(KIND=dp)                                      :: scaled

      DO d = 1, 3
         scaled = (pos(d) - cache%lower(d))/MAX(cache%bin_width(d), rebuild_eps)
         bin_index(d) = INT(scaled) + 1
         bin_index(d) = MAX(1, MIN(cache%nbin(d), bin_index(d)))
      END DO

   END SUBROUTINE nnp_cell_list_bin_from_position

! **************************************************************************************************
!> \brief Flatten 3D bin indices into the head/next storage index.
! **************************************************************************************************
   INTEGER FUNCTION nnp_cell_list_linear_index(cache, bin_index)

      TYPE(nnp_cell_list_cache_type), INTENT(IN)         :: cache
      INTEGER, DIMENSION(3), INTENT(IN)                  :: bin_index

      nnp_cell_list_linear_index = bin_index(1) + cache%nbin(1)* &
                                   ((bin_index(2) - 1) + cache%nbin(2)*(bin_index(3) - 1))

   END FUNCTION nnp_cell_list_linear_index

! **************************************************************************************************
!> \brief Compute the number of PBC image rings needed along each lattice axis to
!>        cover sphere of radius `cutoff` from any atom in the primary cell.
!>        Uses the projection of the cell vectors onto the plane normals so the
!>        result is correct for triclinic boxes; multiplies by cell%perd at the
!>        end so non-periodic axes always return 0 image rings.
!> \param pbc_copies  output: number of image rings on each of the 3 lattice axes.
!> \param cell        cell whose hmat/h_inv/perd/deth are read.
!> \param cutoff      sphere radius (Bohr) to cover.
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_compute_pbc_copies(pbc_copies, cell, cutoff)

      INTEGER, DIMENSION(3), INTENT(OUT)                 :: pbc_copies
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      REAL(KIND=dp)                                      :: proja, projb, projc
      REAL(KIND=dp), DIMENSION(3)                        :: axb, axc, bxc

      ! A degenerate cell (deth <= 0) makes the NORM2 normalisations below
      ! divide by zero and the DO WHILE projection loops never terminate.
      CPASSERT(cell%deth > 0.0_dp)

      axb(1) = cell%hmat(2, 1)*cell%hmat(3, 2) - cell%hmat(3, 1)*cell%hmat(2, 2)
      axb(2) = cell%hmat(3, 1)*cell%hmat(1, 2) - cell%hmat(1, 1)*cell%hmat(3, 2)
      axb(3) = cell%hmat(1, 1)*cell%hmat(2, 2) - cell%hmat(2, 1)*cell%hmat(1, 2)
      axb(:) = axb(:)/NORM2(axb(:))

      axc(1) = cell%hmat(2, 1)*cell%hmat(3, 3) - cell%hmat(3, 1)*cell%hmat(2, 3)
      axc(2) = cell%hmat(3, 1)*cell%hmat(1, 3) - cell%hmat(1, 1)*cell%hmat(3, 3)
      axc(3) = cell%hmat(1, 1)*cell%hmat(2, 3) - cell%hmat(2, 1)*cell%hmat(1, 3)
      axc(:) = axc(:)/NORM2(axc(:))

      bxc(1) = cell%hmat(2, 2)*cell%hmat(3, 3) - cell%hmat(3, 2)*cell%hmat(2, 3)
      bxc(2) = cell%hmat(3, 2)*cell%hmat(1, 3) - cell%hmat(1, 2)*cell%hmat(3, 3)
      bxc(3) = cell%hmat(1, 2)*cell%hmat(2, 3) - cell%hmat(2, 2)*cell%hmat(1, 3)
      bxc(:) = bxc(:)/NORM2(bxc(:))

      proja = ABS(SUM(cell%hmat(:, 1)*bxc(:)))*0.5_dp
      projb = ABS(SUM(cell%hmat(:, 2)*axc(:)))*0.5_dp
      projc = ABS(SUM(cell%hmat(:, 3)*axb(:)))*0.5_dp

      pbc_copies(:) = 0
      DO WHILE ((pbc_copies(1) + 1)*proja <= cutoff)
         pbc_copies(1) = pbc_copies(1) + 1
      END DO
      DO WHILE ((pbc_copies(2) + 1)*projb <= cutoff)
         pbc_copies(2) = pbc_copies(2) + 1
      END DO
      DO WHILE ((pbc_copies(3) + 1)*projc <= cutoff)
         pbc_copies(3) = pbc_copies(3) + 1
      END DO
      pbc_copies(:) = pbc_copies(:)*cell%perd(:)

   END SUBROUTINE nnp_compute_pbc_copies

END MODULE nnp_cell_list

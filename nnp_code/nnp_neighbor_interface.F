!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Per-nnp persistent neighbour-interface state for the NNP hot path.
!>         Mirrors the LAMMPS idea of separating neighbour bookkeeping from the
!>         ACSF/NN loops: species-pair routing is precomputed once and per-element
!>         work buffers are reused.
!>
!>         State lives on the parent nnp_type so independent &NNP
!>         force_evals (committee / MIX) carry independent caches and the
!>         lifetime tracks nnp_env_release. Not safe across an OMP region
!>         around the per-atom loop -- the OMP that exists today is inside
!>         the angular kernel only.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
MODULE nnp_neighbor_interface

   USE kinds,                           ONLY: dp
   USE nnp_environment_types,           ONLY: nnp_dGdr_grp_type,&
                                              nnp_neigh_grp_type,&
                                              nnp_neighbor_interface_state_release,&
                                              nnp_neighbor_workspace_type,&
                                              nnp_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: nnp_neighbor_interface_prepare, &
             nnp_neighbor_interface_reset_neighbor, &
             nnp_grp_grow_dGdr, &
             nnp_neigh_grp_grow, &
             nnp_workspace_grow_caches

CONTAINS

! **************************************************************************************************
!> \brief Ensure pair-routing metadata and reusable workspaces are ready for the current NNP model.
!> \param nnp NNP environment whose neighbor_interface_state will be (re)built if needed.
!> \par History
!>      2026-04 created [dsharma]
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_prepare(nnp)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      LOGICAL                                            :: rebuild

      rebuild = .NOT. nnp%neighbor_interface_state%initialized
      IF (.NOT. rebuild) CALL nnp_neighbor_interface_needs_rebuild(nnp, rebuild)

      IF (rebuild) THEN
         CALL nnp_neighbor_interface_state_release(nnp%neighbor_interface_state)
         CALL nnp_neighbor_interface_build_pair_maps(nnp)
      END IF

   END SUBROUTINE nnp_neighbor_interface_prepare

! **************************************************************************************************
!> \brief Reset per-group neighbour counters for one central element before refilling
!>        the reusable buffers. Zeroes n_rad/n_ang1/n_ang2; the (ind, dist) slabs
!>        are kept allocated and overwritten in place by the next push.
!> \param nnp NNP environment whose neighbor_interface_state will have its counters cleared.
!> \param ind central-element index (1..nnp%n_ele) whose workspace is reset.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_reset_neighbor(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      nnp%neighbor_interface_state%workspace(ind)%neighbor%n_rad(:) = 0
      nnp%neighbor_interface_state%workspace(ind)%neighbor%n_ang1(:) = 0
      nnp%neighbor_interface_state%workspace(ind)%neighbor%n_ang2(:) = 0

   END SUBROUTINE nnp_neighbor_interface_reset_neighbor

! **************************************************************************************************
!> \brief Check whether the persistent state matches the current NNP model.
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_needs_rebuild(nnp, rebuild)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      LOGICAL, INTENT(OUT)                               :: rebuild

      INTEGER                                            :: i

      rebuild = .FALSE.
      IF (nnp%neighbor_interface_state%n_ele /= nnp%n_ele) THEN
         rebuild = .TRUE.
         RETURN
      END IF
      IF (.NOT. ALLOCATED(nnp%neighbor_interface_state%n_rad)) THEN
         rebuild = .TRUE.
         RETURN
      END IF

      DO i = 1, nnp%n_ele
         IF (nnp%neighbor_interface_state%n_rad(i) /= nnp%n_rad(i) .OR. &
             nnp%neighbor_interface_state%n_ang(i) /= nnp%n_ang(i) .OR. &
             nnp%neighbor_interface_state%n_radgrp(i) /= nnp%rad(i)%n_symfgrp .OR. &
             nnp%neighbor_interface_state%n_anggrp(i) /= nnp%ang(i)%n_symfgrp) THEN
            rebuild = .TRUE.
            RETURN
         END IF
      END DO

   END SUBROUTINE nnp_neighbor_interface_needs_rebuild

! **************************************************************************************************
!> \brief Build species-pair routing tables and initialize reusable workspaces.
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_build_pair_maps(nnp)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: i

      ASSOCIATE (state => nnp%neighbor_interface_state)
         state%n_ele = nnp%n_ele
         ALLOCATE (state%n_rad(nnp%n_ele))
         ALLOCATE (state%n_ang(nnp%n_ele))
         ALLOCATE (state%n_radgrp(nnp%n_ele))
         ALLOCATE (state%n_anggrp(nnp%n_ele))
         ALLOCATE (state%pair_map(nnp%n_ele, nnp%n_ele))
         ALLOCATE (state%workspace(nnp%n_ele))

         DO i = 1, nnp%n_ele
            state%n_rad(i) = nnp%n_rad(i)
            state%n_ang(i) = nnp%n_ang(i)
            state%n_radgrp(i) = nnp%rad(i)%n_symfgrp
            state%n_anggrp(i) = nnp%ang(i)%n_symfgrp
         END DO

         DO i = 1, nnp%n_ele
            CALL nnp_neighbor_interface_build_pair_map_for_element(nnp, i)
            CALL nnp_neighbor_interface_init_workspace_metadata(nnp, i)
         END DO

         state%initialized = .TRUE.
      END ASSOCIATE

   END SUBROUTINE nnp_neighbor_interface_build_pair_maps

! **************************************************************************************************
!> \brief Build all species-pair routes for one central element.
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_build_pair_map_for_element(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      INTEGER                                            :: neighbor_ind, idx, s

      DO neighbor_ind = 1, nnp%n_ele
         ASSOCIATE (pair_map => nnp%neighbor_interface_state%pair_map(ind, neighbor_ind))
            pair_map%n_rad = 0
            pair_map%n_ang1 = 0
            pair_map%n_ang2 = 0
            pair_map%max_relevant_cutoff = 0.0_dp

            DO s = 1, nnp%rad(ind)%n_symfgrp
               IF (nnp%rad(ind)%symfgrp(s)%ele_ind(1) == neighbor_ind) pair_map%n_rad = pair_map%n_rad + 1
            END DO
            DO s = 1, nnp%ang(ind)%n_symfgrp
               IF (nnp%ang(ind)%symfgrp(s)%ele_ind(1) == neighbor_ind) pair_map%n_ang1 = pair_map%n_ang1 + 1
               IF (nnp%ang(ind)%symfgrp(s)%ele_ind(2) == neighbor_ind) pair_map%n_ang2 = pair_map%n_ang2 + 1
            END DO

            ALLOCATE (pair_map%rad_groups(MAX(1, pair_map%n_rad)))
            ALLOCATE (pair_map%ang1_groups(MAX(1, pair_map%n_ang1)))
            ALLOCATE (pair_map%ang2_groups(MAX(1, pair_map%n_ang2)))

            idx = 0
            DO s = 1, nnp%rad(ind)%n_symfgrp
               IF (nnp%rad(ind)%symfgrp(s)%ele_ind(1) == neighbor_ind) THEN
                  idx = idx + 1
                  pair_map%rad_groups(idx) = s
                  pair_map%max_relevant_cutoff = MAX(pair_map%max_relevant_cutoff, nnp%rad(ind)%symfgrp(s)%cutoff)
               END IF
            END DO

            idx = 0
            DO s = 1, nnp%ang(ind)%n_symfgrp
               IF (nnp%ang(ind)%symfgrp(s)%ele_ind(1) == neighbor_ind) THEN
                  idx = idx + 1
                  pair_map%ang1_groups(idx) = s
                  pair_map%max_relevant_cutoff = MAX(pair_map%max_relevant_cutoff, nnp%ang(ind)%symfgrp(s)%cutoff)
               END IF
            END DO

            idx = 0
            DO s = 1, nnp%ang(ind)%n_symfgrp
               IF (nnp%ang(ind)%symfgrp(s)%ele_ind(2) == neighbor_ind) THEN
                  idx = idx + 1
                  pair_map%ang2_groups(idx) = s
                  pair_map%max_relevant_cutoff = MAX(pair_map%max_relevant_cutoff, nnp%ang(ind)%symfgrp(s)%cutoff)
               END IF
            END DO
         END ASSOCIATE
      END DO

   END SUBROUTINE nnp_neighbor_interface_build_pair_map_for_element

! **************************************************************************************************
!> \brief Cache max scratch sizes for one central element.
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_init_workspace_metadata(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      INTEGER                                            :: s

      ASSOCIATE (workspace => nnp%neighbor_interface_state%workspace(ind))
         workspace%max_rad_symf = 0
         workspace%max_ang_symf = 0
         workspace%n_input_nodes = nnp%n_rad(ind) + nnp%n_ang(ind)

         DO s = 1, nnp%rad(ind)%n_symfgrp
            workspace%max_rad_symf = MAX(workspace%max_rad_symf, nnp%rad(ind)%symfgrp(s)%n_symf)
         END DO
         DO s = 1, nnp%ang(ind)%n_symfgrp
            workspace%max_ang_symf = MAX(workspace%max_ang_symf, nnp%ang(ind)%symfgrp(s)%n_symf)
         END DO

         ! Per-(SF) scratch reused inside nnp_calc_rad / nnp_calc_ang.
         ! Sized by max_*_symf only -- no dependence on the candidate pool.
         IF (ALLOCATED(workspace%radial_sym)) DEALLOCATE (workspace%radial_sym)
         IF (ALLOCATED(workspace%radial_force)) DEALLOCATE (workspace%radial_force)
         IF (ALLOCATED(workspace%angular_sym)) DEALLOCATE (workspace%angular_sym)
         IF (ALLOCATED(workspace%angular_force)) DEALLOCATE (workspace%angular_force)
         ALLOCATE (workspace%radial_sym(MAX(1, workspace%max_rad_symf)))
         ALLOCATE (workspace%radial_force(3, MAX(1, workspace%max_rad_symf)))
         ALLOCATE (workspace%angular_sym(MAX(1, workspace%max_ang_symf)))
         ALLOCATE (workspace%angular_force(3, 3, MAX(1, workspace%max_ang_symf)))

         ! 1D angular cutoff caches start at a small floor and grow lazily
         ! (see nnp_workspace_grow_caches) to the peak per-element angular
         ! neighbor count, NOT the image-pool capacity.
         IF (ALLOCATED(workspace%fc_cache1)) DEALLOCATE (workspace%fc_cache1)
         IF (ALLOCATED(workspace%dfc_cache1)) DEALLOCATE (workspace%dfc_cache1)
         IF (ALLOCATED(workspace%fc_cache2)) DEALLOCATE (workspace%fc_cache2)
         IF (ALLOCATED(workspace%dfc_cache2)) DEALLOCATE (workspace%dfc_cache2)
         workspace%cache_cap = 8
         ALLOCATE (workspace%fc_cache1(workspace%cache_cap))
         ALLOCATE (workspace%dfc_cache1(workspace%cache_cap))
         ALLOCATE (workspace%fc_cache2(workspace%cache_cap))
         ALLOCATE (workspace%dfc_cache2(workspace%cache_cap))

         ! self_dGdr is sized once at metadata init -- it does not depend on
         ! the linked-cell candidate pool capacity.
         IF (ALLOCATED(workspace%self_dGdr)) DEALLOCATE (workspace%self_dGdr)
         ALLOCATE (workspace%self_dGdr(3, MAX(1, workspace%n_input_nodes)))

         ! Per-group neighbor counters and per-group dense (ind, dist)
         ! containers. n_rad/n_ang1/n_ang2 are sized to n_symfgrp here
         ! (independent of the image-pool capacity); the per-group
         ! ind/dist slabs inside neighbor%rad/ang1/ang2 are grown lazily
         ! by nnp_neigh_grp_grow when neighbors are pushed.
         CALL nnp_release_neighbor_local(workspace)
         ALLOCATE (workspace%neighbor%n_rad(MAX(1, nnp%rad(ind)%n_symfgrp)))
         ALLOCATE (workspace%neighbor%n_ang1(MAX(1, nnp%ang(ind)%n_symfgrp)))
         ALLOCATE (workspace%neighbor%n_ang2(MAX(1, nnp%ang(ind)%n_symfgrp)))
         ALLOCATE (workspace%neighbor%rad(MAX(1, nnp%rad(ind)%n_symfgrp)))
         ALLOCATE (workspace%neighbor%ang1(MAX(1, nnp%ang(ind)%n_symfgrp)))
         ALLOCATE (workspace%neighbor%ang2(MAX(1, nnp%ang(ind)%n_symfgrp)))
         workspace%neighbor%n_rad(:) = 0
         workspace%neighbor%n_ang1(:) = 0
         workspace%neighbor%n_ang2(:) = 0
         workspace%neighbor%pbc_copies = 0

         ! Per-(group) dG/dr buffers. Container arrays are sized once here
         ! to n_symfgrp, and each entry's n_symf is recorded so subsequent
         ! grow calls can size the dense (3, n_symf, cap) slab correctly.
         ! The actual %data is left unallocated -- nnp_grp_grow_dGdr is
         ! called from nnp_calc_acsf the first time a real neighbor count
         ! is observed.
         CALL nnp_release_dGdr_grp_array(workspace%dGdr_rad)
         CALL nnp_release_dGdr_grp_array(workspace%dGdr_ang_jj)
         CALL nnp_release_dGdr_grp_array(workspace%dGdr_ang_kk)
         ALLOCATE (workspace%dGdr_rad(MAX(1, nnp%rad(ind)%n_symfgrp)))
         ALLOCATE (workspace%dGdr_ang_jj(MAX(1, nnp%ang(ind)%n_symfgrp)))
         ALLOCATE (workspace%dGdr_ang_kk(MAX(1, nnp%ang(ind)%n_symfgrp)))
         DO s = 1, nnp%rad(ind)%n_symfgrp
            workspace%dGdr_rad(s)%n_symf = nnp%rad(ind)%symfgrp(s)%n_symf
            workspace%dGdr_rad(s)%cap = 0
         END DO
         DO s = 1, nnp%ang(ind)%n_symfgrp
            workspace%dGdr_ang_jj(s)%n_symf = nnp%ang(ind)%symfgrp(s)%n_symf
            workspace%dGdr_ang_jj(s)%cap = 0
            workspace%dGdr_ang_kk(s)%n_symf = nnp%ang(ind)%symfgrp(s)%n_symf
            workspace%dGdr_ang_kk(s)%cap = 0
         END DO
      END ASSOCIATE

   END SUBROUTINE nnp_neighbor_interface_init_workspace_metadata

! **************************************************************************************************
!> \brief Release just the workspace's neighbor%(rad,ang1,ang2) per-group slabs
!>        before re-allocating them. Smaller-scope than nnp_workspace_release
!>        which also tears down the persistent caches we want to keep across
!>        rebuilds within the same workspace identity.
! **************************************************************************************************
   SUBROUTINE nnp_release_neighbor_local(workspace)

      TYPE(nnp_neighbor_workspace_type), INTENT(INOUT)   :: workspace

      INTEGER                                            :: s

      IF (ALLOCATED(workspace%neighbor%rad)) THEN
         DO s = 1, SIZE(workspace%neighbor%rad)
            IF (ALLOCATED(workspace%neighbor%rad(s)%ind)) DEALLOCATE (workspace%neighbor%rad(s)%ind)
            IF (ALLOCATED(workspace%neighbor%rad(s)%dist)) DEALLOCATE (workspace%neighbor%rad(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%rad)
      END IF
      IF (ALLOCATED(workspace%neighbor%ang1)) THEN
         DO s = 1, SIZE(workspace%neighbor%ang1)
            IF (ALLOCATED(workspace%neighbor%ang1(s)%ind)) DEALLOCATE (workspace%neighbor%ang1(s)%ind)
            IF (ALLOCATED(workspace%neighbor%ang1(s)%dist)) DEALLOCATE (workspace%neighbor%ang1(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%ang1)
      END IF
      IF (ALLOCATED(workspace%neighbor%ang2)) THEN
         DO s = 1, SIZE(workspace%neighbor%ang2)
            IF (ALLOCATED(workspace%neighbor%ang2(s)%ind)) DEALLOCATE (workspace%neighbor%ang2(s)%ind)
            IF (ALLOCATED(workspace%neighbor%ang2(s)%dist)) DEALLOCATE (workspace%neighbor%ang2(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%ang2)
      END IF
      IF (ALLOCATED(workspace%neighbor%n_rad)) DEALLOCATE (workspace%neighbor%n_rad)
      IF (ALLOCATED(workspace%neighbor%n_ang1)) DEALLOCATE (workspace%neighbor%n_ang1)
      IF (ALLOCATED(workspace%neighbor%n_ang2)) DEALLOCATE (workspace%neighbor%n_ang2)
      workspace%neighbor%pbc_copies = -1

   END SUBROUTINE nnp_release_neighbor_local

! **************************************************************************************************
!> \brief Release every per-group dG/dr buffer in a container array, then deallocate the container.
! **************************************************************************************************
   SUBROUTINE nnp_release_dGdr_grp_array(grps)

      TYPE(nnp_dGdr_grp_type), ALLOCATABLE, INTENT(INOUT) :: grps(:)

      INTEGER                                            :: s

      IF (.NOT. ALLOCATED(grps)) RETURN
      DO s = 1, SIZE(grps)
         IF (ALLOCATED(grps(s)%data)) DEALLOCATE (grps(s)%data)
         grps(s)%cap = 0
         grps(s)%n_symf = 0
      END DO
      DEALLOCATE (grps)

   END SUBROUTINE nnp_release_dGdr_grp_array

! **************************************************************************************************
!> \brief Ensure a per-group dG/dr buffer can hold n_needed neighbors.
!>        Uses 1.5x growth (with a small additive floor) so reallocations
!>        are amortized to O(1) over an MD trajectory: after a few steps
!>        the slab settles at the true peak and stops growing.
!>        Guarantees that grp%data is allocated on return (with at least
!>        cap >= 1) so callers can safely ASSOCIATE-bind it inside loops
!>        whose trip count happens to be zero.
!> \param grp        per-element dGdr group whose %data slab will be reallocated if undersized.
!> \param n_needed   minimum required leading-dim capacity (third dim of grp%data).
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_grp_grow_dGdr(grp, n_needed)

      TYPE(nnp_dGdr_grp_type), INTENT(INOUT)             :: grp
      INTEGER, INTENT(IN)                                :: n_needed

      INTEGER                                            :: new_cap

      IF (ALLOCATED(grp%data) .AND. grp%cap >= n_needed) RETURN
      IF (ALLOCATED(grp%data)) DEALLOCATE (grp%data)
      new_cap = MAX(MAX(1, n_needed), INT(grp%cap*1.5_dp) + 8)
      ALLOCATE (grp%data(3, MAX(1, grp%n_symf), new_cap))
      grp%cap = new_cap

   END SUBROUTINE nnp_grp_grow_dGdr

! **************************************************************************************************
!> \brief Ensure a per-group (ind, dist) neighbor buffer can hold n_needed
!>        entries. Unlike nnp_grp_grow_dGdr (which is called once per group
!>        before its accumulator pass), this is called from inside the
!>        linked-cell push loop, so any already-pushed entries MUST be
!>        preserved on grow. We do that with MOVE_ALLOC over a fresh slab.
!>        Uses 1.5x growth so reallocations amortize to O(1) per group per
!>        MD trajectory; the slab settles at the true peak after the first
!>        few steps.
!> \param grp        per-(species-pair) neighbour group whose (ind, dist) buffers grow.
!> \param n_needed   minimum required entry count.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_neigh_grp_grow(grp, n_needed)

      TYPE(nnp_neigh_grp_type), INTENT(INOUT)            :: grp
      INTEGER, INTENT(IN)                                :: n_needed

      INTEGER                                            :: new_cap, n_old
      INTEGER, ALLOCATABLE                               :: new_ind(:)
      REAL(KIND=dp), ALLOCATABLE                         :: new_dist(:, :)

      IF (ALLOCATED(grp%dist) .AND. grp%cap >= n_needed) RETURN

      new_cap = MAX(MAX(8, n_needed), INT(grp%cap*1.5_dp) + 8)

      IF (ALLOCATED(grp%dist)) THEN
         n_old = grp%cap
         ALLOCATE (new_dist(4, new_cap))
         ALLOCATE (new_ind(new_cap))
         IF (n_old > 0) THEN
            new_dist(:, 1:n_old) = grp%dist(:, 1:n_old)
            new_ind(1:n_old) = grp%ind(1:n_old)
         END IF
         CALL MOVE_ALLOC(new_dist, grp%dist)
         CALL MOVE_ALLOC(new_ind, grp%ind)
      ELSE
         ALLOCATE (grp%dist(4, new_cap))
         ALLOCATE (grp%ind(new_cap))
      END IF
      grp%cap = new_cap

   END SUBROUTINE nnp_neigh_grp_grow

! **************************************************************************************************
!> \brief Ensure all four per-element angular cutoff caches can hold at
!>        least n_needed entries. Caches are simple 1D scratch reused
!>        across angular groups within one nnp_calc_acsf call, so the
!>        peak is MAX_s(n_ang1(s)) and MAX_s(n_ang2(s)) for the current
!>        central atom. Grown lazily with 1.5x growth -- after a handful
!>        of MD steps the slab settles at the true element-wide peak.
!> \param workspace  per-element workspace whose fc_cache1/dfc_cache1/fc_cache2/dfc_cache2 grow.
!> \param n_needed   minimum required cache length.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_workspace_grow_caches(workspace, n_needed)

      TYPE(nnp_neighbor_workspace_type), INTENT(INOUT)   :: workspace
      INTEGER, INTENT(IN)                                :: n_needed

      INTEGER                                            :: new_cap

      IF (workspace%cache_cap >= n_needed) RETURN

      new_cap = MAX(MAX(8, n_needed), INT(workspace%cache_cap*1.5_dp) + 8)
      IF (ALLOCATED(workspace%fc_cache1)) DEALLOCATE (workspace%fc_cache1)
      IF (ALLOCATED(workspace%dfc_cache1)) DEALLOCATE (workspace%dfc_cache1)
      IF (ALLOCATED(workspace%fc_cache2)) DEALLOCATE (workspace%fc_cache2)
      IF (ALLOCATED(workspace%dfc_cache2)) DEALLOCATE (workspace%dfc_cache2)
      ALLOCATE (workspace%fc_cache1(new_cap))
      ALLOCATE (workspace%dfc_cache1(new_cap))
      ALLOCATE (workspace%fc_cache2(new_cap))
      ALLOCATE (workspace%dfc_cache2(new_cap))
      workspace%cache_cap = new_cap

   END SUBROUTINE nnp_workspace_grow_caches

END MODULE nnp_neighbor_interface

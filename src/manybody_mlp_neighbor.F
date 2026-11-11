!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Edge-list construction shared by the message-passing MLP drivers.
!>
!>        A message-passing model needs a full (not half) neighbour list expressed as
!>        edges: pairs of atoms within the model cutoff, each tagged with the periodic
!>        cell shift that maps the neighbour into the central atom's image. Every such
!>        driver walks the same FIST neighbour list in the same way, so the traversal
!>        and the cross-rank exchange live here; the callers keep their own index
!>        conventions and integer kinds, which is where they genuinely differ (the
!>        TorchScript path in manybody_e3nn needs 64-bit global indices and real-valued
!>        shifts, the symmetrix path in manybody_symmetrix needs 32-bit packed indices).
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE manybody_mlp_neighbor

   USE cell_types,                      ONLY: cell_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: pos_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                              pair_potential_single_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mlp_collect_edges, mlp_allgather_edges

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_mlp_neighbor'

CONTAINS

! **************************************************************************************************
!> \brief Collect this rank's edges from the FIST neighbour list.
!>
!>        Returns global atom indices and integer cell shifts; the caller converts to
!>        whatever its model interface expects. Pairs are kept when both kinds are
!>        flagged in kind_included, the kind pair carries target_pot_type, and the
!>        distance is within the matching entry of cutsq.
!> \param nonbonded the FIST neighbour list
!> \param potparm the pair-potential table, used to confirm the kind pair
!> \param cell the simulation cell, for the shift vectors
!> \param r_pbc wrapped positions
!> \param target_pot_type the potential type identifying this model
!> \param kind_included per-kind flag: is this kind handled by the model
!> \param cutsq per-kind-pair squared cutoff, already resolved by the caller
!> \param edges (2, nedges) global atom indices, allocated here
!> \param shifts (3, nedges) integer cell shifts, allocated here
!> \param nedges number of edges found
! **************************************************************************************************
   SUBROUTINE mlp_collect_edges(nonbonded, potparm, cell, r_pbc, target_pot_type, &
                                kind_included, cutsq, edges, shifts, nedges)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_pbc
      INTEGER, INTENT(IN)                                :: target_pot_type
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: kind_included
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN)         :: cutsq
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: edges(:, :), shifts(:, :)
      INTEGER, INTENT(OUT)                               :: nedges

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mlp_collect_edges'

      INTEGER                                            :: atom_a, atom_b, handle, i, iend, igrp, &
                                                            ikind, ilist, ipair, istart, jkind, &
                                                            n_max_edges
      INTEGER, DIMENSION(3)                              :: cvi
      LOGICAL                                            :: is_target
      REAL(kind=dp)                                      :: cutsq_ij, drij
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, rij
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      CALL timeset(routineN, handle)

      n_max_edges = 0
      DO ilist = 1, nonbonded%nlists
         n_max_edges = n_max_edges + nonbonded%neighbor_kind_pairs(ilist)%npairs
      END DO

      ALLOCATE (edges(2, n_max_edges), shifts(3, n_max_edges))
      nedges = 0

      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         IF (neighbor_kind_pair%npairs == 0) CYCLE

         Kind_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)

            IF (.NOT. kind_included(ikind)) CYCLE Kind_Loop
            IF (.NOT. kind_included(jkind)) CYCLE Kind_Loop

            pot => potparm%pot(ikind, jkind)%pot
            is_target = .FALSE.
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == target_pot_type) THEN
                  is_target = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. is_target) CYCLE Kind_Loop

            cutsq_ij = cutsq(ikind, jkind)
            cvi = neighbor_kind_pair%cell_vector
            cell_v = MATMUL(cell%hmat, REAL(cvi, dp))

            DO ipair = istart, iend
               atom_a = neighbor_kind_pair%list(1, ipair)
               atom_b = neighbor_kind_pair%list(2, ipair)

               rij(:) = r_pbc(atom_b)%r(:) - r_pbc(atom_a)%r(:) + cell_v
               drij = DOT_PRODUCT(rij, rij)

               IF (drij <= cutsq_ij) THEN
                  nedges = nedges + 1
                  edges(:, nedges) = [atom_a, atom_b]
                  shifts(:, nedges) = cvi
               END IF
            END DO
         END DO Kind_Loop
      END DO

      IF (nedges < n_max_edges) THEN
         BLOCK
            INTEGER, ALLOCATABLE :: tmp_e(:, :), tmp_s(:, :)
            ALLOCATE (tmp_e(2, nedges), tmp_s(3, nedges))
            tmp_e(:, :) = edges(:, 1:nedges)
            tmp_s(:, :) = shifts(:, 1:nedges)
            CALL MOVE_ALLOC(tmp_e, edges)
            CALL MOVE_ALLOC(tmp_s, shifts)
         END BLOCK
      END IF

      CALL timestop(handle)
   END SUBROUTINE mlp_collect_edges

! **************************************************************************************************
!> \brief Replicate every rank's edges so all ranks hold the complete list.
!>
!>        Each rank finds only the edges of the atoms in its own neighbour-list slice,
!>        but a message-passing model needs the whole graph, so the per-rank lists are
!>        concatenated. The serial case is a plain copy.
!> \param para_env the communicator
!> \param local_edges (2, n) this rank's edges
!> \param local_shifts (3, n) this rank's shifts
!> \param all_edges (2, ntot) the concatenated edges, allocated here
!> \param all_shifts (3, ntot) the concatenated shifts, allocated here
! **************************************************************************************************
   SUBROUTINE mlp_allgather_edges(para_env, local_edges, local_shifts, all_edges, all_shifts)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, INTENT(IN)                   :: local_edges(:, :), local_shifts(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: all_edges(:, :), all_shifts(:, :)

      CHARACTER(len=*), PARAMETER :: routineN = 'mlp_allgather_edges'

      INTEGER                                            :: handle, ip, nedges, nedges_tot
      INTEGER, ALLOCATABLE                               :: cnt_e(:), cnt_s(:), disp_e(:), &
                                                            disp_s(:), edge_count(:)

      CALL timeset(routineN, handle)

      nedges = SIZE(local_edges, 2)

      IF (para_env%num_pe == 1) THEN
         ALLOCATE (all_edges, SOURCE=local_edges)
         ALLOCATE (all_shifts, SOURCE=local_shifts)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (edge_count(para_env%num_pe))
      CALL para_env%allgather(nedges, edge_count)
      nedges_tot = SUM(edge_count)

      ALLOCATE (cnt_e(para_env%num_pe), cnt_s(para_env%num_pe))
      ALLOCATE (disp_e(para_env%num_pe), disp_s(para_env%num_pe))
      cnt_e(:) = edge_count*2
      cnt_s(:) = edge_count*3
      disp_e(1) = 0
      disp_s(1) = 0
      DO ip = 2, para_env%num_pe
         disp_e(ip) = disp_e(ip - 1) + cnt_e(ip - 1)
         disp_s(ip) = disp_s(ip - 1) + cnt_s(ip - 1)
      END DO

      ALLOCATE (all_edges(2, nedges_tot), all_shifts(3, nedges_tot))
      CALL para_env%allgatherv(local_edges, all_edges, cnt_e, disp_e)
      CALL para_env%allgatherv(local_shifts, all_shifts, cnt_s, disp_s)

      DEALLOCATE (edge_count, cnt_e, cnt_s, disp_e, disp_s)

      CALL timestop(handle)
   END SUBROUTINE mlp_allgather_edges

END MODULE manybody_mlp_neighbor

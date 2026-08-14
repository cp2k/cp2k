!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Spatial atom decomposition for the message-passing MLP drivers.
!>
!>        Used by the symmetrix MACE backend (manybody_symmetrix) so the ownership
!>        rule lives in one place and can be shared by further MLP drivers.
!>
!>        WHY THIS EXISTS. CP2K's own particle distribution (local_particles) is
!>        built by bin-packing MOLECULES onto ranks by workload - it is index and
!>        load based, never spatial. Using it to decide which atoms a rank's MLP
!>        evaluates is correct but slow: a rank's atoms end up scattered through
!>        the box, so the set of atoms within a cutoff of them (the ghosts, which
!>        a message-passing model needs to evaluate its owned atoms exactly) is
!>        essentially the whole system. Every rank then still does the first
!>        layer for every atom, which caps the speedup at a system-size
!>        independent ~2x and leaves per-rank memory untouched.
!>
!>        Assigning each rank a contiguous SLAB of the cell instead makes the
!>        ghost shell thin: it is the ~2*cutoff-deep skin of the slab, a fraction
!>        that shrinks as the cell grows. That is the standard MD decomposition
!>        (LAMMPS et al.) and it is what makes the work per rank actually fall.
!>
!>        WHAT THIS DELIBERATELY DOES NOT DO. It does not touch, replace or
!>        reorder any CP2K data structure: particle_set, r_last_update_pbc,
!>        local_particles, the FIST neighbour lists and the fist_force reductions
!>        are all left exactly as they are, and every other potential keeps using
!>        them unchanged. FIST is a replicated-data code by design and stays one.
!>        The partition computed here is private to one evaluation of one MLP: it
!>        is derived from the coordinates FIST already replicates on every rank,
!>        used to pick that rank's nodes, and discarded. Nothing is communicated -
!>        every rank runs identical arithmetic over identical replicated input in
!>        the same binary, so all ranks agree on the partition without exchanging
!>        a byte. Forces are still scattered to GLOBAL atom indices and summed by
!>        the existing para_env%sum in fist_force.F, which is the same contract
!>        the Allegro path already relies on.
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE manybody_mlp_partition

   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE fist_nonbond_env_types,          ONLY: pos_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mlp_slab_owner

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_mlp_partition'

CONTAINS

! **************************************************************************************************
!> \brief Assigns every candidate atom to a rank by a deterministic spatial slab.
!>
!>        Slabs are cut along the cell axis with the largest perpendicular width
!>        (the one that gives the fattest slabs, hence the thinnest ghost skin
!>        for a given rank count), and hold an EQUAL NUMBER OF ATOMS rather than
!>        an equal width. Equal counts balance the owned work exactly and stay
!>        sane for inhomogeneous systems - a slab of vacuum would otherwise get a
!>        rank to itself.
!>
!>        Determinism is the whole basis of the no-communication claim: every
!>        rank sorts the same replicated coordinates with the same tie-break and
!>        therefore lands on a bit-identical partition. The global atom index is
!>        the tie-break so that atoms sharing a coordinate cannot order
!>        differently between ranks.
!> \param candidates global indices of the atoms the model covers
!> \param r_pbc replicated PBC-wrapped positions (all atoms)
!> \param cell the simulation cell
!> \param para_env supplies num_pe and mepos
!> \param owned on output, .TRUE. for the candidates this rank owns (indexed like candidates)
! **************************************************************************************************
   SUBROUTINE mlp_slab_owner(candidates, r_pbc, cell, para_env, owned)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: candidates
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_pbc
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: owned

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mlp_slab_owner'

      INTEGER                                            :: axis, handle, i, ic, n_cand, rank, &
                                                            slab_hi, slab_lo
      INTEGER, ALLOCATABLE                               :: order(:)
      REAL(KIND=dp)                                      :: best_width, s(3), width
      REAL(KIND=dp), ALLOCATABLE                         :: key(:)

      CALL timeset(routineN, handle)

      n_cand = SIZE(candidates)
      CPASSERT(SIZE(owned) == n_cand)
      owned = .FALSE.

      ! One rank owns everything: the partition is the identity. Taking this
      ! shortcut keeps the serial path bit-identical to the undecomposed one.
      IF (para_env%num_pe == 1 .OR. n_cand == 0) THEN
         owned = .TRUE.
         CALL timestop(handle)
         RETURN
      END IF

      ! Cut along the axis with the largest perpendicular width. In scaled
      ! coordinates a slab of thickness ds along axis a is ds/|b_a| deep in real
      ! space, where b_a is the corresponding reciprocal (row of h_inv), so the
      ! widest axis is the one with the smallest reciprocal norm. Ties break on
      ! the lower axis index, which keeps the choice deterministic.
      axis = 1
      best_width = -1.0_dp
      DO i = 1, 3
         width = 1.0_dp/MAX(NORM2(cell%h_inv(i, :)), EPSILON(1.0_dp))
         IF (width > best_width + 1.0e-12_dp) THEN
            best_width = width
            axis = i
         END IF
      END DO

      ! Sort the candidates on their scaled coordinate along that axis. The
      ! global index enters as a scaled tie-break so equal coordinates cannot
      ! order differently on different ranks.
      ALLOCATE (key(n_cand), order(n_cand))
      DO ic = 1, n_cand
         CALL real_to_scaled(s, r_pbc(candidates(ic))%r(:), cell)
         key(ic) = s(axis) + REAL(candidates(ic), dp)*1.0e-12_dp
      END DO
      CALL sort(key, n_cand, order)

      ! Hand each rank a contiguous, equal-count block of the sorted order; the
      ! first MOD(n_cand, num_pe) ranks take one extra atom.
      CALL mlp_block_bounds(n_cand, para_env%num_pe, para_env%mepos, slab_lo, slab_hi)
      DO i = slab_lo, slab_hi
         owned(order(i)) = .TRUE.
      END DO

      ! Every candidate must be owned exactly once. Cheap enough to always check:
      ! a partition bug here is silent (plausible forces, wrong total) and this
      ! is the only place that can catch it.
      rank = COUNT(owned)
      CALL para_env%sum(rank)
      IF (rank /= n_cand) THEN
         CALL cp_abort(__LOCATION__, &
                       "MLP spatial partition does not cover the atoms exactly once. "// &
                       "The ranks disagree on the partition, which must never happen: it is "// &
                       "computed from replicated coordinates by identical arithmetic.")
      END IF

      DEALLOCATE (key, order)
      CALL timestop(handle)
   END SUBROUTINE mlp_slab_owner

! **************************************************************************************************
!> \brief Contiguous equal-count block of [1, n] for one rank (1-based, inclusive).
!> \param n number of items
!> \param num_pe number of ranks
!> \param mepos this rank (0-based)
!> \param lo first index of this rank's block
!> \param hi last index of this rank's block (hi < lo when the block is empty)
! **************************************************************************************************
   SUBROUTINE mlp_block_bounds(n, num_pe, mepos, lo, hi)
      INTEGER, INTENT(IN)                                :: n, num_pe, mepos
      INTEGER, INTENT(OUT)                               :: lo, hi

      INTEGER                                            :: base, extra

      base = n/num_pe
      extra = MOD(n, num_pe)
      lo = mepos*base + MIN(mepos, extra) + 1
      hi = lo + base - 1
      IF (mepos < extra) hi = hi + 1
   END SUBROUTINE mlp_block_bounds

END MODULE manybody_mlp_partition

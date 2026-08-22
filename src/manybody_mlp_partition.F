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
!>        Assigning each rank a compact REGION of the cell instead makes the
!>        ghost shell thin: it is the ~2*cutoff-deep skin of the region, a
!>        fraction that shrinks as the cell grows. That is the standard MD
!>        decomposition (LAMMPS et al.) and it is what makes the work per rank
!>        actually fall.
!>
!>        WHY THE REGIONS ARE 3D. Cutting along one axis only is simpler, but its
!>        speedup has a ceiling of L/(2*cutoff) that no amount of hardware moves:
!>        once a slab is thinner than its own ghost skin, every rank evaluates
!>        nearly the whole box. Bisecting recursively keeps the regions
!>        box-shaped, so their surface grows as the 2/3 power of their volume
!>        rather than staying constant, and the ceiling goes away.
!>
!>        WHY 3D IS NOT ALWAYS RIGHT. A slab that spans the whole cell in two
!>        directions pays no ghost in those directions, because there the shell
!>        wraps onto the slab itself. A compact box pays on all six faces. So for
!>        few enough domains the slab is genuinely the better shape, and only
!>        once a domain is comfortably wider than 2*cutoff does the box win.
!>
!>        Rather than pick a shape up front, each bisection cuts whichever axis
!>        leaves the smallest ghost volume, capped per axis by how far a shell can
!>        actually reach. That keeps single-axis cuts while domains are thin and
!>        turns into boxes as they fatten, with no rank-count threshold.
!>
!>        WHAT IT BUYS. The shape does what it is meant to. Mean layer-0 nodes
!>        per domain on 1024 waters, every domain captured and the owned counts
!>        summing to the system exactly:
!>
!>            domains      4      8     16     32
!>            compact   48.0%  34.4%  22.9%  14.8%
!>            one axis  54.7%  41.3%  32.7%  28.9%
!>
!>        So by 32 domains a compact region carries about half the layer-0 work
!>        of a slab, and the gap widens with domain count, which is the whole
!>        point of not cutting a single axis.
!>
!>        WHETHER THAT MAKES IT FASTER IS NOT YET SETTLED, and the honest state
!>        of the evidence is that two measurements disagree. A back-to-back
!>        comparison on a 10-core machine, both variants interleaved and best of
!>        three, gives 0.040 s against 0.099 s at 32 domains, which is rather
!>        better than the layer-0 ratio alone predicts. A run on 64 real cores
!>        found no difference at all, but that machine carried a load average
!>        between 13 and 21 while it ran and its ranks were not pinned, so it
!>        cannot settle anything. The small-machine run was oversubscribed,
!>        which is its own distortion.
!>
!>        Do not quote a speedup for the decomposition until a pinned run on a
!>        quiet machine resolves this. The layer-0 counts above are solid; the
!>        time is not.
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
   PUBLIC :: mlp_domain_owner

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_mlp_partition'

CONTAINS

! **************************************************************************************************
!> \brief Assigns every candidate atom to a rank and a thread by recursive
!>        coordinate bisection.
!>
!>        A domain range is halved by splitting the atoms along the axis whose
!>        bisection leaves the smallest ghost shell, and recursing until one
!>        domain is left. Each cut is by ATOM COUNT rather than by position, so
!>        the owned work balances exactly and inhomogeneous systems stay sane: a
!>        region of vacuum would otherwise get a rank to itself.
!>
!>        Ranks and threads are decomposed in two separate passes rather than as
!>        one flat num_pe*n_threads range. That keeps a rank's threads inside the
!>        rank's own box (so both shells stay compact for any thread count, not
!>        only powers of two), and it is what makes the cost acceptable: the rank
!>        pass descends only into the branch holding this rank, sorting n, then
!>        n/2, then n/4, which is the same O(n log n) a single sort costs. The
!>        thread pass then only ever touches this rank's share.
!>
!>        Determinism is the whole basis of the no-communication claim: every
!>        rank runs identical arithmetic over identical replicated coordinates
!>        and therefore lands on a bit-identical partition. Pruning does not
!>        weaken that, because the cuts on the path down are computed from the
!>        same subsets on every rank. Ties are resolved by holding every subset
!>        in ascending global-index order and relying on the sort being stable,
!>        rather than by perturbing the sort key: an index-proportional epsilon
!>        would grow with system size and, at a few hundred thousand atoms, smear
!>        the key by an appreciable fraction of an Angstrom.
!>
!>        BALANCE. With the proportional cut used below, FLOOR(n/D) is invariant
!>        down the recursion and domain d ends up with one extra atom exactly
!>        when d < MOD(n, D). That is the same assignment mlp_block_bounds makes
!>        in closed form, which is what the assertions below check against.
!> \param candidates global indices of the atoms the model covers
!> \param r_pbc replicated PBC-wrapped positions (all atoms)
!> \param cell the simulation cell
!> \param para_env supplies num_pe and mepos
!> \param n_threads domains per rank (1 for the unthreaded path)
!> \param rcut model interaction cutoff, same length unit as r_pbc
!> \param owned on output, .TRUE. for the candidates this rank owns (indexed like candidates)
!> \param thread_of on output, 0-based domain within this rank, -1 where not owned
! **************************************************************************************************
   SUBROUTINE mlp_domain_owner(candidates, r_pbc, cell, para_env, n_threads, rcut, owned, thread_of)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: candidates
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_pbc
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: n_threads
      REAL(KIND=dp), INTENT(IN)                          :: rcut
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: owned
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: thread_of

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mlp_domain_owner'

      INTEGER                                            :: d, handle, i, ic, it, n_cand, n_mine, &
                                                            span_hi, span_lo
      INTEGER, ALLOCATABLE                               :: dom_count(:), dom_of(:), mine(:), &
                                                            n_owners(:)
      REAL(KIND=dp)                                      :: cap(3), s(3), width(3)
      REAL(KIND=dp), ALLOCATABLE                         :: scaled(:, :)

      CALL timeset(routineN, handle)

      n_cand = SIZE(candidates)
      CPASSERT(SIZE(owned) == n_cand)
      CPASSERT(SIZE(thread_of) == n_cand)
      CPASSERT(n_threads >= 1)
      owned = .FALSE.
      thread_of = -1

      ! One rank on one thread owns everything: the partition is the identity.
      ! Taking this shortcut keeps the serial path bit-identical to the
      ! undecomposed one. It also covers the case of nothing to partition.
      IF ((para_env%num_pe == 1 .AND. n_threads == 1) .OR. n_cand == 0) THEN
         owned = .TRUE.
         thread_of = 0
         CALL timestop(handle)
         RETURN
      END IF

      ! Scaled coordinates once, and the real-space depth of a unit of scaled
      ! coordinate along each axis. A displacement ds along axis d spans
      ! ds/|b_d| in real space, b_d being the corresponding reciprocal vector
      ! (row of h_inv), so widths convert subset extents into comparable
      ! real-space lengths on a triclinic cell just as well as on a cubic one.
      ALLOCATE (scaled(3, n_cand))
      DO ic = 1, n_cand
         CALL real_to_scaled(s, r_pbc(candidates(ic))%r(:), cell)
         scaled(:, ic) = s(:)
      END DO
      DO d = 1, 3
         width(d) = 1.0_dp/MAX(NORM2(cell%h_inv(d, :)), EPSILON(1.0_dp))
      END DO

      ! Ceiling on how far a ghost shell can reach along each axis. On a periodic
      ! axis the shell can wrap right round, so the whole cell width is in reach.
      ! On an aperiodic one it can only reach as far as the atoms go, which is
      ! what keeps a thin slab in a tall vacuum box from being scored as though
      ! the vacuum were full of neighbours.
      DO d = 1, 3
         IF (cell%perd(d) == 1) THEN
            cap(d) = width(d)
         ELSE
            cap(d) = (MAXVAL(scaled(d, :)) - MINVAL(scaled(d, :)))*width(d)
         END IF
      END DO

      ALLOCATE (mine(n_cand), dom_of(n_cand))
      DO ic = 1, n_cand
         mine(ic) = ic
      END DO

      ! Rank pass, pruned to the branch that contains this rank.
      dom_of = -1
      CALL rcb_split(mine, n_cand, 0, para_env%num_pe - 1, para_env%mepos, &
                     scaled, width, cap, rcut, dom_of)
      DO ic = 1, n_cand
         owned(ic) = (dom_of(ic) == para_env%mepos)
      END DO
      n_mine = COUNT(owned)

      ! This rank must hold exactly the block the closed form assigns it. Local
      ! and exact, so it pins the recursion against its specification without a
      ! collective.
      CALL mlp_block_bounds(n_cand, para_env%num_pe, para_env%mepos, span_lo, span_hi)
      CPASSERT(n_mine == MAX(0, span_hi - span_lo + 1))

      ! Thread pass over this rank's own atoms, not pruned: a rank needs every
      ! one of its domains, not just one of them.
      IF (n_threads > 1 .AND. n_mine > 0) THEN
         i = 0
         DO ic = 1, n_cand
            IF (owned(ic)) THEN
               i = i + 1
               mine(i) = ic
            END IF
         END DO
         dom_of = -1
         CALL rcb_split(mine, n_mine, 0, n_threads - 1, -1, scaled, width, cap, rcut, dom_of)
         DO ic = 1, n_cand
            IF (owned(ic)) thread_of(ic) = dom_of(ic)
         END DO
         ! Tallied in one sweep rather than one COUNT per domain: the latter is
         ! n_threads passes over the candidates on the critical path of every
         ! force evaluation, which at high thread counts is real work for an
         ! assertion.
         ALLOCATE (dom_count(0:n_threads - 1))
         dom_count = 0
         DO ic = 1, n_cand
            IF (thread_of(ic) >= 0) dom_count(thread_of(ic)) = dom_count(thread_of(ic)) + 1
         END DO
         DO it = 0, n_threads - 1
            CALL mlp_block_bounds(n_mine, n_threads, it, span_lo, span_hi)
            CPASSERT(dom_count(it) == MAX(0, span_hi - span_lo + 1))
         END DO
         DEALLOCATE (dom_count)
      ELSE
         DO ic = 1, n_cand
            IF (owned(ic)) thread_of(ic) = 0
         END DO
      END IF

      ! Every candidate must be owned exactly once. Checked PER ATOM, not by
      ! comparing totals: a partition that double-claims one atom and orphans
      ! another has the right total and the wrong energy, which is precisely the
      ! silent failure this guard exists to catch. It is also the one test that
      ! sees the ranks disagreeing, which the local assertions above cannot.
      ! One integer per covered atom is negligible next to the edge data the
      ! caller already exchanges.
      ALLOCATE (n_owners(n_cand))
      n_owners = 0
      DO i = 1, n_cand
         IF (owned(i)) n_owners(i) = 1
      END DO
      CALL para_env%sum(n_owners)
      IF (ANY(n_owners /= 1)) THEN
         CALL cp_abort(__LOCATION__, &
                       "MLP spatial partition does not cover the atoms exactly once. "// &
                       "The ranks disagree on the partition, which must never happen: it is "// &
                       "computed from replicated coordinates by identical arithmetic.")
      END IF

      DEALLOCATE (scaled, mine, dom_of, n_owners)
      CALL timestop(handle)
   END SUBROUTINE mlp_domain_owner

! **************************************************************************************************
!> \brief Splits a subset of atoms across a range of domain indices by recursive
!>        coordinate bisection, writing the domain index of each atom it settles.
!>
!>        The subset is held in ascending order of its own indices, and both
!>        halves are rebuilt by sweeping it in that order, so the ascending
!>        invariant survives every level. That is what lets a stable sort resolve
!>        equal coordinates identically on every rank with no tie-break term.
!> \param sub positions into the candidate list, ASCENDING; only sub(1:n) is read
!> \param n number of atoms in the subset
!> \param dom_lo first domain index of the range
!> \param dom_hi last domain index of the range
!> \param my_dom recurse only into the branch holding this domain; negative for all
!> \param scaled scaled coordinates of every candidate, (3, n_cand)
!> \param width real-space depth of a unit of scaled coordinate, per axis
!> \param cap furthest a ghost shell can reach along each axis
!> \param rcut model interaction cutoff
!> \param dom_of on output, the domain index of every atom this call settles
! **************************************************************************************************
   RECURSIVE SUBROUTINE rcb_split(sub, n, dom_lo, dom_hi, my_dom, scaled, width, cap, rcut, dom_of)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sub
      INTEGER, INTENT(IN)                                :: n, dom_lo, dom_hi, my_dom
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scaled
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: width, cap
      REAL(KIND=dp), INTENT(IN)                          :: rcut
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: dom_of

      INTEGER                                            :: axis, base, d, d_lo, d_width, dom_mid, &
                                                            e, extra, i, n_hi, n_lo, p_hi, p_lo
      INTEGER, ALLOCATABLE                               :: order(:), sub_hi(:), sub_lo(:)
      LOGICAL, ALLOCATABLE                               :: take_lo(:)
      REAL(KIND=dp)                                      :: ext(3), ghost(3), hi_s, lo_s, side, &
                                                            vol_eps
      REAL(KIND=dp), ALLOCATABLE                         :: key(:)

      d_width = dom_hi - dom_lo + 1

      ! A single domain takes what is left. An empty subset settles nothing,
      ! which is how domains legitimately end up empty once there are more of
      ! them than there are atoms.
      IF (d_width == 1) THEN
         DO i = 1, n
            dom_of(sub(i)) = dom_lo
         END DO
         RETURN
      END IF
      IF (n == 0) RETURN

      ! Real-space extent of this subset along each axis.
      DO d = 1, 3
         lo_s = scaled(d, sub(1))
         hi_s = lo_s
         DO i = 2, n
            lo_s = MIN(lo_s, scaled(d, sub(i)))
            hi_s = MAX(hi_s, scaled(d, sub(i)))
         END DO
         ext(d) = (hi_s - lo_s)*width(d)
      END DO

      ! Score each axis by the ghost volume that cutting it would leave: halve
      ! this subset along that axis, grow every axis by the cutoff on both sides,
      ! and clamp to how far a shell can actually reach. Halving an axis that is
      ! already no wider than its own shell buys nothing, and the clamp is what
      ! says so, which is why this reproduces slabs while domains are thin and
      ! turns into boxes once they are fat. Recomputing per subset is also what
      ! lets an interface or a droplet be cut on its own shape rather than the
      ! cell's.
      DO d = 1, 3
         ghost(d) = 1.0_dp
         DO e = 1, 3
            side = ext(e)
            IF (e == d) side = 0.5_dp*side
            ghost(d) = ghost(d)*MIN(cap(e), side + 2.0_dp*rcut)
         END DO
      END DO

      ! Smallest ghost volume wins; on a tie the widest axis does, which is what
      ! keeps the choice sensible once every axis is clamped and the volumes stop
      ! discriminating. Both comparisons are one-sided, so the lowest axis index
      ! survives a genuine tie and the choice stays deterministic across ranks.
      vol_eps = 1.0e-12_dp*MAXVAL(ghost)
      axis = 1
      DO d = 2, 3
         IF (ghost(d) < ghost(axis) - vol_eps) THEN
            axis = d
         ELSE IF (ghost(d) <= ghost(axis) + vol_eps .AND. ext(d) > ext(axis) + 1.0e-12_dp) THEN
            axis = d
         END IF
      END DO

      ! Give the lower half of the domain range its exact share of the atoms.
      ! Splitting the count this way, rather than cutting at a coordinate, is
      ! what keeps every leaf at FLOOR(n/D) or one more.
      d_lo = (d_width + 1)/2
      base = n/d_width
      extra = MOD(n, d_width)
      n_lo = d_lo*base + MIN(d_lo, extra)
      n_hi = n - n_lo

      ALLOCATE (key(n), order(n))
      DO i = 1, n
         key(i) = scaled(axis, sub(i))
      END DO
      CALL sort(key, n, order)

      ALLOCATE (take_lo(n))
      take_lo = .FALSE.
      DO i = 1, n_lo
         take_lo(order(i)) = .TRUE.
      END DO

      ALLOCATE (sub_lo(n_lo), sub_hi(n_hi))
      p_lo = 0
      p_hi = 0
      DO i = 1, n
         IF (take_lo(i)) THEN
            p_lo = p_lo + 1
            sub_lo(p_lo) = sub(i)
         ELSE
            p_hi = p_hi + 1
            sub_hi(p_hi) = sub(i)
         END IF
      END DO
      DEALLOCATE (key, order, take_lo)

      dom_mid = dom_lo + d_lo - 1
      IF (my_dom < 0 .OR. my_dom <= dom_mid) THEN
         CALL rcb_split(sub_lo, n_lo, dom_lo, dom_mid, my_dom, scaled, width, cap, rcut, dom_of)
      END IF
      IF (my_dom < 0 .OR. my_dom > dom_mid) THEN
         CALL rcb_split(sub_hi, n_hi, dom_mid + 1, dom_hi, my_dom, scaled, width, cap, rcut, dom_of)
      END IF
      DEALLOCATE (sub_lo, sub_hi)
   END SUBROUTINE rcb_split

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

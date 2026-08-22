!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Linked-cell neighbour search for the message-passing MLP drivers.
!>
!>        WHY THIS EXISTS. A message-passing model evaluates its owned atoms from
!>        a two-hop environment, so a decomposed driver needs the edges centred on
!>        its own atoms and on their one-hop ghosts, and no others. Taking those
!>        out of CP2K's FIST neighbour list means every rank first has to learn
!>        the whole system's edge list, because the FIST list is distributed by
!>        CP2K's own particle map rather than spatially. That is what the
!>        symmetrix driver used to do, and the all-gather behind it grows without
!>        bound: half a gigabyte per rank at half a million atoms.
!>
!>        Building the edges here instead costs nothing extra, because FIST
!>        replicates every coordinate on every rank already. A rank bins all the
!>        atoms once and then walks out only from the ones it needs, so both the
!>        work and the memory follow the domain rather than the system.
!>
!>        IMAGES ARE NOT REPLICATED. The obvious way to handle periodicity is to
!>        pad the atom set with copies of itself, one ring per periodic axis. At
!>        half a million atoms that is thirteen million image records and about
!>        the same half gigabyte the all-gather cost, so it would trade one wall
!>        for another. Instead the stencil wraps: stepping off the end of the bin
!>        grid returns the bin on the far side together with the number of cells
!>        crossed, and that integer count becomes the edge's cell shift. Memory
!>        stays proportional to the atom count, boxes shorter than twice the
!>        cutoff work by letting the stencil reach more than one bin, and an
!>        aperiodic axis is handled by refusing any step that would wrap it.
!>
!>        ORDER IS REPRODUCIBLE WITHOUT SORTING. Every rank bins the same
!>        replicated coordinates into the same grid by the same arithmetic, and
!>        the chains are built in one deterministic sweep, so an atom's
!>        neighbours come out in the same order on every rank and every thread.
!>        A domain that holds an atom as a ghost therefore builds exactly the
!>        edge list, in exactly the order, that the domain owning it builds. That
!>        is what makes a decomposed result independent of how it was decomposed,
!>        and it is why no canonical re-sort of the edges is needed.
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE manybody_mlp_celllist

   USE ISO_C_BINDING,                   ONLY: C_INT32_T
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_1,&
                                              int_8
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mlp_cell_list_type, mlp_edge_buffer_type
   PUBLIC :: mlp_cell_list_build, mlp_cell_list_release
   PUBLIC :: mlp_edge_buffer_reserve, mlp_edge_buffer_release
   PUBLIC :: mlp_expand_tier, mlp_tier_collect
   PUBLIC :: mlp_tier_none, mlp_tier_owned, mlp_tier_ghost, mlp_tier_outer

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_mlp_celllist'

   ! Node tiers. Owned atoms get the readouts, ghosts are the one-hop shell that
   ! makes the owned atoms' first layer exact, and outer atoms are the two-hop
   ! shell that makes the GHOSTS' first layer exact. Outer atoms are edge targets
   ! only: they carry a position and a type and they receive forces, but no
   ! feature is ever computed for them.
   INTEGER(KIND=int_1), PARAMETER :: mlp_tier_none = -1_int_1, &
                                     mlp_tier_owned = 0_int_1, &
                                     mlp_tier_ghost = 1_int_1, &
                                     mlp_tier_outer = 2_int_1

   ! A wrapped stencil can only ever produce a handful of distinct cell shifts.
   ! Refusing an absurd number of them turns a nonsensical cell into a clear
   ! abort rather than a silent allocation of gigabytes.
   INTEGER, PARAMETER, PRIVATE :: max_shift_codes = 50000

   ! Ceiling on the stencil reach. Checked against the unrounded ratio, because
   ! rounding it first is what turns a degenerate axis into a negative reach and
   ! an empty stencil. Checked PER AXIS rather than through the shift-code count,
   ! which only notices when two or more axes are degenerate at once.
   !
   ! On a periodic axis this is currently belt and braces: CP2K's own neighbour
   ! list sizes itself as (2*cutoff/plane_distance + 1)**3 lists and exhausts
   ! memory on such a cell long before the search here is reached. It is kept
   ! because that list is now built only for the cross-check, and because the
   ! aperiodic branch below relies on the same reasoning to justify clamping
   ! rather than aborting.
   INTEGER, PARAMETER, PRIVATE :: max_span = 1000

   ! Ceilings on the grid itself. Both are generous next to any real cell: at a
   ! 4 Angstrom cutoff the per-axis limit is a 20 micron box, and the total is
   ! reached long before the memory is. They exist because the products overflow
   ! a default integer well before they become absurd, and because ALLOCATE with
   ! a negative extent yields a zero-size array instead of failing, so the walk
   ! would then write outside it.
   INTEGER, PARAMETER, PRIVATE :: max_bins_per_axis = 50000
   INTEGER, PARAMETER, PRIVATE :: max_cells = 200000000

   TYPE mlp_cell_list_type
      LOGICAL                        :: initialized = .FALSE.
      INTEGER                        :: n_atoms = 0
      INTEGER                        :: n_cells = 1
      INTEGER                        :: n_shift_codes = 1
      INTEGER, DIMENSION(3)          :: nbin = 1
      INTEGER, DIMENSION(3)          :: span = 0        ! stencil reach, in bins
      INTEGER, DIMENSION(3)          :: perd = 0
      INTEGER, DIMENSION(3)          :: shift_lo = 0
      INTEGER, DIMENSION(3)          :: shift_n = 1     ! shifts per axis, for the code packing
      REAL(KIND=dp)                  :: rcutsq = -1.0_dp
      REAL(KIND=dp), DIMENSION(3)    :: s_lower = -0.5_dp
      REAL(KIND=dp), DIMENSION(3)    :: s_span = 1.0_dp
      ! Positions verbatim as handed in, never re-wrapped: they have to stay the
      ! ones the rest of the driver uses, or an edge's shift stops matching its
      ! displacement.
      REAL(KIND=dp), ALLOCATABLE     :: pos(:, :)         ! (3, n_atoms)
      ! Cartesian translation per packed shift code, so the walk can add a cell
      ! vector with one lookup and no matrix product.
      REAL(KIND=dp), ALLOCATABLE     :: shift_trans(:, :) ! (3, n_shift_codes)
      INTEGER, ALLOCATABLE           :: home_bin(:, :)    ! (3, n_atoms), 1-based
      INTEGER, ALLOCATABLE           :: home_off(:, :)    ! (3, n_atoms), cells removed to reach the grid
      INTEGER, ALLOCATABLE           :: head(:)           ! (n_cells)
      INTEGER, ALLOCATABLE           :: next(:)           ! (n_atoms)
   END TYPE mlp_cell_list_type

   ! Edges for one domain, grouped by centre in ascending centre order. Indices
   ! are the driver's canonical packed numbering, zero based, so they can go
   ! straight into the evaluator's arrays after the local relabelling.
   TYPE mlp_edge_buffer_type
      INTEGER                         :: n_edges = 0
      INTEGER                         :: capacity = 0
      INTEGER(C_INT32_T), ALLOCATABLE :: pairs(:, :)   ! (2, capacity): centre, neighbour
      INTEGER(C_INT32_T), ALLOCATABLE :: shifts(:, :)  ! (3, capacity)
   END TYPE mlp_edge_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Bins the atoms and works out the stencil the cutoff needs.
!> \param cl the cell list to fill
!> \param pos (3, n_atoms) positions, taken verbatim from the driver
!> \param cell the simulation cell, for the metric and the periodicity
!> \param rcutsq squared interaction cutoff, in the units of pos
! **************************************************************************************************
   SUBROUTINE mlp_cell_list_build(cl, pos, cell, rcutsq)
      TYPE(mlp_cell_list_type), INTENT(INOUT)            :: cl
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pos
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: rcutsq

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mlp_cell_list_build'

      CHARACTER(LEN=default_string_length*8)             :: msg
      INTEGER                                            :: b, d, handle, i, ia, ib, ic, n_atoms, &
                                                            off_hi(3), off_lo(3), w_hi(3), w_lo(3)
      INTEGER(KIND=int_8)                                :: n_cells_8, n_codes_8
      INTEGER, DIMENSION(3)                              :: sh
      REAL(KIND=dp)                                      :: depth, rcut, s(3), shifted, width(3)
      REAL(KIND=dp), ALLOCATABLE                         :: scal(:, :)

      CALL timeset(routineN, handle)

      n_atoms = SIZE(pos, 2)
      rcut = SQRT(rcutsq)
      cl%n_atoms = n_atoms
      cl%rcutsq = rcutsq
      cl%perd(:) = cell%perd(:)

      IF (ALLOCATED(cl%pos)) DEALLOCATE (cl%pos)
      ALLOCATE (cl%pos(3, n_atoms))
      cl%pos(:, :) = pos(:, :)

      IF (n_atoms == 0) THEN
         cl%initialized = .TRUE.
         CALL timestop(handle)
         RETURN
      END IF

      ! Real-space depth of one unit of scaled coordinate along each axis. Row d
      ! of h_inv is the reciprocal vector normal to the d faces, so its inverse
      ! norm is the interplanar spacing: identically what plane_distance returns
      ! for (1,0,0), (0,1,0) and (0,0,1), and what the FIST neighbour list uses
      ! for the same purpose. Written out here rather than called because this
      ! form stays finite on a degenerate cell, where plane_distance divides by
      ! zero, and because it needs no round trip through the cell angles.
      DO d = 1, 3
         width(d) = 1.0_dp/MAX(NORM2(cell%h_inv(d, :)), EPSILON(1.0_dp))
      END DO

      ALLOCATE (scal(3, n_atoms))
      DO i = 1, n_atoms
         CALL real_to_scaled(s, pos(:, i), cell)
         scal(:, i) = s(:)
      END DO

      ! A periodic axis is gridded over the whole cell, since the shell wraps
      ! round it. An aperiodic one is gridded over the atoms only, so a thin slab
      ! in a tall vacuum box gets bins where its atoms are instead of one bin
      ! holding everything.
      DO d = 1, 3
         IF (cl%perd(d) /= 0) THEN
            cl%s_lower(d) = -0.5_dp
            cl%s_span(d) = 1.0_dp
         ELSE
            cl%s_lower(d) = MINVAL(scal(d, :))
            cl%s_span(d) = MAX(MAXVAL(scal(d, :)) - cl%s_lower(d), 1.0E-12_dp)
            ! a hair of margin so the topmost atom cannot land one bin past the end
            cl%s_lower(d) = cl%s_lower(d) - 1.0E-12_dp
            cl%s_span(d) = cl%s_span(d) + 2.0E-12_dp
         END IF
      END DO

      ! Bins at least one cutoff deep, so the stencil is the immediate ring and
      ! no more. A finer grid would only add empty bins to walk. When the cell is
      ! shorter than the cutoff the count clamps at one and the stencil reaches
      ! further instead, which is what makes small boxes work.
      DO d = 1, 3
         ! Bound the RATIO, not the rounded count. FLOOR of something past the
         ! reach of a default integer is undefined, and in practice comes back
         ! as a small or negative number, so a guard placed after it would let
         ! through exactly the degenerate input it exists to refuse.
         !
         ! On a periodic axis the grid spans the cell, so this ratio is the cell
         ! width over the cutoff and is bounded by the cell. On an aperiodic one
         ! it spans the ATOMS, and PERIODIC NONE does not require the atoms to
         ! sit inside the box, so a wide enough cluster would otherwise grow the
         ! grid without limit.
         IF (cl%s_span(d)*width(d) > REAL(max_bins_per_axis, dp)*rcut) THEN
            WRITE (msg, '(A,I0,A)') &
               "The MLP neighbour search would need more than ", max_bins_per_axis, &
               " bins along one cell axis. Either the cell is enormous compared with the model "// &
               "cutoff, or, on a non-periodic axis, the atoms span far more than the cell does."
            CALL cp_abort(__LOCATION__, TRIM(msg))
         END IF
         cl%nbin(d) = MAX(1, FLOOR(cl%s_span(d)*width(d)/rcut))
         IF (cl%nbin(d) > max_bins_per_axis) THEN
            WRITE (msg, '(A,I0,A,I0,A)') &
               "The MLP neighbour search would need ", cl%nbin(d), " bins along one cell axis, "// &
               "more than the limit of ", max_bins_per_axis, ". Either the cell is enormous "// &
               "compared with the model cutoff, or, on a non-periodic axis, the atoms span far "// &
               "more than the cell does."
            CALL cp_abort(__LOCATION__, TRIM(msg))
         END IF
         depth = cl%s_span(d)*width(d)/REAL(cl%nbin(d), dp)
         IF (cl%perd(d) == 0) THEN
            ! An aperiodic axis is gridded over the atoms and never walked past
            ! its own ends, so reaching every bin always covers every atom on
            ! it. That holds however thin the axis is, which matters because a
            ! flat slab has an aperiodic extent of essentially nothing and the
            ! unclamped reach below would be enormous.
            cl%span(d) = MIN(CEILING(MIN(rcut/MAX(depth, EPSILON(1.0_dp)), &
                                         REAL(max_span, dp))), cl%nbin(d))
         ELSE
            ! A periodic axis has to be reached across in full, so here the
            ! ratio is real and is checked BEFORE rounding: an axis vanishingly
            ! thin against the cutoff gives a value that overflows a default
            ! integer and comes back negative, which empties the stencil rather
            ! than filling it and silently finds no neighbours at all.
            IF (rcut > REAL(max_span, dp)*depth) THEN
               CALL cp_abort(__LOCATION__, &
                             "The MLP neighbour search would have to reach across more periodic "// &
                             "images than its stencil allows along one cell axis. That axis is "// &
                             "vanishingly thin compared with the model cutoff; check the CELL "// &
                             "section against the model's interaction range.")
            END IF
            cl%span(d) = CEILING(rcut/MAX(depth, EPSILON(1.0_dp)))
         END IF
      END DO
      ! Accumulated wide: three counts each within the per-axis limit still
      ! multiply past a default integer, and a wrapped negative extent makes
      ! ALLOCATE produce a zero-size array rather than fail, after which the
      ! walk writes outside it.
      n_cells_8 = INT(cl%nbin(1), int_8)*INT(cl%nbin(2), int_8)*INT(cl%nbin(3), int_8)
      IF (n_cells_8 > INT(max_cells, int_8)) THEN
         WRITE (msg, '(A,I0,A,I0,A)') &
            "The MLP neighbour search would need ", n_cells_8, " bins in total, more than the "// &
            "limit of ", max_cells, ". The cell is very large compared with the model cutoff."
         CALL cp_abort(__LOCATION__, TRIM(msg))
      END IF
      cl%n_cells = INT(n_cells_8)

      ! Home bin, and how many whole cells had to come off to get there.
      IF (ALLOCATED(cl%home_bin)) DEALLOCATE (cl%home_bin)
      IF (ALLOCATED(cl%home_off)) DEALLOCATE (cl%home_off)
      ALLOCATE (cl%home_bin(3, n_atoms), cl%home_off(3, n_atoms))
      DO i = 1, n_atoms
         DO d = 1, 3
            IF (cl%perd(d) == 1) THEN
               cl%home_off(d, i) = FLOOR(scal(d, i) + 0.5_dp)
            ELSE
               cl%home_off(d, i) = 0
            END IF
            shifted = scal(d, i) - REAL(cl%home_off(d, i), dp)
            b = FLOOR((shifted - cl%s_lower(d))/cl%s_span(d)*REAL(cl%nbin(d), dp)) + 1
            cl%home_bin(d, i) = MIN(MAX(b, 1), cl%nbin(d))
         END DO
      END DO
      DEALLOCATE (scal)

      ! Chains, swept from the last atom down so each chain comes out in
      ! ascending atom order. Deterministic, and identical on every rank.
      IF (ALLOCATED(cl%head)) DEALLOCATE (cl%head)
      IF (ALLOCATED(cl%next)) DEALLOCATE (cl%next)
      ALLOCATE (cl%head(cl%n_cells), cl%next(n_atoms))
      cl%head(:) = 0
      DO i = n_atoms, 1, -1
         ic = linear_bin(cl, cl%home_bin(:, i))
         cl%next(i) = cl%head(ic)
         cl%head(ic) = i
      END DO

      ! Every shift the walk can emit, tabulated once. The wrap contributes the
      ! number of cells the stencil crosses; the two home offsets contribute the
      ! rest.
      DO d = 1, 3
         w_lo(d) = wrap_count(1 - cl%span(d), cl%nbin(d))
         w_hi(d) = wrap_count(cl%nbin(d) + cl%span(d), cl%nbin(d))
         off_lo(d) = MINVAL(cl%home_off(d, :))
         off_hi(d) = MAXVAL(cl%home_off(d, :))
         cl%shift_lo(d) = w_lo(d) + off_lo(d) - off_hi(d)
         cl%shift_n(d) = (w_hi(d) + off_hi(d) - off_lo(d)) - cl%shift_lo(d) + 1
      END DO
      ! Accumulated wide for the same reason as the bin count: the product of
      ! three plausible per-axis counts can still wrap, and a wrapped value that
      ! happens to land in range would pass a narrow test.
      n_codes_8 = INT(cl%shift_n(1), int_8)*INT(cl%shift_n(2), int_8)*INT(cl%shift_n(3), int_8)
      IF (n_codes_8 > INT(max_shift_codes, int_8) .OR. n_codes_8 < 1_int_8) THEN
         CALL cp_abort(__LOCATION__, &
                       "The MLP neighbour search would need an implausible number of periodic "// &
                       "images. The cell is either far smaller than the model cutoff or badly "// &
                       "distorted; check the CELL section against the model's interaction range.")
      END IF
      cl%n_shift_codes = INT(n_codes_8)

      IF (ALLOCATED(cl%shift_trans)) DEALLOCATE (cl%shift_trans)
      ALLOCATE (cl%shift_trans(3, cl%n_shift_codes))
      DO ic = 1, cl%shift_n(3)
         DO ib = 1, cl%shift_n(2)
            DO ia = 1, cl%shift_n(1)
               sh(1) = cl%shift_lo(1) + ia - 1
               sh(2) = cl%shift_lo(2) + ib - 1
               sh(3) = cl%shift_lo(3) + ic - 1
               cl%shift_trans(:, shift_code(cl, sh)) = MATMUL(cell%hmat, REAL(sh, dp))
            END DO
         END DO
      END DO

      cl%initialized = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE mlp_cell_list_build

! **************************************************************************************************
!> \brief Frees a cell list.
!> \param cl ...
! **************************************************************************************************
   SUBROUTINE mlp_cell_list_release(cl)
      TYPE(mlp_cell_list_type), INTENT(INOUT)            :: cl

      IF (ALLOCATED(cl%pos)) DEALLOCATE (cl%pos)
      IF (ALLOCATED(cl%shift_trans)) DEALLOCATE (cl%shift_trans)
      IF (ALLOCATED(cl%home_bin)) DEALLOCATE (cl%home_bin)
      IF (ALLOCATED(cl%home_off)) DEALLOCATE (cl%home_off)
      IF (ALLOCATED(cl%head)) DEALLOCATE (cl%head)
      IF (ALLOCATED(cl%next)) DEALLOCATE (cl%next)
      cl%initialized = .FALSE.
      cl%n_atoms = 0
   END SUBROUTINE mlp_cell_list_release

! **************************************************************************************************
!> \brief Makes room for at least n edges, without ever shrinking.
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE mlp_edge_buffer_reserve(buf, n)
      TYPE(mlp_edge_buffer_type), INTENT(INOUT)          :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (buf%capacity >= n) RETURN
      IF (ALLOCATED(buf%pairs)) DEALLOCATE (buf%pairs)
      IF (ALLOCATED(buf%shifts)) DEALLOCATE (buf%shifts)
      ALLOCATE (buf%pairs(2, n), buf%shifts(3, n))
      buf%capacity = n
   END SUBROUTINE mlp_edge_buffer_reserve

! **************************************************************************************************
!> \brief Frees an edge buffer.
!> \param buf ...
! **************************************************************************************************
   SUBROUTINE mlp_edge_buffer_release(buf)
      TYPE(mlp_edge_buffer_type), INTENT(INOUT)          :: buf

      IF (ALLOCATED(buf%pairs)) DEALLOCATE (buf%pairs)
      IF (ALLOCATED(buf%shifts)) DEALLOCATE (buf%shifts)
      buf%capacity = 0
      buf%n_edges = 0
   END SUBROUTINE mlp_edge_buffer_release

! **************************************************************************************************
!> \brief Walks the stencil around every given centre, tagging the neighbours it
!>        reaches and, if asked, recording the edges.
!>
!>        Called twice per tier: once to count and tag, then once to fill a
!>        buffer sized exactly from that count. Tagging is idempotent, so the
!>        second pass reproduces the first edge for edge. Counting first is not
!>        wasted work, because the per-centre counts it returns are the CSR's
!>        neighbour counts, which would otherwise need their own sweep.
!> \param cl the cell list, read only and safe to share between threads
!> \param centres canonical packed 0-based centres, ASCENDING
!> \param tier per-atom tier, updated in place for anything still untagged
!> \param new_tier tier to give the neighbours this call reaches
!> \param buf edges land here while they fit; n_edges counts them all regardless
!> \param do_emit .FALSE. to count and tag only, .TRUE. to also write what fits
!> \param counts optional per-centre edge count, indexed like centres
! **************************************************************************************************
   SUBROUTINE mlp_expand_tier(cl, centres, tier, new_tier, buf, do_emit, counts)
      TYPE(mlp_cell_list_type), INTENT(IN)               :: cl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: centres
      INTEGER(KIND=int_1), DIMENSION(:), INTENT(INOUT)   :: tier
      INTEGER(KIND=int_1), INTENT(IN)                    :: new_tier
      TYPE(mlp_edge_buffer_type), INTENT(INOUT)          :: buf
      LOGICAL, INTENT(IN)                                :: do_emit
      INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL       :: counts

      INTEGER                                            :: bx, by, bz, code, i, ic, j, kx, ky, kz, &
                                                            n_here, wx, wy, wz
      INTEGER, DIMENSION(3)                              :: bi, oi, sh
      REAL(KIND=dp)                                      :: d2, ri(3), rij(3)

      DO ic = 1, SIZE(centres)
         i = centres(ic) + 1
         ri(:) = cl%pos(:, i)
         bi(:) = cl%home_bin(:, i)
         oi(:) = cl%home_off(:, i)
         n_here = 0

         DO kx = bi(1) - cl%span(1), bi(1) + cl%span(1)
            bx = MODULO(kx - 1, cl%nbin(1)) + 1
            wx = (kx - bx)/cl%nbin(1)
            IF (cl%perd(1) == 0 .AND. wx /= 0) CYCLE
            DO ky = bi(2) - cl%span(2), bi(2) + cl%span(2)
               by = MODULO(ky - 1, cl%nbin(2)) + 1
               wy = (ky - by)/cl%nbin(2)
               IF (cl%perd(2) == 0 .AND. wy /= 0) CYCLE
               DO kz = bi(3) - cl%span(3), bi(3) + cl%span(3)
                  bz = MODULO(kz - 1, cl%nbin(3)) + 1
                  wz = (kz - bz)/cl%nbin(3)
                  IF (cl%perd(3) == 0 .AND. wz /= 0) CYCLE

                  j = cl%head(bz + cl%nbin(3)*(by - 1 + cl%nbin(2)*(bx - 1)))
                  DO WHILE (j > 0)
                     sh(1) = wx + oi(1) - cl%home_off(1, j)
                     sh(2) = wy + oi(2) - cl%home_off(2, j)
                     sh(3) = wz + oi(3) - cl%home_off(3, j)
                     ! An atom is its own neighbour through any image but the
                     ! identity one, which matters as soon as the cell is
                     ! shorter than twice the cutoff.
                     IF (j == i .AND. sh(1) == 0 .AND. sh(2) == 0 .AND. sh(3) == 0) THEN
                        j = cl%next(j)
                        CYCLE
                     END IF
                     code = shift_code(cl, sh)
                     rij(:) = cl%pos(:, j) - ri(:) + cl%shift_trans(:, code)
                     d2 = DOT_PRODUCT(rij, rij)
                     IF (d2 <= cl%rcutsq) THEN
                        IF (tier(j) == mlp_tier_none) tier(j) = new_tier
                        buf%n_edges = buf%n_edges + 1
                        n_here = n_here + 1
                        IF (do_emit .AND. buf%n_edges <= buf%capacity) THEN
                           buf%pairs(1, buf%n_edges) = INT(i - 1, C_INT32_T)
                           buf%pairs(2, buf%n_edges) = INT(j - 1, C_INT32_T)
                           buf%shifts(1, buf%n_edges) = INT(sh(1), C_INT32_T)
                           buf%shifts(2, buf%n_edges) = INT(sh(2), C_INT32_T)
                           buf%shifts(3, buf%n_edges) = INT(sh(3), C_INT32_T)
                        END IF
                     END IF
                     j = cl%next(j)
                  END DO
               END DO
            END DO
         END DO

         IF (PRESENT(counts)) counts(ic) = n_here
      END DO
   END SUBROUTINE mlp_expand_tier

! **************************************************************************************************
!> \brief Lists, in ascending order, every atom currently carrying a given tier.
!> \param tier ...
!> \param want ...
!> \param list 0-based packed indices, must be large enough
!> \param n_list ...
! **************************************************************************************************
   SUBROUTINE mlp_tier_collect(tier, want, list, n_list)
      INTEGER(KIND=int_1), DIMENSION(:), INTENT(IN)      :: tier
      INTEGER(KIND=int_1), INTENT(IN)                    :: want
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: list
      INTEGER, INTENT(OUT)                               :: n_list

      INTEGER                                            :: i

      n_list = 0
      DO i = 1, SIZE(tier)
         IF (tier(i) == want) THEN
            n_list = n_list + 1
            list(n_list) = i - 1
         END IF
      END DO
   END SUBROUTINE mlp_tier_collect

! **************************************************************************************************
!> \brief Number of whole grids crossed to bring a stencil index back in range.
!> \param k stencil index, may sit outside [1, nbin]
!> \param nbin ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION wrap_count(k, nbin) RESULT(w)
      INTEGER, INTENT(IN)                                :: k, nbin
      INTEGER                                            :: w

      w = (k - (MODULO(k - 1, nbin) + 1))/nbin
   END FUNCTION wrap_count

! **************************************************************************************************
!> \brief Packs a cell shift into its index in the translation table.
!> \param cl ...
!> \param sh ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION shift_code(cl, sh) RESULT(code)
      TYPE(mlp_cell_list_type), INTENT(IN)               :: cl
      INTEGER, DIMENSION(3), INTENT(IN)                  :: sh
      INTEGER                                            :: code

      code = 1 + (sh(1) - cl%shift_lo(1)) &
             + cl%shift_n(1)*((sh(2) - cl%shift_lo(2)) &
                              + cl%shift_n(2)*(sh(3) - cl%shift_lo(3)))
   END FUNCTION shift_code

! **************************************************************************************************
!> \brief Flattens a 1-based bin triple.
!> \param cl ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION linear_bin(cl, b) RESULT(ic)
      TYPE(mlp_cell_list_type), INTENT(IN)               :: cl
      INTEGER, DIMENSION(3), INTENT(IN)                  :: b
      INTEGER                                            :: ic

      ic = b(3) + cl%nbin(3)*(b(2) - 1 + cl%nbin(2)*(b(1) - 1))
   END FUNCTION linear_bin

END MODULE manybody_mlp_celllist

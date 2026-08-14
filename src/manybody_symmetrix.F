!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Driver for symmetrix MACE potentials in the FIST path.
!>
!>        The FIST neighbour list is converted into symmetrix's CSR contract (a
!>        per-atom neighbour-count array plus flattened neighbour indices/types,
!>        displacement vectors and distances grouped by central atom) and handed
!>        to symmetrix_eval, which returns the total energy and per-atom forces.
!>        symmetrix works in the model's native units (Angstrom / eV); positions
!>        are converted from Bohr and the returned energy/forces back to Hartree
!>        and Hartree/Bohr here.
!>
!>        The model is ATOM DECOMPOSED across MPI ranks: each rank owns a
!>        contiguous SLAB of the cell (manybody_mlp_partition) and evaluates the
!>        model only for those atoms, with the atoms around them present as
!>        ghosts so the owned atoms' 2-layer message-passing environments stay
!>        exact. Layer 0 runs over owned+ghosts, layer 1 and the readouts over
!>        the owned nodes only. Because FIST replicates every atom's position on
!>        every rank the ghosts cost no communication at all, and CP2K's own data
!>        structures - including its particle distribution - are left untouched.
!>        Each rank's energy/forces/virial are genuine partial contributions and
!>        are NOT scaled by 1/num_pe; the existing cross-rank sums in
!>        fist_force.F assemble the exact totals. This differs from the TorchScript
!>        MACE and NequIP paths in manybody_e3nn.F, which evaluate the whole system
!>        on every rank and divide by num_pe; it matches the masking the Allegro
!>        branch there applies to the energy.
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE manybody_symmetrix

   USE ISO_C_BINDING,                   ONLY: C_INT32_T
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Batatia2022,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type,&
                                              pos_type,&
                                              symmetrix_data_type
   USE kinds,                           ONLY: dp
   USE manybody_mlp_neighbor,           ONLY: mlp_allgather_edges,&
                                              mlp_collect_edges
   USE manybody_mlp_partition,          ONLY: mlp_slab_owner
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                              pair_potential_single_type,&
                                              symmetrix_pot_type,&
                                              symmetrix_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE symmetrix_api,                   ONLY: symmetrix_eval,&
                                              symmetrix_model_atomic_types,&
                                              symmetrix_model_load,&
                                              symmetrix_model_num_types
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: symmetrix_energy_store_force_virial, &
             symmetrix_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_symmetrix'

   TYPE symmetrix_work_type
      INTEGER                               :: n_atoms_use = 0
      INTEGER                               :: n_atoms_local = 0    ! owned nodes:  owned-first [0, n_atoms_local)
      INTEGER                               :: n_atoms_layer0 = 0   ! owned+ghosts: [0, n_atoms_layer0), get layer 0
      INTEGER                               :: n_atoms_node = 0     ! owned+ghosts+outer: the node arrays' extent
      LOGICAL                               :: use_virial = .FALSE.
      TYPE(cell_type), POINTER              :: cell => NULL()
      TYPE(pos_type), DIMENSION(:), POINTER :: r_pbc => NULL()
      TYPE(particle_type), POINTER          :: particle_set(:) => NULL()
      TYPE(mp_para_env_type), POINTER       :: para_env => NULL()
      LOGICAL, ALLOCATABLE                  :: owned(:)             ! global-indexed: owned by this rank
      INTEGER, ALLOCATABLE                  :: kind_z(:)            ! Z per kind, -1 if not covered
      INTEGER, ALLOCATABLE                  :: global_to_packed(:)  ! global -> 0-based CANONICAL packed, -1 if unused
      INTEGER, ALLOCATABLE                  :: use_indices(:)       ! canonical packed -> global
      ! Canonical packing is ownership-INDEPENDENT and therefore identical on
      ! every rank - it has to be, because the packed edge indices are exchanged
      ! between ranks in gather_edges. The owned-first node order the evaluator
      ! needs is a purely local relabelling applied afterwards:
      INTEGER, ALLOCATABLE                  :: node_perm(:)         ! canonical 0-based -> owned-first 0-based
      INTEGER, ALLOCATABLE                  :: node_perm_inv(:)     ! owned-first 0-based -> canonical 0-based
      ! Per-rank edges stay in default INTEGER so the shared gather in
      ! manybody_mlp_neighbor can be used; the C kinds appear on the
      ! replicated arrays below, which are what crosses the C boundary.
      INTEGER, ALLOCATABLE                  :: local_pairs(:, :)    ! (2, nedges) packed 0-based (centre, neigh)
      INTEGER, ALLOCATABLE                  :: local_shifts(:, :)   ! (3, nedges) integer cell shift
      INTEGER(C_INT32_T), ALLOCATABLE       :: pairs(:, :)          ! all-rank gathered
      INTEGER(C_INT32_T), ALLOCATABLE       :: shifts(:, :)
   END TYPE symmetrix_work_type

CONTAINS

! **************************************************************************************************
!> \brief Evaluate the symmetrix potential, storing forces on the env.
!> \param nonbonded ...
!> \param particle_set ...
!> \param local_particles ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param r_last_update_pbc ...
!> \param pot_total ...
!> \param fist_nonbond_env ...
!> \param para_env ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE symmetrix_energy_store_force_virial(nonbonded, particle_set, local_particles, cell, &
                                                  atomic_kind_set, potparm, r_last_update_pbc, &
                                                  pot_total, fist_nonbond_env, para_env, use_virial)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      REAL(kind=dp), INTENT(OUT)                         :: pot_total
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: use_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'symmetrix_energy_store_force_virial'

      INTEGER                                            :: handle
      TYPE(symmetrix_data_type), POINTER                 :: sym_data
      TYPE(symmetrix_pot_type), POINTER                  :: sym_pot
      TYPE(symmetrix_work_type)                          :: work

      CALL timeset(routineN, handle)
      ! Ownership for the MLP evaluation is decided spatially (see
      ! build_owned_mask); CP2K's own particle distribution is left untouched and
      ! is not used here. The argument stays for interface compatibility with the
      ! other manybody drivers.
      MARK_USED(local_particles)

      pot_total = 0.0_dp
      CALL get_potential_config(atomic_kind_set, potparm, sym_pot)
      IF (.NOT. ASSOCIATED(sym_pot)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      work%use_virial = use_virial
      work%cell => cell
      work%r_pbc => r_last_update_pbc
      work%particle_set => particle_set
      work%para_env => para_env

      CALL build_kind_coverage(atomic_kind_set, potparm, work)
      CALL build_packed_atoms(work)
      CALL build_owned_mask(work)
      CALL build_local_edges(nonbonded, potparm, work)
      CALL gather_edges(work)
      ! must follow gather_edges: the ghost/outer tiers are discovered by walking
      ! the complete (gathered) edge list out from the owned atoms
      CALL build_node_permutation(work)
      CALL setup_symmetrix_data(fist_nonbond_env, sym_data, sym_pot, work)
      CALL run_symmetrix_model(sym_data, work, pot_total)

      CALL release_work(work)
      CALL timestop(handle)
   END SUBROUTINE symmetrix_energy_store_force_virial

! **************************************************************************************************
!> \brief Locate the symmetrix potential among the registered pair potentials.
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param sym_pot ...
! **************************************************************************************************
   SUBROUTINE get_potential_config(atomic_kind_set, potparm, sym_pot)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(symmetrix_pot_type), POINTER                  :: sym_pot

      INTEGER                                            :: i, ikind, jkind
      TYPE(pair_potential_single_type), POINTER          :: pot

      NULLIFY (sym_pot)
      OuterLoop: DO ikind = 1, SIZE(atomic_kind_set)
         DO jkind = ikind, SIZE(atomic_kind_set)
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == symmetrix_type) THEN
                  sym_pot => pot%set(i)%symmetrix
                  EXIT OuterLoop
               END IF
            END DO
         END DO
      END DO OuterLoop
   END SUBROUTINE get_potential_config

! **************************************************************************************************
!> \brief Mark which kinds the model covers and record their atomic numbers.
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_kind_coverage(atomic_kind_set, potparm, work)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER                                            :: i, ikind, jkind, n_kinds, z
      TYPE(pair_potential_single_type), POINTER          :: pot

      n_kinds = SIZE(atomic_kind_set)
      ALLOCATE (work%kind_z(n_kinds))
      work%kind_z = -1

      DO ikind = 1, n_kinds
         DO jkind = 1, n_kinds
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == symmetrix_type) THEN
                  CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
                  work%kind_z(ikind) = z
                  CALL get_atomic_kind(atomic_kind_set(jkind), z=z)
                  work%kind_z(jkind) = z
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE build_kind_coverage

! **************************************************************************************************
!> \brief Flags the atoms this rank owns, using a spatial slab decomposition.
!>
!>        Atom decomposition: each rank evaluates the model only for the atoms
!>        it owns (the other ranks' atoms still take part as ghosts, completing
!>        the owned atoms' environments). Ownership must partition the atoms
!>        exactly once across ranks, which mlp_slab_owner guarantees (one
!>        contiguous spatial slab per rank; see manybody_mlp_partition).
!>
!>        With no distribution (serial, or FIST without one) every atom is
!>        owned, which reduces the evaluation to the plain whole-system one.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_owned_mask(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER                                            :: ic, n_cand
      INTEGER, ALLOCATABLE                               :: candidates(:)
      LOGICAL, ALLOCATABLE                               :: cand_owned(:)

      ALLOCATE (work%owned(SIZE(work%particle_set)))
      work%owned = .FALSE.

      ! Ownership is decided SPATIALLY, not from CP2K's local_particles. The
      ! latter bin-packs molecules by workload, so a rank's atoms are scattered
      ! and their ghost shell is the whole box - correct but with no work saved
      ! in the first layer. A contiguous slab keeps the ghost shell thin. This
      ! map is private to this evaluation and never touches CP2K's own
      ! distribution; see manybody_mlp_partition.
      n_cand = work%n_atoms_use
      ALLOCATE (candidates(n_cand), cand_owned(n_cand))
      DO ic = 1, n_cand
         candidates(ic) = work%use_indices(ic)
      END DO
      CALL mlp_slab_owner(candidates, work%r_pbc, work%cell, work%para_env, cand_owned)
      DO ic = 1, n_cand
         work%owned(candidates(ic)) = cand_owned(ic)
      END DO
      DEALLOCATE (candidates, cand_owned)
   END SUBROUTINE build_owned_mask

! **************************************************************************************************
!> \brief Maps the covered atoms onto CANONICAL packed 0-based node indices.
!>
!>        Canonical = ascending global atom order, independent of which rank owns
!>        what. This is REQUIRED: gather_edges exchanges packed edge indices
!>        between ranks, so every rank must agree on the numbering. (Making the
!>        packing owned-first here instead silently corrupts the exchange - each
!>        rank would interpret its peers' indices in its own numbering.)
!>        The owned-first order the evaluator needs is a local relabelling built
!>        by build_node_permutation and applied after the gather.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_packed_atoms(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER                                            :: iat, iat_use, ikind, n_atoms

      n_atoms = SIZE(work%particle_set)
      ALLOCATE (work%global_to_packed(n_atoms))
      work%global_to_packed = -1

      iat_use = 0
      DO iat = 1, n_atoms
         ikind = work%particle_set(iat)%atomic_kind%kind_number
         IF (work%kind_z(ikind) > 0) THEN
            iat_use = iat_use + 1
            work%global_to_packed(iat) = iat_use - 1
         END IF
      END DO
      work%n_atoms_use = iat_use

      ALLOCATE (work%use_indices(work%n_atoms_use))
      iat_use = 0
      DO iat = 1, n_atoms
         IF (work%global_to_packed(iat) >= 0) THEN
            iat_use = iat_use + 1
            work%use_indices(iat_use) = iat
         END IF
      END DO
   END SUBROUTINE build_packed_atoms

! **************************************************************************************************
!> \brief Builds the local canonical -> owned-first node relabelling.
!>
!>        The evaluator wants this rank's owned nodes contiguous at the front, so
!>        that layer 1 and the readouts can run over just the first
!>        n_atoms_local nodes while layer 0 runs over all of them. This
!>        permutation is applied when the CSR is built (run_symmetrix_model) and
!>        undone when forces are scattered back - it never crosses a rank
!>        boundary, so it is free to differ between ranks.
!>        At P=1 (or with no distribution) every node is owned and this is the
!>        identity, which is what makes the serial path bit-identical.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_node_permutation(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER, PARAMETER                                 :: tier_ghost = 1, tier_outer = 2, &
                                                            tier_owned = 0

      INTEGER                                            :: c0, e, iu, j0, next
      INTEGER, ALLOCATABLE                               :: tier(:)

      ! Classify every canonical node into one of three tiers:
      !   owned  - this rank's slab; gets layer 1 and the readouts
      !   ghost  - within one hop of an owned atom; needs layer 0 so that the
      !            owned atoms' layer-1 features are exact
      !   outer  - only reachable in two hops. A ghost needs its OWN complete
      !            edge list for a correct H1, and those edges land here. Outer
      !            atoms are edge targets only: they carry a position and a type
      !            and they receive forces, but no features are ever computed for
      !            them. (LAMMPS buys these with a 2*r_cut communication cutoff;
      !            here they are free, because FIST replicates every coordinate.)
      ! Atoms in none of the tiers are dropped entirely - that is the whole point
      ! of the spatial partition, and where the work saving comes from.
      ALLOCATE (tier(work%n_atoms_use))
      tier = -1
      DO iu = 1, work%n_atoms_use
         IF (work%owned(work%use_indices(iu))) tier(iu) = tier_owned
      END DO
      DO e = 1, SIZE(work%pairs, 2)
         c0 = INT(work%pairs(1, e)) + 1
         j0 = INT(work%pairs(2, e)) + 1
         IF (tier(c0) == tier_owned .AND. tier(j0) < 0) tier(j0) = tier_ghost
      END DO
      DO e = 1, SIZE(work%pairs, 2)
         c0 = INT(work%pairs(1, e)) + 1
         j0 = INT(work%pairs(2, e)) + 1
         IF (tier(c0) == tier_ghost .AND. tier(j0) < 0) tier(j0) = tier_outer
      END DO

      ALLOCATE (work%node_perm(work%n_atoms_use))
      work%node_perm = -1
      next = 0
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == tier_owned) THEN
            work%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      work%n_atoms_local = next
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == tier_ghost) THEN
            work%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      work%n_atoms_layer0 = next
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == tier_outer) THEN
            work%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      work%n_atoms_node = next

      ALLOCATE (work%node_perm_inv(work%n_atoms_node))
      DO iu = 1, work%n_atoms_use
         IF (work%node_perm(iu) >= 0) work%node_perm_inv(work%node_perm(iu) + 1) = iu - 1
      END DO
      DEALLOCATE (tier)
   END SUBROUTINE build_node_permutation

! **************************************************************************************************
!> \brief Walk this rank's neighbour list into packed 0-based edges + cell shifts.
!> \param nonbonded ...
!> \param potparm ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_local_edges(nonbonded, potparm, work)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER                                            :: e, i, ikind, jkind, nedges, nkeep, &
                                                            nkinds, pa, pb
      INTEGER, ALLOCATABLE                               :: g_edges(:, :), g_shifts(:, :)
      LOGICAL, ALLOCATABLE                               :: kind_included(:)
      REAL(kind=dp), ALLOCATABLE                         :: cutsq(:, :)
      TYPE(pair_potential_single_type), POINTER          :: pot

      ! Resolve the cutoff from the symmetrix set of each kind pair, not from the
      ! pair's overall rcutsq: with a mixed force field the latter is the maximum
      ! over every potential on that pair and would over-collect edges.
      nkinds = SIZE(work%kind_z)
      ALLOCATE (kind_included(nkinds), cutsq(nkinds, nkinds))
      kind_included(:) = work%kind_z(:) >= 0
      cutsq(:, :) = 0.0_dp
      DO ikind = 1, nkinds
         DO jkind = 1, nkinds
            pot => potparm%pot(ikind, jkind)%pot
            IF (.NOT. ASSOCIATED(pot)) CYCLE
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == symmetrix_type) THEN
                  cutsq(ikind, jkind) = pot%set(i)%symmetrix%rcutsq
                  EXIT
               END IF
            END DO
         END DO
      END DO

      CALL mlp_collect_edges(nonbonded, potparm, work%cell, work%r_pbc, symmetrix_type, &
                             kind_included, cutsq, g_edges, g_shifts, nedges)

      ! Translate the global atom indices into the canonical packed numbering the
      ! ranks agreed on; drop edges touching an atom the model does not cover.
      ALLOCATE (work%local_pairs(2, nedges), work%local_shifts(3, nedges))
      nkeep = 0
      DO e = 1, nedges
         pa = work%global_to_packed(g_edges(1, e))
         pb = work%global_to_packed(g_edges(2, e))
         IF (pa < 0 .OR. pb < 0) CYCLE
         nkeep = nkeep + 1
         work%local_pairs(:, nkeep) = [pa, pb]
         work%local_shifts(:, nkeep) = g_shifts(:, e)
      END DO

      IF (nkeep < nedges) THEN
         BLOCK
            INTEGER, ALLOCATABLE :: tp(:, :), ts(:, :)
            ALLOCATE (tp(2, nkeep), ts(3, nkeep))
            tp(:, :) = work%local_pairs(:, 1:nkeep)
            ts(:, :) = work%local_shifts(:, 1:nkeep)
            CALL MOVE_ALLOC(tp, work%local_pairs)
            CALL MOVE_ALLOC(ts, work%local_shifts)
         END BLOCK
      END IF
   END SUBROUTINE build_local_edges

! **************************************************************************************************
!> \brief Gather every rank's edges so all ranks evaluate the full system.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE gather_edges(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER, ALLOCATABLE                               :: all_pairs(:, :), all_shifts(:, :)

      CALL mlp_allgather_edges(work%para_env, work%local_pairs, work%local_shifts, &
                               all_pairs, all_shifts)

      ! The replicated lists are what crosses the C boundary, so narrow them to
      ! the interface's 32-bit kind here (the guard in symmetrix_eval checks the
      ! edge count against HUGE(1_C_INT)).
      ALLOCATE (work%pairs(2, SIZE(all_pairs, 2)), work%shifts(3, SIZE(all_shifts, 2)))
      work%pairs(:, :) = INT(all_pairs, C_INT32_T)
      work%shifts(:, :) = INT(all_shifts, C_INT32_T)
   END SUBROUTINE gather_edges

! **************************************************************************************************
!> \brief Load (once) and cache the model on the env; allocate force scratch.
!> \param fist_nonbond_env ...
!> \param sym_data ...
!> \param sym_pot ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE setup_symmetrix_data(fist_nonbond_env, sym_data, sym_pot, work)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(symmetrix_data_type), POINTER                 :: sym_data
      TYPE(symmetrix_pot_type), POINTER                  :: sym_pot
      TYPE(symmetrix_work_type), INTENT(IN)              :: work

      INTEGER                                            :: e, num_types
      INTEGER(C_INT32_T), ALLOCATABLE                    :: model_z(:)

      CALL fist_nonbond_env_get(fist_nonbond_env, symmetrix_data=sym_data)

      IF (.NOT. ASSOCIATED(sym_data)) THEN
         ALLOCATE (sym_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, symmetrix_data=sym_data)
         NULLIFY (sym_data%use_indices, sym_data%force)
         CALL cite_reference(Batatia2022)
         CALL symmetrix_model_load(sym_data%model, sym_pot%pot_file_name)
         ! map atomic number Z -> symmetrix 0-based type index, once per model.
         ! The lookup covers both the model's elements and any element a covered
         ! CP2K kind carries (a mismatch is caught at evaluation, not by
         ! overflow); the kind coverage is fixed for the lifetime of the env.
         num_types = symmetrix_model_num_types(sym_data%model)
         ALLOCATE (model_z(num_types))
         CALL symmetrix_model_atomic_types(sym_data%model, model_z)
         ALLOCATE (sym_data%idx_of_z(0:MAX(MAXVAL(INT(model_z)), MAXVAL(work%kind_z))))
         sym_data%idx_of_z = -1
         DO e = 1, num_types
            sym_data%idx_of_z(INT(model_z(e))) = e - 1
         END DO
         DEALLOCATE (model_z)
      END IF

      IF (ASSOCIATED(sym_data%force)) THEN
         IF (SIZE(sym_data%force, 2) /= work%n_atoms_use) THEN
            DEALLOCATE (sym_data%force, sym_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(sym_data%force)) THEN
         ALLOCATE (sym_data%force(3, work%n_atoms_use))
         ALLOCATE (sym_data%use_indices(work%n_atoms_use))
      END IF
      sym_data%use_indices(:) = work%use_indices(:)
   END SUBROUTINE setup_symmetrix_data

! **************************************************************************************************
!> \brief Build the CSR neighbour list, evaluate, and store scaled energy/forces.
!> \param sym_data ...
!> \param work ...
!> \param pot_total ...
! **************************************************************************************************
   SUBROUTINE run_symmetrix_model(sym_data, work, pot_total)
      TYPE(symmetrix_data_type), POINTER                 :: sym_data
      TYPE(symmetrix_work_type), INTENT(IN)              :: work
      REAL(kind=dp)                                      :: pot_total

      INTEGER                                            :: c0, e, iat, iu, j0, n_edges, n_use, &
                                                            slot, z
      INTEGER(C_INT32_T), ALLOCATABLE                    :: edge_centre(:), neigh_indices(:), &
                                                            neigh_types(:), node_types(:), &
                                                            num_neigh(:)
      INTEGER, ALLOCATABLE                               :: cursor(:)
      LOGICAL, ALLOCATABLE                               :: keep(:)
      REAL(kind=dp)                                      :: cell_v(3), disp(3), energy, &
                                                            hmat_ang(3, 3), virial(3, 3)
      REAL(kind=dp), ALLOCATABLE                         :: forces(:, :), pos_ang(:, :), r(:), &
                                                            xyz(:, :)

      pot_total = 0.0_dp
      ! n_atoms_use counts the covered atoms globally, n_atoms_node what THIS
      ! rank holds: with more ranks than atoms a rank owns nothing and tiers to
      ! an empty node set. Its stored partials are still scattered by
      ! symmetrix_add_force_virial, so zero them rather than leaving the
      ! previous evaluation's values behind.
      IF (work%n_atoms_use == 0 .OR. work%n_atoms_node == 0) THEN
         IF (ASSOCIATED(sym_data%force)) sym_data%force(:, :) = 0.0_dp
         sym_data%virial(:, :) = 0.0_dp
         RETURN
      END IF

      n_use = work%n_atoms_node
      n_edges = SIZE(work%pairs, 2)

      ! packed positions (Angstrom) and per-atom symmetrix type index; the
      ! Z -> type-index map was cached on sym_data at model load
      ! Node arrays are built in OWNED-FIRST order: node p holds canonical atom
      ! node_perm_inv(p+1). The evaluator needs this rank's owned nodes
      ! contiguous at the front (see symmetrix_c_eval_mp); the canonical order
      ! is what the ranks agreed on in gather_edges, so every canonical index
      ! read out of work%pairs is translated through node_perm below.
      ALLOCATE (pos_ang(3, n_use), node_types(n_use))
      DO iu = 1, n_use
         iat = work%use_indices(work%node_perm_inv(iu) + 1)
         pos_ang(:, iu) = work%r_pbc(iat)%r(:)*angstrom
         z = work%kind_z(work%particle_set(iat)%atomic_kind%kind_number)
         IF (sym_data%idx_of_z(z) < 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "An atom element is not among the atomic numbers supported by the symmetrix model.")
         END IF
         node_types(iu) = INT(sym_data%idx_of_z(z), C_INT32_T)
      END DO
      hmat_ang(:, :) = work%cell%hmat*angstrom

      ! Keep only the edges this rank needs: those centred on a node that gets
      ! layer 0 (owned or ghost). Edges centred on outer atoms, or on atoms
      ! dropped entirely, belong to other ranks. n_edges is recounted here, so
      ! the arrays below are sized to the kept edges only.
      ALLOCATE (keep(n_edges))
      keep = .FALSE.
      DO e = 1, SIZE(work%pairs, 2)
         c0 = work%node_perm(INT(work%pairs(1, e)) + 1)
         IF (c0 >= 0 .AND. c0 < work%n_atoms_layer0) THEN
            IF (work%node_perm(INT(work%pairs(2, e)) + 1) >= 0) keep(e) = .TRUE.
         END IF
      END DO
      n_edges = COUNT(keep)

      ! count neighbours per central atom for the CSR layout
      ALLOCATE (num_neigh(n_use))
      num_neigh = 0_C_INT32_T
      DO e = 1, SIZE(work%pairs, 2)
         IF (.NOT. keep(e)) CYCLE
         c0 = work%node_perm(INT(work%pairs(1, e)) + 1)
         num_neigh(c0 + 1) = num_neigh(c0 + 1) + 1_C_INT32_T
      END DO

      ! prefix-sum cursor (0-based offsets) for the counting sort
      ALLOCATE (cursor(n_use))
      cursor(1) = 0
      DO iu = 2, n_use
         cursor(iu) = cursor(iu - 1) + INT(num_neigh(iu - 1))
      END DO

      ALLOCATE (neigh_indices(n_edges), neigh_types(n_edges), edge_centre(n_edges))
      ALLOCATE (xyz(3, n_edges), r(n_edges))
      DO e = 1, SIZE(work%pairs, 2)
         IF (.NOT. keep(e)) CYCLE
         c0 = work%node_perm(INT(work%pairs(1, e)) + 1)
         j0 = work%node_perm(INT(work%pairs(2, e)) + 1)
         slot = cursor(c0 + 1) + 1
         cursor(c0 + 1) = cursor(c0 + 1) + 1

         cell_v = MATMUL(hmat_ang, REAL(work%shifts(:, e), dp))
         disp(:) = pos_ang(:, j0 + 1) - pos_ang(:, c0 + 1) + cell_v

         neigh_indices(slot) = INT(j0, C_INT32_T)
         edge_centre(slot) = INT(c0, C_INT32_T)
         neigh_types(slot) = node_types(j0 + 1)
         xyz(:, slot) = disp(:)
         r(slot) = NORM2(disp)
      END DO

      ALLOCATE (forces(3, n_use))
      ! Atom decomposition: this rank evaluates the model only for the atoms it
      ! owns - packed nodes [0, n_atoms_local) - with the rest present as ghosts
      ! so the owned atoms' message-passing environments are complete.
      CALL symmetrix_eval(sym_data%model, node_types, num_neigh, neigh_indices, neigh_types, &
                          edge_centre, xyz, r, work%use_virial, energy, forces, virial, &
                          n_local=work%n_atoms_local, n_layer0=work%n_atoms_layer0)

      ! energy eV->Hartree, forces eV/Angstrom -> Hartree/Bohr, virial eV->Hartree.
      ! NO 1/num_pe: the energy already covers only this rank's owned atoms, and
      ! the forces/virial are already genuine partial derivatives of that partial
      ! energy. CP2K's cross-rank sums in fist_force.F assemble the exact totals.
      ! The virial lands as r (outer) f in CP2K's memory order; this is symmetric
      ! for a rotation-invariant total energy.
      pot_total = energy/evolt
      ! forces come back in owned-first node order; sym_data%force is indexed by
      ! canonical packed node (that is what use_indices/add_force_virial expect)
      sym_data%force(:, :) = 0.0_dp
      DO iu = 1, n_use
         sym_data%force(:, work%node_perm_inv(iu) + 1) = forces(:, iu)/evolt*angstrom
      END DO
      IF (work%use_virial) sym_data%virial(:, :) = virial(:, :)/evolt

      DEALLOCATE (pos_ang, node_types, num_neigh, cursor, keep, &
                  neigh_indices, neigh_types, edge_centre, xyz, r, forces)
   END SUBROUTINE run_symmetrix_model

! **************************************************************************************************
!> \brief Accumulate stored forces into the nonbonded environment.
!> \param fist_nonbond_env ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE symmetrix_add_force_virial(fist_nonbond_env, f_nonbond, pv_nonbond, use_virial)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial

      INTEGER                                            :: iat, iu
      TYPE(symmetrix_data_type), POINTER                 :: sym_data

      CALL fist_nonbond_env_get(fist_nonbond_env, symmetrix_data=sym_data)
      IF (.NOT. ASSOCIATED(sym_data)) RETURN

      IF (use_virial) pv_nonbond = pv_nonbond + sym_data%virial

      DO iu = 1, SIZE(sym_data%use_indices)
         iat = sym_data%use_indices(iu)
         f_nonbond(1:3, iat) = f_nonbond(1:3, iat) + sym_data%force(1:3, iu)
      END DO
   END SUBROUTINE symmetrix_add_force_virial

! **************************************************************************************************
!> \brief Release the per-call work scratch.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE release_work(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      IF (ALLOCATED(work%owned)) DEALLOCATE (work%owned)
      IF (ALLOCATED(work%kind_z)) DEALLOCATE (work%kind_z)
      IF (ALLOCATED(work%global_to_packed)) DEALLOCATE (work%global_to_packed)
      IF (ALLOCATED(work%use_indices)) DEALLOCATE (work%use_indices)
      IF (ALLOCATED(work%node_perm)) DEALLOCATE (work%node_perm)
      IF (ALLOCATED(work%node_perm_inv)) DEALLOCATE (work%node_perm_inv)
      IF (ALLOCATED(work%local_pairs)) DEALLOCATE (work%local_pairs)
      IF (ALLOCATED(work%local_shifts)) DEALLOCATE (work%local_shifts)
      IF (ALLOCATED(work%pairs)) DEALLOCATE (work%pairs)
      IF (ALLOCATED(work%shifts)) DEALLOCATE (work%shifts)
      NULLIFY (work%cell, work%r_pbc, work%particle_set, work%para_env)
   END SUBROUTINE release_work

END MODULE manybody_symmetrix

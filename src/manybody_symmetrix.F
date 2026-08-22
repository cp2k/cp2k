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
!>        The model is ATOM DECOMPOSED across MPI ranks and threads: each domain
!>        owns a compact region of the cell (manybody_mlp_partition) and
!>        evaluates the model only for those atoms, with the atoms around them
!>        present as ghosts so the owned atoms' 2-layer message-passing
!>        environments stay exact. The edges come from a linked-cell search over
!>        the replicated coordinates (manybody_mlp_celllist), so a domain builds
!>        only what it needs and no rank ever holds the whole system's edges. Layer 0 runs over owned+ghosts, layer 1 and the readouts over
!>        the owned nodes only. Because FIST replicates every atom's position on
!>        every rank the ghosts cost no communication at all, and CP2K's own data
!>        structures - including its particle distribution - are left untouched.
!>        Each rank's energy/forces/virial are genuine partial contributions and
!>        are NOT scaled by 1/num_pe; the existing cross-rank sums in
!>        fist_force.F assemble the exact totals. Same contract as the Allegro
!>        path in manybody_nequip.F.
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE manybody_symmetrix

   USE ISO_C_BINDING,                   ONLY: C_INT32_T
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type,&
                                              pos_type,&
                                              symmetrix_data_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_1,&
                                              int_8
   USE manybody_mlp_celllist,           ONLY: &
        mlp_cell_list_build, mlp_cell_list_release, mlp_cell_list_type, mlp_edge_buffer_release, &
        mlp_edge_buffer_reserve, mlp_edge_buffer_type, mlp_expand_tier, mlp_tier_collect, &
        mlp_tier_ghost, mlp_tier_none, mlp_tier_outer, mlp_tier_owned
   USE manybody_mlp_partition,          ONLY: mlp_domain_owner
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                              pair_potential_single_type,&
                                              symmetrix_pot_type,&
                                              symmetrix_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE symmetrix_api,                   ONLY: symmetrix_eval,&
                                              symmetrix_is_kokkos,&
                                              symmetrix_model_atomic_types,&
                                              symmetrix_model_load,&
                                              symmetrix_model_num_types
   USE util,                            ONLY: sort

!$ USE OMP_LIB, ONLY: omp_get_max_threads, &
!$                    omp_get_num_threads, &
!$                    omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: symmetrix_energy_store_force_virial, &
             symmetrix_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_symmetrix'

   ! Cross-checks the locally built edge set against the one derived from CP2K's
   ! own neighbour list. Off by default: it builds the whole system's edges on
   ! every rank every step, which is exactly the cost the local build exists to
   ! avoid. Turn it on when touching the neighbour search.
   ! NOTE: flipping this needs the object file removed, not just a rebuild.
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   ! The node tiers and the owned-first relabelling are per SLAB, and with more
   ! than one thread a rank holds several. They are kept apart from the shared
   ! work state so that a thread's copy is plainly thread-private.
   TYPE symmetrix_thread_type
      INTEGER                               :: n_atoms_local = 0    ! owned nodes:  owned-first [0, n_atoms_local)
      INTEGER                               :: n_atoms_layer0 = 0   ! owned+ghosts: [0, n_atoms_layer0), get layer 0
      INTEGER                               :: n_atoms_node = 0     ! owned+ghosts+outer: the node arrays' extent
      INTEGER, ALLOCATABLE                  :: node_perm(:)         ! canonical 0-based -> owned-first 0-based
      INTEGER, ALLOCATABLE                  :: node_perm_inv(:)     ! owned-first 0-based -> canonical 0-based
      INTEGER(C_INT32_T), ALLOCATABLE       :: num_neigh(:)         ! per node, CSR counts
      REAL(KIND=dp), ALLOCATABLE            :: force_node(:, :)     ! (3, n_atoms_node), this domain's order
   END TYPE symmetrix_thread_type

   TYPE symmetrix_work_type
      INTEGER                               :: n_atoms_use = 0
      LOGICAL                               :: use_virial = .FALSE.
      TYPE(cell_type), POINTER              :: cell => NULL()
      TYPE(pos_type), DIMENSION(:), POINTER :: r_pbc => NULL()
      TYPE(particle_type), POINTER          :: particle_set(:) => NULL()
      TYPE(mp_para_env_type), POINTER       :: para_env => NULL()
      INTEGER                               :: n_threads = 1        ! domains this rank splits into
      REAL(KIND=dp)                         :: rcut = 0.0_dp        ! model cutoff, bohr
      REAL(KIND=dp)                         :: rcutsq = 0.0_dp      ! verbatim from the potential
      INTEGER, ALLOCATABLE                  :: thread_of(:)         ! global-indexed: which of this rank's slabs, -1 if none
      INTEGER, ALLOCATABLE                  :: kind_z(:)            ! Z per kind, -1 if not covered
      INTEGER, ALLOCATABLE                  :: global_to_packed(:)  ! global -> 0-based CANONICAL packed, -1 if unused
      INTEGER, ALLOCATABLE                  :: use_indices(:)       ! canonical packed -> global
      ! Canonical packing is ownership-INDEPENDENT and therefore identical on
      ! every rank, which is what lets a ghost's edge list match the one built
      ! where that atom is owned. The owned-first node order the evaluator needs
      ! is a purely local relabelling, held per domain in symmetrix_thread_type.
      REAL(KIND=dp), ALLOCATABLE            :: packed_pos(:, :)     ! (3, n_atoms_use), bohr
      ! Built once per rank per evaluation and read by every thread. Shared and
      ! never written inside the parallel region.
      TYPE(mlp_cell_list_type)              :: cl
   END TYPE symmetrix_work_type

CONTAINS

! **************************************************************************************************
!> \brief Evaluate the symmetrix potential, storing forces on the env.
!> \param nonbonded ...
!> \param particle_set ...
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

      INTEGER                                            :: handle, handle_stage, ikind, it, iu, &
                                                            my_thread
      REAL(kind=dp), ALLOCATABLE                         :: pot_thr(:), virial_thr(:, :, :)
      TYPE(symmetrix_data_type), POINTER                 :: sym_data
      TYPE(symmetrix_pot_type), POINTER                  :: sym_pot
      TYPE(symmetrix_thread_type), ALLOCATABLE           :: thr_all(:)
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

      ! With a core/shell model the FIST force path takes a branch that never
      ! reaches force_nonbond_manybody, so the manybody forces stored here would
      ! be dropped while the energy still counted - a silent inconsistency
      ! rather than an error. Refuse the combination instead of returning
      ! forces that do not match the reported energy.
      DO ikind = 1, SIZE(atomic_kind_set)
         IF (atomic_kind_set(ikind)%shell_active) THEN
            CALL cp_abort(__LOCATION__, &
                          "The symmetrix MACE backend does not support the core/shell model: "// &
                          "the manybody force contribution is not accumulated on that path, so "// &
                          "the forces would not correspond to the reported energy.")
         END IF
      END DO

      work%use_virial = use_virial
      work%cell => cell
      work%r_pbc => r_last_update_pbc
      work%particle_set => particle_set
      work%para_env => para_env
      ! The decomposition needs the cutoff to weigh one domain shape against
      ! another: it is the depth of the ghost shell every domain has to carry.
      ! The squared form is kept verbatim: it is what the edge cutoff test
      ! compares against, and a round trip through SQRT would move the boundary.
      work%rcutsq = sym_pot%rcutsq
      work%rcut = SQRT(work%rcutsq)

      ! One slab per thread, decomposed exactly as ranks are. The evaluator keeps
      ! its scratch inside the model object, so a thread needs a model of its
      ! own; the count is fixed here and the clones are made once at load.
      work%n_threads = 1
!$    work%n_threads = omp_get_max_threads()
      ! On a GPU build the evaluator shares one execution space, so running the
      ! domains concurrently would have every thread reallocate and launch into
      ! it at once, and would put a full model clone per thread in device
      ! memory. One domain per rank there; the GPU is the parallelism.
      IF (symmetrix_is_kokkos()) work%n_threads = 1
      ! The per-thread state on the environment is sized once and kept, so the
      ! count must not move between evaluations. It cannot within one run, but
      ! indexing past those arrays would be silent, so it is checked.
      IF (ASSOCIATED(fist_nonbond_env%symmetrix_data)) THEN
         IF (ALLOCATED(fist_nonbond_env%symmetrix_data%model)) THEN
            IF (SIZE(fist_nonbond_env%symmetrix_data%model) /= work%n_threads) THEN
               CALL cp_abort(__LOCATION__, &
                             "The OpenMP thread count changed between force evaluations. The "// &
                             "symmetrix model clones and edge buffers are allocated once, per "// &
                             "thread, and cannot follow it.")
            END IF
         END IF
      END IF

      ! Each stage carries its own timer. The stage breakdown is what tells a
      ! rebuild cost apart from the model cost: any stage whose self time does
      ! not fall with rank count, or does not grow with the atom count, is doing
      ! work that does not belong on the per-step path.
      CALL timeset('symmetrix_setup', handle_stage)
      CALL build_kind_coverage(atomic_kind_set, potparm, work)
      CALL build_packed_atoms(work)
      CALL timestop(handle_stage)

      CALL timeset('symmetrix_partition', handle_stage)
      CALL build_owned_mask(work)
      CALL timestop(handle_stage)

      ! One grid per rank, shared by every thread. The domains then walk out of
      ! it only from the atoms they own, so no rank ever holds, or exchanges,
      ! the whole system's edge list.
      CALL timeset('symmetrix_cells', handle_stage)
      ALLOCATE (work%packed_pos(3, work%n_atoms_use))
      DO ikind = 1, work%n_atoms_use
         work%packed_pos(:, ikind) = work%r_pbc(work%use_indices(ikind))%r(:)
      END DO
      CALL mlp_cell_list_build(work%cl, work%packed_pos, work%cell, work%rcutsq)
      CALL timestop(handle_stage)

      IF (debug_this_module) CALL check_edges_against_cell_list(work, nonbonded, potparm)

      CALL setup_symmetrix_data(fist_nonbond_env, sym_data, sym_pot, work)

      ! Evaluate the slabs. Each thread owns one, exactly as a rank does, and
      ! carries its own tiers, CSR and model clone; the gathered edge list and
      ! the packed atom maps are read-only here and shared. Results land in
      ! per-thread accumulators and are folded below, so no thread writes to
      ! anything another thread reads.
      ! Per-domain state lives in a shared array indexed by thread rather than in
      ! a private variable. A private copy of a derived type whose components are
      ! themselves allocatable is exactly the case where the per-thread copy is
      ! easy to get wrong, and every entry here is touched by one thread only,
      ! so sharing the array costs nothing and removes the question.
      !
      ! Forces stay on the domain that computed them, in ITS node order, and are
      ! scattered once during the fold. A per-thread copy of the whole packed
      ! force array instead would be 3*N*T reals, which at half a million atoms
      ! on 64 threads is 783 MB per rank and would be the largest allocation
      ! left now that the gathered edge list is gone.
      ALLOCATE (thr_all(work%n_threads))
      ALLOCATE (pot_thr(work%n_threads))
      ALLOCATE (virial_thr(3, 3, work%n_threads))
      ! Zeroed because the fold below reads every entry. A team smaller than the
      ! domain count would otherwise leave some of them untouched, and the fold
      ! would add whatever the allocation happened to contain.
      pot_thr(:) = 0.0_dp
      virial_thr(:, :, :) = 0.0_dp

      CALL timeset('symmetrix_domains', handle_stage)
      ! The cell was cut into exactly work%n_threads domains, so exactly that
      ! many have to be evaluated. Asking for the count is not the same as
      ! getting it: OpenMP may hand back a smaller team, and a domain nobody
      ! runs is a set of atoms silently missing from the energy rather than an
      ! error. Request the count, then check it was granted.
!$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(work%n_threads) &
!$OMP          SHARED(work, sym_data, pot_thr, virial_thr, thr_all) &
!$OMP          PRIVATE(my_thread)
      my_thread = 0
!$    my_thread = omp_get_thread_num()
!$OMP MASTER
!$    IF (omp_get_num_threads() /= work%n_threads) THEN
!$       CALL cp_abort(__LOCATION__, &
!$                     "The OpenMP runtime provided fewer threads than the atom decomposition "// &
!$                     "was built for, which would leave whole domains unevaluated.")
!$    END IF
!$OMP END MASTER
!$OMP BARRIER
      CALL build_domain_edges(work, my_thread, thr_all(my_thread + 1), &
                              sym_data%edge_buf(my_thread + 1))
      CALL run_symmetrix_model(sym_data, work, thr_all(my_thread + 1), my_thread, &
                               pot_thr(my_thread + 1), &
                               virial_thr(:, :, my_thread + 1))
!$OMP END PARALLEL
      CALL timestop(handle_stage)

      ! Fold in ascending thread order, so a given thread count reproduces
      ! bitwise from run to run rather than in whatever order the threads finish.
      pot_total = 0.0_dp
      sym_data%force(:, :) = 0.0_dp
      sym_data%virial(:, :) = 0.0_dp
      DO it = 1, work%n_threads
         pot_total = pot_total + pot_thr(it)
         DO iu = 1, thr_all(it)%n_atoms_node
            sym_data%force(:, thr_all(it)%node_perm_inv(iu) + 1) = &
               sym_data%force(:, thr_all(it)%node_perm_inv(iu) + 1) + thr_all(it)%force_node(:, iu)
         END DO
         IF (use_virial) sym_data%virial(:, :) = sym_data%virial(:, :) + virial_thr(:, :, it)
      END DO
      DO it = 1, work%n_threads
         CALL release_thread(thr_all(it))
      END DO
      DEALLOCATE (thr_all, pot_thr, virial_thr)

      CALL release_work(work)
      CALL timestop(handle)
   END SUBROUTINE symmetrix_energy_store_force_virial

! **************************************************************************************************
!> \brief Locate the symmetrix potential among the registered pair potentials.
!>
!>        The whole sweep is walked rather than stopping at the first hit,
!>        because ONE model and ONE cutoff are used for every covered atom: the
!>        cutoff sizes the neighbour grid and the stencil, and the model is
!>        loaded once and cached. A second, different model would therefore be
!>        silently ignored and its atoms evaluated with the first model at the
!>        first model's range. Two routes reach that state - &SYMMETRIX and
!>        &MACE are read into separate potential slots, and MULTIPLE_POTENTIAL
!>        adds rather than replaces - so it is refused here rather than trusted
!>        not to happen.
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
      DO ikind = 1, SIZE(atomic_kind_set)
         DO jkind = ikind, SIZE(atomic_kind_set)
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= symmetrix_type) CYCLE
               IF (.NOT. ASSOCIATED(sym_pot)) THEN
                  sym_pot => pot%set(i)%symmetrix
               ELSE IF (pot%set(i)%symmetrix%pot_file_name /= sym_pot%pot_file_name .OR. &
                        pot%set(i)%symmetrix%rcutsq /= sym_pot%rcutsq) THEN
                  CALL cp_abort(__LOCATION__, &
                                "More than one symmetrix model is defined in this force field: '"// &
                                TRIM(sym_pot%pot_file_name)//"' and '"// &
                                TRIM(pot%set(i)%symmetrix%pot_file_name)//"'. The backend "// &
                                "evaluates a single model over every atom it covers, so only "// &
                                "the first would be used. Give one model per force field, in "// &
                                "one &MACE section.")
               END IF
            END DO
         END DO
      END DO
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

      CHARACTER(LEN=default_string_length)               :: kind_name
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

      ! Under QM/MM the potential is stripped from pairs whose kinds are both QM,
      ! but a pair with one QM and one MM kind keeps it, and the loop above then
      ! marks the QM kind covered. The model would evaluate over QM atoms whose
      ! energy the electronic structure already accounts for, and would build
      ! QM-QM edges that the FIST neighbour list deliberately leaves out. Refuse
      ! the combination rather than return that.
      DO ikind = 1, n_kinds
         IF (work%kind_z(ikind) < 0) CYCLE
         CALL get_atomic_kind(atomic_kind_set(ikind), name=kind_name)
         IF (INDEX(kind_name, "_QM_") /= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "The symmetrix MACE backend does not support QM/MM: kind '"// &
                          TRIM(kind_name)//"' is quantum, and the model would evaluate it "// &
                          "alongside the electronic structure that already describes it.")
         END IF
      END DO
   END SUBROUTINE build_kind_coverage

! **************************************************************************************************
!> \brief Build the packed list of model atoms (global<->packed index maps).
!> \param work ...
! **************************************************************************************************
! **************************************************************************************************
!> \brief Flags the atoms this rank owns, from CP2K's particle distribution.
!>
!>        Atom decomposition: each rank evaluates the model only for the atoms
!>        it owns (the other ranks' atoms still take part as ghosts, completing
!>        the owned atoms' environments). Ownership must partition the atoms
!>        exactly once across ranks, which is precisely what local_particles
!>        guarantees. Same construction as the Allegro path in manybody_nequip.F.
!>
!>        With no distribution (serial, or FIST without one) every atom is
!>        owned, which reduces the evaluation to the plain whole-system one.
!> \param local_particles the atoms owned by this rank, per kind
!> \param work ...
! **************************************************************************************************
   SUBROUTINE build_owned_mask(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      INTEGER                                            :: ic, n_cand
      INTEGER, ALLOCATABLE                               :: cand_thread(:), candidates(:)
      LOGICAL, ALLOCATABLE                               :: cand_owned(:)

      ALLOCATE (work%thread_of(SIZE(work%particle_set)))
      work%thread_of = -1

      ! Ownership is decided SPATIALLY, not from CP2K's local_particles. The
      ! latter bin-packs molecules by workload, so a rank's atoms are scattered
      ! and their ghost shell is the whole box - correct but with no work saved
      ! in the first layer. A compact region keeps the ghost shell thin. This
      ! map is private to this evaluation and never touches CP2K's own
      ! distribution; see manybody_mlp_partition.
      n_cand = work%n_atoms_use
      ALLOCATE (candidates(n_cand), cand_owned(n_cand), cand_thread(n_cand))
      DO ic = 1, n_cand
         candidates(ic) = work%use_indices(ic)
      END DO
      ! Ranks and threads are partitioned together here, in one call, because the
      ! MPI consistency check inside has to happen outside any parallel region.
      CALL mlp_domain_owner(candidates, work%r_pbc, work%cell, work%para_env, &
                            work%n_threads, work%rcut, cand_owned, cand_thread)
      DO ic = 1, n_cand
         work%thread_of(candidates(ic)) = cand_thread(ic)
      END DO
      DEALLOCATE (candidates, cand_owned, cand_thread)
   END SUBROUTINE build_owned_mask

! **************************************************************************************************
!> \brief Maps the covered atoms onto CANONICAL packed 0-based node indices.
!>
!>        Canonical = ascending global atom order, independent of which rank owns
!>        what. This is REQUIRED: the edge list a domain builds for one of its
!>        ghosts has to come out identical to the one built where that atom is
!>        owned, and that only holds if every domain numbers the atoms the same
!>        way. The owned-first order the evaluator needs is a purely local
!>        relabelling, built by build_domain_edges once the tiers are known.
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
!> \param my_thread the slab to build the relabelling for (0-based within the rank)
!> \param thr per-slab tiers and relabelling, filled here
! **************************************************************************************************
   SUBROUTINE build_domain_edges(work, my_thread, thr, buf)
      TYPE(symmetrix_work_type), INTENT(IN)              :: work
      INTEGER, INTENT(IN)                                :: my_thread
      TYPE(symmetrix_thread_type), INTENT(INOUT)         :: thr
      TYPE(mlp_edge_buffer_type), INTENT(INOUT)          :: buf

      INTEGER                                            :: e, iu, n_ghost, n_owned, next
      INTEGER(KIND=int_1), ALLOCATABLE                   :: tier(:)
      INTEGER, ALLOCATABLE                               :: ghosts(:), owned(:)

      ! Classify every canonical node into one of three tiers:
      !   owned  - this domain's atoms; get layer 1 and the readouts
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
      !
      ! The tiers are exactly a two-step breadth-first walk out from the owned
      ! set, so they are found by walking the neighbour search rather than by
      ! sweeping a list of the whole system's edges. Everything here is
      ! thread-private; the cell list it reads is shared and never written.
      ALLOCATE (tier(work%n_atoms_use))
      ALLOCATE (owned(work%n_atoms_use), ghosts(work%n_atoms_use))

      ! The two walks below both discover tiers and record edges. A walk that
      ! runs out of room still counts correctly, so an undersized buffer costs
      ! one repeat and nothing else. Buffers persist across evaluations and are
      ! never shrunk, so the repeat happens on the first step and then stops:
      ! in the steady state each domain walks its stencil once, not twice.
      DO
         tier = mlp_tier_none
         n_owned = 0
         DO iu = 1, work%n_atoms_use
            IF (work%thread_of(work%use_indices(iu)) == my_thread) THEN
               tier(iu) = mlp_tier_owned
               n_owned = n_owned + 1
               owned(n_owned) = iu - 1
            END IF
         END DO

         buf%n_edges = 0
         CALL mlp_expand_tier(work%cl, owned(1:n_owned), tier, mlp_tier_ghost, buf, .TRUE.)
         CALL mlp_tier_collect(tier, mlp_tier_ghost, ghosts, n_ghost)
         CALL mlp_expand_tier(work%cl, ghosts(1:n_ghost), tier, mlp_tier_outer, buf, .TRUE.)

         IF (buf%n_edges <= buf%capacity) EXIT
         ! A little headroom, because the edge count breathes from step to step
         ! and reserving exactly what was just seen would make the next small
         ! rise walk the stencil twice again. Under NPT compression, or simply
         ! as a liquid rearranges, that is not a first-step event.
         CALL mlp_edge_buffer_reserve(buf, buf%n_edges + buf%n_edges/8 + 64)
      END DO

      ! Owned-first relabelling. Both tier lists are in ascending canonical
      ! order and the node numbering below follows the same order, so owned
      ! entry k is node k and ghost entry k is node n_owned + k. That is what
      ! lets the emitted edges double as the CSR: they come out grouped by
      ! centre, in ascending centre order, with no sort and no counting pass.
      ALLOCATE (thr%node_perm(work%n_atoms_use))
      thr%node_perm = -1
      next = 0
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == mlp_tier_owned) THEN
            thr%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      thr%n_atoms_local = next
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == mlp_tier_ghost) THEN
            thr%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      thr%n_atoms_layer0 = next
      DO iu = 1, work%n_atoms_use
         IF (tier(iu) == mlp_tier_outer) THEN
            thr%node_perm(iu) = next
            next = next + 1
         END IF
      END DO
      thr%n_atoms_node = next
      CPASSERT(n_owned + n_ghost == thr%n_atoms_layer0)

      ALLOCATE (thr%node_perm_inv(thr%n_atoms_node))
      DO iu = 1, work%n_atoms_use
         IF (thr%node_perm(iu) >= 0) thr%node_perm_inv(thr%node_perm(iu) + 1) = iu - 1
      END DO

      ! Neighbour counts per node, recovered from the edges themselves in one
      ! sweep. The centres are canonical packed indices, so they go through the
      ! relabelling: num_neigh is indexed by NODE, which is what the evaluator
      ! reads it as, and the two numberings coincide only when a domain owns
      ! everything.
      ALLOCATE (thr%num_neigh(thr%n_atoms_node))
      thr%num_neigh = 0_C_INT32_T
      DO e = 1, buf%n_edges
         next = thr%node_perm(INT(buf%pairs(1, e)) + 1)
         thr%num_neigh(next + 1) = thr%num_neigh(next + 1) + 1_C_INT32_T
      END DO

      DEALLOCATE (tier, owned, ghosts)
   END SUBROUTINE build_domain_edges

! **************************************************************************************************
!> \brief Checks the cell-list edge set against the one CP2K's neighbour list
!>        gives, as an exact multiset of (centre, neighbour, shift) triples.
!>
!>        This is the acceptance test for the neighbour search, and it is an
!>        integer comparison rather than a numerical one: the two searches share
!>        no code, so agreeing on every edge pins wrapping, aperiodic axes and
!>        boxes shorter than twice the cutoff all at once. CP2K's own list is the
!>        oracle because it is independent, long established, and handles
!>        periodicity correctly.
!>
!>        Both sides are sorted before comparing, since the orders differ by
!>        construction and only the SET has to agree.
!>
!>        NOT VALID WHILE THE CELL MOVES. CP2K's list rebuilds on particle drift
!>        alone and never on a change of cell, so under NPT or CELL_OPT the
!>        reference goes stale against a deforming box while this backend keeps
!>        up with it. The disagreement that follows is the oracle being wrong,
!>        not the search. Compare on a fixed cell.
!> \param work must already hold the gathered whole-system edge list
! **************************************************************************************************
   SUBROUTINE check_edges_against_cell_list(work, nonbonded, potparm)
      TYPE(symmetrix_work_type), INTENT(IN)              :: work
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_edges_against_cell_list'

      CHARACTER(LEN=default_string_length*2)             :: msg
      INTEGER                                            :: atom_a, atom_b, e, handle, i, iend, &
                                                            igrp, ikind, ilist, ipair, istart, &
                                                            jkind, n_cent, n_ref, n_use, npairs, &
                                                            pa, pb
      INTEGER(KIND=int_1), ALLOCATABLE                   :: tier(:)
      INTEGER(KIND=int_8), ALLOCATABLE                   :: key_new(:), key_ref(:)
      INTEGER, ALLOCATABLE                               :: centres(:), perm(:)
      INTEGER, DIMENSION(3)                              :: cvi
      LOGICAL                                            :: do_symmetrix
      LOGICAL, ALLOCATABLE                               :: is_centre(:)
      REAL(kind=dp)                                      :: cell_v(3), cutsq_ij, drij, rij(3)
      TYPE(mlp_edge_buffer_type)                         :: buf
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      CALL timeset(routineN, handle)
      n_use = work%n_atoms_use

      ! Reference side: this rank's slice of CP2K's neighbour list, walked
      ! exactly as the driver used to walk it. Under a full list every edge
      ! centred on an atom is produced by the one rank holding that atom, so the
      ! slice is complete for its own centres and can be compared on its own,
      ! with no exchange.
      ALLOCATE (key_ref(0))
      ALLOCATE (is_centre(n_use))
      is_centre = .FALSE.
      n_ref = 0
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            IF (work%kind_z(ikind) < 0 .OR. work%kind_z(jkind) < 0) CYCLE Kind_Loop
            pot => potparm%pot(ikind, jkind)%pot
            do_symmetrix = .FALSE.
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == symmetrix_type) THEN
                  do_symmetrix = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. do_symmetrix) CYCLE Kind_Loop
            cutsq_ij = pot%set(i)%symmetrix%rcutsq
            cvi = neighbor_kind_pair%cell_vector
            cell_v = MATMUL(work%cell%hmat, REAL(cvi, dp))
            DO ipair = istart, iend
               atom_a = neighbor_kind_pair%list(1, ipair)
               atom_b = neighbor_kind_pair%list(2, ipair)
               rij(:) = work%r_pbc(atom_b)%r(:) - work%r_pbc(atom_a)%r(:) + cell_v
               drij = DOT_PRODUCT(rij, rij)
               IF (drij <= cutsq_ij) THEN
                  pa = work%global_to_packed(atom_a)
                  pb = work%global_to_packed(atom_b)
                  n_ref = n_ref + 1
                  CALL push_key(key_ref, n_ref, edge_key(pa, pb, cvi, n_use))
                  is_centre(pa + 1) = .TRUE.
               END IF
            END DO
         END DO Kind_Loop
      END DO

      ! Comparison side: the same centres, walked out of the cell list.
      ALLOCATE (centres(n_use))
      n_cent = 0
      DO i = 1, n_use
         IF (is_centre(i)) THEN
            n_cent = n_cent + 1
            centres(n_cent) = i - 1
         END IF
      END DO
      ALLOCATE (tier(n_use))
      tier = mlp_tier_owned
      buf%n_edges = 0
      CALL mlp_expand_tier(work%cl, centres(1:n_cent), tier, mlp_tier_ghost, buf, .FALSE.)
      CALL mlp_edge_buffer_reserve(buf, MAX(buf%n_edges, 1))
      buf%n_edges = 0
      CALL mlp_expand_tier(work%cl, centres(1:n_cent), tier, mlp_tier_ghost, buf, .TRUE.)

      IF (buf%n_edges /= n_ref) THEN
         WRITE (msg, '(A,I0,A,I0,A,3(1X,I0),A,3(1X,I0),A,3(1X,I0))') &
            "The MLP cell list found ", buf%n_edges, " edges where CP2K's neighbour list has ", &
            n_ref, ". Grid:", work%cl%nbin, "  stencil:", work%cl%span, &
            "  periodicity:", work%cl%perd
         CALL cp_abort(__LOCATION__, TRIM(msg))
      END IF

      ! One integer per edge, ordered so equal triples collide and nothing else
      ! does. Both sides are sorted before comparing, since their orders differ
      ! by construction and only the SET has to agree.
      ALLOCATE (key_new(n_ref), perm(MAX(n_ref, 1)))
      DO e = 1, n_ref
         key_new(e) = edge_key(INT(buf%pairs(1, e)), INT(buf%pairs(2, e)), &
                               INT(buf%shifts(:, e)), n_use)
      END DO
      IF (n_ref > 0) THEN
         CALL sort(key_ref, n_ref, perm)
         CALL sort(key_new, n_ref, perm)
         DO e = 1, n_ref
            IF (key_ref(e) /= key_new(e)) THEN
               CALL cp_abort(__LOCATION__, &
                             "The MLP cell list and CP2K's neighbour list disagree on which "// &
                             "pairs are within the model cutoff.")
            END IF
         END DO
      END IF

      DEALLOCATE (tier, centres, is_centre, key_ref, key_new, perm)
      CALL mlp_edge_buffer_release(buf)
      CALL timestop(handle)
   END SUBROUTINE check_edges_against_cell_list

! **************************************************************************************************
!> \brief Appends one key, growing the array geometrically.
!> \param a ...
!> \param n position to write, one past the old end
!> \param v ...
! **************************************************************************************************
   SUBROUTINE push_key(a, n, v)
      INTEGER(KIND=int_8), ALLOCATABLE, INTENT(INOUT)    :: a(:)
      INTEGER, INTENT(IN)                                :: n
      INTEGER(KIND=int_8), INTENT(IN)                    :: v

      INTEGER(KIND=int_8), ALLOCATABLE                   :: tmp(:)

      IF (n > SIZE(a)) THEN
         ALLOCATE (tmp(MAX(2*SIZE(a), 1024)))
         tmp(1:SIZE(a)) = a(:)
         CALL MOVE_ALLOC(tmp, a)
      END IF
      a(n) = v
   END SUBROUTINE push_key

! **************************************************************************************************
!> \brief Collapses one edge into a single sortable integer.
!> \param c centre, 0-based packed
!> \param j neighbour, 0-based packed
!> \param sh cell shift
!> \param n_use number of packed atoms
!> \return ...
! **************************************************************************************************
   PURE FUNCTION edge_key(c, j, sh, n_use) RESULT(key)
      INTEGER, INTENT(IN)                                :: c, j, sh(3), n_use
      INTEGER(KIND=int_8)                                :: key

      INTEGER, PARAMETER                                 :: sh_max = 16, sh_base = 2*sh_max + 1

      INTEGER                                            :: code

! A fixed, deliberately generous shift base, so both sides encode a shift
! the same way whatever range each of them happened to need.

      code = (sh(1) + sh_max) + sh_base*((sh(2) + sh_max) + sh_base*(sh(3) + sh_max))
      key = (INT(c, int_8)*INT(n_use, int_8) + INT(j, int_8))*INT(sh_base**3, int_8) &
            + INT(code, int_8)
   END FUNCTION edge_key

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
         ! One model per thread, built here rather than per force evaluation.
         ! Each clone re-reads the file, which is why this is done once: the
         ! evaluator holds its scratch inside the model, so threads cannot share
         ! one.
         ALLOCATE (sym_data%model(work%n_threads))
         ALLOCATE (sym_data%edge_buf(work%n_threads))
         DO e = 1, work%n_threads
            CALL symmetrix_model_load(sym_data%model(e), sym_pot%pot_file_name)
         END DO
         ! map atomic number Z -> symmetrix 0-based type index, once per model.
         ! The lookup covers both the model's elements and any element a covered
         ! CP2K kind carries (a mismatch is caught at evaluation, not by
         ! overflow); the kind coverage is fixed for the lifetime of the env.
         num_types = symmetrix_model_num_types(sym_data%model(1))
         ALLOCATE (model_z(num_types))
         CALL symmetrix_model_atomic_types(sym_data%model(1), model_z)
         ALLOCATE (sym_data%idx_of_z(0:MAX(MAXVAL(INT(model_z)), MAXVAL(work%kind_z))))
         sym_data%idx_of_z = -1
         DO e = 1, num_types
            sym_data%idx_of_z(INT(model_z(e))) = e - 1
         END DO
         DEALLOCATE (model_z)
      END IF

      IF (ASSOCIATED(sym_data%force)) THEN
         IF (SIZE(sym_data%force, 2) /= work%n_atoms_use) &
            DEALLOCATE (sym_data%force, sym_data%use_indices)
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
   SUBROUTINE run_symmetrix_model(sym_data, work, thr, my_thread, pot_total, virial_out)
      TYPE(symmetrix_data_type), POINTER                 :: sym_data
      TYPE(symmetrix_work_type), INTENT(IN)              :: work
      TYPE(symmetrix_thread_type), INTENT(INOUT)         :: thr
      INTEGER, INTENT(IN)                                :: my_thread
      REAL(kind=dp), INTENT(OUT)                         :: pot_total
      REAL(kind=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_out

      CHARACTER(LEN=*), PARAMETER :: routineN = 'run_symmetrix_model'

      INTEGER                                            :: c0, e, handle_stage, iat, iu, j0, &
                                                            n_edges, n_use, z
      INTEGER(C_INT32_T), ALLOCATABLE                    :: edge_centre(:), neigh_indices(:), &
                                                            neigh_types(:), node_types(:)
      REAL(kind=dp)                                      :: cell_v(3), disp(3), energy, &
                                                            hmat_ang(3, 3), virial(3, 3)
      REAL(kind=dp), ALLOCATABLE                         :: forces(:, :), pos_ang(:, :), r(:), &
                                                            xyz(:, :)

      pot_total = 0.0_dp
      virial_out(:, :) = 0.0_dp
      ALLOCATE (thr%force_node(3, MAX(thr%n_atoms_node, 1)))
      thr%force_node(:, :) = 0.0_dp
      IF (work%n_atoms_use == 0 .OR. thr%n_atoms_node == 0) RETURN

      n_use = thr%n_atoms_node
      n_edges = sym_data%edge_buf(my_thread + 1)%n_edges

      CALL timeset('symmetrix_csr', handle_stage)

      ! packed positions (Angstrom) and per-atom symmetrix type index; the
      ! Z -> type-index map was cached on sym_data at model load
      ! Node arrays are built in OWNED-FIRST order: node p holds canonical atom
      ! node_perm_inv(p+1). The evaluator needs this domain's owned nodes
      ! contiguous at the front (see symmetrix_c_eval_mp).
      ALLOCATE (pos_ang(3, n_use), node_types(n_use))
      DO iu = 1, n_use
         iat = work%use_indices(thr%node_perm_inv(iu) + 1)
         pos_ang(:, iu) = work%r_pbc(iat)%r(:)*angstrom
         z = work%kind_z(work%particle_set(iat)%atomic_kind%kind_number)
         IF (sym_data%idx_of_z(z) < 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "An atom element is not among the atomic numbers supported by the symmetrix model.")
         END IF
         node_types(iu) = INT(sym_data%idx_of_z(z), C_INT32_T)
      END DO
      hmat_ang(:, :) = work%cell%hmat*angstrom

      ! The neighbour search already emitted exactly the edges this domain
      ! needs, grouped by centre and in ascending centre order, and returned the
      ! per-centre counts. That IS the CSR layout, so the arrays below are one
      ! linear pass: no filtering of other domains' edges, no counting sweep and
      ! no prefix-sum scatter, all of which used to be needed because the edge
      ! list arrived as an unordered pile covering the whole system.
      ALLOCATE (neigh_indices(n_edges), neigh_types(n_edges), edge_centre(n_edges))
      ALLOCATE (xyz(3, n_edges), r(n_edges))
      DO e = 1, n_edges
         c0 = thr%node_perm(INT(sym_data%edge_buf(my_thread + 1)%pairs(1, e)) + 1)
         j0 = thr%node_perm(INT(sym_data%edge_buf(my_thread + 1)%pairs(2, e)) + 1)

         cell_v = MATMUL(hmat_ang, REAL(sym_data%edge_buf(my_thread + 1)%shifts(:, e), dp))
         disp(:) = pos_ang(:, j0 + 1) - pos_ang(:, c0 + 1) + cell_v

         neigh_indices(e) = INT(j0, C_INT32_T)
         edge_centre(e) = INT(c0, C_INT32_T)
         neigh_types(e) = node_types(j0 + 1)
         xyz(:, e) = disp(:)
         r(e) = NORM2(disp)
      END DO

      ALLOCATE (forces(3, n_use))
      CALL timestop(handle_stage)

      ! Atom decomposition: this domain evaluates the model only for the atoms it
      ! owns - packed nodes [0, n_atoms_local) - with the rest present as ghosts
      ! so the owned atoms' message-passing environments are complete.
      CALL timeset('symmetrix_eval', handle_stage)
      CALL symmetrix_eval(sym_data%model(my_thread + 1), node_types, thr%num_neigh, neigh_indices, neigh_types, &
                          edge_centre, xyz, r, work%use_virial, energy, forces, virial, &
                          n_local=thr%n_atoms_local, n_layer0=thr%n_atoms_layer0)
      CALL timestop(handle_stage)

      ! energy eV->Hartree, forces eV/Angstrom -> Hartree/Bohr, virial eV->Hartree.
      ! NO 1/num_pe: the energy already covers only this rank's owned atoms, and
      ! the forces/virial are already genuine partial derivatives of that partial
      ! energy. CP2K's cross-rank sums in fist_force.F assemble the exact totals.
      ! The virial lands as r (outer) f in CP2K's memory order; this is symmetric
      ! for a rotation-invariant total energy.
      pot_total = energy/evolt
      ! forces stay in this domain's own node order and are scattered onto the
      ! canonical packed numbering during the fold
      DO iu = 1, n_use
         thr%force_node(:, iu) = forces(:, iu)/evolt*angstrom
      END DO
      IF (work%use_virial) virial_out(:, :) = virial(:, :)/evolt

      DEALLOCATE (pos_ang, node_types, neigh_indices, neigh_types, edge_centre, xyz, r, forces)
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
!> \brief Release one slab's tier arrays.
!> \param thr ...
! **************************************************************************************************
   SUBROUTINE release_thread(thr)
      TYPE(symmetrix_thread_type), INTENT(INOUT)         :: thr

      IF (ALLOCATED(thr%num_neigh)) DEALLOCATE (thr%num_neigh)
      IF (ALLOCATED(thr%force_node)) DEALLOCATE (thr%force_node)

      IF (ALLOCATED(thr%node_perm)) DEALLOCATE (thr%node_perm)
      IF (ALLOCATED(thr%node_perm_inv)) DEALLOCATE (thr%node_perm_inv)
   END SUBROUTINE release_thread

! **************************************************************************************************
!> \brief Release the per-call work scratch.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE release_work(work)
      TYPE(symmetrix_work_type), INTENT(INOUT)           :: work

      CALL mlp_cell_list_release(work%cl)
      IF (ALLOCATED(work%thread_of)) DEALLOCATE (work%thread_of)
      IF (ALLOCATED(work%kind_z)) DEALLOCATE (work%kind_z)
      IF (ALLOCATED(work%global_to_packed)) DEALLOCATE (work%global_to_packed)
      IF (ALLOCATED(work%use_indices)) DEALLOCATE (work%use_indices)
      IF (ALLOCATED(work%packed_pos)) DEALLOCATE (work%packed_pos)
      NULLIFY (work%cell, work%r_pbc, work%particle_set, work%para_env)
   END SUBROUTINE release_work

END MODULE manybody_symmetrix

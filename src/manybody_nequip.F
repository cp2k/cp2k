!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      Implementation of NequIP and Allegro potentials - [gtocci] 2022
!>      Index mapping of atoms from .xyz to Allegro config.yaml file - [mbilichenko] 2024
!>      Refactoring and update to NequIP version >= v0.7.0 - [gtocci] 2026
!> \author Gabriele Tocci
! **************************************************************************************************
MODULE manybody_nequip

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type,&
                                              nequip_data_type,&
                                              pos_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: nequip_pot_type,&
                                              nequip_type,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_type
   USE particle_types,                  ONLY: particle_type
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: &
        torch_dict_create, torch_dict_get, torch_dict_insert, torch_dict_release, torch_dict_type, &
        torch_model_forward, torch_model_freeze, torch_model_load, torch_tensor_data_ptr, &
        torch_tensor_from_array, torch_tensor_release, torch_tensor_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: nequip_energy_store_force_virial, &
             nequip_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_nequip'

   TYPE, PRIVATE :: nequip_work_type
      INTEGER                       :: target_pot_type
      INTEGER                       :: n_atoms_use
      LOGICAL                       :: use_virial

      TYPE(cell_type), POINTER              :: cell => NULL()
      TYPE(pos_type), DIMENSION(:), POINTER :: r_pbc => NULL()
      TYPE(distribution_1d_type), POINTER   :: local_particles => NULL()
      TYPE(particle_type), POINTER          :: particle_set(:) => NULL()
      TYPE(mp_para_env_type), POINTER       :: para_env => NULL()

      LOGICAL, ALLOCATABLE                  :: use_atom(:)
      INTEGER(kind=int_8), ALLOCATABLE      :: local_edges(:, :)
      REAL(kind=dp), ALLOCATABLE            :: local_shifts(:, :)
      INTEGER(kind=int_8), ALLOCATABLE      :: final_edges(:, :)
      REAL(kind=dp), ALLOCATABLE            :: final_shifts(:, :)
      INTEGER, DIMENSION(:), ALLOCATABLE    :: kind_mapper
      LOGICAL, ALLOCATABLE                  :: sum_energy(:)
   END TYPE nequip_work_type

CONTAINS

! **************************************************************************************************
!> \brief ...
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
!> \param target_pot_type ...
!> \par History
!>      Implementation of the nequip potential - [gtocci] 2022
!>      Refactoring and unifying NequIP and Allegro - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE nequip_energy_store_force_virial(nonbonded, particle_set, local_particles, cell, &
                                               atomic_kind_set, potparm, r_last_update_pbc, &
                                               pot_total, fist_nonbond_env, para_env, use_virial, &
                                               target_pot_type)

      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      REAL(kind=dp)                                      :: pot_total
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: use_virial
      INTEGER, INTENT(IN)                                :: target_pot_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nequip_energy_store_force_virial'

      INTEGER                                            :: handle
      TYPE(nequip_data_type), POINTER                    :: neq_data
      TYPE(nequip_pot_type), POINTER                     :: neq_pot
      TYPE(nequip_work_type)                             :: nequip_work
      TYPE(torch_dict_type)                              :: outputs

      CALL timeset(routineN, handle)

      CALL nequip_work_create(nequip_work, atomic_kind_set, particle_set, local_particles, cell, &
                              r_last_update_pbc, para_env, potparm, target_pot_type, use_virial, &
                              neq_pot)

      IF (.NOT. ASSOCIATED(neq_pot)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL build_local_edges_shifts(nonbonded, potparm, nequip_work)

      CALL build_torch_edge_indexes(nequip_work)

      CALL setup_neq_data(fist_nonbond_env, neq_data, neq_pot, nequip_work)

      IF (nequip_work%target_pot_type == nequip_type) THEN
         CALL prepare_edges_shifts_nequip(nequip_work)
      ELSE
         CALL prepare_edges_shifts_allegro(nequip_work)
      END IF

      CALL run_torch_model(neq_data, neq_pot, nequip_work, outputs)

      CALL process_outputs(outputs, neq_data, neq_pot, pot_total, nequip_work)

      CALL torch_dict_release(outputs)
      CALL release_nequip_work(nequip_work)

      CALL timestop(handle)
   END SUBROUTINE nequip_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param cell ...
!> \param r_pbc ...
!> \param para_env ...
!> \param potparm ...
!> \param target_pot_type ...
!> \param use_virial ...
!> \param neq_pot ...
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE nequip_work_create(nequip_work, atomic_kind_set, particle_set, local_particles, cell, &
                                 r_pbc, para_env, potparm, target_pot_type, use_virial, neq_pot)
      TYPE(nequip_work_type), INTENT(OUT)                :: nequip_work
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_pbc
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, INTENT(IN)                                :: target_pot_type
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(nequip_pot_type), INTENT(OUT), POINTER        :: neq_pot

      nequip_work%target_pot_type = target_pot_type
      nequip_work%use_virial = use_virial
      nequip_work%cell => cell
      nequip_work%r_pbc => r_pbc
      nequip_work%particle_set => particle_set
      nequip_work%para_env => para_env
      nequip_work%local_particles => local_particles

      CALL get_potential_config(atomic_kind_set, potparm, target_pot_type, neq_pot)

      IF (.NOT. ASSOCIATED(neq_pot)) THEN
         RETURN
      END IF

      CALL build_kind_mapper(atomic_kind_set, neq_pot, nequip_work)

      CALL init_atom_masks(nequip_work)

   END SUBROUTINE nequip_work_create

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE release_nequip_work(nequip_work)
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      IF (ALLOCATED(nequip_work%final_edges)) DEALLOCATE (nequip_work%final_edges)
      IF (ALLOCATED(nequip_work%final_shifts)) DEALLOCATE (nequip_work%final_shifts)
      IF (ALLOCATED(nequip_work%local_edges)) DEALLOCATE (nequip_work%local_edges)
      IF (ALLOCATED(nequip_work%local_shifts)) DEALLOCATE (nequip_work%local_shifts)
      IF (ALLOCATED(nequip_work%use_atom)) DEALLOCATE (nequip_work%use_atom)
      IF (ALLOCATED(nequip_work%kind_mapper)) DEALLOCATE (nequip_work%kind_mapper)
      IF (ALLOCATED(nequip_work%sum_energy)) DEALLOCATE (nequip_work%sum_energy)
      NULLIFY (nequip_work%cell, nequip_work%r_pbc, nequip_work%particle_set, nequip_work%para_env, &
               nequip_work%local_particles)

   END SUBROUTINE release_nequip_work

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param potparm ...
!> \param nequip_work ...
!> \par History
!>      Build edges and cell shifts for the GNN - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE build_local_edges_shifts(nonbonded, potparm, nequip_work)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      INTEGER                                            :: atom_a, atom_b, i, idx_i, idx_j, iend, &
                                                            igrp, ikind, ilist, ipair, istart, &
                                                            jkind, n_max_edges, nedges, npairs
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      LOGICAL                                            :: do_nequip_allegro
      REAL(kind=dp)                                      :: cutsq_ij, drij, rij(3)
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, cvi
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      n_max_edges = 0
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         n_max_edges = n_max_edges + neighbor_kind_pair%npairs
      END DO

      ALLOCATE (nequip_work%local_edges(2, n_max_edges), nequip_work%local_shifts(3, n_max_edges))
      nedges = 0

      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE

         Kind_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)

            idx_i = nequip_work%kind_mapper(ikind)
            idx_j = nequip_work%kind_mapper(jkind)

            IF (idx_i < 1 .OR. idx_j < 1) THEN
               ! pair involving atom not defined in the NequIP model, skipping..
               CYCLE Kind_Loop
            END IF
            pot => potparm%pot(ikind, jkind)%pot
            do_nequip_allegro = .FALSE.
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == nequip_work%target_pot_type) THEN
                  do_nequip_allegro = .TRUE.
                  EXIT
               END IF
            END DO

            IF (.NOT. do_nequip_allegro) CYCLE Kind_Loop

            cutsq_ij = pot%set(i)%nequip%cutoff_matrix(idx_i, idx_j)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            cell_v = MATMUL(nequip_work%cell%hmat, cvi)

            DO ipair = istart, iend
               atom_a = neighbor_kind_pair%list(1, ipair)
               atom_b = neighbor_kind_pair%list(2, ipair)

               rij(:) = nequip_work%r_pbc(atom_b)%r(:) - nequip_work%r_pbc(atom_a)%r(:) + cell_v
               drij = DOT_PRODUCT(rij, rij)

               IF (drij <= cutsq_ij) THEN
                  nedges = nedges + 1
                  nequip_work%local_edges(:, nedges) = [atom_a, atom_b]
                  nequip_work%local_shifts(:, nedges) = cvi
               END IF
            END DO
         END DO Kind_Loop
      END DO

      IF (nedges < n_max_edges) THEN
         BLOCK
            INTEGER(kind=int_8), ALLOCATABLE :: tmp_idx(:, :)
            REAL(kind=dp), ALLOCATABLE :: tmp_sft(:, :)

            ALLOCATE (tmp_idx(2, nedges), tmp_sft(3, nedges))

            tmp_idx(:, :) = nequip_work%local_edges(:, 1:nedges)
            tmp_sft(:, :) = nequip_work%local_shifts(:, 1:nedges)

            CALL MOVE_ALLOC(tmp_idx, nequip_work%local_edges)
            CALL MOVE_ALLOC(tmp_sft, nequip_work%local_shifts)
         END BLOCK
      END IF

   END SUBROUTINE build_local_edges_shifts

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param target_pot_type ...
!> \param neq_pot ...
!> \par History
!>      Get the NequIP or Allegro potential - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE get_potential_config(atomic_kind_set, potparm, target_pot_type, neq_pot)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, INTENT(IN)                                :: target_pot_type
      TYPE(nequip_pot_type), INTENT(OUT), POINTER        :: neq_pot

      INTEGER                                            :: i, ikind, jkind
      TYPE(pair_potential_single_type), POINTER          :: pot

      NULLIFY (neq_pot)
      OuterLoop: DO ikind = 1, SIZE(atomic_kind_set)
         DO jkind = ikind, SIZE(atomic_kind_set)
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == target_pot_type) THEN
                  neq_pot => pot%set(i)%nequip
                  EXIT OuterLoop
               END IF
            END DO
         END DO
      END DO OuterLoop
   END SUBROUTINE get_potential_config

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \par History
!>      Inits masks for torch evaluation (use_atom) and MPI summation (sum_energy) - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE init_atom_masks(nequip_work)
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      INTEGER                                            :: iat, ikind, ilocal, n_atoms, n_local

      IF (.NOT. ALLOCATED(nequip_work%kind_mapper)) THEN
         CPABORT("kind_mapper not initialized before init_atom_masks")
      END IF

      n_atoms = SIZE(nequip_work%particle_set)

      IF (ALLOCATED(nequip_work%use_atom)) DEALLOCATE (nequip_work%use_atom)
      ALLOCATE (nequip_work%use_atom(n_atoms))
      nequip_work%use_atom = .FALSE.

      DO iat = 1, n_atoms
         ikind = nequip_work%particle_set(iat)%atomic_kind%kind_number
         IF (nequip_work%kind_mapper(ikind) > 0) THEN
            nequip_work%use_atom(iat) = .TRUE.
         END IF
      END DO
      nequip_work%n_atoms_use = COUNT(nequip_work%use_atom)

      IF (ALLOCATED(nequip_work%sum_energy)) DEALLOCATE (nequip_work%sum_energy)
      ALLOCATE (nequip_work%sum_energy(n_atoms))
      nequip_work%sum_energy = .FALSE.

      IF (ASSOCIATED(nequip_work%local_particles)) THEN
         DO ikind = 1, SIZE(nequip_work%local_particles%n_el)
            IF (nequip_work%kind_mapper(ikind) > 0) THEN
               n_local = nequip_work%local_particles%n_el(ikind)
               DO ilocal = 1, n_local
                  iat = nequip_work%local_particles%list(ikind)%array(ilocal)
                  nequip_work%sum_energy(iat) = .TRUE.
               END DO
            END IF
         END DO
      ELSE
         nequip_work%sum_energy(:) = nequip_work%use_atom(:)
      END IF

   END SUBROUTINE init_atom_masks

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param neq_pot ...
!> \param nequip_work ...
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE build_kind_mapper(atomic_kind_set, neq_pot, nequip_work)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(nequip_pot_type), POINTER                     :: neq_pot
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      CHARACTER(LEN=100)                                 :: model_sym
      CHARACTER(LEN=default_string_length)               :: kind_sym
      INTEGER                                            :: i, ikind, n_kinds

      n_kinds = SIZE(atomic_kind_set)

      IF (ALLOCATED(nequip_work%kind_mapper)) DEALLOCATE (nequip_work%kind_mapper)
      ALLOCATE (nequip_work%kind_mapper(n_kinds))
      nequip_work%kind_mapper = -1

      DO ikind = 1, n_kinds
         kind_sym = atomic_kind_set(ikind)%element_symbol
         CALL uppercase(kind_sym)

         DO i = 1, neq_pot%num_types
            model_sym = neq_pot%type_names_torch(i)
            CALL uppercase(model_sym)
            IF (TRIM(kind_sym) == TRIM(model_sym)) THEN
               nequip_work%kind_mapper(ikind) = i
               EXIT
            END IF
         END DO
      END DO
   END SUBROUTINE build_kind_mapper

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param neq_data ...
!> \param pot ...
!> \param nequip_work ...
!> \par History
!>      load the NequIP/Allegro model, initialize forces, positions  - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE setup_neq_data(fist_nonbond_env, neq_data, pot, nequip_work)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(nequip_data_type), POINTER                    :: neq_data
      TYPE(nequip_pot_type), POINTER                     :: pot
      TYPE(nequip_work_type), INTENT(IN)                 :: nequip_work

      INTEGER                                            :: iat, iat_use, n_atoms

      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=neq_data)

      IF (.NOT. ASSOCIATED(neq_data)) THEN
         ALLOCATE (neq_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, nequip_data=neq_data)
         NULLIFY (neq_data%use_indices, neq_data%force)

         CALL torch_model_load(neq_data%model, pot%pot_file_name)
         CALL torch_model_freeze(neq_data%model)
      END IF

      IF (ASSOCIATED(neq_data%force)) THEN
         IF (SIZE(neq_data%force, 2) /= nequip_work%n_atoms_use) &
            DEALLOCATE (neq_data%force, neq_data%use_indices)
      END IF

      IF (.NOT. ASSOCIATED(neq_data%force)) THEN
         ALLOCATE (neq_data%force(3, nequip_work%n_atoms_use))
         ALLOCATE (neq_data%use_indices(nequip_work%n_atoms_use))
      END IF

      n_atoms = SIZE(nequip_work%use_atom)
      iat_use = 0
      DO iat = 1, n_atoms
         IF (nequip_work%use_atom(iat)) THEN
            iat_use = iat_use + 1
            neq_data%use_indices(iat_use) = iat
         END IF
      END DO
   END SUBROUTINE setup_neq_data

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \par History
!>      Prepare edges and cell shifts for NequIP  - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE prepare_edges_shifts_nequip(nequip_work)
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      INTEGER                                            :: ipair, nedges, nedges_tot
      INTEGER(kind=int_8), ALLOCATABLE                   :: temp_edge_index(:, :)
      INTEGER, ALLOCATABLE                               :: displ(:), displ_cell(:), edge_count(:), &
                                                            edge_count_cell(:)

      nedges = SIZE(nequip_work%local_edges, 2)

      ALLOCATE (edge_count(nequip_work%para_env%num_pe), edge_count_cell(nequip_work%para_env%num_pe))
      ALLOCATE (displ_cell(nequip_work%para_env%num_pe), displ(nequip_work%para_env%num_pe))

      CALL nequip_work%para_env%allgather(nedges, edge_count)
      nedges_tot = SUM(edge_count)

      ALLOCATE (temp_edge_index(2, nedges_tot))
      ALLOCATE (nequip_work%final_shifts(3, nedges_tot))

      edge_count_cell(:) = edge_count*3
      edge_count = edge_count*2
      displ(1) = 0
      displ_cell(1) = 0
      DO ipair = 2, nequip_work%para_env%num_pe
         displ(ipair) = displ(ipair - 1) + edge_count(ipair - 1)
         displ_cell(ipair) = displ_cell(ipair - 1) + edge_count_cell(ipair - 1)
      END DO

      CALL nequip_work%para_env%allgatherv(nequip_work%local_shifts, nequip_work%final_shifts, edge_count_cell, displ_cell)
      CALL nequip_work%para_env%allgatherv(nequip_work%local_edges, temp_edge_index, edge_count, displ)

      ALLOCATE (nequip_work%final_edges(nedges_tot, 2))
      nequip_work%final_edges(:, :) = TRANSPOSE(temp_edge_index)

      DEALLOCATE (edge_count, edge_count_cell, displ, displ_cell, temp_edge_index)

   END SUBROUTINE prepare_edges_shifts_nequip

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \par History
!>      Prepare edges and cell shifts for Allegro  - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE prepare_edges_shifts_allegro(nequip_work)
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      ALLOCATE (nequip_work%final_shifts, SOURCE=nequip_work%local_shifts)
      ALLOCATE (nequip_work%final_edges(SIZE(nequip_work%local_edges, 2), 2))
      nequip_work%final_edges(:, :) = TRANSPOSE(nequip_work%local_edges)
   END SUBROUTINE prepare_edges_shifts_allegro

! **************************************************************************************************
!> \brief ...
!> \param nequip_work ...
!> \par History
!>      Build edges from cp2k global neigh lists to local/packed ones for torch - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE build_torch_edge_indexes(nequip_work)
      TYPE(nequip_work_type), INTENT(INOUT)              :: nequip_work

      INTEGER                                            :: atom_a, atom_b, i, iat, iat_use, n_atoms
      INTEGER, ALLOCATABLE                               :: global_to_packed(:)

      n_atoms = SIZE(nequip_work%particle_set)

      ! for allegro ensure ghost atoms are included in the evaluation
      IF (nequip_work%target_pot_type /= nequip_type) THEN
         ! label atoms in the local edges
         DO i = 1, SIZE(nequip_work%local_edges, 2)
            atom_a = INT(nequip_work%local_edges(1, i))
            atom_b = INT(nequip_work%local_edges(2, i))
            nequip_work%use_atom(atom_a) = .TRUE.
            nequip_work%use_atom(atom_b) = .TRUE.
         END DO
         nequip_work%n_atoms_use = COUNT(nequip_work%use_atom)
      END IF

      ! mapping from global CP2K index to packed/local Torch index
      ALLOCATE (global_to_packed(n_atoms))
      global_to_packed = 0
      iat_use = 0
      DO iat = 1, n_atoms
         IF (nequip_work%use_atom(iat)) THEN
            iat_use = iat_use + 1
            global_to_packed(iat) = iat_use
         END IF
      END DO

      ! remap local_edges to use 0-based dense indices for torch
      DO i = 1, SIZE(nequip_work%local_edges, 2)
         atom_a = INT(nequip_work%local_edges(1, i))
         atom_b = INT(nequip_work%local_edges(2, i))

         nequip_work%local_edges(1, i) = INT(global_to_packed(atom_a) - 1, kind=int_8)
         nequip_work%local_edges(2, i) = INT(global_to_packed(atom_b) - 1, kind=int_8)
      END DO

      DEALLOCATE (global_to_packed)

   END SUBROUTINE build_torch_edge_indexes

! **************************************************************************************************
!> \brief ...
!> \param neq_data ...
!> \param pot ...
!> \param nequip_work ...
!> \param outputs ...
!> \par History
!>      Run forward pass using torch api  - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE run_torch_model(neq_data, pot, nequip_work, outputs)
      TYPE(nequip_data_type), POINTER                    :: neq_data
      TYPE(nequip_pot_type), POINTER                     :: pot
      TYPE(nequip_work_type), INTENT(IN)                 :: nequip_work
      TYPE(torch_dict_type), INTENT(OUT)                 :: outputs

      INTEGER                                            :: iat, iat_use, ikind
      INTEGER(kind=int_8), ALLOCATABLE                   :: atom_types(:)
      REAL(kind=dp), ALLOCATABLE                         :: lattice(:, :), pos(:, :)
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_tensor_type)                            :: cell_t, idx_t, pos_t, shift_t, types_t

      ALLOCATE (lattice(3, 3))
      lattice(:, :) = nequip_work%cell%hmat/pot%unit_length_val

      ALLOCATE (pos(3, nequip_work%n_atoms_use), atom_types(nequip_work%n_atoms_use))
      iat_use = 0
      DO iat = 1, SIZE(nequip_work%particle_set)
         IF (.NOT. nequip_work%use_atom(iat)) CYCLE
         iat_use = iat_use + 1

         ikind = nequip_work%particle_set(iat)%atomic_kind%kind_number
         IF (nequip_work%kind_mapper(ikind) < 1) THEN
            CALL cp_abort(__LOCATION__, "Atom symbol not found in NequIP model!")
         END IF

         ! Convert 1-based Fortran index to 0-based PyTorch index
         atom_types(iat_use) = nequip_work%kind_mapper(ikind) - 1
         pos(:, iat_use) = nequip_work%r_pbc(iat)%r(:)/pot%unit_length_val
      END DO

      CALL torch_dict_create(inputs)

      CALL torch_tensor_from_array(pos_t, pos)
      CALL torch_tensor_from_array(shift_t, nequip_work%final_shifts)
      CALL torch_tensor_from_array(cell_t, lattice)

      CALL torch_dict_insert(inputs, "pos", pos_t)
      CALL torch_dict_insert(inputs, "edge_cell_shift", shift_t)
      CALL torch_dict_insert(inputs, "cell", cell_t)
      CALL torch_tensor_release(pos_t)
      CALL torch_tensor_release(shift_t)
      CALL torch_tensor_release(cell_t)

      CALL torch_tensor_from_array(idx_t, nequip_work%final_edges)
      CALL torch_dict_insert(inputs, "edge_index", idx_t)
      CALL torch_tensor_release(idx_t)

      CALL torch_tensor_from_array(types_t, atom_types)
      CALL torch_dict_insert(inputs, "atom_types", types_t)
      CALL torch_tensor_release(types_t)

      CALL torch_dict_create(outputs)
      CALL torch_model_forward(neq_data%model, inputs, outputs)

      CALL torch_dict_release(inputs)

      IF (ALLOCATED(pos)) DEALLOCATE (pos)
      IF (ALLOCATED(lattice)) DEALLOCATE (lattice)
      IF (ALLOCATED(atom_types)) DEALLOCATE (atom_types)

   END SUBROUTINE run_torch_model

! **************************************************************************************************
!> \brief ...
!> \param outputs ...
!> \param neq_data ...
!> \param pot ...
!> \param pot_total ...
!> \param nequip_work ...
!> \par History
!>      Collect potential, forces, virial  - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE process_outputs(outputs, neq_data, pot, pot_total, nequip_work)
      TYPE(torch_dict_type), INTENT(IN)                  :: outputs
      TYPE(nequip_data_type), POINTER                    :: neq_data
      TYPE(nequip_pot_type), POINTER                     :: pot
      REAL(kind=dp), INTENT(OUT)                         :: pot_total
      TYPE(nequip_work_type), INTENT(IN)                 :: nequip_work

      INTEGER                                            :: iat, iat_use
      REAL(kind=dp), POINTER                             :: e_ptr(:, :), f_ptr(:, :), v_ptr(:, :, :)
      TYPE(torch_tensor_type)                            :: t_energy, t_forces, t_virial

      NULLIFY (f_ptr, e_ptr, v_ptr)

      CALL torch_dict_get(outputs, "forces", t_forces)
      CALL torch_tensor_data_ptr(t_forces, f_ptr)

      neq_data%force = f_ptr*pot%unit_forces_val
      CALL torch_tensor_release(t_forces)
      CALL torch_dict_get(outputs, "atomic_energy", t_energy)
      CALL torch_tensor_data_ptr(t_energy, e_ptr)

      pot_total = 0.0_dp
      DO iat_use = 1, SIZE(neq_data%use_indices)
         iat = neq_data%use_indices(iat_use)
         ! Only apply the local mask for Allegro models
         IF (nequip_work%target_pot_type /= nequip_type) THEN
            IF (.NOT. nequip_work%sum_energy(iat)) CYCLE
         END IF

         pot_total = pot_total + e_ptr(1, iat_use)
      END DO
      CALL torch_tensor_release(t_energy)
      pot_total = pot_total*pot%unit_energy_val

      IF (nequip_work%target_pot_type == nequip_type) THEN
         neq_data%force = neq_data%force/REAL(nequip_work%para_env%num_pe, dp)
         pot_total = pot_total/REAL(nequip_work%para_env%num_pe, dp)
      END IF

      IF (nequip_work%use_virial) THEN
         CALL torch_dict_get(outputs, "virial", t_virial)
         CALL torch_tensor_data_ptr(t_virial, v_ptr)

         neq_data%virial(:, :) = RESHAPE(v_ptr, [3, 3])*pot%unit_energy_val
         CALL torch_tensor_release(t_virial)
         IF (nequip_work%target_pot_type == nequip_type) THEN
            neq_data%virial = neq_data%virial/REAL(nequip_work%para_env%num_pe, dp)
         END IF
      END IF

   END SUBROUTINE process_outputs

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
!> \par History
!>      Sum forces, virial to nonbond - [gtocci] 2026
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE nequip_add_force_virial(fist_nonbond_env, f_nonbond, pv_nonbond, use_virial)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial

      INTEGER                                            :: iat, iat_use
      TYPE(nequip_data_type), POINTER                    :: neq_data

      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=neq_data)

      IF (use_virial) THEN
         pv_nonbond = pv_nonbond + neq_data%virial
      END IF

      DO iat_use = 1, SIZE(neq_data%use_indices)
         iat = neq_data%use_indices(iat_use)
         f_nonbond(1:3, iat) = f_nonbond(1:3, iat) + neq_data%force(1:3, iat_use)
      END DO

   END SUBROUTINE nequip_add_force_virial

END MODULE manybody_nequip

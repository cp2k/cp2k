!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      allegro implementation
!> \author Gabriele Tocci
! **************************************************************************************************
MODULE manybody_allegro

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: allegro_data_type,&
                                              fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type,&
                                              pos_type
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              sp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: allegro_pot_type,&
                                              allegro_type,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_type
   USE particle_types,                  ONLY: particle_type
   USE torch_api,                       ONLY: &
        torch_dict_create, torch_dict_get, torch_dict_insert, torch_dict_release, torch_dict_type, &
        torch_model_forward, torch_model_freeze, torch_model_load, torch_tensor_data_ptr, &
        torch_tensor_from_array, torch_tensor_release, torch_tensor_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: setup_allegro_arrays, destroy_allegro_arrays, &
             allegro_energy_store_force_virial, allegro_add_force_virial
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_allegro'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param potparm ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param unique_list_a ...
!> \param cell ...
!> \par History
!>      Implementation of the allegro potential - [gtocci] 2023
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE setup_allegro_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, &
                                   unique_list_a, cell)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a, unique_list_a
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_allegro_arrays'

      INTEGER                                            :: handle, i, iend, igrp, ikind, ilist, &
                                                            ipair, istart, jkind, nkinds, nlocal, &
                                                            npairs, npairs_tot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: temp_unique_list_a, work_list, work_list2
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rwork_list
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      CPASSERT(.NOT. ASSOCIATED(glob_loc_list))
      CPASSERT(.NOT. ASSOCIATED(glob_loc_list_a))
      CPASSERT(.NOT. ASSOCIATED(unique_list_a))
      CPASSERT(.NOT. ASSOCIATED(glob_cell_v))
      CALL timeset(routineN, handle)
      npairs_tot = 0
      nkinds = SIZE(potparm%pot, 1)
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop1: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            pot => potparm%pot(ikind, jkind)%pot
            npairs = iend - istart + 1
            IF (pot%no_mb) CYCLE
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) == allegro_type) npairs_tot = npairs_tot + npairs
            END DO
         END DO Kind_Group_Loop1
      END DO
      ALLOCATE (work_list(npairs_tot))
      ALLOCATE (work_list2(npairs_tot))
      ALLOCATE (glob_loc_list(2, npairs_tot))
      ALLOCATE (glob_cell_v(3, npairs_tot))
      ! Fill arrays with data
      npairs_tot = 0
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop2: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            npairs = iend - istart + 1
            IF (pot%no_mb) CYCLE
            cell_v = MATMUL(cell%hmat, cvi)
            DO i = 1, SIZE(pot%type)
               ! ALLEGRO
               IF (pot%type(i) == allegro_type) THEN
                  DO ipair = 1, npairs
                     glob_loc_list(:, npairs_tot + ipair) = list(:, istart - 1 + ipair)
                     glob_cell_v(1:3, npairs_tot + ipair) = cell_v(1:3)
                  END DO
                  npairs_tot = npairs_tot + npairs
               END IF
            END DO
         END DO Kind_Group_Loop2
      END DO
      ! Order the arrays w.r.t. the first index of glob_loc_list
      CALL sort(glob_loc_list(1, :), npairs_tot, work_list)
      DO ipair = 1, npairs_tot
         work_list2(ipair) = glob_loc_list(2, work_list(ipair))
      END DO
      glob_loc_list(2, :) = work_list2
      DEALLOCATE (work_list2)
      ALLOCATE (rwork_list(3, npairs_tot))
      DO ipair = 1, npairs_tot
         rwork_list(:, ipair) = glob_cell_v(:, work_list(ipair))
      END DO
      glob_cell_v = rwork_list
      DEALLOCATE (rwork_list)
      DEALLOCATE (work_list)
      ALLOCATE (glob_loc_list_a(npairs_tot))
      glob_loc_list_a = glob_loc_list(1, :)
      ALLOCATE (temp_unique_list_a(npairs_tot))
      nlocal = 1
      temp_unique_list_a(1) = glob_loc_list_a(1)
      DO ipair = 2, npairs_tot
         IF (glob_loc_list_a(ipair - 1) /= glob_loc_list_a(ipair)) THEN
            nlocal = nlocal + 1
            temp_unique_list_a(nlocal) = glob_loc_list_a(ipair)
         END IF
      END DO
      ALLOCATE (unique_list_a(nlocal))
      unique_list_a(:) = temp_unique_list_a(:nlocal)
      DEALLOCATE (temp_unique_list_a)
      CALL timestop(handle)
   END SUBROUTINE setup_allegro_arrays

! **************************************************************************************************
!> \brief ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param unique_list_a ...
!> \par History
!>      Implementation of the allegro potential - [gtocci] 2023
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE destroy_allegro_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a, unique_list_a)
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a, unique_list_a

      IF (ASSOCIATED(glob_loc_list)) THEN
         DEALLOCATE (glob_loc_list)
      END IF
      IF (ASSOCIATED(glob_loc_list_a)) THEN
         DEALLOCATE (glob_loc_list_a)
      END IF
      IF (ASSOCIATED(glob_cell_v)) THEN
         DEALLOCATE (glob_cell_v)
      END IF
      IF (ASSOCIATED(unique_list_a)) THEN
         DEALLOCATE (unique_list_a)
      END IF

   END SUBROUTINE destroy_allegro_arrays

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param allegro ...
!> \param glob_loc_list_a ...
!> \param r_last_update_pbc ...
!> \param pot_allegro ...
!> \param fist_nonbond_env ...
!> \param unique_list_a ...
!> \param para_env ...
!> \param use_virial ...
!> \par History
!>      Implementation of the allegro potential - [gtocci] 2023
!>      Index mapping of atoms from .xyz to Allegro config.yaml file - [mbilichenko] 2024
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE allegro_energy_store_force_virial(nonbonded, particle_set, cell, atomic_kind_set, &
                                                potparm, allegro, glob_loc_list_a, r_last_update_pbc, &
                                                pot_allegro, fist_nonbond_env, unique_list_a, para_env, use_virial)

      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(allegro_pot_type), POINTER                    :: allegro
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      REAL(kind=dp)                                      :: pot_allegro
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      INTEGER, DIMENSION(:), POINTER                     :: unique_list_a
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: use_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allegro_energy_store_force_virial'

      INTEGER :: atom_a, atom_b, atom_idx, handle, i, iat, iat_use, iend, ifirst, igrp, ikind, &
         ilast, ilist, ipair, istart, iunique, jkind, junique, mpair, n_atoms, n_atoms_use, &
         nedges, nloc_size, npairs, nunique
      INTEGER(kind=int_8), ALLOCATABLE                   :: atom_types(:), temp_atom_types(:)
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:, :)  :: edge_index, t_edge_index, temp_edge_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: work_list
      INTEGER, DIMENSION(:, :), POINTER                  :: list, sort_list
      LOGICAL, ALLOCATABLE                               :: use_atom(:)
      REAL(kind=dp)                                      :: drij, rab2_max, rij(3)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: edge_cell_shifts, lattice, &
                                                            new_edge_cell_shifts, pos
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, cvi
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: atomic_energy, forces, virial
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: virial3d
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:, :)        :: lattice_sp, new_edge_cell_shifts_sp, &
                                                            pos_sp
      REAL(kind=sp), DIMENSION(:, :), POINTER            :: atomic_energy_sp, forces_sp
      TYPE(allegro_data_type), POINTER                   :: allegro_data
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot
      TYPE(torch_dict_type)                              :: inputs, outputs
      TYPE(torch_tensor_type) :: atom_types_tensor, atomic_energy_tensor, forces_tensor, &
         lattice_tensor, new_edge_cell_shifts_tensor, pos_tensor, t_edge_index_tensor, &
         virial_tensor

      CALL timeset(routineN, handle)

      NULLIFY (atomic_energy, forces, atomic_energy_sp, forces_sp, virial3d, virial)
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      use_atom = .FALSE.

      DO ikind = 1, SIZE(atomic_kind_set)
         DO jkind = 1, SIZE(atomic_kind_set)
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= allegro_type) CYCLE
               DO iat = 1, n_atoms
                  IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
                      particle_set(iat)%atomic_kind%kind_number == jkind) use_atom(iat) = .TRUE.
               END DO ! iat
            END DO ! i
         END DO ! jkind
      END DO ! ikind
      n_atoms_use = COUNT(use_atom)

      ! get allegro_data to save force, virial info and to load model
      CALL fist_nonbond_env_get(fist_nonbond_env, allegro_data=allegro_data)
      IF (.NOT. ASSOCIATED(allegro_data)) THEN
         ALLOCATE (allegro_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, allegro_data=allegro_data)
         NULLIFY (allegro_data%use_indices, allegro_data%force)
         CALL torch_model_load(allegro_data%model, pot%set(1)%allegro%allegro_file_name)
         CALL torch_model_freeze(allegro_data%model)
      END IF
      IF (ASSOCIATED(allegro_data%force)) THEN
         IF (SIZE(allegro_data%force, 2) /= n_atoms_use) THEN
            DEALLOCATE (allegro_data%force, allegro_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(allegro_data%force)) THEN
         ALLOCATE (allegro_data%force(3, n_atoms_use))
         ALLOCATE (allegro_data%use_indices(n_atoms_use))
      END IF

      iat_use = 0
      DO iat = 1, n_atoms_use
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            allegro_data%use_indices(iat_use) = iat
         END IF
      END DO

      nedges = 0

      ALLOCATE (edge_index(2, SIZE(glob_loc_list_a)))
      ALLOCATE (edge_cell_shifts(3, SIZE(glob_loc_list_a)))
      ALLOCATE (temp_atom_types(SIZE(glob_loc_list_a)))

      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop_Allegro: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= allegro_type) CYCLE
               rab2_max = pot%set(i)%allegro%rcutsq
               cell_v = MATMUL(cell%hmat, cvi)
               pot => potparm%pot(ikind, jkind)%pot
               allegro => pot%set(i)%allegro
               npairs = iend - istart + 1
               IF (npairs /= 0) THEN
                  ALLOCATE (sort_list(2, npairs), work_list(npairs))
                  sort_list = list(:, istart:iend)
                  ! Sort the list of neighbors, this increases the efficiency for single
                  ! potential contributions
                  CALL sort(sort_list(1, :), npairs, work_list)
                  DO ipair = 1, npairs
                     work_list(ipair) = sort_list(2, work_list(ipair))
                  END DO
                  sort_list(2, :) = work_list
                  ! find number of unique elements of array index 1
                  nunique = 1
                  DO ipair = 1, npairs - 1
                     IF (sort_list(1, ipair + 1) /= sort_list(1, ipair)) nunique = nunique + 1
                  END DO
                  ipair = 1
                  junique = sort_list(1, ipair)
                  ifirst = 1
                  DO iunique = 1, nunique
                     atom_a = junique
                     IF (glob_loc_list_a(ifirst) > atom_a) CYCLE
                     DO mpair = ifirst, SIZE(glob_loc_list_a)
                        IF (glob_loc_list_a(mpair) == atom_a) EXIT
                     END DO
                     ifirst = mpair
                     DO mpair = ifirst, SIZE(glob_loc_list_a)
                        IF (glob_loc_list_a(mpair) /= atom_a) EXIT
                     END DO
                     ilast = mpair - 1
                     nloc_size = 0
                     IF (ifirst /= 0) nloc_size = ilast - ifirst + 1
                     DO WHILE (ipair <= npairs)
                        IF (sort_list(1, ipair) /= junique) EXIT
                        atom_b = sort_list(2, ipair)
                        rij(:) = r_last_update_pbc(atom_b)%r(:) - r_last_update_pbc(atom_a)%r(:) + cell_v
                        drij = DOT_PRODUCT(rij, rij)
                        ipair = ipair + 1
                        IF (drij <= rab2_max) THEN
                           nedges = nedges + 1
                           edge_index(:, nedges) = [atom_a - 1, atom_b - 1]
                           edge_cell_shifts(:, nedges) = cvi
                        END IF
                     END DO
                     ifirst = ilast + 1
                     IF (ipair <= npairs) junique = sort_list(1, ipair)
                  END DO
                  DEALLOCATE (sort_list, work_list)
               END IF
            END DO
         END DO Kind_Group_Loop_Allegro
      END DO

      allegro => pot%set(1)%allegro

      ALLOCATE (temp_edge_index(2, nedges))
      temp_edge_index(:, :) = edge_index(:, :nedges)
      ALLOCATE (new_edge_cell_shifts(3, nedges))
      new_edge_cell_shifts(:, :) = edge_cell_shifts(:, :nedges)
      DEALLOCATE (edge_cell_shifts)

      ALLOCATE (t_edge_index(nedges, 2))

      t_edge_index(:, :) = TRANSPOSE(temp_edge_index)
      DEALLOCATE (temp_edge_index, edge_index)
      ALLOCATE (lattice(3, 3), lattice_sp(3, 3))
      lattice(:, :) = cell%hmat/pot%set(1)%allegro%unit_cell_val
      lattice_sp(:, :) = REAL(lattice, kind=sp)
      iat_use = 0
      ALLOCATE (pos(3, n_atoms_use), atom_types(n_atoms_use))
      DO iat = 1, n_atoms_use
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         ! Find index of the element based on its position in config.yaml file to have correct mapping
         DO i = 1, SIZE(allegro%type_names_torch)
            IF (particle_set(iat)%atomic_kind%element_symbol == allegro%type_names_torch(i)) THEN
               atom_idx = i - 1
            END IF
         END DO
         atom_types(iat_use) = atom_idx
         pos(:, iat) = r_last_update_pbc(iat)%r(:)/allegro%unit_coords_val
      END DO

      CALL torch_dict_create(inputs)

      IF (allegro%do_allegro_sp) THEN
         ALLOCATE (new_edge_cell_shifts_sp(3, nedges), pos_sp(3, n_atoms_use))
         new_edge_cell_shifts_sp(:, :) = REAL(new_edge_cell_shifts(:, :), kind=sp)
         pos_sp(:, :) = REAL(pos(:, :), kind=sp)
         DEALLOCATE (pos, new_edge_cell_shifts)
         CALL torch_tensor_from_array(pos_tensor, pos_sp)
         CALL torch_tensor_from_array(new_edge_cell_shifts_tensor, new_edge_cell_shifts_sp)
         CALL torch_tensor_from_array(lattice_tensor, lattice_sp)
      ELSE
         CALL torch_tensor_from_array(pos_tensor, pos)
         CALL torch_tensor_from_array(new_edge_cell_shifts_tensor, new_edge_cell_shifts)
         CALL torch_tensor_from_array(lattice_tensor, lattice)
      END IF

      CALL torch_dict_insert(inputs, "pos", pos_tensor)
      CALL torch_dict_insert(inputs, "edge_cell_shift", new_edge_cell_shifts_tensor)
      CALL torch_dict_insert(inputs, "cell", lattice_tensor)
      CALL torch_tensor_release(pos_tensor)
      CALL torch_tensor_release(new_edge_cell_shifts_tensor)
      CALL torch_tensor_release(lattice_tensor)

      CALL torch_tensor_from_array(t_edge_index_tensor, t_edge_index)
      CALL torch_dict_insert(inputs, "edge_index", t_edge_index_tensor)
      CALL torch_tensor_release(t_edge_index_tensor)

      CALL torch_tensor_from_array(atom_types_tensor, atom_types)
      CALL torch_dict_insert(inputs, "atom_types", atom_types_tensor)
      CALL torch_tensor_release(atom_types_tensor)

      CALL torch_dict_create(outputs)
      CALL torch_model_forward(allegro_data%model, inputs, outputs)
      pot_allegro = 0.0_dp

      CALL torch_dict_get(outputs, "atomic_energy", atomic_energy_tensor)
      CALL torch_dict_get(outputs, "forces", forces_tensor)
      IF (allegro%do_allegro_sp) THEN
         CALL torch_tensor_data_ptr(atomic_energy_tensor, atomic_energy_sp)
         CALL torch_tensor_data_ptr(forces_tensor, forces_sp)
         allegro_data%force(:, :) = REAL(forces_sp(:, :), kind=dp)*allegro%unit_forces_val
         DO iat_use = 1, SIZE(unique_list_a)
            i = unique_list_a(iat_use)
            pot_allegro = pot_allegro + REAL(atomic_energy_sp(1, i), kind=dp)*allegro%unit_energy_val
         END DO
         DEALLOCATE (new_edge_cell_shifts_sp, pos_sp)
      ELSE
         CALL torch_tensor_data_ptr(atomic_energy_tensor, atomic_energy)
         CALL torch_tensor_data_ptr(forces_tensor, forces)

         allegro_data%force(:, :) = forces(:, :)*allegro%unit_forces_val
         DO iat_use = 1, SIZE(unique_list_a)
            i = unique_list_a(iat_use)
            pot_allegro = pot_allegro + atomic_energy(1, i)*allegro%unit_energy_val
         END DO
         DEALLOCATE (pos, new_edge_cell_shifts)
      END IF
      CALL torch_tensor_release(atomic_energy_tensor)
      CALL torch_tensor_release(forces_tensor)

      IF (use_virial) THEN
         CALL torch_dict_get(outputs, "virial", virial_tensor)
         CALL torch_tensor_data_ptr(virial_tensor, virial3d)
         allegro_data%virial(:, :) = RESHAPE(virial3d, (/3, 3/))*allegro%unit_energy_val
         CALL torch_tensor_release(virial_tensor)
      END IF

      CALL torch_dict_release(inputs)
      CALL torch_dict_release(outputs)

      DEALLOCATE (t_edge_index, atom_types)

      IF (use_virial) allegro_data%virial(:, :) = allegro_data%virial/REAL(para_env%num_pe, dp)
      CALL timestop(handle)
   END SUBROUTINE allegro_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE allegro_add_force_virial(fist_nonbond_env, f_nonbond, pv_nonbond, use_virial)

      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial

      INTEGER                                            :: iat, iat_use
      TYPE(allegro_data_type), POINTER                   :: allegro_data

      CALL fist_nonbond_env_get(fist_nonbond_env, allegro_data=allegro_data)

      IF (use_virial) THEN
         pv_nonbond = pv_nonbond + allegro_data%virial
      END IF

      DO iat_use = 1, SIZE(allegro_data%use_indices)
         iat = allegro_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(f_nonbond, 2))
         f_nonbond(1:3, iat) = f_nonbond(1:3, iat) + allegro_data%force(1:3, iat_use)
      END DO

   END SUBROUTINE allegro_add_force_virial
END MODULE manybody_allegro


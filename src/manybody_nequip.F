!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      nequip implementation
!> \author Gabriele Tocci
! **************************************************************************************************
MODULE manybody_nequip

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type,&
                                              nequip_data_type,&
                                              pos_type
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              sp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential_types,            ONLY: nequip_pot_type,&
                                              nequip_type,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_type
   USE particle_types,                  ONLY: particle_type
   USE torch_api,                       ONLY: torch_dict_create,&
                                              torch_dict_get,&
                                              torch_dict_insert,&
                                              torch_dict_release,&
                                              torch_dict_type,&
                                              torch_model_eval,&
                                              torch_model_load
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: setup_nequip_arrays, destroy_nequip_arrays, &
             nequip_energy_store_force_virial, nequip_add_force_virial
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_nequip'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param potparm ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \param cell ...
!> \par History
!>      Implementation of the nequip potential - [gtocci] 2022
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE setup_nequip_arrays(nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_nequip_arrays'

      INTEGER                                            :: handle, i, iend, igrp, ikind, ilist, &
                                                            ipair, istart, jkind, nkinds, npairs, &
                                                            npairs_tot
      INTEGER, DIMENSION(:), POINTER                     :: work_list, work_list2
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rwork_list
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      CPASSERT(.NOT. ASSOCIATED(glob_loc_list))
      CPASSERT(.NOT. ASSOCIATED(glob_loc_list_a))
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
               IF (pot%type(i) == nequip_type) npairs_tot = npairs_tot + npairs
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
               ! NEQUIP
               IF (pot%type(i) == nequip_type) THEN
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
      CALL timestop(handle)
   END SUBROUTINE setup_nequip_arrays

! **************************************************************************************************
!> \brief ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \par History
!>      Implementation of the nequip potential - [gtocci] 2022
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE destroy_nequip_arrays(glob_loc_list, glob_cell_v, glob_loc_list_a)
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a

      IF (ASSOCIATED(glob_loc_list)) THEN
         DEALLOCATE (glob_loc_list)
      END IF
      IF (ASSOCIATED(glob_loc_list_a)) THEN
         DEALLOCATE (glob_loc_list_a)
      END IF
      IF (ASSOCIATED(glob_cell_v)) THEN
         DEALLOCATE (glob_cell_v)
      END IF

   END SUBROUTINE destroy_nequip_arrays

! **************************************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param nequip ...
!> \param glob_loc_list_a ...
!> \param r_last_update_pbc ...
!> \param pot_nequip ...
!> \param fist_nonbond_env ...
!> \param para_env ...
!> \par History
!>      Implementation of the nequip potential - [gtocci] 2022
!> \author Gabriele Tocci - University of Zurich
! **************************************************************************************************
   SUBROUTINE nequip_energy_store_force_virial(nonbonded, particle_set, cell, atomic_kind_set, &
                                               potparm, nequip, glob_loc_list_a, r_last_update_pbc, &
                                               pot_nequip, fist_nonbond_env, para_env)

      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(nequip_pot_type), POINTER                     :: nequip
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      REAL(kind=dp)                                      :: pot_nequip
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nequip_energy_store_force_virial'

      INTEGER :: atom_a, atom_b, handle, i, iat, iat_use, iend, ifirst, igrp, ikind, ilast, ilist, &
         ipair, istart, iunique, jkind, junique, mpair, n_atoms, n_atoms_use, nedges, nedges_tot, &
         nloc_size, npairs, nunique
      INTEGER(kind=int_8), ALLOCATABLE                   :: atom_types(:)
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:, :)  :: edge_index, t_edge_index, temp_edge_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: displ, displ_cell, edge_count, &
                                                            edge_count_cell
      INTEGER, DIMENSION(:), POINTER                     :: work_list
      INTEGER, DIMENSION(:, :), POINTER                  :: list, sort_list
      LOGICAL, ALLOCATABLE                               :: use_atom(:)
      REAL(kind=dp)                                      :: drij, lattice(3, 3), rab2_max, rij(3)
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi
      REAL(kind=sp), ALLOCATABLE                         :: pos(:, :)
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:, :)        :: edge_cell_shifts, temp_edge_cell_shifts
      REAL(sp), DIMENSION(:, :), POINTER                 :: atomic_energy, forces, total_energy
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(nequip_data_type), POINTER                    :: nequip_data
      TYPE(pair_potential_single_type), POINTER          :: pot
      TYPE(torch_dict_type)                              :: inputs, outputs

      CALL timeset(routineN, handle)

      NULLIFY (total_energy, atomic_energy, forces)
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      use_atom = .FALSE.

      DO ikind = 1, SIZE(atomic_kind_set)
      DO jkind = 1, SIZE(atomic_kind_set)
         pot => potparm%pot(ikind, jkind)%pot
         DO i = 1, SIZE(pot%type)
            IF (pot%type(i) /= nequip_type) CYCLE
            IF (.NOT. ASSOCIATED(nequip)) nequip => pot%set(i)%nequip
            DO iat = 1, n_atoms
               IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
                   particle_set(iat)%atomic_kind%kind_number == jkind) use_atom(iat) = .TRUE.
            END DO ! iat
         END DO ! i
      END DO ! jkind
      END DO ! ikind
      n_atoms_use = COUNT(use_atom)

      ! get nequip_data to save force, virial info and to load model
      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=nequip_data)
      IF (.NOT. ASSOCIATED(nequip_data)) THEN
         ALLOCATE (nequip_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, nequip_data=nequip_data)
         NULLIFY (nequip_data%use_indices, nequip_data%force)
         CALL torch_model_load(nequip_data%model, pot%set(1)%nequip%nequip_file_name)
      END IF
      IF (ASSOCIATED(nequip_data%force)) THEN
         IF (SIZE(nequip_data%force, 2) /= n_atoms_use) THEN
            DEALLOCATE (nequip_data%force, nequip_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(nequip_data%force)) THEN
         ALLOCATE (nequip_data%force(3, n_atoms_use))
         ALLOCATE (nequip_data%use_indices(n_atoms_use))
      END IF

      iat_use = 0
      DO iat = 1, n_atoms_use
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            nequip_data%use_indices(iat_use) = iat
         END IF
      END DO

      nedges = 0
      ALLOCATE (edge_index(2, SIZE(glob_loc_list_a)))
      ALLOCATE (edge_cell_shifts(3, SIZE(glob_loc_list_a)))
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         Kind_Group_Loop_Nequip: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= nequip_type) CYCLE
               rab2_max = pot%set(i)%nequip%rcutsq
               cell_v = MATMUL(cell%hmat, cvi)
               pot => potparm%pot(ikind, jkind)%pot
               nequip => pot%set(i)%nequip
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
                           edge_cell_shifts(:, nedges) = REAL(cvi, kind=sp)
                        END IF
                     END DO
                     ifirst = ilast + 1
                     IF (ipair <= npairs) junique = sort_list(1, ipair)
                  END DO
                  DEALLOCATE (sort_list, work_list)
               END IF
            END DO
         END DO Kind_Group_Loop_Nequip
      END DO

      ALLOCATE (edge_count(para_env%num_pe))
      ALLOCATE (edge_count_cell(para_env%num_pe))
      ALLOCATE (displ_cell(para_env%num_pe))
      ALLOCATE (displ(para_env%num_pe))

      CALL para_env%allgather(nedges, edge_count)
      nedges_tot = SUM(edge_count)

      ALLOCATE (temp_edge_index(2, nedges))
      temp_edge_index(:, :) = edge_index(:, :nedges)
      DEALLOCATE (edge_index)
      ALLOCATE (temp_edge_cell_shifts(3, nedges))
      temp_edge_cell_shifts(:, :) = edge_cell_shifts(:, :nedges)
      DEALLOCATE (edge_cell_shifts)

      ALLOCATE (edge_index(2, nedges_tot))
      ALLOCATE (edge_cell_shifts(3, nedges_tot))
      ALLOCATE (t_edge_index(nedges_tot, 2))

      edge_count_cell(:) = edge_count*3
      edge_count = edge_count*2
      displ(1) = 0
      displ_cell(1) = 0
      DO ipair = 2, para_env%num_pe
         displ(ipair) = displ(ipair - 1) + edge_count(ipair - 1)
         displ_cell(ipair) = displ_cell(ipair - 1) + edge_count_cell(ipair - 1)
      END DO

      CALL para_env%allgatherv(temp_edge_cell_shifts, edge_cell_shifts, edge_count_cell, displ_cell)
      CALL para_env%allgatherv(temp_edge_index, edge_index, edge_count, displ)

      t_edge_index(:, :) = TRANSPOSE(edge_index)
      DEALLOCATE (temp_edge_index, temp_edge_cell_shifts, edge_index)

      lattice = cell%hmat/pot%set(1)%nequip%unit_cell_val

      iat_use = 0
      ALLOCATE (pos(3, n_atoms_use), atom_types(n_atoms_use))

      DO iat = 1, n_atoms_use
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         atom_types(iat_use) = particle_set(iat)%atomic_kind%kind_number - 1
         pos(:, iat) = REAL(r_last_update_pbc(iat)%r(:)/pot%set(1)%nequip%unit_coords_val, kind=sp)
      END DO

      CALL torch_dict_create(inputs)
      CALL torch_dict_insert(inputs, "pos", pos)
      CALL torch_dict_insert(inputs, "edge_index", t_edge_index)
      CALL torch_dict_insert(inputs, "edge_cell_shift", edge_cell_shifts)
      CALL torch_dict_insert(inputs, "cell", REAL(lattice, kind=sp))
      CALL torch_dict_insert(inputs, "atom_types", atom_types)

      CALL torch_dict_create(outputs)

      CALL torch_model_eval(nequip_data%model, inputs, outputs)

      CALL torch_dict_get(outputs, "total_energy", total_energy)
      CALL torch_dict_get(outputs, "atomic_energy", atomic_energy)
      CALL torch_dict_get(outputs, "forces", forces)

      pot_nequip = REAL(total_energy(1, 1), kind=dp)*pot%set(1)%nequip%unit_energy_val
      total_energy = total_energy(1, 1)*REAL(pot%set(1)%nequip%unit_energy_val, kind=sp)

      CALL torch_dict_release(inputs)
      CALL torch_dict_release(outputs)

      DEALLOCATE (pos, t_edge_index, edge_cell_shifts, atom_types, total_energy, atomic_energy)

      ! account for double counting from multiple MPI processes
      IF (PRESENT(para_env)) pot_nequip = pot_nequip/REAL(para_env%num_pe, dp)
      IF (PRESENT(para_env)) THEN
         nequip_data%force = REAL(forces, kind=dp)* &
                             pot%set(1)%nequip%unit_forces_val/REAL(para_env%num_pe, dp)
      ELSE
         nequip_data%force = REAL(forces, kind=dp)*pot%set(1)%nequip%unit_forces_val
      END IF

      DEALLOCATE (forces)

      CALL timestop(handle)
   END SUBROUTINE nequip_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE nequip_add_force_virial(fist_nonbond_env, f_nonbond, pv_nonbond, use_virial)

      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: f_nonbond, pv_nonbond
      LOGICAL, INTENT(IN)                                :: use_virial

      INTEGER                                            :: iat, iat_use
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial
      TYPE(nequip_data_type), POINTER                    :: nequip_data

      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=nequip_data)

      IF (use_virial) THEN
         virial(1, 1) = 0.0_dp; virial(1, 2) = 0.0_dp; virial(1, 3) = 0.0_dp
         virial(2, 1) = 0.0_dp; virial(2, 2) = 0.0_dp; virial(2, 3) = 0.0_dp
         virial(3, 1) = 0.0_dp; virial(3, 2) = 0.0_dp; virial(3, 3) = 0.0_dp

         pv_nonbond(1, 1) = pv_nonbond(1, 1) + virial(1, 1)
         pv_nonbond(1, 2) = pv_nonbond(1, 2) + virial(1, 2)
         pv_nonbond(1, 3) = pv_nonbond(1, 3) + virial(1, 3)
         pv_nonbond(2, 1) = pv_nonbond(2, 1) + virial(2, 1)
         pv_nonbond(2, 2) = pv_nonbond(2, 2) + virial(2, 2)
         pv_nonbond(3, 1) = pv_nonbond(3, 1) + virial(3, 1)
         pv_nonbond(3, 2) = pv_nonbond(3, 2) + virial(3, 2)
         pv_nonbond(3, 3) = pv_nonbond(3, 3) + virial(3, 3)

         CPABORT("Stress tensor for NequIP not yet implemented")
      END IF

      DO iat_use = 1, SIZE(nequip_data%use_indices)
         iat = nequip_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(f_nonbond, 2))
         f_nonbond(1:3, iat) = f_nonbond(1:3, iat) + nequip_data%force(1:3, iat_use)
      END DO

   END SUBROUTINE nequip_add_force_virial
END MODULE manybody_nequip


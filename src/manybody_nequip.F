!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      tersoff implementation
!> \author CJM, I-Feng W. Kuo, Teodoro Laino
!>      nequip implementation
!> \author Gabriele Tocci
! **************************************************************************************************
MODULE manybody_nequip

   USE cell_types,                      ONLY: cell_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type                        
   USE fist_nonbond_env_types,          ONLY: pos_type, fist_nonbond_env_type
   USE kinds,                           ONLY: dp,sp
   USE mathconstants,                   ONLY: pi
   USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                              pair_potential_single_type,&
                                              nequip_pot_type,&
                                              nequip_type
   USE particle_types,                  ONLY: particle_type
   USE cp_para_types,                   ONLY: cp_para_env_type

   USE util,                            ONLY: sort
   USE physcon,                         ONLY: angstrom
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: setup_nequip_arrays, destroy_nequip_arrays, &
            nequip_energy_store_force_virial
!             nequip_forces, nequip_energy
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
!>      Fast implementation of the nequip potential - [tlaino] 2007
!> \author Teodoro Laino - University of Zurich
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
      WRITE(*, *) "npairs_tot ", npairs_tot
   END SUBROUTINE setup_nequip_arrays

! **************************************************************************************************
!> \brief ...
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \par History
!>      Fast implementation of the nequip potential - [tlaino] 2007
!> \author Teodoro Laino - University of Zurich
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
!> \param glob_loc_list ...
!> \param glob_cell_v ...
!> \param glob_loc_list_a ...
!> \par History
!>      Fast implementation of the nequip potential - [tlaino] 2007
!> \author Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE nequip_energy_store_force_virial(nonbonded, particle_set, cell, atomic_kind_set, &
                                              potparm, nequip, glob_loc_list_a, glob_loc_list, &
                                              glob_cell_v, r_last_update_pbc)

!                                             (nonbonded, potparm, glob_loc_list, glob_cell_v, glob_loc_list_a, cell)
      TYPE(nequip_pot_type), POINTER                     :: nequip
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      INTEGER, DIMENSION(:, :), POINTER                  :: glob_loc_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: glob_cell_v
      INTEGER, DIMENSION(:), POINTER                     :: glob_loc_list_a
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nequip_energy_store_force_virial'

!      INTEGER                                            :: handle, nkinds, npairs_tot
      INTEGER                                            :: i, iend, igrp, ikind, ilist, &
                                                            ipair, istart, jkind, npairs
      INTEGER, DIMENSION(:), POINTER                     :: work_list!, work_list2
      INTEGER, DIMENSION(:, :), POINTER                  :: list, sort_list
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi
!      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rwork_list
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pair_potential_single_type), POINTER          :: pot

      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
!      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
!      REAL(kind=dp)                                      :: pot_nequip
!      TYPE(cp_para_env_type), OPTIONAL, &
!         POINTER                                         :: para_env
      
      INTEGER                                            :: iat, iat_use, n_atoms_use, atom_a, atom_b, &
                                                            ifirst, ilast, iunique, nunique, junique, &
                                                            mpair, nloc_size, n_atoms, nedges
      LOGICAL, ALLOCATABLE                               :: use_atom(:)

      REAL(kind=dp), ALLOCATABLE                         :: force(:, :), pos(:, :)
      REAL(kind=dp)                                      :: lattice(3, 3), rij(3), drij, rab2_max, pot_loc
      CHARACTER(len=2), ALLOCATABLE                      :: elem_symbol(:)
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc
      INTEGER, ALLOCATABLE                               :: edge_index(:,:), new_edge_index(:,:)
      REAL(kind=sp), ALLOCATABLE                         :: edge_cell_shifts(:,:), new_edge_cell_shifts(:,:) 

!      IF (output_unit > 0) THEN
!         WRITE (output_unit, '(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy [a.u.]: ",T55,F26.15)') &
!            ADJUSTR(TRIM(use_prog_name(force_env%in_use))), e_pot
!      END IF
!         WRITE(*, *) "glob_cell_v ", SIZE(glob_cell_v) , glob_cell_v
!         WRITE(*, *) "cell ", cell

      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      use_atom = .FALSE.

      WRITE(*,*) "entering loop over kind"

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
      ALLOCATE (pos(3, n_atoms_use), force(3, n_atoms_use), elem_symbol(n_atoms_use))

      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         pos(1:3, iat_use) = particle_set(iat)%r*angstrom
         elem_symbol(iat_use) = particle_set(iat)%atomic_kind%element_symbol
      END DO
      lattice = cell%hmat*angstrom
      nedges = 0
      ALLOCATE (edge_index(SIZE(glob_loc_list_a),2))
      ALLOCATE (edge_cell_shifts(3,SIZE(glob_loc_list_a)))
      DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         WRITE(*, *) "npairs, nonbonded%nlists ", npairs, nonbonded%nlists
         IF (npairs == 0) CYCLE
         Kind_Group_Loop_Nequip: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)
            list => neighbor_kind_pair%list
            cvi = neighbor_kind_pair%cell_vector
            pot => potparm%pot(ikind, jkind)%pot
            WRITE(*,*) "istart, iend: ", istart, iend
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= nequip_type) CYCLE
               rab2_max = pot%set(i)%nequip%rcutsq
               WRITE(*, *) "cutoff manybody_potential ", rab2_max
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
                        ! Energy terms
                        pot_loc = 0.0_dp
                        rij(:) = r_last_update_pbc(atom_b)%r(:) - r_last_update_pbc(atom_a)%r(:) + cell_v
                        drij = DOT_PRODUCT(rij, rij)
                        ipair = ipair + 1
                        IF (drij <= rab2_max) THEN
                           WRITE(*, *) "calc rij: ", drij
                           WRITE(*, *) "glob_loc_list ", SIZE(glob_loc_list)
                           WRITE(*, *) "glob_loc_list_a ", SIZE(glob_loc_list_a)
                           WRITE(*, *) "glob_loc_list_a(atom_b), atom_b ", glob_loc_list_a(atom_b), atom_b
                           WRITE(*,*) "r_last_update_pbc: ", r_last_update_pbc(atom_a)%r(:)*angstrom, atom_a
                           WRITE(*,*) "r_last_update_pbc: ", r_last_update_pbc(atom_b)%r(:)*angstrom, atom_b
                           WRITE(*, *) "glob_cell_v ", SIZE(glob_cell_v)
                           nedges = nedges + 1
                           edge_index(nedges, :) = [atom_a - 1, atom_b - 1]
                           edge_cell_shifts(:, nedges) = REAL(cvi, kind=sp) 
                        END IF
!                           IF (drij > rab2_max) CYCLE
!                           drij = SQRT(drij)
!                           CALL tersoff_energy(pot_loc, tersoff, r_last_update_pbc, atom_a, atom_b, nloc_size, &
!                                               glob_loc_list(:, ifirst:ilast), glob_cell_v(:, ifirst:ilast), cell_v, drij)
!                        pot_manybody = pot_manybody + pot_nequip
                     END DO
                     ifirst = ilast + 1
                     IF (ipair <= npairs) junique = sort_list(1, ipair)
                  END DO
                  DEALLOCATE (sort_list, work_list)
               END IF
            END DO
         END DO Kind_Group_Loop_Nequip
      END DO   
      WRITE(*, *) "nequip nedges: ", nedges 
      WRITE(*, *) "nequip npairs: ", npairs
      ALLOCATE (new_edge_index(nedges,2))
      new_edge_index = edge_index(:nedges, :)
      DEALLOCATE (edge_index)
      ALLOCATE (new_edge_cell_shifts(3, nedges))
      new_edge_cell_shifts = edge_cell_shifts(:,:nedges)
      DEALLOCATE (edge_cell_shifts)
      DO ipair = 1, nedges
         atom_a = new_edge_index(ipair,1) + 1
         atom_b = new_edge_index(ipair,2) + 1
         WRITE(*,*) new_edge_index(ipair, :), r_last_update_pbc(atom_a)%r(:)*angstrom, &
                    r_last_update_pbc(atom_b)%r(:)*angstrom
      END DO
      DO ipair = 1, nedges
         WRITE(*,*) new_edge_cell_shifts(:,ipair)
      END DO

   END SUBROUTINE nequip_energy_store_force_virial

END MODULE manybody_nequip


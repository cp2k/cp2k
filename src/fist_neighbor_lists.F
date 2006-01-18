!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_neighbor_lists [1.0] *
!!
!!   NAME
!!     fist_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists for FIST.
!!
!!   AUTHOR
!!     MK (19.11.2002,24.07.2003)
!!
!!   MODIFICATION HISTORY
!!     11.02.2005: build and update merged (MK)
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             init_neighbor_list_set,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE global_types,                    ONLY: LOW,&
                                             MEDIUM,&
                                             global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "fist_neighbor_lists"

! *** Public subroutines ***

  PUBLIC :: build_fist_neighbor_lists

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_fist_neighbor_lists(atomic_kind_set,particle_set,&
                                       local_particles,cell,r_max,nonbonded,&
                                       globenv,build_from_scratch,error)

!   Purpose: Build all the required neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nonbonded
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(IN)                      :: build_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_fist_neighbor_lists"

    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, maxatom, &
      maxatom_local, natom_a, natom_local_a, nkind, output_unit
    LOGICAL                                  :: print_subcell_grid
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("Entering "//routine_name,globenv,error)

    CALL timeset(routine_name,"I","",handle)

    output_unit = globenv%scr

    IF ((globenv%ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > MEDIUM)) THEN
      print_subcell_grid = .TRUE.
    ELSE
      print_subcell_grid = .FALSE.
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxatom=maxatom)

!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "atom",nkind*int_size)

!   *** Initialize the local data structures ***

    maxatom_local = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom(ikind)%list)

      natom_a = SIZE(atom(ikind)%list)
      natom_local_a = local_particles%n_el(ikind)

      maxatom_local = MAX(maxatom_local,natom_local_a)

      IF (natom_local_a > 0) THEN

        ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_local_a_index",&
                                         natom_local_a*int_size)

!       *** Build index vector for mapping ***

        DO iatom_local=1,natom_local_a
          atom_a = local_particles%list(ikind)%array(iatom_local)
          DO iatom=1,natom_a
            IF (atom_a == atom(ikind)%list(iatom)) THEN
              atom(ikind)%list_local_a_index(iatom_local) = iatom
              EXIT
            END IF
          END DO
        END DO

      END IF

      ALLOCATE (atom(ikind)%r_pbc(3,natom_a),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom(ikind)%r_pbc",3*natom_a*dp_size)

      ALLOCATE (atom(ikind)%s_pbc(3,natom_a),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom(ikind)%s_pbc",3*natom_a*dp_size)

!     *** Calculate PBC coordinates ***

      DO iatom=1,natom_a
        atom_a = atom(ikind)%list(iatom)
        r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
        atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
        atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell)
      END DO

    END DO

!   *** Build the nonbonded neighbor lists ***

    CALL build_neighbor_lists(nonbonded,particle_set,atom,cell,maxatom_local,&
                              print_subcell_grid,output_unit,r_max,&
                              build_from_scratch)

    IF (globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(nonbonded,particle_set,cell,globenv,&
                                "NONBONDED NEIGHBOR LISTS")
    END IF

!   *** Print particle distribution ***

    CALL write_distribution(nonbonded,atomic_kind_set,globenv)

!   *** Release work storage ***

    DO ikind=1,nkind
      NULLIFY (atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_a_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE (atom(ikind)%r_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE (atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "atom")

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("Leaving "//routine_name,globenv,error)

  END SUBROUTINE build_fist_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_neighbor_lists(ab,particle_set,atom,cell,maxatom_local,&
                                  print_subcell_grid,output_unit,r_max,&
                                  build_from_scratch)

!   Purpose: Build simple pair neighbor lists.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max
    LOGICAL, INTENT(IN)                      :: build_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom_local, jcell, jkind, k, &
      kcell, natom_local_a, natom_local_b, nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b, ncell, &
                                                nsubcell, periodic
    LOGICAL                                  :: cell000, include_ab
    REAL(dp)                                 :: rab2, rab2_max, rab_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: abc, r, rab, rb, s, sab, &
                                                sab_max, sb, sb_max, sb_min, &
                                                sb_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    nkind = SIZE(atom)

    IF (build_from_scratch) THEN

!     *** Deallocate the old neighbor list structure ***

      IF (ASSOCIATED(ab)) THEN
        DO iab=1,SIZE(ab)
          CALL deallocate_neighbor_list_set(ab(iab)%neighbor_list_set)
        END DO
        DEALLOCATE (ab,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ab")
      END IF

!     *** Allocate and initialize the new neighbor list structure ***

      ALLOCATE (ab(nkind*nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ab",nkind*nkind*int_size)

      DO iab=1,SIZE(ab)
        NULLIFY (ab(iab)%neighbor_list_set)
      END DO

    ELSE

      IF (.NOT.ASSOCIATED(ab)) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "ab is not associated and cannot be updated")
      END IF

    END IF

!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  abc=abc,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***

    IF (print_subcell_grid) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE NONBONDED NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO jkind=1,nkind

        iab = ikind + nkind*(jkind - 1)

        natom_local_b = SIZE(atom(jkind)%list)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = r_max(ikind,jkind)
        rab2_max = rab_max*rab_max

        sab_max(1) = rab_max/plane_distance(1,0,0,cell)
        sab_max(2) = rab_max/plane_distance(0,1,0,cell)
        sab_max(3) = rab_max/plane_distance(0,0,1,cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MIN(MAX(1,NINT(0.5_dp*subcells/sab_max(:))),20)

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = s(:)*abc(:)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        IF (build_from_scratch) THEN
          CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          ab(iab)%neighbor_list_set,&
                                          r_max=rab_max)
          neighbor_list_set => ab(iab)%neighbor_list_set
        ELSE
          neighbor_list_set => ab(iab)%neighbor_list_set
          IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE
          CALL init_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                      r_max=rab_max)
        END IF

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            atom_b = atom(jkind)%list(jatom_local)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom_local)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop1_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_b(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop1_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_b(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop1_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_b(1) = icell

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  rb(:) = scaled_to_real(sb(:),cell)

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (cell000) THEN
                      IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                      END IF
                    ELSE
                      IF (atom_a == atom_b) THEN
                        include_ab = .TRUE.
                      ELSE IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                      END IF
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max) THEN
                        sab(:) = real_to_scaled(rab(:),cell)
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom_local)%neighbor_list,&
                          neighbor=atom_b,&
                          cell=cell_b,&
                          r=sab,&
                          exclusion_list = particle_set(atom_a)%lpi%&
                                           exclusion(jkind)%list_vdw)
                      END IF
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = MIN(INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp),nsubcell(1))
            j = MIN(INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp),nsubcell(2))
            k = MIN(INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp),nsubcell(3))
            subcell(i,j,k)%natom  = subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            atom_b = atom(jkind)%list(jatom_local)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom_local)

            loop2_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop2_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop2_kcell
              cell_b(3) = kcell

              loop2_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop2_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop2_jcell
                cell_b(2) = jcell

                loop2_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop2_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop2_icell
                  cell_b(1) = icell

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  rb(:) = scaled_to_real(sb(:),cell)

                  loop_k: DO k=1,nsubcell(3)
                    loop_j: DO j=1,nsubcell(2)
                      loop_i: DO i=1,nsubcell(1)

                        IF (subcell(i,j,k)%natom == 0) CYCLE
 
                        IF (sb_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                        IF (sb_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k

                        IF (sb_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                        IF (sb_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j

                        IF (sb_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                        IF (sb_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i

                        DO iatom_subcell=1,subcell(i,j,k)%natom
                          iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          atom_a = atom(ikind)%list(iatom)
                          IF (cell000) THEN
                            IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                            END IF
                          ELSE
                            IF (atom_a == atom_b) THEN
                              include_ab = .TRUE.
                            ELSE IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                            END IF
                          END IF
                          IF (include_ab) THEN
                            rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                            IF (rab2 < rab2_max) THEN
                              sab(:) = real_to_scaled(rab(:),cell)
                              CALL add_neighbor_node(&
                                neighbor_list=kind_a(iatom_local)%neighbor_list,&
                                neighbor=atom_b,&
                                cell=cell_b,&
                                r=sab,&
                                exclusion_list=particle_set(atom_a)%lpi%&
                                               exclusion(jkind)%list_vdw)
                            END IF
                          END IF
                        END DO

                      END DO loop_i
                    END DO loop_j
                  END DO loop_k

                END DO loop2_icell
              END DO loop2_jcell
            END DO loop2_kcell

          END DO

          CALL deallocate_subcell(subcell)

        END IF

      END DO
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_neighbor_lists

! *****************************************************************************

  SUBROUTINE write_distribution(ab,atomic_kind_set,globenv)

!   Purpose: Print the distribution of the simple pair neighbor list.

!   History: - Creation (19.06.2003,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      pointer                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_distribution"

    INTEGER :: group, iab, ikind, ilist, ipe, istat, jkind, mype, n, nkind, &
      nlist, nnode, npe, nsgf_a, nsgf_b, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblock, nelement
    LOGICAL                                  :: ionode
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

!   ---------------------------------------------------------------------------
!   *** Quick return, if no printing is requested ***

    IF (.NOT.globenv%print%distribution) RETURN

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

!   *** Allocate work storage ***

    ALLOCATE (nblock(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock",npe*int_size)
    nblock(:) = 0

    ALLOCATE (nelement(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement",npe*int_size)
    nelement(:) = 0

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(orb_basis_set)) THEN
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf_a)
      ELSE
        nsgf_a = 1
      END IF

      DO jkind=1,nkind

        iab = ikind + nkind*(jkind - 1)

        neighbor_list_set => ab(iab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (ASSOCIATED(orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf_b)
        ELSE
          nsgf_b = 1
        END IF

        n = nsgf_a*nsgf_b

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 nnode=nnode)

          nblock(mype) = nblock(mype) + nnode
          nelement(mype) = nelement(mype) + n*nnode

          neighbor_list => next(neighbor_list)

        END DO

      END DO

    END DO

    CALL mp_sum(nblock,group)
    CALL mp_sum(nelement,group)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I6,T27,I10,T55,I10))")&
        "DISTRIBUTION OF THE NEIGHBOR LISTS",&
        "Process   Number of particle pairs   Number of matrix elements",&
        (ipe-1,nblock(ipe),nelement(ipe),ipe=1,npe)
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,T55,I10)")&
        "Sum",SUM(nblock),SUM(nelement)
    END IF

!   *** Release work storage ***

    DEALLOCATE (nblock,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock")

    DEALLOCATE (nelement,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement")

  END SUBROUTINE write_distribution

! *****************************************************************************

  SUBROUTINE write_neighbor_lists(ab,particle_set,cell,globenv,name)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002,MK)
!            - Adapted to the new parallelized neighbor list version
!              (26.06.2003,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), INTENT(IN)               :: ab
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      pointer                             :: globenv
    CHARACTER(LEN=*), INTENT(IN)             :: name

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: atom_a, atom_b, group, iab, &
                                                ilist, inode, ipe, mype, &
                                                nlist, nneighbor, nnode, npe, &
                                                output_unit
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b
    LOGICAL                                  :: ionode, print_headline
    REAL(dp)                                 :: dab, unit_of_length
    REAL(dp), DIMENSION(3)                   :: ra, rab, rb
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    CALL mp_sync(group)
    IF (ionode) CALL m_flush(output_unit)
    CALL mp_sync(group)

!   *** Create the headline ***

    string = ""

    IF (npe == 1) THEN
      string = TRIM(name)//" IN "//TRIM(unit_of_length_name)
    ELSE
      WRITE (UNIT=string,FMT="(A,I5,A)")&
        TRIM(name)//" IN "//TRIM(unit_of_length_name)//" (PROCESS",mype-1,")"
        CALL compress(string)
    END IF

    print_headline = .TRUE.

!   *** All processes print now their neighbor lists portion consecutively ***

    nneighbor = 0

    DO ipe=1,npe

      IF (ipe == mype) THEN

!       *** Print first part of headline ***

        if (ionode) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(string)

        DO iab=1,SIZE(ab)

          neighbor_list_set => ab(iab)%neighbor_list_set

          IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

!         *** Loop over all atoms and their corresponding neighbor lists ***

          CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                     nlist=nlist)

          neighbor_list => first_list(neighbor_list_set)

          DO ilist=1,nlist

            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   nnode=nnode)
            nneighbor = nneighbor + nnode

            IF (ionode) THEN
               IF (globenv%print%level > LOW) THEN

                  !Print second part of headline ***

                  IF (print_headline) THEN
                     WRITE (UNIT=output_unit,FMT="(/,T3,A,7X,A,2(11X,A),10X,A)")&
                          "Atom  Neighbors  Cell(i,j,k)","X","Y","Z","Distance"
                     print_headline = .FALSE.
                  END IF
                  
                  ra(:) = pbc(particle_set(atom_a)%r,cell,cell_a)
                  
                  WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,I6,3X,3I4,3F12.6)")&
                       atom_a,nnode,cell_a(1:3),ra(1:3)/unit_of_length
                  
                  !             *** Direct the work pointer to the  ***
                  !             *** start point of the current list ***
                  
                  neighbor_node => first_node(neighbor_list)
                  
                  !             *** Traverse the neighbor list of the current ***
                  !             *** atom and print the stored information     ***
                  
                  IF (globenv%print%level > MEDIUM) THEN
                     
                     DO inode=1,nnode
                        
                        CALL get_neighbor_node(neighbor_node=neighbor_node,&
                             neighbor=atom_b,&
                             cell=cell_b,&
                             r=rab)
                        
                        rb(:) = ra(:) + rab(:)
                        dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
                        
                        WRITE (UNIT=output_unit,&
                             FMT="(T10,I6,3X,3I4,3F12.6,2X,F12.6)")&
                             atom_b,cell_b(1:3),rb(1:3)/unit_of_length,dab/unit_of_length
                        
                        neighbor_node => next(neighbor_node)
                        
                     END DO

                  END IF
                  
               END IF
            END IF

            neighbor_list => next(neighbor_list)

          END DO

        END DO

        string = ""

        IF (npe == 1) THEN
          WRITE (UNIT=string,FMT="(A,I12)")&
            "Total number of neighbor interactions:",nneighbor
        ELSE
          WRITE (UNIT=string,FMT="(A,I12,A,I12)")&
            "Total number of neighbor interactions for process",mype-1,":",&
            nneighbor
        END IF
        CALL compress(string)
        IF (ionode) WRITE (UNIT=output_unit,FMT="(/,T2,A)") TRIM(string)

      END IF

      CALL mp_sync(group)
      IF (ionode) CALL m_flush(output_unit)
      CALL mp_sync(group)

    END DO

  END SUBROUTINE write_neighbor_lists

! *****************************************************************************

END MODULE fist_neighbor_lists

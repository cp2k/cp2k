!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
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
!!     Matthias Krack (19.11.2002,24.07.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_list_types,          ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type,&
                                             write_distribution,&
                                             write_neighbor_lists
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

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
                                       local_particles,cell,r_cut,nonbonded,&
                                       globenv)

!   Purpose: Build all the required neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_cut
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nonbonded
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_fist_neighbor_lists"

    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, maxatom, &
      maxatom_local, natom_a, natom_local_a, nkind, output_unit
    LOGICAL                                  :: print_subcell_grid
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    output_unit = globenv%scr

    IF ((globenv%ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
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
      NULLIFY (atom(ikind)%list_local_b_index)
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

!         *** Build index vector for mapping ***

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

!       *** Calculate PBC coordinates ***

        DO iatom=1,natom_a
          atom_a = atom(ikind)%list(iatom)
          r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
          atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell)
        END DO

    END DO

!   *** Build the nonbonded neighbor lists ***

    CALL build_neighbor_lists(nonbonded,particle_set,atom,cell,maxatom_local,&
                              print_subcell_grid,output_unit,r_cut)

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

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_fist_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_neighbor_lists(ab,particle_set,atom,cell,maxatom_local,&
                                  print_subcell_grid,output_unit,r_cut)

!   Purpose: Build simple pair neighbor lists.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE ( particle_type ), DIMENSION ( : ), POINTER :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_cut

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom_local, jcell, &
      jkind, k, kcell, natom_local_a, natom_local_b, nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b, ncell, &
                                                nsubcell, periodic
    LOGICAL                                  :: include_ab
    REAL(dp)                                 :: rab2, rab2_max, rab_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: r, rab, rb, s, sab, sab_max, &
                                                sb, sb_max, sb_min, sb_pbc, abc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

!   *** Deallocate the old neighbor list structure ***

    IF (ASSOCIATED(ab)) THEN
      DO iab=1,SIZE(ab)
        CALL deallocate_neighbor_list_set(ab(iab)%neighbor_list_set)
      END DO
      DEALLOCATE (ab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ab")
    END IF

!   *** Allocate and initialize the new neighbor list structure ***

    nkind = SIZE(atom)

    ALLOCATE (ab(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ab",nkind*nkind*int_size)

    DO iab=1,SIZE(ab)
      NULLIFY (ab(iab)%neighbor_list_set)
    END DO

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

        rab_max = r_cut(ikind,jkind)
        rab2_max = rab_max*rab_max

        sab_max(1) = rab_max/plane_distance(1,0,0,cell)
        sab_max(2) = rab_max/plane_distance(0,1,0,cell)
        sab_max(3) = rab_max/plane_distance(0,0,1,cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = s(:)*abc(:)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                        ab(iab)%neighbor_list_set,&
                                        r_max=rab_max)
        neighbor_list_set => ab(iab)%neighbor_list_set

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

                  rb(:) = scaled_to_real(sb(:),cell)

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (atom_a > atom_b) THEN
                      include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                    ELSE
                      include_ab = (MODULO(atom_a + atom_b,2) == 0)
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max) THEN
                        IF (rab2 > 1.0E-6_dp) THEN
                          sab(:) = real_to_scaled(rab(:),cell)
                          CALL add_neighbor_node(&
                            neighbor_list=kind_a(iatom_local)%neighbor_list,&
                            neighbor=atom_b,&
                            cell=cell_b,&
                            r=rab,&
                            s=sab,&
                            exclusion_list = particle_set ( atom_a ) % lpi %  &
                                             exclusion ( jkind ) % list )
                         END IF
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
                          IF (atom_a > atom_b) THEN
                            include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                          ELSE
                            include_ab = (MODULO(atom_a + atom_b,2) == 0)
                          END IF
                          IF (include_ab) THEN
                            rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                            IF (rab2 < rab2_max) THEN
                              IF (rab2 > 1.0E-6_dp) THEN
                                sab(:) = real_to_scaled(rab(:),cell)
                                CALL add_neighbor_node(&
                                  neighbor_list=kind_a(iatom_local)%neighbor_list,&
                                  neighbor=atom_b,&
                                  cell=cell_b,&
                                  r=rab,&
                                  s=sab,&
                                  exclusion_list=particle_set ( atom_a ) % lpi %  &
                                                 exclusion ( jkind ) % list )
                              END IF
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

END MODULE fist_neighbor_lists

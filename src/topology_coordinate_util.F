!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collection of subroutine needed for topology related things
!> \par History
!>     jgh (23-05-2004) Last atom of molecule information added
! *****************************************************************************
MODULE topology_coordinate_util
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE exclusion_types,                 ONLY: exclusion_type
  USE external_potential_types,        ONLY: allocate_potential,&
                                             fist_potential_type,&
                                             get_potential,&
                                             set_potential
  USE input_constants,                 ONLY: do_fist,&
                                             do_skip_12,&
                                             do_skip_13,&
                                             do_skip_14
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             set_molecule_kind
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE physcon,                         ONLY: massunit
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: array1_list_type,&
                                             reorder_structure
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_coordinate_util'

  PRIVATE
  PUBLIC :: topology_coordinate_pack

CONTAINS

! *****************************************************************************
!> \brief Take info readin from different file format and stuff it into
!>      compatible data structure in cp2k
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param topology ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param subsys_section ...
!> \param force_env_section ...
!> \param exclusions ...
!> \param ignore_outside_box ...
!> \param error ...
!> \par History
!>      Teodoro Laino - modified in order to optimize the list of molecules
!>                      to build the exclusion lists
! *****************************************************************************
  SUBROUTINE topology_coordinate_pack (particle_set,atomic_kind_set,&
       molecule_kind_set,molecule_set,topology,qmmm,qmmm_env,&
       subsys_section,force_env_section,exclusions,ignore_outside_box,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section, &
                                                force_env_section
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exclusions
    LOGICAL, INTENT(IN), OPTIONAL            :: ignore_outside_box
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_coordinate_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname
    INTEGER :: atom_i, atom_j, counter, dim0, dim1, dim2, dim3, first, &
      handle, handle2, i, iatom, ikind, iw, j, k, last, method_name_id, n, &
      natom, stat
    INTEGER, DIMENSION(:), POINTER           :: iatomlist, id_element, &
                                                id_work, kind_of, list, &
                                                list2, molecule_list, &
                                                natom_of_kind, wlist
    INTEGER, DIMENSION(:, :), POINTER        :: pairs
    LOGICAL :: autogen, check, disable_exclusion_lists, do_center, explicit, &
      failure, found, my_ignore_outside_box, my_qmmm, &
      present_12_excl_ei_list, present_12_excl_vdw_list
    REAL(KIND=dp)                            :: bounds(2,3), cdims(3), &
                                                dims(3), qeff, vec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge, cpoint, mass
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: ex_bend_list, ex_bond_list, &
                                                ex_bond_list_ei, &
                                                ex_bond_list_vdw, ex_onfo_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: exclude_section, &
                                                topology_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    topology_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY",error=error)
    CALL timeset(routineN,handle)

    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm=qmmm
    atom_info => topology%atom_info
    conn_info => topology%conn_info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Determine topology%[natom_type,atom_names] and save mass(natom_type)
    !    and element(natom_type)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_1',handle2)
    counter=0
    NULLIFY(id_work,mass,id_element,charge)
    ALLOCATE(id_work(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mass(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(id_element(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(charge(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    id_work=str2id(s2s(""))
    IF (iw>0) WRITE(iw,*)"molecule_kind_set ::",SIZE(molecule_kind_set)
    DO i = 1, SIZE(molecule_kind_set)
       j=molecule_kind_set(i)%molecule_list(1)
       molecule => molecule_set(j)
       molecule_kind => molecule_set(j)%molecule_kind
       IF (iw>0) WRITE(iw,*)"molecule number ::",j," has molecule kind number ::",i
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            natom=natom,atom_list=atom_list)
       CALL get_molecule(molecule=molecule,&
            first_atom=first,last_atom=last)
       IF (iw>0) WRITE(iw,*)"boundaries of molecules (first, last) ::",first,last
       DO j = 1, natom
          IF(.NOT.ANY(id_work(1:counter).EQ.atom_list(j)%id_name)) THEN
             counter = counter + 1
             id_work(counter)   = atom_list(j)%id_name
             mass(counter)   = atom_info%atm_mass(first+j-1)
             id_element(counter)= atom_info%id_element(first+j-1)
             charge(counter) = atom_info%atm_charge(first+j-1)
             IF(iw>0) WRITE(iw,'(7X,A,1X,A5,F10.5,5X,A2,5X,F10.5)') &
                  "NEW ATOMIC KIND",id2str(id_work(counter)),mass(counter),id2str(id_element(counter)),charge(counter)
          ELSE
             found = .FALSE.
             DO k = 1, counter
                IF((id_work(k)==atom_list(j)%id_name).AND.(charge(k)==atom_info%atm_charge(first+j-1))) THEN
                   found = .TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.found) THEN
                counter = counter + 1
                id_work(counter)   = atom_list(j)%id_name
                mass(counter)   = atom_info%atm_mass(first+j-1)
                id_element(counter)= atom_info%id_element(first+j-1)
                charge(counter) = atom_info%atm_charge(first+j-1)
                IF(iw>0) WRITE(iw,'(7X,A,1X,A5,F10.5,5X,A2,5X,F10.5)') &
                     "NEW ATOMIC KIND",id2str(id_work(counter)),mass(counter),id2str(id_element(counter)),charge(counter)
             END IF
          END IF
       END DO
    END DO
    topology%natom_type = counter
    ALLOCATE(atom_info%id_atom_names(topology%natom_type),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO k=1,counter
       atom_info%id_atom_names(k)=id_work(k)
    ENDDO
    DEALLOCATE(id_work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL reallocate(mass,1,counter)
    CALL reallocate(id_element,1,counter)
    CALL reallocate(charge,1,counter)
    IF(iw>0) &
         WRITE(iw,'(5X,A,I3)') "Total Number of Atomic Kinds = ",topology%natom_type
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Allocate the data structure for the atomic kind information
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_2',handle2)
    NULLIFY(atomic_kind_set)
    ALLOCATE(atomic_kind_set(topology%natom_type))
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3.  Allocate the data structure for the atomic information
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_3',handle2)
    NULLIFY(particle_set)
    CALL allocate_particle_set(particle_set,topology%natoms,error)
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the atomic_kind_set(ikind)%[name,kind_number,mass]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_4',handle2)
    DO i=1,topology%natom_type
       atomic_kind => atomic_kind_set(i)
       mass(i) = mass(i) * massunit
       CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=i,&
            name=id2str(atom_info%id_atom_names(i)),element_symbol=id2str(id_element(i)),&
            mass=mass(i))
       IF (iw>0) THEN
          WRITE(iw,'(A,I5,A,I5,4A)')"Atomic Kind n.:",i," out of:",topology%natom_type,&
               " name:   ",TRIM(id2str(atom_info%id_atom_names(i))),"   element:   ",&
               TRIM(id2str(id_element(i)))
       END IF
    END DO
    DEALLOCATE(mass,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(id_element,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Determine number of atom of each kind (ie natom_of_kind and kind_of)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_5',handle2)
    ALLOCATE(kind_of(topology%natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(natom_of_kind(topology%natom_type),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    kind_of(:) = 0
    natom_of_kind(:) = 0
    DO i=1,topology%natom_type
       DO j=1,topology%natoms
          IF((atom_info%id_atom_names(i)==atom_info%id_atmname(j)).AND.(charge(i)==atom_info%atm_charge(j))) THEN
             natom_of_kind(i)=natom_of_kind(i)+1
             IF (kind_of(j)==0) kind_of(j)=i
          END IF
       END DO
    END DO
    IF (ANY(kind_of==0)) THEN
       DO i = 1, topology%natoms
          IF (kind_of(i)==0) THEN
             WRITE(*,*)i,kind_of(i)
             WRITE(*,*)"Two molecules have been defined as identical molecules but atoms mismatch charges!!"
          END IF
       END DO
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the atom_kind_set(ikind)%[natom,atom_list]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_6',handle2)
    DO i=1,topology%natom_type
       atomic_kind => atomic_kind_set(i)
       NULLIFY (iatomlist)
       ALLOCATE(iatomlist(natom_of_kind(i)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       counter=0
       DO j=1,topology%natoms
          IF(kind_of(j)==i) THEN
             counter=counter+1
             iatomlist(counter)=j
          END IF
       END DO
       IF(iw>0) THEN
          WRITE(iw,'(A,I6,A)') "      Atomic kind ",i," contains particles"
          DO J=1,SIZE(iatomlist)
             IF (MOD(J,5).EQ.0) THEN ! split long lines
                WRITE(iw,'(I12)') iatomlist(J)
             ELSE
                WRITE(iw,'(I12)',ADVANCE="NO") iatomlist(J)
             ENDIF
          ENDDO
          WRITE(iw,*)
       ENDIF
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
            natom=natom_of_kind(i),atom_list=iatomlist)
       DEALLOCATE(iatomlist,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(natom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Possibly center the coordinates and fill in coordinates in particle_set
    !-----------------------------------------------------------------------------
    CALL section_vals_val_get(subsys_section,&
         "TOPOLOGY%CENTER_COORDINATES%_SECTION_PARAMETERS_",l_val=do_center,&
         error=error)
    CALL timeset(routineN//'_7a',handle2)
    bounds(1,1) = MINVAL(atom_info%r(1,:))
    bounds(2,1) = MAXVAL(atom_info%r(1,:))

    bounds(1,2) = MINVAL(atom_info%r(2,:))
    bounds(2,2) = MAXVAL(atom_info%r(2,:))

    bounds(1,3) = MINVAL(atom_info%r(3,:))
    bounds(2,3) = MAXVAL(atom_info%r(3,:))

    dims     = bounds(2,:)-bounds(1,:)
    cdims(1) = topology%cell%hmat(1,1)
    cdims(2) = topology%cell%hmat(2,2)
    cdims(3) = topology%cell%hmat(3,3)
    IF (iw>0) THEN
       WRITE(iw,'(A,3F12.6)')"System sizes: ",dims,"Cell sizes (diagonal): ",cdims
    END IF
    check = .TRUE.
    DO i = 1, 3
       IF (topology%cell%perd(i)==0) THEN
         check = check.AND.(dims(i)<cdims(i))
       END IF
    END DO
    my_ignore_outside_box = .FALSE.
    IF (PRESENT(ignore_outside_box)) my_ignore_outside_box = ignore_outside_box
    CALL cp_assert(my_ignore_outside_box .OR. check,cp_failure_level,cp_assertion_failed,routineP,&
         "A non-periodic calculation has been requested but the system size "//&
         "exceeds the cell size in at least one of the non-periodic directions! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    IF (do_center) THEN
       CALL section_vals_val_get(subsys_section,&
            "TOPOLOGY%CENTER_COORDINATES%CENTER_POINT",explicit=explicit,&
            error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsys_section,&
               "TOPOLOGY%CENTER_COORDINATES%CENTER_POINT",r_vals=cpoint,&
               error=error)
          vec = cpoint
       ELSE
          vec = cdims/2.0_dp
       END IF
       dims = (bounds(2,:)+bounds(1,:))/2.0_dp - vec
    ELSE
       dims = 0.0_dp
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//'_7b',handle2)
    DO i=1,topology%natoms
       ikind = kind_of(i)
       IF (iw>0) THEN
          WRITE(iw,*)"atom number :: ",i,"kind number ::",ikind
       END IF
       particle_set(i)%atomic_kind => atomic_kind_set(ikind)
       particle_set(i)%r(:) = atom_info%r(:,i) - dims
       particle_set(i)%atom_index = i
    END DO
    CALL timestop(handle2)
    DEALLOCATE(kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Fill in the exclusions%list_exclude_vdw
    ! 9. Fill in the exclusions%list_exclude_ei
    ! 10. Fill in the exclusions%list_onfo
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_89',handle2)
    CALL section_vals_val_get(force_env_section,"METHOD",i_val=method_name_id,error=error)
    CALL section_vals_val_get(subsys_section,"TOPOLOGY%DISABLE_EXCLUSION_LISTS",&
                              l_val=disable_exclusion_lists,error=error)
    IF ((method_name_id == do_fist).AND.(.NOT.disable_exclusion_lists)) THEN
       CPPostcondition(PRESENT(exclusions),cp_failure_level,routineP,error,failure)
       natom=topology%natoms
       ! allocate exclusions. Most likely they would only be needed for the local_particles
       ALLOCATE(exclusions(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I=1,natom
          NULLIFY(exclusions(i)%list_exclude_vdw)
          NULLIFY(exclusions(i)%list_exclude_ei)
          NULLIFY(exclusions(i)%list_onfo)
       ENDDO
       ! Reorder bonds
       ALLOCATE(ex_bond_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I=1,natom
          ALLOCATE(ex_bond_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       N = 0
       IF(ASSOCIATED(conn_info%bond_a)) THEN
         N = SIZE(conn_info%bond_a)
         CALL reorder_structure(ex_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)
       END IF

       ! Check if a list of 1-2 exclusion bonds is defined.. if not use all bonds
       NULLIFY(ex_bond_list_vdw, ex_bond_list_ei)
       ! VdW
       exclude_section => section_vals_get_subs_vals(topology_section,"EXCLUDE_VDW_LIST",error=error)
       CALL section_vals_get(exclude_section, explicit=explicit, error=error)
       present_12_excl_vdw_list = .FALSE.
       IF (explicit) present_12_excl_vdw_list = .TRUE.
       IF (present_12_excl_vdw_list) THEN
          ALLOCATE(ex_bond_list_vdw(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO I=1,natom
             ALLOCATE(ex_bond_list_vdw(I)%array1(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          CALL setup_exclusion_list(exclude_section,"BOND", ex_bond_list, ex_bond_list_vdw,&
                                    particle_set, error)
       ELSE
          ex_bond_list_vdw => ex_bond_list
       END IF
       ! EI
       exclude_section => section_vals_get_subs_vals(topology_section,"EXCLUDE_EI_LIST",error=error)
       CALL section_vals_get(exclude_section, explicit=explicit, error=error)
       present_12_excl_ei_list = .FALSE.
       IF (explicit) present_12_excl_ei_list = .TRUE.
       IF (present_12_excl_ei_list) THEN
          ALLOCATE(ex_bond_list_ei(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO I=1,natom
             ALLOCATE(ex_bond_list_ei(I)%array1(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          CALL setup_exclusion_list(exclude_section,"BOND", ex_bond_list, ex_bond_list_ei,&
                                    particle_set, error)
       ELSE
          ex_bond_list_ei => ex_bond_list
       END IF


       CALL section_vals_val_get(topology_section, "AUTOGEN_EXCLUDE_LISTS", &
            l_val=autogen, error=error)
       ! Reorder bends
       ALLOCATE(ex_bend_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I=1,natom
          ALLOCATE(ex_bend_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       IF (autogen) THEN
          ! Construct autogenerated 1-3 pairs, i.e. all possible 1-3 pairs instead
          ! of only the bends that are present in the topology.
          ALLOCATE(pairs(0,2),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          N = 0
          DO iatom = 1, natom
             DO i = 1, SIZE(ex_bond_list(iatom)%array1)
                ! a neighboring atom of iatom:
                atom_i = ex_bond_list(iatom)%array1(i)
                DO j = 1, i-1
                   ! another neighboring atom of iatom
                   atom_j = ex_bond_list(iatom)%array1(j)
                   ! It is only a true bend if there is no shorter path.
                   ! No need to check if i and j correspond to the same atom.
                   ! Check if i and j are not involved in a bond:
                   check = .FALSE.
                   DO counter = 1, SIZE(ex_bond_list(atom_i)%array1)
                     IF (ex_bond_list(atom_i)%array1(counter) == atom_j) THEN
                        check = .TRUE.
                        EXIT
                     END IF
                   END DO
                   IF (check) CYCLE
                   ! Add the genuine 1-3 pair
                   N = N + 1
                   IF (SIZE(pairs,dim=1)<=N) THEN
                      CALL reallocate(pairs, 1, N+5, 1, 2)
                   END IF
                   pairs(N,1) = atom_i
                   pairs(N,2) = atom_j
                END DO
             END DO
          END DO
          CALL reorder_structure(ex_bend_list, pairs(:,1), pairs(:,2), N, error)
          DEALLOCATE(pairs,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          IF (ASSOCIATED(conn_info%theta_a)) THEN
             N = SIZE(conn_info%theta_a)
             CALL reorder_structure(ex_bend_list, conn_info%theta_a, conn_info%theta_c, N, error)
          END IF
       END IF

       ! Reorder onfo
       ALLOCATE(ex_onfo_list(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I=1,natom
          ALLOCATE(ex_onfo_list(I)%array1(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       IF (autogen) THEN
          ! Construct autogenerated 1-4 pairs, i.e. all possible 1-4 pairs instead
          ! of only the onfo's that are present in the topology.
          ALLOCATE(pairs(0,2),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          N = 0
          DO iatom = 1, natom
             DO i = 1, SIZE(ex_bond_list(iatom)%array1)
                ! a neighboring atom of iatom:
                atom_i = ex_bond_list(iatom)%array1(i)
                DO j = 1, SIZE(ex_bend_list(iatom)%array1)
                   ! a next neighboring atom of iatom:
                   atom_j = ex_bend_list(iatom)%array1(j)
                   ! It is only a true onfo if there is no shorter path.
                   ! check if i and j are not the same atom
                   IF (atom_i==atom_j) CYCLE
                   ! check if i and j are not involved in a bond
                   check = .FALSE.
                   DO counter = 1, SIZE(ex_bond_list(atom_i)%array1)
                     IF (ex_bond_list(atom_i)%array1(counter) == atom_j) THEN
                        check = .TRUE.
                        EXIT
                     END IF
                   END DO
                   IF (check) CYCLE
                   ! check if i and j are not involved in a bend
                   check = .FALSE.
                   DO counter = 1, SIZE(ex_bend_list(atom_i)%array1)
                     IF (ex_bend_list(atom_i)%array1(counter) == atom_j) THEN
                        check = .TRUE.
                        EXIT
                     END IF
                   END DO
                   IF (check) CYCLE
                   ! Add the true onfo.
                   N = N + 1
                   IF (SIZE(pairs,dim=1)<=N) THEN
                      CALL reallocate(pairs, 1, N+5, 1, 2)
                   END IF
                   pairs(N,1) = atom_i
                   pairs(N,2) = atom_j
                END DO
             END DO
          END DO
          CALL reorder_structure(ex_onfo_list, pairs(:,1), pairs(:,2), N, error)
          DEALLOCATE(pairs,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          IF(ASSOCIATED(conn_info%onfo_a)) THEN
             N = SIZE(conn_info%onfo_a)
             CALL reorder_structure(ex_onfo_list, conn_info%onfo_a, conn_info%onfo_b, N, error)
          END IF
       END IF

       ! Build the exclusion (and onfo) list per atom.
       DO iatom = 1, SIZE(particle_set)
          ! Setup exclusion list for VDW: always exclude itself
          dim0 = 1
          ! exclude bond-neighbors (only if do_skip_12 .OR. do_skip_13 .OR. do_skip_14)
          dim1 = 0
          IF  (topology%exclude_vdw==do_skip_12.OR. &
               topology%exclude_vdw==do_skip_13.OR. &
               topology%exclude_vdw==do_skip_14) dim1 = SIZE(ex_bond_list_vdw(iatom)%array1)
          dim1 = dim1 + dim0
          dim2 = 0
          IF  (topology%exclude_vdw==do_skip_13.OR. &
               topology%exclude_vdw==do_skip_14) dim2 = SIZE(ex_bend_list(iatom)%array1)
          dim2 = dim1 + dim2
          dim3 = 0
          IF  (topology%exclude_vdw==do_skip_14) dim3 = SIZE(ex_onfo_list(iatom)%array1)
          dim3 = dim2 + dim3
          IF (dim3 /= 0) THEN
             NULLIFY(list, wlist)
             ALLOCATE(wlist(dim3),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             wlist(  dim0:dim0) = iatom
             IF (dim1>dim0) wlist(dim0+1:dim1) = ex_bond_list_vdw(iatom)%array1
             IF (dim2>dim1) wlist(dim1+1:dim2) = ex_bend_list(iatom)%array1
             IF (dim3>dim2) wlist(dim2+1:dim3) = ex_onfo_list(iatom)%array1
             ! Get a unique list
             DO i = 1, SIZE(wlist)-1
                IF (wlist(i)==0)CYCLE
                DO j = i+1, SIZE(wlist)
                   IF(wlist(j)==wlist(i)) wlist(j) = 0
                END DO
             END DO
             dim3 = SIZE(wlist)-COUNT(wlist==0)
             ALLOCATE(list(dim3),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             j = 0
             DO i = 1, SIZE(wlist)
                IF (wlist(i)==0)CYCLE
                j = j + 1
                list(j)=wlist(i)
             END DO
             DEALLOCATE(wlist,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ! Unique list completed
             NULLIFY(list2)
             IF((topology%exclude_vdw==topology%exclude_ei).AND.&
                  (.NOT.present_12_excl_ei_list).AND.(.NOT.present_12_excl_vdw_list)) THEN
                list2 => list
             ELSE
                ! Setup exclusion list for EI : always exclude itself
                dim0 = 1
                ! exclude bond-neighbors (only if do_skip_12 .OR. do_skip_13 .OR. do_skip_14)
                dim1 = 0
                IF  (topology%exclude_ei==do_skip_12.OR. &
                     topology%exclude_ei==do_skip_13.OR. &
                     topology%exclude_ei==do_skip_14) dim1 = SIZE(ex_bond_list_ei(iatom)%array1)
                dim1 = dim1 + dim0
                dim2 = 0
                IF  (topology%exclude_ei==do_skip_13.OR. &
                     topology%exclude_ei==do_skip_14) dim2 = SIZE(ex_bend_list(iatom)%array1)
                dim2 = dim1 + dim2
                dim3 = 0
                IF  (topology%exclude_ei==do_skip_14) dim3 = SIZE(ex_onfo_list(iatom)%array1)
                dim3 = dim2 + dim3

                IF (dim3 /= 0) THEN
                   ALLOCATE(wlist(dim3),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   wlist(  dim0:dim0) = iatom
                   IF (dim1>dim0) wlist(dim0+1:dim1) = ex_bond_list_ei(iatom)%array1
                   IF (dim2>dim1) wlist(dim1+1:dim2) = ex_bend_list(iatom)%array1
                   IF (dim3>dim2) wlist(dim2+1:dim3) = ex_onfo_list(iatom)%array1
                   ! Get a unique list
                   DO i = 1, SIZE(wlist)-1
                      IF (wlist(i)==0)CYCLE
                      DO j = i+1, SIZE(wlist)
                         IF(wlist(j)==wlist(i)) wlist(j) = 0
                      END DO
                   END DO
                   dim3 = SIZE(wlist)-COUNT(wlist==0)
                   ALLOCATE(list2(dim3),stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   j = 0
                   DO i = 1, SIZE(wlist)
                      IF (wlist(i)==0)CYCLE
                      j = j + 1
                      list2(j)=wlist(i)
                   END DO
                   DEALLOCATE(wlist,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ! Unique list completed
                END IF
             END IF
          END IF
          exclusions(iatom)%list_exclude_vdw => list
          exclusions(iatom)%list_exclude_ei  => list2
          ! Keep a list of onfo atoms for proper selection of specialized 1-4
          ! potentials instead of conventional nonbonding potentials.
          ALLOCATE(exclusions(iatom)%list_onfo(SIZE(ex_onfo_list(iatom)%array1)))
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! copy of data, not copy of pointer
          exclusions(iatom)%list_onfo = ex_onfo_list(iatom)%array1
          IF (iw>0) THEN
             IF (ASSOCIATED(list)) &
                  WRITE(iw,*)"exclusion list_vdw :: ",&
                  "atom num :",iatom,"exclusion list ::",&
                  list
             IF(topology%exclude_vdw/=topology%exclude_ei) THEN
                IF (ASSOCIATED(list2)) &
                     WRITE(iw,*)"exclusion list_ei :: ",&
                     "atom num :",iatom,"exclusion list ::",&
                     list2
             END IF
             IF (ASSOCIATED(exclusions(iatom)%list_onfo)) &
                  WRITE(iw,*)"onfo list :: ",&
                  "atom num :",iatom,"onfo list ::",&
                  exclusions(iatom)%list_onfo
          END IF
       END DO
       ! deallocate onfo
       DO I=1,natom
          DEALLOCATE(ex_onfo_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(ex_onfo_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! deallocate bends
       DO I=1,natom
          DEALLOCATE(ex_bend_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(ex_bend_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! deallocate bonds
       IF (present_12_excl_ei_list) THEN
          DO I=1,natom
             DEALLOCATE(ex_bond_list_ei(I)%array1,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          DEALLOCATE(ex_bond_list_ei,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          NULLIFY(ex_bond_list_ei)
       END IF
       IF (present_12_excl_vdw_list) THEN
          DO I=1,natom
             DEALLOCATE(ex_bond_list_vdw(I)%array1,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          DEALLOCATE(ex_bond_list_vdw,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          NULLIFY(ex_bond_list_vdw)
       END IF
       DO I=1,natom
          DEALLOCATE(ex_bond_list(I)%array1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(ex_bond_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Set the atomic_kind_set()%fist_potentail%[qeff] (PART 1)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_10',handle2)
    CALL section_vals_val_get(force_env_section,"METHOD",i_val=method_name_id,error=error)
    IF(method_name_id == do_fist) THEN
       DO i=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(i)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
          qeff = charge(i)
          NULLIFY(fist_potential)
          CALL allocate_potential(fist_potential,error)
          CALL set_potential(potential=fist_potential,qeff=qeff)
          CALL set_atomic_kind(atomic_kind=atomic_kind,fist_potential=fist_potential)
       END DO
    END IF
    DEALLOCATE(charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 12. Set the atom_list for molecule_kind in molecule_kind_set (PART 2)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//'_11',handle2)
    DO i=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(i)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            natom=natom,molecule_list=molecule_list,&
            atom_list=atom_list)
       molecule => molecule_set(molecule_list(1))
       CALL get_molecule(molecule=molecule,&
            first_atom=first,last_atom=last)
       DO j=1,natom
          DO k=1,SIZE(atomic_kind_set)
             atomic_kind => atomic_kind_set(k)
             CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
             IF(method_name_id == do_fist) THEN
                CALL get_atomic_kind(atomic_kind=atomic_kind,fist_potential=fist_potential)
                CALL get_potential(potential=fist_potential,qeff=qeff)
                IF((id2str(atom_list(j)%id_name) == atmname).AND.(qeff==atom_info%atm_charge(first+j-1)) ) THEN
                   atom_list(j)%atomic_kind => atomic_kind_set(k)
                   EXIT
                END IF
             ELSE
                IF(id2str(atom_list(j)%id_name) == atmname) THEN
                   atom_list(j)%atomic_kind => atomic_kind_set(k)
                   EXIT
                END IF
             END IF
          END DO
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,atom_list=atom_list)
    END DO
    CALL timestop(handle2)

    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_coordinate_pack

! *****************************************************************************
!> \brief Builds the exclusion list for VDW and EI if an explicit list of terms
!>        is provided by the user. Otherwise all possibilities are excluded
!> \param exclude_section ...
!> \param keyword ...
!> \param ex_bond_list ...
!> \param ex_bond_list_w ...
!> \param particle_set ...
!> \param error ...
!> \par History
!>      Teodoro Laino [tlaino] - 12.2009
! *****************************************************************************
  SUBROUTINE setup_exclusion_list(exclude_section, keyword, ex_bond_list,&
       ex_bond_list_w, particle_set, error)
    TYPE(section_vals_type), POINTER         :: exclude_section
    CHARACTER(LEN=*), INTENT(IN)             :: keyword
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: ex_bond_list, ex_bond_list_w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_exclusion_list', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: flag1, flag2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: names
    INTEGER                                  :: i, ind, j, k, l, m, n_rep
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ex_bond_list),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(ex_bond_list_w),cp_failure_level,routineP,error,failure)
    SELECT CASE(keyword)
    CASE ("BOND")
       CALL section_vals_val_get(exclude_section,keyword,n_rep_val=n_rep,error=error)
       DO j=1, SIZE(ex_bond_list)
          CPPostcondition(ASSOCIATED(ex_bond_list(j)%array1),cp_failure_level,routineP,error,failure)
          CPPostcondition(ASSOCIATED(ex_bond_list_w(j)%array1),cp_failure_level,routineP,error,failure)

          flag1=particle_set(j)%atomic_kind%name
          m = SIZE(ex_bond_list(j)%array1)
          CALL reallocate(ex_bond_list_w(j)%array1,1,m)

          l = 0
          DO k =1, m
             ind = ex_bond_list(j)%array1(k)
             flag2=particle_set(ind)%atomic_kind%name
             DO i=1, n_rep
                CALL section_vals_val_get(exclude_section,keyword,i_rep_val=i,&
                     c_vals=names,error=error)
                IF  (((TRIM(names(1))==TRIM(flag1)).AND.(TRIM(names(2))==TRIM(flag2))).OR.&
                     ((TRIM(names(1))==TRIM(flag2)).AND.(TRIM(names(2))==TRIM(flag1)))) THEN
                   l = l + 1
                   ex_bond_list_w(j)%array1(l) = ind
                END IF
             END DO
          END DO
          CALL reallocate(ex_bond_list_w(j)%array1,1,l)
       END DO
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE setup_exclusion_list

END MODULE topology_coordinate_util

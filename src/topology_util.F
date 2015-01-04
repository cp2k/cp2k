!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collection of subroutine needed for topology related things
!> \par History
!>     jgh (23-05-2004) Last atom of molecule information added
! *****************************************************************************
MODULE topology_util
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE graphcon,                        ONLY: graph_type,&
                                             hash_molecule,&
                                             reorder_graph,&
                                             vertex
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mm_mapping_library,              ONLY: amber_map,&
                                             charmm_map,&
                                             gromos_map
  USE periodic_table,                  ONLY: get_ptable_info
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             str2id
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_util'

! *****************************************************************************
  TYPE array1_list_type
     INTEGER, POINTER, DIMENSION(:) :: array1 => NULL()
  END TYPE array1_list_type

! *****************************************************************************
  TYPE array2_list_type
     INTEGER, POINTER, DIMENSION(:) :: array1 => NULL()
     INTEGER, POINTER, DIMENSION(:) :: array2 => NULL()
  END TYPE array2_list_type

  PRIVATE
  PUBLIC :: topology_set_atm_mass,&
            topology_reorder_atoms,&
            topology_molecules_check,&
            check_subsys_element,&
            reorder_structure,&
            find_molecule,&
            array1_list_type,&
            array2_list_type,&
            give_back_molecule,&
            reorder_list_array,&
            tag_molecule

  INTERFACE reorder_structure
     MODULE PROCEDURE reorder_structure1d, reorder_structure2d
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param topology ...
!> \param qmmm ...
!> \param qmmm_env_mm ...
!> \param subsys_section ...
!> \param force_env_section ...
!> \param error ...
!> \par History
!>      Teodoro Laino 09.2006 - Rewritten with a graph matching algorithm
! *****************************************************************************
  SUBROUTINE topology_reorder_atoms(topology,qmmm,qmmm_env_mm,subsys_section,force_env_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env_mm
    TYPE(section_vals_type), POINTER         :: subsys_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_reorder_atoms', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: mol_id
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: molname(:), telement(:), &
                                                tlabel_atmname(:), &
                                                tlabel_molname(:), &
                                                tlabel_resname(:)
    INTEGER :: handle, i, iatm, iindex, ikind, imol, imol_ref, iref, iund, &
      iw, j, k, location, max_mol_num, mm_index, n, n_rep, n_var, natom, &
      natom_loc, nkind, nlinks, old_hash, old_mol, output_unit, qm_index, &
      stat, unique_mol
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms, qm_atom_index
    INTEGER, POINTER :: atm_map1(:), atm_map2(:), atm_map3(:), &
      map_atom_type(:), map_mol_hash(:), mm_indexes_n(:), mm_indexes_v(:), &
      mol_bnd(:,:), mol_hash(:), mol_num(:), new_position(:), order(:), &
      tmp_n(:), tmp_v(:), wrk(:)
    LOGICAL                                  :: explicit, failure, matches, &
                                                my_qmmm
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tr
    REAL(KIND=dp), POINTER                   :: tatm_charge(:), tatm_mass(:)
    TYPE(array1_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom_bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(graph_type), DIMENSION(:), POINTER  :: reference_set
    TYPE(section_vals_type), POINTER         :: generate_section, &
                                                isolated_section, qm_kinds, &
                                                qmmm_link_section, &
                                                qmmm_section
    TYPE(vertex), DIMENSION(:), POINTER      :: reference, unordered

    failure = .FALSE.
    NULLIFY(logger,generate_section,isolated_section,tmp_v,tmp_n)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    output_unit= cp_logger_get_default_io_unit(logger)
    IF (output_unit>0) WRITE(output_unit,'(T2,"REORDER |  ")')

    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env_mm)) my_qmmm = qmmm

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    natom = topology%natoms

    NULLIFY( new_position, reference_set)
    NULLIFY(tlabel_atmname, telement, mol_num, tlabel_molname, tlabel_resname)
    NULLIFY(tr, tatm_charge, tatm_mass, atm_map1, atm_map2, atm_map3)
    ! This routine can be called only at a very high level where these structures are still
    ! not even taken into account...
    CPPostcondition(.NOT.ASSOCIATED(atom_info%map_mol_num),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(atom_info%map_mol_typ),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(atom_info%map_mol_res),cp_failure_level,routineP,error,failure)
    !ALLOCATE all the temporary arrays needed for this routine
    ALLOCATE(new_position(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mol_num(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(molname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tlabel_atmname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tlabel_molname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tlabel_resname(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tr(3,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tatm_charge(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tatm_mass(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(telement(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(atm_map1(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(atm_map2(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! The only information we have at this level is the connectivity of the system.
    ! 0. Build a list of mapping atom types
    ALLOCATE(map_atom_type(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! 1. Build a list of bonds
    ALLOCATE(atom_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       map_atom_type(I) = atom_info%id_atmname(i)
       ALLOCATE(atom_bond_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(atom_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)
    ALLOCATE(atom_info%map_mol_num(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    atom_info%map_mol_num = -1
    CALL find_molecule(atom_bond_list,atom_info%map_mol_num, atom_info%id_molname)
    max_mol_num = MAXVAL(atom_info%map_mol_num)
    ! In atom_info%map_mol_num have already been mapped molecules
    ALLOCATE(mol_bnd(2,max_mol_num),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mol_hash(max_mol_num),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_mol_hash(max_mol_num),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! 2. Sort the map_mol_num array.. atm_map1 contains now the mapped index
    !    of the reordered array
    CALL sort(atom_info%map_mol_num, natom, atm_map1)
    old_mol=0
    iindex =0
    imol=0
    DO i=1,natom
       IF (old_mol.NE.atom_info%map_mol_num(I)) THEN
          old_mol=atom_info%map_mol_num(I)
          iindex=0
          IF (imol>0) THEN
             mol_bnd(2,imol) = i-1
          END IF
          imol =imol+1
          mol_bnd(1,imol) = i
       ENDIF
       iindex=iindex+1
       atm_map2(atm_map1(i))=iindex
    ENDDO
    mol_bnd(2,imol) = natom
    ! Indexes of the two molecules to check
    iref = 1
    iund = max_mol_num/2+1
    ! Allocate reference and unordered
    NULLIFY(reference, unordered)
    ! This is the real matching of graphs
    DO j = 1, max_mol_num
       CALL setup_graph(j, reference,  map_atom_type,&
            atom_bond_list, mol_bnd, atm_map1, atm_map2, error=error)

       ALLOCATE(order(SIZE(reference)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL hash_molecule(reference, order, mol_hash(j))

       DEALLOCATE(order,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO I=1,SIZE(reference)
          DEALLOCATE(reference(I)%bonds,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(reference, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    ! Reorder molecules hashes
    CALL sort(mol_hash, max_mol_num, map_mol_hash)
    ! Now find unique molecules and reorder atoms too (if molecules match)
    old_hash=-1
    unique_mol=0
    natom_loc=0
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,"REORDER |  ",A)')&
            "Reordering Molecules. The Reordering of molecules will override all",&
            "information regarding molecule names and residue names.",&
            "New ones will be provided on the basis of the connectivity!"
    END IF
    DO j = 1, max_mol_num
       IF (mol_hash(j).NE.old_hash) THEN
          unique_mol=unique_mol+1
          old_hash=mol_hash(j)
          CALL setup_graph_set(reference_set,unique_mol,map_mol_hash(j),&
               map_atom_type,atom_bond_list, mol_bnd, atm_map1, atm_map2,&
               atm_map3, error)
          ! Reorder Last added reference
          mol_id = TRIM(ADJUSTL(cp_to_string(unique_mol)))
          DO i=1,SIZE(atm_map3)
             natom_loc=natom_loc+1
             new_position(natom_loc) = atm_map3(i)
             molname(natom_loc)      = mol_id
             mol_num(natom_loc)      = unique_mol
          END DO
          DEALLOCATE(atm_map3,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          CALL setup_graph(map_mol_hash(j), unordered,  map_atom_type,&
                      atom_bond_list, mol_bnd, atm_map1, atm_map2, atm_map3,&
                      error)
          DO imol_ref=1,unique_mol
             !
             reference => reference_set(imol_ref)%graph
             ALLOCATE(order(SIZE(reference)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL reorder_graph(reference, unordered, order, matches)
             IF (matches) EXIT
             DEALLOCATE(order,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          IF (matches) THEN
             ! Reorder according the correct index
             ALLOCATE(wrk(SIZE(order)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL sort(order, SIZE(order), wrk)
             DO i=1,SIZE(order)
                natom_loc=natom_loc+1
                new_position(natom_loc) = atm_map3(wrk(i))
                molname(natom_loc)      = mol_id
                mol_num(natom_loc)      = unique_mol
             END DO
             !
             DEALLOCATE(order,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(wrk,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
             unique_mol=unique_mol+1
             CALL setup_graph_set(reference_set,unique_mol,map_mol_hash(j),&
                  map_atom_type,atom_bond_list, mol_bnd, atm_map1, atm_map2,&
                  atm_map3, error)
             ! Reorder Last added reference
             mol_id = TRIM(ADJUSTL(cp_to_string(unique_mol)))
             DO i=1,SIZE(atm_map3)
                natom_loc=natom_loc+1
                new_position(natom_loc) = atm_map3(i)
                molname(natom_loc)      = mol_id
                mol_num(natom_loc)      = unique_mol
             END DO
             DEALLOCATE(atm_map3,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          DO I=1,SIZE(unordered)
             DEALLOCATE(unordered(I)%bonds,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          DEALLOCATE(unordered, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(atm_map3, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDDO
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,"REORDER |  ",A,I7,A)')"Number of unique molecules found:",unique_mol,"."
    END IF
    CPPostcondition(natom_loc==natom,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_atom_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(atm_map1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(atm_map2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mol_bnd,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mol_hash,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_mol_hash,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Deallocate working arrays
    DO I=1,natom
       DEALLOCATE(atom_bond_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(atom_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(atom_info%map_mol_num,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Deallocate reference_set
    DO i = 1,SIZE(reference_set)
       DO j = 1, SIZE(reference_set(i)%graph)
          DEALLOCATE(reference_set(i)%graph(j)%bonds,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(reference_set(i)%graph,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(reference_set,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !Lets swap the atoms now
    DO iatm=1,natom
       location = new_position(iatm)
       tlabel_molname(iatm)       = id2str(atom_info%id_molname(location))
       tlabel_resname(iatm)       = id2str(atom_info%id_resname(location))
       tlabel_atmname(iatm)       = id2str(atom_info%id_atmname(location))
       telement(iatm)             = id2str(atom_info%id_element(location))
       tr(1,iatm)                 = atom_info%r(1,location)
       tr(2,iatm)                 = atom_info%r(2,location)
       tr(3,iatm)                 = atom_info%r(3,location)
       tatm_charge(iatm)          = atom_info%atm_charge(location)
       tatm_mass(iatm)            = atom_info%atm_mass(location)
    END DO
    IF (topology%create_molecules) THEN
       DO iatm=1,natom
          tlabel_molname(iatm)   = "MOL"//TRIM(molname(iatm))
          tlabel_resname(iatm)   = "R"//TRIM(molname(iatm))
       END DO
       topology%create_molecules = .FALSE.
    END IF
    DO iatm=1,natom
       atom_info%id_molname(iatm)   = str2id(tlabel_molname(iatm))
       atom_info%id_resname(iatm)   = str2id(tlabel_resname(iatm))
       atom_info%id_atmname(iatm)   = str2id(tlabel_atmname(iatm))
       atom_info%id_element(iatm)   = str2id(telement(iatm))
       atom_info%resid(iatm)        = mol_num(iatm)
       atom_info%r(1,iatm)          = tr(1,iatm)
       atom_info%r(2,iatm)          = tr(2,iatm)
       atom_info%r(3,iatm)          = tr(3,iatm)
       atom_info%atm_charge(iatm)   = tatm_charge(iatm)
       atom_info%atm_mass(iatm)     = tatm_mass(iatm)
    END DO

    ! Let's reorder all the list provided in the input..
    ALLOCATE(wrk(SIZE(new_position)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL sort(new_position, SIZE(new_position), wrk)

    ! NOTE: In the future the whole list of possible integers should be updated at this level..
    ! Let's update all the integer lists in the qmmm_env_mm section and in the input sections
    !       from where qmmm_env_qm will be read.
    IF (my_qmmm) THEN
       ! Update the qmmm_env_mm
       CPPrecondition(SIZE(qmmm_env_mm%qm_atom_index)/=0,cp_failure_level,routineP,error,failure)
       CPPrecondition(ALL(qmmm_env_mm%qm_atom_index /= 0),cp_failure_level,routineP,error,failure)
       ALLOCATE(qm_atom_index(SIZE(qmmm_env_mm%qm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iatm=1,SIZE(qmmm_env_mm%qm_atom_index)
          qm_atom_index(iatm) = wrk(qmmm_env_mm%qm_atom_index(iatm))
       END DO
       qmmm_env_mm%qm_atom_index = qm_atom_index
       CPPostcondition(ALL(qmmm_env_mm%qm_atom_index /= 0),cp_failure_level,routineP,error,failure)
       DEALLOCATE(qm_atom_index,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Update the qmmm_section: MM_INDEX of the QM_KIND
       qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)
       qm_kinds => section_vals_get_subs_vals(qmmm_section,"QM_KIND",error=error)
       CALL section_vals_get(qm_kinds,n_repetition=nkind,error=error)
       DO ikind=1,nkind
          CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,n_rep_val=n_var,&
               error=error)
          DO k = 1, n_var
             CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,i_rep_val=k,&
                  i_vals=mm_indexes_v,error=error)
             ALLOCATE(mm_indexes_n(SIZE(mm_indexes_v)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO j = 1, SIZE(mm_indexes_v)
                mm_indexes_n(j) = wrk(mm_indexes_v(j))
             END DO
             CALL section_vals_val_set(qm_kinds,"MM_INDEX", i_rep_section=ikind, &
                                       i_vals_ptr=mm_indexes_n, i_rep_val=k, error=error)
          END DO
       END DO
       ! Handle the link atoms
       IF (qmmm_env_mm%qmmm_link) THEN
          ! Update the qmmm_env_mm
          CPPrecondition(SIZE(qmmm_env_mm%mm_link_atoms)>0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mm_link_atoms(SIZE(qmmm_env_mm%mm_link_atoms)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO iatm=1,SIZE(qmmm_env_mm%mm_link_atoms)
             mm_link_atoms(iatm) = wrk(qmmm_env_mm%mm_link_atoms(iatm))
          END DO
          qmmm_env_mm%mm_link_atoms = mm_link_atoms
          CPPostcondition(ALL(qmmm_env_mm%mm_link_atoms /= 0),cp_failure_level,routineP,error,failure)
          DEALLOCATE(mm_link_atoms,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Update the qmmm_section: MM_INDEX,QM_INDEX of the LINK atom list
          qmmm_link_section => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
          CALL section_vals_get(qmmm_link_section,n_repetition=nlinks,error=error)
          CPPostcondition(nlinks /= 0,cp_failure_level,routineP,error,failure)
          DO ikind= 1, nlinks
             CALL section_vals_val_get(qmmm_link_section,"QM_INDEX",i_rep_section=ikind,i_val=qm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,error=error)
             mm_index = wrk(mm_index)
             qm_index = wrk(qm_index)
             CALL section_vals_val_set(qmmm_link_section,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,error=error)
             CALL section_vals_val_set(qmmm_link_section,"QM_INDEX",i_rep_section=ikind,i_val=qm_index,error=error)
          END DO
       END IF
    END IF
    !
    !LIST of ISOLATED atoms
    !
    generate_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%GENERATE",error=error)
    isolated_section => section_vals_get_subs_vals(generate_section,"ISOLATED_ATOMS",error=error)
    CALL section_vals_get(isolated_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(isolated_section,"LIST", n_rep_val=n_rep, error=error)
       DO i = 1, n_rep
          CALL section_vals_val_get(isolated_section,"LIST", i_vals=tmp_v, i_rep_val=i, error=error)
          ALLOCATE(tmp_n(SIZE(tmp_v)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j = 1, SIZE(tmp_v)
             tmp_n(j) = wrk(tmp_v(j))
          END DO
          CALL section_vals_val_set(isolated_section,"LIST", i_vals_ptr=tmp_n, i_rep_val=i, error=error)
       END DO
    END IF
    DEALLOCATE(wrk,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !DEALLOCATE all the temporary arrays needed for this routine
    DEALLOCATE(new_position,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tlabel_atmname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tlabel_molname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tlabel_resname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(telement,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tr,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tatm_charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tatm_mass,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(molname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mol_num,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! DEALLOCATE the bond structures in the connectivity
    DEALLOCATE(conn_info%bond_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conn_info%bond_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (output_unit>0) WRITE(output_unit,'(T2,"REORDER |  ")')
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_reorder_atoms

! *****************************************************************************
!> \brief Set up a SET of graph kind
!> \param graph_set ...
!> \param idim ...
!> \param ind ...
!> \param array2 ...
!> \param atom_bond_list ...
!> \param map_mol ...
!> \param atm_map1 ...
!> \param atm_map2 ...
!> \param atm_map3 ...
!> \param error ...
!> \author Teodoro Laino 10.2006
! *****************************************************************************
  SUBROUTINE setup_graph_set(graph_set, idim, ind, array2, atom_bond_list, map_mol,&
       atm_map1, atm_map2, atm_map3, error)
    TYPE(graph_type), DIMENSION(:), POINTER  :: graph_set
    INTEGER, INTENT(IN)                      :: idim, ind
    INTEGER, DIMENSION(:), POINTER           :: array2
    TYPE(array1_list_type), DIMENSION(:), &
      INTENT(IN)                             :: atom_bond_list
    INTEGER, DIMENSION(:, :), POINTER        :: map_mol
    INTEGER, DIMENSION(:), POINTER           :: atm_map1, atm_map2, atm_map3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_graph_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ldim
    LOGICAL                                  :: failure
    TYPE(graph_type), DIMENSION(:), POINTER  :: tmp_graph_set

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ldim = 0
       NULLIFY(tmp_graph_set)
       IF (ASSOCIATED(graph_set)) THEN
          ldim = SIZE(graph_set)
          CPPostcondition(ldim+1==idim,cp_failure_level,routineP,error,failure)
          NULLIFY(tmp_graph_set)
          CALL allocate_graph_set(graph_set, tmp_graph_set, error=error)
       END IF
       CALL allocate_graph_set(tmp_graph_set, graph_set, ldim, ldim+1, error)
       CALL setup_graph(ind, graph_set(ldim+1)%graph, array2,&
                      atom_bond_list, map_mol, atm_map1, atm_map2, atm_map3,&
                      error=error)
    END IF

  END SUBROUTINE setup_graph_set

! *****************************************************************************
!> \brief Allocate a new graph_set deallocating an old one..
!> \param graph_set_in ...
!> \param graph_set_out ...
!> \param ldim_in ...
!> \param ldim_out ...
!> \param error ...
!> \author Teodoro Laino 10.2006
! *****************************************************************************
  SUBROUTINE allocate_graph_set(graph_set_in, graph_set_out, ldim_in, ldim_out,  error)
    TYPE(graph_type), DIMENSION(:), POINTER  :: graph_set_in, graph_set_out
    INTEGER, INTENT(IN), OPTIONAL            :: ldim_in, ldim_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_graph_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b_dim, i, j, mydim_in, &
                                                mydim_out, stat, v_dim
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(graph_set_out),cp_failure_level,routineP,error,failure)
       mydim_in  = 0
       mydim_out = 0
       IF (ASSOCIATED(graph_set_in)) THEN
          mydim_in  = SIZE(graph_set_in)
          mydim_out = SIZE(graph_set_in)
       END IF
       IF (PRESENT(ldim_in))  mydim_in  = ldim_in
       IF (PRESENT(ldim_out)) mydim_out = ldim_out
       ALLOCATE(graph_set_out(mydim_out),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, mydim_out
          NULLIFY(graph_set_out(i)%graph)
       END DO
       ! Copy graph structure into the temporary array
       DO i = 1, mydim_in
          v_dim = SIZE(graph_set_in(i)%graph)
          ALLOCATE(graph_set_out(i)%graph(v_dim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j = 1, v_dim
             graph_set_out(i)%graph(j)%kind = graph_set_in(i)%graph(j)%kind
             b_dim = SIZE(graph_set_in(i)%graph(j)%bonds)
             ALLOCATE(graph_set_out(i)%graph(j)%bonds(b_dim),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             graph_set_out(i)%graph(j)%bonds = graph_set_in(i)%graph(j)%bonds
             DEALLOCATE(graph_set_in(i)%graph(j)%bonds,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(graph_set_in(i)%graph,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       IF (ASSOCIATED(graph_set_in)) THEN
          DEALLOCATE(graph_set_in,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE allocate_graph_set

! *****************************************************************************
!> \brief Set up a graph kind
!> \param ind ...
!> \param graph ...
!> \param array2 ...
!> \param atom_bond_list ...
!> \param map_mol ...
!> \param atm_map1 ...
!> \param atm_map2 ...
!> \param atm_map3 ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE setup_graph(ind, graph, array2, atom_bond_list, map_mol,&
       atm_map1, atm_map2, atm_map3, error)
    INTEGER, INTENT(IN)                      :: ind
    TYPE(vertex), DIMENSION(:), POINTER      :: graph
    INTEGER, DIMENSION(:), POINTER           :: array2
    TYPE(array1_list_type), DIMENSION(:), &
      INTENT(IN)                             :: atom_bond_list
    INTEGER, DIMENSION(:, :), POINTER        :: map_mol
    INTEGER, DIMENSION(:), POINTER           :: atm_map1, atm_map2
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: atm_map3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_graph', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idim, ifirst, ilast, j, &
                                                nbonds, nelement, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (PRESENT(atm_map3)) THEN
       CPPostcondition(.NOT.ASSOCIATED(atm_map3),cp_failure_level,routineP,error,failure)
    END IF
    CPPostcondition(.NOT.ASSOCIATED(graph),cp_failure_level,routineP,error,failure)
    ! Setup reference graph
    idim = 0
    ifirst = map_mol(1,ind)
    ilast  = map_mol(2,ind)
    nelement = ilast-ifirst+1
    ALLOCATE(graph(nelement),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(atm_map3)) THEN
       ALLOCATE(atm_map3(nelement),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DO i = ifirst, ilast
       idim = idim + 1
       graph(idim)%kind = array2(atm_map1(i))
       nbonds = SIZE(atom_bond_list(atm_map1(i))%array1)
       ALLOCATE(graph(idim)%bonds(nbonds),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO j = 1, nbonds
          graph(idim)%bonds(j) = atm_map2(atom_bond_list(atm_map1(i))%array1(j))
       END DO
       IF (PRESENT(atm_map3)) THEN
          atm_map3(idim) = atm_map1(i)
       END IF
    END DO

  END SUBROUTINE setup_graph

! *****************************************************************************
!> \brief Order arrays of lists
!> \param Ilist1 ...
!> \param Ilist2 ...
!> \param Ilist3 ...
!> \param Ilist4 ...
!> \param nsize ...
!> \param ndim ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE reorder_list_array(Ilist1, Ilist2, Ilist3, Ilist4, nsize, ndim, error)
    INTEGER, DIMENSION(:), POINTER           :: Ilist1
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Ilist2, Ilist3, Ilist4
    INTEGER, INTENT(IN)                      :: nsize, ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorder_list_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, istart, ldim, stat
    INTEGER, DIMENSION(:), POINTER           :: tmp, wrk
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(nsize>0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wrk(Ndim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL sort(Ilist1,Ndim,wrk)
    IF (nsize/=1) THEN
       ALLOCATE(tmp(Ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp = Ilist2(1:Ndim)
       DO i = 1,Ndim
          Ilist2(i)=tmp(wrk(i))
       END DO
       SELECT CASE(nsize)
       CASE(3)
          tmp = Ilist3(1:Ndim)
          DO i = 1,Ndim
             Ilist3(i)=tmp(wrk(i))
          END DO
       CASE(4)
          tmp = Ilist3(1:Ndim)
          DO i = 1,Ndim
             Ilist3(i)=tmp(wrk(i))
          END DO
          tmp = Ilist4(1:Ndim)
          DO i = 1,Ndim
             Ilist4(i)=tmp(wrk(i))
          END DO
       END SELECT
       DEALLOCATE(tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       istart = 1
       DO i = 1,Ndim
          IF (Ilist1(i)/=Ilist1(istart)) THEN
             iend = i-1
             ldim = iend-istart+1
             CALL  reorder_list_array_low(Ilist2, Ilist3, Ilist4, nsize,&
                  ldim, istart, iend, error)
             istart = i
          END IF
       END DO
       ! Last term to sort
       iend = Ndim
       ldim = iend-istart+1
       CALL  reorder_list_array_low(Ilist2, Ilist3, Ilist4, nsize,&
            ldim, istart, iend, error)
    END IF
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE reorder_list_array

! *****************************************************************************
!> \brief Low level routine for ordering arrays of lists
!> \param Ilist2 ...
!> \param Ilist3 ...
!> \param Ilist4 ...
!> \param nsize ...
!> \param ldim ...
!> \param istart ...
!> \param iend ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE reorder_list_array_low(Ilist2, Ilist3, Ilist4, nsize,&
       ldim, istart, iend, error)
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: Ilist2, Ilist3, Ilist4
    INTEGER, INTENT(IN)                      :: nsize, ldim, istart, iend
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorder_list_array_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, DIMENSION(:), POINTER           :: tmp_2, tmp_3, tmp_4
    LOGICAL                                  :: failure

    failure = .FALSE.
    SELECT CASE(nsize)
    CASE(2)
       ALLOCATE(tmp_2(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_2(:) = Ilist2(istart:iend)
       CALL reorder_list_array(tmp_2,nsize=nsize-1,ndim=ldim,&
            error=error)
       Ilist2(istart:iend) = tmp_2(:)
       DEALLOCATE(tmp_2,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE(3)
       ALLOCATE(tmp_2(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_3(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_2(:) = Ilist2(istart:iend)
       tmp_3(:) = Ilist3(istart:iend)
       CALL reorder_list_array(tmp_2,tmp_3,nsize=nsize-1,ndim=ldim,&
            error=error)
       Ilist2(istart:iend) = tmp_2(:)
       Ilist3(istart:iend) = tmp_3(:)
       DEALLOCATE(tmp_2,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_3,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE(4)
       ALLOCATE(tmp_2(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_3(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_4(ldim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_2(:) = Ilist2(istart:iend)
       tmp_3(:) = Ilist3(istart:iend)
       tmp_4(:) = Ilist4(istart:iend)
       CALL reorder_list_array(tmp_2, tmp_3, tmp_4, nsize=nsize-1,ndim=ldim,&
            error=error)
       Ilist2(istart:iend) = tmp_2(:)
       Ilist3(istart:iend) = tmp_3(:)
       Ilist4(istart:iend) = tmp_4(:)
       DEALLOCATE(tmp_2,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_3,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_4,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE reorder_list_array_low

! *****************************************************************************
!> \brief ...
!> \param icheck ...
!> \param bond_list ...
!> \param i ...
!> \param mol_natom ...
!> \param mol_map ...
!> \param my_mol ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE give_back_molecule(icheck, bond_list,i,mol_natom,mol_map,my_mol)
    LOGICAL, DIMENSION(:), POINTER           :: icheck
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    INTEGER, INTENT(IN)                      :: i
    INTEGER, INTENT(INOUT)                   :: mol_natom
    INTEGER, DIMENSION(:), POINTER           :: mol_map
    INTEGER, INTENT(IN)                      :: my_mol

    INTEGER                                  :: j, k

    IF (mol_map(i)==my_mol) THEN
       icheck(i) = .TRUE.
       DO j = 1, SIZE(bond_list(i)%array1)
          k = bond_list(i)%array1(j)
          IF (icheck(k)) CYCLE
          mol_natom = mol_natom + 1
          CALL give_back_molecule(icheck, bond_list, k, mol_natom, mol_map, my_mol)
       END DO
    ELSE
       ! Do nothing means only that bonds were found between molecules
       ! as we defined them.. This could easily be a bond detected but not
       ! physically present.. so just skip this part and go on..
    END IF
  END SUBROUTINE give_back_molecule

! *****************************************************************************
!> \brief gives back a mapping of molecules.. icheck needs to be initialized with -1
!> \param icheck ...
!> \param bond_list ...
!> \param i ...
!> \param my_mol ...
!> \author Teodoro Laino 04.2007 - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE tag_molecule(icheck, bond_list, i, my_mol)
    INTEGER, DIMENSION(:), POINTER           :: icheck
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: bond_list
    INTEGER, INTENT(IN)                      :: i, my_mol

    INTEGER                                  :: j, k

    icheck(i) = my_mol
    DO j = 1, SIZE(bond_list(i)%array1)
       k = bond_list(i)%array1(j)
       IF (k<=i) CYCLE
       CALL tag_molecule(icheck, bond_list, k, my_mol)
    END DO

  END SUBROUTINE tag_molecule

! *****************************************************************************
!> \brief Given two lists of corresponding element a complex type is built in
!>      which each element of the type has a list of mapping elements
!> \param work ...
!> \param list1 ...
!> \param list2 ...
!> \param N ...
!> \param error ...
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE reorder_structure1d(work, list1, list2, N, error)
    TYPE(array1_list_type), DIMENSION(:), &
      INTENT(INOUT)                          :: work
    INTEGER, DIMENSION(:), INTENT(IN)        :: list1, list2
    INTEGER, INTENT(IN)                      :: N
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorder_structure1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, index1, index2, Nsize, stat
    INTEGER, DIMENSION(:), POINTER           :: wrk_tmp
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO I=1,N
       index1=list1(I)
       index2=list2(I)

       wrk_tmp => work(index1)%array1
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index1)%array1(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index1)%array1(1:Nsize) = wrk_tmp
       work(index1)%array1(Nsize+1) = index2
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       wrk_tmp => work(index2)%array1
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index2)%array1(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index2)%array1(1:Nsize) = wrk_tmp
       work(index2)%array1(Nsize+1) = index1
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO

  END SUBROUTINE reorder_structure1d

! *****************************************************************************
!> \brief Given two lists of corresponding element a complex type is built in
!>      which each element of the type has a list of mapping elements
!> \param work ...
!> \param list1 ...
!> \param list2 ...
!> \param list3 ...
!> \param N ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE reorder_structure2d(work, list1, list2, list3, N, error)
    TYPE(array2_list_type), DIMENSION(:), &
      INTENT(INOUT)                          :: work
    INTEGER, DIMENSION(:), INTENT(IN)        :: list1, list2, list3
    INTEGER, INTENT(IN)                      :: N
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorder_structure2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, index1, index2, index3, &
                                                Nsize, stat
    INTEGER, DIMENSION(:), POINTER           :: wrk_tmp
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO I=1,N
       index1=list1(I)
       index2=list2(I)
       index3=list3(I)

       wrk_tmp => work(index1)%array1
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index1)%array1(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index1)%array1(1:Nsize) = wrk_tmp
       work(index1)%array1(Nsize+1) = index2
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       wrk_tmp => work(index2)%array1
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index2)%array1(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index2)%array1(1:Nsize) = wrk_tmp
       work(index2)%array1(Nsize+1) = index1
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       wrk_tmp => work(index1)%array2
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index1)%array2(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index1)%array2(1:Nsize) = wrk_tmp
       work(index1)%array2(Nsize+1) = index3
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       wrk_tmp => work(index2)%array2
       Nsize=SIZE(wrk_tmp)
       ALLOCATE(work(index2)%array2(Nsize+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work(index2)%array2(1:Nsize) = wrk_tmp
       work(index2)%array2(Nsize+1) = -index3
       DEALLOCATE(wrk_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO

  END SUBROUTINE reorder_structure2d

! *****************************************************************************
!> \brief each atom will be assigned a molecule number based on bonded fragments
!>      The array mol_info should be initialized with -1 before calling the
!>      find_molecule routine
!> \param atom_bond_list ...
!> \param mol_info ...
!> \param mol_name ...
!> \author Joost 05.2006
! *****************************************************************************
  SUBROUTINE find_molecule(atom_bond_list,mol_info,mol_name)
    TYPE(array1_list_type), DIMENSION(:), &
      INTENT(IN)                             :: atom_bond_list
    INTEGER, DIMENSION(:), POINTER           :: mol_info, mol_name

    INTEGER                                  :: I, my_mol_name, N, nmol

    N=SIZE(atom_bond_list)
    nmol=0
    DO I=1,N
       IF (mol_info(I)==-1) THEN
          nmol=nmol+1
          my_mol_name = mol_name(I)
          CALL spread_mol(atom_bond_list,mol_info,i,nmol, my_mol_name,&
               mol_name)
       ENDIF
    ENDDO
  END SUBROUTINE find_molecule

! *****************************************************************************
!> \brief spreads the molnumber over the bonded list
!> \param atom_bond_list ...
!> \param mol_info ...
!> \param iatom ...
!> \param imol ...
!> \param my_mol_name ...
!> \param mol_name ...
!> \author Joost 05.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE spread_mol(atom_bond_list,mol_info,iatom,imol, &
       my_mol_name, mol_name)
    TYPE(array1_list_type), DIMENSION(:), &
      INTENT(IN)                             :: atom_bond_list
    INTEGER, DIMENSION(:), POINTER           :: mol_info
    INTEGER, INTENT(IN)                      :: iatom, imol, my_mol_name
    INTEGER, DIMENSION(:), POINTER           :: mol_name

    INTEGER                                  :: atom_b, i

    mol_info(iatom)=imol
    DO I=1,SIZE(atom_bond_list(iatom)%array1)
       atom_b=atom_bond_list(iatom)%array1(I)
       ! In this way we're really sure that all atoms belong to the same
       ! molecule. This should correct possible errors in the generation of
       ! the bond list..
       IF (mol_info(atom_b)==-1.AND.my_mol_name==mol_name(atom_b)) THEN
          CALL spread_mol(atom_bond_list,mol_info,atom_b,imol,my_mol_name, mol_name)
          IF (mol_info(atom_b)/=imol) STOP "internal error"
       END IF
    ENDDO
  END SUBROUTINE spread_mol

! *****************************************************************************
!> \brief Use info from periodic table and set atm_mass
!> \param topology ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE topology_set_atm_mass(topology,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_set_atm_mass', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: upper_sym_1
    CHARACTER(LEN=default_string_length)     :: atmname_upper
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: keyword
    INTEGER                                  :: handle, i, i_rep, iatom, &
                                                ielem_found, iw, n_rep, &
                                                natom, stat
    LOGICAL                                  :: failure, user_defined
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mass
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: kind_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    atom_info => topology%atom_info
    natom = topology%natoms

    ! Available external info
    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    CALL section_vals_get(kind_section,n_repetition=n_rep,error=error)
    ALLOCATE(keyword(n_rep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(mass(n_rep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    mass = HUGE(0.0_dp)
    DO i_rep=1,n_rep
       CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
            c_val=keyword(i_rep),i_rep_section=i_rep,error=error)
       CALL uppercase(keyword(i_rep))
       CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
            keyword_name="MASS",n_rep_val=i,error=error)
       IF (i>0)  CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
            keyword_name="MASS",r_val=mass(i_rep),error=error)
    END DO
    !
    DO iatom=1,natom
       !If we reach this point then we've definitely identified the element..
       !Let's look if an external mass has been defined..
       user_defined = .FALSE.
       DO i= 1, SIZE(keyword)
          atmname_upper = id2str(atom_info%id_atmname(iatom))
          CALL uppercase(atmname_upper)
          IF (TRIM(atmname_upper)==TRIM(keyword(i)) .AND. mass(i)/=HUGE(0.0_dp)) THEN
             atom_info%atm_mass(iatom) = mass(i)
             user_defined = .TRUE.
             EXIT
          END IF
       END DO
       ! If name didn't match let's try with the element
       IF (.NOT.user_defined) THEN
          upper_sym_1=id2str(atom_info%id_element(iatom))
          CALL get_ptable_info(symbol=upper_sym_1,ielement=ielem_found,amass=atom_info%atm_mass(iatom))
       END IF
       IF(iw>0) WRITE(iw,'(7X,A,A5,A,F12.5)') "In topology_set_atm_mass :: element = ",&
            id2str(atom_info%id_element(iatom))," a_mass ", atom_info%atm_mass(iatom)
    END DO
    DEALLOCATE(keyword,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mass,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_set_atm_mass

! *****************************************************************************
!> \brief Check and verify that all molecules of the same kind are bonded the same
!> \param topology ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE topology_molecules_check(topology,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_molecules_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, first, first_loc, &
                                                handle, i, iatom, iw, k, &
                                                loc_counter, mol_num, &
                                                mol_typ, n, natom, stat
    LOGICAL                                  :: failure, icheck_num, &
                                                icheck_typ
    TYPE(array1_list_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom_bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    natom = topology%natoms

    IF (iw>0) WRITE(iw,'(A)')"Start of Molecule_Check",&
         "  Checking consistency between the generated molecules"

    ALLOCATE(atom_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(atom_bond_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(atom_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)

    mol_typ = atom_info%map_mol_typ(1)
    mol_num = atom_info%map_mol_num(1)
    counter     = 1
    loc_counter = 1
    first       = 1
    first_loc   = 1
    DO iatom = 2, natom
       icheck_num = (atom_info%map_mol_num(iatom) == mol_num)
       icheck_typ = (atom_info%map_mol_typ(iatom) == mol_typ)
       IF ((icheck_typ.AND.(.NOT.icheck_num)).OR.(.NOT.icheck_typ)) THEN
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 1. Check each molecule have the same number of atoms
          !-----------------------------------------------------------------------------
          IF (counter /= loc_counter) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "different number of atoms for same molecule kind"//&
                  " molecule type  = "//cp_to_string(mol_typ)//&
                  " molecule number= "//cp_to_string(mol_num)//&
                  " expected number of atoms="//cp_to_string(counter)//" found="//&
                  cp_to_string(loc_counter))
          END IF
       END IF
       IF (.NOT.icheck_typ) THEN
          first       = iatom
          first_loc   = iatom
          counter     = 1
          loc_counter = 1
          mol_typ     = atom_info%map_mol_typ(iatom)
       END IF
       IF (icheck_num) THEN
          IF (icheck_typ) loc_counter = loc_counter + 1
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 2. Check that each molecule has the same atom sequences
          !-----------------------------------------------------------------------------
          IF( atom_info%id_atmname(iatom) /= &
              atom_info%id_atmname(first+loc_counter-1) ) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "different atom name for same molecule kind"//&
                  " atom number    = "//cp_to_string(iatom)//&
                  " molecule type  = "//cp_to_string(mol_typ)//&
                  " molecule number= "//cp_to_string(mol_num)//&
                  " expected atom name="//TRIM(id2str(atom_info%id_atmname(first+loc_counter-1)))//&
                  " found="//TRIM(id2str(atom_info%id_atmname(iatom))))
          END IF
          !-----------------------------------------------------------------------------
          !-----------------------------------------------------------------------------
          ! 3. Check that each molecule have the same bond sequences
          !-----------------------------------------------------------------------------
          IF (SIZE(atom_bond_list(iatom)%array1)/=SIZE(atom_bond_list(first+loc_counter-1)%array1))&
               failure = .TRUE.
          IF (.NOT.failure) THEN
             DO k = 1,SIZE(atom_bond_list(iatom)%array1)
                IF (ALL(atom_bond_list(first+loc_counter-1)%array1-first/=&
                        atom_bond_list(iatom)%array1(k)-first_loc)) THEN
                   CALL stop_program(routineN,moduleN,__LINE__,&
                        "different sequence of bonds for same molecule kind"//&
                        " molecule type  = "//cp_to_string(mol_typ)//&
                        " molecule number= "//cp_to_string(mol_num)//&
                        " NOT FOUND! Check the connectivity of your system.")
                END IF
             END DO
          ELSE
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "different number of bonds for same molecule kind"//&
                  " molecule type  = "//cp_to_string(mol_typ)//&
                  " molecule number= "//cp_to_string(mol_num)//&
                  " expected bonds="//&
                  cp_to_string(SIZE(atom_bond_list(first+loc_counter-1)%array1))//" - "//&
                  cp_to_string(SIZE(atom_bond_list(iatom)%array1))//&
                  " NOT FOUND! Check the connectivity of your system.")
          END IF
       ELSE
          mol_num = atom_info%map_mol_num(iatom)
          loc_counter = 1
          first_loc   = iatom
       END IF
       IF (mol_num==1.AND.icheck_typ) counter = counter + 1
    END DO
    IF (iw>0) WRITE(iw,'(A)')"End of Molecule_Check"

    DO I=1,natom
       DEALLOCATE(atom_bond_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(atom_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)

  END SUBROUTINE topology_molecules_check

! *****************************************************************************
!> \brief Check and returns the ELEMENT label
!> \param element_in ...
!> \param atom_name_in ...
!> \param element_out ...
!> \param subsys_section ...
!> \param use_mm_map_first ...
!> \param error ...
!> \par History
!>      12.2005 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE check_subsys_element(element_in, atom_name_in, element_out, subsys_section, use_mm_map_first, error)
    CHARACTER(len=*), INTENT(IN)             :: element_in, atom_name_in
    CHARACTER(len=default_string_length), &
      INTENT(OUT)                            :: element_out
    TYPE(section_vals_type), POINTER         :: subsys_section
    LOGICAL                                  :: use_mm_map_first
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'check_subsys_element', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: atom_name, element_symbol, &
                                                keyword
    INTEGER                                  :: i, i_rep, n_rep
    LOGICAL                                  :: defined_kind_section, &
                                                failure, found, in_ptable
    TYPE(section_vals_type), POINTER         :: kind_section

    failure = .FALSE.
    found          = .FALSE.
    element_symbol = element_in
    atom_name      = atom_name_in
    element_out    = ""
    defined_kind_section = .FALSE.

    ! First check if a KIND section is overriding the element
    ! definition
    CALL uppercase(atom_name)
    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    CALL section_vals_get(kind_section,n_repetition=n_rep,error=error)
    DO i_rep=1,n_rep
       CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
            c_val=keyword,i_rep_section=i_rep,error=error)
       CALL uppercase(keyword)
       IF (TRIM(keyword)==TRIM(atom_name)) THEN
          CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
               keyword_name="ELEMENT",n_rep_val=i,error=error)
          IF (i>0) THEN
             CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
                  keyword_name="ELEMENT",c_val=element_symbol,error=error)
             defined_kind_section = .TRUE.
             EXIT
          ELSE
             element_symbol=element_in
             defined_kind_section = .TRUE.
          END IF
       END IF
    END DO

    ! Let's check the validity of the element so far stored..
    ! if we are not having a connectivity file, we must first match against the ptable.
    ! this helps to resolve Ca/CA (calcium and Calpha) or Cn/CN7 (Coppernicum (112) CN) conflicts
    ! so, in the presence of connectivity CA will be 'C', while in the absence of connectivity CA will be 'Ca'
    IF (defined_kind_section .OR. .NOT. use_mm_map_first) THEN
       ! lengths larger than 2 should not match, because 'trailing' characters are ignored.
       in_ptable=.FALSE.
       IF (LEN_TRIM(element_symbol)<=2) CALL get_ptable_info(element_symbol,found=in_ptable)
       IF (in_ptable) THEN
          element_out = TRIM(element_symbol)
          found=.TRUE.
       ENDIF
    ENDIF

    ! This is clearly a user error
    CALL cp_assert(found.OR.(.NOT.defined_kind_section),cp_fatal_level,cp_assertion_failed,routineP,&
         "Element <"//TRIM(element_symbol)//"> provided for KIND <"//TRIM(atom_name_in)//"> "//&
         "which cannot be mapped with any standard element label. Please correct your "//&
         "input file!"//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    ! Last chance.. are these atom_kinds of AMBER or CHARMM or GROMOS FF ?
    CALL uppercase(element_symbol)
    IF ((.NOT.found).AND.(ASSOCIATED(amber_map))) THEN
       ! First we go through the AMBER library
       DO i = 1, SIZE(amber_map%kind)
          IF (element_symbol==amber_map%kind(i)) THEN
             found =.TRUE.
             EXIT
          END IF
       END DO
       IF (found) THEN
          element_out = amber_map%element(i)
       END IF
    END IF
    IF ((.NOT.found).AND.(ASSOCIATED(charmm_map))) THEN
       ! Then we go through the CHARMM library
       DO i = 1, SIZE(charmm_map%kind)
          IF (element_symbol==charmm_map%kind(i)) THEN
             found =.TRUE.
             EXIT
          END IF
       END DO
       IF (found) THEN
          element_out = charmm_map%element(i)
       END IF
    END IF
    IF ((.NOT.found).AND.(ASSOCIATED(gromos_map))) THEN
       ! Then we go through the GROMOS library
       DO i = 1, SIZE(gromos_map%kind)
          IF (element_symbol==gromos_map%kind(i)) THEN
             found =.TRUE.
             EXIT
          END IF
       END DO
       IF (found) THEN
          element_out = gromos_map%element(i)
       END IF
    END IF

    ! final check. We has a connectivity, so we first tried to match against the ff_maps, but the element was not there.
    ! Try again the periodic table.
    IF(.NOT.found) THEN
       in_ptable=.FALSE.
       IF (LEN_TRIM(element_symbol)<=2) CALL get_ptable_info(element_symbol,found=in_ptable)
       IF (in_ptable) THEN
          element_out = TRIM(element_symbol)
          found=.TRUE.
       ENDIF
    ENDIF

    ! If no element is found the job stops here.
    IF (.NOT.found) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            " Unknown element for KIND <"//TRIM(atom_name_in)//">."//&
            " This problem can be fixed specifying properly elements in PDB"//&
            " or specifying a KIND section or getting in touch with one of "//&
            " the developers!")
    END IF

  END SUBROUTINE check_subsys_element

END MODULE topology_util

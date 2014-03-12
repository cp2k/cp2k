!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collection of subroutine needed for topology related things
!> \par History
!>     jgh (23-05-2004) Last atom of molecule information added
! *****************************************************************************

MODULE topology_connectivity_util
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_conn_g87,&
                                             do_conn_g96,&
                                             do_conn_user,&
                                             do_ff_charmm,&
                                             do_ff_harmonic
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       allocate_molecule_kind_set, atom_type, bend_type, bond_type, &
       get_molecule_kind, impr_type, molecule_kind_type, opbend_type, &
       set_molecule_kind, torsion_type, ub_type
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             get_molecule,&
                                             local_molecule_type,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE string_table,                    ONLY: id2str
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_connectivity_util'

  PRIVATE
  PUBLIC :: topology_connectivity_pack, topology_conn_multiple

CONTAINS

! *****************************************************************************
!> \brief topology connectivity pack
!> \History 11/2009 (Louis Vanduyhuys): added Out of Plane bends based on
!>                                      impropers in topology
! *****************************************************************************
  SUBROUTINE topology_connectivity_pack(molecule_kind_set,molecule_set,&
       topology,subsys_section,error)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_connectivity_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: counter, first, handle, handle2, i, ibend, ibond, idim, iimpr, &
      ikind, imol, inter_bends, inter_bonds, inter_imprs, inter_torsions, &
      inter_ubs, intra_bends, intra_bonds, intra_imprs, intra_torsions, &
      intra_ubs, inum, ires, istart_mol, istart_typ, itorsion, ityp, iub, iw, &
      j, j1, j2, j3, j4, jind, last, min_index, natom, nelectron, nsgf, &
      nval_tot1, nval_tot2, nvar1, nvar2, output_unit, stat
    INTEGER, DIMENSION(:), POINTER :: c_var_a, c_var_b, c_var_c, c_var_d, &
      c_var_type, first_list, last_list, map_atom_mol, map_atom_type, &
      map_cvar_mol, map_cvars, map_var_mol, map_vars, molecule_list
    INTEGER, DIMENSION(:, :), POINTER        :: bnd_ctype, bnd_type
    LOGICAL                                  :: failure, found, found_last
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list
    TYPE(local_molecule_type), POINTER       :: lmi
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(opbend_type), DIMENSION(:), POINTER :: opbend_list
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    TYPE(ub_type), DIMENSION(:), POINTER     :: ub_list

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    ALLOCATE(map_atom_mol(topology%natoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_atom_type(topology%natoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Set the topology%[nmol_type,nmol,nmol_conn]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_1",handle2)
    natom  = topology%natoms
    topology%nmol      = 1
    topology%nmol_type = 1
    topology%nmol_conn = 0
    map_atom_mol       =-1
    map_atom_type      =-1
    map_atom_mol(1)    = 1
    map_atom_type(1)   = 1
    DO i=1,natom-1
       IF ((atom_info%map_mol_typ(i+1) /= atom_info%map_mol_typ(i)).OR.&
           ((atom_info%map_mol_res(i+1) /= atom_info%map_mol_res(i)).AND.&
            (.NOT.(topology%conn_type == do_conn_user)))) THEN
          topology%nmol_type = topology%nmol_type + 1
       END IF
       map_atom_type(i+1) = topology%nmol_type
       IF ((atom_info%map_mol_typ(i+1) /= atom_info%map_mol_typ(i)).OR.&
           (atom_info%map_mol_num(i+1) /= atom_info%map_mol_num(i)).OR.&
           (atom_info%map_mol_res(i+1) /= atom_info%map_mol_res(i))) THEN
          topology%nmol = topology%nmol + 1
       END IF
       map_atom_mol(i+1) = topology%nmol
       IF ((atom_info%map_mol_typ(i+1) == atom_info%map_mol_typ(i)).AND.&
           (atom_info%map_mol_num(i+1) == atom_info%map_mol_num(i)).AND.&
           (atom_info%map_mol_res(i+1) /= atom_info%map_mol_res(i))) THEN
          topology%nmol_conn = topology%nmol_conn + 1
       END IF
    END DO
    IF (iw>0) WRITE(iw,*)"topology%nmol ::",topology%nmol
    IF (iw>0) WRITE(iw,*)"topology%nmol_type ::",topology%nmol_type
    IF (iw>0) WRITE(iw,*)"topology%nmol_conn ::",topology%nmol_conn
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1.1 Clean the temporary arrays to avoid quadratic loops around
    !     after this fix all topology_pack will be linear scaling
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_1.1",handle2)
    istart_mol = map_atom_mol (1)
    istart_typ = map_atom_type(1)
    DO i = 2, natom
       IF ((map_atom_mol(i)/=istart_mol).AND.(map_atom_type(i)==istart_typ)) THEN
          map_atom_mol(i) = -map_atom_mol(i)
       ELSE IF (map_atom_type(i)/=istart_typ) THEN
          istart_mol      = map_atom_mol(i)
          istart_typ      = map_atom_type(i)
       END IF
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Allocate the molecule_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_2",handle2)
    IF(topology%nmol_type <= 0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"No molecule kind defined")
    ELSE
       NULLIFY(molecule_kind_set)
       i = topology%nmol_type
       CALL allocate_molecule_kind_set(molecule_kind_set,i,error=error)
       IF(iw>0) WRITE(iw,*) "    Allocated molecule_kind_set, Dimenstion of ",&
            SIZE(molecule_kind_set)
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Allocate the molecule_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_3",handle2)
    IF(topology%nmol <= 0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"No molecule defined")
    ELSE
       NULLIFY(molecule_set)
       i = topology%nmol
       CALL allocate_molecule_set(molecule_set,i,error=error)
       IF(iw>0) WRITE(iw,*) "    Allocated molecule_set, dimension of ",&
            topology%nmol
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the molecule_kind_set%[kind_number,name,nsgf,nelectron]
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_4",handle2)
    natom  = topology%natoms
    ikind  = -1
    DO i = 1, natom
       IF (ikind /= map_atom_type(i)) THEN
          ikind = map_atom_type(i)
          molecule_kind => molecule_kind_set(ikind)
          nsgf=0
          nelectron=0
          name = TRIM(id2str(atom_info%id_molname(i)))
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               kind_number=ikind,&
               molname_generated=topology%molname_generated,&
               name=TRIM(name),&
               nsgf=nsgf,&
               nelectron=nelectron)
       END IF
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Set the molecule_list for molecule_kind in molecule_kind_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_5",handle2)
    natom   = topology%natoms
    ikind   = map_atom_type(1)
    imol    = ABS(map_atom_mol(1))
    counter = 0
    DO i = 1, natom-1
       IF (ikind /= map_atom_type(i+1)) THEN
          found      = .TRUE.
          found_last = .FALSE.
          imol       =  ABS(map_atom_mol(i))
       ELSEIF (ikind == topology%nmol_type) THEN
          found      = .TRUE.
          found_last = .TRUE.
          imol       =  ABS(map_atom_mol(natom))
       ELSE
          found      = .FALSE.
          found_last = .FALSE.
       END IF

       IF (found) THEN
          ALLOCATE(molecule_list(imol-counter),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j=1,SIZE(molecule_list)
             molecule_list(j)=j+counter
          END DO
          molecule_kind => molecule_kind_set(ikind)
          CALL set_molecule_kind(molecule_kind=molecule_kind,&
               molecule_list=molecule_list)
          IF(iw>0) WRITE(iw,*) "      molecule_list",ikind,molecule_list(:)
          IF (found_last) EXIT
          counter = imol
          ikind   = map_atom_type(i+1)
       END IF
    END DO
    ! Treat separately the case in which the last atom is also a molecule
    IF (i==natom) THEN
       imol    = ABS(map_atom_mol(natom))
       ! Last atom is also a molecule by itself
       ALLOCATE(molecule_list(imol-counter),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO j=1,SIZE(molecule_list)
          molecule_list(j)=j+counter
       END DO
       molecule_kind => molecule_kind_set(ikind)
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list)
       IF(iw>0) WRITE(iw,*) "      molecule_list",ikind,molecule_list(:)
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the molecule_set(imol)%molecule_kind via set_molecule
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_6",handle2)
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            molecule_list=molecule_list)
       DO i=1,SIZE(molecule_list)
          molecule => molecule_set(molecule_list(i))
          CALL set_molecule(molecule,molecule_kind=molecule_kind)
       END DO
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Set the molecule_set(imol)%[first_atom,last_atom] via set_molecule_set
    !-----------------------------------------------------------------------------
    ALLOCATE(first_list(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(last_list(SIZE(molecule_set)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timeset(routineN//"_7",handle2)
    first_list(:) = 0
    last_list (:) = 0
    ityp=atom_info%map_mol_typ(1)
    inum=atom_info%map_mol_num(1)
    ires=atom_info%map_mol_res(1)
    imol          = 1
    first_list(1) = 1
    DO j = 2, natom
       IF  ((atom_info%map_mol_typ(j) /= ityp) .OR. &
            (atom_info%map_mol_num(j) /= inum) .OR. &
            (atom_info%map_mol_res(j) /= ires)) THEN
          ityp = atom_info%map_mol_typ(j)
          inum = atom_info%map_mol_num(j)
          ires = atom_info%map_mol_res(j)
          imol = imol + 1
          first_list(imol)=j
       END IF
    END DO
    CPPostcondition(imol==topology%nmol,cp_failure_level,routineP,error,failure)
    DO ikind=1,topology%nmol-1
       last_list(ikind)=first_list(ikind+1)-1
    END DO
    last_list(topology%nmol)=topology%natoms
    CALL set_molecule_set(molecule_set,first_list,last_list)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Set and NULLIFY the molecule_set(imol)%lmi via set_molecule_set
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_8",handle2)
    DO imol=1,SIZE(molecule_set)
       molecule => molecule_set(imol)
       NULLIFY(lmi)
       CALL set_molecule(molecule,lmi=lmi)
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Set the atom_list for molecule_kind in molecule_kind_set (PART 1)
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_9",handle2)
    counter = 0
    DO imol=1,SIZE(molecule_set)
       molecule => molecule_set(imol)
       molecule_kind => molecule_set(imol)%molecule_kind
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            kind_number=i)
       IF(counter/=i) THEN
          counter = i
          CALL get_molecule(molecule=molecule,&
               first_atom=first,last_atom=last)
          natom = 0
          IF(first/=0 .AND. last/=0) natom = last-first+1
          ALLOCATE(atom_list(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,natom
             !Atomic kind information will be filled in (PART 2)
             NULLIFY(atom_list(i)%atomic_kind)
             atom_list(i)%id_name=atom_info%id_atmname(i+first-1)
             IF(iw>0) WRITE(iw,'(5X,A,3I5,1X,A5)') "atom_list ",&
                  imol,counter,i,TRIM(id2str(atom_list(i)%id_name))
          END DO
          CALL set_molecule_kind(molecule_kind=molecule_kind,atom_list=atom_list)
       END IF
    END DO
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Set the molecule_kind%[nbond,bond_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_10",handle2)
    ! First map bonds on molecules
    nvar1 = 0
    nvar2 = 0
    NULLIFY(map_vars, map_cvars, bnd_type, bnd_ctype)
    IF(ASSOCIATED(conn_info%bond_a))   nvar1 = SIZE(conn_info%bond_a)
    IF(ASSOCIATED(conn_info%c_bond_a)) nvar2 = SIZE(conn_info%c_bond_a)
    nval_tot1 = nvar1
    nval_tot2 = 0
    ALLOCATE(map_var_mol(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_cvar_mol(nvar2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_var_mol  = -1
    map_cvar_mol = -1
    DO i = 1, nvar1
       j1 = map_atom_mol(conn_info%bond_a(i))
       j2 = map_atom_mol(conn_info%bond_b(i))
       IF (j1==j2) THEN
          IF ( j1 > 0 ) map_var_mol(i) = map_atom_type(conn_info%bond_a(i))
       END IF
    END DO
    DO i = 1, nvar2
       min_index = MIN(conn_info%c_bond_a(i),conn_info%c_bond_b(i))
       j1 = map_atom_mol(min_index)
       IF ( j1 > 0 ) map_cvar_mol(i) = map_atom_type(min_index)
    END DO
    CALL find_bnd_typ(topology%nmol_type, map_vars,  map_var_mol,  bnd_type,  nvar1, error)
    CALL find_bnd_typ(topology%nmol_type, map_cvars, map_cvar_mol, bnd_ctype, nvar2, error)
    DO i = 1, topology%nmol_type
       intra_bonds = 0
       inter_bonds = 0
       IF (ALL(bnd_type(:,i)>0)) THEN
          intra_bonds = bnd_type(2,i)-bnd_type(1,i)+1
       END IF
       IF (ALL(bnd_ctype(:,i)>0)) THEN
          inter_bonds = bnd_ctype(2,i)-bnd_ctype(1,i)+1
       END IF
       ibond = intra_bonds+inter_bonds
       IF(iw>0) THEN
          WRITE(iw,*) "    Total number bonds for molecule type ",i," :",ibond
          WRITE(iw,*) "    intra (bonds inside  molecules) :: ",intra_bonds
          WRITE(iw,*) "    inter (bonds between molecules) :: ",inter_bonds
       END IF
       molecule_kind => molecule_kind_set(i)
       nval_tot2 = nval_tot2 + ibond*SIZE(molecule_kind%molecule_list)

       ALLOCATE(bond_list(ibond),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ibond = 0
       DO j=bnd_type(1,i),bnd_type(2,i)
          IF (j==0) CYCLE
          ibond = ibond + 1
          jind  = map_vars(j)
          first = first_list(map_atom_mol(conn_info%bond_a(jind)))
          bond_list(ibond)%a=conn_info%bond_a(jind)-first+1
          bond_list(ibond)%b=conn_info%bond_b(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          bond_list(ibond)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             bond_list(ibond)%itype   = conn_info%bond_type(jind)
          END IF
          !point this to the right bond_kind_type if using force field
          NULLIFY(bond_list(ibond)%bond_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",&
                  i,"  intra bond",&
                  conn_info%bond_a(jind),&
                  conn_info%bond_b(jind),&
                  "offset number at",&
                  conn_info%bond_a(jind)-first+1,&
                  conn_info%bond_b(jind)-first+1
          END IF
       END DO
       DO j=bnd_ctype(1,i),bnd_ctype(2,i)
          IF (j==0) CYCLE
          ibond = ibond + 1
          jind  = map_cvars(j)
          min_index = MIN(conn_info%c_bond_a(jind),conn_info%c_bond_b(jind))
          first = first_list(map_atom_mol(min_index))
          bond_list(ibond)%a=conn_info%c_bond_a(jind)-first+1
          bond_list(ibond)%b=conn_info%c_bond_b(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          bond_list(ibond)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             bond_list(ibond)%itype   = conn_info%c_bond_type(jind)
          END IF
          !point this to the right bond_kind_type if using force field
          NULLIFY(bond_list(ibond)%bond_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",&
                  i," inter bond",&
                  conn_info%c_bond_a(jind),&
                  conn_info%c_bond_b(jind),&
                  "offset number at",&
                  conn_info%c_bond_a(jind)-first+1,&
                  conn_info%c_bond_b(jind)-first+1
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            nbond=SIZE(bond_list),bond_list=bond_list)
    END DO
    IF ((nval_tot1/=nval_tot2).AND.(output_unit>0)) THEN
       WRITE(output_unit,'(/)')
       WRITE(output_unit,'(T5,A)')"ERROR| Mismatching found between the total number  of atoms"
       WRITE(output_unit,'(T5,A)')"ERROR| and the number of atoms computed multiplying the Nr."
       WRITE(output_unit,'(T5,A)')"ERROR| of molecules by the  number of atoms  building  that"
       WRITE(output_unit,'(T5,A)')"ERROR| kind of molecule."
       WRITE(output_unit,'(T5,A)')"ERROR| This happens when the connectivity is wrongly  built"
       WRITE(output_unit,'(T5,A)')"ERROR| One example could be two same kind of molecules have"
       WRITE(output_unit,'(T5,A)')"ERROR| a different number of atoms. Check the connectivity!"
    END IF
    CPPostcondition(nval_tot1==nval_tot2,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_var_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvar_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_vars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_ctype,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Set the molecule_kind%[nbend,bend_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    ! Allocate c_var_a, c_var_b, c_var_c, c_var_type
    CALL timeset(routineN//"_11_pre",handle2)
    idim = 0
    ALLOCATE(c_var_a(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_b(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_c(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    found = ASSOCIATED(conn_info%theta_type)
    IF (found) THEN
       ALLOCATE(c_var_type(idim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(conn_info%c_bond_a).AND.ASSOCIATED(conn_info%theta_a)) THEN
       DO j=1,SIZE(conn_info%theta_a)
          j1 = map_atom_mol(conn_info%theta_a(j))
          j2 = map_atom_mol(conn_info%theta_b(j))
          j3 = map_atom_mol(conn_info%theta_c(j))
          IF (j1/=j2.OR.j2/=j3) THEN
             idim = idim + 1
          END IF
       END DO
       CALL reallocate(c_var_a,1,idim)
       CALL reallocate(c_var_b,1,idim)
       CALL reallocate(c_var_c,1,idim)
       IF (found) THEN
          CALL reallocate(c_var_type,1,idim)
       END IF
       idim = 0
       DO j=1,SIZE(conn_info%theta_a)
          j1 = map_atom_mol(conn_info%theta_a(j))
          j2 = map_atom_mol(conn_info%theta_b(j))
          j3 = map_atom_mol(conn_info%theta_c(j))
          IF (j1/=j2.OR.j2/=j3) THEN
             idim = idim + 1
             c_var_a(idim)    = conn_info%theta_a(j)
             c_var_b(idim)    = conn_info%theta_b(j)
             c_var_c(idim)    = conn_info%theta_c(j)
             IF (found) THEN
                c_var_type(idim) = conn_info%theta_type(j)
             END IF
          END IF
       END DO
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_11",handle2)
    ! map bends on molecules
    nvar1 = 0
    nvar2 = 0
    NULLIFY(map_vars, map_cvars, bnd_type, bnd_ctype)
    IF(ASSOCIATED(conn_info%theta_a))   nvar1 = SIZE(conn_info%theta_a)
    IF(ASSOCIATED(c_var_a))             nvar2 = SIZE(c_var_a)
    nval_tot1 = nvar1
    nval_tot2 = 0
    ALLOCATE(map_var_mol(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_cvar_mol(nvar2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_var_mol  = -1
    map_cvar_mol = -1
    DO i = 1, nvar1
       j1 = map_atom_mol(conn_info%theta_a(i))
       j2 = map_atom_mol(conn_info%theta_b(i))
       j3 = map_atom_mol(conn_info%theta_c(i))
       IF (j1==j2.AND.j2==j3) THEN
          IF ( j1 > 0 ) map_var_mol(i) = map_atom_type(conn_info%theta_a(i))
       END IF
    END DO
    DO i = 1, nvar2
       min_index = MIN(c_var_a(i),c_var_b(i),c_var_c(i))
       j1 = map_atom_mol(min_index)
       IF ( j1 > 0 ) map_cvar_mol(i) = map_atom_type(min_index)
    END DO
    CALL find_bnd_typ(topology%nmol_type, map_vars,  map_var_mol,  bnd_type,  nvar1, error)
    CALL find_bnd_typ(topology%nmol_type, map_cvars, map_cvar_mol, bnd_ctype, nvar2, error)
    DO i=1,topology%nmol_type
       intra_bends = 0
       inter_bends = 0
       IF (ALL(bnd_type(:,i)>0)) THEN
          intra_bends = bnd_type(2,i)-bnd_type(1,i)+1
       END IF
       IF (ALL(bnd_ctype(:,i)>0)) THEN
          inter_bends = bnd_ctype(2,i)-bnd_ctype(1,i)+1
       END IF
       ibend = intra_bends+inter_bends
       IF(iw>0) THEN
          WRITE(iw,*) "    Total number of angles for molecule type ",i," :",ibend
          WRITE(iw,*) "    intra (angles inside  molecules) :: ",intra_bends
          WRITE(iw,*) "    inter (angles between molecules) :: ",inter_bends
       END IF
       molecule_kind => molecule_kind_set(i)
       nval_tot2 = nval_tot2 + ibend*SIZE(molecule_kind%molecule_list)
       ALLOCATE(bend_list(ibend),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ibend = 0
       DO j=bnd_type(1,i),bnd_type(2,i)
          IF (j==0) CYCLE
          ibend = ibend + 1
          jind  = map_vars(j)
          first = first_list(map_atom_mol(conn_info%theta_a(jind)))
          bend_list(ibend)%a=conn_info%theta_a(jind)-first+1
          bend_list(ibend)%b=conn_info%theta_b(jind)-first+1
          bend_list(ibend)%c=conn_info%theta_c(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          bend_list(ibend)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             bend_list(ibend)%itype   = conn_info%theta_type(jind)
          END IF
          !point this to the right bend_kind_type if using force field
          NULLIFY(bend_list(ibend)%bend_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                  "molecule_kind",ikind,"intra bend",&
                  conn_info%theta_a(jind),&
                  conn_info%theta_b(jind),&
                  conn_info%theta_c(jind),&
                  "offset number at",&
                  conn_info%theta_a(jind)-first+1,&
                  conn_info%theta_b(jind)-first+1,&
                  conn_info%theta_c(jind)-first+1
          END IF
       END DO
       DO j=bnd_ctype(1,i),bnd_ctype(2,i)
          IF (j==0) CYCLE
          ibend = ibend + 1
          jind  = map_cvars(j)
          min_index = MIN(c_var_a(jind),c_var_b(jind),c_var_c(jind))
          first = first_list(map_atom_mol(min_index))
          bend_list(ibend)%a=c_var_a(jind)-first+1
          bend_list(ibend)%b=c_var_b(jind)-first+1
          bend_list(ibend)%c=c_var_c(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          bend_list(ibend)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             bend_list(ibend)%itype   = c_var_type(jind)
          END IF
          !point this to the right bend_kind_type if using force field
          NULLIFY(bend_list(ibend)%bend_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                  "molecule_kind",ikind,"inter bend",&
                  c_var_a(jind),&
                  c_var_b(jind),&
                  c_var_c(jind),&
                  "offset number at",&
                  c_var_a(jind)-first+1,&
                  c_var_b(jind)-first+1,&
                  c_var_c(jind)-first+1
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            nbend=SIZE(bend_list),bend_list=bend_list)
    END DO
    CPPostcondition(nval_tot1==nval_tot2,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_var_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvar_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_vars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_ctype,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_a,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_b,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_c,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (found) THEN
       DEALLOCATE(c_var_type,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 12. Set the molecule_kind%[nub,ub_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_12_pre",handle2)
    idim = 0
    ALLOCATE(c_var_a(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_b(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_c(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(ASSOCIATED(conn_info%c_bond_a).AND.ASSOCIATED(conn_info%ub_a)) THEN
       DO j=1,SIZE(conn_info%ub_a)
          j1 = map_atom_mol(conn_info%ub_a(j))
          j2 = map_atom_mol(conn_info%ub_b(j))
          j3 = map_atom_mol(conn_info%ub_c(j))
          IF (j1/=j2.OR.j2/=j3) THEN
             idim = idim + 1
          END IF
       END DO
       CALL reallocate(c_var_a,1,idim)
       CALL reallocate(c_var_b,1,idim)
       CALL reallocate(c_var_c,1,idim)
       idim = 0
       DO j=1,SIZE(conn_info%ub_a)
          j1 = map_atom_mol(conn_info%ub_a(j))
          j2 = map_atom_mol(conn_info%ub_b(j))
          j3 = map_atom_mol(conn_info%ub_c(j))
          IF (j1/=j2.OR.j2/=j3) THEN
             idim = idim + 1
             c_var_a(idim)    = conn_info%ub_a(j)
             c_var_b(idim)    = conn_info%ub_b(j)
             c_var_c(idim)    = conn_info%ub_c(j)
          END IF
       END DO
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_12",handle2)
    ! map UBs on molecules
    nvar1 = 0
    nvar2 = 0
    NULLIFY(map_vars, map_cvars, bnd_type, bnd_ctype)
    IF(ASSOCIATED(conn_info%ub_a))   nvar1 = SIZE(conn_info%ub_a)
    IF(ASSOCIATED(c_var_a))          nvar2 = SIZE(c_var_a)
    nval_tot1 = nvar1
    nval_tot2 = 0
    ALLOCATE(map_var_mol(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_cvar_mol(nvar2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_var_mol  = -1
    map_cvar_mol = -1
    DO i = 1, nvar1
       j1 = map_atom_mol(conn_info%ub_a(i))
       j2 = map_atom_mol(conn_info%ub_b(i))
       j3 = map_atom_mol(conn_info%ub_c(i))
       IF (j1==j2.AND.j2==j3) THEN
          IF ( j1 > 0 ) map_var_mol(i) = map_atom_type(conn_info%ub_a(i))
       END IF
    END DO
    DO i = 1, nvar2
       min_index = MIN(c_var_a(i),c_var_b(i),c_var_c(i))
       j1 = map_atom_mol(min_index)
       IF ( j1 > 0 ) map_cvar_mol(i) = map_atom_type(min_index)
    END DO
    CALL find_bnd_typ(topology%nmol_type, map_vars,  map_var_mol,  bnd_type,  nvar1, error)
    CALL find_bnd_typ(topology%nmol_type, map_cvars, map_cvar_mol, bnd_ctype, nvar2, error)
    DO i=1,topology%nmol_type
       intra_ubs = 0
       inter_ubs = 0
       IF (ALL(bnd_type(:,i)>0)) THEN
          intra_ubs = bnd_type(2,i)-bnd_type(1,i)+1
       END IF
       IF (ALL(bnd_ctype(:,i)>0)) THEN
          inter_ubs = bnd_ctype(2,i)-bnd_ctype(1,i)+1
       END IF
       iub = intra_ubs+inter_ubs
       IF(iw>0) THEN
          WRITE(iw,*) "    Total number of Urey-Bradley for molecule type ",i," :",iub
          WRITE(iw,*) "    intra (UB inside  molecules) :: ",intra_ubs
          WRITE(iw,*) "    inter (UB between molecules) :: ",inter_ubs
       END IF
       molecule_kind => molecule_kind_set(i)
       nval_tot2 = nval_tot2 + iub*SIZE(molecule_kind%molecule_list)
       ALLOCATE(ub_list(iub),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       iub = 0
       DO j=bnd_type(1,i),bnd_type(2,i)
          IF (j==0) CYCLE
          iub = iub + 1
          jind  = map_vars(j)
          first = first_list(map_atom_mol(conn_info%ub_a(jind)))
          ub_list(iub)%a=conn_info%ub_a(jind)-first+1
          ub_list(iub)%b=conn_info%ub_b(jind)-first+1
          ub_list(iub)%c=conn_info%ub_c(jind)-first+1
          ub_list(iub)%id_type=do_ff_charmm
          !point this to the right ub_kind_type if using force field
          NULLIFY(ub_list(iub)%ub_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                  "molecule_kind",i,"intra UB",&
                  conn_info%ub_a(jind),&
                  conn_info%ub_b(jind),&
                  conn_info%ub_c(jind),&
                  "offset number at",&
                  conn_info%ub_a(jind)-first+1,&
                  conn_info%ub_b(jind)-first+1,&
                  conn_info%ub_c(jind)-first+1
          END IF
       END DO
       DO j=bnd_ctype(1,i),bnd_ctype(2,i)
          IF (j==0) CYCLE
          iub = iub + 1
          jind  = map_cvars(j)
          min_index = MIN(c_var_a(jind),c_var_b(jind),c_var_c(jind))
          first = first_list(map_atom_mol(min_index))
          ub_list(iub)%a=c_var_a(jind)-first+1
          ub_list(iub)%b=c_var_b(jind)-first+1
          ub_list(iub)%c=c_var_c(jind)-first+1
          ub_list(iub)%id_type=do_ff_charmm
          !point this to the right ub_kind_type if using force field
          NULLIFY(ub_list(iub)%ub_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                  "molecule_kind",i,"inter UB",&
                  c_var_a(jind),&
                  c_var_b(jind),&
                  c_var_c(jind),&
                  "offset number at",&
                  c_var_a(jind)-first+1,&
                  c_var_b(jind)-first+1,&
                  c_var_c(jind)-first+1
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            nub=SIZE(ub_list),ub_list=ub_list)
    END DO
    CPPostcondition(nval_tot1==nval_tot2,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_var_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvar_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_vars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_ctype,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_a,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_b,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_c,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 13. Set the molecule_kind%[ntorsion,torsion_list] via set_molecule_kind
    !-----------------------------------------------------------------------------
    ! Allocate c_var_a, c_var_b, c_var_c, c_var_d, c_var_type
    CALL timeset(routineN//"_13_pre",handle2)
    idim = 0
    ALLOCATE(c_var_a(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_b(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_c(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_d(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    found = ASSOCIATED(conn_info%phi_type)
    IF (found) THEN
       ALLOCATE(c_var_type(idim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(conn_info%c_bond_a).AND.ASSOCIATED(conn_info%phi_a)) THEN
       DO j=1,SIZE(conn_info%phi_a)
          j1 = map_atom_mol(conn_info%phi_a(j))
          j2 = map_atom_mol(conn_info%phi_b(j))
          j3 = map_atom_mol(conn_info%phi_c(j))
          j4 = map_atom_mol(conn_info%phi_d(j))
          IF (j1/=j2.OR.j2/=j3.OR.j3/=j4) THEN
             idim = idim + 1
          END IF
       END DO
       CALL reallocate(c_var_a,1,idim)
       CALL reallocate(c_var_b,1,idim)
       CALL reallocate(c_var_c,1,idim)
       CALL reallocate(c_var_d,1,idim)
       IF (found) THEN
          CALL reallocate(c_var_type,1,idim)
       END IF
       idim = 0
       DO j=1,SIZE(conn_info%phi_a)
          j1 = map_atom_mol(conn_info%phi_a(j))
          j2 = map_atom_mol(conn_info%phi_b(j))
          j3 = map_atom_mol(conn_info%phi_c(j))
          j4 = map_atom_mol(conn_info%phi_d(j))
          IF (j1/=j2.OR.j2/=j3.OR.j3/=j4) THEN
             idim = idim + 1
             c_var_a(idim)    = conn_info%phi_a(j)
             c_var_b(idim)    = conn_info%phi_b(j)
             c_var_c(idim)    = conn_info%phi_c(j)
             c_var_d(idim)    = conn_info%phi_d(j)
             IF (found) THEN
                c_var_type(idim) = conn_info%phi_type(j)
             END IF
          END IF
       END DO
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_13",handle2)
    ! map torsions on molecules
    nvar1 = 0
    nvar2 = 0
    NULLIFY(map_vars, map_cvars, bnd_type, bnd_ctype)
    IF(ASSOCIATED(conn_info%phi_a))   nvar1 = SIZE(conn_info%phi_a)
    IF(ASSOCIATED(c_var_a))           nvar2 = SIZE(c_var_a)
    nval_tot1 = nvar1
    nval_tot2 = 0
    ALLOCATE(map_var_mol(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_cvar_mol(nvar2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_var_mol  = -1
    map_cvar_mol = -1
    DO i = 1, nvar1
       j1 = map_atom_mol(conn_info%phi_a(i))
       j2 = map_atom_mol(conn_info%phi_b(i))
       j3 = map_atom_mol(conn_info%phi_c(i))
       j4 = map_atom_mol(conn_info%phi_d(i))
       IF (j1==j2.AND.j2==j3.AND.j3==j4) THEN
          IF ( j1 > 0 ) map_var_mol(i) = map_atom_type(conn_info%phi_a(i))
       END IF
    END DO
    DO i = 1, nvar2
       min_index = MIN(c_var_a(i),c_var_b(i),c_var_c(i),c_var_d(i))
       j1 = map_atom_mol(min_index)
       IF ( j1 > 0 ) map_cvar_mol(i) = map_atom_type(min_index)
    END DO
    CALL find_bnd_typ(topology%nmol_type, map_vars,  map_var_mol,  bnd_type,  nvar1, error)
    CALL find_bnd_typ(topology%nmol_type, map_cvars, map_cvar_mol, bnd_ctype, nvar2, error)
    DO i=1,topology%nmol_type
       intra_torsions = 0
       inter_torsions = 0
       IF (ALL(bnd_type(:,i)>0)) THEN
          intra_torsions = bnd_type(2,i)-bnd_type(1,i)+1
       END IF
       IF (ALL(bnd_ctype(:,i)>0)) THEN
          inter_torsions = bnd_ctype(2,i)-bnd_ctype(1,i)+1
       END IF
       itorsion = intra_torsions+inter_torsions
       IF(iw>0) THEN
          WRITE(iw,*) "    Total number of torsions for molecule type ",i," :",itorsion
          WRITE(iw,*) "    intra (torsions inside  molecules) :: ",intra_torsions
          WRITE(iw,*) "    inter (torsions between molecules) :: ",inter_torsions
       END IF
       molecule_kind => molecule_kind_set(i)
       nval_tot2 = nval_tot2 + itorsion*SIZE(molecule_kind%molecule_list)
       ALLOCATE(torsion_list(itorsion),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       itorsion = 0
       DO j=bnd_type(1,i),bnd_type(2,i)
          IF (j==0) CYCLE
          itorsion = itorsion + 1
          jind  = map_vars(j)
          first = first_list(map_atom_mol(conn_info%phi_a(jind)))
          torsion_list(itorsion)%a=conn_info%phi_a(jind)-first+1
          torsion_list(itorsion)%b=conn_info%phi_b(jind)-first+1
          torsion_list(itorsion)%c=conn_info%phi_c(jind)-first+1
          torsion_list(itorsion)%d=conn_info%phi_d(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          torsion_list(itorsion)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             torsion_list(itorsion)%itype   = conn_info%phi_type(jind)
          END IF
          !point this to the right torsion_kind_type if using force field
          NULLIFY(torsion_list(itorsion)%torsion_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"intra TOR",&
                  conn_info%phi_a(jind),&
                  conn_info%phi_b(jind),&
                  conn_info%phi_c(jind),&
                  conn_info%phi_d(jind),&
                  "offset number at",&
                  conn_info%phi_a(jind)-first+1,&
                  conn_info%phi_b(jind)-first+1,&
                  conn_info%phi_c(jind)-first+1,&
                  conn_info%phi_d(jind)-first+1
          END IF
       END DO
       DO j=bnd_ctype(1,i),bnd_ctype(2,i)
          IF (j==0) CYCLE
          itorsion = itorsion + 1
          jind  = map_cvars(j)
          min_index = MIN(c_var_a(jind),c_var_b(jind),c_var_c(jind),c_var_d(jind))
          first = first_list(map_atom_mol(min_index))
          torsion_list(itorsion)%a=c_var_a(jind)-first+1
          torsion_list(itorsion)%b=c_var_b(jind)-first+1
          torsion_list(itorsion)%c=c_var_c(jind)-first+1
          torsion_list(itorsion)%d=c_var_d(jind)-first+1
          ! Set by default id_type to charmm and modify when handling the forcefield
          torsion_list(itorsion)%id_type = do_ff_charmm
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             torsion_list(itorsion)%itype   = c_var_type(jind)
          END IF
          !point this to the right torsion_kind_type if using force field
          NULLIFY(torsion_list(itorsion)%torsion_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"inter TOR",&
                  c_var_a(jind),&
                  c_var_b(jind),&
                  c_var_c(jind),&
                  c_var_d(jind),&
                  "offset number at",&
                  c_var_a(jind)-first+1,&
                  c_var_b(jind)-first+1,&
                  c_var_c(jind)-first+1,&
                  c_var_d(jind)-first+1
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            ntorsion=SIZE(torsion_list),torsion_list=torsion_list)
    END DO
    CPPostcondition(nval_tot1==nval_tot2,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_var_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvar_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_vars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_ctype,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_a,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_b,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_c,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_d,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (found) THEN
       DEALLOCATE(c_var_type,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 14. Set the molecule_kind%[nimpr,impr_list] via set_molecule_kind
    !     Also set the molecule_kind%[nopbend,opbend_list]
    !-----------------------------------------------------------------------------
    ! Allocate c_var_a, c_var_b, c_var_c, c_var_d, c_var_type
    CALL timeset(routineN//"_14_pre",handle2)
    idim = 0
    ALLOCATE(c_var_a(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_b(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_c(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_var_d(idim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    found = ASSOCIATED(conn_info%impr_type)
    IF (found) THEN
       ALLOCATE(c_var_type(idim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(conn_info%c_bond_a).AND.ASSOCIATED(conn_info%impr_a)) THEN
       DO j=1,SIZE(conn_info%impr_a)
          j1 = map_atom_mol(conn_info%impr_a(j))
          j2 = map_atom_mol(conn_info%impr_b(j))
          j3 = map_atom_mol(conn_info%impr_c(j))
          j4 = map_atom_mol(conn_info%impr_d(j))
          IF (j1/=j2.OR.j2/=j3.OR.j3/=j4) THEN
             idim = idim + 1
          END IF
       END DO
       CALL reallocate(c_var_a,1,idim)
       CALL reallocate(c_var_b,1,idim)
       CALL reallocate(c_var_c,1,idim)
       CALL reallocate(c_var_d,1,idim)
       IF (found) THEN
          CALL reallocate(c_var_type,1,idim)
       END IF
       idim = 0
       DO j=1,SIZE(conn_info%impr_a)
          j1 = map_atom_mol(conn_info%impr_a(j))
          j2 = map_atom_mol(conn_info%impr_b(j))
          j3 = map_atom_mol(conn_info%impr_c(j))
          j4 = map_atom_mol(conn_info%impr_d(j))
          IF (j1/=j2.OR.j2/=j3.OR.j3/=j4) THEN
             idim = idim + 1
             c_var_a(idim)    = conn_info%impr_a(j)
             c_var_b(idim)    = conn_info%impr_b(j)
             c_var_c(idim)    = conn_info%impr_c(j)
             c_var_d(idim)    = conn_info%impr_d(j)
             IF (found) THEN
                c_var_type(idim) = conn_info%impr_type(j)
             END IF
          END IF
       END DO
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_14",handle2)
    ! map imprs on molecules
    nvar1 = 0
    nvar2 = 0
    NULLIFY(map_vars, map_cvars, bnd_type, bnd_ctype)
    IF(ASSOCIATED(conn_info%impr_a))   nvar1 = SIZE(conn_info%impr_a)
    IF(ASSOCIATED(c_var_a))            nvar2 = SIZE(c_var_a)
    nval_tot1 = nvar1
    nval_tot2 = 0
    ALLOCATE(map_var_mol(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_cvar_mol(nvar2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_var_mol  = -1
    map_cvar_mol = -1
    DO i = 1, nvar1
       j1 = map_atom_mol(conn_info%impr_a(i))
       j2 = map_atom_mol(conn_info%impr_b(i))
       j3 = map_atom_mol(conn_info%impr_c(i))
       j4 = map_atom_mol(conn_info%impr_d(i))
       IF (j1==j2.AND.j2==j3.AND.j3==j4) THEN
          IF ( j1 > 0 ) map_var_mol(i) = map_atom_type(conn_info%impr_a(i))
       END IF
    END DO
    DO i = 1, nvar2
       min_index = MIN(c_var_a(i),c_var_b(i),c_var_c(i),c_var_d(i))
       j1 = map_atom_mol(min_index)
       IF ( j1 > 0 ) map_cvar_mol(i) = map_atom_type(min_index)
    END DO
    CALL find_bnd_typ(topology%nmol_type, map_vars,  map_var_mol,  bnd_type,  nvar1, error)
    CALL find_bnd_typ(topology%nmol_type, map_cvars, map_cvar_mol, bnd_ctype, nvar2, error)
    DO i=1,topology%nmol_type
       intra_imprs  = 0
       inter_imprs  = 0
       IF (ALL(bnd_type(:,i)>0)) THEN
          intra_imprs  = bnd_type(2,i)-bnd_type(1,i)+1
       END IF
       IF (ALL(bnd_ctype(:,i)>0)) THEN
          inter_imprs  = bnd_ctype(2,i)-bnd_ctype(1,i)+1
       END IF
       iimpr   = intra_imprs+inter_imprs
       IF(iw>0) THEN
          WRITE(iw,*) "    Total number of imprs for molecule type ",i," :",iimpr
          WRITE(iw,*) "    intra (imprs inside  molecules) :: ",intra_imprs
          WRITE(iw,*) "    inter (imprs between molecules) :: ",inter_imprs
          WRITE(iw,*) "    Total number of opbends for molecule type ",i," :",iimpr
          WRITE(iw,*) "    intra (opbends inside  molecules) :: ",intra_imprs
          WRITE(iw,*) "    inter (opbends between molecules) :: ",inter_imprs
       END IF
       molecule_kind => molecule_kind_set(i)
       nval_tot2 = nval_tot2 + iimpr*SIZE(molecule_kind%molecule_list)
       ALLOCATE(impr_list(iimpr),STAT=stat)
       ALLOCATE(opbend_list(iimpr),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       iimpr  = 0
       DO j=bnd_type(1,i),bnd_type(2,i)
          IF (j==0) CYCLE
          iimpr  = iimpr + 1
          jind  = map_vars(j)
          first = first_list(map_atom_mol(conn_info%impr_a(jind)))
          impr_list(iimpr)%a=conn_info%impr_a(jind)-first+1
          impr_list(iimpr)%b=conn_info%impr_b(jind)-first+1
          impr_list(iimpr)%c=conn_info%impr_c(jind)-first+1
          impr_list(iimpr)%d=conn_info%impr_d(jind)-first+1
          ! Atom sequence for improper is A B C D in which A is central atom,
          ! B is deviating atom and C & D are secondairy atoms. Atom sequence for
          ! opbend is B D C A in which A is central atom, B is deviating. Hence
          ! to create an opbend out of an improper, B and D need to be interchanged.
          opbend_list(iimpr)%a=conn_info%impr_b(jind)-first+1
          opbend_list(iimpr)%b=conn_info%impr_d(jind)-first+1
          opbend_list(iimpr)%c=conn_info%impr_c(jind)-first+1
          opbend_list(iimpr)%d=conn_info%impr_a(jind)-first+1
          ! Set by default id_type of improper to charmm and modify when handling the forcefield
          impr_list(iimpr)%id_type = do_ff_charmm
          ! Set by default id_type of opbend to harmonic and modify when handling the forcefield
          opbend_list(iimpr)%id_type = do_ff_harmonic
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             impr_list(iimpr)%itype   = conn_info%impr_type(jind)
          END IF
          !point this to the right impr_kind_type if using force field
          NULLIFY(impr_list(iimpr)%impr_kind)
          NULLIFY(opbend_list(iimpr)%opbend_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"intra IMPR",&
                  conn_info%impr_a(jind),&
                  conn_info%impr_b(jind),&
                  conn_info%impr_c(jind),&
                  conn_info%impr_d(jind),&
                  "offset number at",&
                  conn_info%impr_a(jind)-first+1,&
                  conn_info%impr_b(jind)-first+1,&
                  conn_info%impr_c(jind)-first+1,&
                  conn_info%impr_d(jind)-first+1
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"intra OPBEND",&
                  conn_info%impr_b(jind),&
                  conn_info%impr_d(jind),&
                  conn_info%impr_c(jind),&
                  conn_info%impr_a(jind),&
                  "offset number at",&
                  conn_info%impr_b(jind)-first+1,&
                  conn_info%impr_d(jind)-first+1,&
                  conn_info%impr_c(jind)-first+1,&
                  conn_info%impr_a(jind)-first+1
          END IF
       END DO
       DO j=bnd_ctype(1,i),bnd_ctype(2,i)
          IF (j==0) CYCLE
          iimpr = iimpr + 1
          jind  = map_cvars(j)
          min_index = MIN(c_var_a(jind),c_var_b(jind),c_var_c(jind),c_var_d(jind))
          first = first_list(map_atom_mol(min_index))
          impr_list(iimpr)%a=c_var_a(jind)-first+1
          impr_list(iimpr)%b=c_var_b(jind)-first+1
          impr_list(iimpr)%c=c_var_c(jind)-first+1
          impr_list(iimpr)%d=c_var_d(jind)-first+1
          opbend_list(iimpr)%a=c_var_b(jind)-first+1
          opbend_list(iimpr)%b=c_var_d(jind)-first+1
          opbend_list(iimpr)%c=c_var_c(jind)-first+1
          opbend_list(iimpr)%d=c_var_a(jind)-first+1
          ! Set by default id_type of improper to charmm and modify when handling the forcefield
          impr_list(iimpr)%id_type = do_ff_charmm
          ! Set by default id_type of opbend to harmonic and modify when handling the forcefield
          opbend_list(iimpr)%id_type = do_ff_harmonic
          IF((topology%conn_type==do_conn_g96).OR.(topology%conn_type==do_conn_g87)) THEN
             impr_list(iimpr)%itype   = c_var_type(jind)
          END IF
          !point this to the right impr_kind_type and opbend_kind_type if using force field
          NULLIFY(impr_list(iimpr)%impr_kind)
          NULLIFY(opbend_list(iimpr)%opbend_kind)
          IF(iw>0) THEN
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"inter IMPR",&
                  c_var_a(jind),&
                  c_var_b(jind),&
                  c_var_c(jind),&
                  c_var_d(jind),&
                  "offset number at",&
                  c_var_a(jind)-first+1,&
                  c_var_b(jind)-first+1,&
                  c_var_c(jind)-first+1,&
                  c_var_d(jind)-first+1
             WRITE(iw,'(7X,A,I3,1X,A,I4,I4,I4,I4,1X,A,I4,I4,I4,I4)') &
                  "molecule_kind",i,"inter OPBEND",&
                  c_var_b(jind),&
                  c_var_d(jind),&
                  c_var_c(jind),&
                  c_var_a(jind),&
                  "offset number at",&
                  c_var_b(jind)-first+1,&
                  c_var_d(jind)-first+1,&
                  c_var_c(jind)-first+1,&
                  c_var_a(jind)-first+1
          END IF
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            nimpr=SIZE(impr_list),impr_list=impr_list)
       CALL set_molecule_kind(molecule_kind=molecule_kind,&
            nopbend=SIZE(opbend_list),opbend_list=opbend_list)
    END DO
    CPPostcondition(nval_tot1==nval_tot2,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_var_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvar_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_vars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_cvars,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bnd_ctype,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_a,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_b,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_c,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(c_var_d,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (found) THEN
       DEALLOCATE(c_var_type,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! Finally deallocate some stuff.
    !-----------------------------------------------------------------------------
    DEALLOCATE(first_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_atom_mol,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(map_atom_type,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
  END SUBROUTINE topology_connectivity_pack

! *****************************************************************************
!> \brief used to achieve linear scaling in the connectivity_pack
!> \par History
!>      Teodoro Laino
! *****************************************************************************
  SUBROUTINE find_bnd_typ(nmol_type, map_vars, map_var_mol, bnd_type, nvar1, error)
    INTEGER, INTENT(IN)                      :: nmol_type
    INTEGER, DIMENSION(:), POINTER           :: map_vars, map_var_mol
    INTEGER, DIMENSION(:, :), POINTER        :: bnd_type
    INTEGER, INTENT(IN)                      :: nvar1
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'find_bnd_typ', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ibond, j, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE(map_vars(nvar1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL sort(map_var_mol,nvar1,map_vars)
    ALLOCATE(bnd_type(2,nmol_type),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    bnd_type = 0
    IF (nvar1==0) RETURN
    DO j = 1, nvar1
       IF (map_var_mol(j)/=-1) EXIT
    END DO
    IF (j==nvar1+1) RETURN
    i = map_var_mol(j)
    bnd_type(1,i) = j
    DO ibond = j, nvar1
       IF (map_var_mol(ibond) /=i) THEN
          bnd_type(2,i) = ibond-1
          i = map_var_mol(ibond)
          bnd_type(1,i) = ibond
       END IF
    END DO
    bnd_type(2,i) = nvar1

  END SUBROUTINE find_bnd_typ

! *****************************************************************************
!> \brief   Handles the multiple unit cell option for the connectivity
!> \author  Teodoro Laino [tlaino] - 06.2009
! *****************************************************************************
  SUBROUTINE topology_conn_multiple(topology, subsys_section, error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_conn_multiple', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, fac, i, ind, j, k, m, &
                                                natoms_orig, nbond, nbond_c, &
                                                nimpr, nonfo, nphi, ntheta, &
                                                nub
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL                                  :: failure
    TYPE(connectivity_info_type), POINTER    :: conn_info

    failure = .FALSE.
    NULLIFY(multiple_unit_cell)
    CALL section_vals_val_get(subsys_section, "TOPOLOGY%MULTIPLE_UNIT_CELL",&
                              i_vals=multiple_unit_cell, error=error)
    IF (ANY(multiple_unit_cell/=1)) THEN
       fac = PRODUCT(multiple_unit_cell)
       conn_info => topology%conn_info

       nbond = 0
       IF (ASSOCIATED(conn_info%bond_a)) THEN
          nbond  = SIZE(conn_info%bond_a)
          CALL reallocate(conn_info%bond_a,1,nbond*fac)
          CALL reallocate(conn_info%bond_b,1,nbond*fac)
       END IF

       ntheta = 0
       IF (ASSOCIATED(conn_info%theta_a)) THEN
          ntheta = SIZE(conn_info%theta_a)
          CALL reallocate(conn_info%theta_a,1,ntheta*fac)
          CALL reallocate(conn_info%theta_b,1,ntheta*fac)
          CALL reallocate(conn_info%theta_c,1,ntheta*fac)
       END IF

       nphi = 0
       IF (ASSOCIATED(conn_info%phi_a)) THEN
          nphi = SIZE(conn_info%phi_a)
          CALL reallocate(conn_info%phi_a,1,nphi*fac)
          CALL reallocate(conn_info%phi_b,1,nphi*fac)
          CALL reallocate(conn_info%phi_c,1,nphi*fac)
          CALL reallocate(conn_info%phi_d,1,nphi*fac)
       END IF

       nimpr = 0
       IF (ASSOCIATED(conn_info%impr_a)) THEN
          nimpr = SIZE(conn_info%impr_a)
          CALL reallocate(conn_info%impr_a,1,nimpr*fac)
          CALL reallocate(conn_info%impr_b,1,nimpr*fac)
          CALL reallocate(conn_info%impr_c,1,nimpr*fac)
          CALL reallocate(conn_info%impr_d,1,nimpr*fac)
       END IF

       nbond_c = 0
       IF (ASSOCIATED(conn_info%c_bond_a)) THEN
          nbond_c = SIZE(conn_info%c_bond_a)
          CALL reallocate(conn_info%c_bond_a,1,nbond_c*fac)
          CALL reallocate(conn_info%c_bond_b,1,nbond_c*fac)
       END IF

       nub = 0
       IF (ASSOCIATED(conn_info%ub_a)) THEN
          nub = SIZE(conn_info%ub_a)
          CALL reallocate(conn_info%ub_a,1,nub*fac)
          CALL reallocate(conn_info%ub_b,1,nub*fac)
          CALL reallocate(conn_info%ub_c,1,nub*fac)
       END IF

       nonfo = 0
       IF (ASSOCIATED(conn_info%onfo_a)) THEN
          nonfo = SIZE(conn_info%onfo_a)
          CALL reallocate(conn_info%onfo_a,1,nonfo*fac)
          CALL reallocate(conn_info%onfo_b,1,nonfo*fac)
       END IF

       natoms_orig = topology%natoms/fac
       ind = 0
       DO k = 1, multiple_unit_cell(3)
          DO j = 1, multiple_unit_cell(2)
             DO i = 1, multiple_unit_cell(1)
                ind = ind + 1
                IF (ind==1) CYCLE
                a = (ind-1)*natoms_orig

                ! Bonds
                IF (nbond>0) THEN
                   m = (ind-1)*nbond
                   conn_info%bond_a(m+1:m+nbond)=conn_info%bond_a(1:nbond)+a
                   conn_info%bond_b(m+1:m+nbond)=conn_info%bond_b(1:nbond)+a
                END IF
                ! Theta
                IF (ntheta>0) THEN
                   m = (ind-1)*ntheta
                   conn_info%theta_a(m+1:m+ntheta)=conn_info%theta_a(1:ntheta)+a
                   conn_info%theta_b(m+1:m+ntheta)=conn_info%theta_b(1:ntheta)+a
                   conn_info%theta_c(m+1:m+ntheta)=conn_info%theta_c(1:ntheta)+a
                END IF
                ! Phi
                IF (nphi>0) THEN
                   m = (ind-1)*nphi
                   conn_info%phi_a(m+1:m+nphi)=conn_info%phi_a(1:nphi)+a
                   conn_info%phi_b(m+1:m+nphi)=conn_info%phi_b(1:nphi)+a
                   conn_info%phi_c(m+1:m+nphi)=conn_info%phi_c(1:nphi)+a
                   conn_info%phi_d(m+1:m+nphi)=conn_info%phi_d(1:nphi)+a
                END IF
                ! Impropers
                IF (nimpr>0) THEN
                   m = (ind-1)*nimpr
                   conn_info%impr_a(m+1:m+nimpr)=conn_info%impr_a(1:nimpr)+a
                   conn_info%impr_b(m+1:m+nimpr)=conn_info%impr_b(1:nimpr)+a
                   conn_info%impr_c(m+1:m+nimpr)=conn_info%impr_c(1:nimpr)+a
                   conn_info%impr_d(m+1:m+nimpr)=conn_info%impr_d(1:nimpr)+a
                END IF
                ! Para_res
                IF (nbond_c>0) THEN
                   m = (ind-1)*nbond_c
                   conn_info%c_bond_a(m+1:m+nbond_c)=conn_info%c_bond_a(1:nbond_c)+a
                   conn_info%c_bond_b(m+1:m+nbond_c)=conn_info%c_bond_b(1:nbond_c)+a
                END IF
                ! Urey-Bradley
                IF (nub>0) THEN
                   m = (ind-1)*nub
                   conn_info%ub_a(m+1:m+nub)=conn_info%ub_a(1:nub)+a
                   conn_info%ub_b(m+1:m+nub)=conn_info%ub_b(1:nub)+a
                   conn_info%ub_c(m+1:m+nub)=conn_info%ub_c(1:nub)+a
                END IF
                ! ONFO
                IF (nonfo>0) THEN
                   m = (ind-1)*nonfo
                   conn_info%onfo_a(m+1:m+nonfo)=conn_info%onfo_a(1:nonfo)+a
                   conn_info%onfo_b(m+1:m+nonfo)=conn_info%onfo_b(1:nonfo)+a
                END IF
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE topology_conn_multiple

END MODULE topology_connectivity_util

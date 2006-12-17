!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/topology_constraint_util [1.0] *
!!
!!   NAME
!!     topology_constraint_util
!!
!!   FUNCTION
!!     Collection of subroutine needed for topology related things
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!    jgh (23-05-2004) Last atom of molecule information added
!!
!!   SOURCE
!******************************************************************************

MODULE topology_constraint_util
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_clone,&
                                             colvar_counters,&
                                             colvar_create,&
                                             colvar_p_reallocate,&
                                             colvar_setup,&
                                             dist_colvar_id,&
                                             torsion_colvar_id
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE input_constants,                 ONLY: do_constr_none, do_constr_atomic, do_constr_molec
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       atom_type, bond_type, colvar_constraint_type, fixd_constraint_type, &
       g3x3_constraint_type, g4x6_constraint_type, get_molecule_kind, &
       molecule_kind_type, set_molecule_kind, setup_colvar_counters
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_colvar_constraint_type,&
                                             local_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type,&
                                             molecule_type,&
                                             set_molecule
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: constraint_info_type,&
                                             topology_parameters_type,&
                                             constr_list_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_constraint_util'

  PRIVATE
  PUBLIC :: topology_constraint_pack


!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_constraint_util/topology_constraint_pack [1.0] *
!!
!!   NAME
!!     topology_constraint_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

  SUBROUTINE topology_constraint_pack ( molecule_kind_set,molecule_set,&
       topology,qmmm_env,particle_set,input_file,para_env,subsys_section,error)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_constraint_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: molname, name
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_typeh
    INTEGER :: first, first_atom, handle, i, ii, iw, j, k, k1loc, k2loc, &
      kdim, kk, last, last_atom, m, n_start_colv, natom, nbond, ncolv_mol, &
      nfixd_restraint, nfixed_atoms, ng3x3, ng3x3_restraint, ng4x6, &
      ng4x6_restraint, nhdist, nmolecule, nrep, offset, stat
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: failure, fix_atom_mm, &
                                                fix_atom_molname, &
                                                fix_atom_qm, fix_atom_qmmm, &
                                                fix_fixed_atom
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: missed_molname
    REAL(KIND=dp)                            :: rmod, rvec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: hdist, r
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(colvar_counters)                    :: ncolv
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    TYPE(local_colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: lcolv
    TYPE(local_constraint_type), POINTER     :: lci
    TYPE(local_g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: lg3x3
    TYPE(local_g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: lg4x6
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: hbonds_section,&
                                                fixd_restr_rest,&
                                                colvar_rest
    TYPE(constr_list_type), POINTER, DIMENSION(:)  :: constr_x_mol
    INTEGER :: output_unit, gind
    LOGICAL :: ionode, ishbond, ldummy
    INTEGER, POINTER, DIMENSION(:) :: inds
    CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: cnds
    LOGICAL :: restart_restraint_pos, restart_restraint_clv

    failure = .FALSE.
    NULLIFY(logger, constr_x_mol)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,'I','',handle)

    cons_info => topology%cons_info
    hbonds_section => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%HBONDS",error=error)
    fixd_restr_rest => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%FIX_ATOM_RESTART",error=error)
    CALL section_vals_get(fixd_restr_rest,explicit=restart_restraint_pos, error=error)
    colvar_rest => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%COLVAR_RESTART",error=error)
    CALL section_vals_get(colvar_rest,explicit=restart_restraint_clv, error=error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. NULLIFY the molecule_set(imol)%lci via set_molecule_set
    !-----------------------------------------------------------------------------
    DO i=1,topology%nmol
       molecule => molecule_set(i)
       NULLIFY(lci)
       ALLOCATE(lci,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(lci%lcolv)
       NULLIFY(lci%lg3x3)
       NULLIFY(lci%lg4x6)
       CALL set_molecule(molecule,lci=lci)
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Add more stuff to COLVAR constraint if constraint hydrogen is on
    !-----------------------------------------------------------------------------
    IF(topology%const_hydr) THEN
       topology%const_colv = .TRUE.
       NULLIFY(atom_typeh, hdist)
       ALLOCATE(constr_x_mol(SIZE(molecule_kind_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, SIZE(molecule_kind_set)
          ALLOCATE(constr_x_mol(i)%constr(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          constr_x_mol(i)%constr(1) = 1
       END DO
       CALL section_vals_val_get(hbonds_section,"MOLECULE",n_rep_val=nrep,error=error)
       IF (nrep/=0) THEN
          NULLIFY(inds)
          DO i = 1, SIZE(molecule_kind_set)
             constr_x_mol(i)%constr(1) = 0
          END DO
          CALL section_vals_val_get(hbonds_section,"MOLECULE",i_vals=inds,error=error)
          DO i = 1, SIZE(inds)
             constr_x_mol(inds(i))%constr(1) = 1
          END DO
       ELSE
          CALL section_vals_val_get(hbonds_section,"MOLNAME",n_rep_val=nrep,error=error)
          IF (nrep/=0) THEN
             NULLIFY(cnds)
             DO i = 1, SIZE(molecule_kind_set)
                constr_x_mol(i)%constr(1) = 0
             END DO
             CALL section_vals_val_get(hbonds_section,"MOLNAME",c_vals=cnds,error=error)
             DO i = 1, SIZE(cnds)
                DO k = 1, SIZE(molecule_kind_set)
                   molecule_kind => molecule_kind_set(k)
                   name = molecule_kind%name
                   ldummy = qmmm_ff_precond_only_qm(id1=name)
                   IF (cnds(i) == name) constr_x_mol(k)%constr(1) = 1
                END DO
             END DO
          END IF
       END IF
       CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",n_rep_val=nrep,error=error)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",c_vals=atom_typeh,error=error)
       CALL section_vals_val_get(hbonds_section,"TARGETS",n_rep_val=nrep,error=error)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"TARGETS",r_vals=hdist,error=error)
       IF (ASSOCIATED(hdist)) THEN
          CPPostcondition(SIZE(hdist)==SIZE(atom_typeh),cp_failure_level,routineP,error,failure)
       END IF
       nhdist = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          IF (constr_x_mol(i)%constr(1) == 0) CYCLE
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bond_list=bond_list, nbond=nbond, atom_list=atom_list,&
               molecule_list=molecule_list)
          ! Let's tag all requested atoms involving Hydrogen
          ! on the first molecule of this kind
          molecule => molecule_set(molecule_list(1))
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
          natom = last_atom - first_atom +1
          DO k=1,nbond
             ishbond = .FALSE.
             j = bond_list(k)%a
             IF (j<1.OR.j>natom)CYCLE
             atomic_kind => atom_list(j)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name,element_symbol=element_symbol)
             ldummy = qmmm_ff_precond_only_qm(id1=name)
             IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H"))    ishbond = .TRUE.
             IF (.NOT.ishbond) THEN
                j = bond_list(k)%b
                IF (j<1.OR.j>natom)CYCLE
                atomic_kind => atom_list(j)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name,element_symbol=element_symbol)
                ldummy = qmmm_ff_precond_only_qm(id1=name)
                IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H")) ishbond = .TRUE.
             END IF
             IF (ishbond) THEN
                nhdist = nhdist + 1
             END IF
          END DO
       END DO
       n_start_colv = cons_info%nconst_colv
       cons_info%nconst_colv = nhdist + n_start_colv
       CALL reallocate(cons_info%const_colv_mol,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%const_colv_molname,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%const_colv_target,1,cons_info%nconst_colv)
       CALL colvar_p_reallocate(cons_info%colvar_set,1,cons_info%nconst_colv,error=error)
       ! Fill in Restraints info
       CALL reallocate(cons_info%colv_restraint,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%colv_k0,1,cons_info%nconst_colv)
       cons_info%colv_restraint(n_start_colv+1:cons_info%nconst_colv) = cons_info%hbonds_restraint
       cons_info%colv_k0(n_start_colv+1:cons_info%nconst_colv)        = cons_info%hbonds_k0
       !
       nhdist = 0
       DO i=1,SIZE(molecule_kind_set)
          IF (constr_x_mol(i)%constr(1) == 0) CYCLE
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               bond_list=bond_list, nbond=nbond, atom_list=atom_list,&
               molecule_list=molecule_list)
          molecule => molecule_set(molecule_list(1))
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
          natom = last_atom  - first_atom + 1
          offset = first_atom - 1
          DO k=1,nbond
             ishbond = .FALSE.
             j = bond_list(k)%a
             IF (j<1.OR.j>natom)CYCLE
             atomic_kind => atom_list(j)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=name,element_symbol=element_symbol)
             ldummy = qmmm_ff_precond_only_qm(id1=name)
             IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H"))    ishbond = .TRUE.
             IF (.NOT.ishbond) THEN
                j = bond_list(k)%b
                IF (j<1.OR.j>natom)CYCLE
                atomic_kind => atom_list(j)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     name=name,element_symbol=element_symbol)
                ldummy = qmmm_ff_precond_only_qm(id1=name)
                IF((name(1:1)=="H").OR.(TRIM(element_symbol)=="H")) ishbond = .TRUE.
             END IF
             IF (ishbond) THEN
                nhdist = nhdist + 1
                rvec = particle_set(offset+bond_list(k)%a)%r-particle_set(offset+bond_list(k)%b)%r
                rmod = SQRT(DOT_PRODUCT(rvec,rvec))
                IF (ASSOCIATED(hdist)) THEN
                   IF (SIZE(hdist)>0) THEN
                      IF (bond_list(k)%a==j) atomic_kind => atom_list(bond_list(k)%b)%atomic_kind
                      IF (bond_list(k)%b==j) atomic_kind => atom_list(bond_list(k)%a)%atomic_kind
                      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name,element_symbol=element_symbol)
                      ldummy = qmmm_ff_precond_only_qm(id1=name)
                      DO m = 1, SIZE(hdist)
                         IF (TRIM(name)==TRIM(atom_typeh(m))) EXIT
                         IF (TRIM(element_symbol)==TRIM(atom_typeh(m))) EXIT
                      END DO
                      IF (m<=SIZE(hdist)) THEN
                         rmod = hdist(m)
                      END IF
                   END IF
                END IF
                cons_info%const_colv_mol(nhdist+n_start_colv)     = i
                cons_info%const_colv_molname(nhdist+n_start_colv) = "UNDEF"
                cons_info%const_colv_target(nhdist+n_start_colv)  = rmod
                CALL colvar_create(cons_info%colvar_set(nhdist+n_start_colv)%colvar,&
                     dist_colvar_id, error)
                cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%i_at = bond_list(k)%a
                cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%j_at = bond_list(k)%b
                CALL colvar_setup(cons_info%colvar_set(nhdist+n_start_colv)%colvar, error)
             END IF
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(constr_x_mol,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Set the COLVAR constraint molecule_kind_set(ikind)%colv_list
    !-----------------------------------------------------------------------------
    IF(topology%const_colv) THEN
       CALL give_constraint_array(cons_info%const_colv_mol,&
                                  cons_info%const_colv_molname,&
                                  constr_x_mol,&
                                  molecule_kind_set,&
                                  error)
       gind = 0
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ncolv_mol = SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(colv_list(ncolv_mol),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ncolv_mol = 0
          DO kk = 1, SIZE(constr_x_mol(ii)%constr)
             j = constr_x_mol(ii)%constr(kk)
             ncolv_mol = ncolv_mol + 1
             gind = gind + 1
             kdim = SIZE(cons_info%colvar_set(j)%colvar%i_atom)
             ALLOCATE(colv_list(ncolv_mol)%i_atoms(kdim),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colv_list(ncolv_mol)%type_id        = cons_info%colvar_set(j)%colvar%type_id
             colv_list(ncolv_mol)%i_atoms        = cons_info%colvar_set(j)%colvar%i_atom
             ! Restraint
             colv_list(ncolv_mol)%restraint%active = cons_info%colv_restraint(j)
             colv_list(ncolv_mol)%restraint%k0     = cons_info%colv_k0(j)
             IF (cons_info%const_colv_target(j)==-HUGE(0.0_dp)) THEN
                CALL colvar_eval_mol_f(cons_info%colvar_set(j)%colvar, topology%cell,&
                     particle_set, error=error)
                colv_list(ncolv_mol)%expected_value = cons_info%colvar_set(j)%colvar%ss
             ELSE
                colv_list(ncolv_mol)%expected_value = cons_info%const_colv_target(j)
             END IF
             ! In case of Restraint let's check for possible restart values
             IF (colv_list(ncolv_mol)%restraint%active) THEN
                IF (restart_restraint_clv) THEN
                   CALL section_vals_val_get(colvar_rest,"_DEFAULT_KEYWORD_",&
                              i_rep_val=gind, r_val=rmod, error=error)
                   colv_list(ncolv_mol)%expected_value = rmod
                ELSE
                   rmod = colv_list(ncolv_mol)%expected_value
                   CALL section_vals_val_set(colvar_rest,"_DEFAULT_KEYWORD_",&
                        i_rep_val=gind, r_val=rmod, error=error)
                END IF
             END IF
             ! Only if torsion let's take into account the singularity in the definition
             ! of the dihedral
             IF (cons_info%colvar_set(j)%colvar%type_id==torsion_colvar_id) THEN
                cons_info%colvar_set(j)%colvar%torsion_param%o0=colv_list(ncolv_mol)%expected_value
             END IF
          END DO
          CALL setup_colvar_counters(colv_list,ncolv)
          CALL set_molecule_kind(molecule_kind,colv_list=colv_list,ncolv=ncolv)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lcolv(ncolv_mol),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ncolv_mol = 0
             DO kk=1,SIZE(constr_x_mol(ii)%constr)
                k = constr_x_mol(ii)%constr(kk)
                ncolv_mol = ncolv_mol + 1
                lcolv(ncolv_mol)%lambda    = 0.0_dp
                lcolv(ncolv_mol)%sigma     = 0.0_dp
                ! Set Up colvar variable
                NULLIFY(lcolv(ncolv_mol)%colvar, lcolv(ncolv_mol)%colvar_old)
                ! Colvar
                CALL colvar_clone(lcolv(ncolv_mol)%colvar, cons_info%colvar_set(k)%colvar,&
                     i_atom_offset=first_atom-1,error=error)
                ! Colvar_old
                CALL colvar_clone(lcolv(ncolv_mol)%colvar_old,lcolv(ncolv_mol)%colvar,error=error)
                ! Check for consistency in the constraint definition
                IF  (ANY(lcolv(ncolv_mol)%colvar%i_atom > last_atom).OR.&
                     ANY(lcolv(ncolv_mol)%colvar%i_atom < first_atom)) THEN
                   WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                   WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                        " but the atoms specified in the constraint and the atoms defined for",&
                        " the molecule DO NOT match!",&
                        "This could be very probable due to a wrong connectivity, or an error",&
                        " in the constraint specification in the input file.",&
                        " Please check it carefully!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lcolv=lcolv)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(constr_x_mol,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the group 3x3 constraint g3x3_list
    !-----------------------------------------------------------------------------
    IF(topology%const_33) THEN
       CALL give_constraint_array(cons_info%const_g33_mol,&
                                  cons_info%const_g33_molname,&
                                  constr_x_mol,&
                                  molecule_kind_set,&
                                  error)
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ng3x3 =  SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(g3x3_list(ng3x3),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ng3x3 = 0
          ng3x3_restraint = 0
          DO kk = 1, SIZE(constr_x_mol(ii)%constr)
             j = constr_x_mol(ii)%constr(kk)
             ng3x3 = ng3x3 + 1
             g3x3_list(ng3x3)%a   = cons_info%const_g33_a(j)
             g3x3_list(ng3x3)%b   = cons_info%const_g33_b(j)
             g3x3_list(ng3x3)%c   = cons_info%const_g33_c(j)
             g3x3_list(ng3x3)%dab = cons_info%const_g33_dab(j)
             g3x3_list(ng3x3)%dac = cons_info%const_g33_dac(j)
             g3x3_list(ng3x3)%dbc = cons_info%const_g33_dbc(j)
             ! Restraint
             g3x3_list(ng3x3)%restraint%active = cons_info%g33_restraint(j)
             g3x3_list(ng3x3)%restraint%k0     = cons_info%g33_k0(j)
             IF (g3x3_list(ng3x3)%restraint%active) ng3x3_restraint = ng3x3_restraint+ 1
          END DO
          CALL set_molecule_kind(molecule_kind,ng3x3=ng3x3,ng3x3_restraint=ng3x3_restraint,g3x3_list=g3x3_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg3x3(ng3x3),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO kk = 1, SIZE(constr_x_mol(ii)%constr)
                k = constr_x_mol(ii)%constr(kk)
                lg3x3(k)%scale     = 0.0_dp
                lg3x3(k)%scale_old = 0.0_dp
                lg3x3(k)%fa        = 0.0_dp
                lg3x3(k)%fb        = 0.0_dp
                lg3x3(k)%fc        = 0.0_dp
                lg3x3(k)%ra_old    = 0.0_dp
                lg3x3(k)%rb_old    = 0.0_dp
                lg3x3(k)%rc_old    = 0.0_dp
                lg3x3(k)%va        = 0.0_dp
                lg3x3(k)%vb        = 0.0_dp
                lg3x3(k)%vc        = 0.0_dp
                lg3x3(k)%lambda    = 0.0_dp
                IF  ((g3x3_list(k)%a+first_atom-1 < first_atom).OR.&
                     (g3x3_list(k)%b+first_atom-1 < first_atom).OR.&
                     (g3x3_list(k)%c+first_atom-1 < first_atom).OR.&
                     (g3x3_list(k)%a+first_atom-1 > last_atom ).OR.&
                     (g3x3_list(k)%b+first_atom-1 > last_atom ).OR.&
                     (g3x3_list(k)%c+first_atom-1 > last_atom )) THEN
                   WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                   WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                        " but the atoms specified in the constraint and the atoms defined for",&
                        " the molecule DO NOT match!",&
                        "This could be very probable due to a wrong connectivity, or an error",&
                        " in the constraint specification in the input file.",&
                        " Please check it carefully!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lg3x3=lg3x3)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(constr_x_mol,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Set the group 4x6 constraint g4x6_list
    !-----------------------------------------------------------------------------
    IF(topology%const_46) THEN
       CALL give_constraint_array(cons_info%const_g46_mol,&
                                  cons_info%const_g46_molname,&
                                  constr_x_mol,&
                                  molecule_kind_set,&
                                  error)
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ng4x6 = SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(g4x6_list(ng4x6),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ng4x6 = 0
          ng4x6_restraint = 0
          DO kk = 1, SIZE(constr_x_mol(ii)%constr)
             j = constr_x_mol(ii)%constr(kk)
             ng4x6 = ng4x6 + 1
             g4x6_list(ng4x6)%a   = cons_info%const_g46_a(j)
             g4x6_list(ng4x6)%b   = cons_info%const_g46_b(j)
             g4x6_list(ng4x6)%c   = cons_info%const_g46_c(j)
             g4x6_list(ng4x6)%d   = cons_info%const_g46_d(j)
             g4x6_list(ng4x6)%dab = cons_info%const_g46_dab(j)
             g4x6_list(ng4x6)%dac = cons_info%const_g46_dac(j)
             g4x6_list(ng4x6)%dbc = cons_info%const_g46_dbc(j)
             g4x6_list(ng4x6)%dad = cons_info%const_g46_dad(j)
             g4x6_list(ng4x6)%dbd = cons_info%const_g46_dbd(j)
             g4x6_list(ng4x6)%dcd = cons_info%const_g46_dcd(j)
             ! Restraint
             g4x6_list(ng4x6)%restraint%active = cons_info%g46_restraint(j)
             g4x6_list(ng4x6)%restraint%k0     = cons_info%g46_k0(j)
             IF (g4x6_list(ng4x6)%restraint%active) ng4x6_restraint = ng4x6_restraint+ 1
          END DO
          CALL set_molecule_kind(molecule_kind,ng4x6=ng4x6,ng4x6_restraint=ng4x6_restraint,g4x6_list=g4x6_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg4x6(ng4x6),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO kk = 1, SIZE(constr_x_mol(ii)%constr)
                k = constr_x_mol(ii)%constr(kk)
                lg4x6(k)%scale     = 0.0_dp
                lg4x6(k)%scale_old = 0.0_dp
                lg4x6(k)%fa        = 0.0_dp
                lg4x6(k)%fb        = 0.0_dp
                lg4x6(k)%fc        = 0.0_dp
                lg4x6(k)%fd        = 0.0_dp
                lg4x6(k)%fe        = 0.0_dp
                lg4x6(k)%ff        = 0.0_dp
                lg4x6(k)%ra_old    = 0.0_dp
                lg4x6(k)%rb_old    = 0.0_dp
                lg4x6(k)%rc_old    = 0.0_dp
                lg4x6(k)%rd_old    = 0.0_dp
                lg4x6(k)%re_old    = 0.0_dp
                lg4x6(k)%rf_old    = 0.0_dp
                lg4x6(k)%va        = 0.0_dp
                lg4x6(k)%vb        = 0.0_dp
                lg4x6(k)%vc        = 0.0_dp
                lg4x6(k)%vd        = 0.0_dp
                lg4x6(k)%ve        = 0.0_dp
                lg4x6(k)%vf        = 0.0_dp
                lg4x6(k)%lambda    = 0.0_dp
                IF  ((g4x6_list(k)%a+first_atom-1 < first_atom).OR.&
                     (g4x6_list(k)%b+first_atom-1 < first_atom).OR.&
                     (g4x6_list(k)%c+first_atom-1 < first_atom).OR.&
                     (g4x6_list(k)%d+first_atom-1 < first_atom).OR.&
                     (g4x6_list(k)%a+first_atom-1 > last_atom ).OR.&
                     (g4x6_list(k)%b+first_atom-1 > last_atom ).OR.&
                     (g4x6_list(k)%c+first_atom-1 > last_atom ).OR.&
                     (g4x6_list(k)%d+first_atom-1 > last_atom )) THEN
                   WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
                   WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
                        " but the atoms specified in the constraint and the atoms defined for",&
                        " the molecule DO NOT match!",&
                        "This could be very probable due to a wrong connectivity, or an error",&
                        " in the constraint specification in the input file.",&
                        " Please check it carefully!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END DO
             CALL set_molecule(molecule=molecule,lg4x6=lg4x6)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(constr_x_mol,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the group fixed_atom constraint fixd_list
    !-----------------------------------------------------------------------------
    IF(topology%const_atom) THEN
       ALLOCATE(missed_molname(SIZE(cons_info%fixed_molnames,1)),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       missed_molname=.TRUE.
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list, name=molname)
          ldummy = qmmm_ff_precond_only_qm(id1=molname)
          WHERE(molname.EQ.cons_info%fixed_molnames)
            missed_molname=.FALSE.
          END WHERE
          ! Try to figure out how many atoms of the list belong to this molecule_kind
          nfixed_atoms = 0
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first, last_atom = last )
             fix_atom_molname=.FALSE.
             IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                IF (ANY(cons_info%fixed_molnames.EQ.molname)) fix_atom_molname=.TRUE.
             ENDIF
             DO k = first, last
                fix_atom_qmmm=.FALSE.
                IF (PRESENT(qmmm_env)) THEN
                   SELECT CASE(cons_info%freeze_qm)
                   CASE(do_constr_atomic)
                      IF (ANY(qmmm_env%qm_atom_index==k)) fix_atom_qmmm=.TRUE.
                   CASE(do_constr_molec)
                      IF (ANY(qmmm_env%qm_molecule_index == molecule_list(j))) fix_atom_qmmm = .TRUE.
                   END SELECT
                   SELECT CASE(cons_info%freeze_mm)
                   CASE(do_constr_atomic)
                      IF (ALL(qmmm_env%qm_atom_index/=k)) fix_atom_qmmm=.TRUE.
                   CASE(do_constr_molec)
                      IF (ALL(qmmm_env%qm_molecule_index /= molecule_list(j))) fix_atom_qmmm = .TRUE.
                   END SELECT
                END IF
                IF  (ANY(cons_info%fixed_atoms == k) .OR. fix_atom_qmmm .OR. fix_atom_molname) THEN
                        nfixed_atoms = nfixed_atoms + 1
                END IF
             END DO
          END DO
          ALLOCATE(fixd_list(nfixed_atoms),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          kk = 0
          nfixd_restraint = 0
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first, last_atom = last )
             fix_atom_molname=.FALSE.
             IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                DO k1loc = 1, SIZE(cons_info%fixed_molnames)
                   IF (cons_info%fixed_molnames(k1loc).EQ.molname) THEN
                      fix_atom_molname=.TRUE.
                      EXIT
                   END IF
                END DO
             ENDIF
             DO k=first,last
                ! FIXED LIST ATOMS
                fix_fixed_atom = .FALSE.
                DO k2loc=1,SIZE(cons_info%fixed_atoms)
                   IF (cons_info%fixed_atoms(k2loc) == k) THEN
                      fix_fixed_atom = .TRUE.
                      EXIT
                   END IF
                END DO
                ! QMMM FIXED ATOMS (QM OR MM)
                fix_atom_qmmm=.FALSE.
                fix_atom_mm  =.FALSE.
                fix_atom_qm  =.FALSE.
                IF (PRESENT(qmmm_env)) THEN
                   SELECT CASE(cons_info%freeze_qm)
                   CASE(do_constr_atomic)
                      IF (ANY(qmmm_env%qm_atom_index==k)) THEN 
                         fix_atom_qmmm=.TRUE.
                         fix_atom_qm = .TRUE.
                      END IF
                   CASE(do_constr_molec)
                      IF (ANY(qmmm_env%qm_molecule_index == molecule_list(j))) THEN 
                         fix_atom_qmmm = .TRUE.
                         fix_atom_qm   = .TRUE.
                      END IF
                   END SELECT
                   SELECT CASE(cons_info%freeze_mm)
                   CASE(do_constr_atomic)
                      IF (ALL(qmmm_env%qm_atom_index/=k)) THEN 
                         fix_atom_qmmm=.TRUE.
                         fix_atom_mm = .TRUE.
                      END IF
                   CASE(do_constr_molec)
                      IF (ALL(qmmm_env%qm_molecule_index /= molecule_list(j))) THEN 
                         fix_atom_qmmm = .TRUE.
                         fix_atom_mm = .TRUE.
                      END IF
                   END SELECT
                   ! We should never reach this point but let's check it anyway
                   IF (fix_atom_qm.AND.fix_atom_mm) THEN
                      CALL stop_program(TRIM(routineN),"Atom number: "//cp_to_string(k)//&
                           " has been defined both QM and MM. General Error!")
                   END IF
                END IF
                ! Check that the fixed atom constraint/restraint is unique
                IF ((fix_fixed_atom.AND.fix_atom_qmmm).OR.(fix_fixed_atom.AND.fix_atom_molname)&
                     .OR.(fix_atom_qmmm.AND.fix_atom_molname)) THEN 
                   CALL stop_program(TRIM(routineN),"Atom number: "//cp_to_string(k)//&
                        " has been constrained/restrained to be fixed in more than one"//&
                        " input section. Check and correct your input file!")
                END IF
                ! Let's store the atom index
                IF (fix_fixed_atom.OR.fix_atom_qmmm.OR.fix_atom_molname) THEN
                   kk = kk + 1
                   fixd_list(kk)%fixd  = k
                   fixd_list(kk)%coord = particle_set(k)%r
                   ! Possibly Restraint
                   IF     (fix_fixed_atom) THEN
                      fixd_list(kk)%restraint%active = cons_info%fixed_restraint(k2loc)
                      fixd_list(kk)%restraint%k0     = cons_info%fixed_k0(k2loc)
                   ELSEIF (fix_atom_qm) THEN
                      fixd_list(kk)%restraint%active = cons_info%fixed_qm_restraint
                      fixd_list(kk)%restraint%k0     = cons_info%fixed_qm_k0
                   ELSEIF (fix_atom_mm) THEN
                      fixd_list(kk)%restraint%active = cons_info%fixed_mm_restraint
                      fixd_list(kk)%restraint%k0     = cons_info%fixed_mm_k0
                   ELSEIF (fix_atom_molname) THEN
                      fixd_list(kk)%restraint%active = cons_info%fixed_mol_restraint(k1loc)
                      fixd_list(kk)%restraint%k0     = cons_info%fixed_mol_k0(k1loc)
                   ELSE
                      ! Should never reach this point
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END IF
                   IF (fixd_list(kk)%restraint%active) THEN 
                      nfixd_restraint = nfixd_restraint+ 1
                      IF (restart_restraint_pos) THEN
                         ! Read  coord0 value for restraint
                         CALL section_vals_val_get(fixd_restr_rest,"_DEFAULT_KEYWORD_",&
                              i_rep_val=nfixd_restraint, r_vals=r, error=error)
                         fixd_list(kk)%coord = r
                      ELSE
                         ! Write coord0 value for restraint
                         ALLOCATE(r(3),stat=stat)
                         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                         r = fixd_list(kk)%coord
                         CALL section_vals_val_set(fixd_restr_rest,"_DEFAULT_KEYWORD_",&
                              i_rep_val=nfixd_restraint, r_vals_ptr=r, error=error)                         
                      END IF
                   END IF
                END IF
             END DO
          END DO
          IF (iw>0) WRITE(iw,*)"FIXD_LIST ::",fixd_list
          CALL set_molecule_kind(molecule_kind, nfixd=nfixed_atoms, nfixd_restraint=nfixd_restraint,&
               fixd_list=fixd_list)
       END DO
       CPPostcondition(COUNT(missed_molname)==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(missed_molname,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE topology_constraint_pack

!******************************************************************************
!!****** topology_constraint_util/give_constraint_array [1.0] *
!!
!!   NAME
!!     give_constraint_array
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Gives back a list of molecule to which apply the constraint
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************
  SUBROUTINE give_constraint_array(const_mol, const_molname,constr_x_mol,&
                                  molecule_kind_set, error)
    
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error
    INTEGER, POINTER, DIMENSION(:)           :: const_mol
    CHARACTER(LEN=default_string_length),&
         POINTER, DIMENSION(:)               :: const_molname
    TYPE(constr_list_type), POINTER, DIMENSION(:)  :: constr_x_mol

    CHARACTER(len=*), PARAMETER :: routineN = 'give_constraint_array', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure, ldummy
    INTEGER :: handle, i, k, isize, stat
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    CHARACTER(LEN=default_string_length)     :: myname, name

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       NULLIFY(molecule_kind)
       ALLOCATE(constr_x_mol(SIZE(molecule_kind_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, SIZE(constr_x_mol)
          NULLIFY(constr_x_mol(i)%constr)
          ALLOCATE(constr_x_mol(i)%constr(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       CPPostcondition(SIZE(const_mol)==SIZE(const_molname),cp_failure_level,routineP,error,failure)
       DO i = 1, SIZE(const_mol)
          IF (const_mol(i)/=0) THEN 
             k = const_mol(i)
             isize = SIZE(constr_x_mol(k)%constr)
             CALL reallocate(constr_x_mol(k)%constr, 1, isize+1)
             constr_x_mol(k)%constr(isize+1) = i
          ELSE
             myname = const_molname(i)
             DO k = 1, SIZE(molecule_kind_set)
                molecule_kind => molecule_kind_set(k)
                name = molecule_kind%name
                ldummy = qmmm_ff_precond_only_qm(id1=name)
                IF (name==myname) THEN
                   isize = SIZE(constr_x_mol(k)%constr)
                   CALL reallocate(constr_x_mol(k)%constr, 1, isize+1)
                   constr_x_mol(k)%constr(isize+1) = i
                END IF
             END DO
          END IF
       END DO
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE give_constraint_array

END MODULE topology_constraint_util

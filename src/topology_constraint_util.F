!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collection of subroutine needed for topology related things
!> \par History
!>     jgh (23-05-2004) Last atom of molecule information added
! *****************************************************************************
MODULE topology_constraint_util
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             is_hydrogen
  USE cell_types,                      ONLY: use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: &
       colvar_clone, colvar_counters, colvar_create, colvar_p_reallocate, &
       colvar_release, colvar_setup, colvar_type, dist_colvar_id, &
       torsion_colvar_id, xyz_diag_colvar_id, xyz_outerdiag_colvar_id
  USE colvar_utils,                    ONLY: post_process_colvar
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: do_constr_atomic,&
                                             do_constr_molec
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       atom_type, bond_type, colvar_constraint_type, fixd_constraint_type, &
       g3x3_constraint_type, g4x6_constraint_type, get_molecule_kind, &
       molecule_kind_type, set_molecule_kind, setup_colvar_counters, &
       vsite_constraint_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_colvar_constraint_type,&
                                             local_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type,&
                                             molecule_type,&
                                             set_molecule
  USE particle_types,                  ONLY: particle_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types_low,                  ONLY: qmmm_env_mm_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: constr_list_type,&
                                             constraint_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_constraint_util'

  PRIVATE
  PUBLIC :: topology_constraint_pack

CONTAINS

! *****************************************************************************
!> \brief Pack in all the information needed for the constraints
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param topology ...
!> \param qmmm_env ...
!> \param particle_set ...
!> \param input_file ...
!> \param subsys_section ...
!> \param gci ...
! *****************************************************************************
  SUBROUTINE topology_constraint_pack ( molecule_kind_set, molecule_set,&
       topology,qmmm_env,particle_set,input_file,subsys_section,gci)
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input_file, subsys_section
    TYPE(global_constraint_type), POINTER    :: gci

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_constraint_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: molname, name
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_typeh, cnds
    INTEGER :: cind, first, first_atom, gind, handle, handle2, i, ii, itype, &
      iw, j, k, k1loc, k2loc, kk, last, last_atom, m, n_start_colv, natom, &
      nbond, ncolv_glob, ncolv_mol, nfixd_list_gci, nfixd_restart, &
      nfixd_restraint, nfixed_atoms, ng3x3, ng3x3_restraint, ng4x6, &
      ng4x6_restraint, nhdist, nmolecule, nrep, nvsite, nvsite_restraint, &
      offset
    INTEGER, DIMENSION(:), POINTER           :: constr_x_glob, inds, &
                                                molecule_list
    LOGICAL :: exclude_mm, exclude_qm, failure, fix_atom_mm, &
      fix_atom_molname, fix_atom_qm, fix_atom_qmmm, fix_fixed_atom, &
      found_molname, is_qm, ishbond, ldummy, restart_restraint_clv, &
      restart_restraint_pos, use_clv_info
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: missed_molname
    REAL(KIND=dp)                            :: rmod, rvec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: hdist, r
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(colvar_counters)                    :: ncolv
    TYPE(constr_list_type), DIMENSION(:), &
      POINTER                                :: constr_x_mol
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list, fixd_list_gci
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
    TYPE(section_vals_type), POINTER         :: colvar_func_info, &
                                                colvar_rest, fixd_restr_rest, &
                                                hbonds_section
    TYPE(vsite_constraint_type), &
      DIMENSION(:), POINTER                  :: vsite_list

    failure = .FALSE.
    NULLIFY(logger, constr_x_mol, constr_x_glob)
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog")
    CALL timeset(routineN,handle)
    CALL timeset(routineN//"_1",handle2)

    cons_info => topology%cons_info
    hbonds_section => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%HBONDS")
    fixd_restr_rest => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%FIX_ATOM_RESTART")
    CALL section_vals_get(fixd_restr_rest,explicit=restart_restraint_pos)
    colvar_rest => section_vals_get_subs_vals(input_file,&
         "MOTION%CONSTRAINT%COLVAR_RESTART")
    CALL section_vals_get(colvar_rest,explicit=restart_restraint_clv)
    colvar_func_info => section_vals_get_subs_vals(subsys_section,&
         "COLVAR%COLVAR_FUNC_INFO")
    CALL section_vals_get(colvar_func_info,explicit=use_clv_info)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. NULLIFY the molecule_set(imol)%lci via set_molecule_set
    !-----------------------------------------------------------------------------
    DO i=1,topology%nmol
       molecule => molecule_set(i)
       NULLIFY(lci)
       ! only allocate the lci if constraints are active. Can this stuff be distributed ?
       IF (topology%const_atom .OR. topology%const_hydr .OR. &
           topology%const_33   .OR. topology%const_46   .OR. &
           topology%const_colv ) THEN
          ALLOCATE(lci)
          NULLIFY(lci%lcolv)
          NULLIFY(lci%lg3x3)
          NULLIFY(lci%lg4x6)
       ENDIF
       CALL set_molecule(molecule,lci=lci)
    END DO
    ALLOCATE(gci)
    NULLIFY(gci%lcolv,&
            gci%lg3x3,&
            gci%lg4x6,&
            gci%fixd_list,&
            gci%colv_list,&
            gci%g3x3_list,&
            gci%g4x6_list,&
            gci%vsite_list)
    gci%ntot   = 0
    gci%ng3x3  = 0
    gci%ng4x6  = 0
    gci%nvsite = 0
    gci%ng3x3_restraint  = 0
    gci%ng4x6_restraint  = 0
    gci%nvsite_restraint  = 0
    CALL setup_colvar_counters(gci%colv_list,gci%ncolv)
    gci%nrestraint = gci%ng3x3_restraint +&
                     gci%ng4x6_restraint +&
                     gci%nvsite_restraint +&
                     gci%ncolv%nrestraint
    CALL timestop(handle2)
    CALL timeset(routineN//"_2",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Add more stuff to COLVAR constraint if constraint hydrogen is on
    !-----------------------------------------------------------------------------
    IF(topology%const_hydr) THEN
       topology%const_colv = .TRUE.
       NULLIFY(atom_typeh, hdist)
       ALLOCATE(constr_x_mol(SIZE(molecule_kind_set)))
       DO i = 1, SIZE(molecule_kind_set)
          ALLOCATE(constr_x_mol(i)%constr(1))
          constr_x_mol(i)%constr(1) = 1
       END DO
       CALL section_vals_val_get(hbonds_section,"MOLECULE",n_rep_val=nrep)
       IF (nrep/=0) THEN
          NULLIFY(inds)
          DO i = 1, SIZE(molecule_kind_set)
             constr_x_mol(i)%constr(1) = 0
          END DO
          CALL section_vals_val_get(hbonds_section,"MOLECULE",i_vals=inds)
          DO i = 1, SIZE(inds)
             constr_x_mol(inds(i))%constr(1) = 1
          END DO
       ELSE
          CALL section_vals_val_get(hbonds_section,"MOLNAME",n_rep_val=nrep)
          IF (nrep/=0) THEN
             NULLIFY(cnds)
             DO i = 1, SIZE(molecule_kind_set)
                constr_x_mol(i)%constr(1) = 0
             END DO
             CALL section_vals_val_get(hbonds_section,"MOLNAME",c_vals=cnds)
             DO i = 1, SIZE(cnds)
                found_molname = .FALSE.
                DO k = 1, SIZE(molecule_kind_set)
                   molecule_kind => molecule_kind_set(k)
                   name = molecule_kind%name
                   ldummy = qmmm_ff_precond_only_qm(id1=name)
                   IF (cnds(i) == name) THEN
                      constr_x_mol(k)%constr(1) = 1
                      found_molname = .TRUE.
                   END IF
                END DO
                CALL print_warning_molname(found_molname, cnds(i))
             END DO
          END IF
       END IF
       CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",n_rep_val=nrep)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"ATOM_TYPE",c_vals=atom_typeh)
       CALL section_vals_val_get(hbonds_section,"TARGETS",n_rep_val=nrep)
       IF (nrep/=0) &
            CALL section_vals_val_get(hbonds_section,"TARGETS",r_vals=hdist)
       IF (ASSOCIATED(hdist)) THEN
          CPASSERT(SIZE(hdist)==SIZE(atom_typeh))
       END IF
       CALL section_vals_val_get(hbonds_section,"exclude_qm",l_val=exclude_qm)
       CALL section_vals_val_get(hbonds_section,"exclude_mm",l_val=exclude_mm)
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
             CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
             is_qm = qmmm_ff_precond_only_qm(id1=name)
             IF((name(1:1)=="H").OR.is_hydrogen(atomic_kind))    ishbond = .TRUE.
             IF(      is_qm .AND.exclude_qm) ishbond = .FALSE.
             IF(.NOT.(is_qm).AND.exclude_mm) ishbond = .FALSE.
             IF (.NOT.ishbond) THEN
                j = bond_list(k)%b
                IF (j<1.OR.j>natom)CYCLE
                atomic_kind => atom_list(j)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
                is_qm = qmmm_ff_precond_only_qm(id1=name)
                IF((name(1:1)=="H").OR.is_hydrogen(atomic_kind)) ishbond = .TRUE.
                IF(      is_qm .AND.exclude_qm) ishbond = .FALSE.
                IF(.NOT.(is_qm).AND.exclude_mm) ishbond = .FALSE.
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
       CALL reallocate(cons_info%const_colv_target_growth,1,cons_info%nconst_colv)
       CALL colvar_p_reallocate(cons_info%colvar_set,1,cons_info%nconst_colv)
       ! Fill in Restraints info
       CALL reallocate(cons_info%colv_intermolecular,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%colv_restraint,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%colv_k0,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%colv_exclude_qm,1,cons_info%nconst_colv)
       CALL reallocate(cons_info%colv_exclude_mm,1,cons_info%nconst_colv)
       ! Bonds involving hydrogens are by their nature only intramolecular
       cons_info%colv_intermolecular(n_start_colv+1:cons_info%nconst_colv) = .FALSE.
       cons_info%colv_exclude_qm(n_start_colv+1:cons_info%nconst_colv) = .FALSE.
       cons_info%colv_exclude_mm(n_start_colv+1:cons_info%nconst_colv) = .FALSE.
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
             CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
             is_qm = qmmm_ff_precond_only_qm(id1=name)
             IF((name(1:1)=="H").OR.is_hydrogen(atomic_kind))    ishbond = .TRUE.
             IF(      is_qm .AND.exclude_qm) ishbond = .FALSE.
             IF(.NOT.(is_qm).AND.exclude_mm) ishbond = .FALSE.
             IF (.NOT.ishbond) THEN
                j = bond_list(k)%b
                IF (j<1.OR.j>natom)CYCLE
                atomic_kind => atom_list(j)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=name)
                is_qm = qmmm_ff_precond_only_qm(id1=name)
                IF((name(1:1)=="H").OR.is_hydrogen(atomic_kind)) ishbond = .TRUE.
                IF(      is_qm .AND.exclude_qm) ishbond = .FALSE.
                IF(.NOT.(is_qm).AND.exclude_mm) ishbond = .FALSE.
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
                cons_info%const_colv_target_growth(nhdist+n_start_colv)  = 0.0_dp
                CALL colvar_create(cons_info%colvar_set(nhdist+n_start_colv)%colvar,&
                     dist_colvar_id)
                cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%i_at = bond_list(k)%a
                cons_info%colvar_set(nhdist+n_start_colv)%colvar%dist_param%j_at = bond_list(k)%b
                CALL colvar_setup(cons_info%colvar_set(nhdist+n_start_colv)%colvar)
             END IF
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr)
       END DO
       DEALLOCATE(constr_x_mol)
    END IF

    CALL timestop(handle2)
    CALL timeset(routineN//"_3",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Set the COLVAR constraint molecule_kind_set(ikind)%colv_list
    !-----------------------------------------------------------------------------
    IF(topology%const_colv) THEN
       ! Post Process of COLVARS..
       DO ii = 1, SIZE(cons_info%colvar_set)
          CALL post_process_colvar(cons_info%colvar_set(ii)%colvar, particle_set)
       END DO
       ! Real constraint/restraint part..
       CALL give_constraint_array(cons_info%const_colv_mol,&
                                  cons_info%const_colv_molname,&
                                  cons_info%colv_intermolecular,&
                                  constr_x_mol,&
                                  constr_x_glob,&
                                  molecule_kind_set,&
                                  cons_info%colv_exclude_qm,&
                                  cons_info%colv_exclude_mm)
       ! Intramolecular constraints
       gind = 0
       cind = 0
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule, molecule_list=molecule_list)
          ncolv_mol = SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(colv_list(ncolv_mol))
          ! Starting index of the first molecule of this kind.
          ! We need the index if no target is provided in the input file
          ! for the collective variable.. The target will be computed on the
          ! first molecule of the kind...
          molecule => molecule_set(molecule_list(1))
          CALL get_molecule ( molecule, first_atom = first_atom)
          CALL setup_colv_list(colv_list, constr_x_mol(ii)%constr, gind,&
               cons_info, topology, particle_set, restart_restraint_clv,&
               colvar_rest, first_atom)
          CALL setup_colvar_counters(colv_list,ncolv)
          CALL set_molecule_kind(molecule_kind,colv_list=colv_list,ncolv=ncolv)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lcolv(ncolv_mol))
             CALL setup_lcolv(lcolv, constr_x_mol(ii)%constr, first_atom, last_atom,&
                  cons_info, particle_set, colvar_func_info, use_clv_info, cind)
             CALL set_molecule(molecule=molecule,lcolv=lcolv)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr)
       END DO
       DEALLOCATE(constr_x_mol)
       ! Intermolecular constraints
       ncolv_glob = 0
       IF (ASSOCIATED(constr_x_glob)) THEN
          ncolv_glob = SIZE(constr_x_glob)
          ALLOCATE(colv_list(ncolv_glob))
          CALL setup_colv_list(colv_list, constr_x_glob, gind, cons_info,&
               topology, particle_set, restart_restraint_clv, colvar_rest,&
               first_atom=1)
          CALL setup_colvar_counters(colv_list,ncolv)
          ALLOCATE(lcolv(ncolv_glob))
          CALL setup_lcolv(lcolv, constr_x_glob, 1, SIZE(particle_set), cons_info,&
                           particle_set, colvar_func_info, use_clv_info, cind)
          gci%colv_list => colv_list
          gci%lcolv     => lcolv
          gci%ncolv     = ncolv
          ! Total number of Intermolecular constraints
          gci%ntot = gci%ncolv%ntot + gci%ntot
          DEALLOCATE(constr_x_glob)
       END IF
    END IF

    CALL timestop(handle2)
    CALL timeset(routineN//"_4",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Set the group 3x3 constraint g3x3_list
    !-----------------------------------------------------------------------------
    IF(topology%const_33) THEN
       CALL give_constraint_array(cons_info%const_g33_mol,&
                                  cons_info%const_g33_molname,&
                                  cons_info%g33_intermolecular,&
                                  constr_x_mol,&
                                  constr_x_glob,&
                                  molecule_kind_set,&
                                  cons_info%g33_exclude_qm,&
                                  cons_info%g33_exclude_mm)
       ! Intramolecular constraints
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
                                 nmolecule=nmolecule,&
                                 molecule_list=molecule_list)
          ng3x3 =  SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(g3x3_list(ng3x3))
          CALL setup_g3x3_list(g3x3_list, constr_x_mol(ii)%constr, cons_info, ng3x3_restraint)
          CALL set_molecule_kind(molecule_kind,ng3x3=ng3x3,ng3x3_restraint=ng3x3_restraint,g3x3_list=g3x3_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg3x3(ng3x3))
             CALL setup_lg3x3(lg3x3, g3x3_list, first_atom, last_atom)
             CALL set_molecule(molecule=molecule,lg3x3=lg3x3)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr)
       END DO
       DEALLOCATE(constr_x_mol)
       ! Intermolecular constraints
       IF (ASSOCIATED(constr_x_glob)) THEN
          ng3x3 =  SIZE(constr_x_glob)
          ALLOCATE(g3x3_list(ng3x3))
          CALL setup_g3x3_list(g3x3_list, constr_x_glob, cons_info, ng3x3_restraint)
          ALLOCATE(lg3x3(ng3x3))
          CALL setup_lg3x3(lg3x3, g3x3_list, first_atom,last_atom)
          gci%g3x3_list => g3x3_list
          gci%lg3x3     => lg3x3
          gci%ng3x3     = ng3x3
          gci%ng3x3_restraint     = ng3x3_restraint
          ! Total number of Intermolecular constraints
          gci%ntot = 3*gci%ng3x3 + gci%ntot
          DEALLOCATE(constr_x_glob)
       END IF
    END IF

    CALL timestop(handle2)
    CALL timeset(routineN//"_5",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Set the group 4x6 constraint g4x6_list
    !-----------------------------------------------------------------------------
    IF(topology%const_46) THEN
       CALL give_constraint_array(cons_info%const_g46_mol,&
                                  cons_info%const_g46_molname,&
                                  cons_info%g46_intermolecular,&
                                  constr_x_mol,&
                                  constr_x_glob,&
                                  molecule_kind_set,&
                                  cons_info%g46_exclude_qm,&
                                  cons_info%g46_exclude_mm)
       ! Intramolecular constraints
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          ng4x6 = SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(g4x6_list(ng4x6))
          CALL setup_g4x6_list(g4x6_list, constr_x_mol(ii)%constr, cons_info, ng4x6_restraint)
          CALL set_molecule_kind(molecule_kind,ng4x6=ng4x6,ng4x6_restraint=ng4x6_restraint,g4x6_list=g4x6_list)
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first_atom, last_atom= last_atom)
             ALLOCATE(lg4x6(ng4x6))
             CALL setup_lg4x6(lg4x6, g4x6_list, first_atom, last_atom)
             CALL set_molecule(molecule=molecule,lg4x6=lg4x6)
          END DO
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr)
       END DO
       DEALLOCATE(constr_x_mol)
       ! Intermolecular constraints
       IF (ASSOCIATED(constr_x_glob)) THEN
          ng4x6 =  SIZE(constr_x_glob)
          ALLOCATE(g4x6_list(ng4x6))
          CALL setup_g4x6_list(g4x6_list, constr_x_glob, cons_info, ng4x6_restraint)
          ALLOCATE(lg4x6(ng4x6))
          CALL setup_lg4x6(lg4x6, g4x6_list, first_atom,last_atom)
          gci%g4x6_list => g4x6_list
          gci%lg4x6     => lg4x6
          gci%ng4x6     = ng4x6
          gci%ng4x6_restraint     = ng4x6_restraint
          ! Total number of Intermolecular constraints
          gci%ntot = 6*gci%ng4x6 + gci%ntot
          DEALLOCATE(constr_x_glob)
       END IF
    END IF

    CALL timestop(handle2)
    CALL timeset(routineN//"_6",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Set the group vsite constraint vsite_list
    !-----------------------------------------------------------------------------
    IF(topology%const_vsite) THEN
       CALL give_constraint_array(cons_info%const_vsite_mol,&
                                  cons_info%const_vsite_molname,&
                                  cons_info%vsite_intermolecular,&
                                  constr_x_mol,&
                                  constr_x_glob,&
                                  molecule_kind_set,&
                                  cons_info%vsite_exclude_qm,&
                                  cons_info%vsite_exclude_mm)
       ! Intramolecular constraints
       DO ii = 1, SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(ii)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list)
          nvsite = SIZE(constr_x_mol(ii)%constr)
          ALLOCATE(vsite_list(nvsite))
          CALL setup_vsite_list(vsite_list, constr_x_mol(ii)%constr, cons_info, nvsite_restraint)
          CALL set_molecule_kind(molecule_kind,nvsite=nvsite,nvsite_restraint=nvsite_restraint,&
               vsite_list=vsite_list)
       END DO
       DO j = 1, SIZE(constr_x_mol)
          DEALLOCATE(constr_x_mol(j)%constr)
       END DO
       DEALLOCATE(constr_x_mol)
       ! Intermolecular constraints
       IF (ASSOCIATED(constr_x_glob)) THEN
          nvsite =  SIZE(constr_x_glob)
          ALLOCATE(vsite_list(nvsite))
          CALL setup_vsite_list(vsite_list, constr_x_glob, cons_info, nvsite_restraint)
          gci%vsite_list => vsite_list
          gci%nvsite     = nvsite
          gci%nvsite_restraint     = nvsite_restraint
          ! Total number of Intermolecular constraints
          gci%ntot = gci%nvsite + gci%ntot
          DEALLOCATE(constr_x_glob)
       END IF
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_7",handle2)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Set the group fixed_atom constraint fixd_list
    !-----------------------------------------------------------------------------
    IF(topology%const_atom) THEN
       ALLOCATE(fixd_list_gci(SIZE(particle_set)))
       nfixd_list_gci = 0
       ALLOCATE(missed_molname(SIZE(cons_info%fixed_molnames,1)))
       missed_molname=.TRUE.
       nfixd_restart = 0
       DO i=1,SIZE(molecule_kind_set)
          molecule_kind => molecule_kind_set(i)
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               nmolecule=nmolecule,molecule_list=molecule_list, name=molname)
          is_qm = qmmm_ff_precond_only_qm(id1=molname)
          WHERE(molname.EQ.cons_info%fixed_molnames)
            missed_molname=.FALSE.
          END WHERE
          ! Try to figure out how many atoms of the list belong to this molecule_kind
          nfixed_atoms    = 0
          DO j=1,nmolecule
             molecule => molecule_set(molecule_list(j))
             CALL get_molecule ( molecule, first_atom = first, last_atom = last )
             fix_atom_molname=.FALSE.
             IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                DO k = 1, SIZE(cons_info%fixed_molnames)
                   IF (cons_info%fixed_molnames(k).EQ.molname) THEN
                      fix_atom_molname=.TRUE.
                      IF (is_qm.AND.cons_info%fixed_exclude_qm(k))        fix_atom_molname=.FALSE.
                      IF ((.NOT.is_qm).AND.cons_info%fixed_exclude_mm(k)) fix_atom_molname=.FALSE.
                   END IF
                END DO
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
          ALLOCATE(fixd_list(nfixed_atoms))
          kk              = 0
          nfixd_restraint = 0
          IF (nfixed_atoms/=0) THEN
             DO j=1,nmolecule
                molecule => molecule_set(molecule_list(j))
                CALL get_molecule ( molecule, first_atom = first, last_atom = last )
                fix_atom_molname=.FALSE.
                IF (ASSOCIATED(cons_info%fixed_molnames)) THEN
                   DO k1loc = 1, SIZE(cons_info%fixed_molnames)
                      IF (cons_info%fixed_molnames(k1loc).EQ.molname) THEN
                         fix_atom_molname=.TRUE.
                         itype = cons_info%fixed_mol_type(k1loc)
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
                         itype = cons_info%fixed_type(k2loc)
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
                            fix_atom_qmmm = .TRUE.
                            fix_atom_qm   = .TRUE.
                            itype         = cons_info%freeze_qm_type
                         END IF
                      CASE(do_constr_molec)
                         IF (ANY(qmmm_env%qm_molecule_index == molecule_list(j))) THEN
                            fix_atom_qmmm = .TRUE.
                            fix_atom_qm   = .TRUE.
                            itype         = cons_info%freeze_qm_type
                         END IF
                      END SELECT
                      SELECT CASE(cons_info%freeze_mm)
                      CASE(do_constr_atomic)
                         IF (ALL(qmmm_env%qm_atom_index/=k)) THEN
                            fix_atom_qmmm = .TRUE.
                            fix_atom_mm   = .TRUE.
                            itype         = cons_info%freeze_mm_type
                         END IF
                      CASE(do_constr_molec)
                         IF (ALL(qmmm_env%qm_molecule_index /= molecule_list(j))) THEN
                            fix_atom_qmmm = .TRUE.
                            fix_atom_mm   = .TRUE.
                            itype         = cons_info%freeze_mm_type
                         END IF
                      END SELECT
                      ! We should never reach this point but let's check it anyway
                      IF (fix_atom_qm.AND.fix_atom_mm) THEN
                         CALL stop_program(routineN,moduleN,__LINE__,&
                                           "Atom number: "//cp_to_string(k)//&
                                          " has been defined both QM and MM. General Error!")
                      END IF
                   END IF
                   ! Check that the fixed atom constraint/restraint is unique
                   IF ((fix_fixed_atom.AND.fix_atom_qmmm).OR.(fix_fixed_atom.AND.fix_atom_molname)&
                        .OR.(fix_atom_qmmm.AND.fix_atom_molname)) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                                        "Atom number: "//cp_to_string(k)//&
                                        " has been constrained/restrained to be fixed in more than one"//&
                                        " input section. Check and correct your input file!")
                   END IF
                   ! Let's store the atom index
                   IF (fix_fixed_atom.OR.fix_atom_qmmm.OR.fix_atom_molname) THEN
                      kk = kk + 1
                      fixd_list(kk)%fixd  = k
                      fixd_list(kk)%coord = particle_set(k)%r
                      fixd_list(kk)%itype = itype
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
                         CPASSERT(.FALSE.)
                      END IF
                      IF (fixd_list(kk)%restraint%active) THEN
                         nfixd_restraint = nfixd_restraint+ 1
                         nfixd_restart   = nfixd_restart  + 1
                         ! Check that we use the components that we really want..
                         SELECT CASE(itype)
                         CASE(use_perd_x)
                            fixd_list(kk)%coord(2) = HUGE(0.0_dp)
                            fixd_list(kk)%coord(3) = HUGE(0.0_dp)
                         CASE(use_perd_y)
                            fixd_list(kk)%coord(1) = HUGE(0.0_dp)
                            fixd_list(kk)%coord(3) = HUGE(0.0_dp)
                         CASE(use_perd_z)
                            fixd_list(kk)%coord(1) = HUGE(0.0_dp)
                            fixd_list(kk)%coord(2) = HUGE(0.0_dp)
                         CASE(use_perd_xy)
                            fixd_list(kk)%coord(3) = HUGE(0.0_dp)
                         CASE(use_perd_xz)
                            fixd_list(kk)%coord(2) = HUGE(0.0_dp)
                         CASE(use_perd_yz)
                            fixd_list(kk)%coord(1) = HUGE(0.0_dp)
                         END SELECT
                         IF (restart_restraint_pos) THEN
                            ! Read  coord0 value for restraint
                            CALL section_vals_val_get(fixd_restr_rest,"_DEFAULT_KEYWORD_",&
                                 i_rep_val=nfixd_restart, r_vals=r)
                            SELECT CASE(itype)
                            CASE(use_perd_x)
                               CPASSERT(SIZE(r)==1)
                               fixd_list(kk)%coord(1) = r(1)
                            CASE(use_perd_y)
                               CPASSERT(SIZE(r)==1)
                               fixd_list(kk)%coord(2) = r(1)
                            CASE(use_perd_z)
                               CPASSERT(SIZE(r)==1)
                               fixd_list(kk)%coord(3) = r(1)
                            CASE(use_perd_xy)
                               CPASSERT(SIZE(r)==2)
                               fixd_list(kk)%coord(1) = r(1)
                               fixd_list(kk)%coord(2) = r(2)
                            CASE(use_perd_xz)
                               CPASSERT(SIZE(r)==2)
                               fixd_list(kk)%coord(1) = r(1)
                               fixd_list(kk)%coord(3) = r(2)
                            CASE(use_perd_yz)
                               CPASSERT(SIZE(r)==2)
                               fixd_list(kk)%coord(2) = r(1)
                               fixd_list(kk)%coord(3) = r(2)
                            CASE(use_perd_xyz)
                               CPASSERT(SIZE(r)==3)
                               fixd_list(kk)%coord(1:3) = r(1:3)
                            END SELECT
                         ELSE
                            ! Write coord0 value for restraint
                            SELECT CASE(itype)
                            CASE(use_perd_x)
                               ALLOCATE(r(1))
                               r(1) = fixd_list(kk)%coord(1)
                            CASE(use_perd_y)
                               ALLOCATE(r(1))
                               r(1) = fixd_list(kk)%coord(2)
                            CASE(use_perd_z)
                               ALLOCATE(r(1))
                               r(1) = fixd_list(kk)%coord(3)
                            CASE(use_perd_xy)
                               ALLOCATE(r(2))
                               r(1) = fixd_list(kk)%coord(1)
                               r(2) = fixd_list(kk)%coord(2)
                            CASE(use_perd_xz)
                               ALLOCATE(r(2))
                               r(1) = fixd_list(kk)%coord(1)
                               r(2) = fixd_list(kk)%coord(3)
                            CASE(use_perd_yz)
                               ALLOCATE(r(2))
                               r(1) = fixd_list(kk)%coord(1)
                               r(2) = fixd_list(kk)%coord(3)
                            CASE(use_perd_xyz)
                               ALLOCATE(r(3))
                               r(1:3) = fixd_list(kk)%coord(1:3)
                            END SELECT
                            CALL section_vals_val_set(fixd_restr_rest,"_DEFAULT_KEYWORD_",&
                                 i_rep_val=nfixd_restart, r_vals_ptr=r)
                         END IF
                      END IF
                   END IF
                END DO
             END DO
          END IF
          IF (iw>0) THEN
             WRITE(iw,*)"MOLECULE KIND:",i," NR. FIXED ATOMS:",SIZE(fixd_list(:)%fixd)," LIST::",fixd_list(:)%fixd
          END IF
          CALL set_molecule_kind(molecule_kind, nfixd=nfixed_atoms, nfixd_restraint=nfixd_restraint,&
               fixd_list=fixd_list)
          fixd_list_gci(nfixd_list_gci+1:nfixd_list_gci+nfixed_atoms) = fixd_list
          nfixd_list_gci = nfixd_list_gci + nfixed_atoms
       END DO
       IF (iw>0) THEN
          WRITE(iw,*)"TOTAL NUMBER OF FIXED ATOMS:",nfixd_list_gci
       END IF
       CPASSERT(COUNT(missed_molname)==0)
       DEALLOCATE(missed_molname)
       ! Intermolecular constraints
       IF (gci%ntot /=0) THEN
          ALLOCATE(fixd_list(nfixd_list_gci))
          fixd_list(1:nfixd_list_gci)=fixd_list_gci(1:nfixd_list_gci)
          gci%fixd_list=>fixd_list
       END IF
       DEALLOCATE(fixd_list_gci)
    END IF
    ! Final setup of the number of possible restraints
    gci%nrestraint = gci%ng3x3_restraint +&
                     gci%ng4x6_restraint +&
                     gci%nvsite_restraint +&
                     gci%ncolv%nrestraint
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO")
    CALL timestop(handle2)
    CALL timestop(handle)
  END SUBROUTINE topology_constraint_pack

! *****************************************************************************
!> \brief Setup the colv_list for the packing of constraints
!> \param colv_list ...
!> \param ilist ...
!> \param gind ...
!> \param cons_info ...
!> \param topology ...
!> \param particle_set ...
!> \param restart_restraint_clv ...
!> \param colvar_rest ...
!> \param first_atom ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_colv_list(colv_list, ilist, gind, cons_info, topology,&
       particle_set, restart_restraint_clv, colvar_rest, first_atom)

    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    INTEGER, DIMENSION(:), POINTER           :: ilist
    INTEGER, INTENT(INOUT)                   :: gind
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    LOGICAL, INTENT(IN)                      :: restart_restraint_clv
    TYPE(section_vals_type), POINTER         :: colvar_rest
    INTEGER, INTENT(IN)                      :: first_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_colv_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, kdim, kk, ncolv_mol
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rmod
    TYPE(colvar_type), POINTER               :: local_colvar

    failure = .FALSE.
    ncolv_mol = 0
    DO kk = 1, SIZE(ilist)
       j = ilist(kk)
       ncolv_mol = ncolv_mol + 1
       kdim = SIZE(cons_info%colvar_set(j)%colvar%i_atom)
       ALLOCATE(colv_list(ncolv_mol)%i_atoms(kdim))
       colv_list(ncolv_mol)%inp_seq_num    = j
       colv_list(ncolv_mol)%type_id        = cons_info%colvar_set(j)%colvar%type_id
       colv_list(ncolv_mol)%i_atoms        = cons_info%colvar_set(j)%colvar%i_atom
       colv_list(ncolv_mol)%use_points     = cons_info%colvar_set(j)%colvar%use_points
       ! Restraint
       colv_list(ncolv_mol)%restraint%active = cons_info%colv_restraint(j)
       colv_list(ncolv_mol)%restraint%k0     = cons_info%colv_k0(j)
       IF (cons_info%const_colv_target(j)==-HUGE(0.0_dp)) THEN
          ! Let's compute the value..
          NULLIFY(local_colvar)
          CALL colvar_clone(local_colvar, cons_info%colvar_set(j)%colvar,&
               i_atom_offset=first_atom-1)
          CALL colvar_eval_mol_f(local_colvar, topology%cell, particle_set)
          colv_list(ncolv_mol)%expected_value = local_colvar%ss
          CALL colvar_release(local_colvar)
       ELSE
          colv_list(ncolv_mol)%expected_value = cons_info%const_colv_target(j)
       END IF
       colv_list(ncolv_mol)%expected_value_growth_speed = cons_info%const_colv_target_growth(j)
       ! In case of Restraint let's check for possible restart values
       IF (colv_list(ncolv_mol)%restraint%active.AND.&
           (colv_list(ncolv_mol)%expected_value_growth_speed == 0.0_dp)) THEN
          gind = gind + 1
          IF (restart_restraint_clv) THEN
             CALL section_vals_val_get(colvar_rest,"_DEFAULT_KEYWORD_",&
                  i_rep_val=gind, r_val=rmod)
             colv_list(ncolv_mol)%expected_value = rmod
          ELSE
             rmod = colv_list(ncolv_mol)%expected_value
             CALL section_vals_val_set(colvar_rest,"_DEFAULT_KEYWORD_",&
                  i_rep_val=gind, r_val=rmod)
          END IF
       END IF
       ! Only if torsion let's take into account the singularity in the definition
       ! of the dihedral
       IF (cons_info%colvar_set(j)%colvar%type_id==torsion_colvar_id) THEN
          cons_info%colvar_set(j)%colvar%torsion_param%o0=colv_list(ncolv_mol)%expected_value
       END IF
    END DO
  END SUBROUTINE setup_colv_list

! *****************************************************************************
!> \brief Setup the g3x3_list for the packing of constraints
!> \param g3x3_list ...
!> \param ilist ...
!> \param cons_info ...
!> \param ng3x3_restraint ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_g3x3_list(g3x3_list, ilist, cons_info, ng3x3_restraint)
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    INTEGER, DIMENSION(:), POINTER           :: ilist
    TYPE(constraint_info_type), POINTER      :: cons_info
    INTEGER, INTENT(OUT)                     :: ng3x3_restraint

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_g3x3_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ng3x3
    LOGICAL                                  :: failure

    failure = .FALSE.
    ng3x3_restraint = 0
    DO ng3x3 = 1, SIZE(ilist)
       j = ilist(ng3x3)
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

  END SUBROUTINE setup_g3x3_list

! *****************************************************************************
!> \brief Setup the g4x6_list for the packing of constraints
!> \param g4x6_list ...
!> \param ilist ...
!> \param cons_info ...
!> \param ng4x6_restraint ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_g4x6_list(g4x6_list, ilist, cons_info, ng4x6_restraint)
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    INTEGER, DIMENSION(:), POINTER           :: ilist
    TYPE(constraint_info_type), POINTER      :: cons_info
    INTEGER, INTENT(OUT)                     :: ng4x6_restraint

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_g4x6_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ng4x6
    LOGICAL                                  :: failure

    failure = .FALSE.
    ng4x6 = 0
    ng4x6_restraint = 0
    DO ng4x6 = 1, SIZE(ilist)
       j = ilist(ng4x6)
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

  END SUBROUTINE setup_g4x6_list

! *****************************************************************************
!> \brief Setup the vsite_list for the packing of constraints
!> \param vsite_list ...
!> \param ilist ...
!> \param cons_info ...
!> \param nvsite_restraint ...
!> \par History
!> \author Marcel Baer [2008]
! *****************************************************************************
  SUBROUTINE setup_vsite_list(vsite_list, ilist, cons_info, nvsite_restraint)
    TYPE(vsite_constraint_type), &
      DIMENSION(:), POINTER                  :: vsite_list
    INTEGER, DIMENSION(:), POINTER           :: ilist
    TYPE(constraint_info_type), POINTER      :: cons_info
    INTEGER, INTENT(OUT)                     :: nvsite_restraint

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_vsite_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, nvsite
    LOGICAL                                  :: failure

    failure = .FALSE.
    nvsite = 0
    nvsite_restraint = 0
    DO nvsite = 1, SIZE(ilist)
       j = ilist(nvsite)
       vsite_list(nvsite)%a   = cons_info%const_vsite_a(j)
       vsite_list(nvsite)%b   = cons_info%const_vsite_b(j)
       vsite_list(nvsite)%c   = cons_info%const_vsite_c(j)
       vsite_list(nvsite)%d   = cons_info%const_vsite_d(j)
       vsite_list(nvsite)%wbc = cons_info%const_vsite_wbc(j)
       vsite_list(nvsite)%wdc = cons_info%const_vsite_wdc(j)
       ! Restraint
       vsite_list(nvsite)%restraint%active = cons_info%vsite_restraint(j)
       vsite_list(nvsite)%restraint%k0     = cons_info%vsite_k0(j)
       IF (vsite_list(nvsite)%restraint%active) nvsite_restraint = nvsite_restraint+ 1
    END DO

  END SUBROUTINE setup_vsite_list
! *****************************************************************************
!> \brief Setup the lcolv for the packing of constraints
!> \param lcolv ...
!> \param ilist ...
!> \param first_atom ...
!> \param last_atom ...
!> \param cons_info ...
!> \param particle_set ...
!> \param colvar_func_info ...
!> \param use_clv_info ...
!> \param cind ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_lcolv(lcolv, ilist, first_atom, last_atom, cons_info,&
                         particle_set, colvar_func_info, use_clv_info,&
                         cind)
    TYPE(local_colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: lcolv
    INTEGER, DIMENSION(:), POINTER           :: ilist
    INTEGER, INTENT(IN)                      :: first_atom, last_atom
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(section_vals_type), POINTER         :: colvar_func_info
    LOGICAL, INTENT(IN)                      :: use_clv_info
    INTEGER, INTENT(INOUT)                   :: cind

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_lcolv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind, k, kk
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_vals

    failure = .FALSE.
    DO kk=1,SIZE(ilist)
       k = ilist(kk)
       lcolv(kk)%init      = .FALSE.
       lcolv(kk)%lambda    = 0.0_dp
       lcolv(kk)%sigma     = 0.0_dp

       ! Set Up colvar variable
       NULLIFY(lcolv(kk)%colvar, lcolv(kk)%colvar_old)
       ! Colvar
       CALL colvar_clone(lcolv(kk)%colvar, cons_info%colvar_set(k)%colvar,&
            i_atom_offset=first_atom-1)

       ! Some COLVARS may need additional information for evaluating the
       ! functional form: this is the case for COLVARS which depend on the
       ! initial position of the atoms: This information is stored in a proper
       ! container in the COLVAR_RESTART section..
       IF  ((lcolv(kk)%colvar%type_id==xyz_diag_colvar_id).OR.&
            (lcolv(kk)%colvar%type_id==xyz_outerdiag_colvar_id)) THEN
          cind = cind + 1
          IF (use_clv_info) THEN
             CALL section_vals_val_get(colvar_func_info,"_DEFAULT_KEYWORD_",&
                  i_rep_val=cind, r_vals=r_vals)
             SELECT CASE(lcolv(kk)%colvar%type_id)
             CASE(xyz_diag_colvar_id)
                CPASSERT(SIZE(r_vals)==3)
                lcolv(kk)%colvar%xyz_diag_param%r0=r_vals
             CASE(xyz_outerdiag_colvar_id)
                CPASSERT(SIZE(r_vals)==6)
                lcolv(kk)%colvar%xyz_outerdiag_param%r0(:,1)=r_vals(1:3)
                lcolv(kk)%colvar%xyz_outerdiag_param%r0(:,2)=r_vals(4:6)
             END SELECT
          ELSE
             SELECT CASE(lcolv(kk)%colvar%type_id)
             CASE(xyz_diag_colvar_id)
                ALLOCATE(r_vals(3))
                ind=first_atom-1+lcolv(kk)%colvar%xyz_diag_param%i_atom
                r_vals=particle_set(ind)%r
                lcolv(kk)%colvar%xyz_diag_param%r0=r_vals
             CASE(xyz_outerdiag_colvar_id)
                ALLOCATE(r_vals(6))
                ind=first_atom-1+lcolv(kk)%colvar%xyz_outerdiag_param%i_atoms(1)
                r_vals(1:3)=particle_set(ind)%r
                ind=first_atom-1+lcolv(kk)%colvar%xyz_outerdiag_param%i_atoms(2)
                r_vals(4:6)=particle_set(ind)%r
                lcolv(kk)%colvar%xyz_outerdiag_param%r0(:,1)=r_vals(1:3)
                lcolv(kk)%colvar%xyz_outerdiag_param%r0(:,2)=r_vals(4:6)
             END SELECT
             CALL section_vals_val_set(colvar_func_info,"_DEFAULT_KEYWORD_",&
                  i_rep_val=cind, r_vals_ptr=r_vals)
          END IF
       END IF

       ! Setup Colvar_old
       CALL colvar_clone(lcolv(kk)%colvar_old,lcolv(kk)%colvar)

       ! Check for consistency in the constraint definition
       IF  (ANY(lcolv(kk)%colvar%i_atom > last_atom).OR.&
            ANY(lcolv(kk)%colvar%i_atom < first_atom)) THEN
          WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
          WRITE(*,'(T5,"|",T8,A)')"A constraint has been defined for a molecule type",&
               " but the atoms specified in the constraint and the atoms defined for",&
               " the molecule DO NOT match!",&
               "This could be very probable due to a wrong connectivity, or an error",&
               " in the constraint specification in the input file.",&
               " Please check it carefully!"
          CPASSERT(.FALSE.)
       END IF
    END DO
  END SUBROUTINE setup_lcolv

! *****************************************************************************
!> \brief Setup the lg3x3 for the packing of constraints
!> \param lg3x3 ...
!> \param g3x3_list ...
!> \param first_atom ...
!> \param last_atom ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_lg3x3(lg3x3, g3x3_list, first_atom, last_atom)
    TYPE(local_g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: lg3x3
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    INTEGER, INTENT(IN)                      :: first_atom, last_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_lg3x3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: kk
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO kk = 1, SIZE(lg3x3)
       lg3x3(kk)%init      = .FALSE.
       lg3x3(kk)%scale     = 0.0_dp
       lg3x3(kk)%scale_old = 0.0_dp
       lg3x3(kk)%fa        = 0.0_dp
       lg3x3(kk)%fb        = 0.0_dp
       lg3x3(kk)%fc        = 0.0_dp
       lg3x3(kk)%ra_old    = 0.0_dp
       lg3x3(kk)%rb_old    = 0.0_dp
       lg3x3(kk)%rc_old    = 0.0_dp
       lg3x3(kk)%va        = 0.0_dp
       lg3x3(kk)%vb        = 0.0_dp
       lg3x3(kk)%vc        = 0.0_dp
       lg3x3(kk)%lambda    = 0.0_dp
       IF  ((g3x3_list(kk)%a+first_atom-1 < first_atom).OR.&
            (g3x3_list(kk)%b+first_atom-1 < first_atom).OR.&
            (g3x3_list(kk)%c+first_atom-1 < first_atom).OR.&
            (g3x3_list(kk)%a+first_atom-1 > last_atom ).OR.&
            (g3x3_list(kk)%b+first_atom-1 > last_atom ).OR.&
            (g3x3_list(kk)%c+first_atom-1 > last_atom )) THEN
          WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
          WRITE(*,'(T5,"|",T8,A)')"A constraint has been defined for a molecule type",&
               " but the atoms specified in the constraint and the atoms defined for",&
               " the molecule DO NOT match!",&
               "This could be very probable due to a wrong connectivity, or an error",&
               " in the constraint specification in the input file.",&
               " Please check it carefully!"
          CPASSERT(.FALSE.)
       END IF
    END DO

  END SUBROUTINE setup_lg3x3

! *****************************************************************************
!> \brief Setup the lg4x6 for the packing of constraints
!> \param lg4x6 ...
!> \param g4x6_list ...
!> \param first_atom ...
!> \param last_atom ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2007]
! *****************************************************************************
  SUBROUTINE setup_lg4x6(lg4x6, g4x6_list, first_atom, last_atom)
    TYPE(local_g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: lg4x6
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    INTEGER, INTENT(IN)                      :: first_atom, last_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_lg4x6', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: kk
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO kk = 1, SIZE(lg4x6)
       lg4x6(kk)%init      = .FALSE.
       lg4x6(kk)%scale     = 0.0_dp
       lg4x6(kk)%scale_old = 0.0_dp
       lg4x6(kk)%fa        = 0.0_dp
       lg4x6(kk)%fb        = 0.0_dp
       lg4x6(kk)%fc        = 0.0_dp
       lg4x6(kk)%fd        = 0.0_dp
       lg4x6(kk)%fe        = 0.0_dp
       lg4x6(kk)%ff        = 0.0_dp
       lg4x6(kk)%ra_old    = 0.0_dp
       lg4x6(kk)%rb_old    = 0.0_dp
       lg4x6(kk)%rc_old    = 0.0_dp
       lg4x6(kk)%rd_old    = 0.0_dp
       lg4x6(kk)%re_old    = 0.0_dp
       lg4x6(kk)%rf_old    = 0.0_dp
       lg4x6(kk)%va        = 0.0_dp
       lg4x6(kk)%vb        = 0.0_dp
       lg4x6(kk)%vc        = 0.0_dp
       lg4x6(kk)%vd        = 0.0_dp
       lg4x6(kk)%ve        = 0.0_dp
       lg4x6(kk)%vf        = 0.0_dp
       lg4x6(kk)%lambda    = 0.0_dp
       IF  ((g4x6_list(kk)%a+first_atom-1 < first_atom).OR.&
            (g4x6_list(kk)%b+first_atom-1 < first_atom).OR.&
            (g4x6_list(kk)%c+first_atom-1 < first_atom).OR.&
            (g4x6_list(kk)%d+first_atom-1 < first_atom).OR.&
            (g4x6_list(kk)%a+first_atom-1 > last_atom ).OR.&
            (g4x6_list(kk)%b+first_atom-1 > last_atom ).OR.&
            (g4x6_list(kk)%c+first_atom-1 > last_atom ).OR.&
            (g4x6_list(kk)%d+first_atom-1 > last_atom )) THEN
          WRITE(*,'(T5,"|",T8,A)')"Error in constraints setup!"
          WRITE(*,'(T5,"|",T8,A)')"A constrained has been defined for a molecule type",&
               " but the atoms specified in the constraint and the atoms defined for",&
               " the molecule DO NOT match!",&
               "This could be very probable due to a wrong connectivity, or an error",&
               " in the constraint specification in the input file.",&
               " Please check it carefully!"
          CPASSERT(.FALSE.)
       END IF
    END DO

  END SUBROUTINE setup_lg4x6

! *****************************************************************************
!> \brief Gives back a list of molecule to which apply the constraint
!> \param const_mol ...
!> \param const_molname ...
!> \param const_intermolecular ...
!> \param constr_x_mol ...
!> \param constr_x_glob ...
!> \param molecule_kind_set ...
!> \param exclude_qm ...
!> \param exclude_mm ...
!> \par History
!>      Updated 2007 for intermolecular constraints
!> \author Teodoro Laino [2006]
! *****************************************************************************
  SUBROUTINE give_constraint_array(const_mol, const_molname, const_intermolecular,&
       constr_x_mol, constr_x_glob, molecule_kind_set, exclude_qm, exclude_mm)

    INTEGER, DIMENSION(:), POINTER           :: const_mol
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: const_molname
    LOGICAL, DIMENSION(:), POINTER           :: const_intermolecular
    TYPE(constr_list_type), DIMENSION(:), &
      POINTER                                :: constr_x_mol
    INTEGER, DIMENSION(:), POINTER           :: constr_x_glob
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    LOGICAL, DIMENSION(:), POINTER           :: exclude_qm, exclude_mm

    CHARACTER(len=*), PARAMETER :: routineN = 'give_constraint_array', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: myname, name
    INTEGER                                  :: handle, i, iglob, isize, k
    LOGICAL                                  :: failure, found_molname, is_qm
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(molecule_kind)
    ALLOCATE(constr_x_mol(SIZE(molecule_kind_set)))
    DO i = 1, SIZE(constr_x_mol)
       NULLIFY(constr_x_mol(i)%constr)
       ALLOCATE(constr_x_mol(i)%constr(0))
    END DO
    CPASSERT(SIZE(const_mol)==SIZE(const_molname))
    iglob = 0
    DO i = 1, SIZE(const_mol)
       IF (const_intermolecular(i)) THEN
          ! Intermolecular constraint
          iglob = iglob + 1
          CALL reallocate(constr_x_glob, 1, iglob)
          constr_x_glob(iglob) = i
       ELSE
          ! Intramolecular constraint
          IF (const_mol(i)/=0) THEN
             k = const_mol(i)
             CALL cp_assert(k<=SIZE(molecule_kind_set),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "A constraint has been specified providing the molecule index. But the "//&
                  " molecule index ("//cp_to_string(k)//") is out of range of the possible"//&
                  " molecule kinds ("//cp_to_string(SIZE(molecule_kind_set))//")."//&
CPSourceFileRef,&
                  failure=failure)
             isize = SIZE(constr_x_mol(k)%constr)
             CALL reallocate(constr_x_mol(k)%constr, 1, isize+1)
             constr_x_mol(k)%constr(isize+1) = i
          ELSE
             myname = const_molname(i)
             found_molname = .FALSE.
             DO k = 1, SIZE(molecule_kind_set)
                molecule_kind => molecule_kind_set(k)
                name = molecule_kind%name
                is_qm = qmmm_ff_precond_only_qm(id1=name)
                IF (is_qm.AND.exclude_qm(i))      CYCLE
                IF (.NOT.is_qm.AND.exclude_mm(i)) CYCLE
                IF (name==myname) THEN
                   isize = SIZE(constr_x_mol(k)%constr)
                   CALL reallocate(constr_x_mol(k)%constr, 1, isize+1)
                   constr_x_mol(k)%constr(isize+1) = i
                   found_molname = .TRUE.
                END IF
             END DO
             CALL print_warning_molname(found_molname, myname)
          END IF
       END IF
    END DO
    CALL timestop(handle)
  END SUBROUTINE give_constraint_array

! *****************************************************************************
!> \brief Prints a warning message if undefined molnames are used to define constraints
!> \param found ...
!> \param name ...
!> \author Teodoro Laino [2007] - Zurich University
! *****************************************************************************
  SUBROUTINE print_warning_molname(found, name)
    LOGICAL, INTENT(IN)                      :: found
    CHARACTER(LEN=*), INTENT(IN)             :: name

    CHARACTER(len=*), PARAMETER :: routineN = 'print_warning_molname', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_assert(found,cp_warning_level,cp_assertion_failed,routineP,&
         " MOLNAME ("//TRIM(name)//") was defined for constraints, but this molecule name "//&
         "is not defined. Please check carefully your PDB, PSF (has priority over PDB) or "//&
         "input driven CP2K coordinates. In case you may not find the reason for this warning "//&
         "it may be a good idea to print all molecule information (including kind name) activating "//&
         "the print_key MOLECULES specific of the SUBSYS%PRINT section. "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

  END SUBROUTINE print_warning_molname

END MODULE topology_constraint_util

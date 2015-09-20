!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Reads the input sections "topology"
!> \par History
!>      JGH (26-01-2002) Added read_topology_section
!> \author JGH
! *****************************************************************************
MODULE topology_input
  USE colvar_types,                    ONLY: colvar_clone,&
                                             colvar_p_type
  USE input_constants,                 ONLY: do_conn_generate,&
                                             do_conn_mol_set,&
                                             do_conn_off,&
                                             do_conn_user,&
                                             do_constr_none,&
                                             do_coord_off
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_unset
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: constraint_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_input'

  PRIVATE
  PUBLIC :: read_topology_section, read_constraints_section

CONTAINS

! *****************************************************************************
!> \brief reads the input section topology
!> \param topology ...
!> \param topology_section ...
!> \par History
!>      none
!> \author JGH (26-01-2002)
! *****************************************************************************
SUBROUTINE read_topology_section ( topology, topology_section)
    TYPE(topology_parameters_type)           :: topology
    TYPE(section_vals_type), POINTER         :: topology_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ival
    LOGICAL                                  :: failure

  CALL timeset(routineN,handle)
  failure = .FALSE.
  CALL section_vals_val_get(topology_section,"CHARGE_OCCUP",l_val=topology%charge_occup)
  CALL section_vals_val_get(topology_section,"CHARGE_BETA",l_val=topology%charge_beta)
  CALL section_vals_val_get(topology_section,"CHARGE_EXTENDED",l_val=topology%charge_extended)
  ival = COUNT((/topology%charge_occup,topology%charge_beta,topology%charge_extended/))
  CALL cp_assert(ival<=1,cp_failure_level,cp_assertion_failed,routineP,&
       "Only one between <CHARGE_OCCUP,CHARGE_BETA,CHARGE_EXTENDED> can be defined! "//&
CPSourceFileRef,&
       only_ionode=.TRUE.)
  CALL section_vals_val_get(topology_section,"PARA_RES",l_val=topology%para_res)
  CALL section_vals_val_get(topology_section,"GENERATE%REORDER",l_val=topology%reorder_atom)
  CALL section_vals_val_get(topology_section,"GENERATE%CREATE_MOLECULES",l_val=topology%create_molecules)
  CALL section_vals_val_get(topology_section,"MOL_CHECK",l_val=topology%molecules_check)
  CALL section_vals_val_get(topology_section,"USE_G96_VELOCITY",l_val=topology%use_g96_velocity)
  CALL section_vals_val_get(topology_section,"COORD_FILE_FORMAT",i_val=topology%coord_type)
  SELECT CASE(topology%coord_type)
  CASE(do_coord_off)
     ! Do Nothing
  CASE DEFAULT
     topology%coordinate = .TRUE.
     CALL section_vals_val_get(topology_section,"COORD_FILE_NAME",c_val=topology%coord_file_name)
  END SELECT
  CALL section_vals_val_get(topology_section,"CONN_FILE_FORMAT",i_val=topology%conn_type)
  SELECT CASE(topology%conn_type)
  CASE(do_conn_off,do_conn_generate,do_conn_mol_set,do_conn_user)
     ! Do Nothing
  CASE DEFAULT
     CALL section_vals_val_get(topology_section,"CONN_FILE_NAME",c_val=topology%conn_file_name)
  END SELECT
  CALL section_vals_val_get(topology_section,"EXCLUDE_VDW",i_val=topology%exclude_vdw)
  CALL section_vals_val_get(topology_section,"EXCLUDE_EI",i_val=topology%exclude_ei)
  CALL section_vals_val_get(topology_section,"GENERATE%BONDPARM",i_val=topology%bondparm_type)
  CALL section_vals_val_get(topology_section,"GENERATE%BONDPARM_FACTOR",r_val=topology%bondparm_factor)
  CALL timestop(handle)
END SUBROUTINE read_topology_section

! *****************************************************************************
!> \brief Read all the distance parameters. Put them in the
!>      constraint_distance array.
!> \param topology ...
!> \param colvar_p ...
!> \param constraint_section ...
!> \par History
!>      JGH (26-01-2002) Distance parameters are now stored in tables. The position
!>         within the table is used as handle for the topology
!>      teo Read the CONSTRAINT section within the new input style
!> \author teo
! *****************************************************************************
SUBROUTINE read_constraints_section(topology, colvar_p, constraint_section)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(section_vals_type), POINTER         :: constraint_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_constraints_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: icolvar, ig, isize, &
                                                isize_old, itype, jg, msize, &
                                                msize_old, n_rep, ncons, &
                                                nrep, stat
    INTEGER, DIMENSION(:), POINTER           :: ilist, tmplist
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rlist
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(section_vals_type), POINTER :: collective_section, fix_atom_section, &
      g3x3_section, g4x6_section, hbonds_section, vsite_section

  cons_info => topology%cons_info
  failure   = .FALSE.
  IF (ASSOCIATED(constraint_section)) THEN
     hbonds_section     => section_vals_get_subs_vals(constraint_section,"HBONDS")
     g3x3_section       => section_vals_get_subs_vals(constraint_section,"G3X3")
     g4x6_section       => section_vals_get_subs_vals(constraint_section,"G4X6")
     vsite_section      => section_vals_get_subs_vals(constraint_section,"VIRTUAL_SITE")
     fix_atom_section   => section_vals_get_subs_vals(constraint_section,"FIXED_ATOMS")
     collective_section => section_vals_get_subs_vals(constraint_section,"COLLECTIVE")
     ! HBONDS
     CALL section_vals_get(hbonds_section, explicit=topology%const_hydr)
     CALL check_restraint(hbonds_section,&
                          is_restraint=cons_info%hbonds_restraint, &
                          k0=cons_info%hbonds_k0,&
                          label="HBONDS")
     ! G3X3
     CALL section_vals_get(g3x3_section,explicit=explicit, n_repetition=ncons)
     IF (explicit) THEN
        topology%const_33    = .TRUE.
        cons_info%nconst_g33 = ncons
        !
        ALLOCATE(cons_info%const_g33_mol(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_molname(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_a(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_b(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_c(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_dab(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_dac(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g33_dbc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g33_intermolecular(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g33_restraint(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g33_k0(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g33_exclude_qm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g33_exclude_mm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        DO ig = 1, ncons
           CALL check_restraint(g3x3_section,&
                                is_restraint=cons_info%g33_restraint(ig), &
                                k0=cons_info%g33_k0(ig),&
                                i_rep_section=ig,&
                                label="G3X3")
           cons_info%const_g33_mol(ig)     = 0
           cons_info%const_g33_molname(ig) = "UNDEF"
           ! Exclude QM or MM
           CALL section_vals_val_get(g3x3_section,"EXCLUDE_QM",i_rep_section=ig,&
                l_val=cons_info%g33_exclude_qm(ig))
           CALL section_vals_val_get(g3x3_section,"EXCLUDE_MM",i_rep_section=ig,&
                l_val=cons_info%g33_exclude_mm(ig))
           ! Intramolecular restraint
           CALL section_vals_val_get(g3x3_section,"INTERMOLECULAR",i_rep_section=ig,&
                l_val=cons_info%g33_intermolecular(ig))
           ! If it is intramolecular let's unset (in case user did it)
           ! the molecule and molname field
           IF (cons_info%g33_intermolecular(ig)) THEN
              CALL section_vals_val_unset(g3x3_section,"MOLECULE",i_rep_section=ig)
              CALL section_vals_val_unset(g3x3_section,"MOLNAME",i_rep_section=ig)
           END IF
           ! Let's tag to which molecule we want to apply constraints
           CALL section_vals_val_get(g3x3_section,"MOLECULE",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(g3x3_section,"MOLECULE",i_rep_section=ig,&
                   i_val=cons_info%const_g33_mol(ig))
           END IF
           CALL section_vals_val_get(g3x3_section,"MOLNAME",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(g3x3_section,"MOLNAME",i_rep_section=ig,&
                   c_val=cons_info%const_g33_molname(ig))
           END IF
           IF  ((cons_info%const_g33_mol(ig)/=0).AND.(cons_info%const_g33_molname(ig)/="UNDEF")) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           IF  ((cons_info%const_g33_mol(ig)==0).AND.(cons_info%const_g33_molname(ig)=="UNDEF").AND.&
                (.NOT.cons_info%g33_intermolecular(ig))) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           CALL section_vals_val_get(g3x3_section,"ATOMS",i_rep_section=ig,&
                i_vals=ilist)
           CALL section_vals_val_get(g3x3_section,"DISTANCES",i_rep_section=ig,&
                r_vals=rlist)
           cons_info%const_g33_a(ig) = ilist(1)
           cons_info%const_g33_b(ig) = ilist(2)
           cons_info%const_g33_c(ig) = ilist(3)

           cons_info%const_g33_dab(ig) = rlist(1)
           cons_info%const_g33_dac(ig) = rlist(2)
           cons_info%const_g33_dbc(ig) = rlist(3)
        END DO
     END IF
     ! G4X6
     CALL section_vals_get(g4x6_section,explicit=explicit, n_repetition=ncons)
     IF (explicit) THEN
        topology%const_46   =.TRUE.
        cons_info%nconst_g46 = ncons
        !
        ALLOCATE(cons_info%const_g46_mol(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_molname(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_a(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_b(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_c(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_d(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dab(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dac(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dbc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dad(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dbd(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_g46_dcd(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g46_intermolecular(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g46_restraint(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g46_k0(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g46_exclude_qm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%g46_exclude_mm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        DO ig = 1, ncons
           CALL check_restraint(g4x6_section,&
                                is_restraint=cons_info%g46_restraint(ig), &
                                k0=cons_info%g46_k0(ig),&
                                i_rep_section=ig,&
                                label="G4X6")
           cons_info%const_g46_mol(ig)     = 0
           cons_info%const_g46_molname(ig) = "UNDEF"
           ! Exclude QM or MM
           CALL section_vals_val_get(g4x6_section,"EXCLUDE_QM",i_rep_section=ig,&
                l_val=cons_info%g46_exclude_qm(ig))
           CALL section_vals_val_get(g4x6_section,"EXCLUDE_MM",i_rep_section=ig,&
                l_val=cons_info%g46_exclude_mm(ig))
           ! Intramolecular restraint
           CALL section_vals_val_get(g4x6_section,"INTERMOLECULAR",i_rep_section=ig,&
                l_val=cons_info%g46_intermolecular(ig))
           ! If it is intramolecular let's unset (in case user did it)
           ! the molecule and molname field
           IF (cons_info%g46_intermolecular(ig)) THEN
              CALL section_vals_val_unset(g4x6_section,"MOLECULE",i_rep_section=ig)
              CALL section_vals_val_unset(g4x6_section,"MOLNAME",i_rep_section=ig)
           END IF
           ! Let's tag to which molecule we want to apply constraints
           CALL section_vals_val_get(g4x6_section,"MOLECULE",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(g4x6_section,"MOLECULE",i_rep_section=ig,&
                   i_val=cons_info%const_g46_mol(ig))
           END IF
           CALL section_vals_val_get(g4x6_section,"MOLNAME",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(g4x6_section,"MOLNAME",i_rep_section=ig,&
                   c_val=cons_info%const_g46_molname(ig))
           END IF
           IF  ((cons_info%const_g46_mol(ig)/=0).AND.(cons_info%const_g46_molname(ig)/="UNDEF")) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           IF  ((cons_info%const_g46_mol(ig)==0).AND.(cons_info%const_g46_molname(ig)=="UNDEF").AND.&
                (.NOT.cons_info%g46_intermolecular(ig))) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           CALL section_vals_val_get(g4x6_section,"ATOMS",i_rep_section=ig,&
                i_vals=ilist)
           CALL section_vals_val_get(g4x6_section,"DISTANCES",i_rep_section=ig,&
                r_vals=rlist)
           cons_info%const_g46_a(ig)   = ilist(1)
           cons_info%const_g46_b(ig)   = ilist(2)
           cons_info%const_g46_c(ig)   = ilist(3)
           cons_info%const_g46_d(ig)   = ilist(4)
           cons_info%const_g46_dab(ig) = rlist(1)
           cons_info%const_g46_dac(ig) = rlist(2)
           cons_info%const_g46_dad(ig) = rlist(3)
           cons_info%const_g46_dbc(ig) = rlist(4)
           cons_info%const_g46_dbd(ig) = rlist(5)
           cons_info%const_g46_dcd(ig) = rlist(6)
        END DO
     END IF
     ! virtual
     CALL section_vals_get(vsite_section,explicit=explicit, n_repetition=ncons)
     IF (explicit) THEN
        topology%const_vsite   =.TRUE.
        cons_info%nconst_vsite = ncons
        !
        ALLOCATE(cons_info%const_vsite_mol(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_molname(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_a(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_b(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_c(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_d(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_wbc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_vsite_wdc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%vsite_intermolecular(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%vsite_restraint(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%vsite_k0(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%vsite_exclude_qm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%vsite_exclude_mm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        DO ig = 1, ncons
           CALL check_restraint(vsite_section,&
                                is_restraint=cons_info%vsite_restraint(ig), &
                                k0=cons_info%vsite_k0(ig),&
                                i_rep_section=ig,&
                                label="Virtual_SITE")
           cons_info%const_vsite_mol(ig)     = 0
           cons_info%const_vsite_molname(ig) = "UNDEF"
           ! Exclude QM or MM
           CALL section_vals_val_get(vsite_section,"EXCLUDE_QM",i_rep_section=ig,&
                l_val=cons_info%vsite_exclude_qm(ig))
           CALL section_vals_val_get(vsite_section,"EXCLUDE_MM",i_rep_section=ig,&
                l_val=cons_info%vsite_exclude_mm(ig))
           ! Intramolecular restraint
           CALL section_vals_val_get(vsite_section,"INTERMOLECULAR",i_rep_section=ig,&
                l_val=cons_info%vsite_intermolecular(ig))
           ! If it is intramolecular let's unset (in case user did it)
           ! the molecule and molname field
           IF (cons_info%vsite_intermolecular(ig)) THEN
              CALL section_vals_val_unset(vsite_section,"MOLECULE",i_rep_section=ig)
              CALL section_vals_val_unset(vsite_section,"MOLNAME",i_rep_section=ig)
           END IF
           ! Let's tag to which molecule we want to apply constraints
           CALL section_vals_val_get(vsite_section,"MOLECULE",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(vsite_section,"MOLECULE",i_rep_section=ig,&
                   i_val=cons_info%const_vsite_mol(ig))
           END IF
           CALL section_vals_val_get(vsite_section,"MOLNAME",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(vsite_section,"MOLNAME",i_rep_section=ig,&
                   c_val=cons_info%const_vsite_molname(ig))
           END IF
           IF  ((cons_info%const_vsite_mol(ig)/=0).AND.(cons_info%const_vsite_molname(ig)/="UNDEF")) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           IF  ((cons_info%const_vsite_mol(ig)==0).AND.(cons_info%const_vsite_molname(ig)=="UNDEF").AND.&
                (.NOT.cons_info%vsite_intermolecular(ig))) THEN
              CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
           END IF
           CALL section_vals_val_get(vsite_section,"ATOMS",i_rep_section=ig,&
                i_vals=ilist)
           CALL section_vals_val_get(vsite_section,"PARAMETERS",i_rep_section=ig,&
                r_vals=rlist)
           cons_info%const_vsite_a(ig)   = ilist(1)
           cons_info%const_vsite_b(ig)   = ilist(2)
           cons_info%const_vsite_c(ig)   = ilist(3)
           cons_info%const_vsite_d(ig)   = ilist(4)
           cons_info%const_vsite_wbc(ig) = rlist(1)
           cons_info%const_vsite_wdc(ig) = rlist(2)
        END DO
     END IF
     ! FIXED ATOMS
     CALL section_vals_get(fix_atom_section,explicit=explicit, n_repetition=ncons)
     IF (explicit) THEN
        NULLIFY(tmplist,tmpstringlist)
        isize  = 0
        msize  = 0
        ALLOCATE(cons_info%fixed_atoms(isize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_type(isize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_restraint(isize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_k0(isize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_molnames(msize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_mol_type(isize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_mol_restraint(msize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_mol_k0(msize),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_exclude_qm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%fixed_exclude_mm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        DO ig = 1, ncons
           isize_old  = isize
           msize_old  = msize
           CALL section_vals_val_get(fix_atom_section,"COMPONENTS_TO_FIX", i_rep_section=ig,&
                i_val=itype)
           CALL section_vals_val_get(fix_atom_section,"LIST", i_rep_section=ig,&
                n_rep_val=n_rep)
           DO jg = 1, n_rep
              CALL section_vals_val_get(fix_atom_section,"LIST", i_rep_section=ig, &
                   i_rep_val=jg,i_vals=tmplist)
              CALL reallocate(cons_info%fixed_atoms,1,isize+SIZE(tmplist))
              cons_info%fixed_atoms(isize+1:isize+SIZE(tmplist)) = tmplist
              CALL reallocate(cons_info%fixed_restraint,1,isize+SIZE(tmplist))
              CALL reallocate(cons_info%fixed_k0,1,isize+SIZE(tmplist))
              CALL reallocate(cons_info%fixed_type,1,isize+SIZE(tmplist))
              cons_info%fixed_type(isize+1:isize+SIZE(tmplist)) = itype
              isize = SIZE(cons_info%fixed_atoms)
           END DO
           !Check for restraints
           IF ((isize-isize_old)>0) THEN
              CALL check_restraint(fix_atom_section,&
                                   is_restraint=cons_info%fixed_restraint(isize_old+1), &
                                   k0=cons_info%fixed_k0(isize_old+1),&
                                   i_rep_section=ig,&
                                   label="FIXED ATOM")
              cons_info%fixed_restraint(isize_old+1:isize)=cons_info%fixed_restraint(isize_old+1)
              cons_info%fixed_k0(isize_old+1:isize)=cons_info%fixed_k0(isize_old+1)
           END IF
           CALL section_vals_val_get(fix_atom_section,"MOLNAME", i_rep_section=ig,&
                n_rep_val=n_rep)
           IF (n_rep /= 0) THEN
              DO jg = 1, n_rep
                 CALL section_vals_val_get(fix_atom_section,"MOLNAME", i_rep_section=ig, &
                      i_rep_val=jg, c_vals=tmpstringlist)
                 CALL reallocate(cons_info%fixed_molnames,1,msize+SIZE(tmpstringlist,1))
                 CALL reallocate(cons_info%fixed_mol_type,1,msize+SIZE(tmpstringlist,1))
                 CALL reallocate(cons_info%fixed_mol_restraint,1,msize+SIZE(tmpstringlist,1))
                 CALL reallocate(cons_info%fixed_mol_k0,1,msize+SIZE(tmpstringlist,1))
                 cons_info%fixed_molnames(msize+1:msize+SIZE(tmpstringlist,1))=tmpstringlist
                 cons_info%fixed_mol_type(msize+1:msize+SIZE(tmpstringlist,1))=itype
                 msize  = SIZE(cons_info%fixed_molnames)
              END DO
              ! Exclude QM or MM work only if defined MOLNAME
              CALL reallocate(cons_info%fixed_exclude_qm,1,msize)
              CALL reallocate(cons_info%fixed_exclude_mm,1,msize)
              CALL section_vals_val_get(fix_atom_section,"EXCLUDE_QM",i_rep_section=ig,&
                   l_val=cons_info%fixed_exclude_qm(msize_old+1))
              CALL section_vals_val_get(fix_atom_section,"EXCLUDE_MM",i_rep_section=ig,&
                   l_val=cons_info%fixed_exclude_mm(msize_old+1))
              cons_info%fixed_exclude_qm(msize_old+1:msize) = cons_info%fixed_exclude_qm(msize_old+1)
              cons_info%fixed_exclude_mm(msize_old+1:msize) = cons_info%fixed_exclude_mm(msize_old+1)
           END IF
           !Check for restraints
           IF (n_rep /=0) THEN
              CALL check_restraint(fix_atom_section,&
                                   is_restraint=cons_info%fixed_mol_restraint(msize_old+1), &
                                   k0=cons_info%fixed_mol_k0(msize_old+1),&
                                   i_rep_section=ig,&
                                   label="FIXED ATOM")
              cons_info%fixed_mol_restraint(msize_old+1:msize)=cons_info%fixed_mol_restraint(msize_old+1)
              cons_info%fixed_mol_k0(msize_old+1:msize)=cons_info%fixed_mol_k0(msize_old+1)
           END IF
           CALL section_vals_val_get(fix_atom_section,"MM_SUBSYS",i_rep_section=ig,&
                                     n_rep_val=nrep,explicit=explicit)
           IF (nrep==1 .AND. explicit) THEN
              CPPostcondition(cons_info%freeze_mm==do_constr_none,cp_failure_level,routineP,failure)
              CALL section_vals_val_get(fix_atom_section,"MM_SUBSYS",i_val=cons_info%freeze_mm,&
                   i_rep_section=ig)
              cons_info%freeze_mm_type = itype
           END IF
           CALL section_vals_val_get(fix_atom_section,"QM_SUBSYS",i_rep_section=ig,&
                n_rep_val=nrep,explicit=explicit)
           IF (nrep==1 .AND. explicit) THEN
              CPPostcondition(cons_info%freeze_qm==do_constr_none,cp_failure_level,routineP,failure)
              CALL section_vals_val_get(fix_atom_section,"QM_SUBSYS",i_val=cons_info%freeze_qm,&
                   i_rep_section=ig)
              cons_info%freeze_qm_type = itype
           END IF
           IF (cons_info%freeze_mm/=do_constr_none) THEN
              CALL check_restraint(fix_atom_section,&
                                   is_restraint=cons_info%fixed_mm_restraint, &
                                   k0=cons_info%fixed_mm_k0,&
                                   i_rep_section=ig,&
                                   label="FIXED ATOM")
           END IF
           IF (cons_info%freeze_qm/=do_constr_none) THEN
              CALL check_restraint(fix_atom_section,&
                                   is_restraint=cons_info%fixed_qm_restraint, &
                                   k0=cons_info%fixed_qm_k0,&
                                   i_rep_section=ig,&
                                   label="FIXED ATOM")
           END IF

        END DO
        IF ( (isize/=0).OR.(msize/=0).OR.&
             (cons_info%freeze_mm/=do_constr_none).OR.&
             (cons_info%freeze_qm/=do_constr_none) ) THEN
           topology%const_atom    = .TRUE.
        END IF
     END IF
     ! Collective Constraints
     CALL section_vals_get(collective_section,explicit=explicit, n_repetition=ncons)
     IF (explicit) THEN
        topology%const_colv = .TRUE.
        DO ig = 1, ncons
           CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig, i_val=icolvar)
           CPPostcondition(icolvar<=SIZE(colvar_p),cp_failure_level,routineP,failure)
        END DO
        cons_info%nconst_colv = ncons
        ALLOCATE(cons_info%const_colv_mol(ncons), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_colv_molname(ncons), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_colv_target(ncons), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%const_colv_target_growth(ncons), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colvar_set(ncons), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colv_intermolecular(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colv_restraint(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colv_k0(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colv_exclude_qm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        ALLOCATE(cons_info%colv_exclude_mm(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        DO ig = 1, ncons
           CALL check_restraint(collective_section,&
                                is_restraint=cons_info%colv_restraint(ig), &
                                k0=cons_info%colv_k0(ig),&
                                i_rep_section=ig,&
                                label="COLLECTIVE")
           cons_info%const_colv_mol(ig)     = 0
           cons_info%const_colv_molname(ig) = "UNDEF"
           ! Exclude QM or MM
           CALL section_vals_val_get(collective_section,"EXCLUDE_QM",i_rep_section=ig,&
                l_val=cons_info%colv_exclude_qm(ig))
           CALL section_vals_val_get(collective_section,"EXCLUDE_MM",i_rep_section=ig,&
                l_val=cons_info%colv_exclude_mm(ig))
           ! Intramolecular restraint
           CALL section_vals_val_get(collective_section,"INTERMOLECULAR",i_rep_section=ig,&
                l_val=cons_info%colv_intermolecular(ig))
           ! If it is intramolecular let's unset (in case user did it)
           ! the molecule and molname field
           IF (cons_info%colv_intermolecular(ig)) THEN
              CALL section_vals_val_unset(collective_section,"MOLECULE",i_rep_section=ig)
              CALL section_vals_val_unset(collective_section,"MOLNAME",i_rep_section=ig)
           END IF
           ! Let's tag to which molecule we want to apply constraints
           CALL section_vals_val_get(collective_section,"MOLECULE",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(collective_section,"MOLECULE",i_rep_section=ig,&
                   i_val=cons_info%const_colv_mol(ig))
           END IF
           CALL section_vals_val_get(collective_section,"MOLNAME",i_rep_section=ig,&
                n_rep_val=nrep)
           IF (nrep/=0) THEN
              CALL section_vals_val_get(collective_section,"MOLNAME",i_rep_section=ig,&
                   c_val=cons_info%const_colv_molname(ig))
           END IF
           IF  (((cons_info%const_colv_mol(ig)/=0).AND.(cons_info%const_colv_molname(ig)/="UNDEF"))) THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                   "Both MOLNAME and MOLECULE specified for CONSTRAINT section. "//&
CPSourceFileRef,&
                   only_ionode=.TRUE.)
           END IF
           IF  ((cons_info%const_colv_mol(ig)==0).AND.(cons_info%const_colv_molname(ig)=="UNDEF").AND.&
                (.NOT.cons_info%colv_intermolecular(ig))) THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                   "Constraint section error: you have to specify at least one of the "//&
                   "following keywords: MOLECULE, MOLNAME or INTERMOLECULAR! "//&
CPSourceFileRef,&
                   only_ionode=.TRUE.)
           END IF
           NULLIFY(cons_info%colvar_set(ig)%colvar)
           CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig,&
                i_val=icolvar)
           CALL colvar_clone(cons_info%colvar_set(ig)%colvar,&
                colvar_p(icolvar)%colvar)
           CALL section_vals_val_get(collective_section,"TARGET",&
                n_rep_val=n_rep,i_rep_section=ig)
           IF (n_rep /=0) THEN
              CALL section_vals_val_get(collective_section,"TARGET",&
                   r_val=cons_info%const_colv_target(ig),i_rep_section=ig)
           ELSE
              cons_info%const_colv_target(ig) = -HUGE(0.0_dp)
           END IF
           CALL section_vals_val_get(collective_section,"TARGET_GROWTH",&
                r_val=cons_info%const_colv_target_growth(ig),i_rep_section=ig)
        END DO
     END IF
  END IF

END SUBROUTINE read_constraints_section

! *****************************************************************************
!> \brief Reads input and decides if apply restraints instead of constraints
!> \param cons_section ...
!> \param is_restraint ...
!> \param k0 ...
!> \param i_rep_section ...
!> \param label ...
!> \author teo
! *****************************************************************************
SUBROUTINE check_restraint(cons_section, is_restraint, k0, i_rep_section, label)
    TYPE(section_vals_type), POINTER         :: cons_section
    LOGICAL, INTENT(OUT)                     :: is_restraint
    REAL(KIND=dp), INTENT(OUT)               :: k0
    INTEGER, INTENT(IN), OPTIONAL            :: i_rep_section
    CHARACTER(LEN=*), INTENT(IN)             :: label

    CHARACTER(len=*), PARAMETER :: routineN = 'check_restraint', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: nlabel
    INTEGER                                  :: output_unit
    LOGICAL                                  :: explicit, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: restraint_section

  failure      = .FALSE.
  is_restraint = .FALSE.
  logger => cp_get_default_logger()
  output_unit= cp_logger_get_default_io_unit(logger)
  CALL section_vals_get(cons_section,explicit=explicit)
  IF (explicit) THEN
     restraint_section => section_vals_get_subs_vals(cons_section,"RESTRAINT",&
          i_rep_section=i_rep_section)
     CALL section_vals_get(restraint_section,explicit=is_restraint)
     IF (is_restraint) THEN
        CALL section_vals_val_get(restraint_section,"K",r_val=k0)
        IF (output_unit>0) THEN
           nlabel = cp_to_string(i_rep_section)
           WRITE(output_unit,FMT='(T2,"RESTRAINT|",1X,A,F9.6)')&
                "Active restraint on "//label//" section Nr."//&
                TRIM(nlabel)//". K [a.u.]=",k0
        END IF
     END IF
  END IF
END SUBROUTINE check_restraint

END MODULE topology_input


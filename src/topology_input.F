!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****** cp2k/topology_input [1.0] *
!!
!!   NAME
!!     topology_input
!!
!!   FUNCTION
!!     Reads the input sections "topology"
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Added read_topology_section
!!  
!!
!!   SOURCE
!******************************************************************************

MODULE topology_input
  USE colvar_types,                    ONLY: angle_colvar_id,&
                                             colvar_type,&
                                             colvar_p_type,&
                                             dist_colvar_id,&
                                             torsion_colvar_id
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_conn_psf,&
                                             do_conn_mol_psf,&
                                             do_coord_off
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: xstring
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: constraint_info_type,&
                                             topology_parameters_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_input'

  PRIVATE
  PUBLIC :: read_topology_section, read_constraints_section


!!*****
!******************************************************************************

CONTAINS


!!*****
!******************************************************************************
!!****** topology_input/read_topology_section [1.0] *
!!
!!   NAME
!!     read_topology_section
!!
!!   SYNOPSIS
!!     Subroutine read_topology_section(topology, globenv)
!!       Implicit None
!!       Type(topology_parameters_type):: topology
!!       Type(global_environment_type), pointer:: globenv
!!     End Subroutine read_topology_section
!!
!!   FUNCTION
!!     reads the input section topology
!!
!!   AUTHOR
!!     JGH (26-01-2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!******************************************************************************
SUBROUTINE read_topology_section ( topology, topology_section, globenv, error )
    TYPE(topology_parameters_type)           :: topology
    TYPE(section_vals_type), POINTER         :: topology_section
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ia, ie, n_rep
    LOGICAL                                  :: failure

  CALL timeset(routineN,handle)
  CALL xstring(globenv % project_name,ia,ie)
  failure = .FALSE.
  topology % coord_file_name = globenv % project_name(ia:ie) // '.dat'
  topology % conn_file_name = globenv % project_name(ia:ie) // '.top'
  
  IF (.NOT. failure) THEN
     ! Logical
     CALL section_vals_val_get(topology_section,"CHARGE_OCCUP",l_val=topology % charge_occup,error=error)
     CALL section_vals_val_get(topology_section,"CHARGE_BETA",l_val=topology % charge_beta,error=error)
     CALL section_vals_val_get(topology_section,"PARA_RES",l_val=topology % para_res,error=error)
     CALL section_vals_val_get(topology_section,"DUMP_TOP",l_val=topology % dump_topology,error=error)
     CALL section_vals_val_get(topology_section,"REORDER",l_val=topology % reorder_atom,error=error)
     CALL section_vals_val_get(topology_section,"MOL_CHECK",l_val=topology % molecules_check,error=error)
     ! Character
     CALL section_vals_val_get(topology_section,"COORD_FILE",n_rep_val=n_rep, error=error)
     IF (n_rep /=0) THEN
        CALL section_vals_val_get(topology_section,"COORD_FILE",c_val=topology % coord_file_name,error=error)
     END IF
     CALL section_vals_val_get(topology_section,"CONN_FILE",n_rep_val=n_rep, error=error)
     IF (n_rep /=0) THEN
        CALL section_vals_val_get(topology_section,"CONN_FILE",c_val=topology % conn_file_name,error=error)
     END IF
     ! Enumeration
     CALL section_vals_val_get(topology_section,"COORDINATE",i_val=topology % coord_type ,error=error)
     IF (topology % coord_type /= do_coord_off) topology%coordinate = .TRUE.
     CALL section_vals_val_get(topology_section,"CONNECTIVITY",i_val=topology % conn_type,error=error)
     IF (topology % conn_type == do_conn_psf) topology % connectivity = .TRUE.
     IF (topology % conn_type == do_conn_mol_psf) topology % connectivity = .TRUE.
     CALL section_vals_val_get(topology_section,"EXCLUDE_VDW",i_val=topology % exclude_vdw,error=error)
     CALL section_vals_val_get(topology_section,"EXCLUDE_EI",i_val=topology % exclude_ei,error=error)
     CALL section_vals_val_get(topology_section,"EXCLUDE_CHECK",i_val=topology % exclude_check,error=error)
     CALL section_vals_val_get(topology_section,"BONDPARM",i_val=topology % bondparm_type,error=error)
     CALL section_vals_val_get(topology_section,"BONDPARM_FACTOR", &
                               r_val=topology % bondparm_factor,error=error)
     CALL section_vals_val_get(topology_section,"NHCOPT",i_val=topology % nhcopt,error=error)

  END IF
  CALL timestop(handle)
END SUBROUTINE read_topology_section

!!*****
!******************************************************************************
!!****** topology_input/input_constraints [1.0] *
!!
!!   NAME
!!     read_constraints_section
!!
!!   FUNCTION
!!     Read all the distance parameters. Put them in the
!!     constraint_distance array.
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Distance parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!     teo Read the CONSTRAINT section within the new input style
!!
!!   NOTES
!!
!******************************************************************************

SUBROUTINE read_constraints_section(topology, globenv, colvar_p, constraint_section, error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(colvar_p_type), DIMENSION(:), POINTER :: colvar_p
    TYPE(section_vals_type), POINTER         :: constraint_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_constraints_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: icolvar, ig, isize, j, jj, &
                                                n_rep, ncons, stat
    INTEGER, DIMENSION(:), POINTER           :: ilist, tmplist
    LOGICAL                                  :: explicit, failure, freeze_mm, &
                                                freeze_qm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rlist
    TYPE(constraint_info_type), POINTER      :: cons_info
    TYPE(section_vals_type), POINTER         :: collective_section, &
                                                fix_atom_section, &
                                                g3x3_section, g4x6_section,&
                                                hbonds_section

  cons_info    => topology%cons_info
  failure = .FALSE.
  IF (.NOT.failure) THEN
     hbonds_section => section_vals_get_subs_vals(constraint_section,"HBONDS",error=error)
     g3x3_section => section_vals_get_subs_vals(constraint_section,"G3X3",error=error)
     g4x6_section => section_vals_get_subs_vals(constraint_section,"G4X6",error=error)
     fix_atom_section => section_vals_get_subs_vals(constraint_section,"FIXED_ATOMS",error=error)
     collective_section => section_vals_get_subs_vals(constraint_section,"COLLECTIVE",error=error)
     ! HBONDS
     CALL section_vals_get(hbonds_section,explicit=topology%const_hydr, error=error)
     ! G3X3
     CALL section_vals_get(g3x3_section,explicit=explicit, n_repetition=ncons, error=error)
     IF (explicit) THEN
        topology%const_33    =.TRUE.
        cons_info%nconst_g33 =ncons
        !
        ALLOCATE(cons_info%const_g33_mol(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_a(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_b(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_c(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_dab(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_dac(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g33_dbc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ig = 1, ncons
           CALL section_vals_val_get(g3x3_section,"MOLECULE",i_rep_section=ig,&
                i_val=cons_info%const_g33_mol(ig),error=error)
           CALL section_vals_val_get(g3x3_section,"ATOMS",i_rep_section=ig,&
                i_vals=ilist,error=error)
           CALL section_vals_val_get(g3x3_section,"DISTANCES",i_rep_section=ig,&
                r_vals=rlist,error=error)
           cons_info%const_g33_a(ig) = ilist(1)
           cons_info%const_g33_b(ig) = ilist(2)
           cons_info%const_g33_c(ig) = ilist(3)

           cons_info%const_g33_dab(ig) = rlist(1)
           cons_info%const_g33_dac(ig) = rlist(2)
           cons_info%const_g33_dbc(ig) = rlist(3)
        END DO
     END IF
     ! G4X6
     CALL section_vals_get(g4x6_section,explicit=explicit, n_repetition=ncons, error=error)
     IF (explicit) THEN
        topology%const_46   =.TRUE.
        cons_info%nconst_g46 = ncons
        !
        ALLOCATE(cons_info%const_g46_mol(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_a(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_b(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_c(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_d(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dab(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dac(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dbc(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dad(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dbd(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(cons_info%const_g46_dcd(ncons),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ig = 1, ncons
           CALL section_vals_val_get(g4x6_section,"MOLECULE",i_rep_section=ig,&
                i_val=cons_info%const_g46_mol(ig),error=error)
           CALL section_vals_val_get(g4x6_section,"ATOMS",i_rep_section=ig,&
                i_vals=ilist,error=error)
           CALL section_vals_val_get(g4x6_section,"DISTANCES",i_rep_section=ig,&
                r_vals=rlist,error=error)
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
     ! FIXED ATOMS
     CALL section_vals_get(fix_atom_section,explicit=explicit, n_repetition=ncons, error=error)
     IF (explicit) THEN
        NULLIFY(tmplist,tmpstringlist)
        isize = 0
        ALLOCATE(cons_info%fixed_atoms(isize),stat=stat)
        ALLOCATE(cons_info%fixed_molnames(0),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ig = 1, ncons
           CALL section_vals_val_get(fix_atom_section,"LIST", i_rep_section=ig, n_rep_val=n_rep, error=error)
           IF (n_rep /= 0) THEN
              CALL section_vals_val_get(fix_atom_section,"LIST", i_rep_section=ig, &
                   i_vals=tmplist, error=error)
              CALL reallocate(cons_info%fixed_atoms,1,isize+SIZE(tmplist))
              cons_info%fixed_atoms(isize+1:isize+SIZE(tmplist)) = tmplist
              isize = SIZE(cons_info%fixed_atoms)
           END IF
           CALL section_vals_val_get(fix_atom_section,"RANGE", i_rep_section=ig, n_rep_val=n_rep, error=error)
           IF (n_rep /= 0) THEN
              CALL section_vals_val_get(fix_atom_section,"RANGE", i_rep_section=ig, &
                   i_vals=tmplist, error=error)
              CALL reallocate(cons_info%fixed_atoms,1,isize+tmplist(2)-tmplist(1)+1)
              jj = 0
              DO j = tmplist(1), tmplist(2)
                 jj = jj + 1
                 cons_info%fixed_atoms(isize+jj:isize+jj) = j
              END DO
              isize = SIZE(cons_info%fixed_atoms)
           END IF
           CALL section_vals_val_get(fix_atom_section,"MOLNAME", i_rep_section=ig, n_rep_val=n_rep, error=error)
           IF (n_rep /= 0) THEN
              CALL section_vals_val_get(fix_atom_section,"MOLNAME", i_rep_section=ig, &
                   c_vals=tmpstringlist, error=error)
              ALLOCATE(cons_info%fixed_molnames(1:SIZE(tmpstringlist,1)))
              cons_info%fixed_molnames(:)=tmpstringlist(:)
           END IF
           CALL section_vals_val_get(fix_atom_section,"MM_SUBSYS",l_val=freeze_mm,i_rep_section=ig,error=error)
           CALL section_vals_val_get(fix_atom_section,"QM_SUBSYS",l_val=freeze_qm,i_rep_section=ig,error=error)
           cons_info%freeze_mm = cons_info%freeze_mm .OR. freeze_mm
           cons_info%freeze_qm = cons_info%freeze_qm .OR. freeze_qm
        END DO
        IF (isize /= 0 .OR. cons_info%freeze_mm .OR. cons_info%freeze_qm &
             .OR. SIZE(cons_info%fixed_molnames,1)>0 ) THEN
           topology%const_atom    = .TRUE.
           cons_info%nfixed_atoms =  isize
        END IF
     END IF
     ! INTERNAL
     CALL section_vals_get(collective_section,explicit=explicit, n_repetition=ncons, error=error)
     IF (explicit) THEN
        CPPostcondition(ncons<=SIZE(colvar_p),cp_failure_level,routineP,error,failure)
        isize = 0
        ! DISTANCE
        DO ig = 1, ncons
           CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig, i_val=icolvar, error=error)
           CPPostcondition(icolvar<=SIZE(colvar_p),cp_failure_level,routineP,error,failure)
           IF (colvar_p(icolvar)%colvar%type_id==dist_colvar_id) isize = isize + 1
        END DO

        IF (isize /= 0) THEN
           topology%const_dist=.TRUE.
           cons_info%nconst_dist = isize
           !
           ALLOCATE(cons_info%const_dist_mol(isize), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(cons_info%const_dist_a(isize), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(cons_info%const_dist_b(isize), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(cons_info%const_dist_dab(isize), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           isize = 0
           DO ig = 1, ncons
              ! DISTANCE
              CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig, i_val=icolvar, error=error)
              IF (colvar_p(icolvar)%colvar%type_id==dist_colvar_id) THEN
                 isize = isize + 1
                 CALL section_vals_val_get(collective_section,"MOLECULE",i_val=cons_info%const_dist_mol(isize),&
                      i_rep_section=ig,error=error)
                 CALL section_vals_val_get(collective_section,"TARGET",&
                      r_val=cons_info%const_dist_dab(isize),i_rep_section=ig, error=error)
                 cons_info%const_dist_a(isize) = colvar_p(icolvar)%colvar%dist_param%i_at
                 cons_info%const_dist_b(isize) = colvar_p(icolvar)%colvar%dist_param%j_at
              END IF
           END DO
        END IF
        ! ANGLE
        isize = 0
        DO ig = 1, ncons
           CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig, i_val=icolvar, error=error)
           CPPostcondition(icolvar<=SIZE(colvar_p),cp_failure_level,routineP,error,failure)
           IF (colvar_p(icolvar)%colvar%type_id==angle_colvar_id) isize = isize + 1
        END DO
        IF (isize /= 0) THEN
           ! Not implemented yet
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
        ! TORSION
        isize = 0
        DO ig = 1, ncons
           CALL section_vals_val_get(collective_section,"COLVAR",i_rep_section=ig, i_val=icolvar, error=error)
           CPPostcondition(icolvar<=SIZE(colvar_p),cp_failure_level,routineP,error,failure)
           IF (colvar_p(icolvar)%colvar%type_id==torsion_colvar_id) isize = isize + 1
        END DO
        IF (isize /= 0) THEN
           ! Not implemented yet
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     END IF
  END IF

END SUBROUTINE read_constraints_section

END MODULE topology_input

!******************************************************************************

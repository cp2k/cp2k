!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_program
  USE topology_types,                  ONLY: constraint_info_type,&
                                             topology_parameters_type
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE kinds, ONLY: dp
  USE input_constants,                 ONLY: do_conn_psf,&
                                             do_conn_generate,&
                                             do_conn_off,&
                                             do_coord_pdb,&
                                             do_coord_xyz,&
                                             do_coord_off,&
                                             do_skip_vdw_12,&
                                             do_skip_vdw_13,&
                                             do_skip_vdw_14,&
                                             do_skip_ei_12,&
                                             do_skip_ei_13,&
                                             do_skip_ei_14,&
                                             do_bondparm_covalent,&
                                             do_bondparm_vdw,&
                                             do_nhcopt_global,&
                                             do_nhcopt_molecule,&
                                             do_nhcopt_massive,&
                                             do_nhcopt_none
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_input"

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
!!       Type(global_environment_type), Intent (IN):: globenv
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
!!I---------------------------------------------------------------------------I
!!I SECTION: &topology ... &end                                               I
!!I                                                                           I
!!I coord_file       "filename"                                               I 
!!I COORDINATE       [PDB,CHM,...]                                            I
!!I topo_file        "filename"                                               I
!!I TOPOLOGY         [PSF,GENERATE,OFF,...]                                   I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!******************************************************************************
SUBROUTINE read_topology_section ( topology, topology_section, globenv, error )
  IMPLICIT NONE
  ! Arguments
  TYPE(topology_parameters_type)                       :: topology
  TYPE(global_environment_type), INTENT(IN),OPTIONAL   :: globenv
  TYPE(section_vals_type), POINTER                     :: topology_section
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_section', &
      routineP = moduleN//':'//routineN
  INTEGER :: ia, ie, n_rep
  LOGICAL :: failure 

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL xstring(globenv % project_name,ia,ie)
  failure = .FALSE.
  topology % coord_file_name = globenv % project_name(ia:ie) // '.dat'
  topology % conn_file_name = globenv % project_name(ia:ie) // '.top'
  
  IF (.NOT. .FALSE.) THEN
     ! Logical
     CALL section_vals_val_get(topology_section,"CHARGE_OCCUP",l_val=topology % charge_occup,error=error)
     CALL section_vals_val_get(topology_section,"CHARGE_BETA",l_val=topology % charge_beta,error=error)
     CALL section_vals_val_get(topology_section,"PARA_RES",l_val=topology % para_res,error=error)
     CALL section_vals_val_get(topology_section,"DUMP_TOP",l_val=topology % dump_topology,error=error)
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
     CALL section_vals_val_get(topology_section,"EXCLUDE_VDW",i_val=topology % exclude_vdw,error=error)
     CALL section_vals_val_get(topology_section,"EXCLUDE_EI",i_val=topology % exclude_ei,error=error)
     CALL section_vals_val_get(topology_section,"EXCLUDE_CHECK",i_val=topology % exclude_check,error=error)
     CALL section_vals_val_get(topology_section,"BONDPARM",i_val=topology % bondparm_type,error=error)
     CALL section_vals_val_get(topology_section,"NHCOPT",i_val=topology % nhcopt,error=error)

  END IF
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

SUBROUTINE read_constraints_section(topology, globenv, constraint_section, error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv
    TYPE(section_vals_type), POINTER                 :: constraint_section  
      TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'read_constraints_section', &
      routineP = moduleN//':'//routineN
  TYPE(constraint_info_type),POINTER       :: cons_info
  TYPE(section_vals_type), POINTER         ::  g3x3_section, g4x6_section, fix_atom_section,&
                                               internal_section, dist_section, angle_section,&
                                               tors_section, comb_section
  LOGICAL :: explicit, failure, freeze_mm, freeze_qm
  INTEGER :: ncons, stat, isize, j, k, jj, icons, ig, n_rep
  INTEGER, DIMENSION(:), POINTER       :: ilist, tmplist
  REAL(KIND=dp), DIMENSION(:), POINTER :: rlist 

  cons_info    => topology%cons_info
  failure = .FALSE.
  IF (.NOT.failure) THEN
     g3x3_section => section_vals_get_subs_vals(constraint_section,"G3X3",error=error)
     g4x6_section => section_vals_get_subs_vals(constraint_section,"G4X6",error=error)
     fix_atom_section => section_vals_get_subs_vals(constraint_section,"FIXED_ATOMS",error=error)
     internal_section => section_vals_get_subs_vals(constraint_section,"INTERNALS",error=error)
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
        NULLIFY(tmplist)
        isize = 0
        ALLOCATE(cons_info%fixed_atoms(isize),stat=stat)
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
           CALL section_vals_val_get(fix_atom_section,"MM_SUBSYS",l_val=freeze_mm,i_rep_section=ig,error=error)
           CALL section_vals_val_get(fix_atom_section,"QM_SUBSYS",l_val=freeze_qm,i_rep_section=ig,error=error)
           cons_info%freeze_mm = cons_info%freeze_mm .OR. freeze_mm
           cons_info%freeze_qm = cons_info%freeze_qm .OR. freeze_qm
        END DO
        IF (isize /= 0 .OR. cons_info%freeze_mm .OR. cons_info%freeze_qm ) THEN
           topology%const_atom    = .TRUE.
           cons_info%nfixed_atoms =  isize
        END IF
     END IF
     ! INTERNAL
     CALL section_vals_get(internal_section,explicit=explicit, n_repetition=ncons, error=error)
     IF (explicit) THEN
        isize = 0
        ! DISTANCE
        DO ig = 1, ncons
           dist_section => section_vals_get_subs_vals(internal_section,"DISTANCE",i_rep_section=ig,error=error)
           CALL section_vals_get(dist_section,explicit=explicit, n_repetition=icons, error=error)
           isize = isize + icons 
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
              dist_section => section_vals_get_subs_vals(internal_section,"DISTANCE",i_rep_section=ig,error=error)
              CALL section_vals_get(dist_section,explicit=explicit, n_repetition=icons, error=error)
              DO j = 1, icons
                 isize = isize + 1
                 CALL section_vals_val_get(dist_section,"MOLECULE", i_val=cons_info%const_dist_mol(isize),&
                      i_rep_section=j,error=error)
                 CALL section_vals_val_get(dist_section,"ATOMS", i_vals=tmplist,&
                      i_rep_section=j,error=error)
                 CALL section_vals_val_get(dist_section,"DISTANCE",&
                      r_val=cons_info%const_dist_dab(isize),i_rep_section=j, error=error)
                 cons_info%const_dist_a(isize) = tmplist(1)
                 cons_info%const_dist_b(isize) = tmplist(2)
              END DO
           END DO
        END IF
        ! ANGLE
        isize = 0
        DO ig = 1, ncons
           angle_section => section_vals_get_subs_vals(internal_section,"ANGLE",i_rep_section=ig,error=error)
           CALL section_vals_get(angle_section,explicit=explicit, n_repetition=icons, error=error)
           isize = isize + icons 
        END DO
        IF (isize /= 0) THEN
           ! Not implemented yet
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
        ! TORSION
        isize = 0
        DO ig = 1, ncons
           tors_section => section_vals_get_subs_vals(internal_section,"TORSION",i_rep_section=ig,error=error)
           CALL section_vals_get(tors_section,explicit=explicit, n_repetition=icons, error=error)
           isize = isize + icons 
        END DO
        IF (isize /= 0) THEN
           ! Not implemented yet
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
        ! COMB
        isize = 0
        DO ig = 1, ncons
           comb_section => section_vals_get_subs_vals(internal_section,"COMB",i_rep_section=ig,error=error)
           CALL section_vals_get(comb_section,explicit=explicit, n_repetition=icons, error=error)
           isize = isize + icons 
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

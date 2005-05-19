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
  USE cp_error_handling,               ONLY: cp_error_type
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
  USE input_cp2k, ONLY: parsed_cp2k_input
  USE kinds, ONLY: dp
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
!!I TOPOLOGY         [PSF,GROMOS,GENERATE,OFF,...]                            I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!******************************************************************************
SUBROUTINE read_topology_section ( topology, globenv )


  TYPE(topology_parameters_type)                       :: topology
  TYPE(global_environment_type), INTENT(IN),OPTIONAL   :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_section', &
      routineP = moduleN//':'//routineN

  CHARACTER(LEN=20)                                    :: str2,string
  CHARACTER(LEN=6)                                     :: label
  INTEGER                                              :: i,ia,ie,ierror,ilen,iw

  LOGICAL                                              :: PM1,PM2,PM3,PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

 
  CALL write_checkpoint_information("entering "//routineN,globenv)

  IF(PM1) WRITE(iw,*) "  Entering read_topology_section"

  CALL xstring(globenv % project_name,ia,ie)

  topology % coord_file_name = globenv % project_name(ia:ie) // '.dat'
  topology % conn_file_name = globenv % project_name(ia:ie) // '.top'

  ! parse the input section
  label = '&TOPOL'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print%level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &TOPOL found on file '
        WRITE ( iw, '( T2, A )' ) globenv % input_file_name
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 40
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE ('CHARGE_O','CHARGE_OCCUP')
           topology % charge_occup = .TRUE.
        CASE ('CHARGE_B','CHARGE_BETA')
           topology % charge_beta = .TRUE.
        CASE ('PARA_RES')
           topology % para_res = .TRUE.
        CASE ('DUMP_TOP')
           topology % dump_topology = .TRUE.
        CASE ('REORDER')
           topology % reorder_atom = .TRUE.
        CASE ('COORD_FI','COORD_FILE','COORD_FILE_NAME')
           ilen = 40
           CALL get_next(topology % coord_file_name,ilen)
        CASE ('COORDINA','COORDINATE')
           topology % coordinate = .TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "PDB") THEN
              topology%coordinate = .TRUE.
              topology%coord_type = "PDB"
           ELSE IF (str2 == "CHM") THEN
              topology%coordinate = .TRUE.
              topology%coord_type = "CHM"
           ELSE IF (str2 == "OFF") THEN
              topology%coordinate = .FALSE.
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option coordinate type' )
           END IF
        CASE ('CONN_FIL','CONN_FILE')
           topology % connectivity = .TRUE.
           ilen = 40
           CALL get_next(topology % conn_file_name,ilen)
        CASE ('CONNECTI','CONNECTIVITY')
           topology % connectivity = .TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "PSF") THEN
              topology % connectivity = .TRUE.
              topology % conn_type = "PSF"
           ELSE IF (str2 == "GROMOS") THEN
              topology % connectivity = .TRUE.
              topology % conn_type = "GROMOS"
           ELSE IF (str2 == "GENERATE") THEN
              topology % connectivity = .FALSE.
              topology % conn_type = "GENERATE"
           ELSE IF (str2 == "OFF") THEN
              topology % connectivity = .FALSE.
              topology % conn_type = "OFF"
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option topology type' )
           END IF
        CASE ('EXCLUDE_VDW')
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           topology % exclude_vdw = str2
        CASE ('EXCLUDE_EI')
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           topology % exclude_ei = str2
        CASE ('EXCLUDE_CH','EXCLUDE_CHECK')
           CALL get_next(topology%exclude_check)
           WRITE(*,*) "exclude_check ",topology%exclude_check
        CASE ('BONDPARM')
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "COVALENT") THEN
              topology % bondparm_type = "COVALENT"
           ELSE IF (str2 == "VDW") THEN
              topology % bondparm_type = "VDW"
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option bondparm type' )
           END IF
           CALL get_next(topology%bondparm_factor)
        CASE ('NHCOPT')
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "GLOBAL") THEN
              topology % nhcopt = "GLOBAL"
           ELSE IF (str2 == "MOLECULE") THEN
              topology % nhcopt = "MOLECULE"
           ELSE IF (str2 == "MASSIVE") THEN
              topology % nhcopt = "MASSIVE"
           ELSE IF (str2 == "NONE") THEN
              topology % nhcopt = "NONE"
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option nhc opt' )
           END IF
        CASE ('CONSTRAINT')
              topology % constraint = .TRUE.
        CASE DEFAULT
           IF(PM1)  WRITE(iw,*) "read_topology_section unknown option"//string
        END SELECT
       CALL read_line
     END DO
  END IF

  CALL parser_end


  IF(PM1)  WRITE(iw,*) "  Exiting  read_topology_section"
  CALL write_checkpoint_information("exiting "//routineN,globenv)

END SUBROUTINE read_topology_section



!!*****
!******************************************************************************
!!****** force_fields/input_constraints [1.0] *
!!
!!   NAME
!!     read_constraints_section
!!
!!   FUNCTION
!!     Read all the distance parameters. Put them in the
!!     constraint_distance array.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Distance parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I constraints                                                               I
!!I   HYDROGEN                                                                I
!!I   ATOM                                                                    I
!!I   DIST MOLKIND ATM_A ATM_B D_AB                                           I
!!I   G3x3 MOLKIND ATM_A ATM_B ATM_C D_AB D_AC D_BC                           I
!!I   G4x6 MOLKIND ATM_A ATM_B ATM_C ATM_D D_AB D_AC D_BC D_AD D_BD D_CD      I
!!I                                                                           I
!!I end constraints                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
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
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string
    INTEGER                                  :: i, ierror, ilen, iw

  LOGICAL                                    :: PM1, PM2, PM3, PM4

  IF (parsed_cp2k_input(constraint_section,check_this_section=.TRUE.,error=error)) THEN
     CALL read_constraints_section_new(topology, constraint_section=constraint_section,&
          error=error)
     RETURN
  END IF
  
!TODEL_START_TL
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routineN,globenv)

  IF(PM1) WRITE(iw,*) "  Entering read_constraints_section"

  cons_info => topology%cons_info

  ! parse the input section
  label = '&CONSTRAINTS'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print%level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &CONSTRAINTS found on file '
        WRITE ( iw, '( T2, A )' ) globenv % input_file_name
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 7
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        IF ( INDEX ( string, 'ATOM' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_atom=.TRUE.
           CALL stop_program ( "read_constraints_section", &
                               "CONATOM not done yet " )
        ELSE IF ( INDEX ( string, 'CONHYDR' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_hydr=.TRUE.
           CALL stop_program ( "read_constraints_section", &
               "CONHYDR will not work until force field is completed " )
        ELSE IF ( INDEX ( string, 'DIST' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_dist=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=cons_info%nconst_dist + 1
           CALL reallocate(cons_info%const_dist_mol,1,i)
           CALL reallocate(cons_info%const_dist_a,1,i)
           CALL reallocate(cons_info%const_dist_b,1,i)
           CALL reallocate(cons_info%const_dist_dab,1,i)
           CALL get_next(cons_info%const_dist_mol(i))
           CALL get_next(cons_info%const_dist_a(i))
           CALL get_next(cons_info%const_dist_b(i))
           CALL get_next(cons_info%const_dist_dab(i))
           cons_info%nconst_dist = i
        ELSE IF ( INDEX ( string, 'G3X3' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_33=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=cons_info%nconst_g33 + 1
           CALL reallocate(cons_info%const_g33_mol,1,i)
           CALL reallocate(cons_info%const_g33_a,1,i)
           CALL reallocate(cons_info%const_g33_b,1,i)
           CALL reallocate(cons_info%const_g33_c,1,i)
           CALL reallocate(cons_info%const_g33_dab,1,i)
           CALL reallocate(cons_info%const_g33_dac,1,i)
           CALL reallocate(cons_info%const_g33_dbc,1,i)
           CALL get_next(cons_info%const_g33_mol(i))
           CALL get_next(cons_info%const_g33_a(i))
           CALL get_next(cons_info%const_g33_b(i))
           CALL get_next(cons_info%const_g33_c(i))
           CALL get_next(cons_info%const_g33_dab(i))
           CALL get_next(cons_info%const_g33_dac(i))
           CALL get_next(cons_info%const_g33_dbc(i))
           cons_info%nconst_g33 = i
        ELSE IF ( INDEX ( string, 'G4X6' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_46=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=cons_info%nconst_g46 + 1
           CALL reallocate(cons_info%const_g46_mol,1,i)
           CALL reallocate(cons_info%const_g46_a,1,i)
           CALL reallocate(cons_info%const_g46_b,1,i)
           CALL reallocate(cons_info%const_g46_c,1,i)
           CALL reallocate(cons_info%const_g46_d,1,i)
           CALL reallocate(cons_info%const_g46_dab,1,i)
           CALL reallocate(cons_info%const_g46_dac,1,i)
           CALL reallocate(cons_info%const_g46_dbc,1,i)
           CALL reallocate(cons_info%const_g46_dad,1,i)
           CALL reallocate(cons_info%const_g46_dbd,1,i)
           CALL reallocate(cons_info%const_g46_dcd,1,i)
           CALL get_next(cons_info%const_g46_mol(i))
           CALL get_next(cons_info%const_g46_a(i))
           CALL get_next(cons_info%const_g46_b(i))
           CALL get_next(cons_info%const_g46_c(i))
           CALL get_next(cons_info%const_g46_d(i))
           CALL get_next(cons_info%const_g46_dab(i))
           CALL get_next(cons_info%const_g46_dac(i))
           CALL get_next(cons_info%const_g46_dad(i))
           CALL get_next(cons_info%const_g46_dbc(i))
           CALL get_next(cons_info%const_g46_dbd(i))
           CALL get_next(cons_info%const_g46_dcd(i))
           cons_info%nconst_g46 = i
         END IF
       CALL read_line
     END DO
  END IF

  CALL parser_end

  IF(PM1)  WRITE(iw,*) "  Exiting  read_constraints_section"
  CALL write_checkpoint_information("exiting "//routineN,globenv)
!TODEL_END_TL

END SUBROUTINE read_constraints_section

!!*****
!******************************************************************************
!!****** force_fields/read_constraints_section_new [1.0] *
!!
!!   NAME
!!     read_constraints_section_new
!!
!!   FUNCTION
!!     Read the CONSTRAINT section within the new input style
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!******************************************************************************

SUBROUTINE read_constraints_section_new(topology, constraint_section, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(topology_parameters_type), INTENT(INOUT)      :: topology
  TYPE(section_vals_type), POINTER                   :: constraint_section  
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  ! Local Variables
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
        topology%constraint  =.TRUE.
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
        topology%constraint =.TRUE.
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
           DO k = 1, n_rep
              CALL section_vals_val_get(fix_atom_section,"LIST", i_rep_section=ig, i_rep_val=k,&
                   i_vals=tmplist, error=error)
              CALL reallocate(cons_info%fixed_atoms,1,isize+SIZE(tmplist))
              cons_info%fixed_atoms(isize+1:isize+SIZE(tmplist)) = tmplist
              isize = SIZE(cons_info%fixed_atoms)
           END DO
           CALL section_vals_val_get(fix_atom_section,"RANGE", i_rep_section=ig, n_rep_val=n_rep, error=error)
           DO k = 1, n_rep
              CALL section_vals_val_get(fix_atom_section,"RANGE", i_rep_section=ig, i_rep_val=k,&
                   i_vals=tmplist, error=error)
              CALL reallocate(cons_info%fixed_atoms,1,isize+tmplist(2)-tmplist(1)+1)
              jj = 0
              DO j = tmplist(1), tmplist(2)
                 jj = jj + 1
                 cons_info%fixed_atoms(isize+jj:isize+jj) = j
              END DO
              isize = SIZE(cons_info%fixed_atoms)
           END DO
           CALL section_vals_val_get(fix_atom_section,"MM_SUBSYS",l_val=freeze_mm,i_rep_section=ig,error=error)
           CALL section_vals_val_get(fix_atom_section,"QM_SUBSYS",l_val=freeze_qm,i_rep_section=ig,error=error)
           cons_info%freeze_mm = cons_info%freeze_mm .OR. freeze_mm
           cons_info%freeze_qm = cons_info%freeze_qm .OR. freeze_qm
        END DO
        IF (isize /= 0 .OR. cons_info%freeze_mm .OR. cons_info%freeze_qm ) THEN
           topology%constraint    = .TRUE.
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
              IF (explicit) THEN
                 DO j = 1, icons
                    isize = isize + 1
                    CALL section_vals_val_get(fix_atom_section,"MOLECULE", i_val=cons_info%const_dist_mol(isize),&
                         error=error)
                    CALL section_vals_val_get(fix_atom_section,"ATOMS", i_vals=tmplist, error=error)
                    CALL section_vals_val_get(fix_atom_section,"DISTANCE",&
                         r_val=cons_info%const_dist_dab(isize), error=error)
                    cons_info%const_dist_a(isize) = tmplist(1)
                    cons_info%const_dist_b(isize) = tmplist(2)
                 END DO
              END IF
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

END SUBROUTINE read_constraints_section_new

END MODULE topology_input

!******************************************************************************

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
           ELSE IF (str2 == "OHH_MC") THEN
              topology % connectivity = .FALSE.
              topology % conn_type = "OHH_MC"
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
!!I   ATOM ??????                                                             I
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

SUBROUTINE read_constraints_section(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'read_constraints_section', &
      routineP = moduleN//':'//routineN

    TYPE(constraint_info_type),POINTER       :: cons_info
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string
    INTEGER                                  :: i, ierror, ilen, iw

  LOGICAL                                  :: PM1, PM2, PM3, PM4

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
           CALL get_next(cons_info%const_g46_dbc(i))
           CALL get_next(cons_info%const_g46_dad(i))
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

END SUBROUTINE read_constraints_section

!!*****
!******************************************************************************

END MODULE topology_input

!******************************************************************************

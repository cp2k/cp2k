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
  USE topology_types,                  ONLY: topology_parameters_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_input"

  PRIVATE
  PUBLIC :: read_topology_section

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


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: str2, string
    CHARACTER(LEN=6)                         :: label
    INTEGER                                  :: i,ia, ie, ierror, ilen, iw

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routineN,globenv)

  iw = globenv % scr
  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Entering read_topology_section"
  END IF

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
        ilen = 8
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE ('DUMP_TOP')
           topology % dump_topology = .TRUE.
        CASE ('REORDER')
           topology % reorder_atom = .TRUE.
        CASE ('COORD_FI')
           ilen = 40
           CALL get_next(topology % coord_file_name,ilen)
        CASE ('COORDINA')
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
        CASE ('CONN_FIL')
           topology % connectivity = .TRUE.
           ilen = 40
           CALL get_next(topology % conn_file_name,ilen)
        CASE ('CONNECTI')
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
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option topology type' )
           END IF
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
        CASE ('CONATOM')
           topology%constraint=.TRUE.
           topology%const_atom=.TRUE.
           CALL stop_program ( "read_topology_section", &
               "CONATOM not done yet " )
        CASE ('CONHYDR')
           topology%constraint=.TRUE.
           topology%const_hydr=.TRUE.
           CALL stop_program ( "read_topology_section", &
               "CONHYDR not implemented yet" )
        CASE ('CONDIST')
           topology%constraint=.TRUE.
           topology%const_dist=.TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "ADD") THEN
             i=topology%nconst_dist + 1
             CALL reallocate(topology%const_dist_mol,1,i)
             CALL reallocate(topology%const_dist_a,1,i)
             CALL reallocate(topology%const_dist_b,1,i)
             CALL reallocate(topology%const_dist_dab,1,i)
             CALL get_next(topology%const_dist_mol(i))
             CALL get_next(topology%const_dist_a(i))
             CALL get_next(topology%const_dist_b(i))
             CALL get_next(topology%const_dist_dab(i))
             topology%nconst_dist = i
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option in conatom' )
           ENDIF
        CASE ('CONG3X3')
           topology%constraint=.TRUE.
           topology%const_33=.TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "ADD") THEN
             i=topology%nconst_g33 + 1
             CALL reallocate(topology%const_g33_mol,1,i)
             CALL reallocate(topology%const_g33_a,1,i)
             CALL reallocate(topology%const_g33_b,1,i)
             CALL reallocate(topology%const_g33_c,1,i)
             CALL reallocate(topology%const_g33_dab,1,i)
             CALL reallocate(topology%const_g33_dac,1,i)
             CALL reallocate(topology%const_g33_dbc,1,i)
             CALL get_next(topology%const_g33_mol(i))
             CALL get_next(topology%const_g33_a(i))
             CALL get_next(topology%const_g33_b(i))
             CALL get_next(topology%const_g33_c(i))
             CALL get_next(topology%const_g33_dab(i))
             CALL get_next(topology%const_g33_dac(i))
             CALL get_next(topology%const_g33_dbc(i))
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_topology_section', &
                                'unknown option in conatom' )
           ENDIF
        CASE DEFAULT
           CALL stop_program ( "read_topology_section", &
               "Unknown option "//string )
        END SELECT
       CALL read_line
     END DO
  END IF

  CALL parser_end


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting  read_topology_section"
  END IF

  CALL write_checkpoint_information("exiting "//routineN,globenv)

END SUBROUTINE read_topology_section

!!*****
!******************************************************************************

END MODULE topology_input

!******************************************************************************

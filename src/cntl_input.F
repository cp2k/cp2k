!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/cntl_input [1.0] *
!!
!!   NAME
!!     cntl_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE cntl_input
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
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
  PRIVATE
  PUBLIC :: read_cntl_section

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** cntl_input/read_cntl_section [1.0] *
!!
!!   NAME
!!     read_cntl_section
!!
!!   SYNOPSIS
!!     Subroutine read_cntl_section(setup, globenv)
!!       Implicit None
!!       Type(setup_parameters_type), Intent (INOUT):: setup
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine read_cntl_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   INPUTS
!!I---------------------------------------------------------------------------I
!!I SECTION: &cntl ... &end                                                   I
!!I                                                                           I
!!I printlevel       globenv%print_level                                      I
!!I set_file         "filename"                                               I
!!I coord_file       "filename"                                               I 
!!I COORDINATE       [PDB,CHM,...]                                            I
!!I topo_file        "filename"                                               I
!!I TOPOLOGY         [PSF,GROMOS,GENERATE,OFF,...]                            I
!!I symmetry         [on,off]                                                 I
!!I path_integrals   [on,off]                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************

SUBROUTINE read_cntl_section ( setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source, group
  CHARACTER ( LEN = 20 ) :: string, str2
  CHARACTER ( LEN = 5 ) :: label
  LOGICAL :: found, setf, datf, topf

!------------------------------------------------------------------------------


!..defaults
  setup % unit_type = 'KELVIN'
  setup % symmetry = .FALSE.
  setup % path_integrals = .FALSE.
  CALL xstring(globenv % project_name,ia,ie)
  setup % set_file_name = globenv % project_name(ia:ie) // '.set'
  setup % coord_file_name = globenv % project_name(ia:ie) // '.dat'
  setup % topology_file_name = globenv % project_name(ia:ie) // '.top'

  setup % coordinate = .FALSE.
  setup % topology = .FALSE.
  setup % topology_type = 'OFF'

  iw = globenv % scr

!..parse the input section
  label = '&CNTL'
  CALL parser_init(globenv % input_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     IF (globenv % ionode) &
        WRITE ( iw, '( a )' ) ' No input section &CNTL found '
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser( 'read_cntl_section','unknown option: '//TRIM(string))
        CASE ( 'PRINTLEV')
           CALL get_next ( globenv % print_level )
        CASE ('SYMMETRY')
           ilen = 3
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF ( str2(1:2) == "ON" ) setup%symmetry=.TRUE.
           IF ( str2(1:3) == "OFF" ) setup%symmetry=.FALSE.
        CASE ('PATH_INT')
           ilen = 3
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF ( str2(1:2) == "ON" ) setup%path_integrals=.TRUE.
           IF ( str2(1:3) == "OFF" ) setup%path_integrals=.FALSE.
        CASE ( 'SET_FILE')
           ilen = 40
           CALL get_next(setup % set_file_name,ilen)
        CASE ( 'INPUT_FI')
           ilen = 40
           CALL get_next(setup % coord_file_name,ilen)
        CASE ( 'COORDINA')
           setup % coordinate = .TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "PDB") THEN
              setup%coordinate = .TRUE.
              setup%coordinate_type = "PDB"
           ELSE IF (str2 == "CHM") THEN
              setup%coordinate = .TRUE.
              setup%coordinate_type = "CHM"
           ELSE IF (str2 == "OFF") THEN
              setup%coordinate = .FALSE.
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_cntl_section', &
                                'unknown option coordinate type' )
           END IF
        CASE ( 'TOPO_FIL')
           setup % topology = .TRUE.
           ilen = 40
           CALL get_next(setup % topology_file_name,ilen)
        CASE ( 'TOPOLOGY')
           setup % topology = .TRUE.
           ilen = 40
           CALL get_next(str2,ilen)
           CALL uppercase(str2)
           IF (str2 == "PSF") THEN
              setup % topology = .TRUE.
              setup % topology_type = "PSF"
           ELSE IF (str2 == "GROMOS") THEN
              setup % topology = .TRUE.
              setup % topology_type = "GROMOS"
           ELSE IF (str2 == "GENERATE") THEN
              setup % topology = .TRUE.
              setup % topology_type = "GENERATE"
           ELSE IF (str2 == "OFF") THEN
              setup % topology = .FALSE.
           ELSE
              CALL p_error()
              CALL stop_parser( 'read_cntl_section', &
                                'unknown option topology type' )
           END IF
        END SELECT

! check for trailing rubbish
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

!..check for files
  setf =.TRUE.
  INQUIRE ( FILE = setup % set_file_name, EXIST = found )
  IF ( .NOT. found ) THEN
     setup % set_file_name = globenv % input_file_name
     setf =.FALSE.
  END IF

  datf =.TRUE.
!..can not check for coord_files in PIMD case, since there are several
!..of those and how many there should be is determined later
  IF (.NOT. setup%path_integrals ) THEN
     INQUIRE ( FILE = setup % coord_file_name, EXIST = found )
     IF ( .NOT. found ) THEN
        setup % coord_file_name = globenv % input_file_name
        datf =.FALSE.
     END IF
  END IF

  topf =.TRUE.
  INQUIRE ( FILE = setup % coord_file_name, EXIST = found )
  IF ( .NOT. found ) topf =.FALSE.

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 ) THEN
        WRITE ( iw, '( A, T61, A )' ) &
             ' CONTROL| Unit type ', ADJUSTR ( setup % unit_type )
        IF ( setup%symmetry ) THEN
           WRITE (iw,'(A,T78,A)') ' CONTROL| Use of symmetry ','Yes'
        ELSE
           WRITE (iw,'(A,T78,A)') ' CONTROL| Use of symmetry ',' No'
        END IF
        IF ( setup%path_integrals ) THEN
           WRITE (iw,'(A,T78,A)') ' CONTROL| Path Integrals ','Yes'
        ELSE
           WRITE (iw,'(A,T78,A)') ' CONTROL| Path Integrals ',' No'
        END IF
        IF ( setf ) &
           WRITE ( iw, '( A, T41, A )' ) ' CONTROL| Set file name', &
                ADJUSTR ( setup % set_file_name )
        IF ( datf ) &
           WRITE ( iw, '( A, T41, A )' ) ' CONTROL| Data file name', &
                ADJUSTR ( setup % coord_file_name )
        IF ( topf ) &
           WRITE ( iw, '( A, T41, A )' ) ' CONTROL| Topology file name', &
                ADJUSTR ( setup % topology_file_name )
        WRITE ( iw, '( A, T76, I5 )' ) &
             ' CONTROL| Print level ', globenv % print_level
        WRITE ( iw, '( )' )
     END IF
  END IF

END SUBROUTINE read_cntl_section

!!*****
!******************************************************************************

END MODULE cntl_input

!******************************************************************************

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

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE global_types, ONLY : global_environment_type
  USE input_types, ONLY : setup_parameters_type
  USE kinds, ONLY : dbl
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, stop_parser
  USE string_utilities, ONLY : uppercase, xstring

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
!!I simulation       [md,pimd,debug]                                          I
!!I printlevel       globenv%print_level                                      I
!!I units            [kelvin,atomic]                                          I
!!I periodic         [0,1][0,1][0,1]                                          I
!!I Ewald_type       [ewald,pme,spme]                                         I
!!I Ewald_param      alpha[,gmax,ns_max,epsilon]                              I
!!I set_file         "filename"                                               I
!!I input_file       "filename"                                               I
!!I symmetry         [on,off]                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************

SUBROUTINE read_cntl_section ( setup, ewald_param, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source, group
  CHARACTER ( LEN = 20 ) :: string, str2
  CHARACTER ( LEN = 5 ) :: label
  CHARACTER ( LEN = 3 ), PARAMETER :: yn ( 0:1 ) = (/ ' NO', 'YES' /)

!------------------------------------------------------------------------------

!..defaults
  setup % run_type = 'MD'
  setup % unit_type = 'KELVIN'
  setup % perd = 1
  setup % symmetry = .false.
  ewald_param % alpha = 0.4_dbl
  ewald_param % gmax = 10
  ewald_param % o_spline = 5
  ewald_param % ns_max = 10
  ewald_param % epsilon = 1.e-6_dbl
  ewald_param % ewald_type = 'NONE'
  CALL xstring(globenv % project_name,ia,ie)
  setup % set_file_name = globenv % project_name(ia:ie) // '.set'
  setup % input_file_name = globenv % project_name(ia:ie) // '.dat'

  iw = globenv % scr

!..parse the input section
  label = '&CNTL'
  CALL parser_init(globenv % input_file_name,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF (globenv % ionode) &
        WRITE ( iw, '( a )' ) ' No input section &CNTL found '
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser( 'read_cntl_section','unknown option')
        CASE ( 'SIMULATI')
           ilen = 20
           CALL cfield(setup % run_type,ilen)
           CALL uppercase(setup % run_type )
        CASE ( 'PRINTLEV')
           globenv % print_level = get_int()
        CASE ( 'UNITS')
           ilen = 20
           CALL cfield(setup % unit_type,ilen)
           CALL uppercase(setup % unit_type )
        CASE ( 'PERIODIC')
           setup % perd(1) = get_int()
           setup % perd(2) = get_int()
           setup % perd(3) = get_int()
        CASE ( 'EWALD_TY')
           ilen=20
           CALL cfield(string,ILEN)
           CALL uppercase ( string )
           SELECT CASE(string)
           CASE DEFAULT
              CALL p_error()
              CALL stop_parser( 'read_cntl_section', &
                   'unknown option for Ewald type' )
           CASE( 'EWALD')
              ewald_param % ewald_type = 'ewald'
              CALL uppercase(ewald_param % ewald_type )
           CASE( 'PME')
              ewald_param % ewald_type = 'pme'
              CALL uppercase(ewald_param % ewald_type )
           CASE( 'SPME')
              ewald_param % ewald_type = 'spme'
              CALL uppercase(ewald_param % ewald_type )
           END SELECT

! if no type specified, assume ewald
        CASE ( 'EWALD_PA')
           ewald_param % alpha = get_real()
           SELECT CASE (ewald_param % ewald_type )
           CASE DEFAULT
              ewald_param % gmax = get_int()
           CASE ( 'PME')
              ewald_param % ns_max = get_int()
              IF ( test_next() == 'N' ) THEN
                 ewald_param % epsilon = get_real()
              END IF
           CASE ( 'SPME')
              ewald_param % gmax = get_int()
              ewald_param % o_spline = get_int()
           END SELECT
        CASE ('SYMMETRY')
           ilen = 3
           CALL cfield(str2,ilen)
           CALL uppercase(str2)
           IF ( str2(1:2) == "ON" ) setup%symmetry=.true.
           IF ( str2(1:3) == "OFF" ) setup%symmetry=.false.
        CASE ( 'SET_FILE')
           ilen = 40
           CALL cfield(setup % set_file_name,ilen)
        CASE ( 'INPUT_FI')
           ilen = 40
           CALL cfield(setup % input_file_name,ilen)
        END SELECT

! check for trailing rubbish
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 ) THEN
        WRITE ( iw, '( A, T61, A )' ) &
             ' CONTROL| Run type ', ADJUSTR ( setup % run_type )
        WRITE ( iw, '( A, T61, A )' ) &
             ' CONTROL| Unit type ', ADJUSTR ( setup % unit_type )
        WRITE ( iw, '( A, T78, A )' ) ' CONTROL| Periodic in X direction ', &
             yn(setup % perd(1))
        WRITE ( iw, '( A, T78, A )' ) ' CONTROL| Periodic in Y direction ', &
             yn(setup % perd(2))
        WRITE ( iw, '( A, T78, A )' ) ' CONTROL| Periodic in Z direction ', &
             yn(setup % perd(3))
        IF ( setup%symmetry ) THEN
           WRITE (iw,'(A,T78,A)') ' CONTROL| Use of symmetry ','Yes'
        ELSE
           WRITE (iw,'(A,T78,A)') ' CONTROL| Use of symmetry ',' No'
        END IF
        WRITE ( iw, '( A, T41, A )' ) ' CONTROL| Set file name', &
             ADJUSTR ( setup % set_file_name )
        WRITE ( iw, '( A, T41, A )' ) ' CONTROL| Input file name', &
             ADJUSTR ( setup % input_file_name )
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

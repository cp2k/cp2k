!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE fist_input
  
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE global_types, ONLY : global_environment_type
  USE input_types, ONLY : setup_parameters_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, stop_parser
  USE string_utilities, ONLY : uppercase, xstring
  
  PRIVATE
  PUBLIC :: read_fist_section
  
CONTAINS

!!>---------------------------------------------------------------------------!
!!  SECTION: &fist ... &end                                                   !
!!                                                                            !
!!  simulation       [md,debug]                                               !
!!  printlevel       fistpar%print_level                                      !
!!  units            [kelvin,atomic]                                          !
!!  periodic         [0,1][0,1][0,1]                                          !
!!  Ewald_type       [pme_gauss,ewald_gauss]                                  !
!!  Ewald_param      alpha,[gmax,ns_max,epsilon]                              !
!!  set_file         "filename"                                               !
!!  input_file       "filename"                                               !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_fist_section ( setup, ewald_param, fistpar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: fistpar
  
! Locals
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source, group
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 5 ) :: label
  CHARACTER ( LEN = 3 ), PARAMETER :: yn ( 0:1 ) = (/ ' NO', 'YES' /)
  
!------------------------------------------------------------------------------
  
!..defaults
  setup % run_type = 'MD'
  setup % unit_type = 'KELVIN'
  setup % perd = 1
  ewald_param % alpha = 0.4_dbl
  ewald_param % gmax = 10
  ewald_param % ns_max = 10
  ewald_param % epsilon = 1.e-6_dbl
  ewald_param % ewald_type = 'NONE'
  CALL xstring(fistpar % project_name,ia,ie)
  setup % set_file_name = fistpar % project_name(ia:ie) // '.set'
  setup % input_file_name = fistpar % project_name(ia:ie) // '.dat'
  
  iw = fistpar % scr
     
!..parse the input section
  label = '&FIST'
  CALL parser_init(fistpar % input_file_name,label,ierror,fistpar)
  IF (ierror /= 0 ) THEN
     IF (fistpar % ionode) &
        WRITE ( iw, '( a )' ) ' No input section &FIST found '
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser( 'read_fist','unknown option')
        CASE ( 'SIMULATI')
           ilen = 20
           CALL cfield(setup % run_type,ilen)
           CALL uppercase(setup % run_type )
        CASE ( 'PRINTLEV')
           fistpar % print_level = get_int()
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
           CASE( 'EWALD_GAUSS')
              ewald_param % ewald_type = 'ewald_gauss'
              CALL uppercase(ewald_param % ewald_type )
           CASE( 'PME_GAUSS')
              ewald_param % ewald_type = 'pme_gauss'
              CALL uppercase(ewald_param % ewald_type )
           END SELECT
              
! if no type specified, assume ewald_gauss
        CASE ( 'EWALD_PA')
           ewald_param % alpha = get_real()
           SELECT CASE (ewald_param % ewald_TYPE ( 1:3)) 
           CASE DEFAULT
              ewald_param % gmax = get_int()
           CASE ( 'PME')
              ewald_param % ns_max = get_int()
              IF ( test_next() == 'N' ) THEN
                 ewald_param % epsilon = get_real()
              END IF
           END SELECT
        CASE ( 'SET_FILE')
           ilen = 20
           CALL cfield(setup % set_file_name,ilen)
        CASE ( 'INPUT_FI')
           ilen = 20
           CALL cfield(setup % input_file_name,ilen)
        END SELECT
        
! check for trailing rubbish
        CALL read_line
     END DO
     
  END IF
  CALL parser_end
!..end of parsing the input section
     
!..write some information to output
  IF (fistpar % ionode) THEN
     IF ( fistpar % print_level >= 0 ) THEN
        WRITE ( iw, '( A, T71, A )' ) &
             ' FIST| Run type ', ADJUSTR ( setup % run_type )
        WRITE ( iw, '( A, T71, A )' ) &
             ' FIST| Unit type ', ADJUSTR ( setup % unit_type )
        WRITE ( iw, '( A, T78, A )' ) ' FIST| Periodic in X direction ', &
             yn(setup % perd(1))
        WRITE ( iw, '( A, T78, A )' ) ' FIST| Periodic in Y direction ', &
             yn(setup % perd(2))
        WRITE ( iw, '( A, T78, A )' ) ' FIST| Periodic in Z direction ', &
             yn(setup % perd(3))
        
        WRITE ( iw, '( A, T61, A )' ) ' FIST| Set file name', &
             ADJUSTR ( setup % set_file_name )
        WRITE ( iw, '( A, T61, A )' ) ' FIST| Input file name', &
             ADJUSTR ( setup % input_file_name )
        WRITE ( iw, '( A, T76, I5 )' ) &
             ' FIST| Print level ', fistpar % print_level
        WRITE ( iw, '( )' )
     END IF
  END IF
  
END SUBROUTINE read_fist_section

!******************************************************************************

END MODULE fist_input

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!


MODULE energy_input

  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, stop_parser
  USE string_utilities, ONLY : uppercase
  USE termination, ONLY : stop_program 


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_energy_section, energy_min_type, min_info_type

  TYPE min_info_type
    REAL ( dbl ) :: e, de, tol, e_over, e_class, e_thk, e_electro, err
  END TYPE min_info_type

  TYPE energy_min_type
    CHARACTER ( len=8 ) :: mintype
    INTEGER :: iter
    INTEGER :: ndiis
    REAL ( dbl ) :: max_tol
    REAL ( dbl ) :: sd_delta
    TYPE ( min_info_type ) :: min_info
  END TYPE energy_min_type

CONTAINS

!******************************************************************************

!!>---------------------------------------------------------------------------!
!! SECTION: &energy ... &end                                                  !
!!                                                                            !
!!  mintype:  [DIIS, SD]                 METHOD for minimization              !
!!  sd_delta: delta                      delta t for steepest descent         !
!!  max_tol:  max_tol                    maximum tolerance for convergence    !
!!  iter:     max_iter                   maximum number of iterations         !
!!  ndiis:    ndiis                      number of DIIS vectors               !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_energy_section ( energy_min, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( energy_min_type ), INTENT ( OUT ) :: energy_min
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv     

! Locals
  INTEGER :: ierror, ilen, iw
  CHARACTER ( len = 20 ) :: string
  CHARACTER ( len = 5 ) :: label

!------------------------------------------------------------------------------
  iw = globenv % scr

!..defaults
  energy_min % max_tol = 1.0E-6_dbl
  energy_min % iter = 10000
  energy_min % ndiis  = 5
  energy_min % sd_delta = 1._dbl
  energy_min % mintype = 'SD'

! initialize min variables
  energy_min % min_info % e = 0._dbl
  energy_min % min_info % de = 0._dbl
  energy_min % min_info % err = 1.0E+6_dbl
  energy_min % min_info % tol = 1.0E+6_dbl
  energy_min % min_info % e_over = 0._dbl
  energy_min % min_info % e_electro = 0._dbl
  energy_min % min_info % e_thk = 0._dbl
  energy_min % min_info % e_class = 0._dbl

!..parse the input section

  IF ( globenv % ionode ) THEN

     label = '&ENERGY'

     CALL parser_init ( globenv % input_file_name, label, ierror, globenv )
     IF ( ierror /= 0 ) THEN
        CALL stop_program ( 'read_energy_section', 'no input section &ENERGY found' )
     END IF

     CALL read_line()

     DO WHILE ( test_next ( ) /= 'X' )

        ilen = 8
        CALL cfield ( string, ilen )
        CALL uppercase ( string ) 

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error ( )
           CALL stop_program ( 'read_energy_section', 'unknown option' )

        CASE ( 'MINTYPE' )
           ilen = 8
           CALL cfield ( string, ilen )
           CALL uppercase ( string ) 
           energy_min % mintype = string ( 1:ilen )

        CASE ( 'TOL' )
           energy_min % max_tol = get_real ( )

        CASE ( 'ITER' )
           energy_min % iter = get_int ( )

        CASE ( 'NDIIS' )
           energy_min % ndiis = get_int ( )

        CASE ( 'SD_DELTA' )
           energy_min % sd_delta = get_real ( )

        END SELECT

! Get the next line
        CALL read_line

     END DO

  END IF
  CALL parser_end

!..end of parsing the input section
!..write some information to output
  IF ( globenv % print_level >= 0 ) THEN
     WRITE ( iw, '( A )' ) ' ENERGY MINIMIZATION'
     WRITE ( iw, '( A, T71, A )' ) &
          ' ENERGY| Minimization type ', ADJUSTR ( energy_min % mintype )
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' ENERGY| Max Number of iteration ', energy_min % iter
     WRITE ( iw, '( A, T71, G10.4 )' ) &
          ' ENERGY| Max tolerance ', energy_min % max_tol
     WRITE ( iw, '( A, T71, G10.4 )' ) &
          ' ENERGY| delta         ', energy_min % sd_delta
     IF (energy_min % mintype == 'DIIS' ) THEN
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' ENERGY| number of diis vectors ', energy_min % ndiis
     END IF
  END IF

END SUBROUTINE read_energy_section

!******************************************************************************

END MODULE energy_input

!******************************************************************************

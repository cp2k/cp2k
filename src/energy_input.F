!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!


MODULE energy_input

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             test_next
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_energy_section, energy_min_type, min_info_type

  TYPE min_info_type
    REAL (KIND=dp) :: e, de, tol, e_over, e_class, e_thk, e_electro, e_pol, err
    REAL (KIND=dp) :: e_electro_0, e_electro_1, e_electro_2, e_over_ii, e_over_ij             
  END TYPE min_info_type

  TYPE energy_min_type
    CHARACTER ( len=8 ) :: mintype
    INTEGER :: iter
    INTEGER :: ndiis
    REAL (KIND=dp) :: max_tol
    REAL (KIND=dp) :: sd_delta
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


    TYPE(energy_min_type), INTENT(OUT)       :: energy_min
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=20)                        :: string
    CHARACTER(len=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw

!------------------------------------------------------------------------------

  iw = globenv % scr

!..defaults
  energy_min % max_tol = 1.0E-6_dp
  energy_min % iter = 10000
  energy_min % ndiis  = 5
  energy_min % sd_delta = 1.0_dp
  energy_min % mintype = 'SD'

! initialize min variables
  energy_min % min_info % e = 0.0_dp
  energy_min % min_info % de = 0.0_dp
  energy_min % min_info % err = 1.0E+6_dp
  energy_min % min_info % tol = 1.0E+6_dp
  energy_min % min_info % e_over = 0.0_dp
  energy_min % min_info % e_electro = 0.0_dp
  energy_min % min_info % e_electro_0 = 0.0_dp
  energy_min % min_info % e_electro_1 = 0.0_dp
  energy_min % min_info % e_electro_2 = 0.0_dp
  energy_min % min_info % e_thk = 0.0_dp
  energy_min % min_info % e_class = 0.0_dp

!..parse the input section

  label = '&ENERGY'

  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )

  IF ( ierror /= 0 ) THEN

     IF ( globenv % ionode ) THEN
        CALL stop_program ( 'read_energy_section', 'no input section &ENERGY found' )
     END IF

  ELSE

     CALL read_line()

     DO WHILE ( test_next ( ) /= 'X' )

        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string ) 

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error ( )
           CALL stop_program ( 'read_energy_section', 'unknown option' )

        CASE ( 'MINTYPE' )
           ilen = 8
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           energy_min % mintype = string ( 1:ilen )

        CASE ( 'TOL' )
           CALL get_next ( energy_min % max_tol )

        CASE ( 'ITER' )
           CALL get_next ( energy_min % iter )

        CASE ( 'NDIIS' )
           CALL get_next ( energy_min % ndiis )

        CASE ( 'SD_DELTA' )
           CALL get_next ( energy_min % sd_delta )

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

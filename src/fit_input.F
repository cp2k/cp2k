!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!


MODULE fit_input

  USE fit_parameters_types,            ONLY: fit_parameter_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
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
  PUBLIC :: read_fit_section


CONTAINS

!******************************************************************************

!!>---------------------------------------------------------------------------!
!! SECTION: &fit ... &end                                                     !
!!                                                                            !
!!  fit_type: [AMOEBA_RI,AMOEBA_POL,FORCE_POL]  type of fit run
!!  pot_type:    [LJ, WILL]              type of two body potential to fit    !
!!  restart: [INITIAL, ALL]              restart option                       !
!!  data type: [POS_FORCE, POS_VEL]      type of data to read                 !
!!  out unit: [ATOMIC, KELVIN]           output unit                          !
!!  tol:      max_tol                    maximum tolerance (fractional ) !
!!  iter:     max_iter                   maximum number of iterations         !
!!  nsteps:   nsteps                     number of trajectory steps           !
!!  iskip:    iskip                      analyse every iskip step             !
!!  idump:    idump                      write to dump file every idump steps !
!!  iwrite:   iwrite                    write to screen  every iwrite steps !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_fit_section ( fit_param, globenv )


    TYPE(fit_parameter_type), INTENT(OUT)    :: fit_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=20)                        :: string
    CHARACTER(len=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw

!------------------------------------------------------------------------------

  iw = globenv % scr

!..defaults
  fit_param % max_tol = 1.0E-6_dbl
  fit_param % max_iter = 10000
  fit_param % nsteps = 1
  fit_param % iskip = 1
  fit_param % iwrite = 1
  fit_param % pot_type = 'LJ'
  fit_param % fit_type = 'AMOEBA_RI'
  fit_param % datatype = 'POS_VEL'
  fit_param % outunit = 'ATOMIC'
  fit_param % idump = 10
  fit_param % restart = 'INITIAL'

!..parse the input section

  IF ( globenv % ionode ) THEN

     label = '&FIT'

     CALL parser_init ( globenv % input_file_name, globenv )
     CALL search_label(label,ierror,ignore_case=.TRUE.)
     IF ( ierror /= 0 ) THEN
        CALL stop_program ( 'read_fit_section', 'no input section &FIT found' )
     END IF

     CALL read_line()

     DO WHILE ( test_next ( ) /= 'X' )

        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string ) 

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error ( )
           CALL stop_program ( 'read_fit_section', 'unknown option' )

        CASE ( 'FIT_TYPE' )
           ilen = 10
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           fit_param % fit_type = string ( 1:ilen )

        CASE ( 'POT_TYPE' )
           ilen = 6
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           fit_param % pot_type = string ( 1:ilen )

        CASE ( 'DATA' )
           ilen = 10
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           fit_param % datatype = string ( 1:ilen )

        CASE ( 'RESTART' )
           ilen = 7
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           fit_param % restart = string ( 1:ilen )

        CASE ( 'OUTUNIT' )
           ilen = 6
           CALL get_next ( string, ilen )
           CALL uppercase ( string ) 
           fit_param % outunit = string ( 1:ilen )

        CASE ( 'TOL' )
           CALL get_next ( fit_param % max_tol )

        CASE ( 'ITER' )
           CALL get_next ( fit_param % max_iter )

        CASE ( 'NSTEPS' )
           CALL get_next ( fit_param % nsteps )

        CASE ( 'IDUMP' )
           CALL get_next ( fit_param % idump )

        CASE ( 'ISKIP' )
           CALL get_next ( fit_param % iskip )

        CASE ( 'IWRITE' )
           CALL get_next ( fit_param % iwrite )

        END SELECT

! Get the next line
        CALL read_line

     END DO

  END IF
  CALL parser_end

!..end of parsing the input section
!..write some information to output
  IF ( globenv % print_level >= 0 ) THEN
     WRITE ( iw, '( A )' ) ' FIT| Fit Protocol '
     WRITE ( iw, '( A, T71, A )' ) &
          ' FIT| Fit type ', ADJUSTR ( fit_param % fit_type )
     WRITE ( iw, '( A, T71, A )' ) &
          ' FIT| 2 body potential type ', ADJUSTR ( fit_param % pot_type )
     WRITE ( iw, '( A, T71, A )' ) &
          ' FIT| Data type ', ADJUSTR ( fit_param % datatype )
     WRITE ( iw, '( A, T71, A )' ) &
          ' FIT| Restart type', ADJUSTR ( fit_param % restart )
     WRITE ( iw, '( A, T71, A )' ) &
          ' FIT| Out Units ', ADJUSTR ( fit_param % outunit )
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' FIT| Number of steps ', fit_param % nsteps
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' FIT| Iskip ', fit_param % iskip
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' FIT| Idump ', fit_param % idump
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' FIT| iwrite ', fit_param % iwrite
     WRITE ( iw, '( A, T71, G10.4 )' ) &
          ' FIT| Max tolerance ', fit_param % max_tol
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' FIT| Max iterations ', fit_param % max_iter
  END IF

END SUBROUTINE read_fit_section

!******************************************************************************

END MODULE fit_input

!******************************************************************************

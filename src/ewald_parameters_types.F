!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/ewald_parameters_types [1.0] *
!!
!!   NAME
!!     ewald_parameters_types
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

MODULE ewald_parameters_types

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: uppercase

  TYPE ewald_parameters_type
    CHARACTER( LEN = 12 ) :: ewald_type  ! type of ewald
    CHARACTER( LEN = 12 ) :: unit_type   !  unit_type
    REAL ( dbl ) :: alpha                ! ewald alpha
    REAL ( dbl ) :: eps0                 ! permittivity constant
    REAL ( dbl ) :: epsilon              ! tolerance for small grid (PME)
    INTEGER   :: gmax                    ! max Miller index
    INTEGER   :: ns_max                  ! # grid points for small grid (PME)
    INTEGER   :: o_spline                ! order of spline (SPME)
    INTEGER   :: gtot                    ! total number of g-vectors
    INTEGER   :: group                   ! parallel group
  END TYPE ewald_parameters_type

  PRIVATE
  PUBLIC :: ewald_parameters_type, read_ewald_section
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!I---------------------------------------------------------------------------I
!!I SECTION: &ewald ... &end                                                  I
!!I                                                                           I
!!I units            [atomic, angstrom]                                       I
!!I Ewald_type       [ewald,pme,spme]                                         I
!!I Ewald_param      alpha[,gmax,ns_max,epsilon]                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************
SUBROUTINE read_ewald_section ( ewald_param, globenv )
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw, n
    LOGICAL                                  :: found

  iw = globenv % scr
! Defaults
  ewald_param % alpha = 0.4_dbl
  ewald_param % gmax = 10
  ewald_param % o_spline = 5
  ewald_param % ns_max = 10
  ewald_param % epsilon = 1.e-6_dbl
  ewald_param % ewald_type = 'NONE'
  ewald_param % unit_type = 'ANGSTROM'
  iw = globenv % scr
!..parse the input section
  label = '&EWALD'
  CALL parser_init(globenv % input_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     IF (globenv % ionode) &
        WRITE ( iw, '( a )' ) ' No input section &EWALD found '
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser( 'read_ewald_section','unknown option')
        CASE ( 'UNITS')
           ilen = 20
           CALL get_next(ewald_param % unit_type,ILEN)
           CALL uppercase(ewald_param % unit_type )
        CASE ( 'EWALD_TY')
           ilen=20
           CALL get_next(string,ILEN)
           CALL uppercase ( string )
           SELECT CASE(string)
           CASE DEFAULT
              CALL p_error()
              CALL stop_parser( 'read_ewald_section', &
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
           CALL get_next ( ewald_param % alpha )
           SELECT CASE (ewald_param % ewald_type )
           CASE DEFAULT
              CALL get_next ( ewald_param % gmax )
           CASE ( 'PME')
              CALL get_next ( ewald_param % ns_max )
              IF ( test_next() == 'N' ) THEN
                 CALL get_next ( ewald_param % epsilon )
              END IF
           CASE ( 'SPME')
              CALL get_next ( ewald_param % gmax )
              CALL get_next ( ewald_param % o_spline )
           END SELECT
! check for trailing rubbish
        END SELECT
        CALL read_line
     END DO
  END IF
  CALL parser_end

!..end of parsing the input section
  IF ( globenv % ionode ) THEN
    IF ( ewald_param % ewald_type /= 'NONE' ) THEN

       WRITE ( iw, '( A,T67,A14 )' ) ' Ewald| Summation is done by:', &
            ADJUSTR(ewald_param % ewald_type)
       WRITE ( iw, '( A,A18,A,T71,F10.4 )' ) ' Ewald| Alpha parameter [',ewald_param % unit_type,']', &
            ewald_param % alpha

       SELECT CASE ( ewald_param % ewald_type )
       CASE DEFAULT
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_param % gmax
       CASE ( 'PME')
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' EWALD| Max small-grid points (input) ', ewald_param % ns_max
          WRITE ( iw, '( A,T71,E10.4 )' ) &
               ' EWALD| Gaussian tolerance (input) ', ewald_param % epsilon
       CASE ( 'SPME' )
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_param % gmax
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| Spline interpolation order ', ewald_param % o_spline
       END SELECT
       CALL convert_to_cp2k_units ( ewald_param % unit_type,              &
                                    length = ewald_param % alpha,         &
                                    l_power = -1 )

       CALL convert_to_cp2k_units ( eps0 = ewald_param % eps0 )
    ELSE

       WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'

    END IF
  ELSE
    iw = -1
  END IF
END SUBROUTINE read_ewald_section

!******************************************************************************
END MODULE ewald_parameters_types
!!*****

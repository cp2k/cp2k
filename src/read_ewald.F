!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/read_ewald [1.0] *
!!
!!   NAME
!!     read_ewald
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

MODULE read_ewald

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_set
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
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
  PUBLIC :: read_ewald_section
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
SUBROUTINE read_ewald_section ( ewald_env, globenv )
  IMPLICIT NONE
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
! locals
  LOGICAL :: found, setf, datf, topf
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source, group
  CHARACTER ( len = default_sting_length ) :: string, str2, label, &
                                              ewald_type, unit_type
  REAL (KIND=dp) :: eps0, alpha, epsilo
  INTEGER :: gmax, o_spline, ns_max
  iw = globenv % scr
! Defaults
   alpha = 0.4_dp
   gmax = 10
   o_spline = 5
   ns_max = 10
   epsilon = 1.e-6_dp
   ewald_type = 'NONE'
   unit_type = 'ANGSTROM'
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
           CALL get_next( unit_type,ILEN)
           CALL uppercase( unit_type )
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
              ewald_type = 'ewald'
              CALL uppercase( ewald_type )
           CASE( 'PME')
              ewald_type = 'pme'
              CALL uppercase( ewald_type )
           CASE( 'SPME')
              ewald_type = 'spme'
              CALL uppercase( ewald_type )
           END SELECT
! if no type specified, assume ewald
        CASE ( 'EWALD_PA')
           CALL get_next ( alpha )
           SELECT CASE ( ewald_type )
           CASE DEFAULT
              CALL get_next ( gmax )
           CASE ( 'PME')
              CALL get_next ( ns_max )
              IF ( test_next() == 'N' ) THEN
                 CALL get_next ( epsilon )
              END IF
           CASE ( 'SPME')
              CALL get_next ( gmax )
              CALL get_next ( o_spline )
           END SELECT
! check for trailing rubbish
        END SELECT
        CALL read_line
     END DO
  END IF
  CALL parser_end
  

!..end of parsing the input section
  IF ( globenv % ionode ) THEN
    IF ( ewald_env % ewald_type /= 'NONE' ) THEN

       WRITE ( iw, '( A,T67,A14 )' ) ' Ewald| Summation is done by:', &
            ADJUSTR(ewald_env % ewald_type)
       WRITE ( iw, '( A,A18,A,T71,F10.4 )' ) ' Ewald| Alpha parameter [',ewald_env % unit_type,']', &
            ewald_env % alpha

       SELECT CASE ( ewald_env % ewald_type )
       CASE DEFAULT
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_env % gmax
       CASE ( 'PME')
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' EWALD| Max small-grid points (input) ', ewald_env % ns_max
          WRITE ( iw, '( A,T71,E10.4 )' ) &
               ' EWALD| Gaussian tolerance (input) ', ewald_env % epsilon
       CASE ( 'SPME' )
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_env % gmax
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| Spline interpolation order ', ewald_env % o_spline
       END SELECT
       CALL convert_to_cp2k_units ( ewald_env % unit_type,              &
                                    length = ewald_env % alpha,         &
                                    l_power = -1 )

       CALL convert_to_cp2k_units ( eps0 = ewald_env % eps0 )
    ELSE

       WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'

    END IF
  ELSE
    iw = -1
  END IF
  CALL ewald_env_set ( ewald_env, eps0 = eps0, alpha = alpha, &
                       ewald_type = ewald_type, unit_type = unit_type, &
                       o_spline = o_spline, gmax = gmax, ns_max = ns_max, &
                       epsilon = epsilon )
END SUBROUTINE read_ewald_section

!******************************************************************************
END MODULE read_ewald
!!*****

!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/ewald_environment_types [1.0] *
!!
!!   NAME
!!     ewald_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM NOV-30-2003
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE ewald_environment_types
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_retain,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl, default_string_length
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
!
  TYPE ewald_environment_type
    PRIVATE
    INTEGER :: id_nr, ref_count 
    CHARACTER( LEN = default_string_length ) :: ewald_type  ! type of ewald
    CHARACTER( LEN = default_string_length ) :: unit_type   !  unit_type
    REAL ( dbl ) :: alpha               ! ewald alpha
    REAL ( dbl ) :: eps0                ! permittivity constant
    REAL ( dbl ) :: epsilon             ! tolerance for small grid (PME)
    INTEGER   :: gmax                   ! max Miller index
    INTEGER   :: ns_max                 ! # grid points for small grid (PME)
    INTEGER   :: o_spline               ! order of spline (SPME)
    TYPE ( cp_para_env_type ), POINTER :: para_env
  END TYPE ewald_environment_type

! *****************************************************************************

!!****s* ewald_environment_types/ewald_environment_p_type [1.0] *
!!
!!   NAME
!!     ewald_environment_p_type
!!
!!   FUNCTION
!!     to build arrays of pointers
!!
!!   ATTRIBUTES
!!     - ewald_env: the pointer to the ewald_env
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     11/03
!!
!!   SOURCE
  !***************************************************************************
  TYPE ewald_environment_p_type
     TYPE(ewald_environment_type), POINTER :: ewald_env
  END TYPE ewald_environment_p_type
!!***
  !****************************************************************************

! *** Public data types ***

  PUBLIC :: ewald_environment_type,        & 
            ewald_environment_p_type

! *** Public subroutines ***

  PUBLIC :: ewald_env_get,&
            ewald_env_set,&
            ewald_env_create, &
            ewald_env_retain, & 
            ewald_env_release, &
            read_ewald_section

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='ewald_environment_types'
  INTEGER, PRIVATE, SAVE :: last_ewald_env_id_nr=0

!***
! *****************************************************************************

 CONTAINS

! *****************************************************************************

  SUBROUTINE ewald_env_get( ewald_env, ewald_type, unit_type, alpha, eps0, &
                            epsilon, gmax, ns_max, o_spline, group,  &
                            para_env, id_nr, error )
    IMPLICIT NONE

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE ( ewald_environment_type ), POINTER       :: ewald_env
    INTEGER, INTENT(out), OPTIONAL                 :: id_nr
    CHARACTER ( default_string_length ), OPTIONAL  :: ewald_type 
    CHARACTER ( default_string_length ), OPTIONAL  :: unit_type 
    REAL ( dbl ), OPTIONAL                         :: alpha
    REAL ( dbl ), OPTIONAL                         :: eps0
    REAL ( dbl ), OPTIONAL                         :: epsilon
    INTEGER, OPTIONAL                              :: gmax
    INTEGER, OPTIONAL                              :: ns_max
    INTEGER, OPTIONAL                              :: o_spline
    INTEGER, OPTIONAL                              :: group
    TYPE ( cp_para_env_type ), POINTER, OPTIONAL   :: para_env
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='ewald_env_get',&
                                   routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)

    IF ( PRESENT ( ewald_type ) ) ewald_type = ewald_env % ewald_type
    IF ( PRESENT ( unit_type ) ) unit_type = ewald_env % unit_type
    IF ( PRESENT ( alpha ) ) alpha = ewald_env % alpha
    IF ( PRESENT ( eps0 ) ) eps0 = ewald_env % eps0
    IF ( PRESENT ( epsilon ) ) epsilon = ewald_env % epsilon
    IF ( PRESENT ( gmax ) ) gmax = ewald_env % gmax
    IF ( PRESENT ( ns_max ) ) ns_max = ewald_env % ns_max
    IF ( PRESENT ( o_spline ) ) o_spline = ewald_env % o_spline
    IF ( PRESENT ( group ) ) group = ewald_env % para_env % group
    IF ( PRESENT ( para_env ) ) para_env => ewald_env % para_env

  END SUBROUTINE ewald_env_get
!
! *****************************************************************************

  SUBROUTINE init_ewald_env ( ewald_env, para_env, error )

!   Purpose: Initialise the FIST environment.

!   ***************************************************************************

    TYPE( ewald_environment_type ), INTENT(OUT)     :: ewald_env
    TYPE ( cp_para_env_type ), POINTER, OPTIONAL :: para_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   ---------------------------------------------------------------------------

    ewald_env % ref_count=1
    last_ewald_env_id_nr=last_ewald_env_id_nr+1
    ewald_env%id_nr=last_ewald_env_id_nr
! put in some safe defaults
    ewald_env % ewald_type = "NONE"
    ewald_env % unit_type = "AU"
    IF ( PRESENT ( para_env ) ) THEN
      CALL cp_para_env_retain ( para_env )
      ewald_env % para_env => para_env
    END IF 

  END SUBROUTINE init_ewald_env

! *****************************************************************************

  SUBROUTINE ewald_env_set( ewald_env, ewald_type, unit_type, alpha, eps0, &
                            epsilon, gmax, ns_max, o_spline, &
                            para_env, id_nr, error )
    IMPLICIT NONE

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE ( ewald_environment_type ), POINTER       :: ewald_env
    INTEGER, INTENT(out), OPTIONAL                 :: id_nr
    CHARACTER ( default_string_length ), OPTIONAL  :: ewald_type 
    CHARACTER ( default_string_length ), OPTIONAL  :: unit_type 
    REAL ( dbl ), OPTIONAL                         :: alpha
    REAL ( dbl ), OPTIONAL                         :: eps0
    REAL ( dbl ), OPTIONAL                         :: epsilon
    INTEGER, OPTIONAL                              :: gmax
    INTEGER, OPTIONAL                              :: ns_max
    INTEGER, OPTIONAL                              :: o_spline
    TYPE ( cp_para_env_type ), POINTER, OPTIONAL   :: para_env
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='ewald_env_set',&
                                   routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)
    
    IF ( .NOT. failure ) THEN
      IF ( PRESENT ( ewald_type ) ) ewald_env % ewald_type = ewald_type
      IF ( PRESENT ( unit_type ) ) ewald_env % unit_type = unit_type
      IF ( PRESENT ( alpha ) ) ewald_env % alpha = alpha
      IF ( PRESENT ( eps0 ) ) ewald_env % eps0 = eps0
      IF ( PRESENT ( epsilon ) ) ewald_env % epsilon = epsilon
      IF ( PRESENT ( gmax ) ) ewald_env % gmax = gmax
      IF ( PRESENT ( ns_max ) ) ewald_env % ns_max = ns_max
      IF ( PRESENT ( o_spline ) ) ewald_env % o_spline = o_spline
      IF ( PRESENT ( para_env ) )ewald_env %  para_env => para_env
    ENDIF

!***************************************************************************
  END SUBROUTINE ewald_env_set
!***************************************************************************
!
!****f* ewald_environment_types/ewald_env_create [1.0] *
!
!   NAME
!     ewald_env_create
!
!   SYNOPSIS
!     Subroutine ewald_env_create(ewald_env, error)
!       Type(ewald_environment_type), Pointer:: ewald_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine ewald_env_create
!
!   FUNCTION
!     allocates and intitializes a ewald_env
!
!   ARGUMENTS
!     - ewald_env: the object to create
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE ewald_env_create ( ewald_env, para_env, error )
  TYPE(ewald_environment_type), POINTER :: ewald_env
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='ewald_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE( ewald_env, stat=stat )
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_ewald_env ( ewald_env, para_env = para_env, error=error)
  END IF
END SUBROUTINE ewald_env_create
!***************************************************************************

!****f* ewald_environment_types/ewald_env_retain [1.0] *
!
!   NAME
!     ewald_env_retain
!
!   SYNOPSIS
!     Subroutine ewald_env_retain(ewald_env, error)
!       Type(ewald_environment_type), Pointer:: ewald_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine ewald_env_retain
!
!   FUNCTION
!     retains the given ewald_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - ewald_env: the object to retain
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE ewald_env_retain(ewald_env,error)
  TYPE(ewald_environment_type), POINTER :: ewald_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='ewald_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
!!  
  CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ewald_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ewald_env%ref_count=ewald_env%ref_count+1
  END IF
END SUBROUTINE ewald_env_retain
!***************************************************************************

!****f* ewald_environment_types/ewald_env_release [1.0] *
!
!   NAME
!     ewald_env_release
!
!   SYNOPSIS
!     Subroutine ewald_env_release(ewald_env, error)
!       Type(ewald_environment_type), Pointer:: ewald_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine ewald_env_release
!
!   FUNCTION
!     releases the given ewald_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - ewald_env: the object to release
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE ewald_env_release(ewald_env,error)
  TYPE(ewald_environment_type), POINTER :: ewald_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='ewald_env_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: i, stat
  failure=.FALSE.

  IF (ASSOCIATED(ewald_env)) THEN
     CPPrecondition(ewald_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ewald_env%ref_count=ewald_env%ref_count-1
     IF (ewald_env%ref_count<1) THEN
       CALL cp_para_env_release ( ewald_env % para_env, error )
        DEALLOCATE (  ewald_env, stat = stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)       
     ENDIF
  END IF      
  NULLIFY(ewald_env)
END SUBROUTINE ewald_env_release
!***************************************************************************




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
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source
  CHARACTER ( len = default_string_length ) :: string, str2, label
                                             
  iw = globenv % scr
  CALL convert_to_cp2k_units ( eps0 = ewald_env % eps0 )
! Defaults
   ewald_env % alpha = 0.4_dbl
   ewald_env % gmax = 10
   ewald_env % o_spline = 5
   ewald_env % ns_max = 10
   ewald_env % epsilon = 1.e-6_dbl
   ewald_env % ewald_type = 'NONE'
   ewald_env % unit_type = 'ANGSTROM'
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
           CALL stop_parser( 'read_ewald_section','unknown option'//trim(string))
        CASE ( 'UNITS', 'UNIT')
           ilen = 20
           CALL get_next( ewald_env % unit_type,ILEN)
           CALL uppercase( ewald_env % unit_type )
        CASE ( 'EWALD_TY')
           ilen=20
           CALL get_next(string,ILEN)
           CALL uppercase ( string )
           SELECT CASE(string)
           CASE DEFAULT
              CALL p_error()
              CALL stop_parser( 'read_ewald_section', &
                   'unknown option for Ewald type' )
           CASE( 'NONE')
              ewald_env % ewald_type = 'none'
              CALL uppercase( ewald_env % ewald_type )
           CASE( 'EWALD')
              ewald_env % ewald_type = 'ewald'
              CALL uppercase( ewald_env % ewald_type )
           CASE( 'PME')
              ewald_env % ewald_type = 'pme'
              CALL uppercase( ewald_env % ewald_type )
           CASE( 'SPME')
              ewald_env % ewald_type = 'spme'
              CALL uppercase( ewald_env % ewald_type )
           END SELECT
! if no type specified, assume ewald
        CASE ( 'EWALD_PA')
           CALL get_next ( ewald_env % alpha )
           SELECT CASE ( ewald_env % ewald_type )
           CASE DEFAULT
              CALL get_next ( ewald_env % gmax )
           CASE ( 'PME')
              CALL get_next ( ewald_env % ns_max )
              IF ( test_next() == 'N' ) THEN
                 CALL get_next ( ewald_env % epsilon )
              END IF
           CASE ( 'SPME')
              CALL get_next ( ewald_env % gmax )
              CALL get_next ( ewald_env % o_spline )
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

    ELSE

       WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'

    END IF
  ELSE
    iw = -1
  END IF
  CALL convert_to_cp2k_units ( ewald_env % unit_type,              &
                               length = ewald_env % alpha,         &
                               l_power = -1 )
END SUBROUTINE read_ewald_section

!!*****
END MODULE ewald_environment_types

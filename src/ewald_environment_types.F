!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
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
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_para_env,                     ONLY: cp_para_env_retain,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp, default_string_length
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get,&
                                             section_vals_retain,&
                                             section_vals_release
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme,&
                                             ewald_name_types
  USE physcon, only: angstrom
  USE global_types,                    ONLY: DEBUG,&
                                             HIGH,&
                                             LOW,&
                                             MEDIUM,& 
                                             SILENT
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

! *** Public parameters ***

  PUBLIC :: DEBUG,HIGH,LOW,MEDIUM,SILENT

  TYPE ewald_environment_type
    PRIVATE
    INTEGER :: id_nr, ref_count 
    INTEGER :: ewald_type                 ! type of ewald
    CHARACTER( LEN = default_string_length ) :: unit_type   !  unit_type TO BE DELETED when new input will work!
    REAL (KIND=dp) :: alpha               ! ewald alpha
    REAL (KIND=dp) :: eps0                ! permittivity constant
    REAL (KIND=dp) :: epsilon             ! tolerance for small grid (PME)
    INTEGER   :: gmax                   ! max Miller index
    INTEGER   :: ns_max                 ! # grid points for small grid (PME)
    INTEGER   :: o_spline               ! order of spline (SPME)
    TYPE ( cp_para_env_type ), POINTER :: para_env
    TYPE(section_vals_type), pointer :: poisson_section
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
                            para_env, id_nr, poisson_section, error )
    IMPLICIT NONE

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE ( ewald_environment_type ), POINTER       :: ewald_env
    INTEGER, INTENT(OUT), OPTIONAL                 :: id_nr
    INTEGER, OPTIONAL  :: ewald_type 
    CHARACTER ( default_string_length ), OPTIONAL  :: unit_type 
    REAL (KIND=dp), OPTIONAL                         :: alpha
    REAL (KIND=dp), OPTIONAL                         :: eps0
    REAL (KIND=dp), OPTIONAL                         :: epsilon
    INTEGER, OPTIONAL                              :: gmax
    INTEGER, OPTIONAL                              :: ns_max
    INTEGER, OPTIONAL                              :: o_spline
    INTEGER, OPTIONAL                              :: group
    TYPE ( cp_para_env_type ), POINTER, OPTIONAL   :: para_env
    TYPE(section_vals_type), POINTER, optional :: poisson_section
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='ewald_env_get',&
                                   routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)

    IF ( PRESENT ( id_nr ) ) id_nr = ewald_env % id_nr
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
    IF (PRESENT(poisson_section)) poisson_section => ewald_env%poisson_section

  END SUBROUTINE ewald_env_get
!
! *****************************************************************************

  SUBROUTINE ewald_env_set( ewald_env, ewald_type, unit_type, alpha, eps0, &
                            epsilon, gmax, ns_max, o_spline, &
                            para_env, id_nr, poisson_section, error )
    IMPLICIT NONE

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE ( ewald_environment_type ), POINTER       :: ewald_env
    INTEGER, INTENT(IN), OPTIONAL                  :: id_nr
    INTEGER, OPTIONAL  :: ewald_type 
    CHARACTER ( default_string_length ), OPTIONAL  :: unit_type 
    REAL (KIND=dp), OPTIONAL                         :: alpha
    REAL (KIND=dp), OPTIONAL                         :: eps0
    REAL (KIND=dp), OPTIONAL                         :: epsilon
    INTEGER, OPTIONAL                              :: gmax
    INTEGER, OPTIONAL                              :: ns_max
    INTEGER, OPTIONAL                              :: o_spline
    TYPE ( cp_para_env_type ), POINTER, OPTIONAL   :: para_env
    TYPE(section_vals_type), POINTER, optional     :: poisson_section
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error

    LOGICAL :: failure
    TYPE(section_vals_type), POINTER    :: ewald_section
    CHARACTER(len=*), PARAMETER :: routineN='ewald_env_set',&
                                   routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ewald_env),cp_failure_level,routineP,error,failure)
    
    IF ( .NOT. failure ) THEN
      IF ( PRESENT ( id_nr ) ) ewald_env % id_nr = id_nr
      IF ( PRESENT ( ewald_type ) ) ewald_env % ewald_type = ewald_type
      IF ( PRESENT ( unit_type ) ) ewald_env % unit_type = unit_type
      IF ( PRESENT ( alpha ) ) ewald_env % alpha = alpha
      IF ( PRESENT ( eps0 ) ) ewald_env % eps0 = eps0
      IF ( PRESENT ( epsilon ) ) ewald_env % epsilon = epsilon
      IF ( PRESENT ( gmax ) ) ewald_env % gmax = gmax
      IF ( PRESENT ( ns_max ) ) ewald_env % ns_max = ns_max
      IF ( PRESENT ( o_spline ) ) ewald_env % o_spline = o_spline
      IF ( PRESENT ( para_env ) )ewald_env %  para_env => para_env
      IF (PRESENT(poisson_section)) THEN
         CALL section_vals_retain(poisson_section,error=error)
         CALL section_vals_release(ewald_env%poisson_section,error=error)
         ewald_env%poisson_section => poisson_section
         ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD",&
              error=error)
         !FM call read_ewald_section1(ewald_env, ewald_section, error=error)
      END IF
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
  TYPE(unit_convert_type), POINTER         :: units

  failure=.FALSE.
  
  ALLOCATE( ewald_env, stat=stat )
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    NULLIFY(units)
    ewald_env % ref_count=1
    last_ewald_env_id_nr=last_ewald_env_id_nr+1
    ewald_env%id_nr=last_ewald_env_id_nr
    CALL get_cp2k_units ( units )
    ! Defaults
    ewald_env % eps0 = units%eps0
    ewald_env % alpha = 0.4_dp
    ewald_env % gmax = 10
    ewald_env % o_spline = 5
    ewald_env % ns_max = 10
    ewald_env % epsilon = 1.e-6_dp
    ewald_env % ewald_type = do_ewald_none
    ewald_env % unit_type = 'ANGSTROM'
    nullify(ewald_env%poisson_section)
    CALL cp_para_env_retain ( para_env )
    ewald_env % para_env => para_env
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
        CALL section_vals_release(ewald_env%poisson_section,error=error)
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
SUBROUTINE read_ewald_section ( ewald_env, globenv, error )
  IMPLICIT NONE
  ! Arguments
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( global_environment_type ), pointer :: globenv
  TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error 
  ! Local Variables
  LOGICAL :: found, setf, datf, topf
  INTEGER :: ierror, ilen, ia, ie, i, j, n, iw, source
  CHARACTER ( len = default_string_length ) :: string, str2, label
  TYPE(unit_convert_type), POINTER         :: units
  TYPE(section_vals_type), POINTER         :: ewald_section  

  iw = globenv % scr
  NULLIFY(ewald_section, units)
  CALL get_cp2k_units ( units )
  ewald_section => section_vals_get_subs_vals(globenv%input_file,&
       "FORCE_EVAL%MM%POISSON_MM%EWALD",error=error)
  CALL read_ewald_section1(ewald_env, ewald_section, error)

  !..end of parsing the input section
  IF ( globenv % ionode .AND. globenv % print_level .NE. SILENT) THEN
     IF ( ewald_env % ewald_type /= do_ewald_none ) THEN

        WRITE ( iw, '( A,T67,A14 )' ) ' Ewald| Summation is done by:', &
             ADJUSTR(ewald_name_types(ewald_env % ewald_type)(1:14))
        WRITE ( iw, '( A,A18,A,T71,F10.4 )' ) ' Ewald| Alpha parameter [','ANGSTROM',']', &
             ewald_env % alpha /  angstrom

        SELECT CASE ( ewald_env % ewald_type )
        CASE DEFAULT
           WRITE ( iw, '( A,T71,I10 )' ) &
                ' Ewald| G-space max. Miller index', ewald_env % gmax
        CASE ( do_ewald_pme)
           WRITE ( iw, '( A,T71,I10 )' ) &
                ' EWALD| Max small-grid points (input) ', ewald_env % ns_max
           WRITE ( iw, '( A,T71,E10.4 )' ) &
                ' EWALD| Gaussian tolerance (input) ', ewald_env % epsilon
        CASE ( do_ewald_spme )
           WRITE ( iw, '( A,T71,I10 )' ) &
                ' Ewald| G-space max. Miller index', ewald_env % gmax
           WRITE ( iw, '( A,T71,I10 )' ) &
                ' Ewald| Spline interpolation order ', ewald_env % o_spline
        END SELECT

     ELSE

        WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'

     END IF
  END IF

END SUBROUTINE read_ewald_section

!
! New routine for the reading of the EWALD section with the new input style
!

SUBROUTINE read_ewald_section1(ewald_env, ewald_section, error)
  IMPLICIT NONE
  ! Arguments
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE(section_vals_type), POINTER         :: ewald_section  
  TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error 
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'read_ewald_section1', &
       routineP = moduleN//':'//routineN
  LOGICAL :: failure
  INTEGER :: n_rep

  failure = .FALSE.
  IF (.NOT.failure) THEN
     CALL section_vals_val_get(ewald_section,"EWALD_TYPE",i_val=ewald_env%ewald_type,error=error)
     CALL section_vals_val_get(ewald_section,"ALPHA",r_val=ewald_env%alpha,error=error)
     CALL section_vals_val_get(ewald_section,"GMAX",i_val=ewald_env%gmax,error=error)
     CALL section_vals_val_get(ewald_section,"NS_MAX",i_val=ewald_env%ns_max,error=error)
     CALL section_vals_val_get(ewald_section,"O_SPLINE",i_val=ewald_env%o_spline,error=error)
     CALL section_vals_val_get(ewald_section,"EPSILON",r_val=ewald_env%epsilon,error=error)
  END IF

END SUBROUTINE read_ewald_section1

END MODULE ewald_environment_types

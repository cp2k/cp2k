!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
!!     JGH FEB-13-2007 : Distributed/replicated realspace grids
!!
!***
!******************************************************************************

MODULE ewald_environment_types
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE global_types,                    ONLY: DEBUG,&
                                             HIGH,&
                                             LOW,&
                                             MEDIUM,&
                                             SILENT
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme,&
                                             do_ewald_ewald
  USE input_cp2k_poisson,              ONLY: create_ewald_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE physcon,                         ONLY: angstrom
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

! *** Public parameters ***

  PUBLIC :: DEBUG,HIGH,LOW,MEDIUM,SILENT

  TYPE ewald_environment_type
    PRIVATE
    INTEGER :: id_nr, ref_count
    INTEGER :: ewald_type                 ! type of ewald
    !  unit_type TO BE DELETED when new input will work!
    CHARACTER( LEN = default_string_length ) :: unit_type  
    REAL (KIND=dp) :: alpha               ! ewald alpha
    REAL (KIND=dp) :: eps0                ! permittivity constant
    REAL (KIND=dp) :: epsilon             ! tolerance for small grid (PME)
    INTEGER   :: gmax(3)                  ! max Miller index
    INTEGER   :: ns_max                   ! # grid points for small grid (PME)
    INTEGER   :: o_spline                 ! order of spline (SPME)
    INTEGER   :: rs_distributed_id        ! realspace distribution
    TYPE ( cp_para_env_type ), POINTER   :: para_env
    TYPE(section_vals_type), POINTER     :: poisson_section
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
!***
  !***************************************************************************
  TYPE ewald_environment_p_type
     TYPE(ewald_environment_type), POINTER :: ewald_env
  END TYPE ewald_environment_p_type
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_environment_types'
  INTEGER, PRIVATE, SAVE :: last_ewald_env_id_nr=0

! *****************************************************************************

 CONTAINS

! *****************************************************************************

  SUBROUTINE ewald_env_get( ewald_env, ewald_type, unit_type, alpha, eps0, &
                            epsilon, gmax, ns_max, o_spline, rs_distributed_id, &
                            group, para_env, id_nr, poisson_section, error )

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    INTEGER, OPTIONAL                        :: ewald_type
    CHARACTER(default_string_length), &
      OPTIONAL                               :: unit_type
    REAL(KIND=dp), OPTIONAL                  :: alpha, eps0, epsilon
    INTEGER, OPTIONAL                        :: gmax(3), ns_max, o_spline, &
                                                rs_distributed_id, group
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: poisson_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    IF ( PRESENT ( rs_distributed_id ) ) &
       rs_distributed_id = ewald_env % rs_distributed_id
    IF ( PRESENT ( group ) ) group = ewald_env % para_env % group
    IF ( PRESENT ( para_env ) ) para_env => ewald_env % para_env
    IF ( PRESENT ( poisson_section ) ) poisson_section => ewald_env%poisson_section

  END SUBROUTINE ewald_env_get
!
! *****************************************************************************

  SUBROUTINE ewald_env_set( ewald_env, ewald_type, unit_type, alpha, eps0, &
                            epsilon, gmax, ns_max, o_spline, rs_distributed_id, &
                            para_env, id_nr, poisson_section, error )

!   Purpose: Get the EWALD environment.

!   ***************************************************************************

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    INTEGER, OPTIONAL                        :: ewald_type
    CHARACTER(default_string_length), &
      OPTIONAL                               :: unit_type
    REAL(KIND=dp), OPTIONAL                  :: alpha, eps0, epsilon
    INTEGER, OPTIONAL                        :: gmax(3), ns_max, o_spline, &
                                                rs_distributed_id
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: poisson_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
      IF ( PRESENT ( rs_distributed_id ) ) &
         ewald_env % rs_distributed_id = rs_distributed_id
      IF ( PRESENT ( para_env ) )ewald_env %  para_env => para_env
      IF (PRESENT(poisson_section)) THEN
         CALL section_vals_retain(poisson_section,error=error)
         CALL section_vals_release(ewald_env%poisson_section,error=error)
         ewald_env%poisson_section => poisson_section
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
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
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
    ewald_env % gmax = -1
    ewald_env % o_spline = 5
    ewald_env % rs_distributed_id = -1
    ewald_env % ns_max = 10
    ewald_env % epsilon = 1.e-6_dp
    ewald_env % ewald_type = do_ewald_none
    ewald_env % unit_type = 'ANGSTROM'
    NULLIFY(ewald_env%poisson_section)
    CALL cp_para_env_retain ( para_env ,error=error)
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
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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

SUBROUTINE read_ewald_section ( ewald_env, ewald_section, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(section_vals_type), POINTER         :: ewald_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_ewald_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dummy
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(unit_convert_type), POINTER         :: units

    INTEGER, DIMENSION(:), POINTER :: gmax_read

    NULLIFY( units, enum, keyword, section)

    CALL get_cp2k_units ( units )

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(ewald_section,&
            "EWALD_TYPE",i_val=ewald_env%ewald_type,error=error)
       CALL section_vals_val_get(ewald_section,&
            "ALPHA",r_val=ewald_env%alpha,error=error)

       ! we have no defaults for gmax, gmax is only needed for ewald and spme
       SELECT CASE ( ewald_env % ewald_type )
       CASE ( do_ewald_ewald, do_ewald_spme)
         CALL section_vals_val_get(ewald_section, "GMAX",i_vals=gmax_read,error=error)
         SELECT CASE (SIZE(gmax_read,1))
         CASE (1)
           ewald_env%gmax = gmax_read(1)
         CASE (3)
           ewald_env%gmax = gmax_read
         CASE DEFAULT
           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
         END SELECT  
       CASE DEFAULT
         ! this should not be used for do_ewald_none and do_ewald_ewald
         ewald_env%gmax = HUGE(0)
       END SELECT

       CALL section_vals_val_get(ewald_section,&
            "NS_MAX",i_val=ewald_env%ns_max,error=error)
       CALL section_vals_val_get(ewald_section,&
            "O_SPLINE",i_val=ewald_env%o_spline,error=error)
       CALL section_vals_val_get(ewald_section,&
            "EPSILON",r_val=ewald_env%epsilon,error=error)
       CALL section_vals_val_get(ewald_section,&
            "RS_GRID",i_val=ewald_env%rs_distributed_id,error=error)
    END IF

    logger=>cp_error_get_logger(error)

    iw=cp_print_key_unit_nr(logger,ewald_section,"PRINT%PROGRAM_RUN_INFO",&
                                           extension=".log",error=error)

    IF ( iw > 0 ) THEN
       CALL create_ewald_section(section,error=error)
       keyword => section_get_keyword(section,"EWALD_TYPE",error=error)
       CALL keyword_get(keyword,enum=enum,error=error)
       IF ( ewald_env % ewald_type /= do_ewald_none ) THEN
          dummy = cp_unit_from_cp2k(ewald_env % alpha,"angstrom^-1",error=error)
          WRITE ( iw, '( A,T67,A14 )' ) ' Ewald| Summation is done by:', &
               ADJUSTR(TRIM(enum_i2c(enum,ewald_env % ewald_type,error=error)))
          WRITE ( iw, '( A,A18,A,T71,F10.4 )' ) &
               ' Ewald| Alpha parameter [','ANGSTROM^-1',']', dummy

          SELECT CASE ( ewald_env % ewald_type )
          CASE ( do_ewald_ewald)
             WRITE ( iw, '( A,T51,3I10 )' ) &
                  ' Ewald| G-space max. Miller index', ewald_env % gmax
          CASE ( do_ewald_pme)
             WRITE ( iw, '( A,T71,I10 )' ) &
                  ' EWALD| Max small-grid points (input) ', ewald_env % ns_max
             WRITE ( iw, '( A,T71,E10.4 )' ) &
                  ' EWALD| Gaussian tolerance (input) ', ewald_env % epsilon
          CASE ( do_ewald_spme )
             WRITE ( iw, '( A,T51,3I10 )' ) &
                  ' Ewald| G-space max. Miller index', ewald_env % gmax
             WRITE ( iw, '( A,T71,I10 )' ) &
                  ' Ewald| Spline interpolation order ', ewald_env % o_spline
          CASE DEFAULT
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       ELSE
          WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'
       END IF
       CALL section_release(section,error=error)
    END IF

    CALL cp_print_key_finished_output(iw,logger,ewald_section,&
                 "PRINT%PROGRAM_RUN_INFO", error=error)

END SUBROUTINE read_ewald_section

!***************************************************************************

END MODULE ewald_environment_types

!***************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/fist_nonbond_env_types [1.0] *
!!
!!   NAME
!!     fist_nonbond_env_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_nonbond_env_types
  USE global_types,                    ONLY: global_environment_type
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
  USE kinds,                           ONLY: dbl, &
                                             default_string_length
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                             deallocate_neighbor_list_set
  USE splines,                         ONLY: spline_environment_type,&
                                             spline_env_retain, &
                                             spline_env_release
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='fist_nonbond_env_types'
  INTEGER, PRIVATE, SAVE :: last_fist_nonbond_env_id_nr=0
  PUBLIC :: fist_nonbond_env_type, fist_nonbond_env_set, &
            fist_nonbond_env_get, fist_nonbond_env_create, &
            fist_nonbond_env_retain, fist_nonbond_env_release


  TYPE fist_nonbond_env_type
    INTEGER :: ref_count, id_nr, natom_types
    CHARACTER ( len = default_string_length ) :: unit_type
    TYPE ( neighbor_list_set_p_type ), DIMENSION ( : ), POINTER :: nonbonded
    TYPE ( spline_environment_type ), POINTER :: spline_env
    INTEGER :: counter, last_update, num_update
    REAL ( dbl ), DIMENSION ( :, : ), POINTER :: rlist_cutsq
    REAL ( dbl ) :: verlet_skin
    REAL ( dbl ) :: lup, aup
    REAL ( dbl ), DIMENSION ( :, : ), POINTER :: r_last_update
    REAL ( dbl ) :: long_range_correction
  END TYPE fist_nonbond_env_type

!!***
!******************************************************************************
  CONTAINS
!******************************************************************************

!***************************************************************************
!
!****f* fist_nonbond_env_types/fist_nonbond_env_set [1.0] *
!
!   NAME
!     fist_nonbond_env_set
!
!   SYNOPSIS
!     Subroutine fist_nonbond_env_set (fist_nonbond_env )
!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_nonbond_env_create
!
!   FUNCTION
!     sets a fist_nonbond_env
!
!   ARGUMENTS
!     - fist_nonbond_env: the object to create
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
SUBROUTINE fist_nonbond_env_get (fist_nonbond_env, spline_env,  &
                                 rlist_cutsq, verlet_skin, nonbonded, &
                                 aup, lup, r_last_update, num_update, &
                                 last_update, counter, natom_types, &
                                 long_range_correction, error )
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( neighbor_list_set_p_type ), POINTER, &
       DIMENSION ( : ), OPTIONAL               :: nonbonded
  TYPE ( spline_environment_type ), OPTIONAL, POINTER :: spline_env
  REAL ( dbl ), DIMENSION ( :, : ), POINTER, OPTIONAL  :: rlist_cutsq
  REAL ( dbl ), OPTIONAL                      :: verlet_skin
  REAL ( dbl ), POINTER, OPTIONAL             :: r_last_update ( :, : )
  REAL ( dbl ), OPTIONAL                      :: long_range_correction
  INTEGER, OPTIONAL                           :: lup, aup, last_update, &
                                                 num_update, natom_types, &
                                                 counter
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  IF ( PRESENT ( spline_env ) ) spline_env => fist_nonbond_env % spline_env
  IF ( PRESENT ( rlist_cutsq ) ) rlist_cutsq => fist_nonbond_env % rlist_cutsq
  IF ( PRESENT ( verlet_skin ) ) verlet_skin = fist_nonbond_env % verlet_skin
  IF ( PRESENT ( nonbonded ) ) nonbonded => fist_nonbond_env % nonbonded
  IF ( PRESENT ( r_last_update ) ) &
      r_last_update => fist_nonbond_env % r_last_update
  IF ( PRESENT ( lup ) ) lup = fist_nonbond_env % lup 
  IF ( PRESENT ( aup ) ) aup = fist_nonbond_env % aup
  IF ( PRESENT ( natom_types ) ) natom_types = fist_nonbond_env % natom_types
  IF ( PRESENT ( counter ) ) counter = fist_nonbond_env % counter
  IF ( PRESENT ( last_update ) ) last_update = fist_nonbond_env % last_update
  IF ( PRESENT ( num_update ) ) num_update = fist_nonbond_env % num_update
  IF ( PRESENT ( long_range_correction ) ) &
      long_range_correction = fist_nonbond_env % long_range_correction 
END SUBROUTINE fist_nonbond_env_get 
!***************************************************************************
!
!****f* fist_nonbond_env_types/fist_nonbond_env_set [1.0] *
!
!   NAME
!     fist_nonbond_env_set
!
!   SYNOPSIS
!     Subroutine fist_nonbond_env_set (fist_nonbond_env )
!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_nonbond_env_create
!
!   FUNCTION
!     sets a fist_nonbond_env
!
!   ARGUMENTS
!     - fist_nonbond_env: the object to create
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
SUBROUTINE fist_nonbond_env_set (fist_nonbond_env, spline_env,  &
                                 rlist_cutsq, verlet_skin, nonbonded, &
                                 aup, lup, r_last_update, num_update, &
                                 last_update, counter, natom_types, &
                                  long_range_correction, error )
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( neighbor_list_set_p_type ), POINTER, &
       DIMENSION ( : ), OPTIONAL               :: nonbonded
  TYPE ( spline_environment_type ), OPTIONAL, POINTER :: spline_env
  REAL ( dbl ), DIMENSION ( :, : ), POINTER, OPTIONAL  :: rlist_cutsq
  REAL ( dbl ), OPTIONAL                      :: verlet_skin
  REAL ( dbl ), DIMENSION ( :, : ), POINTER, &
                OPTIONAL                      :: r_last_update
  REAL ( dbl ), OPTIONAL                      :: long_range_correction
  INTEGER, OPTIONAL                           :: lup, aup, last_update, &
                                                 num_update, natom_types, &
                                                 counter
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  IF ( PRESENT ( spline_env ) ) THEN
    CALL spline_env_retain ( spline_env, error )
    CALL spline_env_release( fist_nonbond_env % spline_env, error )
    fist_nonbond_env % spline_env => spline_env
  ENDIF
  IF ( PRESENT ( rlist_cutsq ) ) fist_nonbond_env % rlist_cutsq => rlist_cutsq
  IF ( PRESENT ( verlet_skin ) ) fist_nonbond_env % verlet_skin = verlet_skin
  IF ( PRESENT ( nonbonded ) ) fist_nonbond_env % nonbonded => nonbonded
  IF ( PRESENT ( r_last_update ) ) &
      fist_nonbond_env % r_last_update => r_last_update
  IF ( PRESENT ( lup ) ) fist_nonbond_env % lup = lup
  IF ( PRESENT ( aup ) ) fist_nonbond_env % aup = aup
  IF ( PRESENT ( natom_types ) ) fist_nonbond_env % natom_types = natom_types
  IF ( PRESENT ( counter ) ) fist_nonbond_env % counter = counter
  IF ( PRESENT ( last_update ) ) fist_nonbond_env % last_update = last_update
  IF ( PRESENT ( num_update ) ) fist_nonbond_env % num_update = num_update
  IF ( PRESENT ( long_range_correction ) ) &
      fist_nonbond_env % long_range_correction = long_range_correction
END SUBROUTINE fist_nonbond_env_set 
!***************************************************************************
!
!****f* fist_nonbond_env_types/fist_nonbond_env_create [1.0] *
!
!   NAME
!     fist_nonbond_env_create
!
!   SYNOPSIS
!     Subroutine fist_nonbond_env_creat (fist_nonbond_env )
!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_nonbond_env_create
!
!   FUNCTION
!     allocates and intitializes a fist_nonbond_env
!
!   ARGUMENTS
!     - fist_nonbond_env: the object to create
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
SUBROUTINE fist_nonbond_env_create(fist_nonbond_env, spline_env, natom_types, error)
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( spline_environment_type ), POINTER :: spline_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  INTEGER, INTENT ( IN ) :: natom_types 
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_nonbond_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(fist_nonbond_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_fist_nonbond_env( fist_nonbond_env, spline_env,  &
                                 natom_types, error=error )
  END IF
END SUBROUTINE fist_nonbond_env_create
! *****************************************************************************

  SUBROUTINE init_fist_nonbond_env ( fist_nonbond_env, spline_env, natom_types, error )

!   Purpose: Initialise the FIST nonbond environment.

! ***************************************************************************

    TYPE(fist_nonbond_env_type), POINTER      :: fist_nonbond_env
    TYPE(spline_environment_type), POINTER            :: spline_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    INTEGER,  INTENT ( IN ) :: natom_types

!   ---------------------------------------------------------------------------
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='init_fist_nonbond_create',&
            routineP=moduleN//':'//routineN
    INTEGER :: stat, idim, jdim
    REAL ( dbl ) :: rcut


    NULLIFY (fist_nonbond_env%nonbonded)
    NULLIFY (fist_nonbond_env%spline_env)
    NULLIFY (fist_nonbond_env%r_last_update)
    NULLIFY (fist_nonbond_env% rlist_cutsq )
    fist_nonbond_env % unit_type = "ANGSTROM"
    fist_nonbond_env % verlet_skin = 1._dbl
    fist_nonbond_env % lup = 0
    fist_nonbond_env % aup = 0
    fist_nonbond_env % counter = 0
    fist_nonbond_env % last_update = 0
    fist_nonbond_env % num_update = 0
    fist_nonbond_env % long_range_correction = 0
    fist_nonbond_env%ref_count=1
    fist_nonbond_env%id_nr=last_fist_nonbond_env_id_nr
    last_fist_nonbond_env_id_nr=last_fist_nonbond_env_id_nr+1

    fist_nonbond_env%natom_types = natom_types
    ALLOCATE ( fist_nonbond_env % rlist_cutsq ( natom_types, natom_types ), stat = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO idim = 1, natom_types
      DO jdim = 1, natom_types
         rcut = spline_env % potparm ( idim, jdim ) % rcutsq + &
                fist_nonbond_env % verlet_skin
         fist_nonbond_env % rlist_cutsq ( idim, jdim ) = rcut 
      END DO
    END DO
    CALL spline_env_retain ( spline_env, error )
    CALL spline_env_release ( fist_nonbond_env % spline_env, error )
    fist_nonbond_env % spline_env => spline_env
  END SUBROUTINE init_fist_nonbond_env
!***************************************************************************

!****f* fist_nonbond_env_types/fist_nonbond_env_retain [1.0] *
!
!   NAME
!     fist_nonbond_env_retain
!
!   SYNOPSIS
!     Subroutine fist_nonbond_env_retain(fist_nonbond_env, error)
!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_nonbond_env_retain
!
!   FUNCTION
!     retains the given fist_nonbond_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - fist_nonbond_env: the object to retain
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
SUBROUTINE fist_nonbond_env_retain(fist_nonbond_env,error)
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_nonbond_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
!!  
  CPPrecondition(ASSOCIATED(fist_nonbond_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
     fist_nonbond_env%ref_count=fist_nonbond_env%ref_count+1
  END IF
END SUBROUTINE fist_nonbond_env_retain
!***************************************************************************

!****f* fist_nonbond_env_types/fist_nonbond_env_release [1.0] *
!
!   NAME
!     fist_nonbond_env_release
!
!   SYNOPSIS
!     Subroutine fist_nonbond_env_release(fist_nonbond_env, error)
!       Type(fist_nonbond_env_type), Pointer:: fist_nonbond_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_nonbond_env_release
!
!   FUNCTION
!     releases the given fist_nonbond_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - fist_nonbond_env: the object to release
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
SUBROUTINE fist_nonbond_env_release(fist_nonbond_env,error)
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_nonbond_env_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: stat, iab
  failure=.FALSE.

  IF (ASSOCIATED(fist_nonbond_env)) THEN
     CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
     fist_nonbond_env%ref_count=fist_nonbond_env%ref_count-1
     IF (fist_nonbond_env%ref_count<1) THEN
        IF (ASSOCIATED(fist_nonbond_env % nonbonded)) THEN
          DO iab=1,SIZE(fist_nonbond_env % nonbonded)
            CALL deallocate_neighbor_list_set (fist_nonbond_env % nonbonded(iab)%neighbor_list_set)
          END DO
          DEALLOCATE (fist_nonbond_env % nonbonded,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        NULLIFY (fist_nonbond_env%nonbonded)
        CALL spline_env_release ( fist_nonbond_env % spline_env, error )
        IF (ASSOCIATED(fist_nonbond_env % r_last_update)) THEN
          DEALLOCATE (fist_nonbond_env % r_last_update,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        NULLIFY (fist_nonbond_env%r_last_update)
        IF (ASSOCIATED(fist_nonbond_env % rlist_cutsq)) THEN
          DEALLOCATE (fist_nonbond_env % rlist_cutsq,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        NULLIFY (fist_nonbond_env% rlist_cutsq )
     END IF
  END IF
  NULLIFY(fist_nonbond_env)
END SUBROUTINE fist_nonbond_env_release
!***************************************************************************
END MODULE fist_nonbond_env_types
!***************************************************************************

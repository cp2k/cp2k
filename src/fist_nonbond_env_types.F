!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
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
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE fist_neighbor_list_types,        ONLY: neighbor_list_set_p_type, &
                                             deallocate_neighbor_list_set
  USE kinds,                           ONLY: dp, &
                                             default_string_length
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_pp_release
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

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
    TYPE ( pair_potential_pp_type ), POINTER :: potparm
    INTEGER :: counter, last_update, num_update
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: rlist_cutsq
    REAL (KIND=dp) :: verlet_skin
    REAL (KIND=dp) :: lup, aup
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: r_last_update
    REAL (KIND=dp) :: long_range_correction
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
SUBROUTINE fist_nonbond_env_get (fist_nonbond_env, potparm, &
                                 nonbonded, rlist_cutsq, verlet_skin,  &
                                 aup, lup, r_last_update, num_update,  &
                                 last_update, counter, natom_types,    &
                                 long_range_correction, error )
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( neighbor_list_set_p_type ), POINTER, &
       DIMENSION ( : ), OPTIONAL               :: nonbonded
  TYPE ( pair_potential_pp_type ), OPTIONAL, POINTER     :: potparm
  REAL (KIND=dp), DIMENSION ( :, : ), POINTER, OPTIONAL  :: rlist_cutsq
  REAL (KIND=dp), OPTIONAL                      :: verlet_skin
  REAL (KIND=dp), POINTER, OPTIONAL             :: r_last_update ( :, : )
  REAL (KIND=dp), OPTIONAL                      :: long_range_correction
  REAL (KIND=dp), OPTIONAL                           :: lup, aup
  INTEGER, OPTIONAL                           :: last_update, &
                                                 num_update, natom_types, &
                                                 counter
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  IF ( PRESENT ( potparm ) ) potparm => fist_nonbond_env % potparm
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
SUBROUTINE fist_nonbond_env_set ( fist_nonbond_env,  potparm, &
                                  rlist_cutsq, verlet_skin, nonbonded, &
                                  aup, lup, r_last_update, num_update, &
                                  last_update, counter, natom_types, &
                                  long_range_correction, error )
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( pair_potential_pp_type ), POINTER, OPTIONAL :: potparm
  TYPE ( neighbor_list_set_p_type ), POINTER, &
       DIMENSION ( : ), OPTIONAL               :: nonbonded
  REAL (KIND=dp), DIMENSION ( :, : ), POINTER, OPTIONAL  :: rlist_cutsq
  REAL (KIND=dp), OPTIONAL                      :: verlet_skin
  REAL (KIND=dp), DIMENSION ( :, : ), POINTER, &
                OPTIONAL                      :: r_last_update
  REAL (KIND=dp), OPTIONAL                      :: long_range_correction
  REAL (KIND=dp), OPTIONAL                           :: lup, aup
  INTEGER, OPTIONAL                           :: last_update, &
                                                 num_update, natom_types, &
                                                 counter
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  IF ( PRESENT ( potparm ) ) THEN
    fist_nonbond_env % potparm => potparm
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
SUBROUTINE fist_nonbond_env_create(fist_nonbond_env, potparm, error)
  TYPE(fist_nonbond_env_type), POINTER :: fist_nonbond_env
  TYPE ( pair_potential_pp_type ), POINTER :: potparm 
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_nonbond_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(fist_nonbond_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_fist_nonbond_env( fist_nonbond_env, potparm, error=error )
  END IF
END SUBROUTINE fist_nonbond_env_create
! *****************************************************************************

  SUBROUTINE init_fist_nonbond_env ( fist_nonbond_env, potparm, error )

!   Purpose: Initialise the FIST nonbond environment.

! ***************************************************************************

    TYPE(fist_nonbond_env_type), POINTER          :: fist_nonbond_env
    TYPE(pair_potential_pp_type), POINTER         :: potparm 
    TYPE(cp_error_type), INTENT(inout), OPTIONAL  :: error

!   ---------------------------------------------------------------------------
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='init_fist_nonbond_create',&
            routineP=moduleN//':'//routineN
    INTEGER :: stat, idim, jdim, natom_types
    REAL (KIND=dp) :: rcut


    NULLIFY (fist_nonbond_env%nonbonded)
    NULLIFY (fist_nonbond_env%r_last_update)
    NULLIFY (fist_nonbond_env% rlist_cutsq )
    fist_nonbond_env % unit_type = "ANGSTROM"
    fist_nonbond_env % verlet_skin = 1.0_dp
    CALL convert_to_cp2k_units ( "ANGSTROM",  &
         length = fist_nonbond_env % verlet_skin )
    fist_nonbond_env % lup = 0
    fist_nonbond_env % aup = 0
    fist_nonbond_env % counter = 0
    fist_nonbond_env % last_update = 0
    fist_nonbond_env % num_update = 0
    fist_nonbond_env % long_range_correction = 0
    fist_nonbond_env%ref_count=1
    fist_nonbond_env%id_nr=last_fist_nonbond_env_id_nr
    last_fist_nonbond_env_id_nr=last_fist_nonbond_env_id_nr+1

    natom_types = SIZE ( potparm%pot, 1 )
    ALLOCATE ( fist_nonbond_env % rlist_cutsq ( natom_types, natom_types ), stat = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO idim = 1, natom_types
      DO jdim = 1, natom_types
         rcut = SQRT ( potparm%pot ( idim, jdim )%pot % rcutsq ) + &
                fist_nonbond_env % verlet_skin
         fist_nonbond_env % rlist_cutsq ( idim, jdim ) = rcut 
      END DO
    END DO
    fist_nonbond_env % potparm => potparm
    fist_nonbond_env % natom_types = natom_types
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
! release potparm
        CALL pair_potential_pp_release ( fist_nonbond_env % potparm, error )
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
        deallocate(fist_nonbond_env, STAT=stat)
        CPPrecondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(fist_nonbond_env)
END SUBROUTINE fist_nonbond_env_release
!***************************************************************************
END MODULE fist_nonbond_env_types
!***************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/dg_types [1.0] *
!!
!!   NAME
!!     dg_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************
MODULE dg_types

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE dg_rho0_types,                        ONLY: dg_rho0_type, &
                                             dg_rho0_create, &
                                             dg_rho0_retain, &
                                             dg_rho0_release

  USE kinds,                           ONLY: dp,dp_size,default_string_length
  IMPLICIT NONE

! Information on the assignment function for Ewald
  TYPE dg_type
    PRIVATE
    INTEGER :: ref_count, id_nr
    INTEGER :: grid_index
    TYPE ( dg_rho0_type ), POINTER :: dg_rho0
  END TYPE dg_type

  TYPE dg_p_type
    TYPE ( dg_type ), POINTER :: dg
  END TYPE dg_p_type

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='dg_types'
  INTEGER, SAVE, PRIVATE :: last_dg_id=0

  PRIVATE
  PUBLIC :: dg_p_type, dg_type, dg_get,  &
            dg_set, dg_retain, dg_release,  &
            dg_create

  CONTAINS

!-----------------------------------------------------------------------------!
 SUBROUTINE dg_get ( dg, id_nr, grid_index, dg_rho0 ) 

    IMPLICIT NONE
!   Purpose: Get the dg_type

!   ***************************************************************************

    TYPE ( dg_type ), POINTER :: dg
    INTEGER, OPTIONAL :: id_nr, grid_index
    TYPE ( dg_rho0_type ), OPTIONAL, POINTER :: dg_rho0
    
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "dg_get",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( id_nr ) ) id_nr = dg % id_nr
    IF ( PRESENT ( grid_index ) ) grid_index = dg % grid_index
    IF ( PRESENT ( dg_rho0 ) ) dg_rho0 => dg % dg_rho0

  END SUBROUTINE dg_get

! *****************************************************************************

  SUBROUTINE dg_create ( dg, error )

!   Purpose: create the dg structure

!   ***************************************************************************

    TYPE( dg_type ), POINTER :: dg
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "dg_create",&
         routineP=moduleN//':'//routine_name
    INTEGER :: istat
    TYPE ( dg_rho0_type ), POINTER :: dg_rho0

!   ---------------------------------------------------------------------------
    ALLOCATE ( dg, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    NULLIFY ( dg_rho0 )
    CALL dg_rho0_create ( dg_rho0, error )
    dg % dg_rho0 => dg_rho0
    last_dg_id=last_dg_id+1
    dg%id_nr=last_dg_id
    dg%ref_count=1

  END SUBROUTINE dg_create
!***************************************************************************
!!****f* dg_types/dg_retain [1.0] *
!!
!!   NAME
!!     dg_retain
!!
!!   SYNOPSIS
!!     Subroutine dg_retain(dg, error)
!!       Type(dg_type), Pointer:: dg
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dg_retain
!!
!!   FUNCTION
!!     retains the given dg_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - dg_type: the dg_type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dg_retain ( dg, error)
  TYPE(dg_type), POINTER :: dg
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dg_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(dg),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(dg%ref_count>0,cp_failure_level,routineP,error)
     dg%ref_count=dg%ref_count+1
  END IF
END SUBROUTINE dg_retain
!***************************************************************************
!!****f* dg_types/dg_release [1.0] *
!!
!!   NAME
!!     dg_release
!!
!!   SYNOPSIS
!!     Subroutine dg_release(dg, error)
!!       Type(dg_type), Pointer:: dg
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dg_release
!!
!!   FUNCTION
!!     releases the given dg_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - dg: the dg_type to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dg_release(dg, error)
  TYPE(dg_type), POINTER :: dg
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dg_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i,istat

  failure=.FALSE.
  
  IF (ASSOCIATED(dg)) THEN
     CPPreconditionNoFail(dg%ref_count>0,cp_failure_level,routineP,error)
     dg%ref_count=dg%ref_count-1
     IF (dg%ref_count==0) THEN
        CALL dg_rho0_release ( dg % dg_rho0, error = error )
        DEALLOCATE (  dg, stat = stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(dg)
END SUBROUTINE dg_release

! *****************************************************************************
  SUBROUTINE dg_set ( dg, dg_rho0, grid_index )

!   Purpose: Set the double grid environment

! ***************************************************************************

    TYPE ( dg_type ), POINTER :: dg
    TYPE ( dg_rho0_type ), POINTER, OPTIONAL :: dg_rho0
    INTEGER, OPTIONAL :: grid_index

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_dg",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( dg_rho0 ) ) THEN
       CALL dg_rho0_retain ( dg_rho0 ) 
       CALL dg_rho0_release ( dg % dg_rho0 )
       dg % dg_rho0 => dg_rho0
    END IF
    IF ( PRESENT ( grid_index ) ) dg % grid_index = grid_index
!---------------------------------------------------------------------------
  END SUBROUTINE dg_set
!-----------------------------------------------------------------------------!
END MODULE dg_types
!-----------------------------------------------------------------------------!
!!*****

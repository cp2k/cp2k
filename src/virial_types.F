!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/virial_types [1.0] *
!!
!!   NAME
!!     virial_types
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

MODULE virial_types

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
  USE kinds,       ONLY : dbl
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: virial_type, virial_set, &
           virial_create, virial_retain, virial_release, zero_virial

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='virial_types'
  INTEGER, PRIVATE, SAVE :: last_virial_id_nr=0

  TYPE virial_type
     INTEGER :: ref_count, id_nr
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_total, pv_kinetic, &
                                         pv_virial, pv_constraint
     LOGICAL :: pv_availability
  END TYPE virial_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE zero_virial ( virial, virial_in_use )

!   Purpose:  zero the virials

!   ***************************************************************************

    TYPE(virial_type), INTENT(OUT)           :: virial
    LOGICAL, INTENT ( IN ), OPTIONAL         :: virial_in_use

    CHARACTER(LEN=*), PARAMETER :: routine_name = "zero_virial",&
                                   routineP = moduleN//':'//routine_name

    virial % pv_total = 0._dbl 
    virial % pv_kinetic = 0._dbl 
    virial % pv_virial = 0._dbl 
    virial % pv_constraint = 0._dbl 
    IF ( PRESENT ( virial_in_use ) ) THEN
       virial % pv_availability = virial_in_use
    ELSE
       virial % pv_availability = .FALSE.
    ENDIF
! *****************************************************************************
  END SUBROUTINE  zero_virial
! *****************************************************************************

!-----------------------------------------------------------------------------!
 SUBROUTINE virial_set ( virial, pv_total, pv_virial, pv_constraint, &
                         pv_kinetic, pv_availability )

    IMPLICIT NONE
!   Purpose: Set the virial_type

!   ***************************************************************************

    TYPE ( virial_type ), POINTER :: virial
    REAL ( dbl ),  DIMENSION ( 3, 3 ), OPTIONAL :: pv_total, pv_virial, pv_constraint, pv_kinetic
    LOGICAL, OPTIONAL :: pv_availability
    
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "virial_set",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( pv_total ) ) virial % pv_total = pv_total
    IF ( PRESENT ( pv_virial ) ) virial % pv_virial = pv_virial
    IF ( PRESENT ( pv_constraint ) ) virial % pv_constraint = pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) virial % pv_kinetic = pv_kinetic
    IF ( PRESENT ( pv_availability ) ) virial % pv_availability = pv_availability

  END SUBROUTINE virial_set
!-----------------------------------------------------------------------------!
 SUBROUTINE virial_get ( virial, pv_total, pv_virial, pv_constraint, &
                         pv_kinetic, pv_availability )

    IMPLICIT NONE
!   Purpose: Get the virial_type

!   ***************************************************************************

    TYPE ( virial_type ), POINTER :: virial
    REAL ( dbl ),  OPTIONAL, DIMENSION ( 3, 3 ) :: pv_total, pv_virial, pv_constraint, pv_kinetic
    LOGICAL, OPTIONAL :: pv_availability
    
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "virial_get",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( pv_total ) ) pv_total = virial % pv_total
    IF ( PRESENT ( pv_virial ) ) pv_virial = virial % pv_virial
    IF ( PRESENT ( pv_constraint ) ) pv_constraint = virial % pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) pv_kinetic = virial % pv_kinetic 
    IF ( PRESENT ( pv_availability ) ) pv_availability = virial % pv_availability 

  END SUBROUTINE virial_get

! *****************************************************************************

  SUBROUTINE virial_create ( virial, virial_in_use, error )

!   Purpose: create the virial structure

!   ***************************************************************************

    TYPE( virial_type ), POINTER :: virial
    LOGICAL, OPTIONAL :: virial_in_use
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "virial_create",&
         routineP=moduleN//':'//routine_name
    INTEGER :: istat

!   ---------------------------------------------------------------------------
    ALLOCATE ( virial, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    IF ( PRESENT ( virial_in_use ) ) THEN
      CALL zero_virial ( virial, virial_in_use = virial_in_use )
    ELSE
      CALL zero_virial ( virial )
    ENDIF
    last_virial_id_nr=last_virial_id_nr+1
    virial%id_nr=last_virial_id_nr
    virial%ref_count=1

  END SUBROUTINE virial_create
!***************************************************************************

!!****f* virial_types/virial_retain [1.0] *
!!
!!   NAME
!!     virial_retain
!!
!!   SYNOPSIS
!!     Subroutine virial_retain(virial, error)
!!       Type(virial_type), Pointer:: virial
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine virial_retain
!!
!!   FUNCTION
!!     retains the given virial_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - virial_type: the virial_type to retain
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
SUBROUTINE virial_retain ( virial, error)
  TYPE(virial_type), POINTER :: virial
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='virial_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(virial),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(virial%ref_count>0,cp_failure_level,routineP,error)
     virial%ref_count=virial%ref_count+1
  END IF
END SUBROUTINE virial_retain
!***************************************************************************
!!****f* virial_types/virial_release [1.0] *
!!
!!   NAME
!!     virial_release
!!
!!   SYNOPSIS
!!     Subroutine virial_release(virial, error)
!!       Type(virial_type), Pointer:: virial
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine virial_release
!!
!!   FUNCTION
!!     releases the given virial_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - virial: the virial_type to release
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
SUBROUTINE virial_release(virial, error)
  TYPE(virial_type), POINTER :: virial
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='virial_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i,istat

  failure=.FALSE.
  
  IF (ASSOCIATED(virial)) THEN
     CPPreconditionNoFail(virial%ref_count>0,cp_failure_level,routineP,error)
     virial%ref_count=virial%ref_count-1
  END IF
  NULLIFY(virial)
END SUBROUTINE virial_release

!******************************************************************************

END MODULE virial_types

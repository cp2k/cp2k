!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
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

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: virial_type, virial_set, virial_get, cp_virial,&
           virial_create, virial_retain, virial_release, zero_virial

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_types'
  INTEGER, PRIVATE, SAVE :: last_virial_id_nr=0

  TYPE virial_type
     INTEGER :: ref_count, id_nr
     REAL (KIND=dp), DIMENSION ( 3, 3 ) :: pv_total, pv_kinetic, &
                                         pv_virial, pv_constraint
     LOGICAL :: pv_availability,pv_numer
  END TYPE virial_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE cp_virial ( virial_in, virial_out )

!   Purpose:  copy virial_in into virial_out

!   ***************************************************************************

    TYPE(virial_type), POINTER               :: virial_in, virial_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_virial', &
      routineP = moduleN//':'//routineN

    virial_out % pv_total = virial_in % pv_total 
    virial_out % pv_kinetic = virial_in % pv_kinetic  
    virial_out % pv_virial = virial_in % pv_virial 
    virial_out % pv_constraint = virial_in % pv_constraint 
    virial_out % pv_availability = virial_in % pv_availability
    virial_out % pv_numer = virial_in % pv_numer
! *****************************************************************************
  END SUBROUTINE  cp_virial
! *****************************************************************************
!******************************************************************************

  SUBROUTINE zero_virial ( virial, virial_in_use, numer_virial )

!   Purpose:  zero the virials

!   ***************************************************************************

    TYPE(virial_type), INTENT(OUT)           :: virial
    LOGICAL, INTENT(IN), OPTIONAL            :: virial_in_use, numer_virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'zero_virial', &
      routineP = moduleN//':'//routineN

    virial % pv_total = 0.0_dp 
    virial % pv_kinetic = 0.0_dp 
    virial % pv_virial = 0.0_dp 
    virial % pv_constraint = 0.0_dp 
    IF ( PRESENT ( virial_in_use ) ) THEN
       virial % pv_availability = virial_in_use
    ELSE
       virial % pv_availability = .FALSE.
    ENDIF
    IF ( PRESENT ( numer_virial ) ) THEN
       virial % pv_numer = numer_virial
    ELSE
       virial % pv_numer = .FALSE.
    ENDIF
! *****************************************************************************
  END SUBROUTINE  zero_virial
! *****************************************************************************

!-----------------------------------------------------------------------------!
 SUBROUTINE virial_set ( virial, pv_total, pv_virial, pv_constraint, &
                         pv_kinetic, pv_availability, pv_numer )

!   Purpose: Set the virial_type

!   ***************************************************************************

    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_virial, &
                                                pv_constraint, pv_kinetic
    LOGICAL, OPTIONAL                        :: pv_availability, pv_numer

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_set', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( pv_total ) ) virial % pv_total = pv_total
    IF ( PRESENT ( pv_virial ) ) virial % pv_virial = pv_virial
    IF ( PRESENT ( pv_constraint ) ) virial % pv_constraint = pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) virial % pv_kinetic = pv_kinetic
    IF ( PRESENT ( pv_availability ) ) virial % pv_availability = pv_availability
    IF ( PRESENT ( pv_numer ) ) virial % pv_numer = pv_numer

  END SUBROUTINE virial_set
!-----------------------------------------------------------------------------!
 SUBROUTINE virial_get ( virial, pv_total, pv_virial, pv_constraint, &
                         pv_kinetic, pv_availability, pv_numer )

!   Purpose: Get the virial_type

!   ***************************************************************************

    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_virial, &
                                                pv_constraint, pv_kinetic
    LOGICAL, OPTIONAL                        :: pv_availability, pv_numer

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_get', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( pv_total ) ) pv_total = virial % pv_total
    IF ( PRESENT ( pv_virial ) ) pv_virial = virial % pv_virial
    IF ( PRESENT ( pv_constraint ) ) pv_constraint = virial % pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) pv_kinetic = virial % pv_kinetic 
    IF ( PRESENT ( pv_availability ) ) pv_availability = virial % pv_availability 
    IF ( PRESENT ( pv_numer ) ) pv_numer = virial % pv_numer 

  END SUBROUTINE virial_get

! *****************************************************************************

  SUBROUTINE virial_create ( virial, virial_in_use, numer_virial, error )

!   Purpose: create the virial structure

!   ***************************************************************************

    TYPE(virial_type), POINTER               :: virial
    LOGICAL, OPTIONAL                        :: virial_in_use, numer_virial
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE ( virial, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    IF (PRESENT(virial_in_use).AND.PRESENT(numer_virial)) THEN
      CALL zero_virial ( virial, virial_in_use = virial_in_use,&
                         numer_virial = numer_virial )
    ELSE IF (PRESENT(virial_in_use)) THEN
      CALL zero_virial ( virial, virial_in_use = virial_in_use )
    ELSE IF (PRESENT(numer_virial)) THEN
      CALL zero_virial ( virial, numer_virial = numer_virial )
    ELSE
      CALL zero_virial ( virial )
    END IF    
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
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'virial_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'virial_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(virial)) THEN
     CPPreconditionNoFail(virial%ref_count>0,cp_failure_level,routineP,error)
     virial%ref_count=virial%ref_count-1
     IF (virial%ref_count.EQ.0) THEN
        DEALLOCATE(virial)
     ENDIF
  END IF
  NULLIFY(virial)
END SUBROUTINE virial_release

!******************************************************************************

END MODULE virial_types

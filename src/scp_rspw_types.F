!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author CJM 10-3-2007
! *****************************************************************************
MODULE scp_rspw_types
  USE f77_blas
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_rspw_types'
  PUBLIC :: scp_rspw_type, scp_rspw_release, scp_rspw_get, &
            scp_rspw_retain

! *****************************************************************************
  TYPE scp_rspw_type
     INTEGER :: ref_count
     TYPE ( pw_env_type ),         POINTER              :: qs_pw_env
  END TYPE scp_rspw_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! *****************************************************************************
SUBROUTINE scp_rspw_get(scp_rspw, qs_pw_env, error)
    TYPE(scp_rspw_type), POINTER             :: scp_rspw
    TYPE(pw_env_type), OPTIONAL, POINTER     :: qs_pw_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_rspw_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(scp_rspw),cp_failure_level,routineP,error,failure)
  CPPrecondition(scp_rspw%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (PRESENT(qs_pw_env)) qs_pw_env => scp_rspw%qs_pw_env
  END IF
END SUBROUTINE scp_rspw_get
! *****************************************************************************
SUBROUTINE scp_rspw_retain(scp_rspw,error)
    TYPE(scp_rspw_type), POINTER             :: scp_rspw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_rspw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(scp_rspw),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scp_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
     scp_rspw%ref_count=scp_rspw%ref_count+1
  END IF
END SUBROUTINE scp_rspw_retain

! *****************************************************************************
!> \brief releases the memory used by the scp_rspw and scp_rspw
! *****************************************************************************
  SUBROUTINE scp_rspw_release(scp_rspw,error)

    TYPE(scp_rspw_type), POINTER             :: scp_rspw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_rspw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

!---------------------------------------------------------------------------

   CALL timeset(routineN,handle)
   failure=.FALSE.
   IF (ASSOCIATED(scp_rspw)) THEN
     CPPrecondition(scp_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
     scp_rspw%ref_count=scp_rspw%ref_count-1
     IF (scp_rspw%ref_count<1) THEN
      IF (ASSOCIATED(scp_rspw%qs_pw_env)) THEN
        CALL pw_env_release (scp_rspw%qs_pw_env, error=error )
      END IF
      DEALLOCATE(scp_rspw, stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END IF
  NULLIFY(scp_rspw)
  CALL timestop(handle)
  END SUBROUTINE scp_rspw_release

END MODULE scp_rspw_types


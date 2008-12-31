!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief rho_types
!> \par History
!>      JGH (22-Feb-03) PW grid options added
!>      gt 16-nov-03 moved initialization in this new module
!>      12-10-04 contains only the various rho: pw and rs stuff
!>      has moved in pwrs types
!> \author gloria,30.09.2002, previously in kg_force
! *****************************************************************************
MODULE kg_rho_types
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_rho_types'
  PUBLIC :: kg_rho_type, kg_rho_release, kg_rho_get, &
            kg_rho_retain

! *****************************************************************************
  TYPE kg_rho_type
     INTEGER :: ref_count, id_nr
     REAL(KIND=dp) :: total_rho_gspace, total_rho_core_rspace, &
                  total_rho_rspace
     TYPE(pw_p_type), POINTER                     :: rho_g, &
                                                     rhop_g, &
                                                     rho_r, &
                                                     rhop_r,&
                                                     rho_core
  END TYPE kg_rho_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! *****************************************************************************
SUBROUTINE kg_rho_get(kg_rho, rho_r, rho_g, rhop_r, rhop_g, &
                      rho_core, error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(pw_p_type), OPTIONAL, POINTER       :: rho_r, rho_g, rhop_r, rhop_g, &
                                                rho_core
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(kg_rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
    IF (PRESENT(rho_r)) rho_r => kg_rho%rho_r
    IF (PRESENT(rho_g)) rho_g => kg_rho%rho_g
    IF (PRESENT(rhop_r)) rhop_r => kg_rho%rhop_r
    IF (PRESENT(rhop_g)) rhop_g => kg_rho%rhop_g
    IF (PRESENT(rho_core)) rho_core => kg_rho%rho_core
  END IF
END SUBROUTINE kg_rho_get
! *****************************************************************************
SUBROUTINE kg_rho_retain(kg_rho,error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(kg_rho),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_rho%ref_count=kg_rho%ref_count+1
  END IF
END SUBROUTINE kg_rho_retain

! *****************************************************************************
!> \brief releases the memory used by the kg_rho and kg_rho
! *****************************************************************************
  SUBROUTINE kg_rho_release(kg_rho,error)

    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

!---------------------------------------------------------------------------

  CALL timeset(routineN,handle)
  failure=.FALSE.
  IF (ASSOCIATED(kg_rho)) THEN
    CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
    kg_rho%ref_count=kg_rho%ref_count-1
    IF (kg_rho%ref_count<1) THEN
      IF (ASSOCIATED(kg_rho%rho_core)) THEN
        CALL pw_release(kg_rho%rho_core%pw,error=error)
        DEALLOCATE(kg_rho%rho_core,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rho_r)) THEN
        CALL pw_release(kg_rho%rho_r%pw,error=error)
        DEALLOCATE(kg_rho%rho_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rho_g)) THEN
        CALL pw_release(kg_rho%rho_g%pw,error=error)
        DEALLOCATE(kg_rho%rho_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rhop_r)) THEN
        CALL pw_release(kg_rho%rhop_r%pw,error=error)
        DEALLOCATE(kg_rho%rhop_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rhop_g)) THEN
        CALL pw_release(kg_rho%rhop_g%pw,error=error)
        DEALLOCATE(kg_rho%rhop_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      DEALLOCATE(kg_rho,stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END IF
  NULLIFY(kg_rho)
  CALL timestop(handle)
  END SUBROUTINE kg_rho_release

END MODULE kg_rho_types

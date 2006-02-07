!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_rho_types [1.0] *
!!
!!   NAME
!!     kg_rho_types
!!
!!   FUNCTION
!!     rho_types    
!!
!!   AUTHOR
!!     gloria,30.09.2002, previously in kg_force 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!     12-10-04 contains only the various rho: pw and rs stuff
!!     has moved in pwrs types
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rho_types
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_rho_types'
  PUBLIC :: kg_rho_type, kg_rho_release, kg_rho_get, &
            kg_rho_retain

  TYPE kg_rho_type
     INTEGER :: ref_count, id_nr
     REAL(KIND=dp) :: total_rho_gspace, total_rho_core_rspace, &
                  total_rho_rspace  
     TYPE(coeff_type), POINTER                    :: rho_g, &
                                                     rhop_g, &
                                                     rho_r, &
                                                     rhop_r,&
                                                     rho_core 
  END TYPE kg_rho_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!***
!******************************************************************************
SUBROUTINE kg_rho_get(kg_rho, rho_r, rho_g, rhop_r, rhop_g, &
                      rho_core, error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(coeff_type), OPTIONAL, POINTER      :: rho_r, rho_g, rhop_r, rhop_g, &
                                                rho_core
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!!**************************************************************************
SUBROUTINE kg_rho_retain(kg_rho,error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************
!!****f* kg_force/release_kg_rho [1.0] *
!!
!!   NAME
!!     release_kg_rho
!!
!!   FUNCTION
!!     releases the memory used by the kg_rho and kg_rho 
!!
!!   MODIFICATION HISTORY
!!
!!*************************************************************************
  SUBROUTINE kg_rho_release(kg_rho,error)               


    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

!---------------------------------------------------------------------------

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.  
  IF (ASSOCIATED(kg_rho)) THEN
    CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
    kg_rho%ref_count=kg_rho%ref_count-1
    IF (kg_rho%ref_count<1) THEN
      IF (ASSOCIATED(kg_rho%rho_core)) THEN
        CALL coeff_deallocate(kg_rho%rho_core)
        DEALLOCATE(kg_rho%rho_core,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rho_r)) THEN
        CALL coeff_deallocate(kg_rho%rho_r)
        DEALLOCATE(kg_rho%rho_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rho_g)) THEN
        CALL coeff_deallocate(kg_rho%rho_g)
        DEALLOCATE(kg_rho%rho_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rhop_r)) THEN
        CALL coeff_deallocate(kg_rho%rhop_r)
        DEALLOCATE(kg_rho%rhop_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(kg_rho%rhop_g)) THEN
        CALL coeff_deallocate(kg_rho%rhop_g)
        DEALLOCATE(kg_rho%rhop_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      DEALLOCATE(kg_rho,stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END IF
  NULLIFY(kg_rho)
  CALL timestop(0.0_dp,handle)
  END SUBROUTINE kg_rho_release

!***************************************************************************
END MODULE kg_rho_types
!******************************************************************************

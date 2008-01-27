!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_rho_methods [1.0] *
!!
!!   NAME
!!     scp_rho_methods
!!
!!   FUNCTION
!!     rho_methods
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE scp_rho_methods
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_release
  USE scp_rho_types,                   ONLY: scp_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_rho_methods'
  INTEGER, PRIVATE, SAVE :: last_rho_id_nr=0

  PUBLIC :: scp_rho_create

!-----------------------------------------------------------------------------!

CONTAINS

!!*************************************************************************
SUBROUTINE scp_rho_create( scp_rho, auxbas_pw_pool, error )
    TYPE(scp_rho_type), POINTER              :: scp_rho
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_rho_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

!-----------------------------------------------------------------------------!

  CALL timeset(routineN,'I',' ',handle)

  failure=.FALSE.

  IF (.NOT.ASSOCIATED(scp_rho)) THEN
    ALLOCATE(scp_rho, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      scp_rho%ref_count=1
      last_rho_id_nr=last_rho_id_nr+1
      scp_rho%id_nr=last_rho_id_nr
      NULLIFY(scp_rho%rhop_r, scp_rho%rhop_g)
    END IF
  END IF

!rhop_r
  IF(.NOT.failure) THEN
    IF(ASSOCIATED(scp_rho%rhop_r)) THEN
      CALL pw_release (scp_rho%rhop_r%pw,error=error)
      DEALLOCATE(scp_rho%rhop_r,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    ALLOCATE(scp_rho%rhop_r,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
    CALL pw_pool_create_pw(auxbas_pw_pool,scp_rho%rhop_r%pw,&
         use_data=REALDATA3D,in_space=REALSPACE,error=error)
  END IF
!rhop_g
  IF(.NOT.failure) THEN
    IF(ASSOCIATED(scp_rho%rhop_g)) THEN
      CALL pw_release(scp_rho%rhop_g%pw,error=error)
      DEALLOCATE(scp_rho%rhop_g,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    ALLOCATE(scp_rho%rhop_g,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
    CALL pw_pool_create_pw(auxbas_pw_pool,scp_rho%rhop_g%pw,&
         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
  END IF
  CALL timestop(0.0_dp,handle)

END SUBROUTINE scp_rho_create
!!*** **********************************************************************

END MODULE scp_rho_methods

!******************************************************************************

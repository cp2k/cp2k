!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods that change a qs calculation to perform an ep calculation
!> \par History
!>      06.2006 created
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE ep_qs_methods
  USE ep_qs_types,                     ONLY: ep_qs_type
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_qs_methods'
  INTEGER, SAVE, PRIVATE :: last_ep_qs_id=0

  PUBLIC :: ep_qs_create

CONTAINS

! *****************************************************************************
!> \brief creates an ep_qs_env
!> \param ep_qs_env the environemt to create
!> \param qs_env the qs env the ep env will live in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE ep_qs_create(ep_qs_env,qs_env,error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CALL timeset(routineN,handle)
  CPPrecondition(.NOT.ASSOCIATED(ep_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     ALLOCATE(ep_qs_env,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     last_ep_qs_id=last_ep_qs_id+1
     ep_qs_env%id_nr=last_ep_qs_id
     ep_qs_env%ref_count=1
     NULLIFY(ep_qs_env%dH_coeffs,ep_qs_env%dS_coeffs)
     ep_qs_env%tot_rho=0._dp
     ep_qs_env%core_energy=0._dp
     ep_qs_env%ep_active=.FALSE.
  END IF
  CALL timestop(handle)
END SUBROUTINE ep_qs_create

END MODULE ep_qs_methods

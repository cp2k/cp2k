!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/ep_qs_methods [1.0] *
!!
!!   NAME
!!     ep_qs_methods
!!
!!   FUNCTION
!!     methods that change a qs calculation to perform an ep calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2006 created
!!
!!   SOURCE
!****************************************************************************
MODULE ep_qs_methods
  USE kinds, only: dp
  USE timings, ONLY: timeset, timestop
  USE ep_qs_types, only: ep_qs_type
  USE qs_environment_types, only: qs_environment_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='ep_qs_methods'
  INTEGER, SAVE, PRIVATE :: last_ep_qs_id=0

  public :: ep_qs_create
!!***
!****************************************************************************

CONTAINS

!!****f* ep_qs_types/ep_qs_create *
!!
!!   NAME
!!     ep_qs_create
!!
!!   FUNCTION
!!     creates an ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the environemt to create
!!     - qs_env: the qs env the ep env will live in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_create(ep_qs_env,qs_env,error)
  TYPE(ep_qs_type), pointer :: ep_qs_env
  TYPE(qs_environment_type), pointer :: qs_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='ep_qs_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: handle, stat

  failure=.false.

  CALL timeset(routineN,handle)
  CPPrecondition(.NOT.ASSOCIATED(ep_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     ALLOCATE(ep_qs_env,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     last_ep_qs_id=last_ep_qs_id+1
     ep_qs_env%id_nr=last_ep_qs_id
     ep_qs_env%ref_count=1
  END IF
  call timestop(handle)
END SUBROUTINE ep_qs_create
!***************************************************************************

END MODULE ep_qs_methods

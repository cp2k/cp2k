!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/replica_types *
!!
!!   NAME
!!     replica_types
!!
!!   FUNCTION
!!     types used to handle many replica of the same system that differ only
!!     in atom positions, and velocity.
!!     This is useful for things like path integrals or nudged elastic band
!!
!!   NOTES
!!     this is a stupid implementation that replicates all the information
!!     about the replicas, if you really want to do a *lot* of replicas on
!!     a lot of processors you should think about distributiong also that 
!!     information
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE replica_types
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
#include "cp_common_uses.h"

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='replica_types'

  PUBLIC :: replica_env_type
  PUBLIC :: rep_env_retain, rep_env_release
!!***
!****************************************************************************

  !!****s* replica_types/replica_env_type *
  !!
  !!   NAME
  !!     replica_env_type
  !!
  !!   FUNCTION
  !!     keeps replicated information about the replicas
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE replica_env_type
     INTEGER :: ref_count, id_nr, f_env_id
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: x,v,f
     INTEGER, DIMENSION(:), pointer :: replica_owner
     TYPE(cp_para_cart_type), pointer :: para_cart
     TYPE(cp_para_env_type), pointer :: para_env
  END TYPE replica_env_type
  !!***
  !****************************************************************************

contains

!!****f* replica_types/rep_env_retain *
!!
!!   NAME
!!     rep_env_retain
!!
!!   FUNCTION
!!     retains the given replica environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rep_env: the replica environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_retain(rep_env,error)
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='rep_env_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  CPPrecondition(ASSOCIATED(rep_env),cp_fm_types,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     rep_env%ref_count=rep_env%ref_count+1
  END IF
END SUBROUTINE rep_env_retain
!***************************************************************************

!!****f* replica_types/rep_env_release *
!!
!!   NAME
!!     rep_env_release
!!
!!   FUNCTION
!!     releases the given replica environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rep_env: the replica environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_release(rep_env,error)
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='rep_env_release',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  IF (associated(rep_env)) THEN
     CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     rep_env%ref_count=rep_env%ref_count-1
     IF (rep_env%ref_count==0) THEN
        IF (rep_env%f_env_id>0) THEN
           CALL destroy_force_env(rep_env%f_env_id,ierr=ierr)
           CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        END IF
        CALL cp_para_env_release(rep_env%para_env_f,error=error)
        CALL cp_cart_release(rep_env%cart,error=error)
        IF (ASSOCIATED(rep_env%x)) THEN
           DEALLOCATE(rep_env%x,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(rep_env%v)) THEN
           DEALLOCATE(rep_env%v,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(rep_env%f)) THEN
           DEALLOCATE(rep_env%f,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(rep_env%replica_owner)) THEN
           DEALLOCATE(rep_env%replica_owner,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(rep_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(rep_env)
END SUBROUTINE rep_env_release
!***************************************************************************

END MODULE replica_types

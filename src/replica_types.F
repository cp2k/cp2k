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
  USE cp_para_env,                     ONLY: cp_cart_release,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE f77_interface,                   ONLY: destroy_force_env
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_types'

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
  !!     - ref_count: reference count
  !!     - id_nr: identity number (unique or each replica_env)
  !!     - nrep: number of replicas
  !!     - nat: number of atoms (in each replica)
  !!     - ndim: number of dimensions (3*nat)
  !!     - f_env_id: id of the force env that will do the calculations for the
  !!       replicas owned by this processor
  !!     - r,v,f: positions, velocities and forces of the replicas.
  !!       the indexing is as follow (idir,iat,irep)
  !!     - replica_owner: which processor column owns the replica irep
  !!     - cart: 2d distribution of the processors for the replicas,
  !!       a column  work together on the same force_env (i.e. changing the 
  !!       row you stay in the same replica), rows have different replicas
  !!     - para_env: the global para env that contains all the replicas
  !!     - para_env_f: parallel environment of the underlying force 
  !!       environment
  !!     - para_env_inter_rep: parallel environment between processors
  !!       on the same row of cart
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE replica_env_type
     INTEGER :: ref_count, id_nr, f_env_id, nrep,nat,ndim
     REAL(kind=dp), DIMENSION(:,:), POINTER :: r,v,f     
     INTEGER, DIMENSION(:), POINTER :: replica_owner
     TYPE(cp_para_cart_type), POINTER :: cart
     TYPE(cp_para_env_type), POINTER :: para_env, para_env_f,para_env_inter_rep
     LOGICAL :: sync_v
  END TYPE replica_env_type
  !!***
  !****************************************************************************

CONTAINS

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
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
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
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(rep_env)) THEN
     CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     rep_env%ref_count=rep_env%ref_count-1
     IF (rep_env%ref_count==0) THEN
        IF (rep_env%f_env_id>0) THEN
           CALL destroy_force_env(rep_env%f_env_id,ierr=ierr)
           CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        END IF
        CALL cp_para_env_release(rep_env%para_env,error=error)
        CALL cp_cart_release(rep_env%cart,error=error)
        CALL cp_para_env_release(rep_env%para_env_f,error=error)
        IF (ASSOCIATED(rep_env%r)) THEN
           DEALLOCATE(rep_env%r,stat=stat)
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

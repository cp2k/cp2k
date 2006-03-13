!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_ks_types [1.0] *
!!
!!   NAME
!!     qs_ks_types
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_types
  USE coefficient_types,               ONLY: coeff_type
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_types'

  PUBLIC :: qs_ks_env_type, qs_ks_env_p_type
  PUBLIC :: qs_ks_release, qs_ks_retain
!!***
!****************************************************************************

!!****s* qs_ks_methods/qs_ks_env_type [1.0] *
!!
!!   NAME
!!     qs_ks_env_type
!!
!!   FUNCTION
!!     calculation environement to calculate the ks matrix,
!!     holds all the needed vars.
!!     assumes that the core hamiltonian and energy are up to date.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - s_mstruct_changed: if the sparsity structure of s in qs_env where this
!!       ks_env lives has changed
!!     - rho_changed: if the rho stored in the qs_env has changed
!!     - forces_up_to_date: if the 
!!     - potential_changed: if some external potentian changed (such as due to the constraint)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2002 created [fawzi]
!!     10.2002 introduced pools, uses updated rho as input,
!!             removed most temporary variables, renamed may vars,
!!             began conversion to LSD [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_ks_env_type
     LOGICAL :: s_mstruct_changed, &
                rho_changed, potential_changed,&
                forces_up_to_date
     INTEGER :: n_evals, &
                id_nr, ref_count

     ! hartree is supposed to contain the hartree potential (for cube ouput)
     ! ugly to keep it always around only for a cube output...
     TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
     TYPE(coeff_type) :: v_hartree_rspace

  END TYPE qs_ks_env_type
!!***
!****************************************************************************

!!****s* qs_ks_methods/qs_ks_env_p_type [1.0] *
!!
!!   NAME
!!     qs_ks_env_p_type
!!
!!   FUNCTION
!!     type to build arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     ks_env: the ks_env pointer
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ks_env_p_type
     TYPE(qs_ks_env_type), POINTER :: ks_env
  END TYPE qs_ks_env_p_type
!!***
  !****************************************************************************
CONTAINS
!!****f* qs_ks_methods/qs_ks_release [1.0] *
!!
!!   NAME
!!     qs_ks_release
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_release(ks_env, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_release
!!
!!   FUNCTION
!!     releases the ks_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks_env to be released
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_release(ks_env,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ks_env)) THEN
       CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_env%ref_count=ks_env%ref_count-1

       IF (ks_env%ref_count<1) THEN
          CALL pw_pool_give_back_coeff(ks_env%auxbas_pw_pool,ks_env%v_hartree_rspace,&
               error=error)
          CALL pw_pool_release(ks_env%auxbas_pw_pool,error=error)
          DEALLOCATE(ks_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(ks_env)
  END SUBROUTINE qs_ks_release
!***************************************************************************

!!****f* qs_ks_methods/qs_ks_retain [1.0] *
!!
!!   NAME
!!     qs_ks_retain
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_retain(ks_env, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_retain
!!
!!   FUNCTION
!!     retains the given ks_environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the KohnSham environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_retain(ks_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_env%ref_count=ks_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_retain
!***************************************************************************

END MODULE qs_ks_types

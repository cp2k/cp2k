!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_ks_qmmm_types [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_types
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_qmmm_types
  USE coefficient_types,               ONLY: coeff_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cube_utils,                      ONLY: cube_info_type
  USE cube_utils,                      ONLY: destroy_cube_info
  USE pw_env_types, ONLY: pw_env_type, pw_env_get, pw_env_release

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ks_qmmm_types'

  PUBLIC :: qs_ks_qmmm_env_type, qs_ks_qmmm_env_p_type
  PUBLIC :: qs_ks_qmmm_release, qs_ks_qmmm_retain
!!***
!****************************************************************************

!!****s* qs_ks_qmmm_methods/qs_ks_qmmm_env_type [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_env_type
!!
!!   FUNCTION
!!     calculation environement to calculate the ks_qmmm matrix,
!!     holds the QM/MM potential and all the needed variables to
!!     compute the QM/MM electrostatic 1-electron ks matrix
!!     assumes that the core hamiltonian and energy are up to date.
!!
!!   NOTES
!!     -
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_ks_qmmm_env_type
     INTEGER :: n_evals, &
                id_nr, ref_count

     TYPE(pw_env_type), pointer :: pw_env
     TYPE(coeff_type)               :: v_qmmm_rspace
     TYPE(cube_info_type),DIMENSION(:), POINTER  :: cube_info
  END TYPE qs_ks_qmmm_env_type
!!***
!****************************************************************************

!!****s* qs_ks_qmmm_methods/qs_ks_qmmm_env_p_type [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_env_p_type
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 [tlaino]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ks_qmmm_env_p_type
     TYPE(qs_ks_qmmm_env_type), POINTER :: ks_env
  END TYPE qs_ks_qmmm_env_p_type
!!***
  !****************************************************************************
CONTAINS
!!****f* qs_ks_methods/qs_ks_qmmm_release [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_release
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_release(ks_qmmm_env, error)
!!       Type(qs_ks_qmmm_env_type), Pointer:: ks_env_qmmm
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_qmmm_release
!!
!!   FUNCTION
!!     releases the ks_qmmm_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_qmmm_env: the ks_qmmm_env to be released
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_qmmm_release(ks_qmmm_env,error)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_release', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: stat,i
    LOGICAL                                  :: failure
    TYPE(pw_pool_type), pointer :: pool

    failure=.FALSE.

    IF (ASSOCIATED(ks_qmmm_env)) THEN
       CPPrecondition(ks_qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_qmmm_env%ref_count=ks_qmmm_env%ref_count-1

       IF (ks_qmmm_env%ref_count<1) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=pool,error=error)
          CALL pw_pool_give_back_coeff(pool,ks_qmmm_env%v_qmmm_rspace,error=error)
          CALL pw_env_release(ks_qmmm_env%pw_env,error=error)
          IF (ASSOCIATED(ks_qmmm_env%cube_info))THEN
             DO i=1,SIZE(ks_qmmm_env%cube_info)
                CALL destroy_cube_info(ks_qmmm_env%cube_info(i))
             END DO
          END IF
          DEALLOCATE(ks_qmmm_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(ks_qmmm_env)
  END SUBROUTINE qs_ks_qmmm_release
!***************************************************************************

!!****f* qs_ks_qmmm_methods/qs_ks_qmmm_retain [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_retain
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_qmmm_retain(ks_qmmm_env, error)
!!       Type(qs_ks_qmmm_env_type), Pointer:: ks_qmmm_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_qmmm_retain
!!
!!   FUNCTION
!!     retains the given ks_environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_qmmm_env: the KohnSham QM/MM environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_qmmm_retain(ks_qmmm_env, error)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ks_qmmm_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_qmmm_env%ref_count=ks_qmmm_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_qmmm_retain
!***************************************************************************

END MODULE qs_ks_qmmm_types

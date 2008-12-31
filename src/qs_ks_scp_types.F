!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      03.2007 
!> \author CJM
! *****************************************************************************
MODULE qs_ks_scp_types
  USE f77_blas
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_scp_types'

  PUBLIC :: qs_ks_scp_env_type, qs_ks_scp_env_p_type
  PUBLIC :: qs_ks_scp_release, qs_ks_scp_retain

! *****************************************************************************
!> \brief calculation environement to calculate the ks_scp_matrix,
!>      holds the SCP potential and all the needed variables to
!>      compute the SCP electrostatic 1-electron ks matrix
!>      assumes that the core hamiltonian and energy are up to date.
!> \par History
!>      03.2007 created 
!> \author CJM
! *****************************************************************************
  TYPE qs_ks_scp_env_type
     INTEGER :: n_evals, &
                id_nr, ref_count
     LOGICAL :: rho_changed, potential_changed, full_reset
     TYPE ( pw_pool_type ), POINTER        :: auxbas_pw_pool
     TYPE ( pw_p_type )                   :: v_scp_rspace
  END TYPE qs_ks_scp_env_type

! *****************************************************************************
!> \brief type to build arrays of pointers
!> \param ks_scp_env the ks_scp_env pointer
!> \par History
!>      04.2007 
!> \author CJM
! *****************************************************************************
  TYPE qs_ks_scp_env_p_type
     TYPE(qs_ks_scp_env_type), POINTER :: ks_scp_env
  END TYPE qs_ks_scp_env_p_type
CONTAINS

! *****************************************************************************
!> \brief releases the ks_scp_env (see doc/ReferenceCounting.html)
!> \param ks_scp_env the ks_scp_env to be released
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      04.2007 created 
!> \author CJM
! *****************************************************************************
  SUBROUTINE qs_ks_scp_release(ks_scp_env, error)
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_scp_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ks_scp_env)) THEN
       CPPrecondition(ks_scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_scp_env%ref_count=ks_scp_env%ref_count-1

       IF (ks_scp_env%ref_count<1) THEN
          CALL pw_pool_give_back_pw(ks_scp_env%auxbas_pw_pool,ks_scp_env%v_scp_rspace%pw,&
               error=error)
          CALL pw_pool_release(ks_scp_env%auxbas_pw_pool,error=error)
          DEALLOCATE(ks_scp_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(ks_scp_env)

  END SUBROUTINE qs_ks_scp_release

! *****************************************************************************
!> \brief retains the given ks_scp_environment
!> \param ks_scp_env the KohnSham SCP environment to retain
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      03.2007 created 
!> \author CJM
! *****************************************************************************
SUBROUTINE qs_ks_scp_retain(ks_scp_env, error)
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_scp_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ks_scp_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_scp_env%ref_count=ks_scp_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_scp_retain
END MODULE qs_ks_scp_types

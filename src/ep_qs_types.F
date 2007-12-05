!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/ep_qs_types [1.0] *
!!
!!   NAME
!!     ep_qs_types
!!
!!   FUNCTION
!!     types for the modifications of the qs calculation to perform an ep
!!     calculation
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
MODULE ep_qs_types
  USE kinds,                           ONLY: dp
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_qs_types'

  PUBLIC :: ep_qs_type, ep_qs_p_type
  PUBLIC :: ep_qs_retain, ep_qs_release
!!***
!****************************************************************************

!!****s* ep_qs_types/ep_qs_type *
!!
!!   NAME
!!     ep_qs_type
!!
!!   FUNCTION
!!     environement for changes to the normal qs calculation done during an ep
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number for the ep_qs_type
!!     - ref_count: reference count (see cp2k/doc/ReferenceCounting.html)
!!     - ep_active: if ep is active
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
TYPE ep_qs_type
   INTEGER :: id_nr, ref_count
   LOGICAL :: ep_active
   TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: dH_coeffs, dS_coeffs
   REAL(kind=dp) :: core_energy,tot_rho
END TYPE ep_qs_type
!***************************************************************************

!!****s* ep_qs_types/ep_qs_p_type *
!!
!!   NAME
!!     ep_qs_p_type
!!
!!   FUNCTION
!!     pointer to an ep_qs_type for arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ep_qs_env: the pointer to the ep_qs_type
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
TYPE ep_qs_p_type
   TYPE(ep_qs_type), POINTER :: ep_qs_env
END TYPE ep_qs_p_type
!***************************************************************************

CONTAINS

!!****f* ep_qs_types/ep_qs_retain *
!!
!!   NAME
!!     ep_qs_retain
!!
!!   FUNCTION
!!     retains the given ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the object to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_retain(ep_qs_env,error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(ep_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(ep_qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_qs_env%ref_count=ep_qs_env%ref_count+1
  END IF
END SUBROUTINE ep_qs_retain
!***************************************************************************

!!****f* ep_qs_types/ep_qs_release *
!!
!!   NAME
!!     ep_qs_release
!!
!!   FUNCTION
!!     releases the given ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the object to release
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_release(ep_qs_env,error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(ep_qs_env)) THEN
     CPPrecondition(ep_qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_qs_env%ref_count=ep_qs_env%ref_count-1
     IF(ep_qs_env%ref_count==0) THEN
        CALL deallocate_matrix_set(ep_qs_env%dH_coeffs,error=error)
        CALL deallocate_matrix_set(ep_qs_env%dS_coeffs,error=error)
        DEALLOCATE(ep_qs_env,stat=stat)
     END IF
  END IF
  NULLIFY(ep_qs_env)
END SUBROUTINE ep_qs_release
!***************************************************************************


!!****f* ep_qs_types/ep_qs_set *
!!
!!   NAME
!!     ep_qs_set
!!
!!   FUNCTION
!!     sets the components of a ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the env to modify
!!     - dH_coeffs_ptr: new coeffs, might be null, aquires ownership of
!!       the pointer
!!     - dS_coeffs: the new coeffs, might be null, aquires ownership of
!!       the pointer
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_set(ep_qs_env, dH_coeffs_ptr, dS_coeffs_ptr, error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: dH_coeffs_ptr, dS_coeffs_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CALL timeset(routineN,handle)
  CPPrecondition(ASSOCIATED(ep_qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ep_qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (PRESENT(dH_coeffs_ptr)) THEN
        IF (ASSOCIATED(ep_qs_env%dH_coeffs)) THEN
           IF (ASSOCIATED(dH_coeffs_ptr)) THEN
              CPAssert(SIZE(dH_coeffs_ptr)==SIZE(dS_coeffs_ptr),cp_failure_level,routineP,error,failure)
              DO i=1,MIN(SIZE(dH_coeffs_ptr),SIZE(dS_coeffs_ptr))
                 CALL cp_assert(dH_coeffs_ptr(i)%matrix%id_nr/=ep_qs_env%dH_coeffs(i)%matrix%id_nr,&
                      cp_failure_level,cp_assertion_failed,routineP,&
CPSourceFileRef,&
                      error,failure)
              END DO
           END IF
           CALL deallocate_matrix_set(ep_qs_env%dH_coeffs,error=error)
        END IF
        ep_qs_env%dH_coeffs => dH_coeffs_ptr
     END IF
     IF (PRESENT(dS_coeffs_ptr)) THEN
        IF (ASSOCIATED(ep_qs_env%dS_coeffs)) THEN
           IF (ASSOCIATED(dS_coeffs_ptr)) THEN
              CPAssert(SIZE(dS_coeffs_ptr)==SIZE(dS_coeffs_ptr),cp_failure_level,routineP,error,failure)
              DO i=1,MIN(SIZE(dS_coeffs_ptr),SIZE(dS_coeffs_ptr))
                 CALL cp_assert(dS_coeffs_ptr(i)%matrix%id_nr/=ep_qs_env%dS_coeffs(i)%matrix%id_nr,&
                      cp_failure_level,cp_assertion_failed,routineP,&
CPSourceFileRef,&
                      error,failure)
              END DO
           END IF
           CALL deallocate_matrix_set(ep_qs_env%dS_coeffs,error=error)
        END IF
        ep_qs_env%dS_coeffs => dS_coeffs_ptr
     END IF
  END IF
  CALL timestop(handle)
END SUBROUTINE ep_qs_set
!***************************************************************************

END MODULE ep_qs_types
!****************************************************************************

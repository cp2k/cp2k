!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/kg_scf_types [1.0] *
!!
!!   NAME
!!     kg_scf_types
!!
!!   FUNCTION
!!     module that contains the definitions of the scf types
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!****************************************************************************
MODULE kg_scf_types
  USE kg_diis_types,                   ONLY: kg_diis_b_release,&
                                             kg_diis_buffer_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_p_type
  USE qs_ot,                           ONLY: qs_ot_destroy,&
                                             qs_ot_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_scf_types'

  PUBLIC :: kg_scf_env_type, kg_scf_env_p_type,&
            kg_scf_env_retain, kg_scf_env_release

!****************************************************************************

!!****s* kg_scf_types/kg_scf_env_type [1.0] *
!!
!!   NAME
!!     kg_scf_env_type
!!
!!   FUNCTION
!!     wrapper for temporary and cached objects used in the scf iteration
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!***************************************************************************

TYPE kg_scf_env_type
   INTEGER :: id_nr
   INTEGER :: ref_count
   INTEGER :: iter_count
   INTEGER :: print_count
   LOGICAL :: print_iter_line
   INTEGER :: method
   INTEGER :: mixing_method
   REAL(KIND = dp) :: iter_param
   REAL(KIND = dp) :: iter_delta
   REAL(KIND = dp) :: p_mix
   CHARACTER (len=15) :: iter_method
   TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: p_mix_new
   TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: ot_preconditioner
   TYPE(qs_ot_type), POINTER, DIMENSION(:)  :: qs_ot_env
   TYPE(kg_diis_buffer_type), POINTER :: scf_diis_buffer
END TYPE kg_scf_env_type

!****************************************************************************
!!****s* kg_scf_types/kg_scf_env_p_type [1.0] *
!!
!!   NAME
!!     kg_scf_env_p_type
!!
!!   FUNCTION
!!     type to have arrays of pointers
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!***************************************************************************

  TYPE kg_scf_env_p_type
     TYPE(kg_scf_env_type), POINTER :: scf_env
  END TYPE kg_scf_env_p_type

!****************************************************************************

CONTAINS

!!****f* kg_scf_types/kg_scf_env_retain [1.0] *
!!
!!   NAME
!!     kg_scf_env_retain
!!
!!   FUNCTION
!!     retains an scf_env (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE kg_scf_env_retain(scf_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scf_env%ref_count=scf_env%ref_count+1
  END IF
END SUBROUTINE kg_scf_env_retain
!***************************************************************************

!!****f* kg_scf_types/kg_scf_env_release [1.0] *
!!
!!   NAME
!!     kg_scf_env_release
!!
!!   FUNCTION
!!     releases an scf_env (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE kg_scf_env_release(scf_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(scf_env)) THEN
     CPPreconditionNoFail(scf_env%ref_count>0,cp_failure_level,routineP,error)
     scf_env%ref_count=scf_env%ref_count-1
     IF (scf_env%ref_count==0) THEN
        IF (ASSOCIATED(scf_env%p_mix_new)) THEN
           CPPreconditionNoFail(.TRUE.,cp_failure_level,routineP,error) ! we should not end up here, and give back using the pools
           CALL deallocate_matrix_set(scf_env%p_mix_new)
        ENDIF
        IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
           DO i=1,SIZE(scf_env%ot_preconditioner)
                CALL destroy_preconditioner(scf_env%ot_preconditioner(i)%preconditioner)
                DEALLOCATE(scf_env%ot_preconditioner(i)%preconditioner)
           ENDDO
           DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
           DO i=1,SIZE(scf_env%qs_ot_env)
              CALL qs_ot_destroy(scf_env%qs_ot_env(i))
           END DO
           DEALLOCATE(scf_env%qs_ot_env, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL kg_diis_b_release(scf_env%scf_diis_buffer,error=error)
        DEALLOCATE(scf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(scf_env)
END SUBROUTINE kg_scf_env_release
!***************************************************************************

END MODULE kg_scf_types

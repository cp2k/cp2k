!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_scf_types [1.0] *
!!
!!   NAME
!!     qs_scf_types
!!
!!   FUNCTION
!!     module that contains the definitions of the scf types
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_scf_types
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_p_type
  USE qs_diis_types,                   ONLY: qs_diis_b_release,&
                                             qs_diis_buffer_type
  USE qs_ot,                           ONLY: qs_ot_destroy,&
                                             qs_ot_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_types'

  INTEGER, PARAMETER, PUBLIC :: &
                                general_diag_method_nr=1,&
                                special_diag_method_nr=2,&
                                ot_method_nr=3
  INTEGER, PARAMETER, PUBLIC :: no_mixing_nr=0, ao_mixing_nr=1

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

  PUBLIC :: qs_scf_env_type,qs_outer_scf_type
  PUBLIC :: scf_env_retain, scf_env_release, scf_env_create, scf_env_did_change

!!***
!****************************************************************************

!!****s* qs_scf_types/qs_scf_env_type [1.0] *
!!
!!   NAME
!!     qs_scf_env_type
!!
!!   FUNCTION
!!     wrapper for temporary and cached objects used in the scf iteration
!!
!!   NOTES
!!
!!   ATTRIBUTES
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
TYPE qs_outer_scf_type
   INTEGER :: iter_count

   ! these are the variable of outer loop.
   ! right now, we assume that they can be easily written as 
   ! small arrays, but we might want to go the cp_fm_types
   ! at a later stage
   ! also, we just store the full iteration history
   REAL(KIND=dp), DIMENSION(:), POINTER :: energy
   REAL(KIND=dp), DIMENSION(:,:), POINTER :: variables
   REAL(KIND=dp), DIMENSION(:,:), POINTER :: gradient
   
END TYPE qs_outer_scf_type

TYPE qs_scf_env_type
   TYPE(qs_outer_scf_type) :: outer_scf
   INTEGER :: ref_count, id_nr, print_count, iter_count
   INTEGER :: nelectron, method, mixing_method 
   REAL(KIND = dp) :: iter_param, iter_delta, p_mix
   CHARACTER (len=15) :: iter_method
   LOGICAL :: print_iter_line
   TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: scf_work1
   TYPE(cp_fm_type), POINTER  :: scf_work2, ortho, scf_work3
   ! p_mix_new should only be needed for mixing methods (ao_mixing_nr)
   TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: p_mix_new 
   TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: ot_preconditioner
   TYPE(qs_ot_type), POINTER, DIMENSION(:)  :: qs_ot_env
   TYPE(qs_diis_buffer_type), POINTER :: scf_diis_buffer
END TYPE qs_scf_env_type

CONTAINS

!!****f* qs_scf_types/scf_env_create [1.0] *
!!
!!   NAME
!!     scf_env_create
!!
!!   SYNOPSIS
!!     Subroutine scf_env_create(scf_env, qs_env, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_create
!!
!!   FUNCTION
!!     allocates and initialize an scf_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env to initialize
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_create(scf_env,  error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT.failure) THEN

       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
       scf_env%nelectron=0
       scf_env%iter_param=0.0_dp
       scf_env%iter_delta=0.0_dp
       scf_env%iter_method=""
       scf_env%print_iter_line=.TRUE.

       scf_env%method= -1
       scf_env%mixing_method = -1
       scf_env%p_mix=0.0_dp

       scf_env%outer_scf%iter_count=0
       NULLIFY(scf_env%outer_scf%variables, &
            scf_env%outer_scf%gradient,scf_env%outer_scf%energy)

       NULLIFY(scf_env%scf_work1, scf_env%scf_work2, scf_env%ortho,&
            scf_env%p_mix_new, scf_env%ot_preconditioner, scf_env%qs_ot_env,&
            scf_env%scf_diis_buffer)
       NULLIFY (scf_env%scf_work3)

    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_create

!!****f* qs_scf_types/scf_env_retain [1.0] *
!!
!!   NAME
!!     scf_env_retain
!!
!!   SYNOPSIS
!!     Subroutine scf_env_retain(scf_env, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_retain
!!
!!   FUNCTION
!!     retains an scf_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE scf_env_retain(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scf_env%ref_count=scf_env%ref_count+1
  END IF
END SUBROUTINE scf_env_retain
!***************************************************************************

!!****f* qs_scf_types/scf_env_did_change [1.0] *
!!
!!   NAME
!!     scf_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine scf_env_did_change(scf_env, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_did_change
!!
!!   FUNCTION
!!     function to be called to inform the scf_env about changes
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env to inform
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE scf_env_did_change(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (ASSOCIATED(scf_env%p_mix_new)) THEN
        CALL deallocate_matrix_set(scf_env%p_mix_new)
     END IF
  END IF
END SUBROUTINE scf_env_did_change

!!****f* qs_scf_types/scf_env_release [1.0] *
!!
!!   NAME
!!     scf_env_release
!!
!!   SYNOPSIS
!!     Subroutine scf_env_release(scf_env, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_release
!!
!!   FUNCTION
!!     releases an scf_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE scf_env_release(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(scf_env)) THEN
     CPPreconditionNoFail(scf_env%ref_count>0,cp_failure_level,routineP,error)
     scf_env%ref_count=scf_env%ref_count-1
     IF (scf_env%ref_count==0) THEN
        CALL cp_fm_vect_dealloc(scf_env%scf_work1,error=error)
        CALL cp_fm_release(scf_env%scf_work2,error=error)
        CALL cp_fm_release(scf_env%scf_work3,error=error)
        CALL cp_fm_release(scf_env%ortho,error=error)
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
        CALL qs_diis_b_release(scf_env%scf_diis_buffer,error=error)

        IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
           DEALLOCATE(scf_env%outer_scf%variables,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
           DEALLOCATE(scf_env%outer_scf%gradient,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
           DEALLOCATE(scf_env%outer_scf%energy,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        DEALLOCATE(scf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(scf_env)

END SUBROUTINE scf_env_release
!***************************************************************************

END MODULE qs_scf_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_type
  USE qs_diis_types,                   ONLY: qs_diis_b_release,&
                                             qs_diis_buffer_type
  USE qs_ot,                           ONLY: qs_ot_destroy,&
                                             qs_ot_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_scf_types'

  INTEGER, PARAMETER, PUBLIC :: pao_mehod_nr=0, diagonalisation_method_nr=1,&
       ot_method_nr=2
  INTEGER, PARAMETER, PUBLIC :: no_mixing_nr=0, ao_mixing_nr=1

  PUBLIC :: qs_scf_env_type, qs_scf_env_p_type
  PUBLIC :: scf_env_retain, scf_env_release

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
!!     -
!!
!!   ATTRIBUTES
!!     - method: 0: pao, 1: diagonalization (default), 2: orbital transforms
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
TYPE qs_scf_env_type
   INTEGER :: ref_count, id_nr, print_count, iter_count
   INTEGER :: nelectron, method, mixing_method
   REAL(KIND = dp) :: iter_param, iter_delta, p_mix
   CHARACTER (len=15) :: iter_method
   LOGICAL :: print_iter_line
   TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: scf_work1
   TYPE(cp_fm_type), POINTER  :: scf_work2, ortho
   TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: p_new
   TYPE(preconditioner_type), POINTER :: ot_preconditioner
   TYPE(qs_ot_type), POINTER, DIMENSION(:)  :: qs_ot_env
   TYPE(qs_diis_buffer_type), POINTER :: scf_diis_buffer
END TYPE qs_scf_env_type
!!***
  !****************************************************************************

!!****s* qs_scf_types/qs_scf_env_p_type [1.0] *
!!
!!   NAME
!!     qs_scf_env_p_type
!!
!!   FUNCTION
!!     type to have arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - scf_env: the pointer to the scf_env
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_scf_env_p_type
     TYPE(qs_scf_env_type), POINTER :: scf_env
  END TYPE qs_scf_env_p_type
!!***
  !****************************************************************************

CONTAINS

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
  TYPE(qs_scf_env_type), POINTER :: scf_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='scf_env_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scf_env%ref_count=scf_env%ref_count+1
  END IF
END SUBROUTINE scf_env_retain
!***************************************************************************

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
  TYPE(qs_scf_env_type), POINTER :: scf_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='scf_env_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.FALSE.
  
  IF (ASSOCIATED(scf_env)) THEN
     CPPreconditionNoFail(scf_env%ref_count>0,cp_failure_level,routineP,error)
     scf_env%ref_count=scf_env%ref_count-1
     IF (scf_env%ref_count==0) THEN
        CALL cp_fm_vect_dealloc(scf_env%scf_work1,error=error)
        CALL cp_fm_release(scf_env%scf_work2,error=error)
        CALL cp_fm_release(scf_env%ortho,error=error)
        CALL deallocate_matrix_set(scf_env%p_new)
        IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
           CALL destroy_preconditioner(scf_env%ot_preconditioner)
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
        DEALLOCATE(scf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(scf_env)
END SUBROUTINE scf_env_release
!***************************************************************************

END MODULE qs_scf_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ks_types'

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
                rho_changed, &
                forces_up_to_date
     INTEGER :: n_evals, &
                id_nr, ref_count

     ! hartree is supposed to contain the hartree potential (for cube ouput)
     ! ugly to keep it always around only for a cube output...
     TYPE(coeff_type) :: v_hartree_rspace
     TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace

     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: v
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
    TYPE(qs_ks_env_type), POINTER :: ks_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_ks_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, ispin
    failure=.FALSE.

    IF (ASSOCIATED(ks_env)) THEN
       CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_env%ref_count=ks_env%ref_count-1

       IF (ks_env%ref_count<1) THEN
          ! try to give back the grids?
          CALL coeff_deallocate(ks_env%v_hartree_rspace)
          IF (ASSOCIATED(ks_env%v_rspace)) THEN
             DO ispin=1,SIZE(ks_env%v_rspace)
                CALL coeff_deallocate(ks_env%v_rspace(ispin))
             END DO
             DEALLOCATE(ks_env%v_rspace,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(ks_env%v)) THEN
             DO ispin=1,SIZE(ks_env%v)
                IF (ASSOCIATED(ks_env%v(ispin)%matrix)) THEN
                   CALL deallocate_matrix(ks_env%v(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(ks_env%v, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
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
  TYPE(qs_ks_env_type), POINTER :: ks_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ks_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_env%ref_count=ks_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_retain
!***************************************************************************

END MODULE qs_ks_types

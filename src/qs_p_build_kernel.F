!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_p_build_kernel [1.0] *
!!
!!   NAME
!!     qs_p_build_kernel
!!
!!   FUNCTION
!!     module that builds the second order perturbation kernel
!!     K_P_P1 = delta_rho|_P delta_rho|_P E drho(P1) drho
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_p_build_kernel

  USE cp_log_handling, ONLY: cp_to_string, cp_logger_type, cp_failure_level,&
       cp_warning_level, cp_note_level
  USE cp_error_handling, ONLY: cp_error_type, cp_assert, cp_unimplemented_error
  USE kinds, ONLY: wp=>dp
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE qs_environment_types, ONLY: qs_environment_type, get_qs_env
  USE cp_block_matrix, ONLY: cp_block_matrix_type, cp_block_matrix_p_type
  USE qs_p_types, ONLY: qs_K_P_P1_env_type
  USE qs_build_ks_matrix, ONLY: qs_ks_env_type, qs_ks_did_change, &
       qs_ks_update_qs_env

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_build_kernel'

  PUBLIC :: kpp1_init, kpp1_dealloc_ref, kpp1_set_p, kpp1_calculate
!***
!****************************************************************************
  
CONTAINS
  
!!****f* qs_p_build_kernel/kpp1_init [1.0] *
!!
!!   NAME
!!     kpp1_init
!!
!!   FUNCTION
!!     initializes a K_P_P1_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************

SUBROUTINE kpp1_init(kpp1_env,qs_env,ks_env,global_env,error)
  TYPE(qs_K_P_P1_env_type), INTENT(inout) :: kpp1_env
  TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
  TYPE(qs_ks_env_type), INTENT(in), TARGET :: ks_env
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kpp1_init',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(routineP, error=error, error_level=cp_warning_level)
  END IF
END SUBROUTINE kpp1_init

!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_dealloc_ref [1.0] *
!!
!!   NAME
!!     kpp1_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory used by a K_P_P1_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE kpp1_dealloc_ref(kpp1_env,error)
  TYPE(qs_K_P_P1_env_type), INTENT(inout) :: kpp1_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kpp1_dealloc_ref',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(routineP,error=error)
  END IF
END SUBROUTINE kpp1_dealloc_ref
!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_set_p [1.0] *
!!
!!   NAME
!!     kpp1_set_p
!!
!!   FUNCTION
!!     sets the value of P, the density where to evaluate the derivatives
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to update
!!     - p: the new value of P
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE kpp1_set_p(kpp1_env,p,global_env,error)
  TYPE(qs_K_P_P1_env_type), INTENT(inout) :: kpp1_env
  TYPE(cp_block_matrix_p_type), DIMENSION (:), INTENT(in) :: p
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kpp1_set_p',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(routineP,error=error)
  END IF
END SUBROUTINE kpp1_set_p
!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_calculate [1.0] *
!!
!!   NAME
!!     kpp1_calculate
!!
!!   FUNCTION
!!     calculates the kernel matrix composted with P1
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement where co calculate the kernel
!!     - P1: the density to compose the second order tensor to
!!       obtain the K_P_P1 matrix
!!     - kpp1_matrix: matrix in which to store the K_P_P1 matrix
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE kpp1_calculate(kpp1_env, p1, kpp1_matrix,global_env, error)
  TYPE(qs_K_P_P1_env_type), INTENT(inout) :: kpp1_env
  TYPE(cp_block_matrix_p_type), DIMENSION (:), INTENT(in) :: p1
  TYPE(cp_block_matrix_p_type), DIMENSION (:), INTENT(inout) :: kpp1_matrix
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kpp1_calculate',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(routineP,error=error)
  END IF
END SUBROUTINE kpp1_calculate
!***************************************************************************

END MODULE qs_p_build_kernel

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
module qs_p_build_kernel
  use cp_log_handling, only: cp_to_string, cp_logger_type, cp_failure_level,&
       cp_warning_level, cp_note_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_unimplemented_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use cp_block_matrix, only: cp_block_matrix_type, cp_block_matrix_p_type
  use qs_p_types, only: qs_k_p_p1_env_type
  use qs_build_ks_matrix, only: qs_ks_env_type, qs_ks_did_change, &
       qs_ks_update_qs_env
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_p_build_kernel'

  public :: qs_K_P_P1_env_type
  public :: kpp1_init, kpp1_dealloc_ref, kpp1_set_p, kpp1_calculate
!***
!****************************************************************************
  
contains
  
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
subroutine kpp1_init(kpp1_env,qs_env,ks_env,global_env,error)
  type(qs_K_P_P1_env_type), intent(inout) :: kpp1_env
  type(qs_environment_type), intent(in), target :: qs_env
  type(qs_ks_env_type), intent(in), target :: ks_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='kpp1_init',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(routineP, error=error, error_level=cp_warning_level)
  end if
end subroutine kpp1_init
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
subroutine kpp1_dealloc_ref(kpp1_env,error)
  type(qs_K_P_P1_env_type), intent(inout) :: kpp1_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='kpp1_dealloc_ref',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(routineP,error=error)
  end if
end subroutine kpp1_dealloc_ref
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
subroutine kpp1_set_p(kpp1_env,p,global_env,error)
  type(qs_K_P_P1_env_type), intent(inout) :: kpp1_env
  type(cp_block_matrix_p_type), dimension (:), intent(in) :: p
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='kpp1_set_p',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(routineP,error=error)
  end if
end subroutine kpp1_set_p
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
subroutine kpp1_calculate(kpp1_env, p1, kpp1_matrix,global_env, error)
  type(qs_K_P_P1_env_type), intent(inout) :: kpp1_env
  type(cp_block_matrix_p_type), dimension (:), intent(in) :: p1
  type(cp_block_matrix_p_type), dimension (:), intent(inout) :: kpp1_matrix
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='kpp1_calculate',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(routineP,error=error)
  end if
end subroutine kpp1_calculate
!***************************************************************************

end module qs_p_build_kernel

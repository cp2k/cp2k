!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_tddfpt_module *
!!
!!   NAME
!!     qs_tddfpt_module
!!
!!   FUNCTION
!!     Performs density functional perturbation theory (tddfpt) calculations.
!!     Uses the self consistent approach. The tddfpt calculation uses the ground
!!     state of the unperturbed system as the initial state. So the ground state
!!     density has either to be provided to the subroutine or the be calculated
!!     before the tddfpt calculation can be started. This is why there are two 
!!     subroutines.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!
!!****

module qs_tddfpt_module

  use cp_error_handling,    only : cp_error_type, &
                                   cp_error_init
  use cp_log_handling,      only : cp_failure_level
  use global_types,         only : global_environment_type
  use kinds,                only : dbl
  use qs_build_ks_matrix,   only : qs_ks_env_type
  use qs_environment_types, only : qs_environment_type
  use qs_p_build_kernel,    only : kpp1_init, &
                                   qs_K_P_P1_env_type
  use qs_p_types,           only : p_env_init, &
                                   qs_p_env_type
  use sparse_matrix_types,  only : copy_sparse_to_full_matrix, &
                                   get_matrix_info, &
                                   real_matrix_type
  implicit none

  private
  public :: qs_tddfpt_calculation

  character(len=*), parameter :: moduleN='qs_tddfpt_module'

contains

! *****************************************************************************

  subroutine qs_tddfpt_calculation(p_env, ks_env, qs_env, globenv)
!!****f* cp2k/qs_tddfpt_module/qs_tddfpt_calculation
!!
!!   NAME
!!     qs_tddfpt_calculation
!!
!!   FUNCTION
!!     Performs the perturbation calculation
!!
!!   ARGUMENTS
!!     - p_env    - inout - the perturbation environment
!!     - qs_env      - in    - qs environment with ground state
!!     - ks_env      - in    - KS matrix environment
!!     - globenv     - in    - global environment
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****

    ! arguments
    type(qs_p_env_type), intent(inout)        :: p_env
    type(qs_ks_env_type), intent(in)          :: ks_env
    type(qs_environment_type), intent(in)     :: qs_env
    type(global_environment_type), intent(in) :: globenv

    ! locals
    character(len=*), parameter               :: routineN = 'qs_tddfpt_calculation'
    character(len=*), parameter               :: routineP = moduleN//':'//routineN

    type(qs_K_P_P1_env_type)                  :: kpp1_env
    type(cp_error_type)                       :: error

    ! abort if this calculation is not wanted
    if (.NOT.qs_env%dft_control%do_tddfpt_calculation) return

    ! write minibanner
    write (unit=globenv%scr, fmt="(T2,A79)")&
        "",&
        "-------------------------------------------------------------------------------",&
        "-                     Time-dependent calculation started                      -",&
        "-------------------------------------------------------------------------------",&
        ""

    call cp_error_init(error, stop_level=cp_failure_level)

    call kpp1_init(kpp1_env, qs_env, ks_env, globenv, error=error)
    call p_env_init(p_env, qs_env, ks_env, kpp1_env, globenv, error=error)

  end subroutine qs_tddfpt_calculation

! *****************************************************************************

end module qs_tddfpt_module

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

  use global_types,         only : global_environment_type
  use kinds,                only : dbl
  use qs_build_ks_matrix,   only : qs_ks_env_type
  use qs_environment_types, only : qs_environment_type
  use qs_p_build_kernel,    only : qs_K_P_P1_env_type
  use qs_p_types,           only : qs_p_env_type
  use qs_tddfpt_utils

  implicit none

  private
  public :: tddfpt_calculation

  character(len=*), parameter :: moduleN='qs_tddfpt_module'

contains

! *****************************************************************************

  subroutine tddfpt_calculation(p_env, ks_env, qs_env, glob_env)
!!****f* cp2k/tddfpt_module/qs_tddfpt_calculation
!!
!!   NAME
!!     tddfpt_calculation
!!
!!   FUNCTION
!!     Performs the perturbation calculation
!!
!!   ARGUMENTS
!!     - p_env       - inout - the perturbation environment
!!     - qs_env      - in    - qs environment with ground state
!!     - ks_env      - in    - KS matrix environment
!!     - glob_env    - in    - global environment
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
    type(global_environment_type), intent(in) :: glob_env

    ! locals
    character(len=*), parameter               :: routineN = 'qs_tddfpt_calculation'
    character(len=*), parameter               :: routineP = moduleN//':'//routineN

    ! abort if this calculation is not wanted
    if (.NOT.qs_env%dft_control%do_tddfpt_calculation) return

    call tddfpt_write_banner(glob_env)

  end subroutine tddfpt_calculation

! *****************************************************************************

end module qs_tddfpt_module

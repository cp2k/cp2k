!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_tddfpt_utils *
!!
!!   NAME
!!     qs_tddfpt_utils
!!
!!   FUNCTION
!!
!!   SYNOPSIS
!!     
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!****

module qs_tddfpt_utils

  use cp_error_handling,    only : cp_error_type, &
                                   cp_error_init
  use cp_log_handling,      only : cp_failure_level
  use global_types,         only : global_environment_type
  use kinds,                only : dbl
  use qs_build_ks_matrix,   only : qs_ks_env_type
  use qs_environment_types, only : qs_environment_type
  use qs_p_build_kernel,    only : kpp1_init, &
                                   qs_K_P_P1_env_type
  use qs_p_types,           only : p_env_init,&
                                   qs_p_env_type

  implicit none

  private
  public :: &
       tddfpt_init, &
       tddfpt_write_banner

contains

! *****************************************************************************

  subroutine tddfpt_init(p_env, qs_env, glob_env)

    type(qs_p_env_type), intent(inout)        :: p_env
    type(qs_environment_type), intent(in)     :: qs_env
    type(global_environment_type), intent(in) :: glob_env

  end subroutine tddfpt_init

! *****************************************************************************

  subroutine tddfpt_write_banner(globenv)
    type(global_environment_type), intent(in) :: globenv
    write (unit=globenv%scr, fmt="(T2,A79)")&
        "",&
        "-------------------------------------------------------------------------------",&
        "-                     Time-dependent calculation started                      -",&
        "-------------------------------------------------------------------------------",&
        ""    
  end subroutine tddfpt_write_banner

! *****************************************************************************

end module qs_tddfpt_utils

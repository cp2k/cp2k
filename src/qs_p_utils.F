!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_p_utils *
!!
!!   NAME
!!     qs_p_utils
!!
!!   FUNCTION
!!     Utility functions for the perturbation calculations.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     22-08-2002, TCH, started development
!!
!!   SEE ALSO
!!     
!!
!!****
module qs_p_utils

  use cp_error_handling,    only : cp_assert, cp_error_type
  use cp_log_handling
  use global_types,         only : global_environment_type
  use kinds,                only : dbl
  use qs_blacs
  use qs_environment_types, only : qs_environment_type, get_qs_env
  use qs_p_types,           only : qs_p_env_type
  use sparse_matrix_types

  implicit none

  character(len=*), private, parameter :: moduleN='qs_p_utils'

  private
  public :: p_op_l1, p_env_eval_l1

contains

! *****************************************************************************

!!****f* cp2k/qs_p_utils/p_op_l1 *
!!
!!   NAME
!!     p_op_l1
!!
!!   FUNCTION
!!     Evaluates Fv - Sv(epsilon) and stores it in res
!!
!!   ARGUMENTS
!!     res      - the result
!!     p_env    - perturbation calculation environment
!!     v        - the matrix to operate on
!!     glob_env - global environment
!!
!!   MODIFICATION HISTORY
!!
!!   SEE ALSO
!!     
!!
!!****
  subroutine p_op_l1(res, p_env, v, glob_env)
    
    type(blacs_matrix_p_type), dimension(:), pointer   :: res, v
    type(qs_p_env_type), intent(inout)                 :: p_env
    type(global_environment_type), intent(in)          :: glob_env

    integer                                            :: i, n_spins, ncol
    type(real_matrix_type), pointer                    :: k_p

    n_spins = size(v)

    do i=1, n_spins
       k_p => p_env%qs_env%k(i)%matrix
       call get_matrix_info(k_p, ncol=ncol)
       call sparse_times_blacs(k_p, v(i)%blacs_matrix, res(i)%blacs_matrix, ncol, glob_env%para_env)       
       call blacs_symm('R', 'U', p_env%n_ao, p_env%n_mo, &
            1._dbl, p_env%m_epsilon(i)%blacs_matrix, v(i)%blacs_matrix, &
            0._dbl, p_env%tmp_aomo_1)
       call sparse_times_blacs(p_env%qs_env%s(i)%matrix, p_env%tmp_aomo_1, &
            res(i)%blacs_matrix, p_env%n_mo,&
            para_env=glob_env%para_env,alpha=1.0_dbl,beta=1.0_dbl)
    end do

  end subroutine p_op_l1

! *****************************************************************************

!!****f* qs_p_types/p_env_eval_l1 [1.0] *
!!
!!   NAME
!!     p_env_eval_l1
!!
!!   FUNCTION
!!     evaluates the first part of the linear functional:
!!       r = Hrho C - S C epsilon = Hrho C + S C m_epsilon
!!
!!   NOTES
!!     uses tmp_aomo_1, Hrho must already be valid
!!
!!   INPUTS
!!     - p_env: the environment where to perform the calculation
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!****
  subroutine p_env_eval_l1(p_env,global_env,error)

    type(qs_p_env_type), intent(inout)                :: p_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout)      :: error

    logical                                           :: failure
    character(len=*), parameter                       :: routineN = 'p_env_eval_l1'
    character(len=*), parameter                       :: routineP = moduleN//':'//routineN
    type(real_matrix_p_type), dimension(:), pointer   :: s

    failure=.false.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
!!!    if (.not. failure) then
!!!       call get_qs_env(p_env%qs_env,s=s)
!!!       call cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
!!!            p_env%c_ortho(1)%blacs_matrix,p_env%r(1)%blacs_matrix)
!!!       call blacs_symm('R','U',p_env%n_ao,p_env%n_mo,1.0_wp,&
!!!            p_env%m_epsilon(1)%blacs_matrix,&
!!!            p_env%c_ortho(1)%blacs_matrix,0.0_wp,p_env%tmp_aomo_1)
!!!       call sparse_times_blacs(s(1)%matrix,p_env%tmp_aomo_1, &
!!!            p_env%r(1)%blacs_matrix,p_env%n_mo,&
!!!            para_env=global_env%para_env,alpha=1.0_wp,beta=1.0_wp)
!!!    end if

    call p_op_l1(p_env%r, p_env, p_env%c_ortho, global_env)

  end subroutine p_env_eval_l1

! *****************************************************************************

end module qs_p_utils

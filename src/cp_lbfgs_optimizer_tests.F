!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_lbfgs_optimizer_tests [1.0] *
!!
!!   NAME
!!     cp_lbfgs_optimizer_tests
!!
!!   FUNCTION
!!     tests the lbfgs optimizer
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** ************************************************************************
module cp_lbfgs_optimizer_tests
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_log, cp_logger_type, cp_note_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_dealloc_ref, cp_error_get_logger, cp_error_get_print_level
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_lbfgs_optimizer_e_rosen, only: cp_lbfgs_opt_e_rosen_type,&
     cp_init, cp_dealloc_ref, cp_next, cp_get
  use cp_e_rosenbrok_function, only: cp_e_rosenbrok_f_type,&
       cp_ros_f_init, cp_bounds_and_x_ros_f_init, cp_dealloc_ref
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_lbfgs_optimizer_tests'

  public :: cp_lbfgs_e_rosen_test

!****************************************************************************
contains

!!****f* cp_lbfgs_optimizer_tests/cp_lbfgs_e_rosen_test [1.0] *
!!
!!   NAME
!!     cp_lbfgs_e_rosen_test
!!
!!   SYNOPSIS
!!     Function cp_lbfgs_e_rosen_test(global_env, error) Result(res)
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_lbfgs_e_rosen_test
!!
!!   FUNCTION
!!     tests the lbfgs optimizer with the extended Rosenbrock function.
!!     Its minimum should be 0.0. Returns true if the test succeded
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
function cp_lbfgs_e_rosen_test(global_env,error) result(res)
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  logical :: res
  
  logical :: failure, finished
  integer :: handle, i, print_every
  integer, parameter :: n=25
  integer, dimension(n) :: kind_of_bound
  real(kind=wp), dimension(n) :: lower_bound, upper_bound, x0
  real(kind=wp) :: f
  character(len=*), parameter :: routineN='cp_lbfgs_e_rosen_test',&
       routineP=moduleN//':'//routineN
  type(cp_lbfgs_opt_e_rosen_type) :: optimizer
  type(cp_e_rosenbrok_f_type) :: obj_f
  failure=.false.
  
  call timeset(routineN//','//moduleN,'I',"",handle)
  if (.not. failure) then
     call cp_ros_f_init(obj_f, n=n, error=error)
     call cp_bounds_and_x_ros_f_init(lower_bound=lower_bound,&
          upper_bound=upper_bound,kind_of_bound=kind_of_bound,&
          x0=x0,rf_env=obj_f,error=error)
     print_every=cp_error_get_print_level(error)
     if (print_every < cp_warning_level) then
        print_every=1
     else if (print_every==cp_warning_level) then
        print_every=0
     else
        print_every=-1
     end if
     call cp_init(optimizer, global_env, obj_f,x0,&
          print_every=print_every, wanted_relative_f_delta=1.0e7_wp,&
          wanted_projected_gradient=1.0e-7_wp,lower_bound=lower_bound,&
          upper_bound=upper_bound, kind_of_bound=kind_of_bound,&
          m=5,error=error)
     do i=1,100
        if (.not.cp_next(optimizer)) exit
     end do
     call cp_get(optimizer, f=f, at_end=finished)
     CPPostcondition(finished,cp_warning_level,routineP,error,failure)
     CPPostcondition(f<2.0e7*epsilon(0.0_wp),cp_warning_level,routineP,error,failure)
     call cp_dealloc_ref(optimizer,error=error)
     call cp_dealloc_ref(obj_f,error=error)
  end if
  res=.not.failure
  call timestop(0.0_wp,handle)
end function cp_lbfgs_e_rosen_test
!***************************************************************************

end module cp_lbfgs_optimizer_tests

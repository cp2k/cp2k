!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_e_rosenbrok_function [1.0] *
!!
!!   NAME
!!     cp_e_rosenbrok_function
!!
!!   FUNCTION
!!     contains a functional that calculates the extended Rosenbrock function
!!
!!   NOTES
!!     inspired by driver1 of the lbfgs optimizer
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_e_rosenbrok_function
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_e_rosenbrok_function'

  public :: cp_e_rosenbrok_f_type

  public :: cp_init, cp_dealloc_ref, cp_eval_at, cp_ros_f_init, &
       cp_ros_f_dealloc_ref, cp_ros_f_eval_at, cp_bounds_and_x_ros_f_init

  interface cp_init
     module procedure cp_ros_f_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_ros_f_dealloc_ref
  end interface
  interface cp_eval_at
     module procedure cp_ros_f_eval_at
  end interface

!!***
!****************************************************************************
  
!!****s* cp_e_rosenbrok_function/cp_e_rosenbrok_f_type [1.0] *
!!
!!   NAME
!!     cp_e_rosenbrok_f_type
!!
!!   FUNCTION
!!     calculates the extended Rosenbrock function
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  type cp_e_rosenbrok_f_type
     logical :: initialized
     integer :: n
  end type cp_e_rosenbrok_f_type
!!***
  !****************************************************************************

contains
  

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_init [1.0] *
!!
!!   NAME
!!     cp_ros_f_init
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_init(rf_env, n, error)
!!       Type(cp_e_rosenbrok_f_type), Intent (INOUT):: rf_env
!!       Integer, Intent (IN):: n
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_init
!!
!!   FUNCTION
!!     initializes the extended Rosenbrock function computation environement
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     rf_env: the environement to be initialized
!!     n: the dimension of the space we are in
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
subroutine cp_ros_f_init(rf_env, n, error)
  type(cp_e_rosenbrok_f_type), intent(inout) :: rf_env
  integer, intent(in) :: n
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle
  character(len=*), parameter :: routineN='cp_ros_f_init',&
	routineP=moduleN//':'//routineN
  failure=.false.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  rf_env%initialized=.true.
  rf_env%n=n
  !call timestop(0.0_wp,handle)
end subroutine cp_ros_f_init
!***************************************************************************

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_ros_f_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_dealloc_ref(rf_env, error)
!!       Type(cp_e_rosenbrok_f_type), Intent (INOUT):: rf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_dealloc_ref
!!
!!   FUNCTION
!!     releases the momory allocated by the rf_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     rf_env: the extended Rosenbrock function environement
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
subroutine cp_ros_f_dealloc_ref(rf_env, error)
  type(cp_e_rosenbrok_f_type), intent(inout) :: rf_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_ros_f_dealloc_ref',&
	routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  rf_env%initialized=.false.
end subroutine cp_ros_f_dealloc_ref
!***************************************************************************

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_eval_at [1.0] *
!!
!!   NAME
!!     cp_ros_f_eval_at
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_eval_at(rf_env, x, f, gradient, error)
!!       Type(cp_e_rosenbrok_f_type), Intent (IN):: rf_env
!!       Real(Kind=wp), Dimension(:), Intent (IN):: x
!!       Real(Kind=wp), Intent (OUT), Optional:: f
!!       Real(Kind=wp), Intent (OUT), Dimension(:), Optional:: gradient
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_eval_at
!!
!!   FUNCTION
!!     evaluete the extended Rosenbrock function at x
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     c_env: the computation environement
!!     x: the position where the function should be evaluated
!!     f: the function value
!!     gradient: the value of its gradient
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
subroutine cp_ros_f_eval_at(rf_env, x, f, gradient, error)
  type(cp_e_rosenbrok_f_type), intent(in) :: rf_env
  real(kind=wp), dimension(:), intent(in) :: x
  real(kind=wp), intent(out), optional :: f
  real(kind=wp), intent(out), dimension(:), optional :: gradient
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle, i
  character(len=*), parameter :: routineN='cp_ros_f_eval_at',&
	routineP=moduleN//':'//routineN
  real(kind=wp) :: t1,t2
  failure=.false.
  
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  CPPrecondition(size(x)==rf_env%n,cp_warning_level,routineP,error,failure)
  !call timeset(routineN//','//moduleN,'I',"",handle)
  if (present(f)) then
     f = 0.25_wp * (x (1) - 1.0_wp) **2  
     do i = 2, rf_env%n  
        f = f + (x (i) - x (i - 1) **2) **2  
     end do
     f = 4.0_wp * f  
  end if
  if (present(gradient)) then
     ! Compute gradient g for the sample problem.
     t1 = x (2) - x (1) **2  
     if (rf_env%n>0) then ! right if n==1 ??
        gradient (1) = 2.0_wp * (x (1) - 1.0_wp) - 16_wp * x (1) * t1  
        do i = 2, rf_env%n - 1  
           t2 = t1  
           t1 = x (i + 1) - x (i) **2  
           gradient(i) = 8.0_wp * t2 - 16_wp * x (i) * t1  
        end do
        gradient (rf_env%n) = 8.0_wp * t1
     end if
  end if
  !call timestop(0.0_wp,handle)
end subroutine cp_ros_f_eval_at
!***************************************************************************

!!****f* cp_e_rosenbrok_function/cp_bounds_and_x_ros_f_init [1.0] *
!!
!!   NAME
!!     cp_bounds_and_x_ros_f_init
!!
!!   SYNOPSIS
!!     Subroutine cp_bounds_and_x_ros_f_init(lower_bound, upper_bound,&
!!         kind_of_bound, x0, rf_env, error)
!!       Real(Kind=wp), Dimension(:), Intent (OUT):: lower_bound,&
!!         upper_bound, x0
!!       Integer, Dimension(:), Intent (OUT):: kind_of_bound
!!       Type(cp_e_rosenbrok_f_type), Intent (IN):: rf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bounds_and_x_ros_f_init
!!
!!   FUNCTION
!!     initializes the lower and upper bound, and x in a way that the minimum
!!     of the extended Rosenbrock function with these bounds is 0.0
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     lower_bound: the lower bounds
!!     upper_bound: the upper bounds
!!     kind_of_bound: 0: no bound, 1: lower bound, 2: bouth bounds, 
!!          3: upper bound
!!     x0: the initial x
!!     rf_env: the extended Rosenbrock function environement
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
subroutine cp_bounds_and_x_ros_f_init(lower_bound, upper_bound, kind_of_bound,&
     x0,rf_env, error)
  real(kind=wp), dimension(:), intent(out) :: lower_bound, upper_bound, x0
  integer, dimension(:), intent(out) :: kind_of_bound
  type(cp_e_rosenbrok_f_type), intent(in) :: rf_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle, i
  character(len=*), parameter :: routineN='cp_bounds_ros_f_init',&
	routineP=moduleN//':'//routineN
  failure=.false.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  CPPrecondition(size(lower_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(size(upper_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(size(kind_of_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(size(x0)==rf_env%n,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     ! odd numbered variables
     do i = 1, rf_env%n, 2  
        kind_of_bound=2
        lower_bound (i) = 1.0_wp  
        upper_bound (i) = 100.0_wp
     end do
     ! even-numbered variables.
     do i = 2, rf_env%n, 2
        kind_of_bound=2
        lower_bound (i) = - 100.0_wp  
        upper_bound (i) = 100.0_wp
     end do
     x0=3.0_wp
  end if
  !call timestop(0.0_wp,handle)
end subroutine cp_bounds_and_x_ros_f_init
!***************************************************************************

end module cp_e_rosenbrok_function

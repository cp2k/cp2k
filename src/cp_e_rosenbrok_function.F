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
MODULE cp_e_rosenbrok_function
  USE cp_error_handling,               ONLY: cp_assert,&
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
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_e_rosenbrok_function'

  PUBLIC :: cp_e_rosenbrok_f_type

  PUBLIC :: cp_init, cp_dealloc_ref, cp_eval_at, cp_ros_f_init, &
       cp_ros_f_dealloc_ref, cp_ros_f_eval_at, cp_bounds_and_x_ros_f_init

  INTERFACE cp_init
     MODULE PROCEDURE cp_ros_f_init
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_ros_f_dealloc_ref
  END INTERFACE
  INTERFACE cp_eval_at
     MODULE PROCEDURE cp_ros_f_eval_at
  END INTERFACE

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
  TYPE cp_e_rosenbrok_f_type
     LOGICAL :: initialized
     INTEGER :: n
  END TYPE cp_e_rosenbrok_f_type
!!***
  !****************************************************************************

CONTAINS
  

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
SUBROUTINE cp_ros_f_init(rf_env, n, error)
  TYPE(cp_e_rosenbrok_f_type), INTENT(inout) :: rf_env
  INTEGER, INTENT(in) :: n
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: handle
  CHARACTER(len=*), PARAMETER :: routineN='cp_ros_f_init',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  rf_env%initialized=.TRUE.
  rf_env%n=n
  !call timestop(0.0_wp,handle)
END SUBROUTINE cp_ros_f_init
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
SUBROUTINE cp_ros_f_dealloc_ref(rf_env, error)
  TYPE(cp_e_rosenbrok_f_type), INTENT(inout) :: rf_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ros_f_dealloc_ref',&
       routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  rf_env%initialized=.FALSE.
END SUBROUTINE cp_ros_f_dealloc_ref
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
SUBROUTINE cp_ros_f_eval_at(rf_env, x, f, gradient, error)
  TYPE(cp_e_rosenbrok_f_type), INTENT(in) :: rf_env
  REAL(kind=wp), DIMENSION(:), INTENT(in) :: x
  REAL(kind=wp), INTENT(out), OPTIONAL :: f
  REAL(kind=wp), INTENT(out), DIMENSION(:), OPTIONAL :: gradient
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: handle, i
  CHARACTER(len=*), PARAMETER :: routineN='cp_ros_f_eval_at',&
       routineP=moduleN//':'//routineN
  REAL(kind=wp) :: t1,t2
  failure=.FALSE.
  
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(x)==rf_env%n,cp_warning_level,routineP,error,failure)
  !call timeset(routineN//','//moduleN,'I',"",handle)
  IF (PRESENT(f)) THEN
     f = 0.25_wp * (x (1) - 1.0_wp) **2  
     DO i = 2, rf_env%n  
        f = f + (x (i) - x (i - 1) **2) **2  
     END DO
     f = 4.0_wp * f  
  END IF
  IF (PRESENT(gradient)) THEN
     ! Compute gradient g for the sample problem.
     t1 = x (2) - x (1) **2  
     IF (rf_env%n>0) THEN ! right if n==1 ??
        gradient (1) = 2.0_wp * (x (1) - 1.0_wp) - 16_wp * x (1) * t1  
        DO i = 2, rf_env%n - 1  
           t2 = t1  
           t1 = x (i + 1) - x (i) **2  
           gradient(i) = 8.0_wp * t2 - 16_wp * x (i) * t1  
        END DO
        gradient (rf_env%n) = 8.0_wp * t1
     END IF
  END IF
  !call timestop(0.0_wp,handle)
END SUBROUTINE cp_ros_f_eval_at
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
SUBROUTINE cp_bounds_and_x_ros_f_init(lower_bound, upper_bound, kind_of_bound,&
     x0,rf_env, error)
  REAL(kind=wp), DIMENSION(:), INTENT(out) :: lower_bound, upper_bound, x0
  INTEGER, DIMENSION(:), INTENT(out) :: kind_of_bound
  TYPE(cp_e_rosenbrok_f_type), INTENT(in) :: rf_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: handle, i
  CHARACTER(len=*), PARAMETER :: routineN='cp_bounds_ros_f_init',&
       routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(rf_env%initialized,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(lower_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(upper_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(kind_of_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(x0)==rf_env%n,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ! odd numbered variables
     DO i = 1, rf_env%n, 2  
        kind_of_bound=2
        lower_bound (i) = 1.0_wp  
        upper_bound (i) = 100.0_wp
     END DO
     ! even-numbered variables.
     DO i = 2, rf_env%n, 2
        kind_of_bound=2
        lower_bound (i) = - 100.0_wp  
        upper_bound (i) = 100.0_wp
     END DO
     x0=3.0_wp
  END IF
  !call timestop(0.0_wp,handle)
END SUBROUTINE cp_bounds_and_x_ros_f_init
!***************************************************************************

END MODULE cp_e_rosenbrok_function

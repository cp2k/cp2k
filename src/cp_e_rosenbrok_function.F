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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_e_rosenbrok_function'
  INTEGER, SAVE, PRIVATE :: last_e_rosenbrock_id=0

  PUBLIC :: cp_e_rosenbrok_f_type

  PUBLIC :: cp_create, cp_release, cp_retain, cp_eval_at, cp_ros_f_create, &
       cp_ros_f_release, cp_ros_f_retain, cp_ros_f_eval_at, &
       cp_bounds_and_x_ros_f_init

  INTERFACE cp_create
     MODULE PROCEDURE cp_ros_f_create
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE cp_ros_f_release
  END INTERFACE
  INTERFACE cp_retain
     MODULE PROCEDURE cp_ros_f_retain
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
     INTEGER :: n, id_nr, ref_count
  END TYPE cp_e_rosenbrok_f_type
!!***
  !****************************************************************************

CONTAINS
  

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_create [1.0] *
!!
!!   NAME
!!     cp_ros_f_create
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_create(rf_env, n, error)
!!       Type(cp_e_rosenbrok_f_type), Pointer:: rf_env
!!       Integer, Intent (IN):: n
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_create
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
SUBROUTINE cp_ros_f_create(rf_env, n, error)
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: rf_env
    INTEGER, INTENT(in)                      :: n
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ros_f_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(rf_env,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     rf_env%ref_count=1
     rf_env%n=n
     last_e_rosenbrock_id=last_e_rosenbrock_id+1
     rf_env%id_nr=last_e_rosenbrock_id
  END IF
END SUBROUTINE cp_ros_f_create
!***************************************************************************

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_retain [1.0] *
!!
!!   NAME
!!     cp_ros_f_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_retain(rf_env, error)
!!       Type(cp_e_rosenbrok_f_type), Pointer:: rf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_retain
!!
!!   FUNCTION
!!     retains the the given Rosenbrok object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rf_env: the Rosenbrok function to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ros_f_retain(rf_env, error)
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: rf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ros_f_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(rf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(rf_env%ref_count>0,cp_failure_level,routineP,error)
     rf_env%ref_count=rf_env%ref_count+1
  END IF
END SUBROUTINE cp_ros_f_retain
!***************************************************************************

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_release [1.0] *
!!
!!   NAME
!!     cp_ros_f_release
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_release(rf_env, error)
!!       Type(cp_e_rosenbrok_f_type), Pointer:: rf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ros_f_release
!!
!!   FUNCTION
!!     releases the given the rf_env
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
SUBROUTINE cp_ros_f_release(rf_env, error)
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: rf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ros_f_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(rf_env)) THEN
     CPPreconditionNoFail(rf_env%ref_count>0,cp_failure_level,routineP,error)
     rf_env%ref_count=rf_env%ref_count-1
     IF (rf_env%ref_count==0) THEN
        DEALLOCATE(rf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(rf_env)
END SUBROUTINE cp_ros_f_release
!***************************************************************************

!!****f* cp_e_rosenbrok_f_type/cp_ros_f_eval_at [1.0] *
!!
!!   NAME
!!     cp_ros_f_eval_at
!!
!!   SYNOPSIS
!!     Subroutine cp_ros_f_eval_at(rf_env, x, f, gradient, master, para_env,&
!!         error)
!!       Type(cp_e_rosenbrok_f_type), Pointer:: rf_env
!!       Real(KIND = dp), Dimension(:), Intent (IN):: x
!!       Real(KIND = dp), Intent (OUT), Optional:: f
!!       Real(KIND = dp), Pointer, Dimension(:), Optional:: gradient
!!       Integer, Intent (IN):: master
!!       Type(cp_para_env_type), Pointer:: para_env
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
SUBROUTINE cp_ros_f_eval_at(rf_env, x, f, gradient,master,para_env, error)
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: rf_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: x
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ros_f_eval_at', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: t1, t2

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(rf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(rf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(x)==rf_env%n,cp_warning_level,routineP,error,failure)
  !call timeset(routineN//','//moduleN,'I',"",handle)
  IF (PRESENT(f)) THEN
     f = 0.25_dp * (x (1) - 1.0_dp) **2  
     DO i = 2, rf_env%n  
        f = f + (x (i) - x (i - 1) **2) **2  
     END DO
     f = 4.0_dp * f  
  END IF
  IF (PRESENT(gradient)) THEN
     ! Compute gradient g for the sample problem.
     IF (master==para_env%mepos) THEN
        CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
        t1 = x (2) - x (1) **2  
        IF (rf_env%n>0) THEN ! right if n==1 ??
           gradient (1) = 2.0_dp * (x (1) - 1.0_dp) - 16_dp * x (1) * t1  
           DO i = 2, rf_env%n - 1  
              t2 = t1  
              t1 = x (i + 1) - x (i) **2  
              gradient(i) = 8.0_dp * t2 - 16_dp * x (i) * t1  
           END DO
           gradient (rf_env%n) = 8.0_dp * t1
        END IF
     END IF
  END IF
  !call timestop(0.0_dp,handle)
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
!!       Real(KIND = dp), Dimension(:), Intent (OUT):: lower_bound,&
!!         upper_bound, x0
!!       Integer, Dimension(:), Intent (OUT):: kind_of_bound
!!       Type(cp_e_rosenbrok_f_type), Pointer:: rf_env
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
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: lower_bound, upper_bound
    INTEGER, DIMENSION(:), INTENT(out)       :: kind_of_bound
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: x0
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: rf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bounds_ros_f_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(ASSOCIATED(rf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(rf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(lower_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(upper_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(kind_of_bound)==rf_env%n,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(x0)==rf_env%n,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ! odd numbered variables
     DO i = 1, rf_env%n, 2  
        kind_of_bound=2
        lower_bound (i) = 1.0_dp  
        upper_bound (i) = 100.0_dp
     END DO
     ! even-numbered variables.
     DO i = 2, rf_env%n, 2
        kind_of_bound=2
        lower_bound (i) = - 100.0_dp  
        upper_bound (i) = 100.0_dp
     END DO
     x0=3.0_dp
  END IF
  !call timestop(0.0_dp,handle)
END SUBROUTINE cp_bounds_and_x_ros_f_init
!***************************************************************************

END MODULE cp_e_rosenbrok_function

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
MODULE cp_lbfgs_optimizer_tests
  USE cp_e_rosenbrok_function,         ONLY: cp_bounds_and_x_ros_f_init,&
                                             cp_e_rosenbrok_f_type,&
                                             cp_ros_f_create,&
                                             cp_ros_f_release
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_get_print_level,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_lbfgs_optimizer_e_rosen,      ONLY: cp_get,&
                                             cp_lbfgs_opt_e_rosen_type,&
                                             cp_next,&
                                             cp_opt_e_rosen_create,&
                                             cp_opt_e_rosen_release
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_lbfgs_optimizer_tests'

  PUBLIC :: cp_lbfgs_e_rosen_test

!****************************************************************************
CONTAINS

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
FUNCTION cp_lbfgs_e_rosen_test(global_env,error) RESULT(res)
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_lbfgs_e_rosen_test', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: n = 25

    INTEGER                                  :: handle, i, print_every
    INTEGER, DIMENSION(n)                    :: kind_of_bound
    LOGICAL                                  :: failure, finished
    REAL(KIND=dp)                            :: f
    REAL(KIND=dp), DIMENSION(n)              :: lower_bound, upper_bound, x0
    TYPE(cp_e_rosenbrok_f_type), POINTER     :: obj_f
    TYPE(cp_lbfgs_opt_e_rosen_type), POINTER :: optimizer

  failure=.FALSE.
  
  CALL timeset(routineN//','//moduleN,'I',"",handle)
  IF (.NOT. failure) THEN
     CALL cp_ros_f_create(obj_f, n=n, error=error)
     CALL cp_bounds_and_x_ros_f_init(lower_bound=lower_bound,&
          upper_bound=upper_bound,kind_of_bound=kind_of_bound,&
          x0=x0,rf_env=obj_f,error=error)
     print_every=cp_error_get_print_level(error)
     IF (print_every < cp_warning_level) THEN
        print_every=1
     ELSE IF (print_every==cp_warning_level) THEN
        print_every=0
     ELSE
        print_every=-1
     END IF
     CALL cp_opt_e_rosen_create(optimizer, global_env%para_env, obj_f,x0,&
          print_every=print_every, wanted_relative_f_delta=1.0e7_dp,&
          wanted_projected_gradient=1.0e-7_dp,lower_bound=lower_bound,&
          upper_bound=upper_bound, kind_of_bound=kind_of_bound,&
          m=5,error=error)
     DO i=1,100
        IF (.not.cp_next(optimizer)) EXIT
     END DO
     CALL cp_get(optimizer, f=f, at_end=finished)
     CPPostcondition(finished,cp_warning_level,routineP,error,failure)
     CPPostcondition(f<2.0e7*EPSILON(0.0_dp),cp_warning_level,routineP,error,failure)
     CALL cp_opt_e_rosen_release(optimizer,error=error)
     CALL cp_ros_f_release(obj_f,error=error)
  END IF
  res=.not.failure
  CALL timestop(0.0_dp,handle)
END FUNCTION cp_lbfgs_e_rosen_test
!***************************************************************************

END MODULE cp_lbfgs_optimizer_tests

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that optimize a functional using the limited memory bfgs
!>      quasi-newton method.
!>      The process set up so that a master runs the real optimizer and the
!>      others help then to calculate the objective function.
!>      The arguments for the objective function are physicaly present in
!>      every processor (nedeed in the actual implementation of pao).
!>      In the future tha arguments themselves could be distributed.
!> \note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify the .template and .instantition instead
!> \par History
!>      09.2003 globenv->para_env, retain/release, better parallel behaviour
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
MODULE cp_lbfgs_optimizer_gopt
  USE cp_lbfgs,                        ONLY: setulb
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_type
  USE gopt_f_methods,                  ONLY: gopt_f_io
  USE gopt_f_types,                    ONLY: gopt_f_release,&
                                             gopt_f_retain,&
                                             gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_bcast
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_lbfgs_optimizer_gopt'
  INTEGER, PRIVATE, SAVE :: last_lbfgs_optimizer_id=0

  ! types
  PUBLIC :: cp_lbfgs_opt_gopt_type

  ! core methods
  PUBLIC :: cp_create, cp_release, cp_retain, cp_get, cp_next

  ! special methos
  PUBLIC :: cp_step, cp_stop

  ! underlying functions
  PUBLIC :: cp_opt_gopt_create, cp_opt_gopt_release,&
       cp_opt_gopt_get, cp_opt_gopt_next,&
       cp_opt_gopt_step, cp_opt_gopt_stop,&
       cp_opt_gopt_retain

  ! initalize the object
  INTERFACE cp_create
     MODULE PROCEDURE cp_opt_gopt_create
  END INTERFACE
  ! releases the given object
  INTERFACE cp_release
     MODULE PROCEDURE cp_opt_gopt_release
  END INTERFACE
  ! retains the given object
  INTERFACE cp_retain
     MODULE PROCEDURE cp_opt_gopt_retain
  END INTERFACE
  ! returns attributes about the object
  INTERFACE cp_get
     MODULE PROCEDURE cp_opt_gopt_get
  END INTERFACE
  ! goes to the next point, returns true if not converged or error
  INTERFACE cp_next
     MODULE PROCEDURE cp_opt_gopt_next
  END INTERFACE
  ! goes to the next point
  INTERFACE cp_step
     MODULE PROCEDURE cp_opt_gopt_step
  END INTERFACE
  ! stops the iteration
  INTERFACE cp_stop
     MODULE PROCEDURE cp_opt_gopt_stop
  END INTERFACE

! *****************************************************************************
!> \brief info for the optimizer (see the description of this module)
!> \param task the actual task of the optimizer (in the master it is up to
!>        date, in case of error also the slaves one get updated.
!> \param csave internal character string used by the lbfgs optimizer,
!>        meaningful only in the master
!> \param lsave logical array used by the lbfgs optimizer, updated only
!>        in the master
!>        On exit with task = 'NEW_X', the following information is
!>        available:
!>           lsave(1) = .true.  the initial x did not satisfy the bounds;
!>           lsave(2) = .true.  the problem contains bounds;
!>           lsave(3) = .true.  each variable has upper and lower bounds.
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique)
!> \param m the dimension of the subspace used to approximate the second
!>        derivative
!> \param print_every every how many iterations output should be written.
!>        if 0 only at end, if print_every<0 never
!> \param master the pid of the master processor
!> \param max_f_per_iter the maximum number of function evaluations per
!>        iteration
!> \param status 0: just initialized, 1: f g calculation,
!>        2: begin new iteration, 3: ended iteration,
!>        4: normal (converged) exit, 5: abnormal (error) exit,
!>        6: daellocated
!> \param n_iter the actual iteration number
!> \param kind_of_bound an array with 0 (no bound), 1 (lower bound),
!>        2 (both bounds), 3 (upper bound), to describe the bounds
!>        of every variable
!> \param i_work_array an integer workarray of dimension 3*n, present only
!>        in the master
!> \param isave is an INTEGER working array of dimension 44.
!>        On exit with task = 'NEW_X', it contains information that
!>        the user may want to access:
!> \param isave (30) = the current iteration number;
!> \param isave (34) = the total number of function and gradient
!>           evaluations;
!> \param isave (36) = the number of function value or gradient
!>           evaluations in the current iteration;
!> \param isave (38) = the number of free variables in the current
!>           iteration;
!> \param isave (39) = the number of active constraints at the current
!>           iteration;
!> \param f the actual best value of the object function
!> \param wanted_relative_f_delta the wanted relative error on f
!>        (to be multiplied by epsilon), 0.0 -> no check
!> \param wanted_projected_gradient the wanted error on the projected
!>        gradient (hessian times the gradient), 0.0 -> no check
!> \param last_f the value of f in the last iteration
!> \param projected_gradient the value of the sup norm of the projected
!>        gradient
!> \param x the actual evaluation point (best one if converged or stopped)
!> \param lower_bound the lower bounds
!> \param upper_bound the upper bounds
!> \param gradient the actual gradient
!> \param dsave info date for lbfgs (master only)
!> \param work_array a work array for lbfgs (master only)
!> \param para_env the parallel environment for this optimizer
!> \param obj_funct the objective function to be optimized
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  TYPE cp_lbfgs_opt_gopt_type
     CHARACTER (len=60) :: task
     CHARACTER (len=60) :: csave
     LOGICAL :: lsave (4)
     INTEGER :: m, print_every, master, max_f_per_iter, status, n_iter
     INTEGER :: ref_count, id_nr
     INTEGER, DIMENSION(:), POINTER :: kind_of_bound, i_work_array, isave
     REAL(kind=dp) :: f, wanted_relative_f_delta, wanted_projected_gradient,&
          last_f, projected_gradient, eold, emin
     REAL(kind=dp), DIMENSION(:), POINTER :: x,lower_bound,upper_bound,&
          gradient,dsave,work_array
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(gopt_f_type), POINTER :: obj_funct
  END TYPE cp_lbfgs_opt_gopt_type

CONTAINS

! *****************************************************************************
!> \brief initializes the optimizer
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      redirects the lbfgs output the the default unit
!> \par History
!>      02.2002 created [fawzi]
!>      09.2003 refactored (retain/release,para_env) [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_opt_gopt_create(optimizer, para_env, obj_funct, x0, m, print_every, &
       wanted_relative_f_delta, wanted_projected_gradient, lower_bound,upper_bound,&
       kind_of_bound, master, max_f_per_iter, error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gopt_f_type), POINTER               :: obj_funct
    REAL(kind=dp), DIMENSION(:), INTENT(in)  :: x0
    INTEGER, INTENT(in), OPTIONAL            :: m, print_every
    REAL(kind=dp), INTENT(in), OPTIONAL      :: wanted_relative_f_delta, &
                                                wanted_projected_gradient
    REAL(kind=dp), DIMENSION(SIZE(x0)), &
      INTENT(in), OPTIONAL                   :: lower_bound, upper_bound
    INTEGER, DIMENSION(SIZE(x0)), &
      INTENT(in), OPTIONAL                   :: kind_of_bound
    INTEGER, INTENT(in), OPTIONAL            :: master, max_f_per_iter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, lenwa, n, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ALLOCATE(optimizer,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(optimizer%kind_of_bound,&
               optimizer%i_work_array,&
               optimizer%isave,&
               optimizer%x,&
               optimizer%lower_bound,&
               optimizer%upper_bound,&
               optimizer%gradient,&
               optimizer%dsave,&
               optimizer%work_array,&
               optimizer%para_env,&
               optimizer%obj_funct)
       optimizer%ref_count=0
       last_lbfgs_optimizer_id=last_lbfgs_optimizer_id+1
       optimizer%id_nr=last_lbfgs_optimizer_id
       n=SIZE(x0)
       optimizer%m=4
       IF (PRESENT(m)) optimizer%m=m
       optimizer%master= para_env%source
       optimizer%para_env => para_env
       CALL cp_para_env_retain(para_env,error=error)
       optimizer%obj_funct => obj_funct
       CALL gopt_f_retain(obj_funct, error=error)
       optimizer%max_f_per_iter=20
       optimizer%print_every=-1
       optimizer%n_iter=0
       optimizer%f=-1.0_dp
       optimizer%last_f=-1.0_dp
       optimizer%projected_gradient=-1.0_dp
       IF(PRESENT(print_every)) optimizer%print_every=print_every
       IF (PRESENT(master)) optimizer%master=master
       IF (optimizer%master==optimizer%para_env%mepos) THEN
          !MK This has to be adapted for a new L-BFGS version possibly
          lenwa = 2*optimizer%m*n + 5*n + 11*optimizer%m*optimizer%m + 8*optimizer%m
          ALLOCATE (optimizer%kind_of_bound(n), optimizer%i_work_array(3*n),&
                    optimizer%isave(44), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (optimizer%x(n), optimizer%lower_bound(n), &
                    optimizer%upper_bound(n), optimizer%gradient(n), &
                    optimizer%dsave(29), optimizer%work_array(lenwa), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          optimizer%x=x0
          optimizer%task='START'
          optimizer%wanted_relative_f_delta=wanted_relative_f_delta
          optimizer%wanted_projected_gradient=wanted_projected_gradient
          optimizer%kind_of_bound=0
          IF (PRESENT(kind_of_bound)) optimizer%kind_of_bound=kind_of_bound
          IF (PRESENT(lower_bound)) optimizer%lower_bound=lower_bound
          IF (PRESENT(upper_bound)) optimizer%upper_bound=upper_bound

          CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
               optimizer%lower_bound, optimizer%upper_bound, &
               optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
               optimizer%wanted_relative_f_delta, &
               optimizer%wanted_projected_gradient, optimizer%work_array, &
               optimizer%i_work_array, optimizer%task, optimizer%print_every,&
               optimizer%csave, optimizer%lsave, optimizer%isave, &
               optimizer%dsave)
       ELSE
          NULLIFY(&
               optimizer%kind_of_bound, optimizer%i_work_array, optimizer%isave,&
               optimizer%lower_bound, optimizer%upper_bound, optimizer%gradient,&
               optimizer%dsave, optimizer%work_array)
          ALLOCATE (optimizer%x(n),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          optimizer%x(:) = 0.0_dp
          ALLOCATE (optimizer%gradient(n),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          optimizer%gradient(:) = 0.0_dp
       END IF
       CALL mp_bcast(optimizer%x,optimizer%master,optimizer%para_env%group)
       optimizer%status = 0
       optimizer%ref_count = 1
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_opt_gopt_create

! *****************************************************************************
!> \brief retains the given optimizer (see doc/ReferenceCounting.html)
!> \param optimizer the optimizer to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_opt_gopt_retain(optimizer,error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(optimizer%ref_count>0,cp_failure_level,routineP,error)
     optimizer%ref_count=optimizer%ref_count+1
  END IF
END SUBROUTINE cp_opt_gopt_retain

! *****************************************************************************
!> \brief releases the optimizer (see doc/ReferenceCounting.html)
!> \param optimizer the object that should be released
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \par History
!>      02.2002 created [fawzi]
!>      09.2003 dealloc_ref->release [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_opt_gopt_release(optimizer,error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

  CALL timeset(routineN,handle)
  failure=.FALSE.

  IF (ASSOCIATED(optimizer)) THEN
     CPPreconditionNoFail(optimizer%ref_count>0,cp_failure_level,routineP,error)
     optimizer%ref_count=optimizer%ref_count-1
     IF (optimizer%ref_count==0) THEN
        optimizer%status=6
        IF (ASSOCIATED(optimizer%kind_of_bound)) THEN
           DEALLOCATE(optimizer%kind_of_bound, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%i_work_array)) THEN
           DEALLOCATE(optimizer%i_work_array, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%isave)) THEN
           DEALLOCATE(optimizer%isave, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%x)) THEN
           DEALLOCATE(optimizer%x, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%lower_bound)) THEN
           DEALLOCATE(optimizer%lower_bound, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%upper_bound)) THEN
           DEALLOCATE(optimizer%upper_bound, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%gradient)) THEN
           DEALLOCATE(optimizer%gradient, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%dsave)) THEN
           DEALLOCATE(optimizer%dsave, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(optimizer%work_array)) THEN
           DEALLOCATE(optimizer%work_array, stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        CALL cp_para_env_release(optimizer%para_env,error=error)
        CALL gopt_f_release(optimizer%obj_funct,error=error)
        NULLIFY(optimizer%obj_funct)
        DEALLOCATE(optimizer, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(optimizer)
  CALL timestop(handle)
END SUBROUTINE cp_opt_gopt_release

! *****************************************************************************
!> \brief takes different valuse from the optimizer
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_opt_gopt_get(optimizer, para_env, &
       obj_funct, m, print_every, id_nr, &
       wanted_relative_f_delta, wanted_projected_gradient,&
       x, lower_bound, upper_bound, kind_of_bound, master,&
       actual_projected_gradient, &
       n_var, n_iter, status, max_f_per_iter,at_end,&
       is_master, last_f, f, error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(gopt_f_type), OPTIONAL, POINTER     :: obj_funct
    INTEGER, INTENT(out), OPTIONAL           :: m, print_every, id_nr
    REAL(kind=dp), INTENT(out), OPTIONAL     :: wanted_relative_f_delta, &
                                                wanted_projected_gradient
    REAL(kind=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: x, lower_bound, upper_bound
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: kind_of_bound
    INTEGER, INTENT(out), OPTIONAL           :: master
    REAL(kind=dp), INTENT(out), OPTIONAL     :: actual_projected_gradient
    INTEGER, INTENT(out), OPTIONAL           :: n_var, n_iter, status, &
                                                max_f_per_iter
    LOGICAL, INTENT(out), OPTIONAL           :: at_end, is_master
    REAL(kind=dp), INTENT(out), OPTIONAL     :: last_f, f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    !  call timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
    CPPrecondition(optimizer%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(is_master)) is_master=optimizer%master==optimizer%para_env%mepos
       IF (PRESENT(master)) master=optimizer%master
       IF (PRESENT(status)) status=optimizer%status
       IF (PRESENT(para_env)) para_env => optimizer%para_env
       IF (PRESENT(id_nr)) id_nr=optimizer%id_nr
       IF (PRESENT(obj_funct)) obj_funct=optimizer%obj_funct
       IF (PRESENT(m)) m=optimizer%m
       IF (PRESENT(max_f_per_iter)) max_f_per_iter=optimizer%max_f_per_iter
       IF (PRESENT(wanted_projected_gradient)) &
            wanted_projected_gradient=optimizer%wanted_projected_gradient
       IF (PRESENT(wanted_relative_f_delta)) &
            wanted_relative_f_delta=optimizer%wanted_relative_f_delta
       IF (PRESENT(print_every)) print_every=optimizer%print_every
       IF (PRESENT(x)) x => optimizer%x
       IF (PRESENT(n_var)) n_var=SIZE(x)
       IF (PRESENT(lower_bound)) lower_bound => optimizer%lower_bound
       IF (PRESENT(upper_bound)) upper_bound => optimizer%upper_bound
       IF (PRESENT(kind_of_bound)) kind_of_bound => optimizer%kind_of_bound
       IF (PRESENT(n_iter)) n_iter= optimizer%n_iter
       IF (PRESENT(last_f)) last_f = optimizer%last_f
       IF (PRESENT(f)) f= optimizer%f
       IF (PRESENT(at_end)) at_end=optimizer%status>3
       IF (PRESENT(actual_projected_gradient))&
            actual_projected_gradient=optimizer%projected_gradient
       IF (optimizer%master==optimizer%para_env%mepos) THEN
          IF (optimizer%isave(30)>1 .AND. (optimizer%task(1:5)=="NEW_X".OR.&
               optimizer%task(1:4)=="STOP".AND.optimizer%task(7:9)=="CPU")) THEN
             ! nr iterations >1 .and. dsave contains the wanted data
             IF (PRESENT(last_f)) last_f = optimizer%dsave(2)
             IF (PRESENT(actual_projected_gradient))&
                  actual_projected_gradient=optimizer%dsave(13)
          ELSE
             CPPrecondition(.NOT.PRESENT(last_f),cp_warning_level,routineP,error,failure)
             CPPrecondition(.NOT.PRESENT(actual_projected_gradient),cp_warning_level,routineP,error,failure)
          END IF
       ELSE
          CALL cp_assert(.NOT.(PRESENT(lower_bound).OR.&
               PRESENT(upper_bound).OR.PRESENT(kind_of_bound)),&
               cp_warning_level,cp_assertion_failed,&
               routineP, "asked undefined types in "//&
CPSourceFileRef,&
               error,failure)
       END IF

    END IF
    !  call timestop(handle)
  END SUBROUTINE cp_opt_gopt_get

! *****************************************************************************
!> \brief does one optimization step
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      use directly mainlb in place of setulb ??
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  RECURSIVE SUBROUTINE cp_opt_gopt_step(optimizer,n_iter,f,last_f,&
       projected_gradient, converged, geo_section, input, force_env,&
       gopt_param, error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    INTEGER, INTENT(out), OPTIONAL           :: n_iter
    REAL(kind=dp), INTENT(out), OPTIONAL     :: f, last_f, projected_gradient
    LOGICAL, INTENT(out), OPTIONAL           :: converged
    TYPE(section_vals_type), POINTER         :: geo_section, input
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_step', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=5)                         :: wildcard
    INTEGER                                  :: dataunit, handle, its, stat
    LOGICAL                                  :: conv, failure, is_master, &
                                                justEntred
    REAL(KIND=dp)                            :: t_diff, t_now, t_old
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xold
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger,xold)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)
    failure=.FALSE.
    justEntred=.TRUE.
    is_master=optimizer%master==optimizer%para_env%mepos
    IF (PRESENT(converged)) converged = optimizer%status==4
    ALLOCATE(xold(SIZE(optimizer%x)),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    xold = optimizer%x
    t_old = m_walltime()

    CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
    CPPrecondition(optimizer%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (optimizer%status>=4) THEN
       CALL cp_assert(.FALSE.,level=cp_warning_level,&
            error_nr=cp_assertion_failed,fromWhere=routineP,&
            message="status>=4, trying to restart in "//&
CPSourceFileRef,&
            error=error)
       optimizer%status=0
       IF (is_master) THEN
          optimizer%task='START'
          CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
               optimizer%lower_bound, optimizer%upper_bound, &
               optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
               optimizer%wanted_relative_f_delta, &
               optimizer%wanted_projected_gradient, optimizer%work_array, &
               optimizer%i_work_array, optimizer%task, optimizer%print_every,&
               optimizer%csave, optimizer%lsave, optimizer%isave, &
               optimizer%dsave)
       END IF
    END IF
    IF (.NOT. failure) THEN
       DO
          ifMaster: IF (is_master) THEN
             IF (optimizer%task(1:7)=='RESTART') THEN
                ! restart the optimizer
                optimizer%status=0
                optimizer%task='START'
                CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
                     optimizer%lower_bound, optimizer%upper_bound, &
                     optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                     optimizer%wanted_relative_f_delta, &
                     optimizer%wanted_projected_gradient, optimizer%work_array, &
                     optimizer%i_work_array, optimizer%task, optimizer%print_every,&
                     optimizer%csave, optimizer%lsave, optimizer%isave, &
                     optimizer%dsave)
             END IF
             IF (optimizer%task(1:2)=='FG') THEN
                IF (optimizer%isave(36)>optimizer%max_f_per_iter) THEN
                   optimizer%task='STOP: CPU, hit max f eval in iter'
                   optimizer%status=5 ! anormal exit
                   CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
                        optimizer%lower_bound, optimizer%upper_bound, &
                        optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                        optimizer%wanted_relative_f_delta, &
                        optimizer%wanted_projected_gradient, optimizer%work_array, &
                        optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                        optimizer%csave, optimizer%lsave, optimizer%isave, &
                        optimizer%dsave)
                ELSE
                   optimizer%status=1
                END IF
             ELSE IF (optimizer%task(1:5)=='NEW_X') THEN
                IF (justEntred) THEN
                   optimizer%status=2
                   CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
                        optimizer%lower_bound, optimizer%upper_bound, &
                        optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                        optimizer%wanted_relative_f_delta, &
                        optimizer%wanted_projected_gradient, optimizer%work_array, &
                        optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                        optimizer%csave, optimizer%lsave, optimizer%isave, &
                        optimizer%dsave)
                ELSE
                   optimizer%status=3
                END IF
             ELSE IF (optimizer%task(1:4)=='CONV') THEN
                optimizer%status=4
             ELSE IF (optimizer%task(1:4)=='STOP') THEN
                optimizer%status=5
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                     routineP,"task became stop in an unknown way in "//&
CPSourceFileRef,&
                     error)
             ELSE IF (optimizer%task(1:5)=='ERROR') THEN
                optimizer%status=5
             ELSE
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                     routineP,"unknown task '"//optimizer%task//"' in "//&
CPSourceFileRef,&
                     error)
             END IF
          END IF ifMaster
          CALL mp_bcast(optimizer%status,optimizer%master, optimizer%para_env%group)
          ! Dump info
          IF (optimizer%status == 3) THEN
             its = 0
             IF (is_master) THEN
                ! Iteration level is taken into account in the optimizer external loop
                its = optimizer%isave(30)
             END IF
          END IF
          !
          SELECT CASE (optimizer%status)
          CASE (1)
             !op=1 evaluate f and g
             CALL cp_eval_at(optimizer%obj_funct, x=optimizer%x,&
                  f=optimizer%f,&
                  gradient=optimizer%gradient,&
                  final_evaluation=.FALSE.,&
                  master=optimizer%master,para_env=optimizer%para_env,&
                  error=error) ! do not use keywords?
             IF (is_master) THEN
                CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
                     optimizer%lower_bound, optimizer%upper_bound, &
                     optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                     optimizer%wanted_relative_f_delta, &
                     optimizer%wanted_projected_gradient, optimizer%work_array, &
                     optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                     optimizer%csave, optimizer%lsave, optimizer%isave, &
                     optimizer%dsave)
             END IF
             CALL mp_bcast(optimizer%x,optimizer%master,optimizer%para_env%group)
          CASE (2)
             !op=2 begin new iter
             CALL mp_bcast(optimizer%x,optimizer%master,optimizer%para_env%group)
             t_old = m_walltime()
          CASE (3)
             !op=3 ended iter
             wildcard = "LBFGS"
             dataunit = cp_print_key_unit_nr(logger,geo_section,&
                  "PRINT%PROGRAM_RUN_INFO",extension=".geoLog",error=error)
             IF (is_master) its = optimizer%isave(30)
             CALL mp_bcast(its,optimizer%master,optimizer%para_env%group)

             ! Some IO and Convergence check
             t_now  = m_walltime()
             t_diff = t_now - t_old
             t_old  = t_now
             CALL gopt_f_io(optimizer%obj_funct, force_env, force_env%root_section, &
                  its, optimizer%f, dataunit, optimizer%eold, optimizer%emin, wildcard, gopt_param,&
                  SIZE(optimizer%x), optimizer%x-xold, optimizer%gradient, conv, used_time=t_diff,&
                  error=error)
             CALL mp_bcast(conv,optimizer%master,optimizer%para_env%group)
             CALL cp_print_key_finished_output(dataunit,logger,geo_section,&
                  "PRINT%PROGRAM_RUN_INFO", error=error)
             optimizer%eold = optimizer%f
             optimizer%emin = MIN(optimizer%emin,optimizer%eold)
             xold = optimizer%x
             IF (PRESENT(converged)) converged=conv
             EXIT
          CASE (4)
             !op=4 (convergence - normal exit)
             ! Specific L-BFGS convergence criteria.. overrides the convergence criteria on
             ! stepsize and gradients
             dataunit = cp_print_key_unit_nr(logger,geo_section,&
                  "PRINT%PROGRAM_RUN_INFO",extension=".geoLog",error=error)
             IF (dataunit>0) THEN
                WRITE(dataunit,'(T2,A)')""
                WRITE(dataunit,'(T2,A)')"***********************************************"
                WRITE(dataunit,'(T2,A)')"* Specific L-BFGS convergence criteria         "
                WRITE(dataunit,'(T2,A)')"* WANTED_PROJ_GRADIENT and WANTED_REL_F_ERROR  "
                WRITE(dataunit,'(T2,A)')"* satisfied .... run CONVERGED!                "
                WRITE(dataunit,'(T2,A)')"***********************************************"
                WRITE(dataunit,'(T2,A)')""
             END IF
             CALL cp_print_key_finished_output(dataunit,logger,geo_section,&
                  "PRINT%PROGRAM_RUN_INFO", error=error)
             IF (PRESENT(converged)) converged=.TRUE.
             EXIT
          CASE (5)
             ! op=5 abnormal exit (error)
             CALL mp_bcast(optimizer%task,optimizer%master,&
                  optimizer%para_env%group)
          CASE (6)
             ! deallocated
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"step on a deallocated opt structure "//&
CPSourceFileRef,&
                  error,failure)
          CASE default
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"unknown status "//cp_to_string(optimizer%status)//&
                  " in "//&
CPSourceFileRef,&
                  error,failure)
             optimizer%status=5
             EXIT
          END SELECT
          IF (optimizer%status == 1 .AND. justEntred) THEN
             optimizer%eold = optimizer%f
             optimizer%emin = optimizer%eold
          END IF
          justEntred=.FALSE.
       END DO
       CALL mp_bcast(optimizer%x,optimizer%master,&
            optimizer%para_env%group)
       CALL cp_opt_gopt_bcast_res(optimizer,&
            n_iter=optimizer%n_iter,&
            f=optimizer%f, last_f=optimizer%last_f,&
            projected_gradient=optimizer%projected_gradient,error=error)
    END IF
    DEALLOCATE(xold,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(f)) f=optimizer%f
    IF (PRESENT(last_f)) last_f=optimizer%last_f
    IF (PRESENT(projected_gradient)) &
         projected_gradient=optimizer%projected_gradient
    IF (PRESENT(n_iter)) n_iter=optimizer%n_iter
    CALL timestop(handle)

  END SUBROUTINE cp_opt_gopt_step

! *****************************************************************************
!> \brief returns the results (and broadcasts them)
!> \param optimizer the optimizer object the info is taken from
!> \param n_iter the number of iterations
!> \param f the actual value of the objective function (f)
!> \param last_f the last value of f
!> \param projected_gradient the infinity norm of the projected gradient
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      private routine
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_opt_gopt_bcast_res(optimizer, n_iter,f,last_f,&
       projected_gradient,error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    INTEGER, INTENT(out), OPTIONAL           :: n_iter
    REAL(kind=dp), INTENT(inout), OPTIONAL   :: f, last_f, projected_gradient
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_bcast_res', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(4)              :: results

    failure=.FALSE.

    !  call timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
    CPPrecondition(optimizer%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (optimizer%master==optimizer%para_env%mepos) THEN
          results=(/ REAL(optimizer%isave(30), kind=dp), &
               optimizer%f, optimizer%dsave(2), optimizer%dsave(13) /)
       END IF
       CALL mp_bcast(results,optimizer%master, &
            optimizer%para_env%group)
       IF (PRESENT(n_iter)) n_iter= NINT(results(1))
       IF (PRESENT(f)) f=results(2)
       IF (PRESENT(last_f)) last_f=results(3)
       IF (PRESENT(projected_gradient)) &
            projected_gradient=results(4)
    END IF
    !  call timestop(handle)
  END SUBROUTINE cp_opt_gopt_bcast_res

! *****************************************************************************
!> \brief goes to the next optimal point (after an optimizer iteration)
!>      returns true if converged
!> \param optimizer the optimizer that goes to the next point
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      if you deactivate convergence control it returns never false
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  RECURSIVE FUNCTION cp_opt_gopt_next(optimizer,n_iter,f,last_f,&
       projected_gradient, converged, geo_section,input, force_env,&
       gopt_param, error) RESULT(res)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    INTEGER, INTENT(out), OPTIONAL           :: n_iter
    REAL(kind=dp), INTENT(out), OPTIONAL     :: f, last_f, projected_gradient
    LOGICAL, INTENT(out)                     :: converged
    TYPE(section_vals_type), POINTER         :: geo_section, input
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_next', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    !call timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
    CPPrecondition(optimizer%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_opt_gopt_step(optimizer,n_iter=n_iter,f=f,&
            last_f=last_f, projected_gradient=projected_gradient,&
            converged=converged,geo_section=geo_section,&
            input=input,force_env=force_env,gopt_param=gopt_param,&
            error=error)
       res= (optimizer%status<40) .AND. .NOT. converged
    ELSE
       res=.FALSE.
    END IF
    !call timestop(handle)
  END FUNCTION cp_opt_gopt_next

! *****************************************************************************
!> \brief stops the optimization
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      necessary???
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_opt_gopt_stop(optimizer, error)
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_stop', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    !  call timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(optimizer),cp_failure_level,routineP,error,failure)
    CPPrecondition(optimizer%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       optimizer%task='STOPPED on user request'
       optimizer%status=4 ! normal exit
       IF (optimizer%master==optimizer%para_env%mepos) THEN
          CALL setulb (SIZE(optimizer%x), optimizer%m, optimizer%x, &
               optimizer%lower_bound, optimizer%upper_bound, &
               optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
               optimizer%wanted_relative_f_delta, &
               optimizer%wanted_projected_gradient, optimizer%work_array, &
               optimizer%i_work_array, optimizer%task, optimizer%print_every, &
               optimizer%csave, optimizer%lsave, optimizer%isave, &
               optimizer%dsave)
       END IF
    END IF
    ! call timestop(handle)
  END SUBROUTINE cp_opt_gopt_stop

  ! template def put here so that line numbers in template and derived
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,USE,type1_retain,type1_release,type1_eval_at)]
! ARGS:
!  INTERFACE = "#include "gopt_f77_methods.h""
!  USE = "USE gopt_f_types, ONLY: gopt_f_type, gopt_f_release,gopt_f_retain"
!  nametype1 = "gopt"
!  type1 = "type(gopt_f_type)"
!  type1_eval_at = "cp_eval_at"
!  type1_release = "gopt_f_release"
!  type1_retain = "gopt_f_retain"


END MODULE cp_lbfgs_optimizer_gopt

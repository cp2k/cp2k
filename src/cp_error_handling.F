!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_error_handling [1.0] *
!!
!!   NAME
!!     cp_error_handling
!!
!!   FUNCTION
!!     Module that contains the routines for error handling
!!     The idea is that in general routines have an optional 
!!     parameter that can be used to change the error behaviour.
!!
!!   NOTES
!!     In general 0 means no error, a negative error code means an internal
!!     error (in the routine, postcondition failure), and a positive error 
!!     an error on the side of the caller (precondition failure).
!!     Error numbers in the dange -1000..1000 are reserved for generic error
!!     codes.
!!     Inspired from the nag nl90 library.
!!     Some help on this is also available in the high level dev notes on pao.
!!
!!     The following error levels are defined:
!!     0 - note    : CP_NOTE_LEVEL
!!     1 - warning : CP_WARNING_LEVEL
!!     2 - failure : CP_FAILURE_LEVEL
!!     3 - fatal   : CP_FATAL_LEVEL
!!     You shouldn't use the numbers but the constants defined in 
!!     cp_log_handling
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_error_handling
  USE cp_log_handling,                 ONLY: CP_FAILURE_LEVEL,&
                                             CP_FATAL_LEVEL,&
                                             CP_NOTE_LEVEL,&
                                             CP_WARNING_LEVEL,&
                                             cp_get_default_logger,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_release,&
                                             cp_logger_retain,&
                                             cp_logger_type,&
                                             cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_alltoall
  USE output_utilities,                ONLY: print_message
  USE timings,                         ONLY: print_stack

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_error_handling'

  ! types
  PUBLIC :: cp_error_type
  ! parameters var
  PUBLIC :: cp_debug, cp_no_error, cp_caller_error, cp_wrong_args_error,&
       cp_precondition_failed, cp_internal_error, cp_postcondition_failed,&
       cp_invariant_failed, cp_assertion_failed, cp_unimplemented_error_nr
  ! public routines
  PUBLIC :: cp_assert, cp_a_l, cp_simple_assert, cp_unimplemented_error

  PUBLIC :: cp_error_init, cp_error_dealloc_ref, cp_error_set,&
       cp_error_get, cp_error_reset, cp_error_get_level,&
       cp_error_get_print_level, cp_error_get_nr,&
       cp_error_get_logger, cp_error_get_stop_level,&
       cp_error_common_stop, cp_error_handle_error,&
       cp_error_message, cp_error_propagate_error,&
       cp_error_check, cp_error_synchronize_error

#ifdef FD_DEBUG
!! global variable to turn on/off the debugging, as i am lazy I use the 
!! preprocessor variable FD_DEBUG to decide if it should be on or off.
!! this way it is always right (on in debug builds, off in optimized ones)
  LOGICAL, PARAMETER :: cp_debug=.TRUE.
#else
  LOGICAL, PARAMETER :: cp_debug=.FALSE.
#endif

!! error number: no error
  INTEGER, PARAMETER :: cp_no_error = 0
!! error number: generic error on the side of the caller
  INTEGER, PARAMETER :: cp_caller_error = 1
!! error number: one or more arguments have and invalid value
  INTEGER, PARAMETER :: cp_wrong_args_error = 100
!! error number: precondition failed
  INTEGER, PARAMETER :: cp_precondition_failed = 200
!! error number: generic error inside the routine
  INTEGER, PARAMETER :: cp_internal_error = -1
!! error number: postcondition failed
  INTEGER, PARAMETER :: cp_postcondition_failed = -200
!! error number: invariant failed
  INTEGER, PARAMETER :: cp_invariant_failed = -100
!! error number: assertion failure
  INTEGER, PARAMETER :: cp_assertion_failed = -300
!! error number: not implemented
  INTEGER, PARAMETER :: cp_unimplemented_error_nr = -1000

!!***
!****************************************************************************

!!****s* cp_error_handling/cp_error_type [1.0] *
!!
!!   NAME
!!     cp_error_type
!!
!!   FUNCTION
!!     this type represent a structure that can be passed to
!!     the functions to change their error behaviour.
!!     It should always be initialized with cp_error_init before use
!!
!!   NOTES
!!     if one of level or error_nr is different from 0 then both must be
!!     different from 0.
!!
!!   ATTRIBUTES
!!     initialized: true if initialized (just for bug catching)
!!     logger: a pointer to the logger where error logs should go 
!!         (if null they go to cp_default_logger)
!!     stop_level: the level at which the execution stops
!!     print_level: starting at which level something gets printed
!!     level: level of the error (0 = cp_note_level=no error, 
!!          cp_warning_level, cp_failure_level or cp_fatal_level=3).
!!          Definition of these constants is in cp_log_handling.
!!     error_nr: the number of the error (0 = no error)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_error_type
     PRIVATE
     LOGICAL :: initialized
     TYPE(cp_logger_type), POINTER :: logger
     INTEGER :: stop_level
     INTEGER :: print_level
     INTEGER :: level
     INTEGER :: error_nr
  END TYPE cp_error_type
!!***
!****************************************************************************

CONTAINS
! ================== error structure init/get/set ===================

!!****f* cp_error_handling/cp_error_init [1.0] *
!!
!!   NAME
!!    cp_error_init
!!
!!   FUNCTION
!!    Initializes the error structure error.
!!    The optional arguments initialize the values in the structure.
!!    Default values: stop_level=2, logger=null (i.e. default logger)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error that will be initialized
!!     stop_level: the stop level of the newly created error (defaults to
!!           cp_failure_level)
!!     print_level: the print level of the newly created error (defaults
!!           to cp_note_level)
!!     logger: the logger where errors are written (defaults to
!!            the default logger)
!!     template_error: an error that gives the default values for 
!!            the previous variables (useful to initialize an internal error)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_init(error, stop_level, print_level, logger, template_error, error_nr, level)


    TYPE(cp_error_type), INTENT(out)         :: error
    INTEGER, INTENT(in), OPTIONAL            :: stop_level, print_level
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: template_error
    INTEGER, INTENT(in), OPTIONAL            :: error_nr, level

    error%initialized=.TRUE.
    error%error_nr=0
    error%level=0
    error%print_level=CP_NOTE_LEVEL
    error%stop_level=CP_FAILURE_LEVEL
    error%print_level=CP_NOTE_LEVEL
    NULLIFY(error%logger)
    IF (PRESENT(template_error)) THEN
       error%logger => template_error%logger
       error%stop_level = template_error%stop_level
       error%print_level = template_error%print_level
    END IF
    IF (PRESENT(logger)) THEN
       error%logger => logger
    END IF
    IF(PRESENT(stop_level)) THEN
       error%stop_level=stop_level
    END IF
    IF (PRESENT(print_level)) THEN
       error%print_level=print_level
    END IF
    IF (PRESENT(level)) error%level=level
    IF(PRESENT(error_nr)) error%error_nr=error_nr
    IF (ASSOCIATED(error%logger)) THEN
       CALL cp_logger_retain(error%logger)
    END IF
  END SUBROUTINE cp_error_init
!****************************************************************************

!!****f* cp_error_handling/cp_error_reset [1.0] *
!!
!!   NAME
!!     cp_error_reset
!!
!!   FUNCTION
!!     if the optional parameter error is present resets its level and 
!!     error_nr to 0 (is tipically called when you want to discard the 
!!     previous errors stored in error, i.e. not often)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you want to reset
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_reset(error)
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    IF (PRESENT(error)) THEN
       IF (.NOT. error%initialized) THEN
          CALL cp_log(cp_get_default_logger(), CP_FAILURE_LEVEL,  &
               'cp_error_handling:cp_error_reset:1', &
               'error structure not initialized')
          CALL cp_error_common_stop('cp_error_reset',&
               'error structure not initialized')
       END IF
       error%error_nr=0
       error%level=0
    END IF
  END SUBROUTINE cp_error_reset
!****************************************************************************

!!****f* cp_error_handling/cp_error_set [1.0] *
!!
!!   NAME
!!     cp_error_set
!!
!!   FUNCTION
!!     routine to set the various components of the error structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you want to change
!!     template_error: gives the default values of logger, stop_level
!!          and print_level (but not of level and error_nr)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_set(error,level,error_nr,stop_level,print_level,logger,&
       template_error)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(in), OPTIONAL            :: level, error_nr, stop_level, &
                                                print_level
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: template_error

    TYPE(cp_logger_type), POINTER            :: old_logger

    old_logger => error%logger
    IF (PRESENT(template_error)) THEN
       error%logger => template_error%logger
       error%stop_level = template_error%stop_level
       error%print_level = template_error%print_level
    END IF
    IF (PRESENT(level)) error%level=level
    IF (PRESENT(error_nr)) error%error_nr=error_nr
    IF (PRESENT(stop_level)) error%stop_level=stop_level
    IF (PRESENT(print_level)) error%print_level=print_level
    IF (PRESENT(logger)) error%logger=>logger
    IF (ASSOCIATED(error%logger)) THEN
       CALL cp_logger_retain(error%logger)
    END IF
    CALL cp_logger_release(old_logger)
  END SUBROUTINE cp_error_set
!***************************************************************************

!!****f* cp_error_handling/cp_error_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_error_dealloc_ref
!!
!!   FUNCTION
!!     Deallocates the memory allocated by the error
!!
!!   NOTES
!!     actually the error allocates no memory, and it is better if it 
!!     stays so even in the future, but I have written this routine 
!!     to remain consistent with the generic types
!!
!!   INPUTS
!!     err_to_free: the error that gets deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_dealloc_ref(err_to_free, error)
    TYPE(cp_error_type), INTENT(inout)       :: err_to_free
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_error_handling:cp_error_dealloc_ref'

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (.NOT.err_to_free%initialized) THEN
       CALL cp_log(cp_get_default_logger(), CP_FAILURE_LEVEL,  &
            'cp_error_handling:cp_error_dealloc_ref:1', &
            'error structure not initialized')
       CALL cp_error_common_stop('cp_error_reset',&
            'error structure not initialized')
    END IF
    IF (.NOT. failure) THEN
       CALL cp_logger_release(err_to_free%logger)
    END IF
  END SUBROUTINE cp_error_dealloc_ref
!***************************************************************************

!!****f* cp_error_handling/cp_error_get [1.0] *
!!
!!   NAME
!!     cp_error_get
!!
!!   FUNCTION
!!     returns the various attributes of the given error
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
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_get(error,level,error_nr,stop_level,print_level,logger)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(out), OPTIONAL           :: level, error_nr, stop_level, &
                                                print_level
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_error_handling:cp_error_get'

    IF (PRESENT(level)) level=error%level
    IF (PRESENT(error_nr)) error_nr=error%error_nr
    IF (PRESENT(stop_level)) stop_level=error%stop_level
    IF (PRESENT(print_level)) print_level=error%print_level
    IF (PRESENT(logger)) logger=>cp_error_get_logger(error)
  END SUBROUTINE cp_error_get
!***************************************************************************

!!****f* cp_error_handling/cp_error_get_level [1.0] *
!!
!!   NAME
!!     cp_error_get_level
!!
!!   FUNCTION
!!     gets the error level from the given error
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_error_get_level(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    IF (PRESENT(error)) THEN
       res=error%level
    ELSE
       res=0
    END IF
  END FUNCTION cp_error_get_level
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_print_level [1.0] *
!!
!!   NAME
!!     cp_error_get_print_level
!!
!!   FUNCTION
!!     gets the print level (starting with which errorlevel something 
!!     gets printed) from the given error
!!     if the error is not present returns cp_note_level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version .2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  ! pure 
  FUNCTION cp_error_get_print_level(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    IF (PRESENT(error)) THEN
       res=error%print_level
    ELSE
       res=cp_note_level
    END IF
  END FUNCTION cp_error_get_print_level
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_nr [1.0] *
!!
!!   NAME
!!     cp_error_get_nr
!!
!!   FUNCTION
!!     gets the error number from the given error
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  ! pure 
  FUNCTION cp_error_get_nr(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    IF (PRESENT(error)) THEN
       res=error%error_nr
    ELSE
       res=0
    END IF
  END FUNCTION cp_error_get_nr
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_logger [1.0] *
!!
!!   NAME
!!     cp_error_get_logger
!!
!!   FUNCTION
!!     gets the logger of the given error.
!!     The error argument is optional, so you can get the logger
!!     of an optional dummy argument.
!!     If the error is not present or its logger is not associated
!!     returns the default logger.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_error_get_logger(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    TYPE(cp_logger_type), POINTER            :: res

    NULLIFY(res)
    IF (PRESENT(error)) res=>error%logger
    IF (.NOT. ASSOCIATED(res)) res => cp_get_default_logger()
    IF (.NOT.ASSOCIATED(res)) THEN
       PRINT *,'ERROR in cp_error_handling:cp_error_get_logger, cp_default_logger is not associated'
       CALL cp_error_common_stop('cp_error_handling:cp_error_get_logger',' cp_default_logger is not associated')
    END IF
  END FUNCTION cp_error_get_logger
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_stop_level [1.0] *
!!
!!   NAME
!!     cp_error_get_stop_level
!!
!!   FUNCTION
!!     gets the error level at which the error will stop the programm
!!     If the error is not present returns cp_failure_level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  ! pure 
  FUNCTION cp_error_get_stop_level(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    IF (PRESENT(error)) THEN
       res=error%stop_level
    ELSE
       res=cp_failure_level
    END IF
  END FUNCTION cp_error_get_stop_level

! =================== error handling ===================

!!****f* cp_error_handling/cp_error_common_stop [1.0] *
!!
!!   NAME
!!     cp_error_common_stop
!!
!!   FUNCTION
!!     a subroutine that call just stop, useful because all the failures
!!     fails calling this routine: the right point to set a breakpoint
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     fromWhere: 'module:routine' that called this function
!!     message: an error message
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_common_stop(fromWhere,message)
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_unit_nr()
    CALL print_stack(unit_nr)
    CALL print_message("ERROR in "//TRIM(fromWhere),unit_nr,2,2,0)
    CALL print_message(message,unit_nr,1,1,1)
    CALL print_stack(unit_nr)
    CALL mp_abort()
    STOP ! just to be sure...
  END SUBROUTINE cp_error_common_stop
!****************************************************************************

!!****f* cp_error_handling/cp_error_handle_error [1.0] *
!!
!!   NAME
!!     cp_error_handle_error
!!
!!   FUNCTION
!!     all the error situations at the end should call this routine
!!
!!     main function for error logging, most error reporting ends up calling
!!     this function, great place to put a breakpoint if you want to catch
!!     errors that do not stop the program.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error you take the info from (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_handle_error(error, fromWhere, message)
    TYPE(cp_error_type), INTENT(in)          :: error
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message

    INTEGER                                  :: mepos
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF (debug_this_module) PRINT *, "entred in cp_error_handle_error"
    
    ! print stack only when stopping? print to stdout?
    CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))

    CALL cp_error_message(error%level, fromWhere, &
         'err='//ADJUSTL(cp_to_string(error%error_nr))//message, error)
    IF (error%level >= error%stop_level) THEN
       logger => cp_error_get_logger(error)
       mepos=logger%para_env%mepos
       CALL cp_error_common_stop(fromWhere//' processor '//&
            cp_to_string(mepos)//' ',message)
    END IF
  END SUBROUTINE cp_error_handle_error
!***************************************************************************

!!****f* cp_error_handling/cp_error_message [1.0] *
!!
!!   NAME
!!     cp_error_message
!!
!!   FUNCTION
!!     writes an error message of level level, subject to control of error
!!     (i.e. if error%printlevel is >= level)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     level: the level of the error (cp_note_level, cp_warning_level,
!!     fromWhere: a string that contains the module name and routine name
!!           where this test happened
!!     message: the error message
!!     error: the error type that controls error reporting (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_message(level, fromWhere, message, error)
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error

    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF (cp_error_get_print_level(error) <= level) THEN
       CALL cp_log(logger, level, fromWhere, message,local=.TRUE.)
    END IF
  END SUBROUTINE cp_error_message
!***************************************************************************

!!****f* cp_error_handling/cp_assert [1.0] *
!!
!!   NAME
!!     cp_assert
!!
!!   FUNCTION
!!     routine that checks an assertion, if it fails set the optional
!!     argument failure is given sets it to true.
!!     This can be used to do many assertion one fater the other and 
!!     check at the end if one failed.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     condition: the condition to check, if true does nothing
!!     level: the level of the error (cp_note_level, cp_warning_level,
!!     error_nr: the number of the error
!!     fromWhere: a string that contains the module name and routine name
!!           where this test happened
!!     message: the error message
!!     error: the error type that controls error reporting (optional)
!!     failure: is set to true if the condition is false otherwise
!!           it is not changed (optional)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_assert(condition, level, error_nr, fromWhere, &
       message, error,failure)
!! the condition that is checked, if false, an error did happen
    LOGICAL, INTENT(in)                      :: condition
    INTEGER, INTENT(in)                      :: level, error_nr
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_error_handling:cp_assert'

    TYPE(cp_error_type)                      :: err

!! level of the error
!! error number
!! is set to true if condition is false

    IF (.NOT. condition) THEN
       IF (debug_this_module) PRINT *, "assertion failed, in "//fromWhere//" msg="//message
       IF (PRESENT(failure)) THEN
          failure=.TRUE.
       END IF
       CALL cp_error_init(err, level=level,error_nr=error_nr, &
            template_error=error)
       IF (PRESENT(error)) THEN
          IF (error%level < level) THEN
             error%level=level
             error%error_nr=error_nr
          END IF
       END IF
       CALL cp_error_handle_error(err, fromWhere,message)
       CALL cp_error_dealloc_ref(err)
    ELSE
       IF (debug_this_module) PRINT *, "assertion checked, in "//fromWhere//" msg="//message
    END IF
  END SUBROUTINE cp_assert
!***************************************************************************

!!****f* cp_error_handling/cp_error_propagate_error [1.0] *
!!
!!   NAME
!!     cp_error_propagate_error
!!
!!   FUNCTION
!!     if the internal error iError of a subroutine has an error levele 
!!     that is at least failure_level(i.e. if there was a serious error
!!     in the subroutine) propagates the error to the actual error
!!     (i.e. sets error to the same error than iError) and sets failure
!!     to true
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     iError: the error object that comes from the subroutine
!!     error: the error object in the actual routine
!!     failure: if iError has an error it is set to true
!!     failure_level: starting with which error level propagation happens 
!!           (defaults to cp_warning_level)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_error_propagate_error(iError, fromWhere, message, error, &
       failure, failure_level)
    TYPE(cp_error_type), INTENT(in)          :: iError
    CHARACTER(len=*)                         :: fromWhere
    CHARACTER(len=*), OPTIONAL               :: message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure
    INTEGER, INTENT(in), OPTIONAL            :: failure_level

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_error_handling:cp_error_propagate_error'

    INTEGER                                  :: f_level

    f_level=cp_warning_level
    IF (PRESENT(failure_level)) f_level=failure_level
    IF (cp_error_get_level(iError) /= 0 .AND. iError%level>=f_level) THEN
       IF (PRESENT(failure)) &
            failure=.TRUE.
       IF (PRESENT(error)) THEN
          IF (error%level < iError%level) THEN
             error%level=iError%level
             error%error_nr=iError%error_nr
          END IF
       END IF
       IF (iError%level >= cp_error_get_print_level(error)) THEN
          IF (PRESENT(message)) THEN
             CALL cp_error_message(iError%level, fromWhere,"propagated error "//&
                  message, error)
          ELSE
             CALL cp_error_message(iError%level, fromWhere,"propagated error",error)
          END IF
       END IF
    END IF
  END SUBROUTINE cp_error_propagate_error
!***************************************************************************

!!****f* cp_error_handling/cp_unimplemented_error [1.0] *
!!
!!   NAME
!!     cp_unimplemented_error
!!
!!   FUNCTION
!!     comodity call to signat that something is not implemented
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_unimplemented_error(fromWhere, message, error, error_level)
    CHARACTER(len=*), INTENT(in)             :: fromWhere
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, INTENT(in), OPTIONAL            :: error_level

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unimplemented_error', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: err_l

  err_l=cp_failure_level

  IF(PRESENT(error_level)) err_l=error_level
  IF (PRESENT(message)) THEN
     CALL cp_assert(.FALSE.,err_l,cp_unimplemented_error_nr,&
          fromWhere, "UNIMPLEMENTED, "//message, error=error)
  ELSE
     CALL cp_assert(.FALSE.,err_l,cp_unimplemented_error_nr,&
          fromWhere, "UNIMPLEMENTED", error=error)
  END IF
END SUBROUTINE cp_unimplemented_error
!***************************************************************************

!!****f* cp_error_handling/cp_error_synchronize_error [1.0] *
!!
!!   NAME
!!     cp_error_synchronize_error
!!
!!   FUNCTION
!!     synchronizes the error values between different processors
!!     (propagating an eventual failure to each processor)
!!
!!   NOTES
!!     allow for non uniform presence of the error argument?)
!!
!!   ARGUMENTS
!!     - error: the error to synchronize
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_error_synchronize_error(error,para_env)
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_synchronize_error', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: codes
    INTEGER, DIMENSION(2)                    :: my_codes
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: my_para_env

  failure=.FALSE.
  NULLIFY(my_para_env)
  IF (PRESENT(error)) THEN
     logger => cp_error_get_logger(error)
     my_para_env => logger%para_env
     IF (PRESENT(para_env)) my_para_env => para_env
     CALL cp_assert(ASSOCIATED(my_para_env),cp_failure_level,cp_failure_level,&
          routineP,"cannot synchronize without a valid para_env",error,failure)
     IF (my_para_env%num_pe>1) THEN
        ALLOCATE(codes(my_para_env%num_pe*2),stat=stat)
        my_codes(1)=error%level
        my_codes(2)=error%error_nr
        CALL mp_alltoall( my_codes, codes, 2, my_para_env%group )
        error%level=codes(1)
        error%error_nr=codes(2)
        DO i=1,my_para_env%num_pe-1
           IF (codes(i*2+1)> error%level) THEN
              error%level=codes(i*2+1)
              error%error_nr=codes(i*2+2)
           END IF
        END DO
     END IF
  END IF
END SUBROUTINE cp_error_synchronize_error
!***************************************************************************

!!****f* cp_error_handling/cp_a_l [1.0] *
!!
!!   NAME
!!     cp_a_l
!!
!!   FUNCTION
!!     assert with linenumber (to be used with the assertion macros
!!     to keep the line short)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - from_
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_a_l(condition, level, fromWhere, lineNr, error,failure)
    LOGICAL, INTENT(in)                      :: condition
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere
    INTEGER, INTENT(in)                      :: lineNr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure

  CALL cp_assert(condition,level=level,error_nr=cp_assertion_failed,&
       fromWhere=fromWhere,&
       message="condition FAILED at line "//ADJUSTL(cp_to_string(lineNr)),&
       error=error, failure=failure)
END SUBROUTINE cp_a_l
!***************************************************************************

!!****f* cp_error_handling/cp_error_check [1.0] *
!!
!!   NAME
!!     cp_error_check
!!
!!   FUNCTION
!!     if error is present sets failure to true if an error did happen
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: the error to check for failure
!!     - failure: set to true if there was an error
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_error_check(error, failure)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_check', &
      routineP = moduleN//':'//routineN

  failure=.FALSE.
  
  IF (PRESENT(error).and.PRESENT(failure)) THEN
     IF (error%error_nr /=0 .OR. error%level /=0) THEN
        failure=.TRUE.
     END IF
  END IF
END SUBROUTINE cp_error_check
!***************************************************************************

!!****f* cp_error_handling/cp_simple_assert *
!!
!!   NAME
!!     cp_simple_assert
!!
!!   FUNCTION
!!     simple assert
!!
!!   ARGUMENTS
!!     -
!!
!!   NOTES
!!     you should use the CPAssert and similar macros defined in 
!!     cp_prep_globals.h that call cp_a_l
!!
!!   MODIFICATION HISTORY
!!     , TCH, created
!!
!!   SEE ALSO
!!     cp_a_l
!!
!!****
  SUBROUTINE cp_simple_assert(condition, routineP, line)
  

    LOGICAL, INTENT(in)                      :: condition
    CHARACTER(len=*), INTENT(in)             :: routineP
    INTEGER, INTENT(in)                      :: line

! code

    IF (.NOT.condition) THEN
       WRITE (*,*) "======== ASSERTION FAILED ========"
       WRITE (*,*) routineP, ": ", cp_to_string(line)
       STOP 1
    END IF
  
  END SUBROUTINE cp_simple_assert

! *****************************************************************************

END MODULE cp_error_handling

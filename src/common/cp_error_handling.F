!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that contains the routines for error handling
!>      The idea is that in general routines have an optional
!>      parameter that can be used to change the error behaviour.
!> \note
!>      In general 0 means no error, a negative error code means an internal
!>      error (in the routine, postcondition failure), and a positive error
!>      an error on the side of the caller (precondition failure).
!>      Error numbers in the dange -1000..1000 are reserved for generic error
!>      codes.
!>      Inspired from the nag nl90 library.
!>      Some help on this is also available in the high level dev notes on pao.
!>
!>      The following error levels are defined:
!>      0 - note    : CP_NOTE_LEVEL
!>      1 - warning : CP_WARNING_LEVEL
!>      2 - failure : CP_FAILURE_LEVEL
!>      3 - fatal   : CP_FATAL_LEVEL
!>      You shouldn't use the numbers but the constants defined in
!>      cp_log_handling
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
MODULE cp_error_handling
  USE cp_log_handling,                 ONLY: &
       CP_FAILURE_LEVEL, CP_FATAL_LEVEL, CP_NOTE_LEVEL, CP_WARNING_LEVEL, &
       cp_get_default_logger, cp_log, cp_logger_get_default_unit_nr, &
       cp_logger_release, cp_logger_retain, cp_logger_type, cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_max
  USE termination,                     ONLY: print_message
  USE timings,                         ONLY: print_stack

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_error_handling'

  !API types
  PUBLIC :: cp_error_type
  !API parameters var
  PUBLIC :: cp_debug, cp_no_error, cp_caller_error, cp_wrong_args_error,&
       cp_precondition_failed, cp_internal_error, cp_postcondition_failed,&
       cp_invariant_failed, cp_assertion_failed, cp_unimplemented_error_nr
  !API public routines
  PUBLIC :: cp_assert, cp_a_l, cp_simple_assert, cp_unimplemented_error

  !API
  PUBLIC :: cp_error_init, cp_error_dealloc_ref, cp_error_set,&
       cp_error_get, cp_error_reset, cp_error_get_level,&
       cp_error_get_print_level, cp_error_get_nr,&
       cp_error_get_logger, cp_error_get_stop_level,&
       cp_error_handle_error,&
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

! *****************************************************************************
!> \brief this type represent a structure that can be passed to
!>      the functions to change their error behaviour.
!>      It should always be initialized with cp_error_init before use
!> \param initialized true if initialized (just for bug catching)
!> \param logger a pointer to the logger where error logs should go
!>          (if null they go to cp_default_logger)
!> \param stop_level the level at which the execution stops
!> \param print_level starting at which level something gets printed
!> \param level level of the error (0 = cp_note_level=no error,
!>           cp_warning_level, cp_failure_level or cp_fatal_level=3).
!>           Definition of these constants is in cp_log_handling.
!> \param error_nr the number of the error (0 = no error)
!> \param info some additional info to store (error message, string ...)
!> \note
!>      if one of level or error_nr is different from 0 then both must be
!>      different from 0.
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  TYPE cp_error_type
     PRIVATE
     LOGICAL :: initialized
     TYPE(cp_logger_type), POINTER :: logger => Null()
     INTEGER :: stop_level
     INTEGER :: print_level
     INTEGER :: level
     INTEGER :: error_nr
     CHARACTER(LEN=512) :: info
  END TYPE cp_error_type

CONTAINS
! ================== error structure init/get/set ===================

! *****************************************************************************
!> \brief Initializes the error structure error.
!>     The optional arguments initialize the values in the structure.
!>     Default values: stop_level=2, logger=null (i.e. default logger)
!> \param error the error that will be initialized
!> \param stop_level the stop level of the newly created error (defaults to
!>            cp_failure_level)
!> \param print_level the print level of the newly created error (defaults
!>            to cp_note_level)
!> \param logger the logger where errors are written (defaults to
!>             the default logger)
!> \param template_error an error that gives the default values for
!>             the previous variables (useful to initialize an internal error)
!> \param error_nr ...
!> \param level ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
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
    error%info="No info available yet, we will keep you informed."
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

! *****************************************************************************
!> \brief if the optional parameter error is present resets its level and
!>      error_nr to 0 (is tipically called when you want to discard the
!>      previous errors stored in error, i.e. not often)
!> \param error the error you want to reset
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief routine to set the various components of the error structure
!> \param error the error you want to change
!> \param level ...
!> \param error_nr ...
!> \param stop_level ...
!> \param print_level ...
!> \param logger ...
!> \param template_error gives the default values of logger, stop_level
!>           and print_level (but not of level and error_nr)
!> \param info ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_error_set(error,level,error_nr,stop_level,print_level,logger,&
       template_error,info)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(in), OPTIONAL            :: level, error_nr, stop_level, &
                                                print_level
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: template_error
    CHARACTER(LEN=*), OPTIONAL               :: info

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
    IF (PRESENT(info)) error%info=info
    CALL cp_logger_release(old_logger)
  END SUBROUTINE cp_error_set

! *****************************************************************************
!> \brief Deallocates the memory allocated by the error
!> \param err_to_free the error that gets deallocated
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
!> \note
!>      actually the error allocates no memory, and it is better if it
!>      stays so even in the future, but I have written this routine
!>      to remain consistent with the generic types
! *****************************************************************************
  SUBROUTINE cp_error_dealloc_ref(err_to_free, error)
    TYPE(cp_error_type), INTENT(inout)       :: err_to_free
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_dealloc_ref', &
      routineP = moduleN//':'//routineN

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

! *****************************************************************************
!> \brief returns the various attributes of the given error
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \param level ...
!> \param error_nr ...
!> \param stop_level ...
!> \param print_level ...
!> \param logger ...
!> \param info ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  SUBROUTINE cp_error_get(error,level,error_nr,stop_level,print_level,logger,info)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(out), OPTIONAL           :: level, error_nr, stop_level, &
                                                print_level
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    CHARACTER(LEN=*), OPTIONAL               :: info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_get', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(level)) level=error%level
    IF (PRESENT(error_nr)) error_nr=error%error_nr
    IF (PRESENT(stop_level)) stop_level=error%stop_level
    IF (PRESENT(print_level)) print_level=error%print_level
    IF (PRESENT(logger)) logger=>cp_error_get_logger(error)
    IF (PRESENT(info)) info=error%info
  END SUBROUTINE cp_error_get

! *****************************************************************************
!> \brief gets the error level from the given error
!> \param error the error you take the info from
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief gets the print level (starting with which errorlevel something
!>      gets printed) from the given error
!>      if the error is not present returns cp_note_level
!> \param error the error you take the info from (optional)
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version .2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief gets the error number from the given error
!> \param error the error you take the info from
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief gets the logger of the given error.
!>      The error argument is optional, so you can get the logger
!>      of an optional dummy argument.
!>      If the error is not present or its logger is not associated
!>      returns the default logger.
!> \param error the error you take the info from (optional)
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  FUNCTION cp_error_get_logger(error) RESULT(res)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    TYPE(cp_logger_type), POINTER            :: res

    NULLIFY(res)
    IF (PRESENT(error)) res=>error%logger
    IF (.NOT.ASSOCIATED(res)) res => cp_get_default_logger()
    IF (.NOT.ASSOCIATED(res)) THEN
       CALL cp_error_common_stop('cp_error_handling:cp_error_get_logger',&
            ' cp_default_logger is not associated')
    END IF
  END FUNCTION cp_error_get_logger

! *****************************************************************************
!> \brief gets the error level at which the error will stop the programm
!>      If the error is not present returns cp_failure_level
!> \param error the error you take the info from (optional)
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief a subroutine that call just stop, useful because all the failures
!>      fails calling this routine: the right point to set a breakpoint
!> \param fromWhere 'module:routine' that called this function
!> \param message an error message
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_error_common_stop(fromWhere,message)
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_unit_nr()
    CALL print_message("ERROR in "//TRIM(ADJUSTL(fromWhere))//&
                       TRIM(ADJUSTL(message)),unit_nr,2,2,2)
    CALL mp_abort()

  END SUBROUTINE cp_error_common_stop

! *****************************************************************************
!> \brief all the error situations at the end should call this routine
!>
!>      main function for error logging, most error reporting ends up calling
!>      this function, great place to put a breakpoint if you want to catch
!>      errors that do not stop the program.
!> \param error the error you take the info from (optional)
!> \param fromWhere ...
!> \param message ...
!> \param ionode ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_error_handle_error(error, fromWhere, message, ionode)
    TYPE(cp_error_type), INTENT(in)          :: error
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    LOGICAL, INTENT(IN)                      :: ionode

    INTEGER                                  :: mepos
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF (debug_this_module) PRINT *, "entred in cp_error_handle_error"

    IF (ionode.AND.(error%level < error%stop_level)) THEN
       CALL cp_error_message(error%level,fromWhere,' '//TRIM(ADJUSTL(message)),error)
!MK No need to frighten anybody with an error code when printing an informal note
!MK       CALL cp_error_message(error%level,fromWhere,&
!MK                             'err='//TRIM(ADJUSTL(cp_to_string(error%error_nr)))//&
!MK                             ' '//TRIM(ADJUSTL(message)),error)
    END IF

    IF (error%level >= error%stop_level) THEN
       mepos=logger%para_env%mepos
       CALL cp_error_message(error%level,fromWhere//' processor '//&
                             TRIM(ADJUSTL(cp_to_string(mepos)))//' ',&
                             'err='//TRIM(ADJUSTL(cp_to_string(error%error_nr)))//&
                             ' '//TRIM(ADJUSTL(message)),error)
       CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       CALL mp_abort(message)
    END IF

  END SUBROUTINE cp_error_handle_error

! *****************************************************************************
!> \brief writes an error message of level level, subject to control of error
!>      (i.e. if error%printlevel is >= level)
!> \param level the level of the error (cp_note_level, cp_warning_level,
!> \param fromWhere a string that contains the module name and routine name
!>            where this test happened
!> \param message the error message
!> \param error the error type that controls error reporting (optional)
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
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

! *****************************************************************************
!> \brief routine that checks an assertion, if it fails set the optional
!>      argument failure is given sets it to true.
!>      This can be used to do many assertion one fater the other and
!>      check at the end if one failed.
!> \param condition the condition to check, if true does nothing
!> \param level the level of the error (cp_note_level, cp_warning_level,
!> \param error_nr the number of the error
!> \param fromWhere a string that contains the module name and routine name
!>            where this test happened
!> \param message the error message
!> \param error the error type that controls error reporting (optional)
!> \param failure is set to true if the condition is false otherwise
!>            it is not changed (optional)
!> \param info ...
!> \param only_ionode ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  SUBROUTINE cp_assert(condition, level, error_nr, fromWhere, &
       message, error,failure,info,only_ionode)
    !! the condition that is checked, if false, an error did happen
    LOGICAL, INTENT(in)                      :: condition
    INTEGER, INTENT(in)                      :: level, error_nr
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure
    CHARACTER(len=*), OPTIONAL               :: info
    LOGICAL, INTENT(in), OPTIONAL            :: only_ionode

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_assert', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: ionode
    TYPE(cp_error_type)                      :: err
    TYPE(cp_logger_type), POINTER            :: logger

    ionode = .TRUE.
    IF (.NOT. condition) THEN
       IF (debug_this_module) PRINT *, "assertion failed, in "//fromWhere//" msg="//message
       IF (PRESENT(failure)) THEN
          failure=.TRUE.
       END IF
       CALL cp_error_init(err, level=level,error_nr=error_nr, &
            template_error=error)
       IF (PRESENT(only_ionode)) THEN
          IF (only_ionode) THEN
             NULLIFY(logger)
             logger => cp_error_get_logger(err)
             ionode = logger%para_env%mepos==logger%para_env%source
          END IF
       END IF
       IF (PRESENT(error)) THEN
          IF (error%level < level) THEN
             error%level=level
             error%error_nr=error_nr
          END IF
          IF (PRESENT(info)) error%info=info
       END IF
       CALL cp_error_handle_error(err, fromWhere,message,ionode)
       CALL cp_error_dealloc_ref(err)
    ELSE
       IF (debug_this_module) PRINT *, "assertion checked, in "//fromWhere//" msg="//message
    END IF
  END SUBROUTINE cp_assert

! *****************************************************************************
!> \brief if the internal error iError of a subroutine has an error levele
!>      that is at least failure_level(i.e. if there was a serious error
!>      in the subroutine) propagates the error to the actual error
!>      (i.e. sets error to the same error than iError) and sets failure
!>      to true if a failure more severe than warning happened
!> \param iError the error object that comes from the subroutine
!> \param fromWhere ...
!> \param message ...
!> \param error the error object in the actual routine
!> \param failure if iError has an error it is set to true
!> \param failure_level starting with which error level propagation happens
!>            (defaults to cp_warning_level)
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  SUBROUTINE cp_error_propagate_error(iError, fromWhere, message, error, &
       failure, failure_level)
    TYPE(cp_error_type), INTENT(in)          :: iError
    CHARACTER(len=*)                         :: fromWhere
    CHARACTER(len=*), OPTIONAL               :: message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure
    INTEGER, INTENT(in), OPTIONAL            :: failure_level

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_propagate_error', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f_level

    f_level=cp_warning_level
    IF (PRESENT(failure_level)) f_level=failure_level
    IF (cp_error_get_level(iError) /= 0 .AND. iError%level>f_level) THEN
       IF (PRESENT(failure)) &
            failure=.TRUE.
       IF (PRESENT(error)) THEN
          IF (error%level < iError%level) THEN
             error%level=iError%level
             error%error_nr=iError%error_nr
          END IF
          error%info=iError%info
       END IF
       IF (debug_this_module) THEN
          IF (iError%level >= cp_error_get_print_level(error)) THEN
             IF (PRESENT(message)) THEN
                CALL cp_error_message(iError%level, fromWhere,"propagated error "//&
                     message, error)
             ELSE
                CALL cp_error_message(iError%level, fromWhere,"propagated error",error)
             END IF
          END IF
       ENDIF
    END IF
  END SUBROUTINE cp_error_propagate_error

! *****************************************************************************
!> \brief comodity call to signat that something is not implemented
!> \param fromWhere routine where the error happened
!> \param message the message to write out (UNIMPLEMENTED is added to it)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param error_level the error level of this error, defaults to
!>        cp_failure_level
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
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

! *****************************************************************************
!> \brief synchronizes the error values between different processors
!>      (propagating an eventual failure to each processor)
!> \param error the error to synchronize
!> \param para_env ...
!> \author fawzi
!> \note
!>      allow for non uniform presence of the error argument?)
! *****************************************************************************
  SUBROUTINE cp_error_synchronize_error(error,para_env)
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_synchronize_error', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)                    :: my_codes
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: my_para_env

    failure=.FALSE.
    NULLIFY(my_para_env,logger)
    IF (PRESENT(error)) THEN
       IF (PRESENT(para_env)) THEN
          my_para_env => para_env
       ELSE
          logger => cp_error_get_logger(error)
          my_para_env => logger%para_env
       END IF
       CALL cp_assert(ASSOCIATED(my_para_env),cp_failure_level,cp_failure_level,&
            routineP,"cannot synchronize without a valid para_env",error,failure)
       IF (my_para_env%num_pe>1) THEN
          my_codes(1)=ABS(error%level)
          my_codes(2)=ABS(error%error_nr)
          CALL mp_max( my_codes, my_para_env%group )
          error%level=my_codes(1)
          error%error_nr=my_codes(2)
       END IF
    END IF
  END SUBROUTINE cp_error_synchronize_error

! *****************************************************************************
!> \brief assert with linenumber (to be used with the assertion macros
!>      to keep the line short)
!> \param condition ...
!> \param level ...
!> \param fromWhere ...
!> \param lineNr ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param failure ...
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
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
         message="condition FAILED at line "//TRIM(ADJUSTL(cp_to_string(lineNr))),&
         error=error, failure=failure)
  END SUBROUTINE cp_a_l

! *****************************************************************************
!> \brief if error is present sets failure to true if an error did happen
!> \param error the error to check for failure
!> \param failure set to true if there was an error
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_error_check(error, failure)
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_error_check', &
      routineP = moduleN//':'//routineN

    failure=.FALSE.

    IF (PRESENT(error).and.PRESENT(failure)) THEN
       IF (error%error_nr /=0 .AND. error%level > CP_WARNING_LEVEL) THEN
          failure=.TRUE.
       END IF
    END IF
  END SUBROUTINE cp_error_check

! *****************************************************************************
!> \brief simple assert
!> \param condition ...
!> \param routineP ...
!> \param line ...
!> \par History
!>      TCH, created
!> \note
!>      you should use the CPAssert and similar macros defined in
!>      cp_prep_globals.h that call cp_a_l
! *****************************************************************************
  SUBROUTINE cp_simple_assert(condition, routineP, line)

    LOGICAL, INTENT(in)                      :: condition
    CHARACTER(len=*), INTENT(in)             :: routineP
    INTEGER, INTENT(in)                      :: line

    CHARACTER(LEN=100)                       :: str

    IF (.NOT.condition) THEN
       WRITE (*,*) "======== ASSERTION FAILED ========"
       str=cp_to_string(line)
       WRITE (*,*) routineP, ": ", LEN_TRIM(str)
       STOP 1
    END IF

  END SUBROUTINE cp_simple_assert

END MODULE cp_error_handling

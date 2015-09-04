!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE message_passing,                 ONLY: mp_abort
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
  PUBLIC :: cp_assert, cp_a_l, cp_unimplemented_error

  !API
  PUBLIC :: cp_error_init, cp_error_dealloc_ref, &
       cp_error_get_logger,&
       cp_error_message

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
     TYPE(cp_logger_type), POINTER :: logger => Null()
  END TYPE cp_error_type

CONTAINS
! ================== error structure init/get/set ===================

! *****************************************************************************
!> \brief Initializes the error structure error.
!>     The optional arguments initialize the values in the structure.
!>     Default values: stop_level=2, logger=null (i.e. default logger)
!> \param error the error that will be initialized
!> \param logger the logger where errors are written (defaults to
!>             the default logger)
!> \param template_error an error that gives the default values for
!>             the previous variables (useful to initialize an internal error)
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! *****************************************************************************
  SUBROUTINE cp_error_init(error, logger, template_error)

    TYPE(cp_error_type), INTENT(out)         :: error
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    TYPE(cp_error_type), INTENT(in), &
      OPTIONAL                               :: template_error

    NULLIFY(error%logger)
    IF (PRESENT(template_error)) THEN
       error%logger => template_error%logger
    END IF
    IF (PRESENT(logger)) THEN
       error%logger => logger
    END IF
    IF (ASSOCIATED(error%logger)) THEN
       CALL cp_logger_retain(error%logger)
    END IF
  END SUBROUTINE cp_error_init


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
    IF(ASSOCIATED(err_to_free%logger))&
      CALL cp_logger_release(err_to_free%logger)
  END SUBROUTINE cp_error_dealloc_ref

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
    CALL cp_log(logger, level, fromWhere, message,local=.TRUE.)
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
!> \param only_ionode ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  SUBROUTINE cp_assert(condition, level, error_nr, fromWhere, &
       message, error,failure,only_ionode)
    !! the condition that is checked, if false, an error did happen
    LOGICAL, INTENT(in)                      :: condition
    INTEGER, INTENT(in)                      :: level, error_nr
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure
    LOGICAL, INTENT(in), OPTIONAL            :: only_ionode

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_assert', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: iam_ionode, myonly_ionode
    TYPE(cp_logger_type), POINTER            :: logger

    IF (.NOT. condition) THEN
       IF (PRESENT(failure)) failure=.TRUE.
       logger => cp_error_get_logger(error)

       myonly_ionode = .FALSE.
       IF(PRESENT(only_ionode)) myonly_ionode = only_ionode
       iam_ionode = .TRUE.
       IF(myonly_ionode) iam_ionode = logger%para_env%mepos==logger%para_env%source
       IF(level > CP_WARNING_LEVEL) iam_ionode = .TRUE.

       IF (iam_ionode) THEN
          CALL cp_log(logger,level,fromWhere//' processor '//&
                      TRIM(ADJUSTL(cp_to_string(logger%para_env%mepos)))//' ',&
                      'err='//TRIM(ADJUSTL(cp_to_string(error_nr)))//&
                      ' '//TRIM(ADJUSTL(message)), local=.TRUE.)
       ENDIF

       IF(level > CP_WARNING_LEVEL) THEN
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
          CALL mp_abort(message)
       ENDIF
    END IF
  END SUBROUTINE cp_assert

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
!> \brief assert with linenumber (to be used with the assertion macros
!>      to keep the line short)
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
  SUBROUTINE cp_a_l(level, fromWhere, lineNr, error, failure)
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere
    INTEGER, INTENT(in)                      :: lineNr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(inout), OPTIONAL         :: failure

    CALL cp_assert(.FALSE.,level=level,error_nr=cp_assertion_failed,&
         fromWhere=fromWhere,&
         message="condition FAILED at line "//TRIM(ADJUSTL(cp_to_string(lineNr))),&
         error=error, failure=failure)
  END SUBROUTINE cp_a_l

END MODULE cp_error_handling

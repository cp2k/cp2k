!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_error_handling [1.0] *
!!
!!   NAME
!!     qs_error_handling
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
!!     There are the following error levels: 0 note, 1 warning, 2 failure,
!!     3 fatal. You should use the constants defined here, not directly 
!!     the numbers
!!     @see qs_log_handling
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
module qs_error_handling
  use qs_log_handling
  use termination, only: stop_program
  implicit none
  logical, private, parameter :: debug_this_module=.false.
  public

#ifdef FD_DEBUG
!! global variable to turn on/off the debugging, as i am lazy I use the 
!! preprocessor variable FD_DEBUG to decide if it should be on or off.
!! this way it is always right (on in debug builds, off in optimized ones)
  logical, parameter :: qs_debug=.true.
#else
  logical, parameter :: qs_debug=.false.
#endif

!! error number: no error
  integer, parameter :: qs_no_error = 0
!! error number: generic error on the side of the caller
  integer, parameter :: qs_caller_error = 1
!! error number: one or more arguments have and invalid value
  integer, parameter :: qs_wrong_args_error = 100
!! error number: precondition failed
  integer, parameter :: qs_precondition_failed = 200
!! error number: generic error inside the routine
  integer, parameter :: qs_internal_error = -1
!! error number: postcondition failed
  integer, parameter :: qs_postcondition_failed = -200
!! error number: invariant failed
  integer, parameter :: qs_invariant_failed = -100
!! error number: assertion failure
  integer, parameter :: qs_assertion_failed = -300

!! some of the core methods are avaliable for the error:
  interface qs_init
     module procedure qs_init_error
  end interface
  interface qs_dealloc_ref
     module procedure qs_dealloc_error_ref
  end interface
  interface qs_set
     module procedure qs_set_error
  end interface
  interface qs_get
     module procedure qs_error_get
  end interface
!!***
!****************************************************************************

!!****s* qs_error_handling/qs_error_type [1.0] *
!!
!!   NAME
!!     qs_error_type
!!
!!   FUNCTION
!!     this type represent a structure that can be passed to
!!     the functions to change their error behaviour.
!!     It should always be initialized with qs_init_error before use
!!
!!   NOTES
!!     if one of level or error_nr is different from 0 then both must be
!!     different from 0.
!!
!!   ATTRIBUTES
!!     initialized: true if initialized (just for bug catching)
!!     logger: a pointer to the logger where error logs should go 
!!         (if null they go to qs_default_logger)
!!     stop_level: the level at which the execution stops
!!     print_level: starting at which level something gets printed
!!     level: level of the error (0 = qs_note_level=no error, 
!!          qs_warning_level, qs_failure_level or qs_fatal_level=3).
!!          Definition of these constants is in qs_log_handling.
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
  type qs_error_type
     private
     logical :: initialized
     type(qs_logger), pointer :: logger
     integer :: stop_level
     integer :: print_level
     integer :: level
     integer :: error_nr
  end type qs_error_type
!!***
!****************************************************************************

contains
! ================== error structure init/get/set ===================

!!****f* qs_error_handling/qs_init_error [1.0] *
!!
!!   NAME
!!    qs_init_error
!!
!!   FUNCTION
!!    Initializes the error structure error.
!!    The optional arguments initialize the values in the structure.
!!    Default values: stop_level=2, logger=null (i.e. default logger
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: the error that will be initialized
!!     stop_level: the stop level of the newly created error (defaults to
!!           qs_failure_level)
!!     print_level: the print level of the newly created error (defaults
!!           to qs_note_level)
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
!!   SOURCE
!***************************************************************************
  subroutine qs_init_error(error, stop_level, print_level, logger, template_error, error_nr, level)
    type(qs_error_type), intent(out) :: error
    integer, optional, intent(in) :: stop_level, print_level, &
         error_nr, level
    type(qs_logger), optional, target, intent(in) :: logger
    type(qs_error_type), intent(in), optional :: template_error

    error%initialized=.TRUE.
    error%error_nr=0
    error%level=0
    error%print_level=0
    error%stop_level=2
    error%print_level=0
    nullify(error%logger)
    if (present(template_error)) then
       error%logger => template_error%logger
       error%stop_level = template_error%stop_level
       error%print_level = template_error%print_level
    end if
    if (present(logger)) then
       error%logger => logger
    end if
    if(present(stop_level)) then
       error%stop_level=stop_level
    end if
    if (present(print_level)) then
       error%print_level=print_level
    end if
    if (present(level)) error%level=level
    if(present(error_nr)) error%error_nr=error_nr
  end subroutine qs_init_error
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_reset_error [1.0] *
!!
!!   NAME
!!     qs_reset_error
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
!!   SOURCE
!***************************************************************************
  subroutine qs_reset_error(error)
    type(qs_error_type), optional, intent(inout) :: error

    if (present(error)) then
       if (.NOT. error%initialized) then
          call qs_log(qs_local_logger, qs_failure_level,  &
               'qs_error_handling:qs_reset_error:1', &
               'error structure not initialized')
          call qs_common_error_stop('qs_reset_error',&
               'error structure not initialized')
       end if
       error%error_nr=0
       error%level=0
    end if
  end subroutine qs_reset_error
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_set_error [1.0] *
!!
!!   NAME
!!     qs_set_error
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
!!   SOURCE
!***************************************************************************
  subroutine qs_set_error(error,level,error_nr,stop_level,print_level,logger,&
       template_error)
    type(qs_error_type), intent(inout) :: error
    integer, optional, intent(in) :: level, error_nr, stop_level,print_level
    type(qs_logger), optional, target, intent(in):: logger
    type(qs_error_type), optional, intent(in) :: template_error

    if (present(template_error)) then
       error%logger => template_error%logger
       error%stop_level = template_error%stop_level
       error%print_level = template_error%print_level
    end if
    if (present(level)) error%level=level
    if (present(error_nr)) error%error_nr=error_nr
    if (present(stop_level)) error%stop_level=stop_level
    if (present(print_level)) error%print_level=print_level
    if (present(logger)) error%logger=>logger
  end subroutine qs_set_error
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_dealloc_error_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_error_ref
!!
!!   FUNCTION
!!     Deallocates the memory allocated by the error
!!
!!   NOTES
!!     actually the error allocates no memeory, and it is better if it 
!!     stays so even in the future, but I have written this routine 
!!     to keep the consistence of the generic types
!!
!!   INPUTS
!!     err_to_free: the error that gets deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_error_ref(err_to_free, error)
    type(qs_error_type), intent(inout) :: err_to_free
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_error_handling:qs_dealloc_error_ref'
    failure=.false.
    call qs_assert(err_to_free%initialized,qs_warning_level,&
         qs_precondition_failed, routineN, "uninitialized error",&
         error,failure)
    if (.not. failure) then
! do nothing
    end if
  end subroutine qs_dealloc_error_ref
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_error_get [1.0] *
!!
!!   NAME
!!     qs_error_get
!!
!!   FUNCTION
!!     returns the various attributes of the given error
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_error_get(error,level,error_nr,stop_level,print_level,logger)
    type(qs_error_type), intent(inout) :: error
    integer, optional, intent(out) :: level, error_nr, stop_level,print_level
    type(qs_logger), optional, pointer :: logger

    logical :: failure
    character(len=*), parameter :: routineN='qs_error_handling:qs_error_get'

    if (present(level)) level=error%level
    if (present(error_nr)) error_nr=error%error_nr
    if (present(stop_level)) stop_level=error%stop_level
    if (present(print_level)) print_level=error%print_level
    if (present(logger)) logger=>qs_get_error_logger(error)
  end subroutine qs_error_get
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_get_error_level [1.0] *
!!
!!   NAME
!!     qs_get_error_level
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
!!   SOURCE
!***************************************************************************
  function qs_get_error_level(error) result(res)
    integer :: res
    type(qs_error_type), intent(in) :: error

    res=error%level
  end function qs_get_error_level
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_get_error_print_level [1.0] *
!!
!!   NAME
!!     qs_get_error_print_level
!!
!!   FUNCTION
!!     gets the print level (starting with which errorlevel something 
!!     gets printed) from the given error
!!     if the error is not present returns qs_note_level
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
!!   SOURCE
!***************************************************************************
  ! pure 
  function qs_get_error_print_level(error) result(res)
    integer :: res
    type(qs_error_type), optional, intent(in) :: error

    if (present(error)) then
       res=error%print_level
    else
       res=qs_note_level
    end if
  end function qs_get_error_print_level
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_get_error_nr [1.0] *
!!
!!   NAME
!!     qs_get_error_nr
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
!!   SOURCE
!***************************************************************************
  ! pure 
  function qs_get_error_nr(error) result(res)
    integer :: res
    type(qs_error_type), intent(in) :: error

    res=error%error_nr
  end function qs_get_error_nr
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_get_error_logger [1.0] *
!!
!!   NAME
!!     qs_get_error_logger
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
!!   SOURCE
!***************************************************************************
  function qs_get_error_logger(error) result(res)
    type(qs_logger), pointer :: res
    type(qs_error_type), intent(in), optional :: error

    nullify(res)
    if (present(error)) res=>error%logger
    if (.not. associated(res)) res=>qs_default_logger
    if (.not.associated(res)) then
       print *,'ERROR in qs_error_handling:qs_get_error_logger, qs_default_logger is not associated'
       call qs_common_error_stop('qs_error_handling:qs_get_error_logger',' qs_default_logger is not associated')
    end if
  end function qs_get_error_logger
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_get_error_stop_level [1.0] *
!!
!!   NAME
!!     qs_get_error_stop_level
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
!!   SOURCE
!***************************************************************************
  ! pure 
  function qs_get_error_stop_level(error) result(res)
    integer :: res
    type(qs_error_type), intent(in), optional :: error

    if (present(error)) then
       res=error%stop_level
    else
       res=qs_failure_level
    end if
  end function qs_get_error_stop_level

! =================== error handling ===================

!!****f* qs_error_handling/qs_common_error_stop [1.0] *
!!
!!   NAME
!!     qs_common_error_stop
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
!!   SOURCE
!***************************************************************************
  subroutine qs_common_error_stop(fromWhere,message)
    character(len=*), intent(in) :: fromWhere, message 
    call stop_program(fromWhere,message)
    stop ! just to be sure...
  end subroutine qs_common_error_stop
!!***
!****************************************************************************

!!****f* qs_error_handling/qs_handle_error [1.0] *
!!
!!   NAME
!!     qs_handle_error
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
!!   SOURCE
!***************************************************************************
  subroutine qs_handle_error(error, fromWhere, message)
    type(qs_error_type), intent(in) :: error
    character(len=*),intent(in) :: fromWhere,message

    if (debug_this_module) print *, "entred in qs_handle_error"
    call qs_error_message(error%level, fromWhere, &
         'err='//adjustl(qs_to_string(error%error_nr))//message, error)

    if (error%level >= error%stop_level) then
       call qs_common_error_stop(fromWhere,message)
    end if
  end subroutine qs_handle_error
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_error_message [1.0] *
!!
!!   NAME
!!     qs_error_message
!!
!!   FUNCTION
!!     writes an error message of level level, subject to control of error
!!     (i.e. if error%printlevel is >= level)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     level: the level of the error (qs_note_level, qs_warning_level,
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
!!   SOURCE
!***************************************************************************
  subroutine qs_error_message(level, fromWhere, message, error)
    integer, intent(in) :: level
    character(len=*), intent(in) :: fromWhere, message
    type(qs_error_type), optional, intent(in) :: error

    type(qs_logger), pointer :: logger

    logger => qs_get_error_logger(error)
    if (qs_get_error_print_level(error) <= level) then
       call qs_log(logger, level, fromWhere, message)
    end if
  end subroutine qs_error_message
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_assert [1.0] *
!!
!!   NAME
!!     qs_assert
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
!!     level: the level of the error (qs_note_level, qs_warning_level,
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
!!   SOURCE
!***************************************************************************
  subroutine qs_assert(condition, level, error_nr, fromWhere, &
       message, error,failure)
!! the condition that is checked, if false, an error did happen
    logical, intent(in) :: condition
!! level of the error
    integer, intent(in) :: level
!! error number
    integer, intent(in) :: error_nr
    character(len=*), intent(in) :: fromWhere, message
    type(qs_error_type), optional, intent(inout) :: error
!! is set to true if condition is false
    logical, optional, intent(inout) :: failure

    logical :: oldFailure
    type(qs_error_type) :: err
    character(len=*),parameter :: routineN='qs_error_handling:qs_assert'
    if (.NOT. condition) then
       if (debug_this_module) print *, "assertion failed, in "//fromWhere//" msg="//message
       if (present(failure)) then
          oldFailure=failure
          failure=.TRUE.
       end if
       call qs_init_error(err, level=level,error_nr=error_nr, &
            template_error=error)
       if (present(error)) then
          if (error%level <= level) then
             error%level=level
             error%error_nr=error_nr
          else if (present(failure) .AND. .NOT. oldFailure) then
             call qs_error_message(qs_warning_level,routineN,&
                  'error level bigger than level, but no failure appened before')
          end if
       end if
       call qs_handle_error(err, fromWhere,message)
    else
       if (debug_this_module) print *, "assertion checked, in "//fromWhere//" msg="//message
    end if
  end subroutine qs_assert
!!***
!***************************************************************************

!!****f* qs_error_handling/qs_propagate_error [1.0] *
!!
!!   NAME
!!     qs_propagate_error
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
!!           (defaults to qs_warning_level
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_propagate_error(iError, fromWhere, message, error, &
       failure, failure_level)
    type(qs_error_type), intent(in) :: iError
    character(len=*) :: fromWhere
    character(len=*), optional :: message
    type(qs_error_type), optional, intent(inout) :: error
    logical, optional, intent(inout):: failure
    integer, optional, intent(in):: failure_level

    character(len=*), parameter :: routineN='qs_error_handling:qs_propagate_error'
    integer :: f_level
    f_level=qs_warning_level
    if (present(failure_level)) f_level=failure_level
    if (qs_get_error_level(iError) /= 0 .and. iError%level>=f_level) then
       if (present(failure)) &
            failure=.true.
       if (present(error)) then
          if (error%level < iError%level) then
             error%level=iError%level
             error%error_nr=iError%error_nr
          end if
       end if
       if (iError%level >= qs_get_error_print_level(error)) then
          if (present(message)) then
             call qs_error_message(iError%level, fromWhere,"propagated error "//&
                  message, error)
          else
             call qs_error_message(iError%level, fromWhere,"propagated error",error)
          end if
       end if
    end if
  end subroutine qs_propagate_error
!!***
!***************************************************************************

end module qs_error_handling

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
!!     There are the following error levels: 0 note, 1 warning, 2 failure,
!!     3 fatal. You should use the constants defined here, not directly 
!!     the numbers
!!     @see cp_log_handling
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
module cp_error_handling
  use cp_log_handling
  use termination, only: stop_program
  implicit none
  private

  logical, private, parameter :: debug_this_module=.false.
  character(len=*), private, parameter :: moduleN='cp_error_handling'

  ! types
  public :: cp_error_type
  ! parameters var
  public :: cp_debug, cp_no_error, cp_caller_error, cp_wrong_args_error,&
       cp_precondition_failed, cp_internal_error, cp_postcondition_failed,&
       cp_invariant_failed, cp_assertion_failed, cp_unimplemented_error_nr
  ! public routines
  public :: cp_assert, cp_unimplemented_error
  ! generic routines
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_reset, cp_get_level,&
       cp_get_print_level, cp_get_nr, cp_get_logger, cp_get_stop_level,&
       cp_common_stop, cp_handle_error, cp_message, cp_propagate_error
  ! underlying procedures
  public :: cp_error_init, cp_error_dealloc_ref, cp_error_set,&
       cp_error_get, cp_error_reset, cp_error_get_level,&
       cp_error_get_print_level, cp_error_get_nr,&
       cp_error_get_logger, cp_error_get_stop_level,&
       cp_error_common_stop, cp_error_handle_error,&
       cp_error_message, cp_error_propagate_error

#ifdef FD_DEBUG
!! global variable to turn on/off the debugging, as i am lazy I use the 
!! preprocessor variable FD_DEBUG to decide if it should be on or off.
!! this way it is always right (on in debug builds, off in optimized ones)
  logical, parameter :: cp_debug=.true.
#else
  logical, parameter :: cp_debug=.false.
#endif

!! error number: no error
  integer, parameter :: cp_no_error = 0
!! error number: generic error on the side of the caller
  integer, parameter :: cp_caller_error = 1
!! error number: one or more arguments have and invalid value
  integer, parameter :: cp_wrong_args_error = 100
!! error number: precondition failed
  integer, parameter :: cp_precondition_failed = 200
!! error number: generic error inside the routine
  integer, parameter :: cp_internal_error = -1
!! error number: postcondition failed
  integer, parameter :: cp_postcondition_failed = -200
!! error number: invariant failed
  integer, parameter :: cp_invariant_failed = -100
!! error number: assertion failure
  integer, parameter :: cp_assertion_failed = -300
!! error number: not implemented
  integer, parameter :: cp_unimplemented_error_nr = -1000

  interface cp_init
     module procedure cp_error_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_error_dealloc_ref
  end interface
  interface cp_set
     module procedure cp_error_set
  end interface
  interface cp_get
     module procedure cp_error_get
  end interface
  interface cp_reset
     module procedure cp_error_reset
  end interface
  interface cp_get_level
     module procedure cp_error_get_level
  end interface
  interface cp_get_print_level
     module procedure cp_error_get_print_level
  end interface
  interface cp_get_nr
     module procedure cp_error_get_nr
  end interface
  interface cp_get_logger
     module procedure cp_error_get_logger
  end interface
  interface cp_get_stop_level
     module procedure cp_error_get_stop_level
  end interface
  interface cp_common_stop
     module procedure cp_error_common_stop
  end interface
  interface cp_handle_error
     module procedure cp_error_handle_error
  end interface
  interface cp_message
     module procedure cp_error_message
  end interface
  interface cp_propagate_error
     module procedure cp_error_propagate_error
  end interface
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
  type cp_error_type
     private
     logical :: initialized
     type(cp_logger_type), pointer :: logger
     integer :: stop_level
     integer :: print_level
     integer :: level
     integer :: error_nr
  end type cp_error_type
!!***
!****************************************************************************

contains
! ================== error structure init/get/set ===================

!!****f* cp_error_handling/cp_error_init [1.0] *
!!
!!   NAME
!!    cp_error_init
!!
!!   SYNOPSIS
!!     Subroutine cp_error_init(error, stop_level, print_level, logger,&
!!         template_error, error_nr, level)
!!       Type(cp_error_type), Intent (OUT):: error
!!       Integer, Optional, Intent (IN):: stop_level, print_level, error_nr,&
!!         level
!!       Type(cp_logger_type), Optional, Target, Intent (IN):: logger
!!       Type(cp_error_type), Intent (IN), Optional:: template_error
!!     End Subroutine cp_error_init
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
  subroutine cp_error_init(error, stop_level, print_level, logger, template_error, error_nr, level)
    type(cp_error_type), intent(out) :: error
    integer, optional, intent(in) :: stop_level, print_level, &
         error_nr, level
    type(cp_logger_type), optional, target, intent(in) :: logger
    type(cp_error_type), intent(in), optional :: template_error

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
  end subroutine cp_error_init
!****************************************************************************

!!****f* cp_error_handling/cp_error_reset [1.0] *
!!
!!   NAME
!!     cp_error_reset
!!
!!   SYNOPSIS
!!     Subroutine cp_error_reset(error)
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_error_reset
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
  subroutine cp_error_reset(error)
    type(cp_error_type), optional, intent(inout) :: error

    if (present(error)) then
       if (.NOT. error%initialized) then
          call cp_log(cp_get_default_logger(), cp_failure_level,  &
               'cp_error_handling:cp_error_reset:1', &
               'error structure not initialized')
          call cp_error_common_stop('cp_error_reset',&
               'error structure not initialized')
       end if
       error%error_nr=0
       error%level=0
    end if
  end subroutine cp_error_reset
!****************************************************************************

!!****f* cp_error_handling/cp_error_set [1.0] *
!!
!!   NAME
!!     cp_error_set
!!
!!   SYNOPSIS
!!     Subroutine cp_error_set(error, level, error_nr, stop_level,&
!!         print_level, logger, template_error)
!!       Type(cp_error_type), Intent (INOUT):: error
!!       Integer, Optional, Intent (IN):: level, error_nr, stop_level,&
!!         print_level
!!       Type(cp_logger_type), Optional, Target, Intent (IN):: logger
!!       Type(cp_error_type), Optional, Intent (IN):: template_error
!!     End Subroutine cp_error_set
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
  subroutine cp_error_set(error,level,error_nr,stop_level,print_level,logger,&
       template_error)
    type(cp_error_type), intent(inout) :: error
    integer, optional, intent(in) :: level, error_nr, stop_level,print_level
    type(cp_logger_type), optional, target, intent(in):: logger
    type(cp_error_type), optional, intent(in) :: template_error

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
  end subroutine cp_error_set
!***************************************************************************

!!****f* cp_error_handling/cp_error_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_error_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_error_dealloc_ref(err_to_free, error)
!!       Type(cp_error_type), Intent (INOUT):: err_to_free
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_error_dealloc_ref
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
  subroutine cp_error_dealloc_ref(err_to_free, error)
    type(cp_error_type), intent(inout) :: err_to_free
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_error_handling:cp_error_dealloc_ref'
    failure=.false.
    call cp_assert(err_to_free%initialized,cp_warning_level,&
         cp_precondition_failed, routineP, "uninitialized error",&
         error,failure)
    if (.not. failure) then
! do nothing
    end if
  end subroutine cp_error_dealloc_ref
!***************************************************************************

!!****f* cp_error_handling/cp_error_get [1.0] *
!!
!!   NAME
!!     cp_error_get
!!
!!   SYNOPSIS
!!     Subroutine cp_error_get(error, level, error_nr, stop_level,&
!!         print_level, logger)
!!       Type(cp_error_type), Intent (INOUT):: error
!!       Integer, Optional, Intent (OUT):: level, error_nr, stop_level,&
!!         print_level
!!       Type(cp_logger_type), Optional, Pointer:: logger
!!     End Subroutine cp_error_get
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
  subroutine cp_error_get(error,level,error_nr,stop_level,print_level,logger)
    type(cp_error_type), intent(inout) :: error
    integer, optional, intent(out) :: level, error_nr, stop_level,print_level
    type(cp_logger_type), optional, pointer :: logger

    logical :: failure
    character(len=*), parameter :: routineP='cp_error_handling:cp_error_get'

    if (present(level)) level=error%level
    if (present(error_nr)) error_nr=error%error_nr
    if (present(stop_level)) stop_level=error%stop_level
    if (present(print_level)) print_level=error%print_level
    if (present(logger)) logger=>cp_error_get_logger(error)
  end subroutine cp_error_get
!***************************************************************************

!!****f* cp_error_handling/cp_error_get_level [1.0] *
!!
!!   NAME
!!     cp_error_get_level
!!
!!   SYNOPSIS
!!     Function cp_error_get_level(error) Result(res)
!!       Integer:: res
!!       Type(cp_error_type), Intent (IN):: error
!!     End Function cp_error_get_level
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
  function cp_error_get_level(error) result(res)
    integer :: res
    type(cp_error_type), intent(in) :: error

    res=error%level
  end function cp_error_get_level
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_print_level [1.0] *
!!
!!   NAME
!!     cp_error_get_print_level
!!
!!   SYNOPSIS
!!     Function cp_error_get_print_level(error) Result(res)
!!       Integer:: res
!!       Type(cp_error_type), Optional, Intent (IN):: error
!!     End Function cp_error_get_print_level
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
  function cp_error_get_print_level(error) result(res)
    integer :: res
    type(cp_error_type), optional, intent(in) :: error

    if (present(error)) then
       res=error%print_level
    else
       res=cp_note_level
    end if
  end function cp_error_get_print_level
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_nr [1.0] *
!!
!!   NAME
!!     cp_error_get_nr
!!
!!   SYNOPSIS
!!     Function cp_error_get_nr(error) Result(res)
!!       Integer:: res
!!       Type(cp_error_type), Intent (IN):: error
!!     End Function cp_error_get_nr
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
  function cp_error_get_nr(error) result(res)
    integer :: res
    type(cp_error_type), intent(in) :: error

    res=error%error_nr
  end function cp_error_get_nr
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_logger [1.0] *
!!
!!   NAME
!!     cp_error_get_logger
!!
!!   SYNOPSIS
!!     Function cp_error_get_logger(error) Result(res)
!!       Type(cp_logger_type), Pointer:: res
!!       Type(cp_error_type), Intent (IN), Optional:: error
!!     End Function cp_error_get_logger
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
  function cp_error_get_logger(error) result(res)
    type(cp_logger_type), pointer :: res
    type(cp_error_type), intent(in), optional :: error

    nullify(res)
    if (present(error)) res=>error%logger
    if (.not. associated(res)) res=>cp_get_default_logger()
    if (.not.associated(res)) then
       print *,'ERROR in cp_error_handling:cp_error_get_logger, cp_default_logger is not associated'
       call cp_error_common_stop('cp_error_handling:cp_error_get_logger',' cp_default_logger is not associated')
    end if
  end function cp_error_get_logger
!****************************************************************************

!!****f* cp_error_handling/cp_error_get_stop_level [1.0] *
!!
!!   NAME
!!     cp_error_get_stop_level
!!
!!   SYNOPSIS
!!     Function cp_error_get_stop_level(error) Result(res)
!!       Integer:: res
!!       Type(cp_error_type), Intent (IN), Optional:: error
!!     End Function cp_error_get_stop_level
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
  function cp_error_get_stop_level(error) result(res)
    integer :: res
    type(cp_error_type), intent(in), optional :: error

    if (present(error)) then
       res=error%stop_level
    else
       res=cp_failure_level
    end if
  end function cp_error_get_stop_level

! =================== error handling ===================

!!****f* cp_error_handling/cp_error_common_stop [1.0] *
!!
!!   NAME
!!     cp_error_common_stop
!!
!!   SYNOPSIS
!!     Subroutine cp_error_common_stop(fromwhere, message)
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!     End Subroutine cp_error_common_stop
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
  subroutine cp_error_common_stop(fromWhere,message)
    character(len=*), intent(in) :: fromWhere, message 
    call stop_program(fromWhere,message)
    stop ! just to be sure...
  end subroutine cp_error_common_stop
!****************************************************************************

!!****f* cp_error_handling/cp_error_handle_error [1.0] *
!!
!!   NAME
!!     cp_error_handle_error
!!
!!   SYNOPSIS
!!     Subroutine cp_error_handle_error(error, fromwhere, message)
!!       Type(cp_error_type), Intent (IN):: error
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!     End Subroutine cp_error_handle_error
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
  subroutine cp_error_handle_error(error, fromWhere, message)
    type(cp_error_type), intent(in) :: error
    character(len=*),intent(in) :: fromWhere,message

    integer :: mepos
    type(cp_logger_type), pointer :: logger


    if (debug_this_module) print *, "entred in cp_error_handle_error"
    call cp_error_message(error%level, fromWhere, &
         'err='//adjustl(cp_to_string(error%error_nr))//message, error)

    if (error%level >= error%stop_level) then
       logger => cp_error_get_logger(error) 
       mepos=logger%para_env%mepos
       call cp_error_common_stop(fromWhere//' processor '//&
            cp_to_string(mepos)//' ',message)
    end if
  end subroutine cp_error_handle_error
!***************************************************************************

!!****f* cp_error_handling/cp_error_message [1.0] *
!!
!!   NAME
!!     cp_error_message
!!
!!   SYNOPSIS
!!     Subroutine cp_error_message(level, fromwhere, message, error)
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!       Type(cp_error_type), Optional, Intent (IN):: error
!!     End Subroutine cp_error_message
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
  subroutine cp_error_message(level, fromWhere, message, error)
    integer, intent(in) :: level
    character(len=*), intent(in) :: fromWhere, message
    type(cp_error_type), optional, intent(in) :: error

    type(cp_logger_type), pointer :: logger

    logger => cp_error_get_logger(error)
    if (cp_error_get_print_level(error) <= level) then
       call cp_log(logger, level, fromWhere, message,local=.true.)
    end if
  end subroutine cp_error_message
!***************************************************************************

!!****f* cp_error_handling/cp_assert [1.0] *
!!
!!   NAME
!!     cp_assert
!!
!!   SYNOPSIS
!!     Subroutine cp_assert(condition, level, error_nr, fromwhere, message,&
!!         error, failure)
!!       Logical, Intent (IN):: condition
!!       Integer, Intent (IN):: level
!!       Integer, Intent (IN):: error_nr
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical, Optional, Intent (INOUT):: failure
!!     End Subroutine cp_assert
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
  subroutine cp_assert(condition, level, error_nr, fromWhere, &
       message, error,failure)
!! the condition that is checked, if false, an error did happen
    logical, intent(in) :: condition
!! level of the error
    integer, intent(in) :: level
!! error number
    integer, intent(in) :: error_nr
    character(len=*), intent(in) :: fromWhere, message
    type(cp_error_type), optional, intent(inout) :: error
!! is set to true if condition is false
    logical, optional, intent(inout) :: failure

    logical :: oldFailure
    type(cp_error_type) :: err
    character(len=*),parameter :: routineP='cp_error_handling:cp_assert'
    if (.NOT. condition) then
       if (debug_this_module) print *, "assertion failed, in "//fromWhere//" msg="//message
       if (present(failure)) then
          oldFailure=failure
          failure=.TRUE.
       end if
       call cp_error_init(err, level=level,error_nr=error_nr, &
            template_error=error)
       if (present(error)) then
          if (error%level < level) then
             error%level=level
             error%error_nr=error_nr
          else if (present(failure) .AND. .NOT. oldFailure) then
             call cp_error_message(cp_warning_level,routineP,&
                  'error level bigger than level, but no failure appened before')
          end if
       end if
       call cp_error_handle_error(err, fromWhere,message)
       call cp_error_dealloc_ref(err)
    else
       if (debug_this_module) print *, "assertion checked, in "//fromWhere//" msg="//message
    end if
  end subroutine cp_assert
!***************************************************************************

!!****f* cp_error_handling/cp_error_propagate_error [1.0] *
!!
!!   NAME
!!     cp_error_propagate_error
!!
!!   SYNOPSIS
!!     Subroutine cp_error_propagate_error(ierror, fromwhere, message,&
!!         error, failure, failure_level)
!!       Type(cp_error_type), Intent (IN):: ierror
!!       Character(Len=*):: fromwhere
!!       Character(Len=*), Optional:: message
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical, Optional, Intent (INOUT):: failure
!!       Integer, Optional, Intent (IN):: failure_level
!!     End Subroutine cp_error_propagate_error
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
  subroutine cp_error_propagate_error(iError, fromWhere, message, error, &
       failure, failure_level)
    type(cp_error_type), intent(in) :: iError
    character(len=*) :: fromWhere
    character(len=*), optional :: message
    type(cp_error_type), optional, intent(inout) :: error
    logical, optional, intent(inout):: failure
    integer, optional, intent(in):: failure_level

    character(len=*), parameter :: routineP='cp_error_handling:cp_error_propagate_error'
    integer :: f_level
    f_level=cp_warning_level
    if (present(failure_level)) f_level=failure_level
    if (cp_error_get_level(iError) /= 0 .and. iError%level>=f_level) then
       if (present(failure)) &
            failure=.true.
       if (present(error)) then
          if (error%level < iError%level) then
             error%level=iError%level
             error%error_nr=iError%error_nr
          end if
       end if
       if (iError%level >= cp_error_get_print_level(error)) then
          if (present(message)) then
             call cp_error_message(iError%level, fromWhere,"propagated error "//&
                  message, error)
          else
             call cp_error_message(iError%level, fromWhere,"propagated error",error)
          end if
       end if
    end if
  end subroutine cp_error_propagate_error
!***************************************************************************

!!****f* cp_error_handling/cp_unimplemented_error [1.0] *
!!
!!   NAME
!!     cp_unimplemented_error
!!
!!   SYNOPSIS
!!     Subroutine cp_unimplemented_error(fromwhere, message, error,&
!!         error_level)
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Character(Len=*), Intent (IN), Optional:: message
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Integer, Intent (IN), Optional:: error_level
!!     End Subroutine cp_unimplemented_error
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
subroutine cp_unimplemented_error(fromWhere, message, error, error_level)
  character(len=*), intent(in) :: fromWhere
  character(len=*), intent(in), optional :: message
  type(cp_error_type), optional, intent(inout) :: error
  integer, intent(in),optional :: error_level
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_unimplemented_error',&
        routineP=moduleN//':'//routineN
  integer :: err_l
  err_l=cp_failure_level

  if(present(error_level)) err_l=error_level
  if (present(message)) then
     call cp_assert(.false.,err_l,cp_unimplemented_error_nr,&
          fromWhere, "UNIMPLEMENTED, "//message, error=error)
  else
     call cp_assert(.false.,err_l,cp_unimplemented_error_nr,&
          fromWhere, "UNIMPLEMENTED", error=error)
  end if
end subroutine cp_unimplemented_error
!***************************************************************************

end module cp_error_handling

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_log_handling [1.0] *
!!
!!   NAME
!!     cp_log_handling
!!
!!   FUNCTION
!!     various routines to log and control the output.
!!     The idea is that decisions about where to log should not be done in
!!     the code that generates the log, but should be globally changeable
!!     a central place.
!!     So some care has been taken to have enough information about the 
!!     place from where the log comes so that in the future intelligent and
!!     flexible decisions can be taken by the logger, without having to change
!!     other code.
!!
!!   NOTES
!!     contains also routines to convert to a string.
!!     in my idea they should have been with variable length,
!!     (i.e. they should have returned a trim(adjustl(actual_result)))
!!     As a logger should be robust, at the moment I have given up.
!!
!!     At the moment logging and output refer to the same object 
!!     (cp_logger_type)
!!     as these are actually different it might be better to separate them
!!     (they have already separate routines in a separate module 
!!     @see cp_output_handling).
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!     @see cp_error_handling
!!
!!   SOURCE
!******************************************************************************

module cp_log_handling
  use output_utilities, only: print_message
  use qs_parser, only: open_file, close_file
  implicit none
  private

! types
  public :: cp_logger_type
! parameter vars
  public :: cp_note_level, cp_warning_level, cp_failure_level,&
       cp_fatal_level, cp_std_in_unit_nr, cp_std_out_unit_nr
! default loggers
  public :: cp_get_local_logger, cp_get_global_logger, cp_get_default_logger
! generic routines
  public :: cp_init, cp_dealloc_ref, cp_would_log, cp_get_unit_nr,&
       cp_log, cp_set_log_level, cp_get_default_unit_nr
  public :: cp_to_string
! underlying routines
  public :: cp_logger_init, cp_logger_dealloc_ref, cp_logger_would_log,&
       cp_logger_get_unit_nr, cp_logger_log, cp_logger_set_log_level,&
       cp_logger_get_default_unit_nr
  public :: cp_int_to_string, cp_real_sp_to_string, &
       cp_real_dp_to_string

  character(len=*), private, parameter :: moduleN='cp_log_handling'
  logical, private, parameter :: debug_this_module=.true.
  integer, parameter :: cp_std_in_unit_nr=5
  integer, parameter :: cp_std_out_unit_nr=6

!! level of an error
  integer, parameter  :: cp_fatal_level=3
!! level of a failure
  integer, parameter  :: cp_failure_level=2
!! level of a warning
  integer, parameter  :: cp_warning_level=1
!! level of a note
  integer, parameter  :: cp_note_level=0

!! a generic function to trasform different types to strings
  interface cp_to_string
     module procedure cp_int_to_string, cp_real_sp_to_string, &
          cp_real_dp_to_string
  end interface
  interface cp_init
     module procedure cp_logger_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_logger_dealloc_ref
  end interface
  interface cp_would_log
     module procedure cp_logger_would_log
  end interface
  interface cp_get_unit_nr
     module procedure cp_logger_get_unit_nr
  end interface
  interface cp_log
     module procedure cp_logger_log
  end interface
  interface cp_set_log_level
     module procedure cp_logger_set_log_level
  end interface
  interface cp_get_default_unit_nr
     module procedure cp_logger_get_default_unit_nr
  end interface
!!***
!****************************************************************************

!!****s* cp_log_handling/cp_logger_type [1.0] *
!!
!!   NAME
!!     cp_logger_type
!!
!!   FUNCTION
!!     type of a logger, at the moment it contains just a print level
!!     starting at which level it should be logged 
!!     (0 note, 1 warning, 2 failure, 3 fatal)
!!     it could be expanded with the ability to focus on one or more
!!     module/object/thread/processor
!!
!!   NOTES
!!     This should be private, but as the output functions have been 
!!     moved to another module and there is no "friend" keyword, they
!!     are public.
!!     DO NOT USE THE INTERNAL COMPONENTS DIRECTLY!!!
!!
!!   ATTRIBUTES
!!     print_level: the level starting at which something gets printed
!!                 (defaults to 0)
!!     default_unit_nr: the deafult unit to which evrything is written
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!   SOURCE
!***************************************************************************
  type cp_logger_type
     integer :: print_level
     integer :: default_unit_nr
     logical :: close_unit_on_dealloc
  end type cp_logger_type
!!***
!****************************************************************************

! the default local logger
  type(cp_logger_type), target ::cp_local_logger=cp_logger_type(0,-1,.false.)
! the default global logger
  type(cp_logger_type), target ::cp_global_logger=cp_logger_type(0,-2,.false.)
! the default logger (normally points to the local logger)
  type(cp_logger_type), pointer ::cp_default_logger ! => null()

contains

!!****f* cp_log_handling/cp_init_default_loggers [1.0] *
!!
!!   NAME
!!     cp_init_default_loggers
!!
!!   FUNCTION
!!     performs the initializations of the logger.
!!     MUST be called before logging occours
!!
!!   NOTES
!!     initializes only if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  subroutine cp_init_default_loggers
    logical,save ::firstTime=.true.
    if (firstTime) then
       firstTime=.false.
       call cp_logger_init(cp_local_logger)
       call cp_logger_init(cp_global_logger)
       cp_default_logger => cp_local_logger
    end if
  end subroutine cp_init_default_loggers
!**************************************************************************

!!****f* cp_log_handling/cp_get_default_logger [1.0] *
!!
!!   NAME
!!     cp_get_default_logger
!!
!!   SYNOPSIS
!!     Function cp_get_default_logger() Result(res)
!!       Type(cp_logger_type), Pointer:: res
!!     End Function cp_get_default_logger
!!
!!   FUNCTION
!!     returns the default logger
!!
!!   NOTES
!!     initializes the default loggers if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
function cp_get_default_logger() result(res)
  type(cp_logger_type), pointer ::res
  
  call cp_init_default_loggers
  res => cp_default_logger
end function cp_get_default_logger
!***************************************************************************

!!****f* cp_log_handling/cp_get_local_logger [1.0] *
!!
!!   NAME
!!     cp_get_local_logger
!!
!!   SYNOPSIS
!!     Function cp_get_local_logger() Result(res)
!!       Type(cp_logger_type), Pointer:: res
!!     End Function cp_get_local_logger
!!
!!   FUNCTION
!!     returns the local logger
!!
!!   NOTES
!!     initializes the default loggers if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
function cp_get_local_logger() result(res)
  type(cp_logger_type), pointer ::res
  
  call cp_init_default_loggers
  res => cp_local_logger
end function cp_get_local_logger
!***************************************************************************

!!****f* cp_log_handling/cp_get_global_logger [1.0] *
!!
!!   NAME
!!     cp_get_global_logger
!!
!!   SYNOPSIS
!!     Function cp_get_global_logger() Result(res)
!!       Type(cp_logger_type), Pointer:: res
!!     End Function cp_get_global_logger
!!
!!   FUNCTION
!!     returns the global logger
!!
!!   NOTES
!!     initializes the default loggers if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
function cp_get_global_logger() result(res)
  type(cp_logger_type), pointer ::res
  
  call cp_init_default_loggers
  res => cp_global_logger
end function cp_get_global_logger
!***************************************************************************

! ================== log ==================

!!****f* cp_log_handling/cp_logger_init [1.0] *
!!
!!   NAME
!!     cp_logger_init
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_init(logger, print_level, default_unit_nr,&
!!         filename, close_unit_on_dealloc)
!!       Type(cp_logger_type), Intent (OUT):: logger
!!       Integer, Intent (IN), Optional:: print_level, default_unit_nr
!!       Character(Len=*), Intent (IN), Optional:: filename
!!       Logical, Intent (IN), Optional:: close_unit_on_dealloc
!!     End Subroutine cp_logger_init
!!
!!   FUNCTION
!!     initializes a logger
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger to initialize
!!     - print_level: the level starting with which something is written
!!       (defaults to cp_note_level)
!!     - default_unit_nr: the default unit_nr for output (if not given, and 
!!       no file is given defaults to the standard output)
!!     - filename: a new file to open (can be given instread of the 
!!       default_unit_nr)
!!     - close_unit_on_dealloc: if the unit should be closed when the logger 
!!       is deallocated (defaults to false if a default_unit_nr is given, to
!!       true if a filename is given)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_init(logger, print_level, default_unit_nr, filename, &
       close_unit_on_dealloc)
    type(cp_logger_type), intent(out) :: logger
    integer, intent(in), optional :: print_level, default_unit_nr
    character(len=*), intent(in), optional :: filename
    logical, intent(in), optional :: close_unit_on_dealloc

    character(len=*), parameter :: routineN='cp_logger_init',&
         routineP=moduleN//':'//routineN
    logical :: exists

    logger%default_unit_nr=cp_std_out_unit_nr
    logger%print_level=cp_note_level
    logger%close_unit_on_dealloc=.false.
    if (present(print_level)) logger%print_level=print_level
    if (present(default_unit_nr)) logger%default_unit_nr=default_unit_nr
    if (present(filename)) then
       if (present(default_unit_nr)) then
          print *," WARNING ", routineP, " cannot initialize logger with both filename and default_unit_nr, default_unit_nr ignored"
       end if
       logger%close_unit_on_dealloc=.true.
       INQUIRE (FILE=filename,EXIST=exists)
       if (exists) then
          print *," WARNING ", routineP, " overwriting file ",filename
       end if
       call open_file(filename,file_status="unknown",&
            file_action="write",&
            unit_number=logger%default_unit_nr)
    end if
    if (present(close_unit_on_dealloc)) then
       logger%close_unit_on_dealloc=close_unit_on_dealloc
    end if
  end subroutine cp_logger_init
!***************************************************************************

!!****f* cp_logger_dealloc_ref *
!!
!!   NAME
!!     cp_logger_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_dealloc_ref(logger)
!!       Type(cp_logger_type), Intent (INOUT):: logger
!!     End Subroutine cp_logger_dealloc_ref
!!
!!   FUNCTION
!!     removes the references of this logger, after this it cannot be
!!     used anymore.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger to dealloc
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_dealloc_ref(logger)
    type(cp_logger_type), intent(inout) :: logger

    character(len=*), parameter :: routineN='cp_logger_type',&
         routineP=moduleN//':'//routineN

    if (logger%close_unit_on_dealloc) then
       call close_file(logger%default_unit_nr)
       logger%close_unit_on_dealloc=.false.
       logger%default_unit_nr=-1
    end if
  end subroutine cp_logger_dealloc_ref
!***************************************************************************

!!****f* cp_log_handling/cp_logger_would_log [1.0] *
!!
!!   NAME
!!     cp_logger_would_log
!!
!!   SYNOPSIS
!!     Function cp_logger_would_log(logger, level, fromwhere)
!!       Logical:: cp_logger_would_log
!!       Type(cp_logger_type), Intent (IN):: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!     End Function cp_logger_would_log
!!
!!   FUNCTION
!!     this function can be called to check if the logger would log
!!     a message with the given level from the given source
!!     you should use this function if you do direct logging
!!     (without using cp_logger_log), or if you want to know if the generation 
!!     of some costly log info is necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!  - logger: the logger you want to log in
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure
  function cp_logger_would_log(logger, level, fromWhere)
    logical :: cp_logger_would_log
    type(cp_logger_type), intent(in) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere

    cp_logger_would_log=level>=logger%print_level
  end function cp_logger_would_log

!!****f* cp_log_handling/cp_logger_get_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_unit_nr(logger, level, fromwhere)
!!       Integer:: cp_logger_get_unit_nr
!!       Type(cp_logger_type), Intent (IN):: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!     End Function cp_logger_get_unit_nr
!!
!!   FUNCTION
!!     returns the unit nr for the requested kind of log.
!!
!!     Try not to use this function, in future in some cases the log 
!!     could not use units (global log might be sent to the master process)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!  - logger: the logger you want to log in
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  function cp_logger_get_unit_nr(logger, level, fromWhere)
    integer :: cp_logger_get_unit_nr
    type(cp_logger_type), intent(in) :: logger
    integer, intent(in) :: level
    character(len=*), intent(in):: fromWhere

    if (logger%default_unit_nr<=0) then
       cp_logger_get_unit_nr=cp_std_out_unit_nr
    else
       cp_logger_get_unit_nr=logger%default_unit_nr
    end if
  end function cp_logger_get_unit_nr

!!****f* cp_log_handling/cp_logger_log [1.0] *
!!
!!   NAME
!!     cp_logger_log
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_log(logger, level, fromwhere, message)
!!       Type(cp_logger_type):: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!     End Subroutine cp_logger_log
!!
!!   FUNCTION
!!     Logs in a standard format
!!
!!   NOTES
!!     In multithreaded situations it should be modified to print also 
!!     the thread (and maybe the processor), and proper locking should be 
!!     performed.
!!     At the moment it just goes to the deafault unit of the logger
!!
!!   INPUTS
!!  - logger is a type that describes where the logging must go.
!!    There are three standard loggers: cp_local_logger, cp_global_logger,
!!    and cp_default_logger, that normally points to one of the previous two
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - message is the message to print
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  subroutine cp_logger_log(logger, level, fromWhere , message)
    type(cp_logger_type) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere, message

    character(7):: label
    integer ::oErr, unitNr
    if (cp_logger_would_log(logger,level,fromWhere)) then
       select case (level)
       case (cp_fatal_level)
          label='ERRORL3'
       case (cp_failure_level)
          label='ERRORL2'
       case (cp_warning_level)
          label='WARNING'
       case (cp_note_level)
          label='NOTE'
       case default
          label='L'//adjustl(cp_int_to_string(level))
       end select
       unitNr=cp_logger_get_unit_nr(logger,level,fromWhere)

       call print_message(label//' in '//fromWhere//', '//message,unitNr,&
            level, level,level/2)
! old (more failsafe, no trim, compress,...) code:
!write (unit=unitNr,fmt='(a)', iostat=oErr) ' '
!write (unit=unitNr,fmt='(5a)', iostat=oErr) label,' in ', fromWhere,',', message
!if (oErr /= 0) then
!   print *, 'ERROR  ',oErr,' while writing log, ', fromWhere, message
!end if
    end if
  end subroutine cp_logger_log

! *************************** cp_logger_type settings ***************************

!!****f* cp_log_handling/cp_logger_set_log_level [1.0] *
!!
!!   NAME
!!     cp_logger_set_log_level
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_set_log_level(logger, level)
!!       Type(cp_logger_type), Intent (OUT):: logger
!!       Integer, Intent (IN):: level
!!     End Subroutine cp_logger_set_log_level
!!
!!   FUNCTION
!!     changes the logging level. Log messages with a level less than the one
!!     given wo not be printed.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger to change
!!     - level : the new logging level for the logger
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  subroutine cp_logger_set_log_level(logger, level)
    type(cp_logger_type), intent(out) :: logger
    integer, intent(in) :: level
    logger%print_level=level
  end subroutine cp_logger_set_log_level

!!****f* cp_log_handling/cp_logger_get_log_level [1.0] *
!!
!!   NAME
!!     cp_logger_get_log_level
!!
!!   FUNCTION
!!     returns the logging level of the given logger.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger you want info from
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
!pure
  function cp_logger_get_log_level(logger)
    integer :: cp_logger_get_log_level
    type(cp_logger_type), intent(in) :: logger
    cp_logger_get_log_level=logger%print_level
  end function cp_logger_get_log_level

!!****f* cp_log_handling/cp_logger_get_default_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_default_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_default_unit_nr(logger)
!!       Integer:: cp_logger_get_default_unit_nr
!!       Type(cp_logger_type), Intent (IN):: logger
!!     End Function cp_logger_get_default_unit_nr
!!
!!   FUNCTION
!!     asks the unit number of the given logger.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger you want info from
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
!pure 
  function cp_logger_get_default_unit_nr(logger)
    integer :: cp_logger_get_default_unit_nr
    type(cp_logger_type), intent(in) :: logger
    cp_logger_get_default_unit_nr=logger%default_unit_nr
  end function cp_logger_get_default_unit_nr

! ============= to string  ============

!!****f* cp_log_handling/cp_int_to_string [1.0] *
!!
!!   NAME
!!     cp_int_to_string
!!
!!   SYNOPSIS
!!     Function cp_int_to_string(i) Result(res)
!!       Character(Len=6):: res
!!       Integer, Intent (IN):: i
!!     End Function cp_int_to_string
!!
!!   FUNCTION
!!     converts an int to a string
!!     (should be a variable length string, but that does not work with 
!!     all the compilers)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     i: the integer to convert
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure
  function cp_int_to_string(i) result(res)
    character(len=6) :: res
    integer, intent(in):: i
    write (res,'(i6)') i
    if (i>999999 .or. i<-99999) then
       write (res, '(es6.1)') real(i)
    end if
  end function cp_int_to_string

!!****f* cp_log_handling/cp_real_sp_to_string [1.0] *
!!
!!   NAME
!!     cp_real_sp_to_string
!!
!!   SYNOPSIS
!!     Function cp_real_sp_to_string(val) Result(res)
!!       Character(Len=10):: res
!!       Real(Kind=kind(1.0)), Intent (IN):: val
!!     End Function cp_real_sp_to_string
!!
!!   FUNCTION
!!     convert a single precision real in a string
!!     (should be a variable length string, but that does not work with 
!!     all the compilers)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     val: the number to convert
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  function cp_real_sp_to_string(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0)), intent(in):: val
    write (res,'(es10.4)') val
  end function cp_real_sp_to_string

!!****f* cp_log_handling/cp_real_dp_to_string [1.0] *
!!
!!   NAME
!!     cp_real_dp_to_string
!!
!!   SYNOPSIS
!!     Function cp_real_dp_to_string(val) Result(res)
!!       Character(Len=10):: res
!!       Real(Kind=kind(1.0D0)), Intent (IN):: val
!!     End Function cp_real_dp_to_string
!!
!!   FUNCTION
!!     convert a double precision real in a string
!!     (should be a variable length string, but that does not work with 
!!     all the compilers)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     val: the number to convert
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  function cp_real_dp_to_string(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0D0)), intent(in):: val
    write (res,'(es10.4)') val
  end function cp_real_dp_to_string

end module cp_log_handling

!!***
!******************************************************************************

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
!!     some practices (use of print *, stop_program, no cp_error_type,
!!     manual retain release of some objects) are dictated by the need to
!!     have minimal dependency
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     08.2002 major update: retain, release, printkeys, para_env, 
!!             local logging [fawzi]
!!
!!     @see cp_error_handling
!!
!!   SOURCE
!******************************************************************************

module cp_log_handling
  use output_utilities, only: print_message
  use qs_parser, only: open_file, close_file
  use machine, only: m_flush
  use message_passing, only: mp_comm_free, mp_bcast
  use termination, only: stop_program
  use cp_para_types, only: cp_para_env_type
  use global_types, only: print_key_type, init_print_keys
  use string_utilities, only: compress_all
  use timesl, only: walltime
  use kinds, only: wp=>dp
  implicit none
  private

! types
  public :: cp_logger_type
! parameter vars
  public :: cp_note_level, cp_warning_level, cp_failure_level,&
       cp_fatal_level, cp_std_in_unit_nr, cp_std_out_unit_nr
! default loggers
  public :: cp_get_default_logger, cp_init_default_logger
! logger routines
  public :: cp_logger_create, cp_logger_retain, cp_logger_release,&
       cp_logger_would_log,&
       cp_logger_get_default_unit_nr, cp_logger_get_unit_nr, &
       cp_logger_set_log_level, cp_logger_generate_filename
  public :: cp_to_string, cp_log
! underlying routines
!!FM  public :: cp_int_to_string, cp_real_sp_to_string, &
!!FM       cp_real_dp_to_string

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
     module procedure cp_int_to_string, cp_real_dp_to_string
!!FM          cp_real_sp_to_string
  end interface
  interface cp_log
     module procedure cp_logger_log
  end interface

  logical, save, private :: default_logger_initialized=.false.
  integer, save, private :: last_logger_id_nr=0
  real(kind=wp), save, private :: start_time

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
!!     - to do
!!     - print_level: the level starting at which something gets printed
!!     - para_env: the parallel environement for the output.
!!       this might be a super environment of your computation environment
!!       i.e. be very careful not to do global operations like broadcast
!!       with a subset of its processors (use your computation environment
!!       instread
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2002 revised [fawzi]
!!     08.2002 major update: retain, release, printkeys, para_env, 
!!             local logging [fawzi]
!!
!!   SOURCE
!***************************************************************************
  type cp_logger_type
     integer :: id_nr, ref_count
     integer :: print_level
     integer :: default_local_unit_nr
     integer :: default_global_unit_nr
     logical :: close_local_unit_on_dealloc, close_global_unit_on_dealloc
     character(len=10) :: suffix
     character(len=50) :: local_filename, global_filename
     type(print_key_type), pointer :: print_keys
     type(cp_para_env_type), pointer :: para_env
  end type cp_logger_type
!!***
!****************************************************************************

! the default logger
  type(cp_logger_type), pointer, private ::cp_default_logger ! => null()

contains

!!****f* cp_log_handling/cp_init_default_logger [1.0] *
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
  subroutine cp_init_default_logger(para_env,print_level,&
       default_global_unit_nr,default_local_unit_nr, global_filename,&
       local_filename, close_global_unit_on_dealloc,&
       close_local_unit_on_dealloc, print_keys, suffix,&
       template_logger)
    type(cp_para_env_type), pointer, optional :: para_env
    integer, intent(in), optional     :: print_level, default_global_unit_nr,&
         default_local_unit_nr
    character(len=*), intent(in), optional :: global_filename,local_filename,&
         suffix
    logical, intent(in), optional :: close_global_unit_on_dealloc,&
         close_local_unit_on_dealloc
    type(print_key_type), target, optional  :: print_keys
    type(cp_logger_type), target, optional :: template_logger

    character(len=*), parameter :: routineN='cp_init_default_logger',&
         routineP=moduleN//':'//routineN

    if (default_logger_initialized) then
       print *,routineP," WARNING called more than once"
       call cp_logger_release(cp_default_logger)
    else
       nullify(cp_default_logger)
    end if
    ! cases because of bug with optional pointers arg in pgf90 3.3
    if (present(para_env)) then
       call cp_logger_create(cp_default_logger,&
            para_env=para_env,print_level=print_level,&
            default_global_unit_nr=default_global_unit_nr,&
            default_local_unit_nr=default_local_unit_nr,&
            global_filename=global_filename,&
            local_filename=local_filename,&
            close_global_unit_on_dealloc=close_global_unit_on_dealloc,&
            close_local_unit_on_dealloc=close_local_unit_on_dealloc,&
            suffix=suffix, print_keys=print_keys,&
            template_logger=template_logger)
    else
       call cp_logger_create(cp_default_logger,&
            print_level=print_level,&
            default_global_unit_nr=default_global_unit_nr,&
            default_local_unit_nr=default_local_unit_nr,&
            global_filename=global_filename,&
            local_filename=local_filename,&
            close_global_unit_on_dealloc=close_global_unit_on_dealloc,&
            close_local_unit_on_dealloc=close_local_unit_on_dealloc,&
            suffix=suffix, print_keys=print_keys,&
            template_logger=template_logger)
    end if
    start_time=walltime()
    default_logger_initialized=.true.
  end subroutine cp_init_default_logger
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

    if (.not. default_logger_initialized) then
       call stop_program("cp_log_handling:cp_get_default_logger",&
            "default logger not yet initialized (call cp_init_default_logger")
    end if
    if (.not.associated(cp_default_logger)) then
       call stop_program("cp_log_handling:cp_get_default_logger",&
            "default logger was null (released one time too much??)")
    end if
    res => cp_default_logger
  end function cp_get_default_logger
!***************************************************************************

! ================== log ==================

!!****f* cp_log_handling/cp_logger_create [1.0] *
!!
!!   NAME
!!     cp_logger_create
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_create(logger, para_env, print_level,&
!!         default_global_unit_nr, default_local_unit_nr, global_filename,&
!!         local_filename, close_global_unit_on_dealloc,&
!!         close_local_unit_on_dealloc, print_keys, suffix,&
!!         template_logger)
!!       Type(cp_logger_type), Pointer:: logger
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Integer, Intent (IN), Optional:: print_level,&
!!         default_global_unit_nr, default_local_unit_nr
!!       Character(Len=*), Intent (IN), Optional:: global_filename,&
!!         local_filename
!!       Logical, Intent (IN), Optional:: close_global_unit_on_dealloc,&
!!         close_local_unit_on_dealloc
!!       Type(print_key_type), Target, Optional:: print_keys
!!       Character(Len=*), Intent (IN), Optional:: suffix
!!       Type(cp_logger_type), Target, Optional:: template_logger
!!     End Subroutine cp_logger_create
!!
!!   FUNCTION
!!     initializes a logger
!!
!!   NOTES
!!     the handling of *_filename, default_*_unit_nr, close_*_unit_on_dealloc
!!     tries to take the right decision with different inputs, and thus is a
!!     little complex.
!!
!!   INPUTS
!!     - logger: the logger to initialize
!!     - para_env: the parallel environment (this is most likely the global
!!       parallel environment
!!     - print_level: the level starting with which something is written
!!       (defaults to cp_note_level)
!!     - default_global_unit_nr: the default unit_nr for output
!!       (if not given, and no file is given defaults to the standard output)
!!     - default_local_unit_nr: the default unit number for local (i.e. task)
!!       output. If not given defaults to a out.taskid file created upon 
!!     - global_filename: a new file to open (can be given instread of the 
!!       global_unit_nr)
!!     - local_filename: a new file to open (with suffix and para_env%mepos 
!!       appended). Can be given instread of the default_local_unit_nr).
!!       the file is created only upon the first local logging request
!!     - close_global_unit_on_dealloc: if the unit should be closed when the
!!       logger is deallocated (defaults to true if a local_filename is given,
!!       to false otherwise)
!!     - close_local_unit_on_dealloc: if the unit should be closed when the
!!       logger is deallocated (defaults to true)
!!     - suffix: the suffix that should be added to all the generated filenames
!!     - template_logger: a logger from where to take the unspecified things
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_create(logger, para_env,print_level,&
       default_global_unit_nr,default_local_unit_nr, global_filename,&
       local_filename, close_global_unit_on_dealloc,&
       close_local_unit_on_dealloc, print_keys, suffix, template_logger)
    type(cp_logger_type), pointer           :: logger
    type(cp_para_env_type),pointer,optional :: para_env
    integer, intent(in), optional           :: print_level,&
         default_global_unit_nr, default_local_unit_nr
    character(len=*), intent(in), optional  :: global_filename,local_filename
    logical, intent(in), optional           :: close_global_unit_on_dealloc,&
         close_local_unit_on_dealloc
    type(print_key_type), target, optional  :: print_keys
    character(len=*), intent(in),optional   :: suffix
    type(cp_logger_type),target,optional    :: template_logger

    character(len=*), parameter :: routineN='cp_logger_create',&
         routineP=moduleN//':'//routineN
    
    integer :: stat
    logical :: exists

    allocate(logger, stat=stat)
    if (stat/=0) call stop_program(routineP,"could not allocate a logger")

    nullify(logger%para_env,logger%print_keys)
    logger%ref_count=1
    last_logger_id_nr=last_logger_id_nr+1
    logger%id_nr=last_logger_id_nr

    if (present(template_logger)) then
      if (template_logger%ref_count<1) &
            call stop_program(routineP,"template_logger%ref_count<1")
      logger%print_level=template_logger%print_level
      logger%default_global_unit_nr=template_logger%default_global_unit_nr
      logger%close_local_unit_on_dealloc=template_logger%close_local_unit_on_dealloc
      if (logger%close_local_unit_on_dealloc) then
         logger%default_local_unit_nr=-1
      else
         logger%default_local_unit_nr=template_logger%default_local_unit_nr
      end if
      logger%close_global_unit_on_dealloc=template_logger%close_global_unit_on_dealloc
      if (logger%close_global_unit_on_dealloc) then
         logger%default_global_unit_nr=-1
      else
         logger%default_global_unit_nr=template_logger%default_global_unit_nr
      end if
      logger%local_filename=template_logger%local_filename
      logger%global_filename=template_logger%global_filename
      logger%print_keys => template_logger%print_keys
      logger%para_env => template_logger%para_env
      logger%suffix=template_logger%suffix
   else
    logger%default_global_unit_nr=cp_std_out_unit_nr
    logger%local_filename="localLog"
    logger%global_filename="mainLog"
    logger%print_level=cp_note_level
    logger%close_global_unit_on_dealloc=.false.
#ifdef __parallel
    logger%default_local_unit_nr=-1
    logger%close_local_unit_on_dealloc=.true.
#else
    logger%default_local_unit_nr=cp_std_out_unit_nr
    logger%close_local_unit_on_dealloc=.false.
#endif
    logger%suffix=""
   end if
    if (present(para_env)) logger%para_env => para_env
    if (.not.associated(logger%para_env)) call stop_program(routineP,&
         "para env not associated")
    if (logger%para_env%ref_count<1) call stop_program(routineP,&
         "para_env%ref_count<1")
    logger%para_env%ref_count = logger%para_env%ref_count+1

    if (present(print_level)) logger%print_level=print_level

    if (present(default_global_unit_nr)) &
         logger%default_global_unit_nr=default_global_unit_nr
    if (present(global_filename)) then
       logger%global_filename=global_filename
       logger%close_global_unit_on_dealloc=.true.
       logger%default_global_unit_nr=-1
    end if
    if (present(close_global_unit_on_dealloc)) then
       logger%close_global_unit_on_dealloc=close_global_unit_on_dealloc
       if (present(default_global_unit_nr).and.present(global_filename).and.&
            .not.close_global_unit_on_dealloc) then
         logger%default_global_unit_nr=default_global_unit_nr
       end if
    end if
    
    if (present(default_local_unit_nr)) &
         logger%default_local_unit_nr=default_local_unit_nr
    if (present(local_filename)) then
       logger%local_filename=local_filename
       logger%close_local_unit_on_dealloc=.true.
       logger%default_local_unit_nr=-1
    end if
    if (present(close_local_unit_on_dealloc)) then
       logger%close_local_unit_on_dealloc=close_local_unit_on_dealloc
       if (present(default_local_unit_nr).and.present(local_filename).and.&
            .not.close_local_unit_on_dealloc) then
         logger%default_local_unit_nr=default_local_unit_nr
       end if
    end if

    if (present(print_keys)) logger%print_keys => print_keys
    if (associated(logger%print_keys)) then
       if (logger%print_keys%ref_count<1) then
          call stop_program(routineP,"print_keys%ref_count<1")
       end if
       logger%print_keys%ref_count=print_keys%ref_count+1
    else
       allocate(logger%print_keys,stat=stat)
       if (stat/=0) then
          call stop_program(routineP,"print_keys allocation stat="//&
               cp_to_string(stat))
          call init_print_keys(logger%print_keys,level=3-logger%print_level)
       end if
    end if
  end subroutine cp_logger_create
!***************************************************************************

!!****f* cp_log_handling/cp_logger_retain [1.0] *
!!
!!   NAME
!!     cp_logger_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_retain(logger)
!!       Type(cp_logger_type), Pointer:: logger
!!     End Subroutine cp_logger_retain
!!
!!   FUNCTION
!!     retains the given logger (to be called to keep a shared copy of 
!!     the logger)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_retain(logger)
    type(cp_logger_type), pointer :: logger

    logical :: failure
    character(len=*), parameter :: routineN='cp_logger_retain',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (.not.associated(logger)) call stop_program(routineP,&
         "logger not associated")
    if (logger%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")
    logger%ref_count=logger%ref_count+1
  end subroutine cp_logger_retain
!***************************************************************************

!!****f* cp_logger_release *
!!
!!   NAME
!!     cp_logger_release
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_release(logger)
!!       Type(cp_logger_type), Pointer:: logger
!!     End Subroutine cp_logger_release
!!
!!   FUNCTION
!!     releases this logger
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger to release
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_release(logger)
    type(cp_logger_type), pointer :: logger

    character(len=*), parameter :: routineN='cp_logger_release',&
         routineP=moduleN//':'//routineN
    integer :: stat

    if (associated(logger)) then
       if (logger%ref_count<1) call stop_program(routineP,"logger%ref_count<1")
       logger%ref_count=logger%ref_count-1
       if (logger%ref_count<1) then
          if (logger%close_global_unit_on_dealloc .and. &
               logger%default_global_unit_nr >= 0) then
             call close_file(logger%default_global_unit_nr)
             logger%close_global_unit_on_dealloc=.false.
             logger%default_global_unit_nr=-1
          end if
          if (logger%close_local_unit_on_dealloc .and. &
               logger%default_local_unit_nr >= 0) then
             call close_file(logger%default_local_unit_nr)
             logger%close_local_unit_on_dealloc=.false.
             logger%default_local_unit_nr=-1
          end if
          if (associated(logger%print_keys)) then
             if (logger%print_keys%ref_count<1) then
                call stop_program(routineP," logger%print_keys%ref_count<1")
             end if
             logger%print_keys%ref_count=logger%print_keys%ref_count-1
             if (logger%print_keys%ref_count<1) then
                deallocate(logger%print_keys,stat=stat)
                if (stat/=0) print *, routineP,&
                     " WARNING print_keys dealloc stat=",stat
             end if
          end if
          call my_cp_para_env_release(logger%para_env)
          deallocate(logger,stat=stat)
          if (stat/=0) print *, routineP," WARNING logger dealloc stat=",stat
       end if
    end if
    nullify(logger)
  end subroutine cp_logger_release
!***************************************************************************

!!****f* cp_log_handling/cp_logger_would_log [1.0] *
!!
!!   NAME
!!     cp_logger_would_log
!!
!!   SYNOPSIS
!!     Function cp_logger_would_log(logger, level, fromwhere, local)&
!!         Result(res)
!!       Logical:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Logical, Intent (IN), Optional:: local
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
!!  - local: if the output is local (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure
  function cp_logger_would_log(logger, level, fromWhere, local) result(res)
    logical :: res
    type(cp_logger_type), pointer ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere
    logical, intent(in), optional :: local

    character(len=*), parameter :: routineN='cp_logger_would_log',&
         routineP=moduleN//':'//routineN
    type(cp_logger_type), pointer ::lggr

    lggr => logger
    if (.not.associated(lggr)) lggr => cp_get_default_logger()
    if (lggr%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")

    res = level>=lggr%print_level
  end function cp_logger_would_log

!!****f* cp_log_handling/cp_logger_get_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_unit_nr(logger, level, fromwhere, local)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Logical, Intent (IN), Optional:: local
!!     End Function cp_logger_get_unit_nr
!!
!!   FUNCTION
!!     returns the unit nr for the requested kind of log.
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
!!  - local: if true returns a local logger (one per task), otherwise
!!    returns a global logger (only the process with para_env%mepos==
!!    para_env%source should write to the global logger). Defaults to
!!    false
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  function cp_logger_get_unit_nr(logger, level, fromWhere, local) result(res)
    integer :: res
    type(cp_logger_type), pointer :: logger
    integer, intent(in) :: level
    character(len=*), intent(in):: fromWhere
    logical, intent(in), optional :: local

    logical :: loc=.false.

    res=cp_logger_get_default_unit_nr(logger,local=local)
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
!!     At the moment it just goes to the default unit of the logger
!!
!!   INPUTS
!!  - logger is a type that describes where the logging must go.
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - message is the message to print
!!  - local: if true this goes to the local logger, otherwise (the default)
!!    it goes to the global logger (and only the master process writes it)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  subroutine cp_logger_log(logger, level, fromWhere , message, local)
    type(cp_logger_type), pointer ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere, message
    logical, intent(in), optional :: local

    character(len=*), parameter :: routineN='cp_logger_log',&
         routineP=moduleN//':'//routineN
    character(7):: label
    integer ::oErr, unitNr
    logical :: loc

    type(cp_logger_type), pointer ::lggr
    lggr => logger
    if (.not.associated(lggr)) lggr => cp_get_default_logger()
    if (lggr%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")

    loc=.false.
    if (present(local)) loc=local
    if (cp_logger_would_log(lggr,level,fromWhere,local=loc)) then
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
       unitNr=cp_logger_get_unit_nr(lggr,level,fromWhere,local=loc)
       if (loc .or. lggr%para_env%mepos==lggr%para_env%source) then
          call print_message(label//' in '//fromWhere//',t'//&
               trim(adjustl(cp_to_string(walltime()-start_time)))//&
               ' '//message,unitNr,&
               level, level,level/2)
! old (more failsafe, no trim, compress,...) code:
!write (unit=unitNr,fmt='(a)', iostat=oErr) ' '
!write (unit=unitNr,fmt='(5a)', iostat=oErr) label,' in ', fromWhere,',', message
!if (oErr /= 0) then
!   print *, 'ERROR  ',oErr,' while writing log, ', fromWhere, message
!end if
          call m_flush(unitNr)
       end if
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
!!       Type(cp_logger_type), Pointer:: logger
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
    type(cp_logger_type), pointer :: logger
    integer, intent(in) :: level

    character(len=*), parameter :: routineN='cp_logger_set_log_level',&
         routineP=moduleN//':'//routineN

    if (.not.associated(logger)) call stop_program(routineP,&
         "logger not associated")
    if (logger%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")
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
  function cp_logger_get_log_level(logger) result(res)
    integer :: res
    type(cp_logger_type), pointer :: logger

    character(len=*), parameter :: routineN='cp_logger_get_log_level',&
         routineP=moduleN//':'//routineN
    type(cp_logger_type), pointer ::lggr

    lggr => logger
    if (.not.associated(lggr)) lggr => cp_get_default_logger()
    if (lggr%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")

    res=lggr%print_level
  end function cp_logger_get_log_level

!!****f* cp_log_handling/cp_logger_get_default_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_default_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_default_unit_nr(logger, local) Result(res)
!!       Integer:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Logical, Intent (IN), Optional:: local
!!     End Function cp_logger_get_default_unit_nr
!!
!!   FUNCTION
!!     asks the default unit number of the given logger.
!!     try to use cp_logger_get_unit_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger you want info from
!!     - local: if you want the local unit nr (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
!pure 
  function cp_logger_get_default_unit_nr(logger, local) result(res)
    integer :: res
    type(cp_logger_type), pointer :: logger
    logical, intent(in), optional :: local

    logical :: loc,exists
    character(len=*), parameter :: routineN='cp_logger_get_default_unit_nr',&
         routineP=moduleN//':'//routineN
    character(len=80) :: filename
    type(cp_logger_type), pointer ::lggr
    loc=.false.

    lggr => logger
    if (.not.associated(lggr)) lggr => cp_get_default_logger()
    if (lggr%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")

    if (present(local)) loc=local
    if (.not.local) then
       if (lggr%default_global_unit_nr<=0) then
          if (lggr%para_env%mepos==lggr%para_env%source) then
             call cp_logger_generate_filename(lggr,filename,lggr%global_filename,&
                  ".out",local=.false.)
             INQUIRE (FILE=trim(filename),EXIST=exists)
             if (exists) then
                print *," WARNING ", routineP, " overwriting file ",filename
             end if
             call open_file(trim(filename),file_status="unknown",&
                  file_action="write",&
                  unit_number=lggr%default_global_unit_nr)
          end if
          ! bcast? dangerous if called only in source processor
          !call mp_bcast(lggr%default_global_unit_nr,lggr%para_env%source,&
          !     lggr%para_env%group)
       end if
       res=lggr%default_global_unit_nr
    else
       if (lggr%default_local_unit_nr<0) then
          call cp_logger_generate_filename(lggr,filename,lggr%local_filename,&
               ".out",local=.true.)
          INQUIRE (FILE=trim(filename),EXIST=exists)
          if (exists) then
             print *," WARNING ", routineP, " overwriting file ",filename
          end if
          call open_file(trim(filename),file_status="unknown",&
               file_action="write",&
               unit_number=lggr%default_local_unit_nr)
       end if
       res=lggr%default_local_unit_nr
    end if
  end function cp_logger_get_default_unit_nr

!!****f* cp_log_handling/my_cp_para_env_release [1.0] *
!!
!!   NAME
!!     my_cp_para_env_release
!!
!!   FUNCTION
!!     releases the para object (to be called when you don't want anymore
!!     the shared copy of this object)
!!
!!   NOTES
!!     local private version to avoid circular dependencies
!!     (see cp_para_env:cp_para_env_release)!
!!
!!   INPUTS
!!     - para_env: the new group
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine my_cp_para_env_release(para_env)
    type(cp_para_env_type), pointer :: para_env

    character(len=*), parameter :: routineN='my_cp_para_env_release',&
         routineP=moduleN//':'//routineN
    integer :: stat

    if (associated(para_env)) then
       if (para_env%ref_count<1) then
          call stop_program(routineP,"para_env%ref_count<1")
       end if
       para_env%ref_count=para_env%ref_count-1
       if (para_env%ref_count<1) then
          if (para_env%owns_group) then
             call mp_comm_free(para_env%group)
          end if
          deallocate(para_env,stat=stat)
          if(stat/=0)print *,routineP," WARNING para_env deallocate stat=",stat
       end if
    end if
    nullify(para_env)
  end subroutine my_cp_para_env_release
!***************************************************************************

!!****f* cp_log_handling/cp_logger_generate_filename [1.0] *
!!
!!   NAME
!!     cp_logger_generate_filename
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_generate_filename(logger, res, root, postfix,&
!!         local)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (INOUT):: res
!!       Character(Len=*), Intent (IN):: root, postfix
!!       Logical, Intent (IN), Optional:: local
!!     End Subroutine cp_logger_generate_filename
!!
!!   FUNCTION
!!     generates a unique filename (ie adding eventual suffixes and 
!!     process ids)
!!
!!   NOTES
!!     this should be a function returning a variable length string.
!!     All spaces are moved to the end of the string.
!!     Not fully optimized: result must be a littel longer than the 
!!     resulting compressed filename
!!
!!   INPUTS
!!     - result: the resulting string
!!     - root: the start of filename
!!     - postfix: the end of the name
!!     - local: if the name should be local to this task (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_logger_generate_filename(logger,res, root, postfix,&
       local)
    type(cp_logger_type), pointer   :: logger
    character(len=*), intent(inout) :: res
    character(len=*), intent(in)    :: root, postfix
    logical, intent(in), optional   :: local

    logical :: failure, loc
    character(len=*), parameter :: routineN='cp_logger_generate_filename',&
         routineP=moduleN//':'//routineN
    type(cp_logger_type), pointer :: lggr
    failure=.false.
    loc=.false.

    lggr => logger
    if (.not.associated(lggr)) lggr => cp_get_default_logger()
    if (lggr%ref_count<1) call stop_program(routineP,&
         "logger%ref_count<1")
    if (present(local)) loc=local
    if (loc) then
       res=trim(root)//trim(lggr%suffix)//'_p'//&
            cp_to_string(lggr%para_env%mepos)//postfix
    else
       res=trim(root)//trim(lggr%suffix)//postfix
    end if
    call compress_all(res)
  end subroutine cp_logger_generate_filename
!***************************************************************************

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
!!     desactivated because unused, and on cray single==double
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
!!FM  function cp_real_sp_to_string(val) result(res)
!!FM    character(len=10) :: res
!!FM    real(kind=kind(1.0)), intent(in):: val
!!FM    write (res,'(es10.4)') val
!!FM  end function cp_real_sp_to_string

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
    character(len=11) :: res
    real(kind=kind(1.0D0)), intent(in):: val
    write (res,'(es11.4)') val
  end function cp_real_dp_to_string

end module cp_log_handling

!******************************************************************************

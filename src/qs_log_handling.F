!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_log_handling [1.0] *
!!
!!   NAME
!!     qs_log_handling
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
!!     At the moment logging and output refer to the same object (qs_logger)
!!     as these are actually different it might be better to separate them
!!     (they have already separate routines in a separate module 
!!     @see qs_output_handling).
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!     @see qs_error_handling
!!
!!   SOURCE
!******************************************************************************

module qs_log_handling
  use output_utilities, only: print_message
  implicit none
  public

  logical, private, parameter :: debug_this_module=.true.
  integer, parameter :: qs_std_in_unit_nr=5
  integer, parameter :: qs_std_out_unit_nr=6

!! type of a logger, at the moment it contains just a print level
!! starting at which level it should be logged 
!! (0 note, 1 warning, 2 failure, 3 fatal)
!! it could be expanded with the ability to focus on one or more
!! module/object/thread/processor
!!
!! This should be private, but as the output functions have been moved to 
!! another module and there is no "friend" keyword, it is not.
!! DO NOT USE THE INTERNAL COMPONENTS DIRECTLY!!!
  type qs_logger
     integer ::printLevel
     integer ::unitNr
  end type qs_logger

!! the default local logger
  type(qs_logger), target ::qs_local_logger=qs_logger(0,-1)
!! the default global logger
  type(qs_logger), target ::qs_global_logger=qs_logger(0,-2)
!! the default logger (normally points to the local logger)
  type(qs_logger), pointer ::qs_default_logger ! => null()

!! level of an error
  integer, parameter  :: qs_fatal_level=3
!! level of a failure
  integer, parameter  :: qs_failure_level=2
!! level of a warning
  integer, parameter  :: qs_warning_level=1
!! level of a note
  integer, parameter  :: qs_note_level=0

!! a generic function to trasform different types to strings
  interface qs_to_string
     module procedure qs_int_to_string, qs_real_to_string_sp, qs_real_to_string_dp
  end interface

contains
! I should convert my f90 doc comments to robodoc comments
! (I decided to use robodoc after having written most of this)


!! performs the initializations of the logger.
!! MUST be called before logging occours
  subroutine qs_init_default_loggers ! add optional initialization param?
    logical,save ::firstTime=.true.
    if (firstTime) then
       firstTime=.false.
       qs_local_logger%unitNr=qs_std_out_unit_nr
       qs_global_logger%unitNr=qs_std_out_unit_nr
       qs_default_logger => qs_local_logger
    end if
  end subroutine qs_init_default_loggers

  ! **************************** to string *******************************

!! converts an int to a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
 ! pure
  function qs_int_to_string(i) result(res)
    character(len=6) :: res
    integer, intent(in):: i
    write (res,'(i6)') i
    if (i>999999 .or. i<-99999) then
       write (res, '(es6.1)') real(i)
    end if
  end function qs_int_to_string

!! convert a single precision real in a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
 ! pure 
  function qs_real_to_string_sp(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0)), intent(in):: val
    write (res,'(es10.4)') val
  end function qs_real_to_string_sp

!! convert a double precision real in a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
  ! pure 
  function qs_real_to_string_dp(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0D0)), intent(in):: val
    write (res,'(es10.4)') val
  end function qs_real_to_string_dp

  ! ******************************** log *********************************

!! this function can be called to check if the logger would log
!! a message with the given level from the given source
!! you should use this function if you do direct logging
!! (without using qs_log), or if you want to know if the generation 
!! of some costly log info is necessary
  ! pure
  function qs_would_log(logger, level, fromWhere)
    logical :: qs_would_log
    type(qs_logger), intent(in) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere

    qs_would_log=level>=logger%printLevel
  end function qs_would_log

!! returns the unit nr for the requested kind of log.
!!
!! Try not to use this function, in future in some cases the log 
!! could not use units (global log might be sent to the master process)
  ! pure 
  function qs_unitnr_for_log(logger, level, fromWhere)
    integer :: qs_unitnr_for_log
    type(qs_logger), intent(in) :: logger
    integer, intent(in) :: level
    character(len=*), intent(in):: fromWhere

    if (logger%unitNr<=0) then
       qs_unitnr_for_log=qs_std_out_unit_nr
    else
       qs_unitnr_for_log=logger%unitNr
    end if
  end function qs_unitnr_for_log

!! Logs in a standard format
!!  - logger is a type that describes where the logging must go.
!!    There are three standard loggers: qs_local_logger, qs_global_logger,
!!    and qs_default_logger, that normally points to one of the previous two
!!  - level describes the of the message: qs_fatal_level, 
!!    qs_failure_level(2), qs_warning_level(1), qs_note_level.
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - message is the message to print
!! In multithreaded situations it should be modified to print also 
!! the thread (and maybe the processor), and proper locking should be 
!! performed.
!! At the moment it just goes to the standard output
  subroutine qs_log(logger, level, fromWhere , message)
    type(qs_logger) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere, message

    character(7):: label
    integer ::oErr, unitNr
    if (.not. associated(qs_default_logger)) then
       print *, 'ERROR: qs_default_logger is not associated'
       print *, 'ERROR: did you call qs_init_default_loggers at startup??'
    end if
    if (qs_would_log(logger,level,fromWhere)) then
       select case (level)
       case (qs_fatal_level)
          label='ERRORL3'
       case (qs_failure_level)
          label='ERRORL2'
       case (qs_warning_level)
          label='WARNING'
       case (qs_note_level)
          label='NOTE'
       case default
          label='L'//adjustl(qs_int_to_string(level))
       end select
       unitNr=qs_unitnr_for_log(logger,level,fromWhere)
       
       call print_message(label//' in '//fromWhere//', '//message,unitNr,&
            level, level,level/2)
       ! old (more failsafe, no trim, compress,...) code:
       !write (unit=unitNr,fmt='(a)', iostat=oErr) ' '
       !write (unit=unitNr,fmt='(5a)', iostat=oErr) label,' in ', fromWhere,',', message
       !if (oErr /= 0) then
       !   print *, 'ERROR  ',oErr,' while writing log, ', fromWhere, message
       !end if
    end if
  end subroutine qs_log

  ! *************************** qs_logger settings ***************************

!! changes the logging level. Log messages with a level less than the one
!! given wo not be printed.
!!
!! This routine might change/disappear!!
  subroutine qs_set_log_level(logger, level)
    type(qs_logger), intent(out) :: logger
    integer, intent(in) :: level
    logger%printLevel=level
  end subroutine qs_set_log_level

!! asks the logging level of the given logger.
!!
!! This routine might change/disappear!!
  pure function qs_get_log_level(logger)
    integer :: qs_get_log_level
    type(qs_logger), intent(in) :: logger
    qs_get_log_level=logger%printLevel
  end function qs_get_log_level

!! changes the output unit of the given logger.
!!
!! This routine might change/disappear!!
  subroutine qs_set_log_unitnr(logger, unitNr)
    type(qs_logger), intent(out) :: logger
    integer, intent(in) :: unitNr
    logger%unitNr=unitNr
  end subroutine qs_set_log_unitnr

!! asks the unit number of the given logger.
!!
!! This routine might change/disappear!!
  pure function qs_get_log_unitnr(logger)
    integer :: qs_get_log_unitnr
    type(qs_logger), intent(in) :: logger
    qs_get_log_unitnr=logger%unitNr
  end function qs_get_log_unitnr

end module qs_log_handling

!!***
!******************************************************************************

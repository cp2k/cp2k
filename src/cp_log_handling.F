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
!!     At the moment logging and output refer to the same object (cp_logger)
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
  implicit none
  public

  logical, private, parameter :: debug_this_module=.true.
  integer, parameter :: cp_std_in_unit_nr=5
  integer, parameter :: cp_std_out_unit_nr=6

!! type of a logger, at the moment it contains just a print level
!! starting at which level it should be logged 
!! (0 note, 1 warning, 2 failure, 3 fatal)
!! it could be expanded with the ability to focus on one or more
!! module/object/thread/processor
!!
!! This should be private, but as the output functions have been moved to 
!! another module and there is no "friend" keyword, it is not.
!! DO NOT USE THE INTERNAL COMPONENTS DIRECTLY!!!
  type cp_logger
     integer ::printLevel
     integer ::unitNr
  end type cp_logger

!! the default local logger
  type(cp_logger), target ::cp_local_logger=cp_logger(0,-1)
!! the default global logger
  type(cp_logger), target ::cp_global_logger=cp_logger(0,-2)
!! the default logger (normally points to the local logger)
  type(cp_logger), pointer ::cp_default_logger ! => null()

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
     module procedure cp_int_to_string, cp_real_to_string_sp, cp_real_to_string_dp
  end interface

contains
! I should convert my f90 doc comments to robodoc comments
! (I decided to use robodoc after having written most of this)


!! performs the initializations of the logger.
!! MUST be called before logging occours
  subroutine cp_init_default_loggers ! add optional initialization param?
    logical,save ::firstTime=.true.
    if (firstTime) then
       firstTime=.false.
       cp_local_logger%unitNr=cp_std_out_unit_nr
       cp_global_logger%unitNr=cp_std_out_unit_nr
       cp_default_logger => cp_local_logger
    end if
  end subroutine cp_init_default_loggers

  ! **************************** to string *******************************

!! converts an int to a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
 ! pure
  function cp_int_to_string(i) result(res)
    character(len=6) :: res
    integer, intent(in):: i
    write (res,'(i6)') i
    if (i>999999 .or. i<-99999) then
       write (res, '(es6.1)') real(i)
    end if
  end function cp_int_to_string

!! convert a single precision real in a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
 ! pure 
  function cp_real_to_string_sp(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0)), intent(in):: val
    write (res,'(es10.4)') val
  end function cp_real_to_string_sp

!! convert a double precision real in a string
!! (should be a variable length string, but does not work with 
!! all the compilers)
  ! pure 
  function cp_real_to_string_dp(val) result(res)
    character(len=10) :: res
    real(kind=kind(1.0D0)), intent(in):: val
    write (res,'(es10.4)') val
  end function cp_real_to_string_dp

  ! ******************************** log *********************************

!! this function can be called to check if the logger would log
!! a message with the given level from the given source
!! you should use this function if you do direct logging
!! (without using cp_log), or if you want to know if the generation 
!! of some costly log info is necessary
  ! pure
  function cp_would_log(logger, level, fromWhere)
    logical :: cp_would_log
    type(cp_logger), intent(in) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere

    cp_would_log=level>=logger%printLevel
  end function cp_would_log

!! returns the unit nr for the requested kind of log.
!!
!! Try not to use this function, in future in some cases the log 
!! could not use units (global log might be sent to the master process)
  ! pure 
  function cp_unitnr_for_log(logger, level, fromWhere)
    integer :: cp_unitnr_for_log
    type(cp_logger), intent(in) :: logger
    integer, intent(in) :: level
    character(len=*), intent(in):: fromWhere

    if (logger%unitNr<=0) then
       cp_unitnr_for_log=cp_std_out_unit_nr
    else
       cp_unitnr_for_log=logger%unitNr
    end if
  end function cp_unitnr_for_log

!! Logs in a standard format
!!  - logger is a type that describes where the logging must go.
!!    There are three standard loggers: cp_local_logger, cp_global_logger,
!!    and cp_default_logger, that normally points to one of the previous two
!!  - level describes the of the message: cp_fatal_level, 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level.
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - message is the message to print
!! In multithreaded situations it should be modified to print also 
!! the thread (and maybe the processor), and proper locking should be 
!! performed.
!! At the moment it just goes to the standard output
  subroutine cp_log(logger, level, fromWhere , message)
    type(cp_logger) ::logger
    integer , intent(in) :: level
    character(len=*), intent(in):: fromWhere, message

    character(7):: label
    integer ::oErr, unitNr
    if (.not. associated(cp_default_logger)) then
       print *, 'ERROR: cp_default_logger is not associated'
       print *, "ERROR: did you call cp_init_default_loggers at startup?"
    end if
    if (cp_would_log(logger,level,fromWhere)) then
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
       unitNr=cp_unitnr_for_log(logger,level,fromWhere)
       
       call print_message(label//' in '//fromWhere//', '//message,unitNr,&
            level, level,level/2)
       ! old (more failsafe, no trim, compress,...) code:
       !write (unit=unitNr,fmt='(a)', iostat=oErr) ' '
       !write (unit=unitNr,fmt='(5a)', iostat=oErr) label,' in ', fromWhere,',', message
       !if (oErr /= 0) then
       !   print *, 'ERROR  ',oErr,' while writing log, ', fromWhere, message
       !end if
    end if
  end subroutine cp_log

  ! *************************** cp_logger settings ***************************

!! changes the logging level. Log messages with a level less than the one
!! given wo not be printed.
!!
!! This routine might change/disappear!!
  subroutine cp_set_log_level(logger, level)
    type(cp_logger), intent(out) :: logger
    integer, intent(in) :: level
    logger%printLevel=level
  end subroutine cp_set_log_level

!! asks the logging level of the given logger.
!!
!! This routine might change/disappear!!
  !pure
  function cp_get_log_level(logger)
    integer :: cp_get_log_level
    type(cp_logger), intent(in) :: logger
    cp_get_log_level=logger%printLevel
  end function cp_get_log_level

!! changes the output unit of the given logger.
!!
!! This routine might change/disappear!!
  subroutine cp_set_log_unitnr(logger, unitNr)
    type(cp_logger), intent(out) :: logger
    integer, intent(in) :: unitNr
    logger%unitNr=unitNr
  end subroutine cp_set_log_unitnr

!! asks the unit number of the given logger.
!!
!! This routine might change/disappear!!
  !pure 
  function cp_get_log_unitnr(logger)
    integer :: cp_get_log_unitnr
    type(cp_logger), intent(in) :: logger
    cp_get_log_unitnr=logger%unitNr
  end function cp_get_log_unitnr

end module cp_log_handling

!!***
!******************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_output_handling [1.0] *
!!
!!   NAME
!!     cp_output_handling
!!
!!   FUNCTION
!!     routines to handle the output, The idea is to remove the 
!!     decision of wheter to output and what to output from the code
!!     that does the output, and centralize it here.
!!
!!   NOTES
!!     These were originally together with the log handling routines,
!!     but have been spawned off. Some dependencies are still there,
!!     and some of the comments about log handling also applies to output
!!     handling: @see cp_log_handling
!!
!!     A section with output name, output every (-1: last, 0: never, 1: always)
!!     separate File, discriminate fromWhere,... should be implemented.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_output_handling
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  implicit none

  logical, private, parameter :: debug_this_module=.true.
  
  public :: cp_would_output, cp_unitnr_for_output,&
       cp_write_output

  !****************************************************************************
contains

!!****f* cp_output_handling/cp_would_output [1.0] *
!!
!!   NAME
!!     cp_would_output
!!
!!   SYNOPSIS
!!     Function cp_would_output(logger, outputname, fromwhere, iter)
!!       Logical:: cp_would_output
!!       Type(cp_logger_type), Intent (IN):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!     End Function cp_would_output
!!
!!   FUNCTION
!!    this function can be called to check if the logger would
!!    write the output. This can be used to know if the generation 
!!    of some costly output is necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  ! pure 
  function cp_would_output(logger, outputName , fromWhere, iter)
    logical :: cp_would_output
    type(cp_logger_type), intent(in) ::logger
!! name of the output, for example 'hamiltonian'
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter

    cp_would_output=.true.
  end function cp_would_output
!!***
  !****************************************************************************

!!****f* cp_output_handling/cp_unitnr_for_output [1.0] *
!!
!!   NAME
!!     cp_unitnr_for_output
!!
!!   SYNOPSIS
!!     Function cp_unitnr_for_output(logger, outputname, fromwhere, iter)
!!       Integer:: cp_unitnr_for_output
!!       Type(cp_logger_type), Intent (IN):: logger
!!       Character(Len=*), Intent (IN):: outputname, fromwhere
!!       Integer, Intent (IN):: iter
!!     End Function cp_unitnr_for_output
!!
!!   FUNCTION
!!    returns the unit nr for the requested kind of output.
!!
!!    Try not to use this function, in future in some cases the log 
!!    could not use units (global log might be sent to the master process)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  !pure 
  function cp_unitnr_for_output(logger, outputName, fromWhere, iter)
    integer :: cp_unitnr_for_output
    type(cp_logger_type), intent(in) :: logger
    character(len=*), intent(in):: outputName, fromWhere
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter

    if (logger%default_unit_nr<=0) then
       cp_unitnr_for_output=cp_std_out_unit_nr
    else
       cp_unitnr_for_output=logger%default_unit_nr
    end if
  end function cp_unitnr_for_output
!!***
  !****************************************************************************

!!****f* cp_output_handling/cp_write_output [1.0] *
!!
!!   NAME
!!     cp_write_output
!!
!!   SYNOPSIS
!!     Subroutine cp_write_output(logger, outputname, fromwhere, iter,&
!!         message, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!       Integer, Intent (IN):: iter
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_write_output
!!
!!   FUNCTION
!!     Used to write the result (output) of something
!!     no frills, no label, just the message is written.
!!     In the future some discrimination based on outputName could be 
!!     introduced.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  subroutine cp_write_output(logger, outputName , fromWhere , iter, message,&
       error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere, message
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineP='cp_output_handling:cp_write_output'
    logical :: failure
    integer ::unitNr,oErr
    failure=.false.
    if (cp_would_output(logger,outputName,fromWhere,iter)) then
       unitNr=cp_unitnr_for_output(logger,outputName,fromWhere,iter) 
       write (unit=unitNr,fmt='(a)', iostat=oErr) message
       CPPostcondition(oErr == 0,cp_failure_level,routineP,error,failure)
    end if
  end subroutine cp_write_output
!!***
  !****************************************************************************

end module cp_output_handling


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
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_output_handling'

  public :: cp_would_output, cp_unitnr_for_output,&
       cp_write_output, cp_finished_output

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
  !!     - logger: the logger that decides where the logging should go
  !!     - outputName: the name of the output, for example 'scf_hamiltonian'
  !!       this should be unique (i.e. no other routine should have the same
  !!       name, because it is used to build an unique identifier)
  !!     - fromWhere: string of the form module:function or file:lineNr
  !!       that gives the origin of the output
  !!     - iter: the iteration number, or some other number that permit to 
  !!       differentiate between different instances of the same object
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 12.2001
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 updated to new logger [fawzi]
  !!
  !!*** *********************************************************************
  ! pure 
  function cp_would_output(logger, outputName , fromWhere, iter, local, error)
    logical :: cp_would_output
    type(cp_logger_type), pointer ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    logical, intent(in), optional :: local
    type(cp_error_type), intent(inout), optional :: error

    cp_would_output=.true.
  end function cp_would_output
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
  !!    you should call cp_finished_output when an iteration output is
  !!    finished (to immediately close the file that might have been opened)
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
  !!     08.2002 updated to new logger [fawzi]
  !!
  !!*** *********************************************************************
  !pure 
  function cp_unitnr_for_output(logger, outputName, fromWhere, iter, local,&
       error) result(res)
    integer :: res
    type(cp_logger_type), pointer :: logger
    character(len=*), intent(in):: outputName, fromWhere
    integer, intent(in) :: iter
    logical, intent(in), optional :: local
    type(cp_error_type), intent(inout), optional :: error

    res=cp_logger_get_default_unit_nr(logger,local=local)
  end function cp_unitnr_for_output
  !****************************************************************************

  !!****f* cp_output_handling/cp_finished_output [1.0] *
  !!
  !!   NAME
  !!     cp_finished_output
  !!
  !!   FUNCTION
  !!     should be called after you finish working with a unit obtained with
  !!     cp_unitnr_for_output, so that the file that might have been opened
  !!     can be closed.
  !!
  !!   NOTES
  !!     at the moment does nothing, but this should change in the future
  !!
  !!   INPUTS
  !!     the inputs should be exactly the same of the corresponding
  !!     cp_unitnr_for_output
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 created [fawzi]
  !!
  !!*** **********************************************************************
  subroutine cp_finished_output(logger, outputName, fromWhere, iter, &
       local,error)
    integer :: cp_unitnr_for_output
    type(cp_logger_type), pointer :: logger
    character(len=*), intent(in):: outputName, fromWhere
    integer, intent(in) :: iter
    logical, intent(in), optional :: local
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_finished_output',&
         routineP=moduleN//':'//routineN
    failure=.false.

  end subroutine cp_finished_output
  !***************************************************************************

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
  !!     - logger: the logger that decides where the logging should go
  !!     - outputName: the name of the output, for example 'scf_hamiltonian'.
  !!       this should be unique (i.e. no other routine should have the same
  !!       name, because it is used to build an unique identifier)
  !!     - fromWhere: string of the form module:function or file:lineNr
  !!       that gives the origin of the output
  !!     - iter: the iteration number, or some other number that permit to 
  !!       differentiate between different instances of the same object
  !!     - local: it the output is local to the task or not. defaults to false
  !!     
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 12.2001
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 updated to new logger [fawzi]
  !!
  !!*** *********************************************************************
  subroutine cp_write_output(logger, outputName , fromWhere , iter, message,&
       local,error)
    type(cp_logger_type), pointer :: logger
    character(len=*), intent(in)  :: outputName
    character(len=*), intent(in)  :: fromWhere, message
    integer, intent(in)           :: iter
    logical, intent(in), optional :: local
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineP='cp_output_handling:cp_write_output'
    logical :: failure
    integer ::unitNr,oErr
    failure=.false.
    if (cp_would_output(logger,outputName,fromWhere,iter,&
         local=local,error=error)) then
       unitNr=cp_unitnr_for_output(logger,outputName,fromWhere,iter,&
            local=local,error=error) 
       write (unit=unitNr,fmt='(a)', iostat=oErr) message
       CPPostcondition(oErr == 0,cp_failure_level,routineP,error,failure)
       call cp_finished_output(logger,outputName,fromWhere,iter,&
            local=local,error=error)
    end if
  end subroutine cp_write_output
  !****************************************************************************

end module cp_output_handling


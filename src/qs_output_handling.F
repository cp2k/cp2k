!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_output_handling [1.0] *
!!
!!   NAME
!!     qs_output_handling
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
!!     handling: @see qs_log_handling
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
module qs_output_handling
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp=>dp
  implicit none

  logical, private, parameter :: debug_this_module=.true.
  
!! outputs the array
  interface qs_output_array
     module procedure qs_output_2d_r_array, qs_output_1d_r_array
  end interface
!!***
  !****************************************************************************
contains

!!****f* qs_output_handling/qs_would_output [1.0] *
!!
!!   NAME
!!     qs_would_output
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
  pure function qs_would_output(logger, outputName , fromWhere, iter)
    logical :: qs_would_output
    type(qs_logger), intent(in) ::logger
!! name of the output, for example 'hamiltonian'
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter

    qs_would_output=.true.
  end function qs_would_output
!!***
  !****************************************************************************

!!****f* qs_output_handling/qs_unitnr_for_output [1.0] *
!!
!!   NAME
!!     f
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
!!            see module qs_error_handling 
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
  pure function qs_unitnr_for_output(logger, outputName, fromWhere, iter)
    integer :: qs_unitnr_for_output
    type(qs_logger), intent(in) :: logger
    character(len=*), intent(in):: outputName, fromWhere
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter

    if (logger%unitNr<=0) then
       qs_unitnr_for_output=qs_std_out_unit_nr
    else
       qs_unitnr_for_output=logger%unitNr
    end if
  end function qs_unitnr_for_output
!!***
  !****************************************************************************

!!****f* qs_output_handling/qs_write_output [1.0] *
!!
!!   NAME
!!     qs_write_output
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
!!            see module qs_error_handling 
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
  subroutine qs_write_output(logger, outputName , fromWhere , iter, message,&
       error)
    type(qs_logger) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere, message
!! the iteration number, or some other number that permit do 
!! differentiate between different instances of the same object
    integer, intent(in) :: iter
    type(qs_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineN='qs_output_handling:qs_write_output'
    logical :: failure
    integer ::unitNr,oErr
    failure=.false.
    if (.not. associated(qs_default_logger)) then
       print *, 'ERROR: qs_default_logger is not associated, did you call qs_init_default_loggers at startup??'
    end if
    if (qs_would_output(logger,outputName,fromWhere,iter)) then
       unitNr=qs_unitnr_for_output(logger,outputName,fromWhere,iter) 
       write (unit=unitNr,fmt='(a)', iostat=oErr) message
       QSPostcondition(oErr == 0,qs_failure_level,routineN,error,failure)
    end if
  end subroutine qs_write_output
!!***
  !****************************************************************************

!!****f* qs_pao_utils/qs_output_2d_r_array [1.0] *
!!
!!   NAME
!!     qs_output_2d_r_array
!!
!!   FUNCTION
!!     writes a 2d array of reals to the output
!!
!!   NOTES
!!     At the moment outputs the comment precede by a # and then the whole
!!     with spaces between the entries (with the default format) and newlines
!!     between the rows of the matrix. 
!!     Useful for matlab, maybe I will change it for Mathematica:
!!     { { array(1,1), array(1,2), ...},
!!       { array(2,1), array(2,2), ...},...}
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere:
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     array: the array to output
!!     rformat: the (optional) format of a single number
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_output_2d_r_array(logger, outputName, fromWhere, iter, array,&
       comment, rformat, error)
    type(qs_logger) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    real(kind=wp), intent(in), dimension(:,:) :: array
    character(len=*), intent(in),optional :: comment, rformat
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_output_2d_array'
    integer :: unitNr, lenFormat, iostat, i
    character(len=*), parameter :: defaultFormat='(es10.3," ")'
    character(len=10) :: nRiga
    failure=.FALSE.

    if (qs_would_output(logger, outputName , fromWhere, iter)) then
       unitNr=qs_unitnr_for_output(logger, outputName, fromWhere, iter)
       nRiga=qs_to_string(size(array,2))
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          QSInvariant(iostat==0,qs_failure_level,routineN,error,failure) ! sets failure to true if it fails
       end if
       do i=1,size(array,1)
          if (present(rformat)) then
             write(unit=unitNr,fmt='(" ",'//nRiga//rFormat//')',iostat=iostat) array(i,:)
          else
             write(unit=unitNr,fmt='(" ",'//nRiga//defaultFormat//')',iostat=iostat) array(i,:)
          end if
          QSInvariant(iostat==0,qs_failure_level,routineN,error,failure) ! sets failure to true if it fails
          if (failure) then
             exit
          end if
       end do
    end if
  end subroutine qs_output_2d_r_array
!!***
  !***************************************************************************

!!****f* qs_pao_utils/qs_output_1d_r_array [1.0] *
!!
!!   NAME
!!     qs_output_1d_r_array
!!
!!   FUNCTION
!!     writes a 1d array of reals to the output
!!
!!   NOTES
!!     At the moment outputs the comment precede by a # and then the whole
!!     with newlines between the entries, i.e. as a (n x 1) matrix.
!!     Useful for matlab, maybe I will change it for Mathematica:
!!     { array(1,1), array(1,2), ...}
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     array: the array to output
!!     rformat: the (optional) format of a single number
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_output_1d_r_array(logger, outputName, fromWhere, iter, array,&
       comment, rformat, error)
    type(qs_logger) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    real(kind=wp), intent(in), dimension(:) :: array
    character(len=*), intent(in),optional :: comment, rformat
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_output_2d_array'
    integer :: unitNr, lenFormat, iostat, i
    character(len=*), parameter :: defaultFormat='(es10.3)'
    failure=.FALSE.

    if (qs_would_output(logger, outputName , fromWhere, iter)) then
       unitNr=qs_unitnr_for_output(logger, outputName, fromWhere, iter)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          QSInvariant(iostat==0,qs_failure_level,routineN,error,failure)
       end if
       ! do i=1,size(array)
       if (present(rformat)) then
          write(unit=unitNr,fmt=rFormat,iostat=iostat) array(:)
       else
          write(unit=unitNr,fmt=defaultFormat,iostat=iostat) array(:)
       end if
       QSInvariant(iostat==0,qs_failure_level,routineN,error,failure)
       !   if (failure) then
       !     exit
       !   end if
       ! end do
    end if
  end subroutine qs_output_1d_r_array
!!***
  !***************************************************************************

end module qs_output_handling


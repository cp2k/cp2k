!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_array_i_utils [1.0] *
!!
!!   NAME
!!     cp_array_i_utils
!!
!!   FUNCTION
!!     various utilities that regard array of different kinds:
!!     output, allocation,...
!!
!!   NOTES
!!     maybe it is not a good idea mixing output and memeory utils...
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     3.2002 templatized [fawzi]
!!
!!   SOURCE
!****************************************************************************
module cp_array_i_utils
  use kinds, only: wp=>dp
  use cp_log_handling
  use cp_error_handling
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_output_handling
  
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_array_i_utils'

  ! generic interfaces
  public :: cp_output, cp_guarantee_size, cp_write
  
  ! the underlying procedures
  public :: cp_1d_i_output, cp_2d_i_output,&
       cp_1d_i_guarantee_size, cp_i_write,&
       cp_1d_i_write, cp_2d_i_write

  interface cp_output
     module procedure cp_1d_i_output, &
          cp_2d_i_output
  end interface
  interface cp_guarantee_size
     module procedure cp_1d_i_guarantee_size
  end interface
  interface cp_write
     module procedure cp_i_write, cp_1d_i_write,&
       cp_2d_i_write
  end interface
!***
!****************************************************************************
contains

!!****f* cp_array_i_utils/cp_i_write *
!!
!!   NAME
!!     cp_i_write
!!
!!   SYNOPSIS
!!     Subroutine cp_i_write(element, unit_nr, el_format, error)
!!       Integer, Intent (IN):: element
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_i_write
!!
!!   FUNCTION
!!     writes an element to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     element: the element to write
!!     unit_nr: the unit to write to (defaults to the standard out)
!!     el_format: the format of a single element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine cp_i_write(element, unit_nr, el_format, error)
  integer, intent(in) :: element
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat
  character(len=*), parameter :: routineN='cp_i_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(i6)"
  failure=.false.
  
  if (present(el_format)) then
     write(unit=unit_nr,fmt=el_format,iostat=iostat) element
  else
     write(unit=unit_nr,fmt=defaultFormat,iostat=iostat) element
  end if
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
end subroutine cp_i_write
!***************************************************************************


!!****f* cp_array_i_utils/cp_1d_i_output [1.0] *
!!
!!   NAME
!!     cp_1d_i_output
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_i_output(logger, outputname, fromwhere, iter, array,&
!!         comment, el_format, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Integer, Intent (IN), Dimension(:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_i_output
!!
!!   FUNCTION
!!     writes a 1d array to the output
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
!!     comment: an (optional) comment that is written before the array
!!     el_format: the (optional) format of a single number
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     3.2002 templatized [fawzi]
!!     4.2002 separated write routine [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_1d_i_output(logger, outputName,&
       fromWhere, iter, array,&
       comment, el_format, error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    integer, intent(in), dimension(:) :: array
    character(len=*), intent(in),optional :: comment, el_format
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_1d_i_output',&
         routineP=moduleN//':'//routineN
    integer :: handle, unitNr, lenFormat, iostat, i
    failure=.false.

    call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_would_output(logger, outputName , fromWhere, iter)) then
       unitNr=cp_unitnr_for_output(logger, outputName, fromWhere, iter)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
       end if
       call cp_1d_i_write(array,unit_nr=unitNr,&
            el_format=el_format,error=error)
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_1d_i_output
!***************************************************************************

!!****f* cp_array_i_utils/cp_1d_i_write [1.0] *
!!
!!   NAME
!!     cp_1d_i_write
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_i_write(array, unit_nr, el_format, error)
!!       Integer, Dimension(:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_i_write
!!
!!   FUNCTION
!!     writes an array to the given unit
!!
!!   NOTES
!!     maybe I will move to a comma separated paretized list
!!
!!   INPUTS
!!     array: the array to write
!!     unit_nr: the unit to write to (defaults to the standard out)
!!     el_format: the format of a single element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine cp_1d_i_write(array, unit_nr, el_format, error)
  integer, dimension(:), intent(in) :: array
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat
  character(len=*), parameter :: routineN='cp_1d_i_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(i6)"
  
  ! do i=1,size(array)
  if (present(el_format)) then
     write(unit=unit_nr,fmt=el_format,iostat=iostat) array(:)
  else
     write(unit=unit_nr,fmt=defaultFormat,iostat=iostat) array(:)
  end if
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  !   if (failure) then
  !     exit
  !   end if
  ! end do
end subroutine cp_1d_i_write
!***************************************************************************

!!****f* cp_array_i_utils/cp_2d_i_output *
!!
!!   NAME
!!     cp_2d_i_output
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_i_output(logger, outputname, fromwhere, iter, array,&
!!         comment, el_format, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Integer, Intent (IN), Dimension(:,:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_i_output
!!
!!   FUNCTION
!!     writes a 2d array to the output
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
!!     el_format: the (optional) format of a single number
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     3.2002 templatized [fawzi]
!!     4.2002 separated write routine [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_2d_i_output(logger, outputName, &
       fromWhere, iter, array,&
       comment, el_format, error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    integer, intent(in), dimension(:,:) :: array
    character(len=*), intent(in),optional :: comment, el_format
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_2d_i_output',&
         routineP=moduleN//':'//routineN
    integer :: handle, unitNr, lenFormat, iostat
    failure=.FALSE.

    call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_would_output(logger, outputName , fromWhere, iter)) then
       unitNr=cp_unitnr_for_output(logger, outputName, fromWhere, iter)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_2d_i_write(array,unit_nr=unitNr,&
               el_format=el_format,error=error)
       end if
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_2d_i_output
!***************************************************************************

!!****f* cp_array_i_utils/cp_2d_i_write *
!!
!!   NAME
!!     cp_2d_i_write
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_i_write(array, unit_nr, el_format, error)
!!       Integer, Dimension(:,:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_i_write
!!
!!   FUNCTION
!!     writes an array to the given unit
!!
!!   NOTES
!!     maybe I will move to a comma separated parentized list
!!
!!   INPUTS
!!     array: the array to write
!!     unit_nr: the unit to write to (defaults to the standard out)
!!     el_format: the format of a single element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine cp_2d_i_write(array, unit_nr, el_format, error)
  integer, dimension(:,:), intent(in) :: array
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat, i
  character(len=*), parameter :: routineN='cp_2d_i_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(i6)"
  character(len=10) :: nRiga
  failure=.false.
  
  nRiga=cp_to_string(size(array,2))
  do i=1,size(array,1)
     if (present(el_format)) then
        write(unit=unit_nr,fmt='(" ",'//nRiga//el_format//')',iostat=iostat) array(i,:)
     else
        write(unit=unit_nr,fmt='(" ",'//nRiga//defaultFormat//')',iostat=iostat) array(i,:)
     end if
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     if (failure) then
        exit
     end if
  end do
end subroutine cp_2d_i_write
!***************************************************************************

!!****f* cp_array_i_utils/cp_1d_i_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_1d_i_guarantee_size
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_i_guarantee_size(array, n, init, error)
!!       Integer, Dimension(:), Pointer:: array
!!       Integer, Intent (IN):: n
!!       Integer, Intent (IN), Optional:: init
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_i_guarantee_size
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     Issues a warning when the size changes (but not on allocation
!!     and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     3.2002 templatized [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_1d_i_guarantee_size(array, n, init, error)
    integer, dimension(:), pointer :: array
    integer, intent(in) :: n
    integer, intent(in), optional :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_1d_i_guarantee_size',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             deallocate(array, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine cp_1d_i_guarantee_size
!***************************************************************************

  ! template def put here so that line numbers in template and derived 
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,use,defaultFormatType1)]
! ARGS:
!  type1 = "integer"
!  nametype1 = "i"
!  use = ""
!  defaultFormatType1 = ""(i6)""


end module cp_array_i_utils

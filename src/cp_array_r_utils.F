!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_array_r_utils [1.0] *
!!
!!   NAME
!!     cp_array_r_utils
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
module cp_array_r_utils
  use kinds, only: wp=>dp
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level, cp_log, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message,cp_assertion_failed, cp_internal_error,&
       cp_error_init, cp_error_dealloc_ref
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_output_handling
  use qs_parser, only: open_file, close_file
  
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_array_r_utils'

  ! generic interfaces
  public :: cp_output, cp_guarantee_size, cp_write
  
  ! the underlying procedures
  public :: cp_1d_r_output, cp_2d_r_output,&
       cp_1d_r_guarantee_size, cp_r_write,&
       cp_1d_r_write, cp_2d_r_write,&
       cp_2d_r_guarantee_size

  interface cp_output
     module procedure cp_1d_r_output, &
          cp_2d_r_output
  end interface
  interface cp_guarantee_size
     module procedure cp_1d_r_guarantee_size,&
          cp_2d_r_guarantee_size
  end interface
  interface cp_write
     module procedure cp_r_write, cp_1d_r_write,&
       cp_2d_r_write
  end interface
!***
!****************************************************************************
contains

!!****f* cp_array_r_utils/cp_r_write *
!!
!!   NAME
!!     cp_r_write
!!
!!   SYNOPSIS
!!     Subroutine cp_r_write(element, unit_nr, el_format, error)
!!       Real(Kind=wp), Intent (IN):: element
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_r_write
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
subroutine cp_r_write(element, unit_nr, el_format, error)
  real(kind=wp), intent(in) :: element
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat
  character(len=*), parameter :: routineN='cp_r_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(es14.6)"
  failure=.false.
  
  if (present(el_format)) then
     write(unit=unit_nr,fmt=el_format,iostat=iostat) element
  else
     write(unit=unit_nr,fmt=defaultFormat,iostat=iostat) element
  end if
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
end subroutine cp_r_write
!***************************************************************************


!!****f* cp_array_r_utils/cp_1d_r_output [1.0] *
!!
!!   NAME
!!     cp_1d_r_output
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_r_output(logger, outputname, fromwhere, iter, array,&
!!         comment, el_format, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Real(Kind=wp), Intent (IN), Dimension(:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_r_output
!!
!!   FUNCTION
!!     writes a 1d array to an aptly named file
!!
!!   NOTES
!!     The filename is written as note to the log.
!!     At the moment outputs the comment precede by a # and then the whole
!!     with newlines between the entries, i.e. as a (n x 1) matrix.
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
!!     6.2002 write to a separate file [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_1d_r_output(logger, outputName,&
       fromWhere, iter, array,&
       comment, el_format, error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    real(kind=wp), intent(in), dimension(:) :: array
    character(len=*), intent(in),optional :: comment, el_format
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, exists
    character(len=*), parameter :: routineN='cp_1d_r_output',&
         routineP=moduleN//':'//routineN
    character(len=60) :: filename
    integer :: handle, unitNr, lenFormat, iostat, i
    failure=.false.

    call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_would_output(logger, outputName , fromWhere, iter)) then
       filename=trim(outputName)//'_'//trim(adjustl(cp_to_string(iter)))//'.dat'
       call cp_log(logger, cp_note_level, fromWhere, "writing output to file'"//&
            trim(filename)) 
       inquire (FILE=trim(filename),EXIST=exists)
       call cp_assert(.not.exists, cp_warning_level,cp_assertion_failed,&
            fromWhere=routineP,message=" overwriting file '"//filename//"'",&
            error=error)
       call open_file(trim(filename),file_status="unknown",&
            file_action="write",&
            unit_number=unitNr)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_1d_r_write(array,unit_nr=unitNr,&
               el_format=el_format,error=error)
       end if
       call close_file(unitNr)
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_1d_r_output
!***************************************************************************

!!****f* cp_array_r_utils/cp_1d_r_write [1.0] *
!!
!!   NAME
!!     cp_1d_r_write
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_r_write(array, unit_nr, el_format, error)
!!       Real(Kind=wp), Dimension(:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_r_write
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
subroutine cp_1d_r_write(array, unit_nr, el_format, error)
  real(kind=wp), dimension(:), intent(in) :: array
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat
  character(len=*), parameter :: routineN='cp_1d_r_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(es14.6)"
  
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
end subroutine cp_1d_r_write
!***************************************************************************

!!****f* cp_array_r_utils/cp_2d_r_output *
!!
!!   NAME
!!     cp_2d_r_output
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_r_output(logger, outputname, fromwhere, iter, array,&
!!         comment, el_format, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Real(Kind=wp), Intent (IN), Dimension(:,:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_r_output
!!
!!   FUNCTION
!!     writes a 2d array to an aptly named file
!!
!!   NOTES
!!     The filename is written as note to the log.
!!     At the moment outputs the comment precede by a # and then the whole
!!     with spaces between the entries (with the default format) and newlines
!!     between the rows of the matrix. 
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
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
!!     6.2002 write to a separate file [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_2d_r_output(logger, outputName, &
       fromWhere, iter, array,&
       comment, el_format, error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    real(kind=wp), intent(in), dimension(:,:) :: array
    character(len=*), intent(in),optional :: comment, el_format
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, exists
    character(len=*), parameter :: routineN='cp_2d_r_output',&
         routineP=moduleN//':'//routineN
    character(len=60) :: filename
    integer :: handle, unitNr, lenFormat, iostat
    failure=.FALSE.

    call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_would_output(logger, outputName , fromWhere, iter)) then
       filename=trim(outputName)//'_'//trim(adjustl(cp_to_string(iter)))//'.dat'
       call cp_log(logger, cp_note_level, fromWhere, "writing output to file'"//&
            trim(filename)) 
       inquire (FILE=trim(filename),EXIST=exists)
       call cp_assert(.not.exists, cp_warning_level,cp_assertion_failed,&
            fromWhere=routineP,message=" overwriting file '"//filename//"'",&
            error=error)
       call open_file(trim(filename),file_status="unknown",&
            file_action="write",&
            unit_number=unitNr)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_2d_r_write(array,unit_nr=unitNr,&
               el_format=el_format,error=error)
       end if
       call close_file(unitNr)
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_2d_r_output
!***************************************************************************

!!****f* cp_array_r_utils/cp_2d_r_write *
!!
!!   NAME
!!     cp_2d_r_write
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_r_write(array, unit_nr, el_format, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_r_write
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
subroutine cp_2d_r_write(array, unit_nr, el_format, error)
  real(kind=wp), dimension(:,:), intent(in) :: array
  integer, intent(in) :: unit_nr
  character(len=*), intent(in), optional :: el_format
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: iostat, i
  character(len=*), parameter :: routineN='cp_2d_r_write',&
        routineP=moduleN//':'//routineN
  character(len=*), parameter :: defaultFormat="(es14.6)"
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
end subroutine cp_2d_r_write
!***************************************************************************

!!****f* cp_array_r_utils/cp_1d_r_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_1d_r_guarantee_size
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_r_guarantee_size(array, n, init, error)
!!       Real(Kind=wp), Dimension(:), Pointer:: array
!!       Integer, Intent (IN):: n
!!       Real(Kind=wp), Intent (IN), Optional:: init
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_r_guarantee_size
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
  subroutine cp_1d_r_guarantee_size(array, n, init, error)
    real(kind=wp), dimension(:), pointer :: array
    integer, intent(in) :: n
    real(kind=wp), intent(in), optional :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_1d_r_guarantee_size',&
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
  end subroutine cp_1d_r_guarantee_size
!***************************************************************************

!!****f* cp_array_r_utils/cp_2d_r_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_2d_r_guarantee_size
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_r_guarantee_size(array, n_rows, n_cols, init,&
!!         error)
!!       Real(Kind=wp), Dimension(:,:), Pointer:: array
!!       Integer, Intent (IN):: n_rows, n_cols
!!       Real(Kind=wp), Intent (IN), Optional:: init
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_r_guarantee_size
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
!!     - array: the array to reallocate if necessary
!!     - n_rows: the wanted number of rows
!!     - n_cols: the wanted number of cols
!!     - init: the initial value of the elements when resized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     5.2001 first version [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_2d_r_guarantee_size(array, n_rows, n_cols ,&
       init, error)
    real(kind=wp), dimension(:,:), pointer :: array
    integer, intent(in) :: n_rows, n_cols
    real(kind=wp), intent(in), optional :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_2d_r_guarantee_size',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(n_cols>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(n_rows>=0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array,1) /= n_rows .or. size(array,2) /= n_cols) then
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             deallocate(array, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             nullify(array)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n_rows,n_cols), stat=stat)
          CPPostconditionNoFail(stat==0, cp_failure_level,routineP,error)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine cp_2d_r_guarantee_size
!***************************************************************************

  ! template def put here so that line numbers in template and derived 
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,use,defaultFormatType1)]
! ARGS:
!  type1 = "real(kind=wp)"
!  nametype1 = "r"
!  use = ""
!  defaultFormatType1 = ""(es14.6)""


end module cp_array_r_utils

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_array_logical_utils [1.0] *
!!
!!   NAME
!!     cp_array_logical_utils
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
MODULE cp_array_logical_utils
  USE kinds, ONLY: wp=>dp
  USE cp_log_handling, ONLY: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level, cp_log, cp_logger_type, cp_logger_generate_filename
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_assert, &
       cp_error_message,cp_assertion_failed, cp_internal_error,&
       cp_error_init, cp_error_dealloc_ref
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE cp_output_handling
  USE qs_parser, ONLY: open_file, close_file
  [USE]
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_array_logical_utils'

  ! generic interfaces
  PUBLIC :: cp_output, cp_guarantee_size, cp_write
  
  ! the underlying procedures
  PUBLIC :: cp_1d_logical_output, cp_2d_logical_output,&
       cp_1d_logical_guarantee_size, cp_logical_write,&
       cp_1d_logical_write, cp_2d_logical_write,&
       cp_2d_logical_guarantee_size

  INTERFACE cp_output
     MODULE PROCEDURE cp_1d_logical_output, &
          cp_2d_logical_output
  END INTERFACE
  INTERFACE cp_guarantee_size
     MODULE PROCEDURE cp_1d_logical_guarantee_size,&
          cp_2d_logical_guarantee_size
  END INTERFACE
  INTERFACE cp_write
     MODULE PROCEDURE cp_logical_write, cp_1d_logical_write,&
       cp_2d_logical_write
  END INTERFACE
!***
!****************************************************************************
CONTAINS

!!****f* cp_array_logical_utils/cp_logical_write *
!!
!!   NAME
!!     cp_logical_write
!!
!!   SYNOPSIS
!!     Subroutine cp_logical_write(element, unit_nr, el_format, error)
!!       Logical, Intent (IN):: element
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_logical_write
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
SUBROUTINE cp_logical_write(element, unit_nr, el_format, error)
  logical, INTENT(in) :: element
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat
  CHARACTER(len=*), PARAMETER :: routineN='cp_logical_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(t1)"
  failure=.FALSE.
  
  IF (PRESENT(el_format)) THEN
     WRITE(unit=unit_nr,fmt=el_format,iostat=iostat) element
  ELSE
     WRITE(unit=unit_nr,fmt=defaultFormat,iostat=iostat) element
  END IF
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
END SUBROUTINE cp_logical_write
!***************************************************************************


!!****f* cp_array_logical_utils/cp_1d_logical_output [1.0] *
!!
!!   NAME
!!     cp_1d_logical_output
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_logical_output(logger, outputname, fromwhere, iter,&
!!         array, comment, local, el_format, error)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Logical, Intent (IN), Dimension(:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_logical_output
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
!!     - logger: the logger that decides where the logging should go
!!     - outputName: the name of the output, for example 'hamiltonian'
!!     - fromWhere: string of the form module:function or file:lineNr
!!       that says where the error happend
!!     - iter: the iteration number, or some other number that permit to 
!!       differentiate between different instances of the same object
!!     - array: the array to output
!!     - comment: an (optional) comment that is written before the array
!!     - local: if true each task writes its own file, otherwise (the default)
!!       only the source writes its data
!!     - el_format: the (optional) format of a single number
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     03.2002 templatized [fawzi]
!!     04.2002 separated write routine [fawzi]
!!     06.2002 write to a separate file [fawzi]
!!     08.2002 added local/global output [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_1d_logical_output(logger, outputName,&
     fromWhere, iter, array,&
     comment, local, el_format, error)
  TYPE(cp_logger_type), POINTER ::logger
  CHARACTER(len=*) , INTENT(in) :: outputName
  CHARACTER(len=*), INTENT(in):: fromWhere
  INTEGER, INTENT(in) :: iter
  logical, INTENT(in), DIMENSION(:) :: array
  CHARACTER(len=*), INTENT(in),OPTIONAL :: comment, el_format
  LOGICAL, INTENT(in), OPTIONAL :: local
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure, exists, loc
  CHARACTER(len=*), PARAMETER :: routineN='cp_1d_logical_output',&
       routineP=moduleN//':'//routineN
  CHARACTER(len=60) :: filename
  INTEGER :: handle, unitNr, lenFormat, iostat, i
  failure=.FALSE.; loc=.FALSE.

  CALL timeset(routineN//','//moduleN,'I',"",handle)
  IF (PRESENT(local)) loc=local
  IF (loc .OR. logger%para_env%mepos==logger%para_env%source) THEN
     IF (cp_would_output(logger, outputName , fromWhere, iter,&
          local=loc,error=error)) THEN

        CALL cp_logger_generate_filename(logger,filename,&
             TRIM(outputName)//'-'//TRIM(ADJUSTL(cp_to_string(iter))),&
             '.dat',local=loc)
        CALL cp_log(logger, cp_note_level, fromWhere,&
             "writing output to file'"//TRIM(filename)//"'") 
        INQUIRE (FILE=TRIM(filename),EXIST=exists)
        CALL cp_assert(.NOT.exists, cp_warning_level,cp_assertion_failed,&
             fromWhere=routineP,message=" overwriting file '"//filename//"'",&
             error=error)
        CALL open_file(TRIM(filename),file_status="unknown",&
             file_action="write",&
             unit_number=unitNr)

        IF (PRESENT(comment)) THEN
           WRITE(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
           CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           CALL cp_1d_logical_write(array,unit_nr=unitNr,&
                el_format=el_format,error=error)
        END IF

        CALL close_file(unitNr)
     END IF
  END IF
  CALL timestop(0.0_wp,handle)
END SUBROUTINE cp_1d_logical_output
!***************************************************************************

!!****f* cp_array_logical_utils/cp_1d_logical_write [1.0] *
!!
!!   NAME
!!     cp_1d_logical_write
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_logical_write(array, unit_nr, el_format, error)
!!       Logical, Dimension(:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_logical_write
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
SUBROUTINE cp_1d_logical_write(array, unit_nr, el_format, error)
  logical, DIMENSION(:), INTENT(in) :: array
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat
  CHARACTER(len=*), PARAMETER :: routineN='cp_1d_logical_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(t1)"
  
  ! do i=1,size(array)
  IF (PRESENT(el_format)) THEN
     WRITE(unit=unit_nr,fmt=el_format,iostat=iostat) array(:)
  ELSE
     WRITE(unit=unit_nr,fmt=defaultFormat,iostat=iostat) array(:)
  END IF
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  !   if (failure) then
  !     exit
  !   end if
  ! end do
END SUBROUTINE cp_1d_logical_write
!***************************************************************************

!!****f* cp_array_logical_utils/cp_2d_logical_output *
!!
!!   NAME
!!     cp_2d_logical_output
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_logical_output(logger, outputname, fromwhere, iter,&
!!         array, comment, local, el_format, error)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Logical, Intent (IN), Dimension(:,:):: array
!!       Character(Len=*), Intent (IN), Optional:: comment, el_format
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_logical_output
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
!!     local: if the output is local to each task
!!     el_format: the (optional) format of a single number
!!     local: if the output is local to this task (defaults to false)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 first version [fawzi]
!!     03.2002 templatized [fawzi]
!!     04.2002 separated write routine [fawzi]
!!     06.2002 write to a separate file [fawzi]
!!     08.2002 added local/global output [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_2d_logical_output(logger, outputName, &
     fromWhere, iter, array,&
     comment, local, el_format, error)
  TYPE(cp_logger_type), POINTER ::logger
  CHARACTER(len=*) , INTENT(in) :: outputName
  CHARACTER(len=*), INTENT(in):: fromWhere
  INTEGER, INTENT(in) :: iter
  logical, INTENT(in), DIMENSION(:,:) :: array
  CHARACTER(len=*), INTENT(in),OPTIONAL :: comment, el_format
  LOGICAL, INTENT(in), OPTIONAL :: local
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure, exists, loc
  CHARACTER(len=*), PARAMETER :: routineN='cp_2d_logical_output',&
       routineP=moduleN//':'//routineN
  CHARACTER(len=60) :: filename
  INTEGER :: handle, unitNr, lenFormat, iostat
  failure=.FALSE.; loc=.FALSE.

  CALL timeset(routineN//','//moduleN,'I',"",handle)
  IF (PRESENT(local)) loc=local
  IF (loc .OR. logger%para_env%mepos==logger%para_env%source) THEN
     IF (cp_would_output(logger, outputName , fromWhere, iter,&
          local=loc,error=error)) THEN

        CALL cp_logger_generate_filename(logger,filename,&
             TRIM(outputName)//'-'//TRIM(ADJUSTL(cp_to_string(iter))),&
             '.dat',local=loc)
        CALL cp_log(logger, cp_note_level, fromWhere, &
             "writing output to file'"//TRIM(filename)//"'") 
        INQUIRE (FILE=TRIM(filename),EXIST=exists)
        CALL cp_assert(.NOT.exists, cp_warning_level,cp_assertion_failed,&
             fromWhere=routineP,message=" overwriting file '"//filename//"'",&
             error=error)
        CALL open_file(TRIM(filename),file_status="unknown",&
             file_action="write",&
             unit_number=unitNr)

        IF (PRESENT(comment)) THEN
           WRITE(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
           CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           CALL cp_2d_logical_write(array,unit_nr=unitNr,&
                el_format=el_format,error=error)
        END IF

        CALL close_file(unitNr)
     END IF
  END IF
  CALL timestop(0.0_wp,handle)
END SUBROUTINE cp_2d_logical_output
!***************************************************************************

!!****f* cp_array_logical_utils/cp_2d_logical_write *
!!
!!   NAME
!!     cp_2d_logical_write
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_logical_write(array, unit_nr, el_format, error)
!!       Logical, Dimension(:,:), Intent (IN):: array
!!       Integer, Intent (IN):: unit_nr
!!       Character(Len=*), Intent (IN), Optional:: el_format
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_logical_write
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
SUBROUTINE cp_2d_logical_write(array, unit_nr, el_format, error)
  logical, DIMENSION(:,:), INTENT(in) :: array
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat, i
  CHARACTER(len=*), PARAMETER :: routineN='cp_2d_logical_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(t1)"
  CHARACTER(len=10) :: nRiga
  failure=.FALSE.
  
  nRiga=cp_to_string(SIZE(array,2))
  DO i=1,SIZE(array,1)
     IF (PRESENT(el_format)) THEN
        WRITE(unit=unit_nr,fmt='(" ",'//nRiga//el_format//')',iostat=iostat) array(i,:)
     ELSE
        WRITE(unit=unit_nr,fmt='(" ",'//nRiga//defaultFormat//')',iostat=iostat) array(i,:)
     END IF
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     IF (failure) THEN
        EXIT
     END IF
  END DO
END SUBROUTINE cp_2d_logical_write
!***************************************************************************

!!****f* cp_array_logical_utils/cp_1d_logical_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_1d_logical_guarantee_size
!!
!!   SYNOPSIS
!!     Subroutine cp_1d_logical_guarantee_size(array, n, init, error)
!!       Logical, Dimension(:), Pointer:: array
!!       Integer, Intent (IN):: n
!!       Logical, Intent (IN), Optional:: init
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_1d_logical_guarantee_size
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
  SUBROUTINE cp_1d_logical_guarantee_size(array, n, init, error)
    logical, DIMENSION(:), POINTER :: array
    INTEGER, INTENT(in) :: n
    logical, INTENT(in), OPTIONAL :: init
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='cp_1d_logical_guarantee_size',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (ASSOCIATED(array)) THEN
          IF (SIZE(array) /= n) THEN
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             DEALLOCATE(array, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
       END IF
       IF (.NOT.ASSOCIATED(array)) THEN
          ALLOCATE(array(n), stat=stat)
          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
          IF (PRESENT(init).AND..NOT.failure) array=init
       END IF
    END IF failureIf
  END SUBROUTINE cp_1d_logical_guarantee_size
!***************************************************************************

!!****f* cp_array_logical_utils/cp_2d_logical_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_2d_logical_guarantee_size
!!
!!   SYNOPSIS
!!     Subroutine cp_2d_logical_guarantee_size(array, n_rows, n_cols, init,&
!!         error)
!!       Logical, Dimension(:,:), Pointer:: array
!!       Integer, Intent (IN):: n_rows, n_cols
!!       Logical, Intent (IN), Optional:: init
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_2d_logical_guarantee_size
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
  SUBROUTINE cp_2d_logical_guarantee_size(array, n_rows, n_cols ,&
       init, error)
    logical, DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(in) :: n_rows, n_cols
    logical, INTENT(in), OPTIONAL :: init
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='cp_2d_logical_guarantee_size',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(n_cols>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(n_rows>=0,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (ASSOCIATED(array)) THEN
          IF (SIZE(array,1) /= n_rows .OR. SIZE(array,2) /= n_cols) THEN
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             DEALLOCATE(array, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY(array)
          END IF
       END IF
       IF (.NOT.ASSOCIATED(array)) THEN
          ALLOCATE(array(n_rows,n_cols), stat=stat)
          CPPostconditionNoFail(stat==0, cp_failure_level,routineP,error)
          IF (PRESENT(init).AND..NOT.failure) array=init
       END IF
    END IF failureIf
  END SUBROUTINE cp_2d_logical_guarantee_size
!***************************************************************************

  ! template def put here so that line numbers in template and derived 
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,use,defaultFormatType1)]
! ARGS:
!  type1 = "logical"
!  nametype1 = "logical"
!  use = ""
!  defaultFormatType1 = ""(t1)""


END MODULE cp_array_logical_utils

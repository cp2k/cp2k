!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"
#define CP_R_LESS_Q(el1,el2) ( el1 < el2 )

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
!!     This is a template
!!
!!     **** DO NOT MODIFY THE .F FILES ****
!!     modify the .template and .instantition files
!!
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
MODULE cp_array_r_utils
  USE kinds, ONLY: dp
  USE cp_log_handling, ONLY: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level, cp_log, cp_logger_type, cp_logger_generate_filename
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_assert, cp_a_l, &
       cp_error_message,cp_assertion_failed, cp_internal_error,&
       cp_error_init, cp_error_dealloc_ref
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE qs_parser, ONLY: open_file, close_file
  USE machine, ONLY: m_flush
  
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_array_r_utils'

  ! types
  PUBLIC :: cp_1d_r_p_type, cp_2d_r_p_type,&
       cp_3d_r_p_type

  ! generic interfaces
  PUBLIC :: cp_guarantee_size, cp_write
  
  ! the underlying procedures
  PUBLIC :: cp_1d_r_guarantee_size, cp_r_write,&
       cp_1d_r_write, cp_2d_r_write,&
       cp_3d_r_write,&
       cp_2d_r_guarantee_size,&
       cp_1d_r_bsearch

  INTERFACE cp_guarantee_size
     MODULE PROCEDURE cp_1d_r_guarantee_size,&
          cp_2d_r_guarantee_size
  END INTERFACE
  INTERFACE cp_write
     MODULE PROCEDURE cp_r_write, cp_1d_r_write,&
       cp_2d_r_write
  END INTERFACE
!***
!****************************************************************************
  
  !!****s* cp_array_r_utils/cp_1d_r_p_type [1.0] *
  !!
  !!   NAME
  !!     cp_1d_r_p_type
  !!
  !!   FUNCTION
  !!     represent a pointer to a 1d array
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     02.2003 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type cp_1d_r_p_type
     real(kind=dp), dimension(:), pointer :: array
  end type cp_1d_r_p_type
  !!***
  !****************************************************************************

  !!****s* cp_array_r_utils/cp_2d_r_p_type [1.0] *
  !!
  !!   NAME
  !!     cp_2d_r_p_type
  !!
  !!   FUNCTION
  !!     represent a pointer to a 2d array
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     02.2003 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type cp_2d_r_p_type
     real(kind=dp), dimension(:,:), pointer :: array
  end type cp_2d_r_p_type
  !!***
  !****************************************************************************

  !!****s* cp_array_r_utils/cp_3d_r_p_type [1.0] *
  !!
  !!   NAME
  !!     cp_3d_r_p_type
  !!
  !!   FUNCTION
  !!     represent a pointer to a 3d array
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     02.2003 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type cp_3d_r_p_type
     real(kind=dp), dimension(:,:,:), pointer :: array
  end type cp_3d_r_p_type
  !!***
  !****************************************************************************

CONTAINS

!!****f* cp_array_r_utils/cp_r_write *
!!
!!   NAME
!!     cp_r_write
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
SUBROUTINE cp_r_write(element, unit_nr, el_format, error)
  real(kind=dp), INTENT(in) :: element
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat
  CHARACTER(len=*), PARAMETER :: routineN='cp_r_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(es14.6)"
  failure=.FALSE.
  
  IF (PRESENT(el_format)) THEN
     WRITE(unit=unit_nr,fmt=el_format,iostat=iostat) element
  ELSE
     WRITE(unit=unit_nr,fmt=defaultFormat,iostat=iostat) element
  END IF
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
END SUBROUTINE cp_r_write
!***************************************************************************


!!****f* cp_array_r_utils/cp_1d_r_write [1.0] *
!!
!!   NAME
!!     cp_1d_r_write
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
SUBROUTINE cp_1d_r_write(array, unit_nr, el_format, error)
  real(kind=dp), DIMENSION(:), INTENT(in) :: array
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat, i
  CHARACTER(len=*), PARAMETER :: routineN='cp_1d_r_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(es14.6)"
  
  failure=.false.
  WRITE(unit=unit_nr,fmt="('( ')",advance="no", iostat=iostat)
  CPPrecondition(iostat==0,cp_failure_level,routineP,error,failure)  
  IF (.NOT.failure) THEN
     IF (PRESENT(el_format)) THEN
        DO i=1,SIZE(array)-1
           WRITE(unit=unit_nr,fmt=el_format,advance="no") array(i)
           IF (MOD(i,5).EQ.0) THEN  ! only a few elements per line 
              WRITE(unit=unit_nr,fmt="(',')")
           ELSE
              WRITE(unit=unit_nr,fmt="(',')",advance="no")
           ENDIF
        END DO
        IF (SIZE(array)>0) &
             WRITE(unit=unit_nr,fmt=el_format,advance="no") array(SIZE(array))
     ELSE
        DO i=1,SIZE(array)-1
           WRITE(unit=unit_nr,fmt=defaultFormat,advance="no") array(i)
           IF (MOD(i,5).EQ.0) THEN  ! only a few elements per line 
              WRITE(unit=unit_nr,fmt="(',')")
           ELSE
              WRITE(unit=unit_nr,fmt="(',')",advance="no")
           ENDIF
        END DO
        IF (SIZE(array)>0) &
             WRITE(unit=unit_nr,fmt=defaultFormat,advance="no") array(SIZE(array))
     END IF
     WRITE(unit=unit_nr,fmt="(' )')")
  END IF
  call m_flush(unit_nr)

END SUBROUTINE cp_1d_r_write
!***************************************************************************

!!****f* cp_array_r_utils/cp_2d_r_write *
!!
!!   NAME
!!     cp_2d_r_write
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
SUBROUTINE cp_2d_r_write(array, unit_nr, el_format, error)
  real(kind=dp), DIMENSION(:,:), INTENT(in) :: array
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat, i
  CHARACTER(len=*), PARAMETER :: routineN='cp_2d_r_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(es14.6)"
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
  call m_flush(unit_nr)
END SUBROUTINE cp_2d_r_write
!***************************************************************************


!!****f* cp_array_r_utils/cp_3d_r_write *
!!
!!   NAME
!!     cp_3d_r_write
!!
!!   FUNCTION
!!     writes an array to the given unit
!!
!!   NOTES
!!     maybe I will move to a comma separated parentized list
!!     the order of the elements is different from 2d
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
!!     05.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_3d_r_write(array, unit_nr, el_format, error)
  real(kind=dp), DIMENSION(:,:,:), INTENT(in) :: array
  INTEGER, INTENT(in) :: unit_nr
  CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: iostat, j, k
  CHARACTER(len=*), PARAMETER :: routineN='cp_3d_r_write',&
        routineP=moduleN//':'//routineN
  CHARACTER(len=*), PARAMETER :: defaultFormat="(es14.6)"
  CHARACTER(len=10) :: nRiga
  failure=.FALSE.
  
  WRITE(unit=unit_nr,fmt='(i10," ",i10," ",i10,/)',iostat=iostat) &
       SIZE(array,1),SIZE(array,2), SIZE(array,3)
  CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)

  nRiga=cp_to_string(SIZE(array,1))
  DO k=1,SIZE(array,3)
     DO j=1,SIZE(array,2)
        IF (PRESENT(el_format)) THEN
           WRITE(unit=unit_nr,fmt='(" ",'//nRiga//el_format//')',&
                iostat=iostat) array(:,j,k)
        ELSE
           WRITE(unit=unit_nr,fmt='(" ",'//nRiga//defaultFormat//')',&
                iostat=iostat) array(:,j,k)
        END IF
        CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
        IF (failure) THEN
           EXIT
        END IF
     END DO
     WRITE(unit=unit_nr,fmt='()',iostat=iostat)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  END DO
END SUBROUTINE cp_3d_r_write
!***************************************************************************

!!****f* cp_array_r_utils/cp_1d_r_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_1d_r_guarantee_size
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
  SUBROUTINE cp_1d_r_guarantee_size(array, n, init, error)
    real(kind=dp), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(in) :: n
    real(kind=dp), INTENT(in), OPTIONAL :: init
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='cp_1d_r_guarantee_size',&
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
  END SUBROUTINE cp_1d_r_guarantee_size
!***************************************************************************

!!****f* cp_array_r_utils/cp_2d_r_guarantee_size [1.0] *
!!
!!   NAME
!!     cp_2d_r_guarantee_size
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
  SUBROUTINE cp_2d_r_guarantee_size(array, n_rows, n_cols ,&
       init, error)
    real(kind=dp), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(in) :: n_rows, n_cols
    real(kind=dp), INTENT(in), OPTIONAL :: init
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='cp_2d_r_guarantee_size',&
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
  END SUBROUTINE cp_2d_r_guarantee_size
!***************************************************************************

!!****f* cp_array_r_utils/cp_1d_r_bsearch [1.0] *
!!
!!   NAME
!!     cp_1d_r_bsearch
!!
!!   FUNCTION
!!     returns the index at which the element el should be inserted in the
!!     array to keep it ordered (array(i)>=el).
!!     If the element is bigger than all the elements in the array returns
!!     the last index+1.
!!
!!   NOTES
!!     the array should be ordered in growing order
!!
!!   INPUTS
!!     - array: the array to search
!!     - el: the element to look for
!!     - l_index: the lower index for binary search (defaults to 1)
!!     - u_index: the upper index for binary search (defaults to size(array))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION cp_1d_r_bsearch(array, el, l_index, u_index, error)&
     result(res)
  real(kind=dp), DIMENSION(:), intent(in) :: array
  real(kind=dp), intent(in) :: el
  INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
  type(cp_error_type), optional, intent(inout) :: error
  integer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_1d_r_bsearch',&
        routineP=moduleN//':'//routineN
  INTEGER :: lindex,uindex,aindex
  failure=.false.
  
  lindex=1
  uindex=size(array)
  if (present(l_index)) lindex=l_index
  if (present(u_index)) uindex=u_index
  DO WHILE (lindex <= uindex)
     aindex=(lindex+uindex)/2
     IF (CP_R_LESS_Q(array(aindex),el)) THEN
        lindex=aindex+1
     ELSE
        uindex=aindex-1
     END IF
  END DO
  res=lindex
END FUNCTION cp_1d_r_bsearch
!***************************************************************************

  ! template def put here so that line numbers in template and derived 
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,USE,defaultFormatType1,lessQ,defines)]
! ARGS:
!  USE = ""
!  defaultFormatType1 = ""(es14.6)""
!  defines = 
!    "#define CP_R_LESS_Q(el1,el2) ( el1 < el2 )
!     "
!  lessQ = "CP_R_LESS_Q"
!  nametype1 = "r"
!  type1 = "real(kind=dp)"


END MODULE cp_array_r_utils

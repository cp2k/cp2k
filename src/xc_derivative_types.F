!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/xc_derivative_types *
!!
!!   NAME
!!     xc_derivative_types
!!
!!   FUNCTION
!!     Provides types for the management of the xc-functionals and
!!     their derivatives.
!!
!!   NOTES
!!     
!!     
!!****
MODULE xc_derivative_types
  USE cp_array_r_utils,                ONLY: cp_3d_r_write
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_cr3d,&
                                             pw_pool_type
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             MAX_LABEL_LENGTH,&
                                             create_split_derivative_desc,&
                                             standardize_dervitative_desc
  IMPLICIT NONE

  PRIVATE

  INTEGER, SAVE :: derivative_id_nr = 1

  CHARACTER(len=*), PARAMETER :: moduleN='xc_derivative_types'

   PUBLIC :: xc_derivative_type, xc_derivative_p_type
   PUBLIC :: xc_derivative_create, xc_derivative_release, xc_derivative_retain,&
        xc_derivative_get, xc_derivative_write
   
!!*** *************************************************************************

!!****s* xc_derivative_types/xc_derivative_type  *
!!
!!   NAME
!!     xc_derivative_type
!!
!!   FUNCTION
!!     represent a derivative of a functional
!!
!!****
  TYPE xc_derivative_type
     INTEGER                                   :: ref_count, id_nr
     CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: desc
     CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), POINTER :: split_desc
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER    :: deriv_data
  END TYPE xc_derivative_type
!!*** **********************************************************************

!!****s* xc_derivative_types/xc_derivative_p_type [1.0] *
!!
!!   NAME
!!     xc_derivative_p_type
!!
!!   FUNCTION
!!     represent a pointer to a derivative (to have arrays of derivatives)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - deriv: the pointer to the derivative
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE xc_derivative_p_type
     TYPE(xc_derivative_type), POINTER :: deriv
  END TYPE xc_derivative_p_type
!!***
  !****************************************************************************


CONTAINS

!!****f* xc_derivative_types/xc_derivative_create [1.0] *
!!
!!   NAME
!!     xc_derivative_create
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_create(derivative, desc, cr3d_ptr, error)
!!       Type(xc_derivative_type), Pointer:: derivative
!!       Character(Len=*), Intent (IN):: desc
!!       Real(Kind=dp), Dimension(:,:,:), Pointer:: cr3d_ptr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_create
!!
!!   FUNCTION
!!     allocates and initializes a derivative type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - derivative  : the object to create
!!     - desc : the derivative description
!!     - cr3d_ptr    : the data array (the ownership of it passes to the 
!!                     derivative type), the array is not zeroed
!!     - error       : variable to control error logging, stopping,... 
!!                     see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_create(derivative, desc, cr3d_ptr, error)

    ! arguments
    TYPE(xc_derivative_type), POINTER              :: derivative
    CHARACTER(len=*), INTENT(in)                   :: desc
    REAL(kind=dp), DIMENSION(:,:,:), POINTER       :: cr3d_ptr
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error

    ! locals
    LOGICAL                                   :: failure
    INTEGER                                   :: stat
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_desc

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_create',&
         routineP = moduleN//':'//routineN

    failure=.FALSE.

    ALLOCATE(derivative, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       derivative%ref_count = 1
       derivative%id_nr     = derivative_id_nr
       derivative_id_nr     = derivative_id_nr + 1
       my_desc = standardize_dervitative_desc(desc)
       derivative%split_desc => &
            create_split_derivative_desc(my_desc)
       derivative%desc = my_desc
       derivative%deriv_data => cr3d_ptr
    END IF

  END SUBROUTINE xc_derivative_create
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_retain [1.0] *
!!
!!   NAME
!!     xc_derivative_retain
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_retain(deriv, error)
!!       Type(xc_derivative_type), Pointer:: deriv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_retain
!!
!!   FUNCTION
!!     retains the given derivative (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - deriv: the object you want to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xc_derivative_retain(deriv,error)
  TYPE(xc_derivative_type), POINTER :: deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xc_derivative_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(deriv%ref_count>0,cp_failure_level,routineP,error)
     deriv%ref_count=deriv%ref_count+1
  END IF
END SUBROUTINE xc_derivative_retain
!***************************************************************************

!!****f* xc_derivative_types/xc_derivative_release [1.0] *
!!
!!   NAME
!!     xc_derivative_release
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_release(derivative, pw_pool, error)
!!       Type(xc_derivative_type), Pointer:: derivative
!!       Type(pw_pool_type), Pointer, Optional:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_release
!!
!!   FUNCTION
!!     allocates and initializes a derivative type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - derivative  : the object to create
!!     - pw_pool     : if given gives back the cr3d array %deriv_data back to it
!!                     instead of deallocating it
!!     - error       : variable to control error logging, stopping,... 
!!                     see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_release(derivative, pw_pool, error)

    ! arguments
    TYPE(xc_derivative_type), POINTER              :: derivative
    TYPE(pw_pool_type), POINTER, OPTIONAL          :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error

    ! locals
    LOGICAL                                   :: failure
    INTEGER                                   :: stat
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_desc

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_release',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(derivative),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative%ref_count>=1,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       derivative%ref_count = derivative%ref_count - 1
       IF (derivative%ref_count == 0) THEN
          IF (PRESENT(pw_pool)) THEN
             IF (ASSOCIATED(pw_pool)) THEN
                CALL pw_pool_give_back_cr3d(pw_pool, derivative%deriv_data,&
                     accept_non_compatible=.TRUE.,error=error)
             END IF
          END IF
          IF (ASSOCIATED(derivative%deriv_data)) THEN
             DEALLOCATE(derivative%deriv_data, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(derivative%split_desc, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

          DEALLOCATE(derivative, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(derivative)

  END SUBROUTINE xc_derivative_release
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_get [1.0] *
!!
!!   NAME
!!     xc_derivative_get
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_get(deriv, desc, split_desc, order,&
!!         deriv_data, error)
!!       Type(xc_derivative_type), Pointer:: deriv
!!       Character(Len=max_derivative_desc_length), Pointer, Optional::&
!!         desc
!!       Character(Len=max_label_length), Dimension(:), Pointer, Optional::&
!!         split_desc
!!       Integer, Intent (OUT), Optional:: order
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: deriv_data
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_get
!!
!!   FUNCTION
!!     returns various information on the given derivative
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - deriv: thederivative you whant information about
!!     - desc: a string that describes the derivative (empty string 
!!       for the function itself, otherwise a string of the form 
!!       "(rho)(rho)(norm_drhoa)")
!!     - split_desc: the same as desc but with an array 
!!       of stings, and a derivative coordinate in each string (the previous
!!       example would become (/"rho","rho","norm_drhoa"/)
!!     - order: the order of the derivative
!!     - deriv_data: the 3d real array with the derivative
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xc_derivative_get(deriv,desc,split_desc,&
     order,deriv_data,error)
  TYPE(xc_derivative_type), POINTER :: deriv
  CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH), POINTER, OPTIONAL :: desc
  CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), POINTER, OPTIONAL :: &
       split_desc
  INTEGER, INTENT(out), OPTIONAL :: order
  REAL(kind=dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: deriv_data
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xc_derivative_get',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(desc)) desc => deriv%desc
     IF (PRESENT(split_desc)) split_desc => deriv%split_desc
     IF (PRESENT(deriv_data)) deriv_data => deriv%deriv_data
     IF (PRESENT(order)) order=SIZE(deriv%split_desc)
  END IF
END SUBROUTINE xc_derivative_get
!***************************************************************************

!!****f* xc_derivative_types/xc_derivative_write [1.0] *
!!
!!   NAME
!!     xc_derivative_write
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_write(deriv, unit_nr, long_description,&
!!         error)
!!       Type(xc_derivative_type), Pointer:: deriv
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: long_description
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_write
!!
!!   FUNCTION
!!     writes information about the current derivative to the given unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - deriv: the derivative to write out
!!     - unit_nr: the unit where to write
!!     - long_description: if a long description (with the derivative data)
!!       should be written out. Defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xc_derivative_write(deriv,unit_nr,long_description,error)
  TYPE(xc_derivative_type), POINTER :: deriv
  INTEGER, INTENT(in) :: unit_nr
  LOGICAL, INTENT(in), OPTIONAL :: long_description
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, my_long_description
  INTEGER :: i
  CHARACTER(len=*), PARAMETER :: routineN='xc_derivative_write',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  my_long_description=.FALSE.
  IF (PRESENT(long_description)) my_long_description=long_description
  
  IF (ASSOCIATED(deriv)) THEN
     WRITE(unit=unit_nr,fmt="(' <xc_derivative_type>{ id_nr=',i6,' ref_count=',i6,',')")&
          deriv%id_nr,deriv%ref_count
     WRITE(unit=unit_nr,fmt="(' desc=',a,',')") TRIM(deriv%desc)
     IF (ASSOCIATED(deriv%split_desc)) THEN
        WRITE(unit=unit_nr,fmt="(' split_desc=(')", advance="NO")
        DO i=1,SIZE(deriv%split_desc)-1
           WRITE(unit=unit_nr,fmt="(a,',')", advance="NO") &
                TRIM(deriv%split_desc(i))
        END DO
        IF (SIZE(deriv%split_desc)>0) THEN
           WRITE(unit=unit_nr,fmt="(a,')')")&
                TRIM(deriv%split_desc(SIZE(deriv%split_desc)))
        ELSE
           WRITE(unit=unit_nr,fmt="(')')")
        END IF
     ELSE
        WRITE(unit=unit_nr,fmt="(' split_desc=*null*')")
     END IF
     IF (ASSOCIATED(deriv%deriv_data)) THEN
        IF (my_long_description) THEN
           WRITE(unit=unit_nr,fmt="(' deriv_data=')")
           CALL cp_3d_r_write(deriv%deriv_data,unit_nr=unit_nr,error=error)
        ELSE
           WRITE(unit=unit_nr,fmt="(' deriv_data=<array(',i6,':',i6,',',i6,':',i6,',',i6,':',i6,')>')")&
                LBOUND(deriv%deriv_data,1),UBOUND(deriv%deriv_data,1),&
                LBOUND(deriv%deriv_data,2),UBOUND(deriv%deriv_data,2),&
                LBOUND(deriv%deriv_data,3),UBOUND(deriv%deriv_data,3)
        END IF
     ELSE
        WRITE(unit=unit_nr,fmt="(' deriv_data=*null*')")
     END IF
     WRITE(unit=unit_nr,fmt="(' }')")
  ELSE
     WRITE(unit=unit_nr,fmt="(' <xc_derivative_type*null*>')")
  END IF
END SUBROUTINE xc_derivative_write
!***************************************************************************

END MODULE xc_derivative_types


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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             standardize_dervitative_desc
  IMPLICIT NONE

  PRIVATE

  INTEGER, SAVE :: derivative_id_nr = 1
  INTEGER, SAVE :: derivative_set_id_nr  = 1

  CHARACTER(len=*), PARAMETER :: moduleN='xc_derivative_types'

!!****t* cp2k/xc_derivative_type  *
!!
!!   NAME
!!     xc_derivative_type
!!
!!   FUNCTION
!!     A node of a linked list contains the actual derivative together
!!     with a string describing it.
!!
!!****
  TYPE xc_derivative_type
     INTEGER                                   :: ref_count, id_nr
     CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: description
     REAL(KIND = dp), DIMENSION(:), POINTER    :: DATA
     TYPE(xc_derivative_type), POINTER         :: next
  END TYPE xc_derivative_type

!!****t* cp2k/xc_derivative_set_type  *
!!
!!   NAME
!!     xc_derivative_set_type
!!
!!   FUNCTION
!!     A derivative set contains the different derivatives of a xc-functional
!!     in form of a linked list
!!
!!   NOTE
!!     To all USERS of this type: never reference data directly. data is an
!!     internal variable, which due to fortrans limitations cannot be made
!!     private
!!
!!****
  TYPE xc_derivative_set_type
     INTEGER                           :: ref_count, id_nr, npoints
     TYPE(xc_derivative_type), POINTER :: first
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: DATA
  END TYPE xc_derivative_set_type

   PUBLIC :: xc_derivative_type, &
            xc_derivative_set_type, &
            xc_derivative_set_create, &
            xc_derivative_set_release, &
            xc_dset_get_derivative

CONTAINS

  ! *****************************************************************************

!!****f* xc_derivative_types/xc_dset_get_derivative [1.0] *
!!
!!   NAME
!!     xc_dset_get_derivative
!!
!!   SYNOPSIS
!!     Function xc_dset_get_derivative(derivative_set, description,&
!!         allocate_deriv, error) Result(res)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Character(Len=*), Intent (IN):: description
!!       Logical, Intent (IN), Optional:: allocate_deriv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(xc_derivative_type), Pointer:: res
!!     End Function xc_dset_get_derivative
!!
!!   FUNCTION
!!     returns the requested xc_derivative
!!
!!   NOTES
!!     -
!!     
!!   INPUTS
!!     - derivative_set: the set where to search for the derivative
!!     - description   : the description of the derivative you want to have
!!     - allocate_deriv: if the derivative should be allocated when not present
!!                       Defaults to false.
!!     - error         : variable to control error logging, stopping,... 
!!                       see module cp_error_handling 
!!
!!*** **********************************************************************
  FUNCTION xc_dset_get_derivative(derivative_set, description, allocate_deriv, error) &
       RESULT(res)

    ! arguments
    TYPE(xc_derivative_set_type), POINTER        :: derivative_set
    CHARACTER(len=*), INTENT(in)                 :: description
    LOGICAL, INTENT(in), OPTIONAL                :: allocate_deriv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    TYPE(xc_derivative_type), POINTER            :: res

    ! locals
    LOGICAL                                   :: failure, my_allocate_deriv
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: std_deriv_desc
    TYPE(xc_derivative_type), POINTER         :: pos

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_get_derivative',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)
    my_allocate_deriv=.FALSE.
    IF (PRESENT(allocate_deriv)) my_allocate_deriv=allocate_deriv
    NULLIFY(res)
    std_deriv_desc = standardize_dervitative_desc(description)
    IF (.NOT. failure) THEN
       pos => derivative_set%first
       DO WHILE(ASSOCIATED(pos))
          IF (pos%description == std_deriv_desc) THEN
             res => pos
             EXIT
          END IF
          pos => pos%next
       END DO
       IF (.not.ASSOCIATED(res).and.my_allocate_deriv) THEN
          CALL xc_derivative_create(res, std_deriv_desc, & 
                                    derivative_set%npoints, error=error)
          res%next => derivative_set%first
          derivative_set%first => res
       END IF
    END IF
  END FUNCTION xc_dset_get_derivative
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_set_create [1.0] *
!!
!!   NAME
!!     xc_derivative_set_create
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_set_create(derivative_set, npoints, error)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Integer, Intent (IN):: npoints
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_set_create
!!
!!   FUNCTION
!!     creates a derivative set object
!!
!!   INPUTS
!!     - derivative_set : the set where to search for the derivative
!!     - npoints : the number of points in the derivatives
!!     - error   : variable to control error logging, stopping,... 
!!                 see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_set_create(derivative_set, npoints, error)

    ! arguments
    TYPE(xc_derivative_set_type), POINTER        :: derivative_set
    INTEGER, INTENT(in)                          :: npoints
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL :: failure
    INTEGER :: stat

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_set_create',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.
    
    CPPrecondition(.not.ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(npoints>=0,cp_failure_level,routineP,error,failure)

    ALLOCATE(derivative_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       NULLIFY(derivative_set%first)
       NULLIFY(derivative_set%data)
       derivative_set%ref_count  = 1
       derivative_set%id_nr      = derivative_set_id_nr
       derivative_set_id_nr      = derivative_set_id_nr + 1
       derivative_set%npoints    = npoints
    END IF

  END SUBROUTINE xc_derivative_set_create
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_set_release [1.0] *
!!
!!   NAME
!!     xc_derivative_set_release
!!
!!   SYNOPSIS
!!     Subroutine xc_derivative_set_release(derivative_set, error)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_derivative_set_release
!!
!!   FUNCTION
!!     releases a derivative set
!!
!!   INPUTS
!!     - set     : the set to release
!!     - error   : variable to control error logging, stopping,... 
!!                 see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_set_release(derivative_set, error)

    ! arguments
    TYPE(xc_derivative_set_type), POINTER        :: derivative_set
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL :: failure
    INTEGER :: stat
    TYPE(xc_derivative_type), POINTER :: pos, tmp

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_set_release',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       derivative_set%ref_count = derivative_set%ref_count - 1
       IF (derivative_set%ref_count == 0) THEN
          IF (ASSOCIATED(derivative_set%data)) THEN
             DEALLOCATE(derivative_set%data, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          pos => derivative_set%first
          DO WHILE (ASSOCIATED(pos))
             tmp => pos
             pos => pos%next
             CALL xc_derivative_release(tmp, error=error)
          END DO
          DEALLOCATE(derivative_set, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(derivative_set)

  END SUBROUTINE xc_derivative_set_release
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_create [1.0] *
!!
!!   NAME
!!     xc_derivative_create
!!
!!   FUNCTION
!!     allocates and initializes a derivative type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - derivative  : the object to create
!!     - description : the derivative description
!!     - npoints     : size of the data array
!!     - error       : variable to control error logging, stopping,... 
!!                     see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_create(derivative, description, npoints, error)

    ! arguments
    TYPE(xc_derivative_type), POINTER              :: derivative
    CHARACTER(len=*), INTENT(in)                   :: description
    INTEGER, INTENT(in)                            :: npoints
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error

    ! locals
    LOGICAL                                   :: failure
    INTEGER                                   :: stat
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_description

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_create',&
         routineP = moduleN//':'//routineN

    failure=.FALSE.

    ALLOCATE(derivative, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       NULLIFY(derivative%data)
       NULLIFY(derivative%next)
       derivative%ref_count = 1
       derivative%id_nr     = derivative_id_nr
       derivative_id_nr     = derivative_id_nr + 1
       my_description = standardize_dervitative_desc(description)
       derivative%description = my_description
       ALLOCATE(derivative%data(npoints), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE xc_derivative_create
  !***************************************************************************

!!****f* xc_derivative_types/xc_derivative_release [1.0] *
!!
!!   NAME
!!     xc_derivative_release
!!
!!   FUNCTION
!!     allocates and initializes a derivative type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - derivative  : the object to create
!!     - error       : variable to control error logging, stopping,... 
!!                     see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_derivative_release(derivative, error)

    ! arguments
    TYPE(xc_derivative_type), POINTER              :: derivative
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error

    ! locals
    LOGICAL                                   :: failure
    INTEGER                                   :: stat
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_description

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_release',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(derivative),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative%ref_count>=1,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       derivative%ref_count = derivative%ref_count - 1
       IF (derivative%ref_count == 0) THEN
          DEALLOCATE(derivative%data, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(derivative, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(derivative)

  END SUBROUTINE xc_derivative_release
  !***************************************************************************

END MODULE xc_derivative_types


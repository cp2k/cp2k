!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_derivative_desc [1.0] *
!!
!!   NAME
!!     xc_derivative_desc
!!
!!   FUNCTION
!!     Module with functions to handle derivative descriptors.
!!     derivative description are strings have the following form
!!     "(rho)(rho)(rhoa)" which means that it is a forth order
!!     derivative, twice with respect to rho, once with respect to rhoa
!!     and once with respect to drhoa.
!!     Possible derivatives are:
!!     - rho: total density
!!     - norm_drho: norm of the gradient of the total density
!!     - rhoa, rhob: alpha and beta spin density (with LSD)
!!     - norm_drhoa, norm_drhob: norm of the gradient of the alpha and beta
!!       spin density
!!     - drhoa_drhob: scalar product between the gradient of rhoa and the
!!       grandient of rhob
!!
!!   AUTHOR
!!     thomas & fawzi
!!
!!   NOTES
!!      add kinetic, drhox, drhoy, drhoz, drhoax,...?
!!
!!****
MODULE xc_derivative_desc
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: &
       MAX_LABEL_LENGTH           = 10, &
       MAX_DERIVATIVE             = 4,  &
       MAX_DERIVATIVE_DESC_LENGTH = (MAX_LABEL_LENGTH+2)*MAX_DERIVATIVE

  LOGICAL, PARAMETER :: debug_this_module=.FALSE.

  CHARACTER(len=*), PARAMETER :: moduleN='xc_derivative_desc'

  PUBLIC :: MAX_DERIVATIVE_DESC_LENGTH, MAX_LABEL_LENGTH,&
            standardize_dervitative_desc, &
            create_split_derivative_desc

!****************************************************************************

CONTAINS

!!****f* xc_derivative_types/standardize_dervitative_desc *
!!
!!   NAME
!!     standardize_dervitative_desc
!!
!!   SYNOPSIS
!!     Function standardize_dervitative_desc(deriv_desc, error) Result(res)
!!       Character(Len=*), Intent (IN):: deriv_desc
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!       Character(Len=max_derivative_desc_length):: res
!!     End Function standardize_dervitative_desc
!!
!!   FUNCTION
!!     returns a unique string for a given derivative string
!!
!!   ARGUMENTS
!!     - deriv_desc: input string which describes the derivative and
!!                   should be standardized
!!     - error     : variable to control error logging, stopping,... 
!!                   see module cp_error_handling 
!!
!!****
 FUNCTION standardize_dervitative_desc(deriv_desc, error) RESULT(res)
    CHARACTER(len=*), INTENT(in) :: deriv_desc
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: res

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='standardize_dervitative_desc',&
         routineP=moduleN//':'//routineN
    CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), POINTER :: deriv_array

    INTEGER :: i, pos, l_label, stat
    LOGICAL :: ordered
    CHARACTER(len=MAX_LABEL_LENGTH) :: tmp

    deriv_array => create_split_derivative_desc(deriv_desc, error=error)

    ordered=.FALSE.
    DO WHILE(.not.ordered)
       ordered=.TRUE.
       DO i=1,SIZE(deriv_array)-1
          IF (deriv_array(i)>deriv_array(i+1)) THEN
             tmp=deriv_array(i+1)
             deriv_array(i+1)=deriv_array(i)
             deriv_array(i)=tmp
             ordered=.FALSE.
          END IF
       END DO
    END DO
    
    res=""
    pos=1
    DO i=1,SIZE(deriv_array)
       l_label=LEN_TRIM(deriv_array(i))
       res(pos:pos+l_label+1)='('//deriv_array(i)(1:l_label)//')'
       pos=pos+l_label+2
    END DO

    DEALLOCATE(deriv_array, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  END FUNCTION standardize_dervitative_desc

!!****f* xc_derivative_types/create_split_derivative_desc [1.0] *
!!
!!   NAME
!!     create_split_derivative_desc
!!
!!   SYNOPSIS
!!     Function create_split_derivative_desc(deriv_desc, error) Result(res)
!!       Implicit None
!!       Character(Len=*), Intent (IN):: deriv_desc
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!       Character(Len=max_label_length), Dimension(:), Pointer:: res
!!     End Function create_split_derivative_desc
!!
!!   FUNCTION
!!     returns an array of 1st order derivative descriptions
!!
!!   ARGUMENTS
!!     - deriv_desc: input string which describes the derivative
!!     - error     : variable to control error logging, stopping,... 
!!                   see module cp_error_handling 
!!
!!*** *********************************************************************
  FUNCTION create_split_derivative_desc(deriv_desc, error) RESULT(res)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: deriv_desc
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), POINTER :: res

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_split_derivative_desc',&
         routineP=moduleN//':'//routineN
    INTEGER :: nderiv, i, j, stat

    nderiv=0
    DO i=1,LEN(deriv_desc)
       IF (deriv_desc(i:i)=='(') THEN
          nderiv=nderiv+1
       END IF
    END DO

    ALLOCATE(res(nderiv),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nderiv=0
    DO i=1,LEN(deriv_desc)
       IF (deriv_desc(i:i)=='(') THEN
          nderiv=nderiv+1
          DO j=i+1,LEN(deriv_desc)
             IF (deriv_desc(j:j)==')') EXIT
          END DO
!tc: should we do a check on the derivative descriptions?
          res(nderiv)=deriv_desc(i+1:j-1)
       END IF
    END DO
    
  END FUNCTION create_split_derivative_desc

! *****************************************************************************

END MODULE xc_derivative_desc

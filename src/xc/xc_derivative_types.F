!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides types for the management of the xc-functionals and
!>      their derivatives.
! *****************************************************************************
MODULE xc_derivative_types
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_cr3d,&
                                             pw_pool_type
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             MAX_LABEL_LENGTH,&
                                             create_split_derivative_desc,&
                                             standardize_derivative_desc
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTEGER, SAVE :: derivative_id_nr = 1

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_derivative_types'

  PUBLIC :: xc_derivative_type, xc_derivative_p_type
  PUBLIC :: xc_derivative_create, xc_derivative_release,&
            xc_derivative_get

! *****************************************************************************
!> \brief represent a derivative of a functional
! *****************************************************************************
  TYPE xc_derivative_type
     INTEGER                                   :: ref_count, id_nr
     CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH) :: desc
     CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), POINTER :: split_desc
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER    :: deriv_data
  END TYPE xc_derivative_type

! *****************************************************************************
!> \brief represent a pointer to a derivative (to have arrays of derivatives)
!> \param deriv the pointer to the derivative
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE xc_derivative_p_type
     TYPE(xc_derivative_type), POINTER :: deriv
  END TYPE xc_derivative_p_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes a derivative type
!> \param derivative the object to create
!> \param desc the derivative description
!> \param cr3d_ptr ...
!> \param error variable to control error logging, stopping,...
!>                      see module cp_error_handling
!> \param cr 3d_ptr    : the data array (the ownership of it passes to the
!>                      derivative type), the array is not zeroed UNMATCHED_PROCEDURE_ARGUMENT: please check 
! *****************************************************************************
  SUBROUTINE xc_derivative_create(derivative, desc, cr3d_ptr, error)

    TYPE(xc_derivative_type), POINTER        :: derivative
    CHARACTER(len=*), INTENT(in)             :: desc
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_create', &
      routineP = moduleN//':'//routineN

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: my_desc
    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(derivative, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       derivative%ref_count = 1
       derivative%id_nr     = derivative_id_nr
       derivative_id_nr     = derivative_id_nr + 1
       CALL standardize_derivative_desc(desc,my_desc,error=error)
       CALL create_split_derivative_desc(my_desc,derivative%split_desc,error=error)
       derivative%desc = my_desc
       derivative%deriv_data => cr3d_ptr
    END IF

  END SUBROUTINE xc_derivative_create

! *****************************************************************************
!> \brief retains the given derivative (see doc/ReferenceCounting.html)
!> \param deriv the object you want to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE xc_derivative_retain(deriv,error)
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(deriv%ref_count>0,cp_failure_level,routineP,error)
     deriv%ref_count=deriv%ref_count+1
  END IF
END SUBROUTINE xc_derivative_retain

! *****************************************************************************
!> \brief allocates and initializes a derivative type
!> \param derivative the object to create
!> \param pw_pool if given gives back the cr3d array %deriv_data back to it
!>                      instead of deallocating it
!> \param error variable to control error logging, stopping,...
!>                      see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_derivative_release(derivative, pw_pool, error)

    TYPE(xc_derivative_type), POINTER        :: derivative
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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

! *****************************************************************************
!> \brief returns various information on the given derivative
!> \param deriv thederivative you whant information about
!> \param desc a string that describes the derivative (empty string
!>        for the function itself, otherwise a string of the form
!>        "(rho)(rho)(norm_drhoa)")
!> \param split_desc the same as desc but with an array of strings,
!>        and a derivative coordinate in each string (the previous
!>        example would become (/"rho","rho","norm_drhoa"/)
!> \param order the order of the derivative
!> \param deriv_data the 3d real array with the derivative
!> \param accept_null_data if deriv_data can be unassociated (defaults to no)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
SUBROUTINE xc_derivative_get(deriv,desc,split_desc,&
     order,deriv_data,accept_null_data, error)
    TYPE(xc_derivative_type), POINTER        :: deriv
    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH), &
      INTENT(out), OPTIONAL                  :: desc
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), OPTIONAL, POINTER        :: split_desc
    INTEGER, INTENT(out), OPTIONAL           :: order
    REAL(kind=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: deriv_data
    LOGICAL, INTENT(in), OPTIONAL            :: accept_null_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_derivative_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_accept_null_data

  failure=.FALSE.
  my_accept_null_data=.FALSE.
  IF (PRESENT(accept_null_data)) my_accept_null_data=accept_null_data

  CPPrecondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(desc)) desc = deriv%desc
     IF (PRESENT(split_desc)) split_desc => deriv%split_desc
     IF (PRESENT(deriv_data)) THEN
        deriv_data => deriv%deriv_data
        IF (.NOT.my_accept_null_data) THEN
           CPPrecondition(ASSOCIATED(deriv_data),cp_failure_level,routineP,error,failure)
        END IF
     END IF
     IF (PRESENT(order)) order=SIZE(deriv%split_desc)
  END IF
END SUBROUTINE xc_derivative_get

END MODULE xc_derivative_types


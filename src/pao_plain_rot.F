!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_plain_rot [1.0] *
!!
!!   NAME
!!     pao_plain_rot
!!
!!   FUNCTION
!!     various function to perform plain rotations, like jacobi,...
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pao_plain_rot
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sup_norm_matrix
  USE kinds,                           ONLY: dp
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='pao_plain_rot'

  PUBLIC :: cp_left_compose_p_rot, cp_right_compose_p_rot,&
       cp_right_compose_p_rot_deriv, cp_left_compose_p_rot_deriv

!!***
!****************************************************************************

CONTAINS


!!****f* pao_plain_rot/cp_right_compose_p_rot [1.0] *
!!
!!   NAME
!!     cp_right_compose_p_rot
!!
!!   SYNOPSIS
!!     Subroutine cp_right_compose_p_rot(matrix, angle, i, j, error)
!!       Real(KIND = dp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(KIND = dp), Intent (IN):: angle
!!       Integer, Intent (IN):: i, j
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_right_compose_p_rot
!!
!!   FUNCTION
!!     calculates matrix*J where J is the plain rotation of angle angle
!!     between i1 and j1.
!!v      i1    j1
!!v i1 (cos , -sin)
!!v j1 (sin ,  cos)
!! could be more efficent.
!!     If the angle is small (in -pi/2..pi/2, or at least -pi..pi), as it 
!!     should be if the basis function are ordered so that the minimal basis
!!     goes in the small basis with just the projection (no rotation), then
!!     there are much more efficient ways to compute this using tan in place 
!!     of the angle (see Parlett, The symmetric eigenvalue problem, 6.8).
!!     Hopefully it gets inlined.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that gets modified
!!     angle: the angle of the rotation
!!     i: the first index defining the rotation plane
!!     j: the second index defining the rotation plane
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_right_compose_p_rot(matrix,angle,i,j,error)
    REAL(KIND = dp), DIMENSION(:,:), INTENT(inout) :: matrix 
    REAL(KIND = dp), INTENT(in) :: angle
    INTEGER, INTENT(in) :: i,j
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: failure
    REAL(KIND = dp) :: rcos,rsin
    REAL(KIND = dp), DIMENSION(SIZE(matrix,1)) :: tmpV1, tmpV2
    CHARACTER(len=*), PARAMETER :: routineN='cp_right_compose_p_rot',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.
    IF (cp_debug .AND. debug_this_module) THEN ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       rcos=COS(angle)
       rsin=SIN(angle)
       tmpV1 = rcos*matrix(:,i)+rsin*matrix(:,j)
       tmpV2 = -rsin*matrix(:,i)+rcos*matrix(:,j)
       matrix(:,i)=tmpV1
       matrix(:,j)=tmpV2
    END IF
  END SUBROUTINE cp_right_compose_p_rot
!***************************************************************************


!!****f* pao_plain_rot/cp_left_compose_p_rot [1.0] *
!!
!!   NAME
!!     cp_left_compose_p_rot
!!
!!   SYNOPSIS
!!     Subroutine cp_left_compose_p_rot(matrix, angle, i, j, error)
!!       Real(KIND = dp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(KIND = dp), Intent (IN):: angle
!!       Integer, Intent (IN):: i, j
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_left_compose_p_rot
!!
!!   FUNCTION
!!     calculates J*matrix where J is the plain rotation of angle angle
!!     between i1 and j1 (less efficent than matrix*J because the slice
!!     in continuos in memory).
!!v      i1    j1
!!v i1 (cos , -sin)
!!v j1 (sin ,  cos)
!!     Quite inefficent.
!!     If the angle is small (in -pi/2..pi/2, or at least -pi..pi), as it 
!!     should be if the basis function are ordered so that the minimal basis
!!     goes in the small basis with just the projection (no rotation), then
!!     there are much more efficient ways to compute this using tan in place 
!!     of the angle (see Parlett, The symmetric eigenvalue problem, 6.8).
!!     Hopefully it gets inlined.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that gets modified
!!     angle: the angle of the rotation
!!     i: the first index defining the rotation plane
!!     j: the second index defining the rotation plane
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_left_compose_p_rot(matrix,angle,i,j, error)
    REAL(KIND = dp), DIMENSION(:,:), INTENT(inout) :: matrix 
    REAL(KIND = dp), INTENT(in) :: angle
    INTEGER, INTENT(in) :: i,j
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: failure
    REAL(KIND = dp) :: rcos,rsin
    REAL(KIND = dp), DIMENSION(SIZE(matrix,2)) :: tmpV1,tmpV2
    CHARACTER(len=*), PARAMETER :: routineN='cp_left_compose_p_rot',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.
    IF (debug_this_module) THEN ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       rcos=COS(angle)
       rsin=SIN(angle)
       tmpV1 = rcos*matrix(i,:)-rsin*matrix(j,:)
       tmpV2 = rsin*matrix(i,:)+rcos*matrix(j,:)
       matrix(i,:)=tmpV1
       matrix(j,:)=tmpV2
    END IF
  END SUBROUTINE cp_left_compose_p_rot
!***************************************************************************

!!****f* pao_plain_rot/cp_right_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     cp_right_compose_p_rot_deriv
!!
!!   SYNOPSIS
!!     Subroutine cp_right_compose_p_rot_deriv(matrix, angle, i, j, error)
!!       Real(KIND = dp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(KIND = dp), Intent (IN):: angle
!!       Integer, Intent (IN):: i, j
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_right_compose_p_rot_deriv
!!
!!   FUNCTION
!!     calculates matrix*dJ where dJ is the derivative of the plain 
!!     rotation of angle angle between i1 and j1.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that gets modified
!!     angle: the angle of the rotation
!!     i: the first index defining the rotation plane
!!     j: the second index defining the rotation plane
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_right_compose_p_rot_deriv(matrix,angle,i,j,error)
    REAL(KIND = dp), DIMENSION(:,:), INTENT(inout) :: matrix 
    REAL(KIND = dp), INTENT(in) :: angle
    INTEGER, INTENT(in) :: i,j
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: failure
    REAL(KIND = dp) :: rcos,rsin
    REAL(KIND = dp), DIMENSION(SIZE(matrix,1)) :: tmpV1, tmpV2
    CHARACTER(len=*), PARAMETER :: routineN='cp_right_compose_p_rot_deriv',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.
    IF (cp_debug .AND. debug_this_module) THEN ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       rcos=COS(angle)
       rsin=SIN(angle)
       tmpV1 = -rsin*matrix(:,i)+rcos*matrix(:,j)
       tmpV2 = -rcos*matrix(:,i)-rsin*matrix(:,j)
       matrix=0.0_dp
       matrix(:,i)=tmpV1
       matrix(:,j)=tmpV2
    END IF
  END SUBROUTINE cp_right_compose_p_rot_deriv
!***************************************************************************

!!****f* pao_plain_rot/cp_left_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     cp_left_compose_p_rot_deriv
!!
!!   SYNOPSIS
!!     Subroutine cp_left_compose_p_rot_deriv(matrix, angle, i, j, error)
!!       Real(KIND = dp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(KIND = dp), Intent (IN):: angle
!!       Integer, Intent (IN):: i, j
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_left_compose_p_rot_deriv
!!
!!   FUNCTION
!!     calculates dJ*matrix where dJ is the derivative of the plain 
!!     rotation of angle angle between i1 and j1.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that gets modified
!!     angle: the angle of the rotation
!!     i: the first index defining the rotation plane
!!     j: the second index defining the rotation plane
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_left_compose_p_rot_deriv(matrix,angle,i,j,error)
    REAL(KIND = dp), DIMENSION(:,:), INTENT(inout) :: matrix 
    REAL(KIND = dp), INTENT(in) :: angle
    INTEGER, INTENT(in) :: i,j
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: failure
    REAL(KIND = dp) :: rcos,rsin
    REAL(KIND = dp), DIMENSION(SIZE(matrix,2)) :: tmpV1, tmpV2
    CHARACTER(len=*), PARAMETER :: routineN='cp_left_compose_p_rot_deriv',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.
    IF (cp_debug .AND. debug_this_module) THEN ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       rcos=COS(angle)
       rsin=SIN(angle)
       tmpV1 = -rsin*matrix(i,:)-rcos*matrix(j,:)
       tmpV2 = rcos*matrix(i,:)-rsin*matrix(j,:)
       matrix=0.0_dp
       matrix(i,:)=tmpV1
       matrix(j,:)=tmpV2
    END IF
  END SUBROUTINE cp_left_compose_p_rot_deriv
!***************************************************************************

END MODULE pao_plain_rot




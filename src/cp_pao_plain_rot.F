!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_plain_rot [1.0] *
!!
!!   NAME
!!     cp_pao_plain_rot
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

module cp_pao_plain_rot
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use cp_matrix_utils, only: cp_sup_norm_matrix
  implicit none
  private

  logical, parameter, private :: debug_this_module=.true.
  character(len=*), parameter, private :: moduleN='cp_pao_plain_rot'

  public :: cp_left_compose_p_rot, cp_right_compose_p_rot,&
       cp_right_compose_p_rot_deriv, cp_left_compose_p_rot_deriv

!!***
!****************************************************************************

contains


!!****f* cp_pao_plain_rot/cp_right_compose_p_rot [1.0] *
!!
!!   NAME
!!     cp_right_compose_p_rot
!!
!!   SYNOPSIS
!!     Subroutine cp_right_compose_p_rot(matrix, angle, i, j, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(Kind=wp), Intent (IN):: angle
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
  subroutine cp_right_compose_p_rot(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(cp_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,1)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='cp_right_compose_p_rot',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (cp_debug .and. debug_this_module) then ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = rcos*matrix(:,i)+rsin*matrix(:,j)
       tmpV2 = -rsin*matrix(:,i)+rcos*matrix(:,j)
       matrix(:,i)=tmpV1
       matrix(:,j)=tmpV2
    end if
  end subroutine cp_right_compose_p_rot
!***************************************************************************


!!****f* cp_pao_plain_rot/cp_left_compose_p_rot [1.0] *
!!
!!   NAME
!!     cp_left_compose_p_rot
!!
!!   SYNOPSIS
!!     Subroutine cp_left_compose_p_rot(matrix, angle, i, j, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(Kind=wp), Intent (IN):: angle
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
  subroutine cp_left_compose_p_rot(matrix,angle,i,j, error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(cp_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,2)) :: tmpV1,tmpV2
    character(len=*), parameter :: routineN='cp_left_compose_p_rot',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (debug_this_module) then ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = rcos*matrix(i,:)-rsin*matrix(j,:)
       tmpV2 = rsin*matrix(i,:)+rcos*matrix(j,:)
       matrix(i,:)=tmpV1
       matrix(j,:)=tmpV2
    end if
  end subroutine cp_left_compose_p_rot
!***************************************************************************

!!****f* cp_pao_plain_rot/cp_right_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     cp_right_compose_p_rot_deriv
!!
!!   SYNOPSIS
!!     Subroutine cp_right_compose_p_rot_deriv(matrix, angle, i, j, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(Kind=wp), Intent (IN):: angle
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
  subroutine cp_right_compose_p_rot_deriv(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(cp_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,1)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='cp_right_compose_p_rot_deriv',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (cp_debug .and. debug_this_module) then ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = -rsin*matrix(:,i)+rcos*matrix(:,j)
       tmpV2 = -rcos*matrix(:,i)-rsin*matrix(:,j)
       matrix=0.0_wp
       matrix(:,i)=tmpV1
       matrix(:,j)=tmpV2
    end if
  end subroutine cp_right_compose_p_rot_deriv
!***************************************************************************

!!****f* cp_pao_plain_rot/cp_left_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     cp_left_compose_p_rot_deriv
!!
!!   SYNOPSIS
!!     Subroutine cp_left_compose_p_rot_deriv(matrix, angle, i, j, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Real(Kind=wp), Intent (IN):: angle
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
  subroutine cp_left_compose_p_rot_deriv(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(cp_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,2)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='cp_left_compose_p_rot_deriv',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (cp_debug .and. debug_this_module) then ! this function must be fast
       CPPrecondition(i /= j,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = -rsin*matrix(i,:)-rcos*matrix(j,:)
       tmpV2 = rcos*matrix(i,:)-rsin*matrix(j,:)
       matrix=0.0_wp
       matrix(i,:)=tmpV1
       matrix(j,:)=tmpV2
    end if
  end subroutine cp_left_compose_p_rot_deriv
!***************************************************************************

end module cp_pao_plain_rot




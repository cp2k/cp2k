!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_pao_plain_rot [1.0] *
!!
!!   NAME
!!     qs_pao_plain_rot
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

module qs_pao_plain_rot
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp=>dp
  use qs_pao_utils, only: qs_sup_norm_matrix
  implicit none

  logical, parameter, private :: debug_this_module=.true.

!!***
!****************************************************************************

contains


!! calculates matrix*J where J is the plain rotation of angle angle
!! between i1 and j1.
!!v      i1    j1
!!v i1 (cos , -sin)
!!v j1 (sin ,  cos)
!! could be more efficent.
!! If the angle is small (in -pi/2..pi/2, or at least -pi..pi), as it 
!! should be if the basis function are ordered so that the minimal basis
!! goes in the small basis with just the projection (no rotation), then
!! there are much more efficient ways to compute this using tan in place 
!! of the angle (see Parlett, The symmetric eigenvalue problem, 6.8).
!! Hopefully it gets inlined.
  subroutine qs_right_compose_p_rot(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(qs_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,1)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='qs_block_rot_matrix_types:qs_right_compose_p_rot'
    failure=.false.
    if (qs_debug .and. debug_this_module) then ! this function must be fast
       QSPrecondition(i /= j,qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = rcos*matrix(:,i)+rsin*matrix(:,j)
       tmpV2 = -rsin*matrix(:,i)+rcos*matrix(:,j)
       matrix(:,i)=tmpV1
       matrix(:,j)=tmpV2
    end if
  end subroutine qs_right_compose_p_rot


!! calculates J*matrix where J is the plain rotation of angle angle
!! between i1 and j1 (less efficent than matrix*J because the slice
!! in continuos in memory).
!!v      i1    j1
!!v i1 (cos , -sin)
!!v j1 (sin ,  cos)
!! Quite inefficent.
!! If the angle is small (in -pi/2..pi/2, or at least -pi..pi), as it 
!! should be if the basis function are ordered so that the minimal basis
!! goes in the small basis with just the projection (no rotation), then
!! there are much more efficient ways to compute this using tan in place 
!! of the angle (see Parlett, The symmetric eigenvalue problem, 6.8).
!! Hopefully it gets inlined.
  subroutine qs_left_compose_p_rot(matrix,angle,i,j, error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(qs_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,2)) :: tmpV1,tmpV2
    character(len=*), parameter :: routineN='qs_block_rot_matrix_types:qs_left_compose_p_rot'
    failure=.false.
    if (debug_this_module) then ! this function must be fast
       QSPrecondition(i /= j,qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       rcos=cos(angle)
       rsin=sin(angle)
       tmpV1 = rcos*matrix(i,:)-rsin*matrix(j,:)
       tmpV2 = rsin*matrix(i,:)+rcos*matrix(j,:)
       matrix(i,:)=tmpV1
       matrix(j,:)=tmpV2
    end if
  end subroutine qs_left_compose_p_rot

!!****f* qs_pao_plain_rot/qs_right_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     qs_right_compose_p_rot_deriv
!!
!!   FUNCTION
!!     calculates matrix*dJ where dJ is the derivative of the plain 
!!     rotation of angle angle between i1 and j1.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_right_compose_p_rot_deriv(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(qs_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,1)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='qs_block_rot_matrix_types:qs_right_compose_p_rot_deriv'
    failure=.false.
    if (qs_debug .and. debug_this_module) then ! this function must be fast
       QSPrecondition(i /= j,qs_failure_level,routineN,error,failure)
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
  end subroutine qs_right_compose_p_rot_deriv
!!***
!***************************************************************************

!!****f* qs_pao_plain_rot/qs_left_compose_p_rot_deriv [1.0] *
!!
!!   NAME
!!     qs_left_compose_p_rot_deriv
!!
!!   FUNCTION
!!     calculates dJ*matrix where dJ is the derivative of the plain 
!!     rotation of angle angle between i1 and j1.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_left_compose_p_rot_deriv(matrix,angle,i,j,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix 
    real(kind=wp), intent(in) :: angle
    integer, intent(in) :: i,j
    type(qs_error_type), intent(inout), optional :: error

    logical :: failure
    real(kind=wp) :: rcos,rsin
    real(kind=wp), dimension(size(matrix,2)) :: tmpV1, tmpV2
    character(len=*), parameter :: routineN='qs_block_rot_matrix_types:qs_left_compose_p_rot_deriv'
    failure=.false.
    if (qs_debug .and. debug_this_module) then ! this function must be fast
       QSPrecondition(i /= j,qs_failure_level,routineN,error,failure)
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
  end subroutine qs_left_compose_p_rot_deriv
!!***
!***************************************************************************

end module qs_pao_plain_rot




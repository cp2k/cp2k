!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_plain_rot [1.0] *
!!
!!   NAME
!!     qs_plain_rot
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
  use qs_pao_types
  use qs_pao_utils, only: qs_sup_norm_matrix
  implicit none

  logical, parameter, private :: debug_this_module=.true.

  ! == core methods ==
!! initializes the object
  interface qs_init
     module procedure qs_init_pao_calc_u_du
  end interface

!! deallocates the memory used by the first argument
  interface qs_dealloc_ref
     module procedure qs_dealloc_calc_u_du_ref
  end interface

!! returns various attributes of the first argument 
  interface qs_get
     module procedure qs_calc_u_du_get
  end interface

!! iterator interface: returns true if there is a next element
  interface qs_next
     module procedure qs_calc_next_du
  end interface

!!***
!***************************************************************************

!!****s* qs_pao_plain_rot/qs_pao_calc_u_du_type [1.0] *
!!
!!   NAME
!!     qs_pao_calc_u_du_type
!!
!!   FUNCTION
!!     stores the parameters that should be bassed between the functions
!!     that calculates u du
!!
!!   NOTES
!!     I would have liked to reduce the calculation to one function
!!     that call a user given procedure to ude dU, but the impossibility to
!!     bind an argument, pass a generic argument or use a template forced 
!!     me to this solution.
!!     Both u and du need some calculation to be extracted, try not to
!!     extract them twice
!!
!!   ATTRIBUTES
!!     initialized true if initialized (to catch some nasty bugs)
!!     L: a matrix that contains the rotations u(1,1),u(2,1),... 
!!        u(min_bas-1,excl_bas).
!!     R: a matrix that contains the rotations u(min_bas+1,full_bas),
!!        u(min_bas+1,excl_bas),...,u(size_min_bas,size_excl_bas)
!!     tmp_m: a temporary matrix, used in the calculations
!!     min_bas: the actual position in the minimal basis
!!     excl_bas: the actual position in the full basis
!!     status: the actual status: (-2) after dealloc, (-1) after error,
!!             (0) just initialized, (1) normal status, (2) at end 
!!     local_angles: a pointer to the local angles that define U (SHARED)
!!     
!!     Read only (virtual) attributes
!!     u: the unitary matrix
!!     du: the derivative
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
  type qs_pao_calc_u_du_type
     private
     logical :: initialized
     real(kind=wp), dimension(:,:), pointer :: L,R,tmp_m
     integer :: min_bas, excl_bas, status
     type(qs_pao_local_angles), pointer :: local_angles
  end type qs_pao_calc_u_du_type
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
    type(qs_error), intent(inout), optional :: error

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
    type(qs_error), intent(inout), optional :: error

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
    type(qs_error), intent(inout), optional :: error

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
    type(qs_error), intent(inout), optional :: error

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

  ! ======= pao specific ======

!!****f* qs_pao_plain_rot/qs_pao_l_compose_u [1.0] *
!!
!!   NAME
!!     qs_pao_l_compose_u
!!
!!   FUNCTION
!!     calculates U*matrix where U is the unitary transformation defined by
!!     the given angles. the sequence is: smallbasis then full basis in 
!!     (local) growing order: u(1,1)*u(2,1)*...*u(size_small_bas,1)*...*
!!      u(size_small_bas,size_full_bas)*matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that will be modified
!!     local_angles: the angles that define U
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
  subroutine qs_pao_l_compose_u(matrix,local_angles,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix
    type(qs_pao_local_angles), intent(in) :: local_angles
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_pao_l_compose_u'
    integer :: n_min_bas, n_excl_bas, i_min_bas,i_excl_bas,n_full_bas
    failure=.false.

    call qs_get(local_angles,n_excl_bas=n_excl_bas,n_min_bas=n_min_bas,&
         n_full_bas=n_full_bas)
    QSPrecondition(n_min_bas>0,qs_warning_level,routineN,error,failure)
    QSPrecondition(n_excl_bas>=0,qs_warning_level,routineN,error,failure)
    QSPrecondition(size(matrix,1)==n_full_bas,qs_warning_level,routineN,error,failure)
    QSPrecondition(size(matrix,2)==n_full_bas,qs_warning_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       do i_excl_bas=n_excl_bas,1,-1
          do i_min_bas=n_min_bas,1,-1
             call qs_left_compose_p_rot(matrix=matrix,&
                  angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                  j=qs_get_non_p_index(local_angles,i=i_excl_bas,error=error),&
                  error=error)
          end do
       end do
    end if failureIf
  end subroutine qs_pao_l_compose_u
!!***
  !***************************************************************************

  ! ============== pao_calc_u_du ============

!!****f* qs_pao_plain_rot/qs_init_pao_calc_u_du [1.0] *
!!
!!   NAME
!!     qs_init_pao_calc_u_du
!!
!!   FUNCTION
!!     calculates the orthonormal matrix u and its derivative with respect 
!!     to the angles.
!!
!!   NOTES
!!     as the derivative are matrices, and can need quite some memory
!!     a user given function (take_derivative) is called with the derivative
!!     when it is available.
!!
!!   INPUTS
!!     c_env: the computation environement
!!     min_bas : the start index for the minimal basis (default 1,
!!               otherwise jumps over some values)
!!     excl_bas : the start index for the excluded basis (default 1, 
!!                otherwise jumps over some values)
!!     local_angles: the angles that define U
!!     u: the unitary transformation, if given it uses it and does not 
!!        recompute it
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
  subroutine qs_init_pao_calc_u_du(c_env,local_angles, min_bas, excl_bas,&
       u, error) !to do (correct full to excl)
    type(qs_pao_calc_u_du_type), intent(out) :: c_env
    type(qs_pao_local_angles), intent(in), target :: local_angles
    integer, optional, intent(in) :: min_bas,excl_bas
    real(kind=wp), dimension(:,:), optional, intent(in):: u
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_init_pao_calc_u_du'
    integer :: n_excl_bas,n_min_bas, stat, i_excl_bas, i_min_bas, i, n_full_bas
    failure=.false.
    c_env%initialized=.true.
    c_env%status=0
    c_env%local_angles => local_angles
    nullify(c_env%L, c_env%R, c_env%tmp_m)
    if (present(min_bas)) then
       c_env%min_bas=min_bas
    else
       c_env%min_bas=0
    end if
    if (present(excl_bas)) then
       c_env%excl_bas=excl_bas
    else
       c_env%excl_bas=0
    end if
    call qs_get(local_angles, n_excl_bas=n_excl_bas, n_min_bas=n_min_bas,&
         n_full_bas=n_full_bas, error=error)
    QSPrecondition(n_excl_bas>=0,qs_failure_level,routineN,error,failure)
    allocate(c_env%L(n_full_bas,n_full_bas),stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    allocate(c_env%R(n_full_bas,n_full_bas),stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    allocate(c_env%tmp_m(n_full_bas,n_full_bas),stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(u)) then
          if (present(min_bas)) then
             QSPrecondition(min_bas==1,qs_failure_level,routineN,error,failure)
          end if
          if (present(excl_bas)) then
             QSPrecondition(excl_bas==1,qs_failure_level,routineN,error,failure)
          end if
          c_env%L=0
          forall (i=1:n_excl_bas) c_env%L(i,i)=1
          c_env%R=u
          call qs_left_compose_p_rot(matrix=c_env%R,&
               angle=-qs_get_angle(local_angles,min_bas=n_min_bas,&
               excl_bas=n_excl_bas,error=error),&
               i=qs_get_p_index(local_angles,i=n_min_bas,error=error),&
               j=qs_get_non_p_index(local_angles,i=n_excl_bas,error=error),&
               error=error)
       else
          c_env%L=0
          c_env%R=0
          forall (i=1:n_full_bas)
             c_env%L(i,i)=1
             c_env%R(i,i)=1
          end forall
          ! def R
          do i_excl_bas=n_excl_bas,c_env%excl_bas+1,-1
             do i_min_bas=n_min_bas,1,-1
                call qs_left_compose_p_rot(matrix=c_env%R,&
                     angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                     j=qs_get_non_p_index(local_angles,&
                     i=i_excl_bas,error=error),&
                     error=error)
             end do
          end do
          i_excl_bas=c_env%excl_bas
          do i_min_bas=n_min_bas,c_env%min_bas+1,-1
             call qs_left_compose_p_rot(matrix=c_env%R,&
                  angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                  j=qs_get_non_p_index(local_angles,&
                  i=i_excl_bas, error=error),&
                  error=error)
          end do
          ! def L
          do i_excl_bas=1,c_env%excl_bas-1
             do i_min_bas=1,n_min_bas
                call qs_right_compose_p_rot(matrix=c_env%L,&
                     angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                     j=qs_get_non_p_index(local_angles,&
                     i=i_excl_bas,error=error),&
                     error=error)
             end do
          end do
          i_excl_bas=c_env%excl_bas
          do i_min_bas=1,c_env%min_bas-1
             call qs_right_compose_p_rot(matrix=c_env%L,&
                  angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                  j=qs_get_non_p_index(local_angles,&
                  i=i_excl_bas,error=error),&
                  error=error)
          end do
       end if
    end if failureIf
    if (failure) then
       c_env%status=-1
    end if
  end subroutine qs_init_pao_calc_u_du
!!***
  !***************************************************************************

!!****f* qs_pao_plain_rot/qs_calc_next_du [1.0] *
!!
!!   NAME
!!     qs_calc_next_du
!!
!!   FUNCTION
!!     calculates the next du, returns .true. if thre is a next du
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     c_env: the computation environement
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
  function qs_calc_next_du(c_env,error) result(res)
    logical :: res
    type(qs_pao_calc_u_du_type), intent(inout) :: c_env
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_calc_next_du'
    real(kind=wp) :: err
    integer :: old_min_bas,old_excl_bas, i
    failure=.false.

    QSPrecondition(c_env%initialized,qs_failure_level,routineN,error,failure)
    if (debug_this_module .and. qs_debug) then
       QSPrecondition(associated(c_env%R),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(c_env%L),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       old_min_bas=c_env%min_bas
       old_excl_bas=c_env%excl_bas
       if (c_env%min_bas>=qs_get_min_bas_size(c_env%local_angles)) then
          if (c_env%excl_bas>=qs_get_excl_bas_size(c_env%local_angles)) then
             c_env%status=2
             res=.false.
             return
          else if (qs_get_excl_bas_size(c_env%local_angles)==0) then
             c_env%status=2
             res=.false.
             return
          else
             c_env%min_bas=1
             c_env%excl_bas=c_env%excl_bas+1
          end if
       else
          c_env%min_bas=c_env%min_bas+1
       end if
       select case (c_env%status)
       case (-2) ! status after dealloc
          res=.false.
       case (-1) ! status after an error
          res=.false.
       case (0) ! initial status
          res=.true.
          c_env%status=1
       case (1)
          ! normal case
          res=.true.
          call qs_right_compose_p_rot(matrix=c_env%L,&
               angle=qs_get_angle(c_env%local_angles,min_bas=old_min_bas,&
               excl_bas=old_excl_bas), &
               i=qs_get_p_index(c_env%local_angles,old_min_bas,error=error), &
               j=qs_get_non_p_index(c_env%local_angles,old_excl_bas,&
               error=error), error=error)
          call qs_left_compose_p_rot(matrix=c_env%R,&
               angle= -qs_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
               excl_bas=c_env%excl_bas), &
               i=qs_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=qs_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error))
       case (2)
          res=.false.
          if (debug_this_module .and. qs_debug) then
             forall (i=1:size(c_env%R)) c_env%R(i,i)=c_env%R(i,i)-1
             err=qs_sup_norm_matrix(c_env%R)
             if (err>2*size(c_env%R)*epsilon(0.0_wp)) then
                call qs_error_message(qs_warning_level,routineN,&
                     message="accumulation of the error bigger than tolerance:"//&
                     qs_to_string(err)//" in "//QSSourceFileRef(),&
                     error=error)
             end if
          end if
       case default
          res=.false.
          call qs_error_message(qs_warning_level,routineN,&
               message="unknown status "//qs_to_string(c_env%status)//&
               " in "//QSSourceFileRef(),&
               error=error)
          c_env%status=-1
       end select
    end if
  end function qs_calc_next_du
!!***
  !***************************************************************************

!!****f* qs_pao_plain_rot/qs_calc_u_du_get [1.0] *
!!
!!   NAME
!!     qs_calc_u_du_get
!!
!!   FUNCTION
!!     returns various attributes of the calculation
!!
!!   NOTES
!!     most attributes require some calculation, try not to get them
!!     twice
!!
!!   INPUTS
!!     c_env: the computation environement
!!     u: the unitary matrix
!!     du: the derivative
!!     min_bas : the minbas index of the angle that has been derived
!!     excl_bas : the exclbas index of the angle that has been derived
!!     local_angles: the angles that define U
!!     finished: is true if the iterator is past the last element
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
  subroutine qs_calc_u_du_get(c_env,u,du,min_bas,excl_bas,local_angles,&
       finished,error)
    type(qs_pao_calc_u_du_type), intent(inout) :: c_env
    type(qs_pao_local_angles), optional, pointer :: local_angles
    integer, optional, intent(out) :: min_bas,excl_bas
    logical,optional, intent(out) :: finished
    real(kind=wp), optional, intent(out), dimension(:,:) :: u,du
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_calc_u_du_get'
    integer :: n_excl_bas, n_full_bas
    real(kind=wp) :: angle
    failure=.false.

    QSPrecondition(c_env%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(c_env%status==1,qs_failure_level,routineN,error,failure)
    if (debug_this_module .and. qs_debug) then
       QSPrecondition(associated(c_env%local_angles),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(c_env%R),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(c_env%L),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(c_env%tmp_m),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       call qs_get(c_env%local_angles,n_excl_bas=n_excl_bas,&
            n_full_bas=n_full_bas)
       angle=qs_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
            excl_bas=c_env%excl_bas,error=error)
       if (present(u)) then
          QSPrecondition(size(u,1)==n_full_bas,qs_failure_level,routineN,error,failure)
          QSPrecondition(size(u,2)==n_full_bas,qs_failure_level,routineN,error,failure)
          c_env%tmp_m=c_env%L
          call qs_right_compose_p_rot(matrix=c_env%tmp_m,&
               angle=angle, i=c_env%min_bas, j=c_env%excl_bas,error=error)
          call dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_wp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_wp,u,n_full_bas)
       end if
       if (present(du)) then
          QSPrecondition(size(du,1)==n_full_bas,qs_failure_level,routineN,error,failure)
          QSPrecondition(size(du,2)==n_full_bas,qs_failure_level,routineN,error,failure)
          c_env%tmp_m=c_env%L
          call qs_right_compose_p_rot_deriv(matrix=c_env%tmp_m,&
               angle=angle, i=c_env%min_bas, j=c_env%excl_bas,error=error)
          call dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_wp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_wp,u,n_full_bas)
       end if
       if (present(min_bas)) &
            min_bas=c_env%min_bas
       if (present(excl_bas)) &
            excl_bas=c_env%excl_bas
       if (present(local_angles)) &
            local_angles => c_env%local_angles
       if (present(finished)) then
          finished=c_env%status /= 1 .and. c_env%status /= 0
          ! catch also status=0 and past end?
       end if
    end if
  end subroutine qs_calc_u_du_get
!!***
  !***************************************************************************

!!****f* qs_pao_plain_rot/qs_dealloc_calc_u_du_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_calc_u_du_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by calc_u_du
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     c_env: the calc_u_du object that should be released
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
  subroutine qs_dealloc_calc_u_du_ref(c_env, error)
    type(qs_pao_calc_u_du_type), intent(inout) :: c_env
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_dealloc_calc_u_du_ref'
    integer :: stat
    failure=.false.
    QSAssert(c_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (associated(c_env%R)) then
          deallocate(c_env%R, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       if (associated(c_env%L)) then
          deallocate(c_env%L, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       if (associated(c_env%tmp_m)) then
          deallocate(c_env%tmp_m, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       nullify(c_env%local_angles)
       c_env%status=-2
    end if
  end subroutine qs_dealloc_calc_u_du_ref
!!***
  !***************************************************************************

end module qs_pao_plain_rot




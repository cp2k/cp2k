!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_unitary_t [1.0] *
!!
!!   NAME
!!     pao_unitary_t
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module pao_unitary_t
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_log, cp_logger_type, cp_note_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_error_dealloc_ref, cp_error_get_logger
  use kinds, only: wp=>dp
  use pao_plain_rot, only: cp_right_compose_p_rot, &
       cp_right_compose_p_rot_deriv, cp_left_compose_p_rot,&
       cp_left_compose_p_rot_deriv
  use pao_types, only: pao_env_type, pao_glob_angles_type
  use pao_l_angles_methods, only: pao_local_angles_type, cp_init, &
       pao_l_ang_dealloc_ref, cp_get, cp_get_angle, cp_get_p_index,&
       cp_get_non_p_index, cp_get_size_min_bas, cp_get_size_excl_bas,&
       cp_get_size_full_bas
  use cp_matrix_utils, only: cp_sup_norm_matrix
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.

  public :: pao_calc_u_du_type

  public :: cp_init, cp_dealloc_ref, cp_get, cp_next

  public :: pao_calc_u_du_init, pao_calc_u_du_dealloc_ref,&
       pao_calc_u_du_get, pao_calc_u_du_next

! == core methods ==
!! initializes the object
  interface cp_init
     module procedure pao_calc_u_du_init
  end interface

!! deallocates the memory used by the first argument
  interface cp_dealloc_ref
     module procedure pao_calc_u_du_dealloc_ref
  end interface

!! returns various attributes of the first argument 
  interface cp_get
     module procedure pao_calc_u_du_get
  end interface

!! iterator interface: returns true if there is a next element
  interface cp_next
     module procedure pao_calc_u_du_next
  end interface

!!***
!***************************************************************************

!!****s* pao_plain_rot/pao_calc_u_du_type [1.0] *
!!
!!   NAME
!!     pao_calc_u_du_type
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
  type pao_calc_u_du_type
     private
     logical :: initialized
     real(kind=wp), dimension(:,:), pointer :: L,R,tmp_m
     integer :: min_bas, excl_bas, status
     type(pao_local_angles_type), pointer :: local_angles
  end type pao_calc_u_du_type
!!***
!****************************************************************************

contains
! ======= pao specific ======

! ============== pao_calc_u_du ============

!!****f* pao_plain_rot/pao_calc_u_du_init [1.0] *
!!
!!   NAME
!!     pao_calc_u_du_init
!!
!!   SYNOPSIS
!!     Subroutine pao_calc_u_du_init(c_env, local_angles, min_bas, excl_bas,&
!!         u, error)
!!       Type(pao_calc_u_du_type), Intent (OUT):: c_env
!!       Type(pao_local_angles_type), Intent (INOUT), Target:: local_angles
!!       Integer, Optional, Intent (IN):: min_bas, excl_bas
!!       Real(Kind=wp), Dimension(:,:), Optional, Intent (IN):: u
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_calc_u_du_init
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
!!            see module cp_error_handling 
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
  subroutine pao_calc_u_du_init(c_env,local_angles, min_bas, excl_bas,&
       u, error) !to do (correct full to excl)
    type(pao_calc_u_du_type), intent(out) :: c_env
    type(pao_local_angles_type), intent(inout), target :: local_angles
    integer, optional, intent(in) :: min_bas,excl_bas
    real(kind=wp), dimension(:,:), optional, intent(in):: u
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot:pao_calc_u_du_init'
    integer :: n_excl_bas,n_min_bas, stat, i_excl_bas, i_min_bas, i, n_full_bas
    integer, dimension(:), pointer :: proj_indexes, non_proj_indexes
    failure=.false.
    c_env%initialized=.true.
    c_env%status=0
    c_env%local_angles => local_angles
    nullify(c_env%L, c_env%R, c_env%tmp_m)
    if (present(min_bas)) then
       c_env%min_bas=min_bas
    else
       c_env%min_bas=1
    end if
    if (present(excl_bas)) then
       c_env%excl_bas=excl_bas
    else
       c_env%excl_bas=1
    end if
    call cp_get(local_angles, size_excl_bas=n_excl_bas, &
         size_min_bas=n_min_bas,&
         size_full_bas=n_full_bas, proj_indexes=proj_indexes,&
         non_proj_indexes=non_proj_indexes, error=error)
    CPPrecondition(n_excl_bas>=0,cp_failure_level,routineP,error,failure)
    allocate(c_env%L(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    allocate(c_env%R(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    allocate(c_env%tmp_m(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(u)) then
          if (present(min_bas)) then
             CPPrecondition(min_bas==1,cp_failure_level,routineP,error,failure)
          end if
          if (present(excl_bas)) then
             CPPrecondition(excl_bas==1,cp_failure_level,routineP,error,failure)
          end if
          c_env%L=0
          forall (i=1:n_excl_bas) c_env%L(i,i)=1
          c_env%R=u
          call cp_left_compose_p_rot(matrix=c_env%R,&
               angle=-cp_get_angle(local_angles,min_bas=n_min_bas,&
               excl_bas=n_excl_bas,error=error),&
               i=cp_get_p_index(local_angles,min_bas=n_min_bas,error=error),&
               j=cp_get_non_p_index(local_angles,excl_bas=n_excl_bas,error=error),&
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
                call cp_left_compose_p_rot(matrix=c_env%R,&
                     angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=proj_indexes(i_min_bas),&
                     j=non_proj_indexes(i_excl_bas),&
                     error=error)
             end do
          end do
          i_excl_bas=c_env%excl_bas
          do i_min_bas=n_min_bas,c_env%min_bas+1,-1
             call cp_left_compose_p_rot(matrix=c_env%R,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          end do
! def L
          do i_excl_bas=1,c_env%excl_bas-1
             do i_min_bas=1,n_min_bas
                call cp_right_compose_p_rot(matrix=c_env%L,&
                     angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=proj_indexes(i_min_bas),&
                     j=non_proj_indexes(i_excl_bas),&
                     error=error)
             end do
          end do
          i_excl_bas=c_env%excl_bas
          do i_min_bas=1,c_env%min_bas-1
             call cp_right_compose_p_rot(matrix=c_env%L,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          end do
       end if
       c_env%min_bas=c_env%min_bas-1 ! to start up goes one before
    end if failureIf
    if (failure) then
       c_env%status=-1
    end if
  end subroutine pao_calc_u_du_init
!!***
!***************************************************************************

!!****f* pao_plain_rot/pao_calc_u_du_next [1.0] *
!!
!!   NAME
!!     pao_calc_u_du_next
!!
!!   SYNOPSIS
!!     Function pao_calc_u_du_next(c_env, error) Result(res)
!!       Logical:: res
!!       Type(pao_calc_u_du_type), Intent (INOUT):: c_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_calc_u_du_next
!!
!!   FUNCTION
!!     calculates the next du, returns .true. if there is a next du
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     c_env: the computation environement
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
!!   SOURCE
!***************************************************************************
  function pao_calc_u_du_next(c_env,error) result(res)
    logical :: res
    type(pao_calc_u_du_type), intent(inout) :: c_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot:pao_calc_u_du_next'
    real(kind=wp) :: err
    integer :: old_min_bas,old_excl_bas, i
    failure=.false.

    CPPrecondition(c_env%initialized,cp_failure_level,routineP,error,failure)
    if (debug_this_module .and. cp_debug) then
       CPPrecondition(associated(c_env%R),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(c_env%L),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       old_min_bas=c_env%min_bas
       old_excl_bas=c_env%excl_bas
       if (c_env%min_bas>=cp_get_size_min_bas(c_env%local_angles)) then
          if (c_env%excl_bas>=cp_get_size_excl_bas(c_env%local_angles)) then
             c_env%status=2
             res=.false.
             return
          else if (cp_get_size_excl_bas(c_env%local_angles)==0) then
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
          call cp_right_compose_p_rot(matrix=c_env%L,&
               angle=cp_get_angle(c_env%local_angles,min_bas=old_min_bas,&
               excl_bas=old_excl_bas), &
               i=cp_get_p_index(c_env%local_angles,old_min_bas,error=error), &
               j=cp_get_non_p_index(c_env%local_angles,old_excl_bas,&
               error=error), error=error)
          call cp_left_compose_p_rot(matrix=c_env%R,&
               angle= -cp_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
               excl_bas=c_env%excl_bas), &
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error))
       case (2)
          res=.false.
          if (debug_this_module .and. cp_debug) then
             do i=1,size(c_env%R)
                c_env%R(i,i)=c_env%R(i,i)-1
             end do
             err=cp_sup_norm_matrix(c_env%R)
             if (err>2*size(c_env%R)*epsilon(0.0_wp)) then
                call cp_error_message(cp_warning_level,routineP,&
                     message="accumulation of the error bigger than tolerance:"//&
                     cp_to_string(err)//" in "//&
                     CPSourceFileRef,&
                     error=error)
             end if
          end if
       case default
          res=.false.
          call cp_error_message(cp_warning_level,routineP,&
               message="unknown status "//cp_to_string(c_env%status)//&
               " in "//CPSourceFileRef,&
               error=error)
          c_env%status=-1
       end select
    end if
  end function pao_calc_u_du_next
!!***
!***************************************************************************

!!****f* pao_plain_rot/pao_calc_u_du_get [1.0] *
!!
!!   NAME
!!     pao_calc_u_du_get
!!
!!   SYNOPSIS
!!     Subroutine pao_calc_u_du_get(c_env, u, du, min_bas, excl_bas,&
!!         local_angles, finished, error)
!!       Type(pao_calc_u_du_type), Intent (INOUT):: c_env
!!       Type(pao_local_angles_type), Optional, Pointer:: local_angles
!!       Integer, Optional, Intent (OUT):: min_bas, excl_bas
!!       Logical, Optional, Intent (OUT):: finished
!!       Real(Kind=wp), Optional, Intent (OUT), Dimension(:,:):: u, du
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_calc_u_du_get
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
!!            see module cp_error_handling 
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
  subroutine pao_calc_u_du_get(c_env,u,du,min_bas,excl_bas,local_angles,&
       finished,error)
    type(pao_calc_u_du_type), intent(inout) :: c_env
    type(pao_local_angles_type), optional, pointer :: local_angles
    integer, optional, intent(out) :: min_bas,excl_bas
    logical,optional, intent(out) :: finished
    real(kind=wp), optional, intent(out), dimension(:,:) :: u,du
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot:pao_calc_u_du_get'
    integer :: n_excl_bas, n_full_bas
    real(kind=wp) :: angle
    failure=.false.

    CPPrecondition(c_env%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(c_env%status==1,cp_failure_level,routineP,error,failure)
    if (debug_this_module .and. cp_debug) then
       CPPrecondition(associated(c_env%local_angles),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(c_env%R),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(c_env%L),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(c_env%tmp_m),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       call cp_get(c_env%local_angles,size_excl_bas=n_excl_bas,&
            size_full_bas=n_full_bas)
       angle=cp_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
            excl_bas=c_env%excl_bas,error=error)
       if (present(u)) then
          CPPrecondition(size(u,1)==n_full_bas,cp_failure_level,routineP,error,failure)
          CPPrecondition(size(u,2)==n_full_bas,cp_failure_level,routineP,error,failure)
          c_env%tmp_m=c_env%L
          call cp_right_compose_p_rot(matrix=c_env%tmp_m,&
               angle=angle,&
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error),error=error)
          call dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_wp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_wp,u,n_full_bas)
       end if
       if (present(du)) then
          CPPrecondition(size(du,1)==n_full_bas,cp_failure_level,routineP,error,failure)
          CPPrecondition(size(du,2)==n_full_bas,cp_failure_level,routineP,error,failure)
          c_env%tmp_m=c_env%L
          call cp_right_compose_p_rot_deriv(matrix=c_env%tmp_m,&
               angle=angle,&
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error),error=error)
          call dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_wp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_wp,du,n_full_bas)
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
  end subroutine pao_calc_u_du_get
!!***
!***************************************************************************

!!****f* pao_plain_rot/pao_calc_u_du_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_calc_u_du_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_calc_u_du_dealloc_ref(c_env, error)
!!       Type(pao_calc_u_du_type), Intent (INOUT):: c_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_calc_u_du_dealloc_ref
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
!!            see module cp_error_handling 
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
  subroutine pao_calc_u_du_dealloc_ref(c_env, error)
    type(pao_calc_u_du_type), intent(inout) :: c_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot:pao_calc_u_du_dealloc_ref'
    integer :: stat
    failure=.false.
    CPAssert(c_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (associated(c_env%R)) then
          deallocate(c_env%R, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       if (associated(c_env%L)) then
          deallocate(c_env%L, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       if (associated(c_env%tmp_m)) then
          deallocate(c_env%tmp_m, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       nullify(c_env%local_angles)
       c_env%status=-2
    end if
  end subroutine pao_calc_u_du_dealloc_ref
!!***
!***************************************************************************


end module pao_unitary_t

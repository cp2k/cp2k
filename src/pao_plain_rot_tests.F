!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_plain_rot_tests [1.0] *
!!
!!   NAME
!!     cp_plain_rot_tests
!!
!!   FUNCTION
!!     test plain rotations (module cp_plain_rot)
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
!*****************************************************************************
module pao_plain_rot_tests
  use cp_log_handling
  use cp_error_handling
  use cp_output_handling
  use kinds, only: wp=>dp
  use pao_utils, only:  cp_get_random_integer
  use cp_matrix_utils, only: cp_sup_norm_matrix
  use mathconstants, only: pi
  use global_types, only: global_environment_type
  use pao_plain_rot
  use pao_types
  use pao_types_tests, only: cp_init_test_pao_env
  use pao_proj_methods
  use pao_env_methods
  use pao_glob_angles_methods
  use pao_l_angles_methods
  use cp_array_utils
  !use pao_matrix_methods
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.

  public :: pao_test_plain_rot, pao_test_u_du
!!***
!******************************************************************************
contains

!!****f* pao_plain_rot_tests/pao_test_plain_rot [1.0] *
!!
!!   NAME
!!     pao_test_plain_rot
!!
!!   SYNOPSIS
!!     Subroutine pao_test_plain_rot(g_env, failed, error)
!!       Type(global_environment_type), Intent (IN), Target:: g_env
!!       Logical, Intent (INOUT):: failed
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_test_plain_rot
!!
!!   FUNCTION
!!     testing the plain rotation module
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     g_env: global environement (for the future)
!!     failure: an in out value that become true if some tests fail
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** ************************************************************************
  subroutine pao_test_plain_rot(g_env, failed, error)
    type(global_environment_type), intent(in), target ::g_env
    logical, intent(inout) :: failed
    type(cp_error_type), optional, intent(inout):: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_plain_rot_tests:cp_test_plain_rot'
    integer, parameter :: n=10
    integer :: i,j
    real(kind=wp) :: m1(n,n),m2(n,n),tmpV1(n), angle, err
    type(cp_logger_type), pointer :: logger
    failure=.FALSE.
    CPPrecondition(n>1,cp_failure_level,routineP,error,failure)
    if (.NOT. failure) then
       logger=>cp_error_get_logger(error)

       i=cp_get_random_integer(from=1,to=n)
       j=cp_get_random_integer(from=1,to=n)
       if (i==j) then
          if (i>1) then
             j=i-1
          else
             j=i+1
          end if
       end if
       call random_number(m1)
       m2=m1

       ! t1r
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test1r: rotation with pi/2 swaps columns and change sign of the second')
       end if
       tmpV1=m2(:,i)
       m2(:,i)=m2(:,j)
       m2(:,j)=-tmpV1
       call cp_right_compose_p_rot(m1,pi/2,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       if (err>10*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test1r (right rot pi/2) with i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test1r: SUCCESS')
       end if

       ! t2r
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test2r: rotation with -angle restores the matrix')
       end if
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call cp_right_compose_p_rot(m1,angle,i,j)
       call cp_right_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test2r (right -angle = rot^-1) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test2r: SUCCESS')
       end if

       !t1l
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test1l: rotation with pi/2 swaps lines and change sign of the first')
       end if
       tmpV1=m2(i,:)
       m2(i,:)=-m2(j,:)
       m2(j,:)=tmpV1
       call cp_left_compose_p_rot(m1,pi/2,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test1l (right rot pi/2) with i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test1l: SUCCESS')
       end if

       ! t2l
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test2l: rotation with -angle restores the matrix')
       end if
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call cp_left_compose_p_rot(m1,angle,i,j)
       call cp_left_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test2l (right -angle = rot^-1) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test2l: SUCCESS')
       end if

       ! t3 
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test3: the unit matrix is the same in all orthonormal basis')
       end if
       m1=0
       forall (i=1:n) m1(i,i)=1
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call cp_left_compose_p_rot(m1,angle,i,j)
       call cp_right_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test3 (unity matrix same in all ONB) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          if (cp_debug .and. debug_this_module) then
             call cp_2d_r_output(logger, 'identity_matrix', routineP,1,m1,'identity matrix')
          end if
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test3: SUCCESS')
       end if

       ! t4
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test4: d/dt (J(t,i,j)J(-t,i,j)M)=0')
       end if
       ! test m1 = identity
       !m1=0
       !forall (i=1:n) m1(i,i)=1
       call random_number(m1)
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call cp_left_compose_p_rot(m1,-angle,i,j)
       !call cp_2d_r_output(logger, 'm1_1', routineP,1,m1,'m1.1')
       call cp_left_compose_p_rot_deriv(m1,angle,i,j)
       !call cp_2d_r_output(logger, 'm1_2', routineP,1,m1,'m1.2')
       call cp_left_compose_p_rot_deriv(m2,-angle,i,j)
       !call cp_2d_r_output(logger, 'm2_1', routineP,1,m2,'m2.1')
       call cp_left_compose_p_rot(m2,angle,i,j)
       !call cp_2d_r_output(logger, 'm2_2', routineP,1,m2,'m2.2')
       err=cp_sup_norm_matrix(m1-m2)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test4 d/dt (J(t,i,j)J(-t,i,j)M)=0 with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          if (cp_debug .and. debug_this_module) then
             call cp_2d_r_output(logger, 'zero_matrix', routineP,1,m1-m2,'zero matrix 1')
          end if
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test4: SUCCESS')
       end if

       ! t5
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test5: d/dt (MJ(t,i,j)J(-t,i,j))=0')
       end if
       call random_number(m1)
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call cp_right_compose_p_rot(m1,-angle,i,j)
       call cp_right_compose_p_rot_deriv(m1,angle,i,j)
       call cp_right_compose_p_rot_deriv(m2,-angle,i,j)
       call cp_right_compose_p_rot(m2,angle,i,j)
       err=cp_sup_norm_matrix(m1-m2)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test5 d/dt (M J(t,i,j)J(-t,i,j))=0 with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          if (cp_debug .and. debug_this_module) then
             call cp_2d_r_output(logger, 'zero_matrix', routineP,2,m1-m2,'zero matrix 2')
          end if
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test5: SUCCESS')
       end if

       ! u_du tests
       if (.not.pao_test_u_du(g_env,error=error)) then
          failure=.true.
       end if

       if (cp_debug .and. debug_this_module .and. .not.failure) then
          call cp_log(logger,cp_warning_level,routineP,&
               'cp_test_plain_rot: all tests succeded, SUCCESS')
       end if
    end if
  end subroutine pao_test_plain_rot
  !******************************************************************************

!!****f* pao_test_plain_rot/pao_test_u_du [1.0] *
!!
!!   NAME
!!     pao_test_u_du
!!
!!   SYNOPSIS
!!     Function pao_test_u_du(glob_env, error) Result(res)
!!       Logical:: res
!!       Type(global_environment_type), Target:: glob_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_test_u_du
!!
!!   FUNCTION
!!     tests u generation and du
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_env: the global environement
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
  function pao_test_u_du(glob_env,error) result(res)
    logical ::res
    type(global_environment_type), target :: glob_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot_tests:pao_test_u_du'
    integer, parameter :: n_atom=3, atom_nr=1
    type(pao_env_type) :: pao_env
    type(pao_glob_angles_type), pointer :: g_angles
    type(pao_local_angles_type) :: l_angles
    real(kind=wp), dimension(:), pointer :: angles
    real(kind=wp),dimension(:,:), allocatable :: m1,m2,m3
    real(kind=wp) :: err
    integer :: stat, full_bas_size, i
    type(cp_logger_type), pointer :: logger
    failure=.false.
    logger=>cp_error_get_logger(error)

    call cp_init_test_pao_env(pao_env,glob_env,n_atom,error=error)
    call cp_get(pao_env, angles_att=g_angles,error=error)
    call cp_get(g_angles,angles=angles)
    call random_number(angles)
    angles=4.0_wp*pi*angles
    call cp_init(l_angles,glob_angles=g_angles, atom_nr=atom_nr)
    CPPostcondition(cp_valid(l_angles),cp_failure_level,routineP,error,failure)
    full_bas_size=cp_get_full_bas_size(l_angles)
    allocate(m1(full_bas_size,full_bas_size),m2(full_bas_size,full_bas_size),&
         m3(full_bas_size,full_bas_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (cp_debug.and.debug_this_module) then
          call cp_log(logger,cp_note_level,routineP, &
               'test u_du 1: l_compose_u, U^* U = 1')
       end if
       m1=0.0_wp
       forall (i=1:size(m1,1)) m1(i,i)=1.0_wp
       call pao_l_compose_u(m1,l_angles,error=error)
       m3=m1
       call dgemm('T','N',full_bas_size,full_bas_size,full_bas_size,&
            1.0_wp,m1,full_bas_size,m3,full_bas_size,&
            0.0_wp,m2,full_bas_size)
       forall (i=1:full_bas_size) m2(i,i)=m2(i,i)-1.0_wp
       err=cp_sup_norm_matrix(m2)
       if (err>2*full_bas_size*epsilon(0.0_wp)) then
          failure=.true.
          call cp_log(logger,cp_failure_level,routineP, &
               'test u_du 1: l_compose_u, U^* U = 1, failed')
          if (cp_debug .and. debug_this_module) then
             call cp_2d_r_output(logger, 'zero_matrix', routineP,3,m2,'zero matrix l_compose_u, U^* U - 1')
             call cp_1d_r_output(logger, 'angles', routineP,3,&
                  angles(1:(full_bas_size*cp_get_min_bas_size(l_angles))),&
                  'angles')
          end if
       else if (cp_debug .and. debug_this_module) then
          call cp_log(logger,cp_warning_level,routineP,&
               'test u_du 1: SUCCESS')
       end if
    end if
    deallocate(m1,m2,stat=stat) ! just to be sure...
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    call cp_dealloc_ref(l_angles,error=error)
    call cp_dealloc_ref(pao_env,error=error)
    res=failure
  end function pao_test_u_du
  !***************************************************************************

end module pao_plain_rot_tests


!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_plain_rot_tests [1.0] *
!!
!!   NAME
!!     qs_plain_rot_tests
!!
!!   FUNCTION
!!     test plain rotations (module qs_plain_rot)
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
module qs_pao_plain_rot_tests
  use qs_log_handling
  use qs_error_handling
  use qs_output_handling
  use kinds, only: wp=>dp
  use qs_pao_utils, only:  qs_get_random_integer, qs_sup_norm_matrix
  use mathconstants, only: pi
  use global_types, only: global_environment_type
  use qs_pao_plain_rot
  use qs_pao_types
  use qs_pao_types_tests, only: qs_init_test_pao_env
  implicit none

  logical, private, parameter :: debug_this_module=.true.

!!***
!******************************************************************************
contains

!!****f* qs_pao_plain_rot_tests/qs_test_plain_rot [1.0] *
!!
!!   NAME
!!     qs_pao_test_plain_rot
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
!!            see module qs_error_handling 
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
  subroutine qs_pao_test_plain_rot(g_env, failed, error)
    type(global_environment_type), intent(in), target ::g_env
    logical, intent(inout) :: failed
    type(qs_error_type), optional, intent(inout):: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_plain_rot_tests:qs_test_plain_rot'
    integer, parameter :: n=10
    integer :: i,j
    real(kind=wp) :: m1(n,n),m2(n,n),tmpV1(n), angle, err
    type(qs_logger), pointer :: logger
    failure=.FALSE.
    QSPrecondition(n>1,qs_failure_level,routineN,error,failure)
    if (.NOT. failure) then
       logger=>qs_get_error_logger(error)

       i=qs_get_random_integer(from=1,to=n)
       j=qs_get_random_integer(from=1,to=n)
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
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test1r: rotation with pi/2 swaps columns and change sign of the second')
       end if
       tmpV1=m2(:,i)
       m2(:,i)=m2(:,j)
       m2(:,j)=-tmpV1
       call qs_right_compose_p_rot(m1,pi/2,i,j)
       err=qs_sup_norm_matrix(m2-m1)
       if (err>10*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test1r (right rot pi/2) with i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test1r: SUCCESS')
       end if

       ! t2r
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test2r: rotation with -angle restores the matrix')
       end if
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call qs_right_compose_p_rot(m1,angle,i,j)
       call qs_right_compose_p_rot(m1,-angle,i,j)
       err=qs_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test2r (right -angle = rot^-1) with angle='//qs_to_string(angle)&
               //' i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test2r: SUCCESS')
       end if

       !t1l
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test1l: rotation with pi/2 swaps lines and change sign of the first')
       end if
       tmpV1=m2(i,:)
       m2(i,:)=-m2(j,:)
       m2(j,:)=tmpV1
       call qs_left_compose_p_rot(m1,pi/2,i,j)
       err=qs_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test1l (right rot pi/2) with i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test1l: SUCCESS')
       end if

       ! t2l
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test2l: rotation with -angle restores the matrix')
       end if
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call qs_left_compose_p_rot(m1,angle,i,j)
       call qs_left_compose_p_rot(m1,-angle,i,j)
       err=qs_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test2l (right -angle = rot^-1) with angle='//qs_to_string(angle)&
               //' i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test2l: SUCCESS')
       end if

       ! t3 
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test3: the unit matrix is the same in all orthonormal basis')
       end if
       m1=0
       forall (i=1:n) m1(i,i)=1
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call qs_left_compose_p_rot(m1,angle,i,j)
       call qs_right_compose_p_rot(m1,-angle,i,j)
       err=qs_sup_norm_matrix(m2-m1)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test3 (unity matrix same in all ONB) with angle='//qs_to_string(angle)&
               //' i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
          if (qs_debug .and. debug_this_module) then
             call qs_output_2d_r_array(logger, 'identity_matrix', routineN,1,m1,'identity matrix')
          end if
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test3: SUCCESS')
       end if

       ! t4
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test4: d/dt (J(t,i,j)J(-t,i,j)M)=0')
       end if
       ! test m1 = identity
       !m1=0
       !forall (i=1:n) m1(i,i)=1
       call random_number(m1)
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call qs_left_compose_p_rot(m1,-angle,i,j)
       !call qs_output_2d_r_array(logger, 'm1_1', routineN,1,m1,'m1.1')
       call qs_left_compose_p_rot_deriv(m1,angle,i,j)
       !call qs_output_2d_r_array(logger, 'm1_2', routineN,1,m1,'m1.2')
       call qs_left_compose_p_rot_deriv(m2,-angle,i,j)
       !call qs_output_2d_r_array(logger, 'm2_1', routineN,1,m2,'m2.1')
       call qs_left_compose_p_rot(m2,angle,i,j)
       !call qs_output_2d_r_array(logger, 'm2_2', routineN,1,m2,'m2.2')
       err=qs_sup_norm_matrix(m1-m2)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test4 d/dt (J(t,i,j)J(-t,i,j)M)=0 with angle='//qs_to_string(angle)&
               //' i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
          if (qs_debug .and. debug_this_module) then
             call qs_output_2d_r_array(logger, 'zero_matrix', routineN,1,m1-m2,'zero matrix 1')
          end if
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test4: SUCCESS')
       end if

       ! t5
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test5: d/dt (MJ(t,i,j)J(-t,i,j))=0')
       end if
       call random_number(m1)
       m2=m1
       call random_number(angle)
       angle=4.0_wp*pi*angle
       call qs_right_compose_p_rot(m1,-angle,i,j)
       call qs_right_compose_p_rot_deriv(m1,angle,i,j)
       call qs_right_compose_p_rot_deriv(m2,-angle,i,j)
       call qs_right_compose_p_rot(m2,angle,i,j)
       err=qs_sup_norm_matrix(m1-m2)
       if (err>2*n*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test5 d/dt (M J(t,i,j)J(-t,i,j))=0 with angle='//qs_to_string(angle)&
               //' i='//qs_to_string(i)// &
               ' j='//qs_to_string(j)//' failed')
          if (qs_debug .and. debug_this_module) then
             call qs_output_2d_r_array(logger, 'zero_matrix', routineN,2,m1-m2,'zero matrix 2')
          end if
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test5: SUCCESS')
       end if

       ! u_du tests
       if (.not.qs_pao_test_u_du(g_env,error=error)) then
          failure=.true.
       end if

       if (qs_debug .and. debug_this_module .and. .not.failure) then
          call qs_log(logger,qs_warning_level,routineN,&
               'qs_test_plain_rot: all tests succeded, SUCCESS')
       end if
    end if
  end subroutine qs_pao_test_plain_rot
!!***
  !******************************************************************************

!!****f* qs_pao_test_plain_rot/qs_pao_test_u_du [1.0] *
!!
!!   NAME
!!     qs_pao_test_u_du
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
  function qs_pao_test_u_du(glob_env,error) result(res)
    logical ::res
    type(global_environment_type), target :: glob_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot_tests:qs_pao_test_u_du'
    integer, parameter :: n_atom=3, atom_nr=1
    type(qs_pao_env_type) :: pao_env
    type(qs_pao_glob_angles_type), pointer :: g_angles
    type(qs_pao_local_angles_type) :: l_angles
    real(kind=wp), dimension(:), pointer :: angles
    real(kind=wp),dimension(:,:), allocatable :: m1,m2,m3
    real(kind=wp) :: err
    integer :: stat, full_bas_size, i
    type(qs_logger), pointer :: logger
    failure=.false.
    logger=>qs_get_error_logger(error)

    call qs_init_test_pao_env(pao_env,glob_env,n_atom,error=error)
    call qs_get(pao_env, angles_att=g_angles,error=error)
    call qs_get(g_angles,angles=angles)
    call random_number(angles)
    angles=4.0_wp*pi*angles
    call qs_init(l_angles,glob_angles=g_angles, atom_nr=atom_nr)
    QSPostcondition(qs_valid(l_angles),qs_failure_level,routineN,error,failure)
    full_bas_size=qs_get_full_bas_size(l_angles)
    allocate(m1(full_bas_size,full_bas_size),m2(full_bas_size,full_bas_size),&
         m3(full_bas_size,full_bas_size),stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (qs_debug.and.debug_this_module) then
          call qs_log(logger,qs_note_level,routineN, &
               'test u_du 1: l_compose_u, U^* U = 1')
       end if
       m1=0.0_wp
       forall (i=1:size(m1,1)) m1(i,i)=1.0_wp
       call qs_pao_l_compose_u(m1,l_angles,error=error)
       m3=m1
       call dgemm('T','N',full_bas_size,full_bas_size,full_bas_size,&
            1.0_wp,m1,full_bas_size,m3,full_bas_size,&
            0.0_wp,m2,full_bas_size)
       forall (i=1:full_bas_size) m2(i,i)=m2(i,i)-1.0_wp
       err=qs_sup_norm_matrix(m2)
       if (err>2*full_bas_size*epsilon(0.0_wp)) then
          failure=.true.
          call qs_log(logger,qs_failure_level,routineN, &
               'test u_du 1: l_compose_u, U^* U = 1, failed')
          if (qs_debug .and. debug_this_module) then
             call qs_output_2d_r_array(logger, 'zero_matrix', routineN,3,m2,'zero matrix l_compose_u, U^* U - 1')
             call qs_output_1d_r_array(logger, 'angles', routineN,3,&
                  angles(1:(full_bas_size*qs_get_min_bas_size(l_angles))),&
                  'angles')
          end if
       else if (qs_debug .and. debug_this_module) then
          call qs_log(logger,qs_warning_level,routineN,&
               'test u_du 1: SUCCESS')
       end if
    end if
    deallocate(m1,m2,stat=stat) ! just to be sure...
    QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
    call qs_dealloc_ref(l_angles,error=error)
    call qs_dealloc_ref(pao_env,error=error)
    res=failure
  end function qs_pao_test_u_du
!!***
  !***************************************************************************

end module qs_pao_plain_rot_tests


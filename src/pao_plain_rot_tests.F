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
MODULE pao_plain_rot_tests
  USE cp_array_utils,                  ONLY: cp_2d_r_output
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sup_norm_matrix
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pao_plain_rot,                   ONLY: cp_left_compose_p_rot,&
                                             cp_left_compose_p_rot_deriv,&
                                             cp_right_compose_p_rot,&
                                             cp_right_compose_p_rot_deriv
  USE pao_utils,                       ONLY: cp_get_random_integer

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.

  PUBLIC :: pao_test_plain_rot, pao_test_u_du
!!***
!******************************************************************************
CONTAINS

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
  SUBROUTINE pao_test_plain_rot(g_env, failed, error)
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: g_env
    LOGICAL, INTENT(inout)                   :: failed
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_plain_rot_tests:cp_test_plain_rot'
    INTEGER, PARAMETER                       :: n = 10

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: angle, err, m1(n,n), m2(n,n), &
                                                tmpV1(n)
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(n>1,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       logger=>cp_error_get_logger(error)

       i=cp_get_random_integer(from=1,to=n)
       j=cp_get_random_integer(from=1,to=n)
       IF (i==j) THEN
          IF (i>1) THEN
             j=i-1
          ELSE
             j=i+1
          END IF
       END IF
       CALL RANDOM_NUMBER(m1)
       m2=m1

       ! t1r
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test1r: rotation with pi/2 swaps columns and change sign of the second')
       END IF
       tmpV1=m2(:,i)
       m2(:,i)=m2(:,j)
       m2(:,j)=-tmpV1
       CALL cp_right_compose_p_rot(m1,pi/2,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       IF (err>10*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test1r (right rot pi/2) with i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test1r: SUCCESS')
       END IF

       ! t2r
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test2r: rotation with -angle restores the matrix')
       END IF
       m2=m1
       CALL RANDOM_NUMBER(angle)
       angle=4.0_dp*pi*angle
       CALL cp_right_compose_p_rot(m1,angle,i,j)
       CALL cp_right_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test2r (right -angle = rot^-1) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test2r: SUCCESS')
       END IF

       !t1l
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test1l: rotation with pi/2 swaps lines and change sign of the first')
       END IF
       tmpV1=m2(i,:)
       m2(i,:)=-m2(j,:)
       m2(j,:)=tmpV1
       CALL cp_left_compose_p_rot(m1,pi/2,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test1l (right rot pi/2) with i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test1l: SUCCESS')
       END IF

       ! t2l
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test2l: rotation with -angle restores the matrix')
       END IF
       m2=m1
       CALL RANDOM_NUMBER(angle)
       angle=4.0_dp*pi*angle
       CALL cp_left_compose_p_rot(m1,angle,i,j)
       CALL cp_left_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test2l (right -angle = rot^-1) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test2l: SUCCESS')
       END IF

       ! t3 
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test3: the unit matrix is the same in all orthonormal basis')
       END IF
       m1=0
       FORALL (i=1:n) m1(i,i)=1
       m2=m1
       CALL RANDOM_NUMBER(angle)
       angle=4.0_dp*pi*angle
       CALL cp_left_compose_p_rot(m1,angle,i,j)
       CALL cp_right_compose_p_rot(m1,-angle,i,j)
       err=cp_sup_norm_matrix(m2-m1)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test3 (unity matrix same in all ONB) with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_2d_r_output(logger, 'identity_matrix', routineP,1,m1,'identity matrix')
          END IF
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test3: SUCCESS')
       END IF

       ! t4
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test4: d/dt (J(t,i,j)J(-t,i,j)M)=0')
       END IF
       ! test m1 = identity
       !m1=0
       !forall (i=1:n) m1(i,i)=1
       CALL RANDOM_NUMBER(m1)
       m2=m1
       CALL RANDOM_NUMBER(angle)
       angle=4.0_dp*pi*angle
       CALL cp_left_compose_p_rot(m1,-angle,i,j)
       !call cp_2d_r_output(logger, 'm1_1', routineP,1,m1,'m1.1')
       CALL cp_left_compose_p_rot_deriv(m1,angle,i,j)
       !call cp_2d_r_output(logger, 'm1_2', routineP,1,m1,'m1.2')
       CALL cp_left_compose_p_rot_deriv(m2,-angle,i,j)
       !call cp_2d_r_output(logger, 'm2_1', routineP,1,m2,'m2.1')
       CALL cp_left_compose_p_rot(m2,angle,i,j)
       !call cp_2d_r_output(logger, 'm2_2', routineP,1,m2,'m2.2')
       err=cp_sup_norm_matrix(m1-m2)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test4 d/dt (J(t,i,j)J(-t,i,j)M)=0 with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_2d_r_output(logger, 'zero_matrix', routineP,1,m1-m2,'zero matrix 1')
          END IF
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test4: SUCCESS')
       END IF

       ! t5
       IF (cp_debug.and.debug_this_module) THEN
          CALL cp_log(logger,cp_note_level,routineP, &
               'test5: d/dt (MJ(t,i,j)J(-t,i,j))=0')
       END IF
       CALL RANDOM_NUMBER(m1)
       m2=m1
       CALL RANDOM_NUMBER(angle)
       angle=4.0_dp*pi*angle
       CALL cp_right_compose_p_rot(m1,-angle,i,j)
       CALL cp_right_compose_p_rot_deriv(m1,angle,i,j)
       CALL cp_right_compose_p_rot_deriv(m2,-angle,i,j)
       CALL cp_right_compose_p_rot(m2,angle,i,j)
       err=cp_sup_norm_matrix(m1-m2)
       IF (err>2*n*EPSILON(0.0_dp)) THEN
          failure=.TRUE.
          CALL cp_log(logger,cp_failure_level,routineP, &
               'test5 d/dt (M J(t,i,j)J(-t,i,j))=0 with angle='//cp_to_string(angle)&
               //' i='//cp_to_string(i)// &
               ' j='//cp_to_string(j)//' failed')
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_2d_r_output(logger, 'zero_matrix', routineP,2,m1-m2,'zero matrix 2')
          END IF
       ELSE IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'test5: SUCCESS')
       END IF

       ! u_du tests
       IF (.not.pao_test_u_du(g_env,error=error)) THEN
          failure=.TRUE.
       END IF

       IF (cp_debug .AND. debug_this_module .AND. .not.failure) THEN
          CALL cp_log(logger,cp_warning_level,routineP,&
               'cp_test_plain_rot: all tests succeded, SUCCESS')
       END IF
    END IF
  END SUBROUTINE pao_test_plain_rot
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
  FUNCTION pao_test_u_du(glob_env,error) RESULT(res)
    TYPE(global_environment_type), TARGET    :: glob_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_plain_rot_tests:pao_test_u_du'
    INTEGER, PARAMETER                       :: atom_nr = 1, n_atom = 3

    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger=>cp_error_get_logger(error)

!!FM    call cp_init_test_pao_env(pao_env,glob_env,n_atom,error=error)
!!FM    call cp_get(pao_env, angles_att=g_angles,error=error)
!!FM    call cp_get(g_angles,angles=angles)
!!FM    call random_number(angles)
!!FM    angles=4.0_dp*pi*angles
!!FM    call cp_init(l_angles,glob_angles=g_angles, atom_nr=atom_nr)
!!FM    CPPostcondition(cp_valid(l_angles),cp_failure_level,routineP,error,failure)
!!FM    full_bas_size=cp_get_full_bas_size(l_angles)
!!FM    allocate(m1(full_bas_size,full_bas_size),m2(full_bas_size,full_bas_size),&
!!FM         m3(full_bas_size,full_bas_size),stat=stat)
!!FM    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!!FM    if (.not. failure) then
!!FM       if (cp_debug.and.debug_this_module) then
!!FM          call cp_log(logger,cp_note_level,routineP, &
!!FM               'test u_du 1: l_compose_u, U^* U = 1')
!!FM       end if
!!FM       m1=0.0_dp
!!FM       forall (i=1:size(m1,1)) m1(i,i)=1.0_dp
!!FM       call pao_l_compose_u(m1,l_angles,error=error)
!!FM       m3=m1
!!FM       call dgemm('T','N',full_bas_size,full_bas_size,full_bas_size,&
!!FM            1.0_dp,m1,full_bas_size,m3,full_bas_size,&
!!FM            0.0_dp,m2,full_bas_size)
!!FM       forall (i=1:full_bas_size) m2(i,i)=m2(i,i)-1.0_dp
!!FM       err=cp_sup_norm_matrix(m2)
!!FM       if (err>2*full_bas_size*epsilon(0.0_dp)) then
!!FM          failure=.true.
!!FM          call cp_log(logger,cp_failure_level,routineP, &
!!FM               'test u_du 1: l_compose_u, U^* U = 1, failed')
!!FM          if (cp_debug .and. debug_this_module) then
!!FM             call cp_2d_r_output(logger, 'zero_matrix', routineP,3,m2,'zero matrix l_compose_u, U^* U - 1')
!!FM             call cp_1d_r_output(logger, 'angles', routineP,3,&
!!FM                  angles(1:(full_bas_size*cp_get_min_bas_size(l_angles))),&
!!FM                  'angles')
!!FM          end if
!!FM       else if (cp_debug .and. debug_this_module) then
!!FM          call cp_log(logger,cp_warning_level,routineP,&
!!FM               'test u_du 1: SUCCESS')
!!FM       end if
!!FM    end if
!!FM    deallocate(m1,m2,stat=stat) ! just to be sure...
!!FM    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!FM    call cp_dealloc_ref(l_angles,error=error)
!!FM    call cp_dealloc_ref(pao_env,error=error)
    res=failure
  END FUNCTION pao_test_u_du
  !***************************************************************************

END MODULE pao_plain_rot_tests


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
MODULE pao_unitary_t
  USE cp_array_utils,                  ONLY: cp_2d_r_output
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sup_norm_matrix
  USE kinds,                           ONLY: dp
  USE pao_l_angles_methods,            ONLY: cp_get_angle,&
                                             cp_get_non_p_index,&
                                             cp_get_p_index,&
                                             cp_get_size_excl_bas,&
                                             cp_get_size_min_bas,&
                                             cp_init,&
                                             pao_l_ang_get,&
                                             pao_local_angles_type
  USE pao_plain_rot,                   ONLY: cp_left_compose_p_rot,&
                                             cp_right_compose_p_rot,&
                                             cp_right_compose_p_rot_deriv

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  PUBLIC :: pao_calc_u_du_type

  PUBLIC :: cp_init, cp_dealloc_ref, cp_get, cp_next

  PUBLIC :: pao_calc_u_du_init, pao_calc_u_du_dealloc_ref,&
       pao_calc_u_du_get, pao_calc_u_du_next

! == core methods ==
!! initializes the object
  INTERFACE cp_init
     MODULE PROCEDURE pao_calc_u_du_init
  END INTERFACE

!! deallocates the memory used by the first argument
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE pao_calc_u_du_dealloc_ref
  END INTERFACE

!! returns various attributes of the first argument 
  INTERFACE cp_get
     MODULE PROCEDURE pao_calc_u_du_get
  END INTERFACE

!! iterator interface: returns true if there is a next element
  INTERFACE cp_next
     MODULE PROCEDURE pao_calc_u_du_next
  END INTERFACE

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
  TYPE pao_calc_u_du_type
     PRIVATE
     LOGICAL :: initialized
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: L,R,tmp_m
     INTEGER :: min_bas, excl_bas, status
     TYPE(pao_local_angles_type), POINTER :: local_angles
  END TYPE pao_calc_u_du_type
!!***
!****************************************************************************

CONTAINS
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
!!       Real(KIND = dp), Dimension(:,:), Optional, Intent (IN):: u
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
  SUBROUTINE pao_calc_u_du_init(c_env,local_angles, min_bas, excl_bas,&
       u, error) !to do (correct full to excl)
    TYPE(pao_calc_u_du_type), INTENT(out)    :: c_env
    TYPE(pao_local_angles_type), &
      INTENT(inout), TARGET                  :: local_angles
    INTEGER, INTENT(in), OPTIONAL            :: min_bas, excl_bas
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL                   :: u
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_plain_rot:pao_calc_u_du_init'

    INTEGER                                  :: i, i_excl_bas, i_min_bas, &
                                                n_excl_bas, n_full_bas, &
                                                n_min_bas, stat
    INTEGER, DIMENSION(:), POINTER           :: non_proj_indexes, proj_indexes
    LOGICAL                                  :: failure

    failure=.FALSE.
    c_env%initialized=.TRUE.
    c_env%status=0
    c_env%local_angles => local_angles
    NULLIFY(c_env%L, c_env%R, c_env%tmp_m)
    IF (PRESENT(min_bas)) THEN
       c_env%min_bas=min_bas
    ELSE
       c_env%min_bas=1
    END IF
    IF (PRESENT(excl_bas)) THEN
       c_env%excl_bas=excl_bas
    ELSE
       c_env%excl_bas=1
    END IF
    CALL pao_l_ang_get(local_angles, size_excl_bas=n_excl_bas, &
         size_min_bas=n_min_bas,&
         size_full_bas=n_full_bas, proj_indexes=proj_indexes,&
         non_proj_indexes=non_proj_indexes, error=error)
    CPPrecondition(n_excl_bas>=0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_env%L(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_env%R(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(c_env%tmp_m(n_full_bas,n_full_bas),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (PRESENT(u)) THEN
          IF (PRESENT(min_bas)) THEN
             CPPrecondition(min_bas==1,cp_failure_level,routineP,error,failure)
          END IF
          IF (PRESENT(excl_bas)) THEN
             CPPrecondition(excl_bas==1,cp_failure_level,routineP,error,failure)
          END IF
          c_env%L=0
          FORALL (i=1:n_excl_bas) c_env%L(i,i)=1
          c_env%R=u
          CALL cp_left_compose_p_rot(matrix=c_env%R,&
               angle=-cp_get_angle(local_angles,min_bas=n_min_bas,&
               excl_bas=n_excl_bas,error=error),&
               i=cp_get_p_index(local_angles,min_bas=n_min_bas,error=error),&
               j=cp_get_non_p_index(local_angles,excl_bas=n_excl_bas,error=error),&
               error=error)
       ELSE
          c_env%L=0
          c_env%R=0
          FORALL (i=1:n_full_bas)
             c_env%L(i,i)=1
             c_env%R(i,i)=1
          END FORALL
! def R
          DO i_excl_bas=n_excl_bas,c_env%excl_bas+1,-1
             DO i_min_bas=n_min_bas,1,-1
                CALL cp_left_compose_p_rot(matrix=c_env%R,&
                     angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=proj_indexes(i_min_bas),&
                     j=non_proj_indexes(i_excl_bas),&
                     error=error)
             END DO
          END DO
          i_excl_bas=c_env%excl_bas
          DO i_min_bas=n_min_bas,c_env%min_bas+1,-1
             CALL cp_left_compose_p_rot(matrix=c_env%R,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          END DO
! def L
          DO i_excl_bas=1,c_env%excl_bas-1
             DO i_min_bas=1,n_min_bas
                CALL cp_right_compose_p_rot(matrix=c_env%L,&
                     angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                     excl_bas=i_excl_bas,error=error),&
                     i=proj_indexes(i_min_bas),&
                     j=non_proj_indexes(i_excl_bas),&
                     error=error)
             END DO
          END DO
          i_excl_bas=c_env%excl_bas
          DO i_min_bas=1,c_env%min_bas-1
             CALL cp_right_compose_p_rot(matrix=c_env%L,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          END DO
       END IF
       c_env%min_bas=c_env%min_bas-1 ! to start up goes one before
    END IF failureIf
    IF (failure) THEN
       c_env%status=-1
    END IF
  END SUBROUTINE pao_calc_u_du_init
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
  FUNCTION pao_calc_u_du_next(c_env,error) RESULT(res)
    TYPE(pao_calc_u_du_type), INTENT(inout)  :: c_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_plain_rot:pao_calc_u_du_next'

    INTEGER                                  :: i, old_excl_bas, old_min_bas
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: err

    failure=.FALSE.

    CPPrecondition(c_env%initialized,cp_failure_level,routineP,error,failure)
    IF (debug_this_module .AND. cp_debug) THEN
       CPPrecondition(ASSOCIATED(c_env%R),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(c_env%L),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       old_min_bas=c_env%min_bas
       old_excl_bas=c_env%excl_bas
       IF (c_env%min_bas>=cp_get_size_min_bas(c_env%local_angles)) THEN
          IF (c_env%excl_bas>=cp_get_size_excl_bas(c_env%local_angles)) THEN
             c_env%status=2
             res=.FALSE.
             RETURN
          ELSE IF (cp_get_size_excl_bas(c_env%local_angles)==0) THEN
             c_env%status=2
             res=.FALSE.
             RETURN
          ELSE
             c_env%min_bas=1
             c_env%excl_bas=c_env%excl_bas+1
          END IF
       ELSE
          c_env%min_bas=c_env%min_bas+1
       END IF
       SELECT CASE (c_env%status)
       CASE (-2) ! status after dealloc
          res=.FALSE.
       CASE (-1) ! status after an error
          res=.FALSE.
       CASE (0) ! initial status
          res=.TRUE.
          c_env%status=1
       CASE (1)
! normal case
          res=.TRUE.
          CALL cp_right_compose_p_rot(matrix=c_env%L,&
               angle=cp_get_angle(c_env%local_angles,min_bas=old_min_bas,&
               excl_bas=old_excl_bas), &
               i=cp_get_p_index(c_env%local_angles,old_min_bas,error=error), &
               j=cp_get_non_p_index(c_env%local_angles,old_excl_bas,&
               error=error), error=error)
          CALL cp_left_compose_p_rot(matrix=c_env%R,&
               angle= -cp_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
               excl_bas=c_env%excl_bas), &
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error))
       CASE (2)
          res=.FALSE.
          IF (debug_this_module .AND. cp_debug) THEN
             DO i=1,SIZE(c_env%R)
                c_env%R(i,i)=c_env%R(i,i)-1
             END DO
             err=cp_sup_norm_matrix(c_env%R)
             IF (err>2*SIZE(c_env%R)*EPSILON(0.0_dp)) THEN
                CALL cp_error_message(cp_warning_level,routineP,&
                     message="accumulation of the error bigger than tolerance:"//&
                     cp_to_string(err)//" in "//&
                     CPSourceFileRef,&
                     error=error)
             END IF
          END IF
       CASE default
          res=.FALSE.
          CALL cp_error_message(cp_warning_level,routineP,&
               message="unknown status "//cp_to_string(c_env%status)//&
               " in "// &
               CPSourceFileRef,&
               error=error)
          c_env%status=-1
       END SELECT
    END IF
  END FUNCTION pao_calc_u_du_next
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
!!       Real(KIND = dp), Optional, Intent (OUT), Dimension(:,:):: u, du
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
  SUBROUTINE pao_calc_u_du_get(c_env,u,du,min_bas,excl_bas,local_angles,&
       finished,error)
    TYPE(pao_calc_u_du_type), INTENT(inout)  :: c_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL                  :: u, du
    INTEGER, INTENT(out), OPTIONAL           :: min_bas, excl_bas
    TYPE(pao_local_angles_type), OPTIONAL, &
      POINTER                                :: local_angles
    LOGICAL, INTENT(out), OPTIONAL           :: finished
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_plain_rot:pao_calc_u_du_get'

    INTEGER                                  :: n_excl_bas, n_full_bas
    INTEGER, SAVE                            :: iter = 0
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: angle
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    iter = iter+1

    CPPrecondition(c_env%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(c_env%status==1,cp_failure_level,routineP,error,failure)
    IF (debug_this_module .AND. cp_debug) THEN
       CPPrecondition(ASSOCIATED(c_env%local_angles),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(c_env%R),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(c_env%L),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(c_env%tmp_m),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       CALL pao_l_ang_get(c_env%local_angles,size_excl_bas=n_excl_bas,&
            size_full_bas=n_full_bas)
       angle=cp_get_angle(c_env%local_angles,min_bas=c_env%min_bas,&
            excl_bas=c_env%excl_bas,error=error)
       IF (PRESENT(u)) THEN
          CPPrecondition(SIZE(u,1)==n_full_bas,cp_failure_level,routineP,error,failure)
          CPPrecondition(SIZE(u,2)==n_full_bas,cp_failure_level,routineP,error,failure)
          c_env%tmp_m=c_env%L
          CALL cp_right_compose_p_rot(matrix=c_env%tmp_m,&
               angle=angle,&
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error),error=error)
          CALL dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_dp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_dp,u,n_full_bas)
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_2d_r_output(logger=logger,&
                  outputName='u',fromWhere=routineP, iter=iter,&
                  array=u,comment=" min_bas="//cp_to_string(c_env%min_bas)//&
                  " excl_bas="//cp_to_string(c_env%excl_bas),error=error)
          END IF
       END IF
       IF (PRESENT(du)) THEN
          CPPrecondition(SIZE(du,1)==n_full_bas,cp_failure_level,routineP,error,failure)
          CPPrecondition(SIZE(du,2)==n_full_bas,cp_failure_level,routineP,error,failure)
          c_env%tmp_m=c_env%L
          CALL cp_right_compose_p_rot_deriv(matrix=c_env%tmp_m,&
               angle=angle,&
               i=cp_get_p_index(c_env%local_angles,c_env%min_bas,error=error),&
               j=cp_get_non_p_index(c_env%local_angles,c_env%excl_bas,&
               error=error),error=error)
          CALL dgemm('N','N',n_full_bas,n_full_bas,n_full_bas,1.0_dp,&
               c_env%tmp_m,n_full_bas,c_env%R,n_full_bas,0.0_dp,du,n_full_bas)
          IF (cp_debug .AND. debug_this_module) THEN
!!FM             nullify(ang)
!!FM             call pao_l_ang_get(c_env%local_angles,angles=ang)
!!FM             call cp_1d_r_output(logger=logger,&
!!FM               outputName='udu_angles',fromWhere=routineP, iter=iter,&
!!FM               array=ang,error=error)
             CALL cp_2d_r_output(logger=logger,&
               outputName='udu_L',fromWhere=routineP, iter=iter,&
               array=c_env%L,comment=" min_bas="//&
               cp_to_string(c_env%min_bas)//&
                  " excl_bas="//cp_to_string(c_env%excl_bas),error=error)
             CALL cp_2d_r_output(logger=logger,&
               outputName='udu_R',fromWhere=routineP, iter=iter,&
               array=c_env%R,comment=" min_bas="//&
               cp_to_string(c_env%min_bas)//&
                  " excl_bas="//cp_to_string(c_env%excl_bas),error=error)
             CALL cp_2d_r_output(logger=logger,&
               outputName='du',fromWhere=routineP, iter=iter,&
               array=du,comment=" min_bas="//cp_to_string(c_env%min_bas)//&
                  " excl_bas="//cp_to_string(c_env%excl_bas),error=error)
          END IF
       END IF
       IF (PRESENT(min_bas)) &
            min_bas=c_env%min_bas
       IF (PRESENT(excl_bas)) &
            excl_bas=c_env%excl_bas
       IF (PRESENT(local_angles)) &
            local_angles => c_env%local_angles
       IF (PRESENT(finished)) THEN
          finished=c_env%status /= 1 .AND. c_env%status /= 0
! catch also status=0 and past end?
       END IF
    END IF
  END SUBROUTINE pao_calc_u_du_get
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
  SUBROUTINE pao_calc_u_du_dealloc_ref(c_env, error)
    TYPE(pao_calc_u_du_type), INTENT(inout)  :: c_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_plain_rot:pao_calc_u_du_dealloc_ref'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPAssert(c_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(c_env%R)) THEN
          DEALLOCATE(c_env%R, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(c_env%L)) THEN
          DEALLOCATE(c_env%L, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(c_env%tmp_m)) THEN
          DEALLOCATE(c_env%tmp_m, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       NULLIFY(c_env%local_angles)
       c_env%status=-2
    END IF
  END SUBROUTINE pao_calc_u_du_dealloc_ref
!!***
!***************************************************************************


END MODULE pao_unitary_t

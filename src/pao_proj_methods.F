!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_proj_methods [1.0] *
!!
!!   NAME
!!     pao_proj_methods
!!
!!   FUNCTION
!!     projection and injection between min bas and full bas
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     2.2002 created [fawzi]
!!     4.2002 uses projection object and qs_env [fawzi]
!!     6.2002 added atomic orthogonalization
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_proj_methods
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_sp_create
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_dealloc_ref,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_propagate_error,&
                                             cp_error_type,&
                                             cp_precondition_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_dealloc_ref,&
                                             cp_get,&
                                             cp_init,&
                                             cp_matrix_block_iterator,&
                                             cp_next
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pao_l_angles_methods,            ONLY: cp_dealloc_ref,&
                                             cp_get,&
                                             cp_init,&
                                             cp_set,&
                                             pao_local_angles_type
  USE pao_qs_env_methods,              ONLY: pao_qs_env_get
  USE pao_types,                       ONLY: pao_env_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_proj_methods'

  PUBLIC :: cp_get

! projection - injection
  PUBLIC :: cp_m2f, cp_f2m, cp_m2f_with_B_M_pi

  PUBLIC :: pao_g_ang_get_min_m, sm_matrix_transf_to_m,&
       sm_matrix_transf_to_f, pao_m_m2f_with_b_m_pi,&
       sm_matrix_transf_p_to_m

  INTERFACE cp_get
     MODULE PROCEDURE pao_g_ang_get_min_m
  END INTERFACE

  ! transfer from full basis to min basis
  INTERFACE cp_f2m
     MODULE PROCEDURE sm_matrix_transf_to_m
  END INTERFACE

  !transfer from min basis to full basis
  INTERFACE cp_m2f
     MODULE PROCEDURE sm_matrix_transf_to_f
  END INTERFACE
  
  ! tranfser to full basis with B matrix pi
  INTERFACE cp_m2f_with_B_M_pi
     MODULE PROCEDURE pao_m_m2f_with_b_m_pi
  END INTERFACE
!!***
!****************************************************************************

CONTAINS

!!****f* pao_types/pao_m_m2f_with_b_m_pi [1.0] *
!!
!!   NAME
!!     pao_m_m2f_with_b_m_pi
!!
!!   SYNOPSIS
!!     Subroutine pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_m_m2f_with_b_m_pi
!!
!!   FUNCTION
!!     Transform the given symmetric matrix to the full basis
!!     with B= N U i and pi=i^*.
!!     Puts N U i source_m pi  into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the min basis
!!     target_m: the result matrix in the full basis
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
  SUBROUTINE pao_m_m2f_with_b_m_pi(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    TYPE(real_matrix_type), POINTER          :: source_m, target_m
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_types:pao_m_m2f_with_b_m_pi'

    INTEGER                                  :: bcol, brow, i, j, max_bsize, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: l_index
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: t1_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val, u_row
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(pao_local_angles_type)              :: l_angles_col, l_angles_row

! make target?
! make target?

    CALL cp_error_init(iError,template_error=error)
    NULLIFY(block_val,l_index,u_row,pao_env)
    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL pao_qs_env_get(qs_env,global_env=global_env, &
            pao_max_full_bas=max_bsize)
       ALLOCATE(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_init(iterator,source_m,error=iError)
       CALL cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       IF (.not.failure) THEN
          blocksDo: DO
             IF ((.not.cp_next(iterator,error=error)).or.failure) EXIT
             CALL cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
             ! get u
             CALL cp_set(l_angles_col, atom_nr=bcol,error=error)
             CALL cp_get(l_angles_col, proj_indexes=l_index)
             CALL cp_set(l_angles_row, atom_nr=brow,error=error)
             CALL cp_get(l_angles_row,NUi_injection=u_row)

             ! t1_m = block_val pi
             t1_m=0.0_dp
             CPPrecondition(ASSOCIATED(l_index),cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                CALL cp_assert(SIZE(l_index)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(block_val,1)==SIZE(u_row,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             IF (.not.failure) THEN ! t1_m = block_val pi
                DO j=1,SIZE(block_val,2)
                   DO i=1,SIZE(block_val,1)
                      t1_m(i,l_index(j))=block_val(i,j)
                   END DO
                END DO
             END IF

             ! get full basis block
             CALL get_block_node(target_m, block_col=bcol,&
                  block_row=brow, BLOCK=block_val)
             IF (.not.ASSOCIATED(block_val)) THEN
                CALL add_block_node(target_m, block_col=bcol,&
                     block_row=brow, BLOCK=block_val)
                CPPostcondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
             END IF
             IF (.not.failure) THEN
                CALL cp_assert(SIZE(block_val,1)==SIZE(u_row,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             ! block_val = NUi t1_m
             IF (.not.failure) THEN
                CALL dgemm('N','N',SIZE(u_row,1),SIZE(block_val,2),&
                     SIZE(u_row,2),1.0_dp,&
                     u_row,SIZE(u_row,1),&
                     t1_m,SIZE(t1_m,1),&
                     0.0_dp,block_val,SIZE(block_val,1))
             END IF

             IF (brow /= bcol) THEN ! do the transposed block
                CALL cp_get(iterator,block_row=bcol,block_col=brow,&
                     block_val=block_val)
                ! get u
                CALL cp_set(l_angles_col, atom_nr=bcol,error=error)
                CALL cp_get(l_angles_col, proj_indexes=l_index)
                CALL cp_set(l_angles_row, atom_nr=brow,error=error)
                CALL cp_get(l_angles_row,NUi_injection=u_row)

                ! t1_m = block_val pi
                t1_m=0.0_dp
                CPPrecondition(ASSOCIATED(l_index),cp_failure_level,routineP,error,failure)
                IF (.not.failure) THEN
                   CALL cp_assert(SIZE(l_index)==SIZE(block_val,1),&
                        cp_failure_level,cp_precondition_failed,routineP,&
                        "size mishmash in "//&
                        CPSourceFileRef,&
                        error,failure)
                   CALL cp_assert(SIZE(block_val,2)==SIZE(u_row,2),&
                        cp_failure_level,cp_precondition_failed,routineP,&
                        "size mishmash in "//&
                        CPSourceFileRef,&
                        error,failure)
                END IF
                IF (.not.failure) THEN ! t1_m = block_val pi
                   DO i=1,SIZE(block_val,2)
                      DO j=1,SIZE(block_val,1)
                         t1_m(i,l_index(j))=block_val(j,i)
                      END DO
                   END DO
                END IF

                ! get full basis block
                CALL get_block_node(target_m, block_col=bcol,&
                     block_row=brow, BLOCK=block_val)
                IF (.not.ASSOCIATED(block_val)) THEN
                   CALL add_block_node(target_m, block_col=bcol,&
                        block_row=brow, BLOCK=block_val)
                   CPPostcondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
                END IF
                IF (.not.failure) THEN
                   CALL cp_assert(SIZE(block_val,1)==SIZE(u_row,1),&
                        cp_failure_level,cp_precondition_failed,routineP,&
                        "size mishmash in "//&
                        CPSourceFileRef,&
                        error,failure)
                END IF
                ! block_val = NUi t1_m
                IF (.not.failure) THEN
                   CALL dgemm('N','N',SIZE(u_row,1),SIZE(block_val,2),&
                        SIZE(u_row,2),1.0_dp,&
                        u_row,SIZE(u_row,1),&
                        t1_m,SIZE(t1_m,1),&
                        0.0_dp,block_val,SIZE(block_val,1))
                END IF
             END IF
             
          END DO blocksDo
       END IF
       CALL cp_dealloc_ref(iterator,error=error)
       CALL cp_dealloc_ref(l_angles_col,error=error)
       CALL cp_dealloc_ref(l_angles_row,error=error)
       DEALLOCATE(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE pao_m_m2f_with_b_m_pi
!***************************************************************************


!!****f* pao_types/sm_matrix_transf_to_f [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_to_f
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_to_f(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_to_f
!!
!!   FUNCTION
!!     Transform the given matrix to the full basis.
!!     Puts N^(-1) U i source_m pi U^* N^(-1*)  into the target_m
!!
!!   NOTES
!!     could be optimized by caching N^(-1) (i.e. atomic S^(1/2))
!!
!!   INPUTS
!!     source_m: matrix in the min basis
!!     target_m: the result matrix in the full basis
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
  SUBROUTINE sm_matrix_transf_to_f(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    TYPE(real_matrix_type), POINTER          :: source_m, target_m
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_matrix_transf_to_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcol, brow, max_bsize, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: t1_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val, u_col, u_row
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(pao_local_angles_type)              :: l_angles_col, l_angles_row

! make target?
! make target ?

    CALL cp_error_init(iError,template_error=error)
    NULLIFY(block_val,u_col, u_row,pao_env)
    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CALL cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       ALLOCATE(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_init(iterator,source_m,error=iError)
       CALL cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       IF (.not.failure) THEN
          blocksDo: DO
             IF ((.not.cp_next(iterator,error=error)).or.failure) EXIT
             CALL cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             CALL cp_set(l_angles_col, atom_nr=bcol,error=error)
             CALL cp_get(l_angles_col,NUi_injection=u_col)
             CALL cp_set(l_angles_row, atom_nr=brow,error=error)
             CALL cp_get(l_angles_row,NUi_injection=u_row)

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_col,2)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_row,2)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF

! t1_m = block_val (N U i)^*
             IF (.not.failure) THEN
                CALL dgemm('N','T',SIZE(block_val,1),SIZE(u_col,1),&
                     SIZE(u_col,2),1.0_dp,&
                     block_val,SIZE(block_val,1),&
                     u_col,SIZE(u_col,1),&
                     0.0_dp,t1_m,SIZE(t1_m,1))

                ! get full basis block
                CALL get_block_node(target_m, block_col=bcol,&
                     block_row=brow, BLOCK=block_val)
                IF (.not.ASSOCIATED(block_val)) THEN
                   CALL add_block_node(target_m, block_col=bcol,&
                        block_row=brow, BLOCK=block_val)
                   CPPostcondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
                END IF
             END IF

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_col,1)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_row,1)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF

! block_val = NUi t1_m
             IF (.not.failure) THEN
                CALL dgemm('N','N',SIZE(u_row,1),SIZE(u_col,1),&
                     SIZE(u_row,2),1.0_dp,&
                     u_row,SIZE(u_row,1),&
                     t1_m,SIZE(t1_m,1),&
                     0.0_dp,block_val,SIZE(block_val,1))
             END IF

          END DO blocksDo
       END IF
       CALL cp_dealloc_ref(iterator,error=error)
       CALL cp_dealloc_ref(l_angles_col,error=error)
       CALL cp_dealloc_ref(l_angles_row,error=error)
       DEALLOCATE(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE sm_matrix_transf_to_f
!***************************************************************************

!!****f* pao_types/sm_matrix_transf_to_m [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_to_m
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_to_m(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* N source_m N U i into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the full basis
!!     target_m: the result matrix in the small basis
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
  SUBROUTINE sm_matrix_transf_to_m(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    TYPE(real_matrix_type), POINTER          :: source_m, target_m
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_matrix_transf_to_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcol, brow, max_bsize, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: t1_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val, u_col, u_row
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(pao_local_angles_type)              :: l_angles_col, l_angles_row

! make target?
! make target?

    CALL cp_error_init(iError,template_error=error)
    NULLIFY(block_val,block_val, u_col, u_row,pao_env)
    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CALL cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       ALLOCATE(t1_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_init(iterator,source_m,error=iError)
       CALL cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       IF (.not.failure) THEN
          blocksDo: DO
             IF ((.not.cp_next(iterator,error=error)).or.failure) EXIT
             CALL cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             CALL cp_set(l_angles_col, atom_nr=bcol,error=error)
             CALL cp_get(l_angles_col,NUi_injection=u_col,error=error)
             CALL cp_set(l_angles_row, atom_nr=brow,error=error)
             CALL cp_get(l_angles_row,NUi_injection=u_row,error=error)

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_row,1)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_col,1)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF

! t1_m = source NUi
             CALL dgemm('N','N',SIZE(block_val,1),SIZE(u_col,2),&
                  SIZE(u_col,1),1.0_dp,&
                  block_val,SIZE(block_val,1),&
                  u_col,SIZE(u_col,1),&
                  0.0_dp,t1_m,SIZE(t1_m,1))
! get small basis block
             CALL get_block_node(target_m, block_col=bcol,&
                  block_row=brow, BLOCK=block_val)
             IF (.not.ASSOCIATED(block_val)) THEN
                CALL add_block_node(target_m, block_col=bcol,&
                     block_row=brow, BLOCK=block_val)
                CPPostcondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
             END IF

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_row,2)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_col,2)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             ! target_m = (NUi)^* t1_m
             IF (.not.failure) THEN 
                CALL dgemm('T','N',SIZE(u_row,2),SIZE(u_col,2),&
                     SIZE(u_row,1),1.0_dp,&
                     u_row,SIZE(u_row,1),&
                     t1_m,SIZE(t1_m,1),&
                     0.0_dp,block_val,SIZE(block_val,1))
             END IF

          END DO blocksDo
       END IF
       CALL cp_dealloc_ref(iterator,error=error)
       CALL cp_dealloc_ref(l_angles_col,error=error)
       CALL cp_dealloc_ref(l_angles_row,error=error)
       DEALLOCATE(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE sm_matrix_transf_to_m
!***************************************************************************

!!****f* pao_types/sm_matrix_transf_p_to_m [1.0] *
!!
!!   NAME
!!     sm_matrix_transf_p_to_m
!!
!!   SYNOPSIS
!!     Subroutine sm_matrix_transf_p_to_m(source_m, target_m, glob_angles,&
!!         qs_env, global_env, error)
!!       Type(real_matrix_type), Pointer:: source_m
!!       Type(real_matrix_type), Pointer:: target_m
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_matrix_transf_p_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* N^(-1) source_m N^(-1) U i into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the full basis
!!     target_m: the result matrix in the small basis
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
  SUBROUTINE sm_matrix_transf_p_to_m(source_m, target_m, glob_angles,&
       qs_env, global_env, error)
    TYPE(real_matrix_type), POINTER          :: source_m, target_m
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_matrix_transf_p_to_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcol, brow, max_bsize, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: t1_m, t2_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val, s_col, s_row, &
                                                u_col, u_row
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(pao_local_angles_type)              :: l_angles_col, l_angles_row
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s

! make target?
! make target?

    CALL cp_error_init(iError,template_error=error)
    NULLIFY(block_val,block_val, u_col, u_row,pao_env)
    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CALL cp_init(l_angles_col,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL cp_init(l_angles_row,atom_nr=1, &
            glob_angles=glob_angles, qs_env=qs_env, global_env=global_env,&
            error=error)
       CALL get_qs_env(qs_env,s=s)
       CALL pao_qs_env_get(qs_env,global_env=global_env,&
            pao_max_full_bas=max_bsize)
       ALLOCATE(t1_m(max_bsize,max_bsize),t2_m(max_bsize,max_bsize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_init(iterator,source_m,error=iError)
       CALL cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
       IF (.not.failure) THEN
          blocksDo: DO
             IF ((.not.cp_next(iterator,error=error)).or.failure) EXIT
             CALL cp_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             CALL cp_set(l_angles_col, atom_nr=bcol,error=error)
             CALL cp_get(l_angles_col,NUi_injection=u_col,error=error)
             CALL cp_set(l_angles_row, atom_nr=brow,error=error)
             CALL cp_get(l_angles_row,NUi_injection=u_row,error=error)

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_row,1)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_col,1)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF

             CALL get_block_node(s(1)%matrix,block_row=bcol,block_col=bcol,&
                  BLOCK=s_col)
             CALL get_block_node(s(1)%matrix,block_row=brow,block_col=brow,&
                  BLOCK=s_row)
             IF (.not.failure) THEN
                CALL cp_assert(SIZE(s_col,1)==SIZE(block_val,2).AND.&
                     SIZE(s_col,1)==SIZE(s_col,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(s_row,2)==SIZE(block_val,1).AND.&
                     SIZE(s_row,1)==SIZE(s_row,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             ! t1_m = block_val S
             IF (.not.failure) THEN
                CALL dgemm('N','N',SIZE(block_val,1),SIZE(block_val,2),&
                     SIZE(block_val,2),1.0_dp,&
                     block_val,SIZE(block_val,1),&
                     s_col,SIZE(s_col,1),&
                     0.0_dp,t1_m,SIZE(t1_m,1))
             END IF
             ! t2_m = S t1_m (could be partially faster if done later)
             IF (.not.failure) THEN
                CALL dgemm('N','N',SIZE(block_val,1),SIZE(block_val,2),&
                     SIZE(block_val,1),1.0_dp,&
                     s_row,SIZE(s_row,1),&
                     t1_m,SIZE(t1_m,1),&
                     0.0_dp,t2_m,SIZE(t2_m,1))
             END IF

! t1_m = t2_m NUi
             CALL dgemm('N','N',SIZE(block_val,1),SIZE(u_col,2),&
                  SIZE(u_col,1),1.0_dp,&
                  t2_m,SIZE(t2_m,1),&
                  u_col,SIZE(u_col,1),&
                  0.0_dp,t1_m,SIZE(t1_m,1))
! get small basis block
             CALL get_block_node(target_m, block_col=bcol,&
                  block_row=brow, BLOCK=block_val)
             IF (.not.ASSOCIATED(block_val)) THEN
                CALL add_block_node(target_m, block_col=bcol,&
                     block_row=brow, BLOCK=block_val)
                CPPostcondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
             END IF

             IF (.not.failure) THEN
                CALL cp_assert(SIZE(u_row,2)==SIZE(block_val,1),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_assert(SIZE(u_col,2)==SIZE(block_val,2),&
                     cp_failure_level,cp_precondition_failed,routineP,&
                     "size mishmash in "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             ! target_m = (NUi)^* t1_m
             IF (.not.failure) THEN 
                CALL dgemm('T','N',SIZE(u_row,2),SIZE(u_col,2),&
                     SIZE(u_row,1),1.0_dp,&
                     u_row,SIZE(u_row,1),&
                     t1_m,SIZE(t1_m,1),&
                     0.0_dp,block_val,SIZE(block_val,1))
             END IF

          END DO blocksDo
       END IF
       CALL cp_dealloc_ref(iterator,error=error)
       CALL cp_dealloc_ref(l_angles_col,error=error)
       CALL cp_dealloc_ref(l_angles_row,error=error)
       DEALLOCATE(t1_m,stat=stat) ! just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE sm_matrix_transf_p_to_m
!***************************************************************************

!!****f* pao_proj_methods/pao_g_ang_get_min_m [1.0] *
!!
!!   NAME
!!     pao_g_ang_get_min_m
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_get_min_m(glob_angles, qs_env, global_env,&
!!         min_overlap_m, min_hamiltonian_m, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_block_matrix_p_type), Dimension(:), Pointer, Optional::&
!!         min_hamiltonian_m
!!       Type(cp_block_matrix_p_type), Dimension(:), Pointer::&
!!         min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_get_min_m
!!
!!   FUNCTION
!!     returns the minimal matrices, calculating them if necessary
!!
!!   NOTES
!!     very ugly: to avoid circular dependencies this is not in the get
!!     of the pao glob angles. to do: fix it!!!
!!
!!   INPUTS
!!     glob_angles: 
!!     min_hamiltonian_m: the hamiltonian in the small basis. 
!!            If not yet calculated calculates it (allocating it if necessary)
!!     min_overlap_m: the overlap matrix in the small basis
!!            If not yet calculated calculates it (allocating it if necessary)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_get_min_m(glob_angles,qs_env,global_env,&
       min_overlap_m,  min_hamiltonian_m, error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), POINTER                  :: min_overlap_m
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: min_hamiltonian_m
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_g_ang_get_min_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_b_matrix_struct_type), POINTER   :: full_bas_m_struct, &
                                                min_bas_m_struct
    TYPE(cp_error_type)                      :: iError
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: full_h, full_s

    CALL cp_error_init(iError, template_error=error)
    failure=.FALSE.
    NULLIFY(full_h, full_s, pao_env)

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,k=full_h,s=full_s,pao_env=pao_env)
       CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       CALL pao_env_get(pao_env,min_bas_m_struct=min_bas_m_struct,&
            full_bas_diag_m_struct=full_bas_m_struct,error=error)
       !       if (present(min_overlap_m)) then
       IF (.not.glob_angles%m_s_m_valid) THEN
          IF (.not.ASSOCIATED(glob_angles%min_overlap_m)) THEN
             ALLOCATE(glob_angles%min_overlap_m(1),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL cp_sp_create(glob_angles%min_overlap_m(1)%matrix,&
                  matrix_struct=min_bas_m_struct,&
                  error=iError)
          END IF
          CALL sm_matrix_transf_to_m(source_m=full_s(1)%matrix, &
               target_m=glob_angles%min_overlap_m(1)%matrix%sm,&
               glob_angles=glob_angles,qs_env=qs_env,&
               global_env=global_env,error=error)
          glob_angles%m_s_m_valid=.TRUE.
       END IF
       CPAssert(ASSOCIATED(glob_angles%min_overlap_m),cp_failure_level,routineP,error,failure)
       IF (glob_angles%m_s_m_valid) THEN
          min_overlap_m => glob_angles%min_overlap_m
       ELSE
          NULLIFY(min_overlap_m)
       END IF
       !       end if
       IF (PRESENT(min_hamiltonian_m)) THEN
          IF (.not.glob_angles%m_h_m_valid) THEN
             IF (.not.ASSOCIATED(glob_angles%min_hamiltonian_m)) THEN
                ALLOCATE(glob_angles%min_hamiltonian_m(1),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                CALL cp_sp_create(glob_angles%min_hamiltonian_m(1)%matrix, &
                     matrix_struct=min_bas_m_struct,&
                     error=iError)
                CALL cp_error_propagate_error(iError,routineP,&
                     CPSourceFileRef,&
                     error,failure)
             END IF
             IF (.not.failure) THEN
                CALL sm_matrix_transf_to_m(source_m=full_h(1)%matrix,&
                     target_m=glob_angles%min_hamiltonian_m(1)%matrix%sm,&
                     glob_angles=glob_angles,qs_env=qs_env,&
                     global_env=global_env,error=error)
                glob_angles%m_h_m_valid=.TRUE.
             END IF
          END IF
          CPAssert(ASSOCIATED(glob_angles%min_overlap_m),cp_failure_level,routineP,error,failure)
          IF (glob_angles%m_h_m_valid) THEN
             min_hamiltonian_m => glob_angles%min_hamiltonian_m
          ELSE
             NULLIFY(min_hamiltonian_m)
          END IF
       END IF
    END IF
    !call timestop(0.0_dp,handle)
  END SUBROUTINE pao_g_ang_get_min_m
!***************************************************************************

END MODULE pao_proj_methods

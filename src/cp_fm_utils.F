!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/cp_fm_utils *
!!
!!   NAME
!!     cp_fm_utils
!!
!!   FUNCTION
!!     A module to cure the fucked up naming conventions in the whole
!!     blacs / full_matrix mess.
!!
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!****
MODULE cp_fm_utils
  USE cp_error_handling,               ONLY: cp_simple_assert
  USE cp_fm_types
  USE kinds,                           ONLY: dbl,&
                                             dp_size
  USE message_passing,                 ONLY: mp_sum
  USE qs_blacs,                        ONLY: blacs_get_element,&
                                             blacs_replicated_copy
  USE termination,                     ONLY: stop_memory
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: moduleN='cp_fm_utils'
  LOGICAL, PARAMETER :: DEBUG_THIS_MODULE=.FALSE.

  PRIVATE
  PUBLIC :: copy_fm_col_to_rep_vec, &
            copy_rep_vec_to_fm_col, &
            copy_fm_to_rep_vec, &
            copy_rep_vec_to_fm, &
            copy_fm_col_to_fm

CONTAINS

! *****************************************************************************

  SUBROUTINE copy_fm_col_to_rep_vec(mat_a, col_a, vec_b)
  
    IMPLICIT NONE

  	! arguments
    TYPE(cp_full_matrix_type), POINTER :: mat_a
    INTEGER, INTENT(in)                :: col_a
    REAL(dbl), DIMENSION(:), POINTER   :: vec_b
  
  	! locals
    INTEGER :: stat
    INTEGER :: a_row_global, a_col_global, &
               a_row_block, a_col_block, &
               a_row_local, a_col_local
    INTEGER :: row, col
    REAL(dbl) :: element
    INTEGER, DIMENSION(:), POINTER :: a_row_indices, a_col_indices
    CHARACTER(len=*), PARAMETER :: routineN='copy_fm_col_to_rep_vec', &
                                   routineP=moduleN//'/'//routineN

    ! code
    CALL cp_fm_get_info(matrix=mat_a, &
                        nrow_global=a_row_global, ncol_global=a_col_global, &
                        nrow_block=a_row_block,   ncol_block=a_col_block, &
                        nrow_local=a_row_local,   ncol_local=a_col_local, &
                        row_indices=a_row_indices, col_indices=a_col_indices)

    CPAAssert(a_row_global == SIZE(vec_b))
    CPAAssert(col_a >=1 .AND. a_col_global>=col_a)

    vec_b(:) = 0.0_dbl

    DO col=1, a_col_local
       IF (a_col_indices(col)==col_a) THEN
          DO row=1, a_row_local
             vec_b(a_row_indices(row)) = mat_a%local_data(row, col)
          END DO
       END IF
    END DO

#if defined(__SCALAPACK)
    CALL mp_sum(vec_b, mat_a%matrix_struct%para_env%group)
#endif

    IF (DEBUG_THIS_MODULE) THEN
       DO row=1, a_row_global
          CALL blacs_get_element(mat_a, row, col_a, element)
#if defined(__SCALAPACK)
          CALL mp_sum(element, mat_a%matrix_struct%para_env%group)
#endif
          IF (element/=vec_b(row)) THEN
             WRITE (*,*) routineP//': inconsistent data'
             STOP 1
          END IF
       END DO
    END IF

  END SUBROUTINE copy_fm_col_to_rep_vec

! *****************************************************************************

  SUBROUTINE copy_rep_vec_to_fm_col(vec_a, mat_b, col_b)

    IMPLICIT NONE

    ! arguments
    REAL(dbl), DIMENSION(:), POINTER   :: vec_a
    TYPE(cp_full_matrix_type), POINTER :: mat_b
    INTEGER, INTENT(in)                :: col_b

    ! locals
    INTEGER :: b_row_global, b_col_global, &
               b_row_block, b_col_block, &
               b_row_local, b_col_local
    INTEGER :: row, col
    INTEGER, DIMENSION(:), POINTER :: b_row_indices, b_col_indices
    REAL(DBL) :: element
    CHARACTER(len=*), PARAMETER :: routineN='copy_rep_vec_to_fm', &
                                   routineP=moduleN//'/'//routineN

    ! code
    CALL cp_fm_get_info(matrix=mat_b, &
                        nrow_global=b_row_global, ncol_global=b_col_global, &
                        nrow_block=b_row_block,   ncol_block=b_col_block, &
                        nrow_local=b_row_local,   ncol_local=b_col_local, &
                        row_indices=b_row_indices, col_indices=b_col_indices)

    CPAAssert(SIZE(vec_a) == b_row_global)
    CPAAssert(col_b >= 1 .AND. col_b <= b_col_global)

    ! very simple code not very efficient
    DO col=1, b_col_local
          IF (b_col_indices(col) == col_b) THEN
             DO row=1, b_row_local
                mat_b%local_data(row,col) = vec_a(b_row_indices(row))
             END DO
          END IF
    END DO

    IF (DEBUG_THIS_MODULE) THEN
          DO row=1, b_row_global 
             CALL blacs_get_element(mat_b, row, col_b, element)
#if defined(__SCALAPACK)
             CALL mp_sum(element, mat_b%matrix_struct%para_env%group)
#endif
             IF (element/=vec_a(row)) THEN
                WRITE (*,*) routineP//': Inconsistent data'
                STOP 1
             END IF
          END DO
    END IF

  END SUBROUTINE copy_rep_vec_to_fm_col

! *****************************************************************************

  SUBROUTINE copy_fm_to_rep_vec(mat_a, vec_b)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_type), POINTER :: mat_a
    REAL(dbl), DIMENSION(:), POINTER   :: vec_b

    ! locals
    INTEGER :: a_row_global, a_col_global, &
               a_row_block, a_col_block, &
               a_row_local, a_col_local
    INTEGER :: row, col, vec_b_pos, i, j
    INTEGER, DIMENSION(:), POINTER :: a_row_indices, a_col_indices
    REAL(DBL) :: element
    CHARACTER(len=*), PARAMETER :: routineN='copy_fm_to_rep_vec', &
                                   routineP=moduleN//'/'//routineN

    ! code
    CALL cp_fm_get_info(matrix=mat_a, &
                        nrow_global=a_row_global, ncol_global=a_col_global, &
                        nrow_block=a_row_block,   ncol_block=a_col_block, &
                        nrow_local=a_row_local,   ncol_local=a_col_local, &
                        row_indices=a_row_indices, col_indices=a_col_indices)

    CPAAssert(a_row_global*a_col_global==SIZE(vec_b)) ! the vector must fit
    
       DO row=1, a_row_local
    DO col=1, a_col_local
          i = a_row_indices(row)
          j = a_col_indices(col)
          vec_b_pos = (a_row_global*(j-1))+i
          CPAAssert(vec_b_pos>=1.AND.vec_b_pos<=SIZE(vec_b))
          vec_b(vec_b_pos) = mat_a%local_data(row,col)
       END DO
    END DO

#if defined(__SCALAPACK)
    CALL mp_sum(vec_b, mat_a%matrix_struct%para_env%group)
#endif
    
    IF (DEBUG_THIS_MODULE) THEN
       DO col=1, a_col_global
          DO row=1, a_row_global
             CALL blacs_get_element(mat_a, row, col, element)
#if defined(__SCALAPACK)
             CALL mp_sum(element, mat_a%matrix_struct%para_env%group)
#endif
             IF (element/=vec_b(((col-1)*a_row_global)+row)) THEN
                WRITE (*,*) routineP//': Inconsistent data'
                STOP 1
             END IF
          END DO
       END DO
    END IF

  END SUBROUTINE copy_fm_to_rep_vec

! *****************************************************************************

  SUBROUTINE copy_rep_vec_to_fm(vec_a, mat_b)

    IMPLICIT NONE

    ! arguments
    REAL(dbl), DIMENSION(:), POINTER   :: vec_a
    TYPE(cp_full_matrix_type), POINTER :: mat_b

    ! locals
    INTEGER :: b_row_global, b_col_global, &
               b_row_block, b_col_block, &
               b_row_local, b_col_local
    INTEGER :: row, col, vec_a_pos, i, j
    INTEGER, DIMENSION(:), POINTER :: b_row_indices, b_col_indices
    REAL(DBL) :: element
    CHARACTER(len=*), PARAMETER :: routineN='copy_rep_vec_to_fm', &
                                   routineP=moduleN//'/'//routineN

    ! code
    CALL cp_fm_get_info(matrix=mat_b, &
                        nrow_global=b_row_global, ncol_global=b_col_global, &
                        nrow_block=b_row_block,   ncol_block=b_col_block, &
                        nrow_local=b_row_local,   ncol_local=b_col_local, &
                        row_indices=b_row_indices, col_indices=b_col_indices)

    CPAAssert(b_row_global*b_col_global==SIZE(vec_a)) ! the matrix must fit

    ! very simple code not very efficient
       DO row=1, b_row_local
    DO col=1, b_col_local
          i = b_row_indices(row)
          j = b_col_indices(col)
          vec_a_pos = (b_row_global*(j-1)) + i
          CPAAssert(vec_a_pos>=1.AND.vec_a_pos<=SIZE(vec_a))
          mat_b%local_data(row,col) = vec_a(vec_a_pos)
       END DO
    END DO

    IF (DEBUG_THIS_MODULE) THEN
       DO col=1, b_col_global
          DO row=1, b_row_global 
             CALL blacs_get_element(mat_b, row, col, element)
#if defined(__SCALAPACK)
             CALL mp_sum(element, mat_b%matrix_struct%para_env%group)
#endif
             IF (element/=vec_a(((col-1)*b_row_global)+row)) THEN
                WRITE (*,*) routineP//': Inconsistent data'
                STOP 1
             END IF
          END DO
       END DO
    END IF

  END SUBROUTINE copy_rep_vec_to_fm

! *****************************************************************************

  SUBROUTINE copy_fm_col_to_fm(mat_a, col_a, mat_b)
    
    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_type), POINTER :: mat_a
    INTEGER, INTENT(in)                :: col_a
    TYPE(cp_full_matrix_type), POINTER :: mat_b

    ! locals
    INTEGER :: a_row_global, stat
    REAL(dbl), DIMENSION(:), POINTER :: temp_vec
    CHARACTER(len=*), PARAMETER :: routineN='copy_fm_col_to_fm', &
                                   routineP=moduleN//'/'//routineN
    

    CALL cp_fm_get_info(matrix=mat_a, nrow_global=a_row_global)

    ALLOCATE(temp_vec(a_row_global), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, 'temp_vec', a_row_global*dp_size)

    CALL copy_fm_col_to_rep_vec(mat_a, col_a, temp_vec)
    CALL copy_rep_vec_to_fm(temp_vec, mat_b)

    DEALLOCATE(temp_vec, stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, 'temp_vec')

  END SUBROUTINE copy_fm_col_to_fm

! *****************************************************************************

END MODULE cp_fm_utils


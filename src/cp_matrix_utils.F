!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_matrix_utils [1.0] *
!!
!!   NAME
!!     cp_matrix_utils
!!
!!   FUNCTION
!!     various utilities to do operations with sparse matrices
!!     and 2d array seen as dense matrices
!!
!!   NOTES
!!     this should depend only on sparse_matrix_types
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
MODULE cp_matrix_utils
  USE cp_array_utils,                  ONLY: cp_2d_i_guarantee_size,&
                                             cp_2d_logical_guarantee_size,&
                                             cp_2d_r_write
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_propagate_error,&
                                             cp_error_type,&
                                             cp_postcondition_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_struct,                    ONLY: cp_sm_struct_type,&
                                             sm_struct_add_local_block
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             copy_local_sm_to_replicated_fm,&
                                             find_block_node,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_matrix_utils'
! types
  PUBLIC :: cp_matrix_block_iterator

! core methods
  PUBLIC :: cp_init, cp_dealloc_ref, cp_get, cp_next, cp_valid

  PUBLIC :: cp_trace

  PUBLIC :: cp_p_norm_matrix, cp_sup_norm_matrix,&
       cp_is_block_diagonal, cp_init_diagonal_blocks,&
       cp_d_struct_create, cp_set_to
  
  ! underliying methods
  PUBLIC :: cp_sm_diag_trace, cp_sm_set_to, cp_sm_b_i_init,&
       cp_sm_b_i_dealloc_ref, cp_sm_b_i_get, &
       cp_sm_b_i_next, cp_sm_b_i_valid, &
       cp_sm_is_block_diagonal, cp_sm_init_diagonal_blocks,&
       cp_sm_set_dist_sparsity_from_sm,cp_sm_alloc_nonsparse_blocks,&
       cp_sm_write, sm_struct_add_sm_local_blocks

  ! return the trace of the given matrices
  INTERFACE cp_trace
     MODULE PROCEDURE cp_sm_diag_trace
  END INTERFACE

!! set all the values to a given value
  INTERFACE cp_set_to
     MODULE PROCEDURE cp_sm_set_to
  END INTERFACE

! = core methods =

!! initializes the object
  INTERFACE cp_init
     MODULE PROCEDURE cp_sm_b_i_init
  END INTERFACE

!! deallocates the memory allocated by the object
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_sm_b_i_dealloc_ref
  END INTERFACE

!! gets various attributes from the object
  INTERFACE cp_get
     MODULE PROCEDURE cp_sm_b_i_get
  END INTERFACE

!! goes to the next element, returns true if valid
  INTERFACE cp_next
     MODULE PROCEDURE cp_sm_b_i_next
  END INTERFACE

!! return true if it is valid (minimal testing)
  INTERFACE cp_valid
     MODULE PROCEDURE cp_sm_b_i_valid
  END INTERFACE

!! returns true if the matrix is block diagonal
  INTERFACE cp_is_block_diagonal
     MODULE PROCEDURE cp_sm_is_block_diagonal
  END INTERFACE

!! initializes the diagonal matrix blocks
  INTERFACE cp_init_diagonal_blocks
     MODULE PROCEDURE cp_sm_init_diagonal_blocks
  END INTERFACE

!!***
!****************************************************************************

!!****s* cp_matrix_utils/cp_matrix_block_iterator [1.0] *
!!
!!   NAME
!!     cp_matrix_block_iterator
!!
!!   FUNCTION
!!     goes through all the local blocks of the given matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if the structure was initialized (bug catching)
!!     matrix: the actual matrix
!!     block_row: the actual block row (-1 if past end, -2 in case
!!            of error, 0 after a dealloc)
!!     block_node: the actual block
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
  TYPE cp_matrix_block_iterator
     PRIVATE
     LOGICAL :: initialized
     TYPE(real_matrix_type), POINTER :: matrix
     INTEGER :: block_row
     TYPE(real_block_node_type), POINTER :: block_node
  END TYPE cp_matrix_block_iterator
!!***
!****************************************************************************

CONTAINS
! ============== real full matrix utilities ============

!!****f* cp_matrix_utils/cp_p_norm_matrix [1.0] *
!!
!!   NAME
!!     cp_p_norm_matrix
!!
!!   SYNOPSIS
!!     Function cp_p_norm_matrix(matrix, p) Result(res)
!!       Real(KIND = dp):: res
!!       Real(KIND = dp), Intent (IN), Dimension(:,:):: matrix
!!       Integer, Optional, Intent (IN):: p
!!     End Function cp_p_norm_matrix
!!
!!   FUNCTION
!!     calculates the p-te power of the lp norm on the given matrix:
!!     ||a_{i j}||_{p}^p=\sum_{i,j} a_{i j}^p
!!     The frobenious norm is the one with p=2
!!
!!   NOTES
!!     could be much more efficient
!!
!!   INPUTS
!!     matrix: the matrix whose norm is wanted
!!     p: exponent whose power is used to calculate the norm,
!!        defaults to 2
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
!!*** **********************************************************************
  !pure
  FUNCTION cp_p_norm_matrix(matrix, p) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: p
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_matrix_utils:cp_froebnius_norm2'

    INTEGER                                  :: i, j, pVal
    LOGICAL                                  :: failure

!type(cp_error_type), optional, intent(inout) :: error ! removed to make it pure

    failure=.FALSE.

    res=0
    IF (.NOT. failure) THEN
       pVal=2
       IF (.NOT. PRESENT(p)) pVal=p
       DO i=1,SIZE(matrix,1)
          DO j=1,SIZE(matrix,2)
             res=res+matrix(i,j)**pVal
          END DO
       END DO
    END IF
  END FUNCTION cp_p_norm_matrix
!***************************************************************************

!!****f* cp_matrix_utils/cp_sup_norm_matrix [1.0] *
!!
!!   NAME
!!     cp_sup_norm_matrix
!!
!!   SYNOPSIS
!!     Function cp_sup_norm_matrix(matrix) Result(res)
!!       Real(KIND = dp):: res
!!       Real(KIND = dp), Dimension(:,:), Intent (IN):: matrix
!!     End Function cp_sup_norm_matrix
!!
!!   FUNCTION
!!     returns the supremum (maximum) norm of the matrix
!!     ||a_{i j}||_{sup}=max_{i,j}(|a_{i j}|)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix whose norm you want to calculate
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
!!*** **********************************************************************
  !pure 
  FUNCTION cp_sup_norm_matrix(matrix) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: matrix
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_matrix_utils:cp_sup_norm_matrix'

    LOGICAL                                  :: failure

!type(cp_error_type), optional, intent(inout) :: error ! removed to make it pure

    failure=.FALSE.

    IF (.NOT. failure) THEN
       res=MAXVAL(ABS(matrix))
    ELSE
       res=-1.0_dp
    END IF
  END FUNCTION cp_sup_norm_matrix
!***************************************************************************

! ================== cp_matrix_block_iterator ================

!!****f* cp_matrix_utils/cp_sm_b_i_init [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_init
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_b_i_init(block_iterator, matrix, block_row,&
!!         block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (OUT):: block_iterator
!!       Type(real_matrix_type), Target, Intent (IN):: matrix
!!       Integer, Intent (IN), Optional:: block_row
!!       Type(real_block_node_type), Target, Optional:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_b_i_init
!!
!!   FUNCTION
!!     initializes the block iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the block iterator to be initialized
!!     matrix: the matrix you iterate on
!!     block_row: the initial block row (defaults to 0, before the first)
!!     block_node: the initial block (defaults to the first of the row)
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
  SUBROUTINE cp_sm_b_i_init(block_iterator, matrix, block_row,&
       block_node, error)
    TYPE(cp_matrix_block_iterator), &
      INTENT(out)                            :: block_iterator
    TYPE(real_matrix_type), INTENT(in), &
      TARGET                                 :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: block_row
    TYPE(real_block_node_type), OPTIONAL, &
      TARGET                                 :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_init'

    LOGICAL                                  :: failure

! well subsequent operations on the iterator could change tha data in the matrix

    failure=.FALSE.

    block_iterator%matrix => matrix
    block_iterator%block_row = 1
    NULLIFY(block_iterator%block_node)
    block_iterator%initialized=.TRUE.
    IF (PRESENT(block_node)) &
         block_iterator%block_node => block_node
    IF (PRESENT(block_row)) &
         block_iterator%block_row = block_row
  END SUBROUTINE cp_sm_b_i_init
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_b_i_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_b_i_dealloc_ref(block_iterator, error)
!!       Type(cp_matrix_block_iterator), Intent (INOUT):: block_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_b_i_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by the matrix block iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the block iterator you want to free
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
  SUBROUTINE cp_sm_b_i_dealloc_ref(block_iterator, error)
    TYPE(cp_matrix_block_iterator), &
      INTENT(inout)                          :: block_iterator
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_matrix_utils:cp_sm_b_i_dealloc_ref'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(block_iterator%matrix, block_iterator%block_node)
       block_iterator%block_row=0
    END IF
  END SUBROUTINE cp_sm_b_i_dealloc_ref
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_b_i_get [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_get
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_b_i_get(block_iterator, block_row, block_col,&
!!         block_val, finished, matrix, block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (IN):: block_iterator
!!       Integer, Optional, Intent (OUT):: block_row, block_col
!!       Real(KIND = dp), Dimension(:,:), Pointer, Optional:: block_val
!!       Logical, Optional, Intent (OUT):: finished
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Type(real_block_node_type), Optional, Pointer:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_b_i_get
!!
!!   FUNCTION
!!     returns the vaious attrubutes of the iterator
!!
!!   NOTES
!!     also att global indexig (firstRow, firstCol)??
!!
!!   INPUTS
!!     block_iterator: the block iterator you get the info from
!!     block_row: the actual block row
!!     block_col: the actual block column
!!     block_val: the values contained in the block (can be modified, 
!!           and the matrix values are modified)
!!     finished: true if the iterator is at end
!!     matrix: the matrix you are iterating on
!!     block_node: the block node (try not to use this)
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
  SUBROUTINE cp_sm_b_i_get(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)
    TYPE(cp_matrix_block_iterator), &
      INTENT(in)                             :: block_iterator
    INTEGER, INTENT(out), OPTIONAL           :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: block_val
    LOGICAL, INTENT(out), OPTIONAL           :: finished
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_get'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(block_row)) &
            block_row=block_iterator%block_row
       IF (PRESENT(block_col)) THEN
          IF (ASSOCIATED(block_iterator%block_node)) THEN
             CALL get_block_node(block_iterator%block_node, block_col=block_col)
          ELSE
             CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
             block_col=-1
          END IF
       END IF
       IF (PRESENT(block_val)) THEN
          CPPrecondition(block_iterator%block_row>0,cp_failure_level,routineP,error,failure)
          IF (.not.failure) THEN
             IF (ASSOCIATED(block_iterator%block_node)) THEN
                CALL get_block_node(block_iterator%block_node,BLOCK=block_val)
             ELSE
                NULLIFY(block_val)
                CPErrorMessage(cp_failure_level,routineP,"invalid BLOCK",error)
             END IF
          END IF
       END IF
       IF (PRESENT(finished)) &
            finished = block_iterator%block_row <= 0
       IF (PRESENT(matrix)) &
            matrix => block_iterator%matrix
       IF (PRESENT(block_node)) &
            block_node => block_iterator%block_node
    END IF
  END SUBROUTINE cp_sm_b_i_get
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_b_i_valid [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_valid
!!
!!   SYNOPSIS
!!     Function cp_sm_b_i_valid(block_iterator, error) Result(res)
!!       Logical:: res
!!       Type(cp_matrix_block_iterator), Intent (IN):: block_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_b_i_valid
!!
!!   FUNCTION
!!     returns true if the iterator is valid (not past end or invalid).
!!     a fast minimal testing is done, and the errors are written as warnings
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the iterator you get the info from
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
  FUNCTION cp_sm_b_i_valid(block_iterator,error) RESULT(res)
    TYPE(cp_matrix_block_iterator), &
      INTENT(in)                             :: block_iterator
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_matrix_utils:cp_sm_b_i_valid'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.(ASSOCIATED(block_iterator%matrix).AND. block_iterator%block_row>0 ))THEN
          CPAssert(block_iterator%block_row > 0,cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(block_iterator%matrix),cp_warning_level,routineP,error,failure)
       END IF
    END IF
    res=.not.failure
  END FUNCTION cp_sm_b_i_valid
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_b_i_next [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_next
!!
!!   SYNOPSIS
!!     Function cp_sm_b_i_next(block_iterator, block_row, block_col,&
!!         block_val, finished, matrix, block_node, error) Result(res)
!!       Type(cp_matrix_block_iterator), Intent (INOUT):: block_iterator
!!       Integer, Optional, Intent (OUT):: block_row, block_col
!!       Real(KIND = dp), Dimension(:,:), Pointer, Optional:: block_val
!!       Logical, Optional, Intent (OUT):: finished
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Type(real_block_node_type), Optional, Pointer:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_sm_b_i_next
!!
!!   FUNCTION
!!     Moves the iterator to the next matrix block.
!!     Returns true if the iterator is valid (not gone past end)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - block_iterator: the block iterator that you want to move
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     the other arguments are the same as the get function, and
!!     are valid only if the result is true
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_sm_b_i_next(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)RESULT(res)
    TYPE(cp_matrix_block_iterator), &
      INTENT(inout)                          :: block_iterator
    INTEGER, INTENT(out), OPTIONAL           :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: block_val
    LOGICAL, INTENT(out), OPTIONAL           :: finished
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_next'

    INTEGER                                  :: i, nblock_row
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (block_iterator%block_row > 0) THEN
          IF (ASSOCIATED(block_iterator%block_node)) THEN
             block_iterator%block_node => next_block_node(block_iterator%block_node)
          ELSE
             block_iterator%block_node => first_block_node(block_iterator%matrix,&
                  block_iterator%block_row)
          END IF
          IF (.not.ASSOCIATED(block_iterator%block_node)) THEN
             CALL get_matrix_info(block_iterator%matrix,nblock_row=nblock_row)
             DO i=block_iterator%block_row+1,nblock_row
                block_iterator%block_node => &
                     first_block_node(block_iterator%matrix,i)
                IF (ASSOCIATED(block_iterator%block_node)) THEN
                   block_iterator%block_row=i
                   EXIT
                END IF
             END DO
          END IF
          IF (.not.ASSOCIATED(block_iterator%block_node)) THEN
             block_iterator%block_row=-1
          ELSE
             res=.TRUE.
             CALL cp_sm_b_i_get(block_iterator, block_row=block_row,&
                  block_col=block_col, block_val=block_val,&
                  finished=finished, matrix=matrix, block_node=block_node,&
                  error=error)
             RETURN
          END IF
       END IF
       CPPostcondition(.not.ASSOCIATED(block_iterator%block_node),cp_warning_level,routineP,error,failure)
       IF (.NOT. block_iterator%block_row <= 0) THEN
          CPErrorMessage(cp_warning_level,routineP,"inconsistent state in iterator",error)
          block_iterator%block_row=-2
       END IF
    ELSE
       block_iterator%block_row=-2
       NULLIFY(block_iterator%block_node)
    END IF
    res=.FALSE.
  END FUNCTION cp_sm_b_i_next
!***************************************************************************

! ============ real block (sparse) matrix utils ===========

!!****f* cp_matrix_utils/cp_sm_diag_trace [1.0] *
!!
!!   NAME
!!     cp_sm_diag_trace
!!
!!   SYNOPSIS
!!     Function cp_sm_diag_trace(diagonal, matrix, para_env, error)&
!!         Result(res)
!!       Real(KIND = dp):: res
!!       Real(KIND = dp), Dimension(:), Intent (IN):: diagonal
!!       Type(real_matrix_type), Intent (IN), Target:: matrix
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_diag_trace
!!
!!   FUNCTION
!!     calculates the trace of the given matriz with the diagonal matrix
!!     given through his diagonal elements. Missing giagonal elements are
!!     assumed to be 0.
!!     Tr(diag(diagonal) matrix)= sum(diagonal(i)*matrix(i,i))
!!
!!   NOTES
!!     -
!!
!!   INPUTS
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
FUNCTION cp_sm_diag_trace(diagonal, matrix, para_env, error) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: diagonal
    TYPE(real_matrix_type), INTENT(in), &
      TARGET                                 :: matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_diag_trace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j
    INTEGER, DIMENSION(:), POINTER           :: first_col, first_row
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(real_matrix_type), POINTER          :: matrix_ptr

  failure=.FALSE.
  matrix_ptr => matrix
  
  CALL timeset(routineN,'I',"",handle)
  CALL cp_init(iterator, matrix=matrix)
  CALL get_matrix_info(matrix_ptr, first_row=first_row, first_col=first_col)
  res=0.0_dp
  DO
     IF (.not.cp_next(iterator,error=error)) EXIT
     CALL cp_get(iterator, block_row=i,block_col=j, block_val=block_val,&
          error=error)
     IF (i==j) THEN
        CPPrecondition(SIZE(block_val,1)==SIZE(block_val,2),cp_failure_level,routineP,error,failure)
        CPPrecondition(first_col(i)==first_row(i),cp_failure_level,routineP,error,failure)
        DO i=1,MIN(SIZE(block_val,1),SIZE(block_val,2))
           res=res+block_val(i,i)+diagonal(first_col(i)+i)
        END DO
     END IF
  END DO
  CALL mp_sum(res,para_env%group)
  CALL timestop(0.0_dp,handle)
END FUNCTION cp_sm_diag_trace
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_set_to [1.0] *
!!
!!   NAME
!!     cp_sm_set_to
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_set_to(matrix, value, error)
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Real(KIND = dp), Optional, Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_set_to
!!
!!   FUNCTION
!!     initializes the allocated blocks of the real matrix with the given
!!     value (defaults to 0)
!!
!!   NOTES
!!     the iteration looks ugly because this function has much error
!!     handling to test the iterator, it should be:
!!     do; if (.not.cp_next(iterator)) exit
!!       call cp_get(..);..
!!     end do
!!
!!   INPUTS
!!     matrix: the matrix to initialize
!!     value: the new value
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
  SUBROUTINE cp_sm_set_to(matrix, value, error)
    TYPE(real_matrix_type), INTENT(inout), &
      TARGET                                 :: matrix
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_set_to'

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: val
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_matrix_block_iterator)           :: iterator

    failure=.FALSE.
    CALL cp_error_init(iError,template_error=error)
    val=0.0_dp

    IF (.NOT. failure) THEN
       CALL cp_init(iterator,matrix,error=iError)
       CALL cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
    END IF
    IF (.not.failure) THEN
       IF (PRESENT(value)) val=value
       DO
          IF (cp_next(iterator,error=iError)) THEN
             CALL cp_get(iterator,block_val=block_val,error=iError)
             CALL cp_error_propagate_error(iError,routineP,&
                  CPSourceFileRef,&
                  error,failure)
             CPPrecondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                block_val=val
             ELSE
                EXIT
             END IF
          ELSE
             CALL cp_error_propagate_error(iError,routineP,&
                  CPSourceFileRef,&
                  error,failure)
             EXIT
          END IF
       END DO
       CALL cp_dealloc_ref(iterator,error=error)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE cp_sm_set_to
!***************************************************************************

!!****f* cp_matrix_utils/cp_build_d_struct [1.0] *
!!
!!   NAME
!!     cp_build_d_struct
!!
!!   FUNCTION
!!     builds the distribution structure of the diagonal blocks
!!     of the given matrix (the array is allocated and you are responsible
!!     for it)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: a pointer that will point to an array with the
!!            distribution structure of the diagonal elements
!!     matrix: the matrix you get the structure from (if not given
!!            distributed in blocks of equal size between the processors)      
!!     n_blocks: the number of diagonal blocks, if no matrix is given 
!!            this must be given
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE cp_d_struct_create(d_struct, matrix, n_blocks, para_env, error)
    INTEGER, DIMENSION(:), POINTER           :: d_struct
    TYPE(real_matrix_type), INTENT(in), &
      OPTIONAL, TARGET                       :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n_blocks
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_build_d_struct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncols, nmin, nrows, &
                                                rest, stat
    LOGICAL                                  :: failure
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(real_matrix_type), POINTER          :: matrix_ptr

  failure=.FALSE.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  IF (PRESENT(matrix)) THEN
     matrix_ptr => matrix
     CALL get_matrix_info(matrix_ptr, nblock_row=nrows, nblock_col=ncols)
     CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
     ncols=MIN(nrows,ncols)
     IF (PRESENT(n_blocks)) THEN
        CPPrecondition(ncols==n_blocks,cp_warning_level,routineP,error,failure)
     END IF
     IF (.not.failure) THEN
        ! assume every diagonal block is on one (and just one) proc
        ALLOCATE(d_struct(ncols),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        d_struct=0
        CALL cp_init(iterator, matrix=matrix_ptr, error=error)
        DO
           IF (.not.cp_next(iterator)) EXIT
           CALL cp_get(iterator,block_row=i,block_col=j,error=error)
           IF (i==j) THEN
              d_struct(i)=para_env%mepos
           END IF
        END DO
        CALL cp_dealloc_ref(iterator,error=error)
        CALL mp_sum(d_struct,para_env%group)
     END IF
  ELSE
     CPPrecondition(PRESENT(n_blocks),cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        ALLOCATE(d_struct(n_blocks),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (.not.failure) THEN
        nmin= n_blocks/para_env%num_pe
        rest=MOD(n_blocks,para_env%num_pe)
        DO j=0,(rest-1)
           DO i=1,nmin+1
              d_struct(j*(nmin+1)+i)=j
           END DO
        END DO
        DO j=rest,(para_env%num_pe-1)
           DO i=1,nmin+1
              d_struct(rest+j*nmin+i)=j
           END DO
        END DO
     END IF
  END IF
  !call timestop(0.0_dp,handle)
END SUBROUTINE cp_d_struct_create
!***************************************************************************


!!****f* cp_matrix_utils/cp_sm_init_diagonal_blocks [1.0] *
!!
!!   NAME
!!     cp_sm_init_diagonal_blocks
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_init_diagonal_blocks(matrix, para_env, d_struct,&
!!         default_value, diagonal_value, error)
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Optional, Dimension(:), Intent (IN), Target:: d_struct
!!       Real(KIND = dp), Intent (IN), Optional:: default_value,&
!!         diagonal_value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_init_diagonal_blocks
!!
!!   FUNCTION
!!     initialize the block diagonal blocks of the given sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to be initialized
!!     d_struct: the distribution structure of the diagonal blocks
!!     default_value: the default value of the elements (if not given
!!             no initialization is performed)
!!     diagonal_value: the value of the elements on the diagonal
!!             (defaults to default_value)
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
  SUBROUTINE cp_sm_init_diagonal_blocks(matrix, para_env, d_struct, &
       default_value, diagonal_value, error)
    TYPE(real_matrix_type), INTENT(inout), &
      TARGET                                 :: matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL, TARGET                       :: d_struct
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: default_value, diagonal_value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_init_diagonal_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncols, nrows
    INTEGER, DIMENSION(:), POINTER           :: my_d_struct
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val
    TYPE(real_block_node_type), POINTER      :: m_block
    TYPE(real_matrix_type), POINTER          :: matrix_ptr

    matrix_ptr => matrix
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix_ptr),cp_failure_level,routineP,error,failure)
    CALL timeset(routineN,'I',"",handle)
    IF (.not.failure) THEN
       CALL get_matrix_info(matrix_ptr, nblock_row=nrows, nblock_col=ncols)
       CPPreconditionNoFail(nrows==ncols,cp_warning_level,routineP,error)
       ncols=MIN(ncols,nrows)
       IF (PRESENT(d_struct)) THEN
          CPPrecondition(ncols==SIZE(d_struct),cp_failure_level,routineP,error,failure)
          my_d_struct => d_struct
       ELSE
          CALL cp_d_struct_create(my_d_struct, n_blocks=ncols, &
               para_env=para_env, error=error)
       END IF
       CPPostcondition(ASSOCIATED(my_d_struct),cp_failure_level,routineP,error,failure)
       IF (.not.failure) THEN
          DO i=1,SIZE(my_d_struct)
             IF (para_env%mepos==my_d_struct(i)) THEN
                m_block => find_block_node(matrix_ptr,i,i)
                NULLIFY(block_val)
                IF (.not.ASSOCIATED(m_block)) THEN
                   CALL add_block_node(matrix_ptr,i,i,block_val)
                ELSE
                   CALL get_block_node(m_block,BLOCK=block_val)
                END IF
                IF (PRESENT(default_value)) block_val=default_value
                IF (PRESENT(diagonal_value)) THEN
                   FORALL (j=1:MIN(SIZE(block_val,1),SIZE(block_val,2)))
                      block_val(j,j)=diagonal_value
                   END FORALL
                END IF
             ELSE
                m_block => find_block_node(matrix_ptr,i,i)
                IF (ASSOCIATED(m_block)) THEN
                   CPErrorMessage(cp_failure_level,routineP,"diag BLOCK on wrong processor",error)
                   ! call remove_real_block_node(matrix_ptr,i,i) ! to do
                END IF
             END IF
          END DO
       END IF
       IF (cp_debug.and.debug_this_module) THEN ! test if diag
          CALL cp_assert(cp_is_block_diagonal(matrix,para_env),&
               cp_failure_level,&
               cp_assertion_failed,routineP, "matrix is not diagonal in "//&
               CPSourceFileRef,&
               error,failure)
       END IF
    END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE cp_sm_init_diagonal_blocks
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_is_block_diagonal [1.0] *
!!
!!   NAME
!!     cp_sm_is_block_diagonal
!!
!!   SYNOPSIS
!!     Function cp_sm_is_block_diagonal(matrix, para_env, error)&
!!         Result(res)
!!       Logical:: res
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_is_block_diagonal
!!
!!   FUNCTION
!!     returns true if the given matrix is block diagonal
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to check
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
FUNCTION cp_sm_is_block_diagonal(matrix, para_env, error) RESULT(res)
    TYPE(real_matrix_type), INTENT(inout), &
      TARGET                                 :: matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_is_block_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncols, nrows, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: myDiag
    LOGICAL                                  :: failure
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(real_matrix_type), POINTER          :: matrix_ptr

  failure=.FALSE.
  matrix_ptr => matrix

  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(ASSOCIATED(matrix_ptr),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_matrix_info(matrix_ptr, nblock_row=nrows, nblock_col=ncols)
     CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
     ncols=MIN(nrows,ncols)
     ALLOCATE(myDiag(ncols),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     myDiag=0
     CALL cp_init(iterator, matrix=matrix, error=error)
     DO
        IF (.not.cp_next(iterator)) EXIT
        CALL cp_get(iterator,block_row=i,block_col=j,error=error)
        IF (i==j) THEN
           myDiag(i)=myDiag(i)+1
        ELSE
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                routineP,"there is an off diagonal BLOCK! in "//&
                CPSourceFileRef,error)
           myDiag(i)=myDiag(i)+100
           myDiag(j)=myDiag(j)+100
        END IF
     END DO
     CALL cp_dealloc_ref(iterator,error=error)
     CALL mp_sum(myDiag,para_env%group)
     ! there are off diagonal blocks
     CPAssertNoFail(ALL(myDiag<100),cp_warning_level,routineP,error)
     DO i=1,SIZE(myDiag)
        myDiag(i)= MOD(myDiag(i),100)
     END DO
     CPAssert(ALL(myDiag==1),cp_failure_level,routineP,error,failure)
     DEALLOCATE(myDiag, stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)     
  END IF
  res=.not.failure
  !call timestop(0.0_dp,handle)
END FUNCTION cp_sm_is_block_diagonal
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_set_dist_sparsity_from_sm [1.0] *
!!
!!   NAME
!!     cp_sm_set_dist_sparsity_from_sm
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_set_dist_sparsity_from_sm(m_struct, para_env,&
!!         sample_matrix, error)
!!       Type(cp_b_matrix_struct_type), Intent (INOUT):: m_struct
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(real_matrix_type), Pointer:: sample_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_set_dist_sparsity_from_sm
!!
!!   FUNCTION
!!     sets the distribution and sparsity taking them from the given
!!     sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - m_struct: the matrix structure to be initialized
!!     - para_env: the mpi environment
!!     - sample_matrix: the matrix you get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_set_dist_sparsity_from_sm(m_struct, para_env,&
     sample_matrix,error)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(inout)                          :: m_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_type), POINTER          :: sample_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_sm_set_dist_sparsity_from_sm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_col, block_row, i, j, &
                                                nblock_cols, nblock_rows
    INTEGER, DIMENSION(:, :), POINTER        :: b_owner
    LOGICAL                                  :: failure
    TYPE(cp_matrix_block_iterator)           :: b_iterator
    TYPE(real_block_node_type), POINTER      :: block_node

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(sample_matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  IF (.not.failure) THEN
     CALL get_matrix_info(sample_matrix,nblock_row=nblock_rows, &
          nblock_col=nblock_cols)
     CALL cp_2d_i_guarantee_size(m_struct%distribution%block_owner,&
          n_rows=nblock_rows, n_cols=nblock_cols, error=error)
     CALL cp_2d_logical_guarantee_size(m_struct%sparsity,&
          n_rows=nblock_rows, n_cols=nblock_cols, error=error)
     b_owner => m_struct%distribution%block_owner
     b_owner=0
     m_struct%sparsity=.TRUE.
     CALL cp_sm_b_i_init(b_iterator, matrix=sample_matrix, error=error)
     DO
        IF (.not.cp_sm_b_i_next(b_iterator,block_row=block_row, &
             block_col=block_col,error=error)) EXIT
        b_owner(block_row,block_col)=para_env%mepos+1
     END DO
     CALL cp_sm_b_i_dealloc_ref(b_iterator,error=error)
     CALL mp_sum(b_owner,para_env%group)
     DO j=1,SIZE(b_owner,2)
      DO i=1,SIZE(b_owner,1)
        IF (b_owner(i,j)==0) THEN
          b_owner(i,j)= MOD(j,para_env%num_pe)
          m_struct%sparsity(i,j)=.FALSE.
        ELSE IF (b_owner(i,j)==para_env%mepos+1) THEN
          b_owner(i,j)=para_env%mepos
          IF (cp_debug .AND. debug_this_module) THEN
            CALL get_block_node(sample_matrix,block_row=i,&
                block_col=j,block_node=block_node)
            CPPostcondition(ASSOCIATED(block_node),cp_failure_level,routineP,error,failure)
          END IF
        ELSE
          b_owner(i,j)=b_owner(i,j)-1
          CPPostcondition(b_owner(i,j)<para_env%num_pe,cp_failure_level,routineP,error,failure)
        END IF
      END DO
     END DO
  END IF
END SUBROUTINE cp_sm_set_dist_sparsity_from_sm
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_alloc_nonsparse_blocks [1.0] *
!!
!!   NAME
!!     cp_sm_init_nonsparse_blocks
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_set_dist_sparsity_from_sm(m_struct, global_env,&
!!         sample_matrix, error)
!!       Type(cp_b_matrix_struct_type), Intent (INOUT):: m_struct
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(real_matrix_type), Pointer:: sample_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_set_dist_sparsity_from_sm
!!
!!   FUNCTION
!!     Allocate the blocks that are not in the sparsity structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - m_struct: the matrix structure to be initialized
!!     - para_env: the mpi environment
!!     - sample_matrix: the matrix you get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     5.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_alloc_nonsparse_blocks(matrix, m_struct, para_env,&
     error)
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_b_matrix_struct_type), &
      INTENT(inout)                          :: m_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_alloc_nonsparse_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nblock_cols, nblock_rows
    INTEGER, DIMENSION(:, :), POINTER        :: b_owner
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(:, :), POINTER        :: sparsity
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  IF (.not.failure)THEN
     CPPrecondition(ASSOCIATED(m_struct%distribution),cp_failure_level,routineP,error,failure)
     CPPrecondition(m_struct%distribution%initialized,cp_failure_level,routineP,error,failure)
  END IF
  IF (.not.failure) THEN
     sparsity => m_struct%sparsity
     b_owner => m_struct%distribution%block_owner
     CALL get_matrix_info(matrix,nblock_row=nblock_rows, &
          nblock_col=nblock_cols)
     IF (.not.ASSOCIATED(sparsity)) THEN
        CALL cp_assert(.FALSE.,cp_warning_level,fromWhere=routineP,&
             error_nr=cp_postcondition_failed,&
             message="no sparsity, allocating a full matrix "//&
             CPSourceFileRef,&
             error=error)
        DO i=1,nblock_rows
           DO j=1,nblock_cols
              IF (.not.ASSOCIATED(find_block_node(matrix,&
                   block_row=i,block_col=j))) THEN
                 NULLIFY(block_val)
                 CALL add_block_node(matrix,block_row=i,&
                      block_col=j,BLOCK=block_val)
                 block_val=0.0_dp
              END IF
           END DO
        END DO
        RETURN
     END IF
  END IF

  IF (.not.failure) THEN
     DO i=1,nblock_rows
        DO j=1,nblock_cols
           IF (sparsity(i,j).and.b_owner(i,j)==para_env%mepos.and.&
                .not.ASSOCIATED(find_block_node(matrix,&
                block_row=i,block_col=j))) THEN
              NULLIFY(block_val)
              CALL add_block_node(matrix,block_row=i,&
                   block_col=j,BLOCK=block_val)
              block_val=0.0_dp
           END IF
        END DO
     END DO
  END IF
END SUBROUTINE cp_sm_alloc_nonsparse_blocks
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_write [1.0] *
!!
!!   NAME
!!     cp_sm_write
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_write(matrix, unit_nr, long_description, local,&
!!         error)
!!       Type(real_matrix_type), Pointer:: matrix
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: long_description, local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_write
!!
!!   FUNCTION
!!     writes a matrix to the given output unit
!!
!!   NOTES
!!     to do
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_write(matrix,unit_nr,long_description,local,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    INTEGER                                  :: iblock_col, iblock_row
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(real_block_node_type), POINTER      :: block_node

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  my_local=.FALSE.
  my_long_description=.FALSE.
  IF (PRESENT(local)) my_local=local
  IF (PRESENT(long_description)) my_long_description=long_description
  IF (.not.my_local) my_local=(logger%para_env%mepos==logger%para_env%source)
  
  IF (ASSOCIATED(matrix)) THEN
     IF (my_local) THEN
        CALL get_matrix_info(matrix, matrix_name=matrix_name, &
             matrix_symmetry=matrix_symmetry)
        WRITE(unit=unit_nr,fmt="(' <real_matrix_type>{ name=',a)")&
             TRIM(matrix_name)
        WRITE(unit=unit_nr,fmt="(' symmetry=',a)")&
             TRIM(matrix_symmetry)
        DO iblock_row=1,matrix%nblock_row
           block_node => first_block_node(matrix,iblock_row)
           DO WHILE (ASSOCIATED(block_node))

              CALL get_block_node(block_node=block_node,&
                BLOCK=matrix_block,block_col=iblock_col)
              IF (ASSOCIATED(matrix_block)) THEN
                 IF (my_long_description) THEN
                    WRITE (unit=unit_nr,fmt="(' block(',i6,',',i6,')=')")&
                         iblock_row,iblock_col
                    CALL cp_2d_r_write(matrix_block,unit_nr=unit_nr,error=error)
                 ELSE
                    WRITE (unit=unit_nr,&
                         fmt="(' block(',i6,',',i6,')=array(',i6,',',i6,')')")&
                         iblock_row,iblock_col,&
                         SIZE(matrix_block,1),SIZE(matrix_block,2)
                 END IF
              ELSE
                 WRITE (unit=unit_nr,fmt="(' block(',i6,',',i6,')=*null*')")&
                      iblock_row,iblock_col
              END IF
              
              block_node => next_block_node(block_node)
           END DO
        END DO
        WRITE(unit=unit_nr,fmt="(' }')")
        CALL m_flush(unit_nr)
     END IF
  ELSE
     IF (my_local) THEN
        WRITE(unit=unit_nr,fmt="(' <real_matrix_type *null*>')")
        CALL m_flush(unit_nr)
     END IF
  END IF
END SUBROUTINE cp_sm_write
!***************************************************************************

!!****f* cp_matrix_utils/sm_struct_add_sm_local_blocks [1.0] *
!!
!!   NAME
!!     sm_struct_add_sm_local_blocks
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_add_sm_local_blocks(sm_struct, matrix, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Type(real_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_add_sm_local_blocks
!!
!!   FUNCTION
!!     Adds the blocks of the given matrix to the local blocks of the
!!     given matrix structure
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the structure to update
!!     - matrix: the matrix from where to get the blocks
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_add_sm_local_blocks(sm_struct,matrix,error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'sm_struct_add_sm_local_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_col, block_row, handle
    LOGICAL                                  :: failure
    TYPE(cp_matrix_block_iterator)           :: m_iter

  failure=.FALSE.
  
  CALL timeset(routineN,"I","",handle)
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL cp_sm_b_i_init(m_iter, matrix=matrix, error=error)
     
     !
     ! appears to be quadratic scaling 
     !  
     DO WHILE( cp_sm_b_i_next(m_iter, block_row=block_row, &
          block_col=block_col,error=error))
        CALL sm_struct_add_local_block(sm_struct, block_row=block_row,&
             block_col=block_col, error=error)
     END DO
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE sm_struct_add_sm_local_blocks
!***************************************************************************

END MODULE cp_matrix_utils
